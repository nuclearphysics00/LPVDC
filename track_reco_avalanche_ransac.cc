/**
 * track_reco_avalanche.cc
 *
 * [最終版: 全機能統合 + 10イベントループ統計出力]
 * 1. L0(t)計算: テーブル補間を行わず、ヒストグラムを直接積分して算出
 * 2. トラッキング: RANSAC初期推定 + 反復接線フィット (Iterative Tangent Fit)
 * 3. 可視化(1): 最終結果に直線の式(y=ax+b)とX切片(y=0でのx)を表示
 * 4. 可視化(2): 波形確認グラフに負の閾値線も表示
 * 5. 統計処理: 10イベント実行し、真値との差分をROOTファイルに出力
 */

// #define DEBUG_MODE 

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <numeric>
#include <filesystem>
#include <type_traits>
#include <utility>

// ROOT includes
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TLatex.h" 
#include "TFitResult.h"
#include "TMatrixDSym.h"
#include "TFile.h"  // 追加
#include "TTree.h"  // 追加

// Garfield++ includes
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh" 

// Local includes
#include "detector_geometry.hh"
#include "geometry/eval_uniform.hh"
#include "eval_maps.hh"

// --- geometry choice --------------------------------------------------------
#if defined(GEOM_PLATE)
  #include "geometry/geom_pap_plate.hh"
#elif defined(GEOM_WIRE)
  #include "geometry/geom_pap_wirecath.hh"
#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

using namespace Garfield;
namespace fs = std::filesystem;

// ==================== Preamplifier Response ====================
double TransferFunction(double t) {
  constexpr double tau = 20.0; 
  constexpr double gain = 1.0; 
  if (t < 0.0) return 0.0;
  return gain * (t / tau) * std::exp(1.0 - t / tau);
}

// ===== Helpers ==========================
template <typename T, typename = void>
struct has_EnableSignalCalculation : std::false_type {};
template <typename T>
struct has_EnableSignalCalculation<T, std::void_t<decltype(std::declval<T&>().EnableSignalCalculation())>>
  : std::true_type {};

static void TryEnableSignalCalculation(Garfield::AvalancheMC& av) {
  if constexpr (has_EnableSignalCalculation<Garfield::AvalancheMC>::value) {
    av.EnableSignalCalculation();
  }
}

struct WireNearest {
  double d_surface_cm;  
  double r_cm;          
  double xw_cm, yw_cm;  
  bool   hasWire;       
};

static void DrawElectrodes(const Detector::Geometry& g) {
  const double xmin = g.xmin, xmax = g.xmax;
  for (const auto& e : g.electrodes) {
    if (e.kind == Detector::ElectrodeKind::PlaneY) {
      auto ln = new TLine(xmin, e.y, xmax, e.y);
      ln->SetLineColor(kGray + 3); ln->SetLineStyle(2); ln->SetLineWidth(3);
      ln->Draw("same");
      continue;
    }
    const double pitch = g.pitchX;
    int nL = int(std::floor((xmin - e.x0) / pitch)) - 1;
    int nR = int(std::ceil ((xmax - e.x0) / pitch)) + 1;

    auto lower = [](std::string s){
      std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s;
    };
    const std::string nm = lower(e.name);
    const bool isAnode = nm.find("anode") != std::string::npos;
    const bool isPW    = nm.find("pw")    != std::string::npos;

    for (int n = nL; n <= nR; ++n) {
      const double x = e.x0 + n * pitch; if (x < xmin || x > xmax) continue;
      auto cir = new TEllipse(x, e.y, e.radius, e.radius);
      cir->SetFillStyle(1001);
      if      (isAnode) cir->SetFillColor(kOrange + 7);
      else if (isPW)    cir->SetFillColor(kAzure  + 2);
      else              cir->SetFillColor(kBlack);
      cir->SetLineColor(kBlack);
      cir->SetLineWidth(1);
      cir->Draw("same");
    }
  }
}

static double DetectGap(const Detector::Geometry& g) {
  double gap = 0.0;
  for (const auto& e : g.electrodes) gap = std::max(gap, std::abs(e.y));
  if (gap <= 0.0) gap = std::min(std::abs(g.ymax), std::abs(g.ymin));
  return gap;
}

static WireNearest NearestWireSurface(const Detector::Geometry& g, double x, double y, const std::string& targetLabel){
  WireNearest out; out.d_surface_cm = 1e300; out.r_cm = 0.0; out.xw_cm = 0.0; out.yw_cm = 0.0; out.hasWire=false;
  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    if (e.name.find(targetLabel) == std::string::npos) continue;

    out.hasWire = true;
    const double k  = std::round((x - e.x0) / g.pitchX);
    const double xw = e.x0 + k * g.pitchX;
    const double d  = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) { out.d_surface_cm = d; out.r_cm = e.radius; out.xw_cm = xw; out.yw_cm = e.y; }
  }
  return out;
}

static void AddSeedsUniformOnLine(const Detector::Geometry& g, double x0, double y0, double ux, double uy, double spacing, std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs) {
  const double L = std::hypot(ux, uy); if (L==0) return;
  ux/=L; uy/=L;
  const double xmin=g.xmin, xmax=g.xmax, ymin=g.ymin, ymax=g.ymax;
  for(double s=-2.0; s<2.0; s+=spacing) {
      double px = x0 + s*ux, py = y0 + s*uy;
      if(px>xmin && px<xmax && py>ymin && py<ymax) {
          xs.push_back(px); ys.push_back(py); zs.push_back(0);
      }
  }
}

// ==================== L0 Calculation Logic ====================
struct TimeHist {
  std::vector<double> t_center_ns;
  std::vector<double> dt_ns;
  std::vector<double> count;
};

static bool LoadTimeHistogram(const char* fname, TimeHist& h) {
  std::ifstream fin(fname);
  if (!fin) { std::perror(fname); return false; }
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    bool hasAlpha = false; for (char c : line) if (std::isalpha((unsigned char)c)) { hasAlpha = true; break; }
    if (hasAlpha) continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line);
    double tcen, dt, cnt;
    if (!(ss >> tcen >> dt >> cnt)) continue;
    h.t_center_ns.push_back(tcen); h.dt_ns.push_back(dt); h.count.push_back(cnt);
  }
  return !h.t_center_ns.empty();
}

// ヒストグラムの生データを走査し、t_hit までの積分(分子)と全体の積分(分母)を毎回計算する
static double CalculateL0Directly(const TimeHist& h, double gap_cm, double t_hit) {
  if (h.t_center_ns.empty()) return 0.0;
  
  double numerator = 0.0;   // 分子: ∫t_min^t N(t') dt'
  double denominator = 0.0; // 分母: ∫t_min^t_max N(t') dt'

  for (size_t i = 0; i < h.t_center_ns.size(); ++i) {
    double t_bin = h.t_center_ns[i];
    double count = h.count[i];
    
    denominator += count;
    if (t_bin <= t_hit) {
      numerator += count;
    }
  }

  if (denominator <= 0.0) return 0.0;
  return (numerator / denominator) * gap_cm;
}

static fs::path GuessRunDirFromCsv(const fs::path& csvPath) {
  fs::path p = csvPath.parent_path();
  std::string leaf = p.filename().string();
  if (leaf == "analysis_L0_prim" || leaf == "analysis_L0" || leaf == "analysis") p = p.parent_path();
  fs::path q = p;
  while (!q.empty()) {
    if (q.filename().string().rfind("run_", 0) == 0) return q;
    if (!q.has_parent_path()) break;
    q = q.parent_path();
  }
  return p;
}

// ---- main ----------------------------------------------------
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const char* timeHistFile = (argc > 1) ? argv[1] : "./root/run_20251225_131147/analysis_L0_prim/t_hist_nt.csv";
  
  #if defined(GEOM_WIRE)
    const char* geomTag = "wirecath";
  #elif defined(GEOM_PLATE)
    const char* geomTag = "plate";
  #endif
  fs::path csvPath(timeHistFile);
  fs::path runDir  = GuessRunDirFromCsv(csvPath);
  std::string figDir = (runDir / ("png_" + std::string(geomTag))).string();
  const char* p_fig = std::getenv("FIG_DIR");
  if (p_fig && std::strlen(p_fig)) figDir = p_fig;
  gSystem->mkdir(figDir.c_str(), kTRUE);

  fs::path waveDir = fs::path(figDir) / "waveform";
  fs::create_directory(waveDir);

  // --- Output for Track Stats (New) ---
  fs::path trackDir = runDir / "track";
  fs::create_directories(trackDir);
  std::string rootOutPath = (trackDir / "track_results.root").string();
  
  TFile* fOut = new TFile(rootOutPath.c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "Track Reconstruction Stats");
  double t_x_true, t_x_reco, t_diff, t_abs_diff;
  double t_a_true, t_b_true, t_a_reco, t_b_reco;
  int    t_event_id;
  
  tree->Branch("event_id", &t_event_id);
  tree->Branch("a_true", &t_a_true);
  tree->Branch("b_true", &t_b_true);
  tree->Branch("x_true", &t_x_true);
  tree->Branch("a_reco", &t_a_reco);
  tree->Branch("b_reco", &t_b_reco);
  tree->Branch("x_reco", &t_x_reco);
  tree->Branch("diff_x", &t_diff);
  tree->Branch("abs_diff", &t_abs_diff);

  auto saveFig = [&](TCanvas& c, const char* fname){
    c.SaveAs((figDir + "/" + fname).c_str());
  };

  std::printf("[env] timeHistFile=%s\n", timeHistFile);
  std::printf("[env] figDir=%s\n", figDir.c_str());
  std::printf("[env] rootOutPath=%s\n", rootOutPath.c_str());

  // Gas & Geometry
  MediumMagboltz gas;
  gas.LoadGasFile("ic4H10_100_0.1atm.gas");
  Detector::Geometry geo;
#if defined(GEOM_PLATE)
  geo = Detector::PAP_PlaneCathode_Periodic(0.20, 0.50);
#elif defined(GEOM_WIRE)
  geo = Detector::PAP_WireCathode_Periodic(0.20, 0.50);
#endif

  auto comp = Detector::BuildField(geo);
  comp->SetMedium(&gas);

  Sensor sensor;
  sensor.AddComponent(comp.get());
  sensor.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, +0.1);

  const char* readoutLabel = std::getenv("READOUT_LABEL");
  if (!readoutLabel || !std::strlen(readoutLabel)) readoutLabel = "Anode";
  sensor.AddElectrode(comp.get(), readoutLabel);
  sensor.SetTransferFunction(TransferFunction);

  // L0 Data Loading
  const double gap_val = DetectGap(geo);
  TimeHist th;
  if (!LoadTimeHistogram(timeHistFile, th)) return 1;

  // Signal Window
  const double t0_sig_ns = 0.0;
  const double dt_sig_ns = 0.5; 
  const double tmax_hist = th.t_center_ns.empty() ? 200.0 : th.t_center_ns.back();
  const double tmax_sig  = tmax_hist + 150.0; 
  const int nSigBins = std::max(200, (int)std::ceil((tmax_sig - t0_sig_ns) / dt_sig_ns));
  sensor.SetTimeWindow(t0_sig_ns, dt_sig_ns, nSigBins);

  std::printf("[setup] Window: 0-%.1fns (%d bins). Label='%s'\n", tmax_sig, nSigBins, readoutLabel);

  // ===== (A) field map (Static, loop independent) =====
  {
    const int NX=200, NY=160;
    TH2D hE("|E|","|E|;X [cm];Y [cm]", NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    TH2D hV("V","V;X [cm];Y [cm]",     NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    double minPos=1e98, maxVal=0.0;
    for (int ix=1; ix<=NX; ++ix) {
      const double x = hE.GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy<=NY; ++iy) {
        const double y = hE.GetYaxis()->GetBinCenter(iy);
        double ex=0,ey=0,ez=0,V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(x,y,0,ex,ey,ez,V,m,status);
        const double E = std::hypot(ex,std::hypot(ey,ez));
        hE.SetBinContent(ix,iy,E); hV.SetBinContent(ix,iy,V);
        if (E>0 && E<minPos) minPos=E;
        if (E>maxVal) maxVal=E;
      }
    }
    TCanvas c("cField","E map",1000,780);
    c.SetRightMargin(0.14); c.SetGrid(); c.SetLogz();
    if (minPos<1e97) hE.SetMinimum(minPos*0.9);
    hE.SetMaximum(maxVal*1.05);
    hE.Draw("COLZ");
    hV.SetContour(40); hV.SetLineColor(kGray+2); hV.SetLineWidth(2);
    hV.Draw("CONT3 SAME");
    ViewField vf; vf.SetCanvas(&c); vf.SetComponent(comp.get());
    vf.SetArea(geo.xmin, geo.ymin, geo.xmax, geo.ymax);
    std::vector<double> xs,ys,zs;
    const double th_angle = +45.0 * M_PI / 180.0;
    AddSeedsUniformOnLine(geo, 0.0, 0.0, std::cos(th_angle), std::sin(th_angle), 0.05, xs,ys,zs);
    if (!xs.empty()) vf.PlotFieldLines(xs,ys,zs,true,true);
    DrawElectrodes(geo);
    saveFig(c, "fieldmap_combo.png");
  }

  // ===== (B) Track Reconstruction (Event Summation Logic) Loop =====
  if (!geo.electrodes.empty()) {
    
    // シード固定またはランダム化
    TRandom3 gen(0); 

    const int nEvents = 10;
    for (int ev = 0; ev < nEvents; ++ev) {
        std::printf("\n>>> Processing Event %d / %d <<<\n", ev + 1, nEvents);
        
        // パラメータ設定 (イベントごとに少し揺らす)
        const double angle_deg = 90.0 - 57.0;
        const double a_true = std::tan(angle_deg * M_PI / 180.0);
        // 切片を0.25中心に±0.05程度揺らす
        const double b_true = 0.25 + gen.Uniform(-0.05, 0.05); 
        
        const double x_min_trk = geo.xmin;
        const double x_max_trk = geo.xmax;
        const int    nClusters = 3062;

        AvalancheMC amcTrack;
        amcTrack.SetSensor(&sensor);
        TryEnableSignalCalculation(amcTrack);
        amcTrack.UseWeightingPotential(true);

        // Drift View はイベントごとに初期化
        ViewDrift driftView;
        driftView.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, 0.1);

        // --- 1. 全クラスター生成 & ワイヤー仕分け ---
        struct ClusterInfo { double x, y; };
        std::map<int, std::vector<ClusterInfo>> wireHitMap;
        std::map<int, std::pair<double, double>> wireCenterMap; 

        std::printf("[Action] Generating %d clusters and grouping by wire...\n", nClusters);

        for (int i = 0; i < nClusters; ++i) {
          const double frac = (i + 0.5) / nClusters;
          const double x = x_min_trk + frac * (x_max_trk - x_min_trk);
          const double y = a_true * x + b_true;

          if (y < geo.ymin || y > geo.ymax) continue;
          
          auto nw = NearestWireSurface(geo, x, y, readoutLabel);
          if (!nw.hasWire) continue;

          int wireIdx = std::round(nw.xw_cm / geo.pitchX);
          wireHitMap[wireIdx].push_back({x, y});
          wireCenterMap[wireIdx] = {nw.xw_cm, nw.yw_cm};
          
          if (i % 50 == 0) {
            amcTrack.EnablePlotting(&driftView);
            amcTrack.AvalancheElectron(x, y, 0, 0); 
            amcTrack.DisablePlotting();
          }
        }

        std::printf("[Action] Processing %zu active wires (Signal Summation)...\n", wireHitMap.size());

        // --- 2. 信号積算 & ヒット生成 ---
        struct HitCand { double xwire, ywire, y_up, y_dn; };
        std::vector<HitCand> hits;
        
        const double signalThreshold = 0.01; 
        TCanvas cWave("cWave", "", 600, 400);

        for (auto const& [wIdx, clusters] : wireHitMap) {
            sensor.ClearSignal();
            
            for (const auto& cl : clusters) {
                double xw = wireCenterMap[wIdx].first;
                double dx = cl.x - xw; 
                amcTrack.AvalancheElectron(dx, cl.y, 0, 0);
            }

            sensor.ConvoluteSignals();

            // Waveform Saving
            {
              cWave.Clear(); cWave.SetGrid();
              std::vector<double> tx(nSigBins), vy(nSigBins);
              for(int k=0; k<nSigBins; ++k) {
                tx[k] = t0_sig_ns + k*dt_sig_ns;
                vy[k] = sensor.GetSignal(readoutLabel, k, true);
              }
              TGraph gW(nSigBins, tx.data(), vy.data());
              std::string title = "Wire Index " + std::to_string(wIdx) + " (Ev " + std::to_string(ev) + ");Time [ns];Signal";
              gW.SetTitle(title.c_str());
              gW.SetLineWidth(1); gW.SetLineColor(kBlue);
              gW.Draw("AL");

              // ★修正: 正負両方の閾値線を描画する
              TLine *lth_p = new TLine(t0_sig_ns, signalThreshold, tmax_sig, signalThreshold);
              lth_p->SetLineColor(kRed); lth_p->SetLineStyle(2); lth_p->Draw("same");

              TLine *lth_n = new TLine(t0_sig_ns, -signalThreshold, tmax_sig, -signalThreshold);
              lth_n->SetLineColor(kRed); lth_n->SetLineStyle(2); lth_n->Draw("same");

              char wname[64]; std::snprintf(wname, 64, "wire_sig_ev%d_%d.png", ev, wIdx);
              cWave.SaveAs((waveDir / wname).string().c_str());
            }

            double t_hit = -1.0;
            for (int k = 0; k < nSigBins; ++k) {
              if (std::abs(sensor.GetSignal(readoutLabel, k, true)) > signalThreshold) {
                t_hit = t0_sig_ns + k * dt_sig_ns;
                break; 
              }
            }

            if (t_hit >= 0.0) {
                // 直接積分で計算
                double L_meas = CalculateL0Directly(th, gap_val, t_hit);
                if (L_meas < 0) L_meas = 0;
                
                double xw = wireCenterMap[wIdx].first;
                double yw = wireCenterMap[wIdx].second;

                HitCand hc;
                hc.xwire = xw; hc.ywire = yw;
                hc.y_up  = yw + L_meas;
                hc.y_dn  = yw - L_meas;
                hits.push_back(hc);
                
                std::printf("   Wire %d: Hit at t=%.2f ns -> L=%.4f cm (Direct Calc)\n", 
                            wIdx, t_hit, L_meas);
            } else {
                std::printf("   Wire %d: No signal above threshold\n", wIdx);
            }
        }

        std::printf("[Result] Total Wires Hit: %zu\n", hits.size());

        // --- 可視化 (Drift Lines) ---
        TCanvas cClEnd("cClEnd", "cluster endpoint tracks", 1100, 850);
        cClEnd.SetGrid();
        TH2D frClEnd("frClEnd", "Electron Drift Lines;X [cm];Y [cm]", 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
        frClEnd.SetStats(0); frClEnd.Draw();
        DrawElectrodes(geo);
        driftView.SetCanvas(&cClEnd);
        driftView.Plot(true, false);
        
        TF1 fTrueDrift("fTrueDrift","[0]*x+[1]", geo.xmin, geo.xmax);
        fTrueDrift.SetParameters(a_true, b_true);
        fTrueDrift.SetLineColor(kBlue); fTrueDrift.SetLineStyle(7); fTrueDrift.SetLineWidth(2);
        fTrueDrift.Draw("same");
        
        saveFig(cClEnd, Form("cluster_endtracks_ev%d.png", ev));

        if (hits.size() < 2) {
          std::printf("[WARN] not enough hits for track fit in event %d\n", ev);
          continue; // スキップ
        }

        // =========================================================================
        //  Line Fitting (RANSAC + Iterative Tangent Fit)
        // =========================================================================

        // --- 1. RANSACによる初期推定 ---
        struct CandPt { double x,y; };
        std::vector<CandPt> cand; 
        for (const auto& h : hits) {
          cand.push_back({h.xwire, h.y_up});
          cand.push_back({h.xwire, h.y_dn});
        }

        auto ssr_for_line = [&](double a, double b)->double{
          const double denom = a*a + 1.0;
          double ssr = 0.0;
          for (const auto& h : hits) {
            double nu = a*h.xwire - h.y_up + b;
            double nd = a*h.xwire - h.y_dn + b;
            ssr += std::min(nu*nu, nd*nd) / denom;
          }
          return ssr;
        };

        double best_a=0, best_b=0, best_ssr=1e300;
        TRandom3 rng(12345);
        int n_ransac = 2000;
        for (int it=0; it<n_ransac; ++it) {
          int i1 = rng.Integer(cand.size());
          int i2 = rng.Integer(cand.size());
          if (std::abs(cand[i1].x - cand[i2].x)<1e-4) continue;
          double a = (cand[i2].y - cand[i1].y)/(cand[i2].x - cand[i1].x);
          double b = cand[i1].y - a*cand[i1].x;
          double ssr = ssr_for_line(a,b);
          if (ssr < best_ssr) { best_ssr=ssr; best_a=a; best_b=b; }
        }
        
        std::printf("[RANSAC] Initial: y = %.3fx + %.3f\n", best_a, best_b);

        // --- 2. Iterative Tangent Fit (接線フィット) ---
        double fit_a = best_a;
        double fit_b = best_b;
        const int max_iter = 10; 

        for (int iter = 0; iter < max_iter; ++iter) {
            std::vector<double> x_tan_list;
            std::vector<double> y_tan_list;

            // 上向き法線ベクトル n = (-a, 1) / sqrt(1+a^2)
            double denom = std::sqrt(1.0 + fit_a * fit_a);
            double nx = -fit_a / denom;
            double ny =  1.0   / denom;

            for (const auto& h : hits) {
                double r = std::abs(h.y_up - h.ywire);
                // 直線の位置判定 (符号付き距離)
                double dist_signed = (fit_a * h.xwire - h.ywire + fit_b) / denom;
                double sign = (dist_signed > 0) ? 1.0 : -1.0;

                // 接点座標 P_tan = P_wire + sign * r * n
                double tx = h.xwire + sign * r * nx;
                double ty = h.ywire + sign * r * ny;

                x_tan_list.push_back(tx);
                y_tan_list.push_back(ty);
            }

            // 接点群に対してフィット
            if (x_tan_list.size() < 2) break;
            TGraph gTemp(x_tan_list.size(), x_tan_list.data(), y_tan_list.data());
            TFitResultPtr res = gTemp.Fit("pol1", "Q S N"); 
            
            if (res.Get()) {
                fit_b = res->Parameter(0);
                fit_a = res->Parameter(1);
            }
        }
        
        // 最終結果を更新
        best_a = fit_a;
        best_b = fit_b;

        std::printf("[true] y = %.3fx + %.3f\n", a_true, b_true);
        std::printf("[reco] y = %.3fx + %.3f (Tangent Fit)\n", best_a, best_b);

        // --- 3. 描画 ---
        TCanvas cTrk("cTrk","track reco",1000,780);
        cTrk.SetGrid();
        TH2D frame("frame", "Reconstructed Track;X;Y", 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
        frame.SetStats(0); frame.Draw();
        DrawElectrodes(geo);

        TF1 fTrue("fTrue","[0]*x+[1]", geo.xmin, geo.xmax);
        fTrue.SetParameters(a_true, b_true);
        fTrue.SetLineColor(kBlue); fTrue.SetLineStyle(2); fTrue.Draw("same");

        // 円描画
        for (const auto& h : hits) {
            double L = std::abs(h.y_up - h.ywire);
            TEllipse *el = new TEllipse(h.xwire, h.ywire, L, L);
            el->SetFillStyle(0);
            el->SetLineColor(kGreen+2);
            el->SetLineWidth(1);
            el->Draw("same");

            TMarker *m = new TMarker(h.xwire, h.ywire, 20);
            m->SetMarkerSize(0.5); m->SetMarkerColor(kBlack);
            m->Draw("same");
        }

        TF1 fReco("fReco","[0]*x+[1]", geo.xmin, geo.xmax);
        fReco.SetParameters(best_a, best_b);
        fReco.SetLineColor(kRed); fReco.SetLineWidth(2);
        fReco.Draw("same");

        // ★ 計算された「接点」の確認プロット (赤点)
        {
            double denom = std::sqrt(1.0 + best_a * best_a);
            double nx = -best_a / denom;
            double ny =  1.0   / denom;

            for (const auto& h : hits) {
                double r = std::abs(h.y_up - h.ywire);
                double dist = (best_a * h.xwire - h.ywire + best_b) / denom;
                double sign = (dist > 0) ? 1.0 : -1.0;
                
                double tx = h.xwire + sign * r * nx;
                double ty = h.ywire + sign * r * ny;

                TMarker *mt = new TMarker(tx, ty, 20);
                mt->SetMarkerColor(kRed); 
                mt->SetMarkerSize(0.8);
                mt->Draw("same");
            }
        }
        
        // 直線の式とX切片を画像内に書き込む
        TLatex latex;
        latex.SetNDC(); 
        latex.SetTextSize(0.035);
        latex.SetTextAlign(13); // 左・上寄せ

        double x_intercept_true = (std::abs(a_true) > 1e-9) ? -b_true / a_true : 0.0;
        double x_intercept_reco = (std::abs(best_a) > 1e-9) ? -best_b / best_a : 0.0;
        double abs_diff_val = std::abs(x_intercept_true - x_intercept_reco);

        // 真の値 (青)
        latex.SetTextColor(kBlue);
        latex.DrawLatex(0.15, 0.85, Form("True: y = %.4fx + %.4f", a_true, b_true));
        latex.DrawLatex(0.15, 0.81, Form("      x(y=0) = %.4f", x_intercept_true)); 

        // 再構成値 (赤)
        latex.SetTextColor(kRed);
        latex.DrawLatex(0.15, 0.75, Form("Reco: y = %.4fx + %.4f", best_a, best_b));
        latex.DrawLatex(0.15, 0.71, Form("      x(y=0) = %.4f", x_intercept_reco)); 
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.15, 0.65, Form("Diff: |x_true - x_reco| = %.4f cm", abs_diff_val));

        saveFig(cTrk, Form("track_reco_geometry_ev%d.png", ev));

        // --- ROOT Tree への書き込み ---
        t_event_id = ev;
        t_a_true = a_true;
        t_b_true = b_true;
        t_x_true = x_intercept_true;
        t_a_reco = best_a;
        t_b_reco = best_b;
        t_x_reco = x_intercept_reco;
        t_diff = x_intercept_true - x_intercept_reco;
        t_abs_diff = abs_diff_val;
        
        tree->Fill();

    } // end event loop

    // ファイル保存
    tree->Write();
    fOut->Close();
    std::printf("\n[Done] All events processed. Root file saved to: %s\n", rootOutPath.c_str());

  }
  return 0;
}