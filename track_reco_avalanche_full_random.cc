/**
 * track_reco_avalanche.cc
 * * [統合・拡張版]
 * 1. 複数イベントループ & パラメータ (a, b) の乱数化
 * 2. TTree へのデータ保存 (diff_x0 は絶対値で記録)
 * 3. 100回に1回の頻度でトラッキング可視化画像を出力
 * 4. RANSAC + 反復接線フィット (Iterative Tangent Fit) の維持
 */

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
#include "TFile.h"
#include "TTree.h"

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

// --- geometry choice ---
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
    auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
    const std::string nm = lower(e.name);
    const bool isAnode = nm.find("anode") != std::string::npos;
    const bool isPW    = nm.find("pw")    != std::string::npos;
    for (int n = nL; n <= nR; ++n) {
      const double x = e.x0 + n * pitch; if (x < xmin || x > xmax) continue;
      auto cir = new TEllipse(x, e.y, e.radius, e.radius);
      cir->SetFillStyle(1001);
      if (isAnode) cir->SetFillColor(kOrange + 7);
      else if (isPW) cir->SetFillColor(kAzure + 2);
      else cir->SetFillColor(kBlack);
      cir->SetLineColor(kBlack); cir->SetLineWidth(1); cir->Draw("same");
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

// ==================== L0 Calculation Logic ====================
struct TimeHist {
  std::vector<double> t_center_ns; std::vector<double> dt_ns; std::vector<double> count;
};

static bool LoadTimeHistogram(const char* fname, TimeHist& h) {
  std::ifstream fin(fname); if (!fin) { std::perror(fname); return false; }
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    bool hasAlpha = false; for (char c : line) if (std::isalpha((unsigned char)c)) { hasAlpha = true; break; }
    if (hasAlpha) continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line); double tcen, dt, cnt;
    if (!(ss >> tcen >> dt >> cnt)) continue;
    h.t_center_ns.push_back(tcen); h.dt_ns.push_back(dt); h.count.push_back(cnt);
  }
  return !h.t_center_ns.empty();
}

static double CalculateL0Directly(const TimeHist& h, double gap_cm, double t_hit) {
  if (h.t_center_ns.empty()) return 0.0;
  double numerator = 0.0, denominator = 0.0;
  for (size_t i = 0; i < h.t_center_ns.size(); ++i) {
    denominator += h.count[i];
    if (h.t_center_ns[i] <= t_hit) numerator += h.count[i];
  }
  return (denominator <= 0.0) ? 0.0 : (numerator / denominator) * gap_cm;
}

static fs::path GuessRunDirFromCsv(const fs::path& csvPath) {
  fs::path p = csvPath.parent_path();
  std::string leaf = p.filename().string();
  if (leaf == "analysis_L0_prim" || leaf == "analysis_L0" || leaf == "analysis") p = p.parent_path();
  fs::path q = p;
  while (!q.empty()) {
    if (q.filename().string().rfind("run_", 0) == 0) return q;
    if (!q.has_parent_path()) break; q = q.parent_path();
  }
  return p;
}

// ---- main ----------------------------------------------------
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(1111);

  // --- 1. Simulation Setup ---
  const int nEvents = 1000; // 実行イベント数 (適宜変更)
  TRandom3 rng(0); 

  const char* timeHistFile = (argc > 1) ? argv[1] : "./root/run_20251225_131147/analysis_L0_prim/t_hist_nt.csv";
  #if defined(GEOM_WIRE)
    const char* geomTag = "wirecath";
  #elif defined(GEOM_PLATE)
    const char* geomTag = "plate";
  #endif

  fs::path csvPath(timeHistFile);
  fs::path runDir = GuessRunDirFromCsv(csvPath);
  std::string figDir = (runDir / ("png_" + std::string(geomTag))).string();
  gSystem->mkdir(figDir.c_str(), kTRUE);
  fs::path waveDir = fs::path(figDir) / "waveform";
  fs::create_directory(waveDir);

  auto saveFig = [&](TCanvas& c, const char* fname){ c.SaveAs((figDir + "/" + fname).c_str()); };

  // --- 2. TTree for Statistical Analysis ---
  TFile* rootOut = new TFile((runDir / "reco_analysis.root").c_str(), "RECREATE");
  TTree* tree = new TTree("tree", "Track Reconstruction Multi-Event");
  double br_true_a, br_true_b, br_true_x0;
  double br_reco_a, br_reco_b, br_reco_x0;
  double br_diff_x0; // ABS(Reco - True)
  int br_ev;
  tree->Branch("ev",      &br_ev,      "ev/I");
  tree->Branch("true_a",  &br_true_a,  "true_a/D");
  tree->Branch("true_b",  &br_true_b,  "true_b/D");
  tree->Branch("true_x0", &br_true_x0, "true_x0/D");
  tree->Branch("reco_a",  &br_reco_a,  "reco_a/D");
  tree->Branch("reco_b",  &br_reco_b,  "reco_b/D");
  tree->Branch("reco_x0", &br_reco_x0, "reco_x0/D");
  tree->Branch("diff_x0", &br_diff_x0, "diff_x0/D");

  // Detector Components
  MediumMagboltz gas; gas.LoadGasFile("ic4H10_100_0.1atm.gas");
  Detector::Geometry geo;
#if defined(GEOM_PLATE)
  geo = Detector::PAP_PlaneCathode_Periodic(0.20, 0.50);
#elif defined(GEOM_WIRE)
  geo = Detector::PAP_WireCathode_Periodic(0.20, 0.50);
#endif

  auto comp = Detector::BuildField(geo); comp->SetMedium(&gas);
  Sensor sensor; sensor.AddComponent(comp.get());
  sensor.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, +0.1);
  const char* readoutLabel = "Anode";
  sensor.AddElectrode(comp.get(), readoutLabel);
  sensor.SetTransferFunction(TransferFunction);

  const double gap_val = DetectGap(geo);
  TimeHist th; if (!LoadTimeHistogram(timeHistFile, th)) return 1;

  const double t0_sig_ns = 0.0, dt_sig_ns = 0.5;
  const double tmax_hist = th.t_center_ns.empty() ? 200.0 : th.t_center_ns.back();
  const double tmax_sig = tmax_hist + 150.0;
  const int nSigBins = std::max(200, (int)std::ceil((tmax_sig - t0_sig_ns) / dt_sig_ns));
  sensor.SetTimeWindow(t0_sig_ns, dt_sig_ns, nSigBins);

  // =========================================================================
  //  EVENT LOOP START
  // =========================================================================
  for (int iEv = 0; iEv < nEvents; ++iEv) {
    std::printf("\n>>> Event %d / %d\n", iEv, nEvents - 1);
    br_ev = iEv;

    // --- 1. Track Parameter Randomization ---
    // angle: 33 deg +/- 10 deg, b: 0.15~0.35 cm
    const double angle_deg = (90.0 - 57.0) + rng.Uniform(-10.0, 10.0);
    const double a_true = std::tan(angle_deg * M_PI / 180.0);
    const double b_true = 0.15 + rng.Uniform(0.20);

    const double x_min_trk = geo.xmin, x_max_trk = geo.xmax;
    const int nClusters = 3062;

    AvalancheMC amcTrack; amcTrack.SetSensor(&sensor);
    TryEnableSignalCalculation(amcTrack); amcTrack.UseWeightingPotential(true);

    struct ClusterInfo { double x, y; };
    std::map<int, std::vector<ClusterInfo>> wireHitMap;
    std::map<int, std::pair<double, double>> wireCenterMap;

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
    }

    // --- 2. Signal Processing & Hit Generation ---
    struct HitCand { double xwire, ywire, y_up, y_dn; };
    std::vector<HitCand> hits;
    const double signalThreshold = 0.01;

    for (auto const& [wIdx, clusters] : wireHitMap) {
      sensor.ClearSignal();
      for (const auto& cl : clusters) {
        double xw = wireCenterMap[wIdx].first;
        amcTrack.AvalancheElectron(cl.x - xw, cl.y, 0, 0);
      }
      sensor.ConvoluteSignals();

      // Debug: Waveform save every 100 events
      if (iEv % 100 == 0) {
        TCanvas cw("cw", "", 600, 400); cw.SetGrid();
        std::vector<double> tx(nSigBins), vy(nSigBins);
        for(int k=0; k<nSigBins; ++k) {
          tx[k] = t0_sig_ns + k*dt_sig_ns; vy[k] = sensor.GetSignal(readoutLabel, k, true);
        }
        TGraph gW(nSigBins, tx.data(), vy.data());
        gW.SetTitle(Form("Ev %d Wire %d;Time [ns];Signal", iEv, wIdx));
        gW.Draw("AL");
        TLine *l1 = new TLine(t0_sig_ns, signalThreshold, tmax_sig, signalThreshold);
        l1->SetLineColor(kRed); l1->SetLineStyle(2); l1->Draw("same");
        TLine *l2 = new TLine(t0_sig_ns, -signalThreshold, tmax_sig, -signalThreshold);
        l2->SetLineColor(kRed); l2->SetLineStyle(2); l2->Draw("same");
        cw.SaveAs((waveDir / Form("ev%d_wire_sig_%d.png", iEv, wIdx)).string().c_str());
      }

      double t_hit = -1.0;
      for (int k = 0; k < nSigBins; ++k) {
        if (std::abs(sensor.GetSignal(readoutLabel, k, true)) > signalThreshold) {
          t_hit = t0_sig_ns + k * dt_sig_ns; break;
        }
      }
      if (t_hit >= 0.0) {
        double L_meas = CalculateL0Directly(th, gap_val, t_hit);
        HitCand hc; hc.xwire = wireCenterMap[wIdx].first; hc.ywire = wireCenterMap[wIdx].second;
        hc.y_up = hc.ywire + L_meas; hc.y_dn = hc.ywire - L_meas;
        hits.push_back(hc);
      }
    }

    if (hits.size() < 2) continue;

    // --- 3. Line Fitting (RANSAC + Tangent Fit) ---
    struct CandPt { double x,y; }; std::vector<CandPt> cand;
    for (const auto& h : hits) { cand.push_back({h.xwire, h.y_up}); cand.push_back({h.xwire, h.y_dn}); }

    auto ssr_for_line = [&](double a, double b){
      double ssr = 0, denom = a*a + 1.0;
      for (const auto& h : hits) {
        double nu = a*h.xwire - h.y_up + b, nd = a*h.xwire - h.y_dn + b;
        ssr += std::min(nu*nu, nd*nd) / denom;
      }
      return ssr;
    };

    double best_a=0, best_b=0, best_ssr=1e300;
    for (int it=0; it<2000; ++it) {
      int i1 = rng.Integer(cand.size()), i2 = rng.Integer(cand.size());
      if (std::abs(cand[i1].x - cand[i2].x)<1e-4) continue;
      double a = (cand[i2].y - cand[i1].y)/(cand[i2].x - cand[i1].x), b = cand[i1].y - a*cand[i1].x;
      double ssr = ssr_for_line(a,b); if (ssr < best_ssr) { best_ssr=ssr; best_a=a; best_b=b; }
    }

    double fit_a = best_a, fit_b = best_b;
    for (int iter = 0; iter < 10; ++iter) {
      std::vector<double> xt, yt; double denom = std::sqrt(1.0 + fit_a * fit_a), nx = -fit_a / denom, ny = 1.0 / denom;
      for (const auto& h : hits) {
        double r = std::abs(h.y_up - h.ywire), dist = (fit_a * h.xwire - h.ywire + fit_b) / denom;
        double s = (dist > 0) ? 1.0 : -1.0;
        xt.push_back(h.xwire + s * r * nx); yt.push_back(h.ywire + s * r * ny);
      }
      TGraph gT(xt.size(), xt.data(), yt.data()); TFitResultPtr res = gT.Fit("pol1", "Q S N");
      if (res.Get()) { fit_b = res->Parameter(0); fit_a = res->Parameter(1); }
    }

    // --- 4. Store Results in TTree ---
    br_true_a = a_true; br_true_b = b_true;
    br_true_x0 = (std::abs(a_true) > 1e-9) ? -b_true / a_true : 0.0;
    br_reco_a = fit_a; br_reco_b = fit_b;
    br_reco_x0 = (std::abs(fit_a) > 1e-9) ? -fit_b / fit_a : 0.0;
    
    // 差分の絶対値を計算する
    br_diff_x0 = std::abs(br_reco_x0 - br_true_x0);
    tree->Fill();

    // デバッグ可視化: 100 イベントごとに描画する
    if (iEv % 100 == 0) {
      TCanvas cTrk("cTrk","",1000,780); cTrk.SetGrid();
      TH2D fr("fr", Form("Event %d: a_true=%.3f, b_true=%.3f;X [cm];Y [cm]", iEv, a_true, b_true), 
              10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
      fr.SetStats(0); fr.Draw(); DrawElectrodes(geo);
      
      TF1 fT("fT","[0]*x+[1]", geo.xmin, geo.xmax); fT.SetParameters(a_true, b_true);
      fT.SetLineColor(kBlue); fT.SetLineStyle(2); fT.Draw("same");
      
      TF1 fR("fR","[0]*x+[1]", geo.xmin, geo.xmax); fR.SetParameters(fit_a, fit_b);
      fR.SetLineColor(kRed); fR.SetLineWidth(2); fR.Draw("same");
      
      for (const auto& h : hits) {
        double L = std::abs(h.y_up - h.ywire); TEllipse *el = new TEllipse(h.xwire, h.ywire, L, L);
        el->SetFillStyle(0); el->SetLineColor(kGreen+2); el->Draw("same");
      }
      TLatex lx; lx.SetNDC(); lx.SetTextSize(0.03);
      lx.DrawLatex(0.15, 0.85, Form("True x0: %.4f", br_true_x0));
      lx.DrawLatex(0.15, 0.81, Form("Reco x0: %.4f", br_reco_x0));
      lx.DrawLatex(0.15, 0.77, Form("Abs Diff: %.4f cm", br_diff_x0));
      saveFig(cTrk, Form("debug_track_reco_ev%d.png", iEv));
    }
  }
  // =========================================================================
  //  ANALYSIS & OUTPUT
  // =========================================================================

  TCanvas cRes("cRes", "Resolution Analysis", 800, 600);
  // Absolute difference distribution (Half-normal shape)
  TH1D* hDiff = new TH1D("hDiff", "Position Error Magnitude; |X_{reco} - X_{true}| [cm]; Entries", 100, 0, 0.05);
  tree->Draw("diff_x0 >> hDiff");
  hDiff->SetFillColor(kCyan-9);
  hDiff->Draw();
  saveFig(cRes, "final_resolution_summary.png");

  tree->Write(); rootOut->Close();
  std::printf("\nDone. Results saved to reco_analysis.root and PNGs.\n");
  return 0;
}