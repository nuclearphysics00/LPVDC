// fieldview_array.cc  (PBS array-friendly; random seeding version)
// - PBS_ARRAY_INDEX または SHARD_ID で指定されたシャードのみ処理
// - 環境変数 OUT_ROOT に出力 ROOT ファイル名（必須推奨）
// - N_SHARDS 未指定時は 1（イベントを shard 分割）
// - 検出器ジオメトリ内でランダムにシード点をばらまき、AvalancheMC で drift
// - 各イベントについて
//      * drift_all     : アバランシェで増えた全 endpoint の到達時刻 tb
//      * drift_primary : そのイベントで最初に到達した電子の時刻 t_first
//   を TTree として保存
// - メタデータは TTree "meta"
// - 電場マップ/δヒートマップ等の画像は MAKE_MAPS=1 かつ SHARD_ID==0 のときだけ作成
// - HALO_MULT（既定 3.0）で「ハロー（ワイヤ表面から）距離閾」を定義（記録のみ）

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>
#include <random>
#include <cstdint>

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
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/FundamentalConstants.hh"  


#include "detector_geometry.hh"
#include "geometry/eval_uniform.hh"
#include "eval_maps.hh"

#if defined(GEOM_PLATE)
  #include "geometry/geom_pap_plate.hh"
#elif defined(GEOM_WIRE)
  #include "geometry/geom_pap_wirecath.hh"
#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

using namespace Garfield;

// ---- helpers -------------------------------------------------

static void DrawElectrodes(const Detector::Geometry& g) {
  const double xmin = g.xmin, xmax = g.xmax;
  for (const auto& e : g.electrodes) {
    if (e.kind == Detector::ElectrodeKind::PlaneY) {
      auto ln = new TLine(xmin, e.y, xmax, e.y);
      ln->SetLineColor(kGray + 3);
      ln->SetLineStyle(2);
      ln->SetLineWidth(3);
      ln->Draw("same");
      continue;
    }
    const double pitch = g.pitchX;
    int nL = int(std::floor((xmin - e.x0) / pitch)) - 1;
    int nR = int(std::ceil ((xmax - e.x0) / pitch)) + 1;

    auto lower = [](std::string s){
      std::transform(s.begin(), s.end(), s.begin(), ::tolower);
      return s;
    };
    const std::string nm = lower(e.name);
    const bool isAnode = nm.find("anode") != std::string::npos;
    const bool isPW    = nm.find("pw")    != std::string::npos;

    for (int n = nL; n <= nR; ++n) {
      const double x = e.x0 + n * pitch;
      if (x < xmin || x > xmax) continue;
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
  for (const auto& e : g.electrodes) {
    gap = std::max(gap, std::abs(e.y));
  }
  if (gap <= 0.0) {
    gap = std::min(std::abs(g.ymax), std::abs(g.ymin));
  }
  return gap;
}

// 任意直線上シード（等間隔）を少数だけ置いてフィールドライン可視化用に
static void AddSeedsUniformOnLine(
    const Detector::Geometry& g,
    double x0, double y0, double ux, double uy,
    double spacing,
    std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs,
    double margin=1.0e-3) {
  const double L = std::hypot(ux, uy);
  if (L == 0) return;
  ux /= L;
  uy /= L;
  const double xmin = g.xmin + margin, xmax = g.xmax - margin;
  const double ymin = g.ymin + margin, ymax = g.ymax - margin;
  auto s_interval = [](double a0,double da,double amin,double amax){
    if (std::abs(da) < 1e-30) {
      if (a0 < amin || a0 > amax) {
        return std::pair<double,double>(+1e300,-1e300);
      }
      return std::pair<double,double>(-1e300,+1e300);
    }
    double s1=(amin-a0)/da, s2=(amax-a0)/da;
    if (s1>s2) std::swap(s1,s2);
    return std::pair<double,double>(s1,s2);
  };
  auto sx = s_interval(x0,ux,xmin,xmax);
  auto sy = s_interval(y0,uy,ymin,ymax);
  double smin = std::max(sx.first, sy.first);
  double smax = std::min(sx.second,sy.second);
  if (!(smin < smax)) return;
  const double usable = smax - smin;
  const int N = std::max(1, (int)std::floor(usable/spacing));
  for (int i=0;i<N;++i){
    double s = smin + (i+0.5)*(usable/N);
    xs.push_back(x0+s*ux);
    ys.push_back(y0+s*uy);
    zs.push_back(0.0);
  }
}

// 名前フィールドを追加した拡張構造体
struct WireNearest {
  double d_surface_cm;   // 表面までの距離（負なら「食い込み」）
  double r_cm;           // そのワイヤ半径
  double xw_cm, yw_cm;   // そのワイヤ中心
  std::string name;      // 名前フィールドを追加
  bool   hasWire;        // 幾何にワイヤがあるか
};

// 最寄りワイヤを特定した際、設計上の名前を格納する
static WireNearest NearestWireSurface(const Detector::Geometry& g, double x, double y){
  WireNearest out;
  out.d_surface_cm = 1e300;
  out.r_cm = 0.0;
  out.xw_cm = out.yw_cm = 0.0;
  out.name = ""; // 名前を初期化
  out.hasWire = false;
  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    out.hasWire = true;
    const double k  = std::round((x - e.x0) / g.pitchX);
    const double xw = e.x0 + k * g.pitchX;
    const double d  = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) {
      out.d_surface_cm = d;
      out.r_cm = e.radius;
      out.xw_cm = xw;
      out.yw_cm = e.y;
      out.name = e.name; // 名前をコピー
    }
  }
  return out;
}

// ---- main ----------------------------------------------------

int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);

  // ROOT style
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetNumberContours(255);
  gROOT->ForceStyle();

  // --- env (array) ---
  const char* p_shards = std::getenv("N_SHARDS");
  const char* p_shard  = std::getenv("PBS_ARRAY_INDEX");
  const char* p_shard2 = std::getenv("SHARD_ID");
  const char* p_out    = std::getenv("OUT_ROOT");
  const char* p_maps   = std::getenv("MAKE_MAPS");
  const char* p_halo   = std::getenv("HALO_MULT");   // ハロー係数（半径×係数）
  const char* p_nevt   = std::getenv("N_EVENTS");    // 全シャード合計イベント数
  const char* p_seed   = std::getenv("SEED_BASE");   // 乱数シードベース

  const int N_SHARDS = p_shards ? std::max(1, std::atoi(p_shards)) : 1;
  int SHARD_ID = -1;
  if (p_shard && std::strlen(p_shard)) SHARD_ID = std::atoi(p_shard);
  else if (p_shard2 && std::strlen(p_shard2)) SHARD_ID = std::atoi(p_shard2);
  else SHARD_ID = 0; // 非アレイ時の既定

  const double HALO_MULT = (p_halo && std::strlen(p_halo)) ? std::atof(p_halo) : 3.0;
  const long long N_EVENTS_TOTAL = p_nevt ? std::atoll(p_nevt) : 10000LL;

  unsigned long long baseSeed = p_seed ? std::strtoull(p_seed, nullptr, 10)
                                       : 0x9e3779b97f4a7c15ULL;
  unsigned long long seed = baseSeed + (unsigned long long)SHARD_ID * 1000003ULL;
  std::mt19937_64 rng(seed);

  char outname[512];
  if (p_out && std::strlen(p_out)) {
    std::snprintf(outname, sizeof(outname), "%s", p_out);
  } else {
    std::snprintf(outname, sizeof(outname), "grid_times.shard%04d.root", SHARD_ID);
  }

  // --- Decide figure directory (next to OUT_ROOT) ---
  std::string outpath = outname;
  auto slashPos = outpath.find_last_of('/');
  std::string baseDir = (slashPos == std::string::npos)
                      ? std::string(".") : outpath.substr(0, slashPos);

#if defined(GEOM_WIRE)
  const char* geomTag = "wirecath";
#elif defined(GEOM_PLATE)
  const char* geomTag = "plate";
#else
  const char* geomTag = "unknown";
#endif

  std::string figDir = baseDir + "/png_" + std::string(geomTag);

  // allow override by FIG_DIR
  const char* p_fig = std::getenv("FIG_DIR");
  if (p_fig && std::strlen(p_fig)) figDir = p_fig;

  // mkdir -p
  gSystem->mkdir(figDir.c_str(), kTRUE);

  // helper to save into figDir
  auto saveFig = [&](TCanvas& c, const char* fname){
    std::string p = figDir + "/" + fname;
    c.SaveAs(p.c_str());
  };

  const bool make_maps = (p_maps && std::atoi(p_maps) != 0) && (SHARD_ID == 0);

  std::printf("[env] N_SHARDS=%d  SHARD_ID=%d  OUT_ROOT=%s  MAKE_MAPS=%d  HALO_MULT=%.3f  N_EVENTS_TOTAL=%lld  seed=0x%llx\n",
              N_SHARDS, SHARD_ID, outname, (int)make_maps, HALO_MULT,
              N_EVENTS_TOTAL, (unsigned long long)seed);

  // --- Gas ---
  MediumMagboltz gas;
  gas.LoadGasFile("ic4H10_100_0.1atm.gas");
  // Garfield のバージョンによっては Torr が無いので自前定義
  constexpr double Torr = 133.322368;  // [Pa] 1 Torr ≒ 133.322 Pa

  std::printf("[gas] P = %.4f atm, T = %.2f K\n",
              gas.GetPressure(),
              gas.GetTemperature());
  // --- Geometry ---
  Detector::Geometry geo;
#if defined(GEOM_PLATE)
  std::puts("[build] geometry = PAP_PlaneCathode_Periodic");
  geo = Detector::PAP_PlaneCathode_Periodic(/*pitch=*/0.60, /*gap=*/0.50);
#elif defined(GEOM_WIRE)
  std::puts("[build] geometry = PAP_WireCathode_Periodic");
  geo = Detector::PAP_WireCathode_Periodic(/*pitch=*/0.20, /*gap=*/0.50);
#endif

  auto comp = Detector::BuildField(geo);
  comp->SetMedium(&gas);

  // --- Sensor for can_drift ---
  Sensor sensor;
  sensor.AddComponent(comp.get());
  sensor.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, +0.1);

  // ===== can_drift(): ドリフト可否 + ハロー判定（記録用） =====
  const bool debug = (gSystem->Getenv("DRIFT_DEBUG")
                      && std::atoi(gSystem->Getenv("DRIFT_DEBUG")) != 0);
  struct DriftCheck {
    bool   driftable;
    bool   in_halo;
    double d_wire_cm;
    double r_wire_cm;
  };
  auto can_drift = [&](double x, double y)->DriftCheck{
    double Ex=0, Ey=0, Ez=0, V=0;
    Garfield::Medium* m=nullptr;
    int status=0;
    comp->ElectricField(x, y, 0.0, Ex, Ey, Ez, V, m, status);
    const bool in_area = sensor.IsInArea(x, y, 0.0);
    WireNearest wn = NearestWireSurface(geo, x, y);
    double dpos = wn.d_surface_cm;
    if (!std::isfinite(dpos)) dpos = 1e300;
    const bool in_halo =
      (wn.hasWire && HALO_MULT > 0.0 && dpos >= 0.0 && dpos < HALO_MULT * wn.r_cm);
    const bool driftable = (m != nullptr) && (status == 0) && in_area;
    if (debug) {
      const double Enorm = std::sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
      std::printf("[chk] (%.4f, %.4f): status=%d m=%p |E|=%.3g  in_area=%d  d_wire=%.4g cm  r=%.4g cm  in_halo=%d\n",
                  x, y, status, (void*)m, Enorm, (int)in_area, dpos, wn.r_cm, (int)in_halo);
    }
    return {driftable, in_halo, dpos, wn.r_cm};
  };

  // ===== (A) 電場マップ（SHARD 0 かつ MAKE_MAPS=1 のときのみ） =====
  if (make_maps) {
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
        hE.SetBinContent(ix,iy,E);
        hV.SetBinContent(ix,iy,V);
        if (E>0 && E<minPos) minPos=E;
        if (E>maxVal) maxVal=E;
      }
    }
    // --- |E| マップ ---
    {
      TCanvas cE("cE","|E| map",900,720);
      cE.SetRightMargin(0.14);
      cE.SetGrid();
      cE.SetLogz();
      if (minPos < 1e97) hE.SetMinimum(minPos * 0.9);
      hE.SetMaximum(maxVal * 1.05);
      hE.SetTitle("|E|;X [cm];Y [cm]");
      hE.Draw("COLZ");
      saveFig(cE, "emap_absE.png");
    }
    // --- ポテンシャルマップ ---
    {
      TCanvas cV("cV","V map",900,720);
      cV.SetRightMargin(0.14);
      cV.SetGrid();
      hV.SetTitle("Potential V;X [cm];Y [cm]");
      hV.Draw("COLZ");
      saveFig(cV, "emap_potential.png");
    }

    // --- E マップ + 等電位 + field lines ---
    TCanvas c("cField","E map + equipotentials + field lines + electrodes",1000,780);
    c.SetRightMargin(0.14);
    c.SetGrid();
    c.SetLogz();
    if (minPos<1e97) hE.SetMinimum(minPos*0.9);
    hE.SetMaximum(maxVal*1.05);
    hE.Draw("COLZ");
    hV.SetContour(40);
    hV.SetLineColor(kGray+2);
    hV.SetLineWidth(2);
    hV.Draw("CONT3 SAME");

    ViewField vf; vf.SetCanvas(&c); vf.SetComponent(comp.get());
    vf.SetArea(geo.xmin, geo.ymin, geo.xmax, geo.ymax);
    std::vector<double> xs,ys,zs;
    const double th = +45.0 * M_PI / 180.0;
    AddSeedsUniformOnLine(geo, 0.0, 0.0, std::cos(th), std::sin(th), 0.05, xs,ys,zs);
    if (!xs.empty()) vf.PlotFieldLines(xs,ys,zs,true,true);
    for (size_t i=0;i<xs.size();++i){
      TMarker m(xs[i],ys[i],20);
      m.SetMarkerColor(kRed);
      m.SetMarkerSize(0.4);
      m.Draw("same");
    }
    DrawElectrodes(geo);
    saveFig(c, "fieldmap_combo.png");

    // ===== δ マップ & ヒスト =====
    auto deltas = MakeDeltaMapsAndHists(*comp, geo,
                                        /*nx=*/220, /*ny=*/160,
                                        /*halo_mult=*/5.0);

    { TCanvas cD("c_delta","delta heatmap (No-Halo)",900,720);
      cD.SetRightMargin(0.14);
      cD.SetLogz();
      deltas.d_noh->SetTitle("#delta (No-Halo);X [cm];Y [cm]");
      deltas.d_noh->GetZaxis()->SetTitle("relative deviation");
      Eval::prepare_logz(deltas.d_noh);
      deltas.d_noh->Draw("COLZ");
      saveFig(cD, "delta_heatmap_nohalo.png"); }

    { TCanvas cDY("c_deltaY","delta_y heatmap (No-Halo)",900,720);
      cDY.SetRightMargin(0.14);
      cDY.SetLogz();
      deltas.dY_noh->SetTitle("#delta_{y} (No-Halo);X [cm];Y [cm]");
      deltas.dY_noh->GetZaxis()->SetTitle("relative deviation");
      Eval::prepare_logz(deltas.dY_noh);
      deltas.dY_noh->Draw("COLZ");
      saveFig(cDY, "deltaY_heatmap_nohalo.png"); }

    { TCanvas ch("c_delta_hist","delta histograms",1100,520);
      ch.Divide(2,1);
      ch.cd(1); gPad->SetLogy();
      deltas.h_dY_all->SetLineColor(kBlue+2);
      deltas.h_dY_all->Draw("HIST");
      deltas.h_dY_noh->SetLineColor(kRed+1);
      deltas.h_dY_noh->Draw("HIST SAME");
      ch.cd(2); gPad->SetLogy();
      deltas.h_d_all->SetLineColor(kBlue+2);
      deltas.h_d_all->Draw("HIST");
      deltas.h_d_noh->SetLineColor(kRed+1);
      deltas.h_d_noh->Draw("HIST SAME");
      saveFig(ch, "delta_hists.png"); }

    // === Pretty field lines (y = ±yfrac*gap) ===
    {
      const int    FL_N    = gSystem->Getenv("FL_N")    ? std::max(8, std::atoi(gSystem->Getenv("FL_N"))) : 64;
      const double yfrac   = gSystem->Getenv("FL_YFRAC")? std::atof(gSystem->Getenv("FL_YFRAC")) : 0.98;
      const bool   both    = gSystem->Getenv("FL_BOTH") ? (std::atoi(gSystem->Getenv("FL_BOTH")) != 0) : true;
      const double edge_eps2 = 1e-3;
      const double gapv     = DetectGap(geo);
      const double yTop     = +yfrac * gapv;
      const double yBot     = -yfrac * gapv;

      auto inMedium = [&](double x, double y){
        double ex=0,ey=0,ez=0,V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(x,y,0.0,ex,ey,ez,V,m,status);
        return (m!=nullptr && status==0);
      };

      auto make_line_seeds = [&](double yseed, std::vector<double>& xs2,
                                 std::vector<double>& ys2, std::vector<double>& zs2){
        const double Lx = (geo.xmax - edge_eps2) - (geo.xmin + edge_eps2);
        for (int i = 0; i < FL_N; ++i) {
          const double x = geo.xmin + edge_eps2 + (i + 0.5) * Lx / FL_N;
          if (!inMedium(x, yseed)) continue;
          xs2.push_back(x);
          ys2.push_back(yseed);
          zs2.push_back(0.0);
        }
      };

      TCanvas cFL("cFieldLines","Field lines from y = ±yfrac*gap", 1000, 780);
      cFL.SetGrid();
      TH2D frame("fl_frame","Field lines; x [cm]; y [cm]",
                 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
      frame.SetStats(0);
      frame.Draw();

      ViewField vfFL; vfFL.SetCanvas(&cFL); vfFL.SetComponent(comp.get());
      vfFL.SetArea(geo.xmin, geo.ymin, geo.xmax, geo.ymax);

      std::vector<double> xs2, ys2, zs2;
      make_line_seeds(yTop, xs2, ys2, zs2);
      if (both) make_line_seeds(yBot, xs2, ys2, zs2);

      if (!xs2.empty()) vfFL.PlotFieldLines(xs2, ys2, zs2, true, true);
      DrawElectrodes(geo);
      char fname[128];
      std::snprintf(fname, sizeof(fname), "fieldlines_yfrac_%03d%s.png",
                    int(std::round(yfrac*100)), both? "_both":"");
      saveFig(cFL, fname);
    }
  }

  // ===== (B-0) テスト: y=0.05 cm 上のドリフトタイムヒスト（SHARD 0 のみ） =====
  if (SHARD_ID == 0) {
    const double ytest = 0.05;
    const int    Ntest = 200;
    const double edge_eps_test = 1e-3;

    std::printf("[test] Drift electrons on line y = %.4f cm (make histogram)\n", ytest);

    TH1D hT_y005("hT_y005",
                 "Drift time at y=0.05 cm; t [ns]; entries",
                 100, 0.0, 400.0);

    const double Lx_all = (geo.xmax - edge_eps_test) - (geo.xmin + edge_eps_test);

    Sensor sensorTest;
    sensorTest.AddComponent(comp.get());
    sensorTest.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, +0.1);

    AvalancheMC amcTest;
    amcTest.SetSensor(&sensorTest);

    auto inGasTest = [&](double x, double y){
      double ex=0,ey=0,ez=0,V=0; Medium* m=nullptr; int status=0;
      comp->ElectricField(x,y,0.0,ex,ey,ez,V,m,status);
      return (m!=nullptr && status==0);
    };
    auto nudgeTest = [&](double& x, double& y)->bool{
      if (inGasTest(x,y)) return true;
      const double dx = 0.25 * (geo.xmax - geo.xmin);
      const double dy = 0.25 * (geo.ymax - geo.ymin);
      const double ox[9] = {0, +dx, -dx, 0, 0, +dx, +dx, -dx, -dx};
      const double oy[9] = {0, 0, 0, +dy, -dy, +dy, -dy, +dy, -dy};
      for (int k=1;k<9;++k){
        const double xx = x + ox[k], yy = y + oy[k];
        if (!inGasTest(xx,yy)) continue;
        x=xx; y=yy; return true;
      }
      return false;
    };

    for (int i = 0; i < Ntest; ++i) {
      double x = geo.xmin + edge_eps_test + (i + 0.5) * Lx_all / Ntest;
      double y = ytest;

      if (!nudgeTest(x, y)) {
        std::printf("[test]  i=%3d  (x=%.4f, y=%.4f) : cannot nudge into gas, skip.\n",
                    i, x, y);
        continue;
      }

      auto chk = can_drift(x, y);
      if (!chk.driftable) {
        std::printf("[test]  i=%3d  (x=%.4f, y=%.4f) : not driftable\n", i, x, y);
        continue;
      }

      amcTest.DriftElectron(x, y, 0.0, 0.0);
      const size_t nEnd = amcTest.GetNumberOfElectronEndpoints();
      if (!nEnd) {
        std::printf("[test]  i=%3d  (x=%.4f, y=%.4f) : no endpoints\n", i, x, y);
        continue;
      }

      double xa, ya, za, ta;
      double xb, yb, zb, tb;
      int status = 0;
      amcTest.GetElectronEndpoint(nEnd - 1, xa, ya, za, ta, xb, yb, zb, tb, status);
      const double t_ns = tb - ta;

      hT_y005.Fill(t_ns);

      auto nw = NearestWireSurface(geo, x, y);

      std::printf("[test]  i=%3d  start=(%.4f, %.4f)cm  "
                  "end=(%.4f, %.4f)cm  t=%.3f ns  status=%d  "
                  "d_wire=%.4g cm  r_wire=%.4g cm  in_halo=%d\n",
                  i, x, y, xb, yb, t_ns, status,
                  nw.d_surface_cm, nw.r_cm, (int)chk.in_halo);
    }

    TCanvas cT("cT_y005","Drift time at y=0.05 cm",800,600);
    cT.SetGrid();
    hT_y005.SetLineWidth(2);
    hT_y005.Draw("HIST");
    saveFig(cT, "test_y005_t_hist.png");
    std::printf("[test] saved histogram PNG: test_y005_t_hist.png\n");
  }

  // ===== (B) ランダムサンプリング → AvalancheMC → endpoint 全部 + t_first =====
  const double edge_eps = 1e-3;
  const double gap_val  = DetectGap(geo);
  const double xminU = geo.xmin + edge_eps;
  const double xmaxU = geo.xmax - edge_eps;
  const double yminU = -gap_val + edge_eps;
  const double ymaxU = +gap_val - edge_eps;

  auto inGas = [&](double x, double y){
    double ex,ey,ez,V; Medium* m=nullptr; int status=0;
    comp->ElectricField(x,y,0.0,ex,ey,ez,V,m,status);
    return (m!=nullptr && status==0);
  };
  auto nudge_to_gas = [&](double& x, double& y)->bool{
    if (inGas(x,y)) return true;
    const double dx = (xmaxU - xminU);
    const double dy = (ymaxU - yminU);
    const double ox[9] = {0, +0.25*dx, -0.25*dx, 0, 0,
                          +0.25*dx, +0.25*dx, -0.25*dx, -0.25*dx};
    const double oy[9] = {0, 0, 0, +0.25*dy, -0.25*dy,
                          +0.25*dy, -0.25*dy, +0.25*dy, -0.25*dy};
    for (int k=1;k<9;++k){
      const double xx = x + ox[k], yy = y + oy[k];
      if (xx<=xminU || xx>=xmaxU || yy<=yminU || yy>=ymaxU) continue;
      if (inGas(xx,yy)) { x=xx; y=yy; return true; }
    }
    return false;
  };

  std::uniform_real_distribution<double> distX(xminU, xmaxU);
  std::uniform_real_distribution<double> distY(yminU, ymaxU);

  Sensor sensorRnd;
  sensorRnd.AddComponent(comp.get());
  sensorRnd.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, +0.1);

  AvalancheMC amcRnd;
  amcRnd.SetSensor(&sensorRnd);

  struct RowAll {
    double xwire;
    double t_ns;    // 各 endpoint の到達時刻
    int    in_halo;
    double dwire_cm;
  };
  struct RowPrim {
    double xwire;
    double t0_ns;   // そのイベントの最初の到達時刻
    int    in_halo;
    double dwire_cm;
  };

  std::vector<RowAll>  rows_all;
  std::vector<RowPrim> rows_prim;
  rows_all.reserve( (size_t)((N_EVENTS_TOTAL + N_SHARDS - 1) / N_SHARDS) + 1024 );
  rows_prim.reserve( rows_all.capacity() );

  long long visited_evt = 0, saved_endpoints = 0;

  std::printf("[rnd] sampling x=[%.4f, %.4f] cm  y=[%.4f, %.4f] cm\n",
              xminU, xmaxU, yminU, ymaxU);
  std::printf("[rnd] N_EVENTS_TOTAL (all shards) = %lld\n", N_EVENTS_TOTAL);

  for (long long iev = 0; iev < N_EVENTS_TOTAL; ++iev) {
    if ((iev % N_SHARDS) != SHARD_ID) continue; // shard 分割
    ++visited_evt;

    double x = distX(rng);
    double y = distY(rng);

    if (!nudge_to_gas(x, y)) continue;

    auto chk = can_drift(x, y);
    if (!chk.driftable) continue;

    amcRnd.DriftElectron(x, y, 0.0, 0.0);
    const size_t nEnd = amcRnd.GetNumberOfElectronEndpoints();
    if (!nEnd) continue;


    if (visited_evt < 20) {
      std::printf("[dbg] iev=%lld  nEnd=%zu\n", iev, nEnd);
    }

    WireNearest nw = NearestWireSurface(geo, x, y);

    // --- (1) このイベントの t_first を求める (Anode到達のみ対象) ---
    double t_first = 1e300;
    for (size_t j = 0; j < nEnd; ++j) {
      double xa,ya,za,ta, xb,yb,zb,tb; int st=0;
      amcRnd.GetElectronEndpoint(j, xa,ya,za,ta, xb,yb,zb,tb, st);
      
      // 終点 (xb, yb) から最寄りワイヤを特定し名前を確認する
      WireNearest nw_end = NearestWireSurface(geo, xb, yb);
      
      // 最寄りワイヤが "Anode" (大文字小文字問わず) である場合のみ記録
      if (nw_end.hasWire && (nw_end.name.find("Anode") != std::string::npos || nw_end.name.find("anode") != std::string::npos)) {
          if (tb < t_first) t_first = tb;
      }
    }
    // アノードに一つも届かなかったイベントはスキップ
    if (!(t_first < 1e299)) continue;

    // --- 1次電子用：このイベントについて1行だけ保存 ---
    rows_prim.emplace_back(RowPrim{
      nw.xw_cm,
      t_first,
      chk.in_halo ? 1 : 0,
      nw.d_surface_cm
    });

    // --- (2) アバランシェ込み：endpoint 全部を保存 (Anode到達のみ対象) ---
    for (size_t j = 0; j < nEnd; ++j) {
      double xa,ya,za,ta, xb,yb,zb,tb; int st=0;
      amcRnd.GetElectronEndpoint(j, xa,ya,za,ta, xb,yb,zb,tb, st);

      // アバランシェ内の各電子についても Anode 到達を確認する
      WireNearest nw_end = NearestWireSurface(geo, xb, yb);
      
      if (nw_end.hasWire && (nw_end.name.find("Anode") != std::string::npos || nw_end.name.find("anode") != std::string::npos)) {
          const double t_ns = tb; // 初期時刻 0 なので tb が到達時刻

          rows_all.emplace_back(RowAll{
            nw.xw_cm,
            t_ns,
            chk.in_halo ? 1 : 0,
            nw.d_surface_cm
          });
          ++saved_endpoints;
      }
    }
  }

  std::printf("[rnd] shard %d: visited=%lld  endpoints_saved=%lld  prim_events=%zu  -> writing %s\n",
              SHARD_ID, visited_evt, saved_endpoints, rows_prim.size(), outname);

  // ===== ROOT output =====
  TFile fout(outname, "RECREATE");
  fout.SetCompressionSettings(207);

  // --- (1) アバランシェ込み: 全 endpoint 用 ---
  TTree t_all("drift_all","random-seed drift endpoints (per shard)");
  int    br_shard = SHARD_ID, br_nshards = N_SHARDS;
  double br_xwire = 0.0, br_tns = 0.0, br_dwire = 0.0;
  int    br_inhalo = 0;

  t_all.Branch("shard_id", &br_shard,   "shard_id/I");
  t_all.Branch("n_shards", &br_nshards, "n_shards/I");
  t_all.Branch("xwire",    &br_xwire,   "xwire/D");
  t_all.Branch("t_ns",     &br_tns,     "t_ns/D");
  t_all.Branch("in_halo",  &br_inhalo,  "in_halo/I");
  t_all.Branch("dwire_cm", &br_dwire,   "dwire_cm/D");

  for (const auto& r : rows_all) {
    br_xwire  = r.xwire;
    br_tns    = r.t_ns;
    br_inhalo = r.in_halo;
    br_dwire  = r.dwire_cm;
    t_all.Fill();
  }

  // --- (2) 1次電子のみ: 各イベント t_first 用 ---
  TTree t_prim("drift_primary","primary-like first-arrival times (per seed)");
  double br_xwire0 = 0.0, br_t0ns = 0.0, br_dwire0 = 0.0;
  int    br_inhalo0 = 0;

  t_prim.Branch("shard_id", &br_shard,    "shard_id/I");
  t_prim.Branch("n_shards", &br_nshards,  "n_shards/I");
  t_prim.Branch("xwire",    &br_xwire0,   "xwire/D");
  t_prim.Branch("t0_ns",    &br_t0ns,     "t0_ns/D");
  t_prim.Branch("in_halo",  &br_inhalo0,  "in_halo/I");
  t_prim.Branch("dwire_cm", &br_dwire0,   "dwire_cm/D");

  for (const auto& r : rows_prim) {
    br_xwire0  = r.xwire;
    br_t0ns    = r.t0_ns;
    br_inhalo0 = r.in_halo;
    br_dwire0  = r.dwire_cm;
    t_prim.Fill();
  }

  // --- メタデータ ---
  TTree meta("meta","run metadata (random sampling region)");
  double meta_xmin=xminU, meta_xmax=xmaxU, meta_ymin=yminU, meta_ymax=ymaxU, meta_Lmax=1.0;
  int    meta_Nx=0, meta_Ny=0;
  meta.Branch("xmin",&meta_xmin,"xmin/D");
  meta.Branch("xmax",&meta_xmax,"xmax/D");
  meta.Branch("ymin",&meta_ymin,"ymin/D");
  meta.Branch("ymax",&meta_ymax,"ymax/D");
  meta.Branch("Nx",&meta_Nx,"Nx/I");
  meta.Branch("Ny",&meta_Ny,"Ny/I");
  meta.Branch("Lmax",&meta_Lmax,"Lmax/D");
  meta.Fill();

  t_all.Write();
  t_prim.Write();
  meta.Write();
  fout.Close();

  std::printf("[done] drift_all: %zu rows, drift_primary: %zu rows written\n",
              rows_all.size(), rows_prim.size());
  return 0;
}