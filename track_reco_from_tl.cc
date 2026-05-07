// track_reco_avalanche.cc
//
// - t_to_Lo_table.csv から L0(t) テーブルを読み込み
// - PAP_WireCathode_Periodic ジオメトリ（detector_geometry.hh）を使用
// - 直線トラック上にクラスター（一次電子）を配置
// - 各クラスターについて AvalancheMC でアバランシェを発生させ、
//   ワイヤに到達した電子のうち最も早い到着時間 t_first を取得
// - t_first → L0(t_first) でドリフト長を再構成し、
//   再構成ヒットから最小二乗直線フィット
// - 検出器構造＋真のトラック＋再構成ヒット＋フィット直線を PNG に保存
//
// コンパイル例：
//  g++ track_reco_avalanche.cc \
//     `root-config --cflags --glibs` \
//     -I$GARFIELD_HOME/include \
//     -L$GARFIELD_HOME/lib64 -lGarfield \
//     -DGEOM_WIRE \
//     -o track_reco_avalanche
//
// 実行： ./track_reco_avalanche
//
// ※ カレントディレクトリに t_to_Lo_table.csv を置くこと。

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>
#include <fstream>
#include <sstream>
#include <memory>

#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TF1.h"
#include "TRandom3.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"

#include "detector_geometry.hh"

#if defined(GEOM_PLATE)
  #include "geometry/geom_pap_plate.hh"
#elif defined(GEOM_WIRE)
  #include "geometry/geom_pap_wirecath.hh"
#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

using namespace Garfield;

// ==================== Geometry helpers (fieldview_array.cc 由来) =============

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

// 最寄りワイヤ表面までの距離 d と、そのワイヤ半径 r, 位置 (xw,yw) を取得
struct WireNearest {
  double d_surface_cm;   // 表面までの距離（負なら「食い込み」）
  double r_cm;           // そのワイヤ半径
  double xw_cm, yw_cm;   // そのワイヤ中心
  bool   hasWire;        // 幾何にワイヤがあるか
};

static WireNearest NearestWireSurface(const Detector::Geometry& g,
                                      double x, double y) {
  WireNearest out;
  out.d_surface_cm = 1e300;
  out.r_cm         = 0.0;
  out.xw_cm        = 0.0;
  out.yw_cm        = 0.0;
  out.hasWire      = false;

  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    out.hasWire = true;
    const double k  = std::round((x - e.x0) / g.pitchX);
    const double xw = e.x0 + k * g.pitchX;
    const double d  = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) {
      out.d_surface_cm = d;
      out.r_cm         = e.radius;
      out.xw_cm        = xw;
      out.yw_cm        = e.y;
    }
  }
  return out;
}

// ==================== t_to_L0 テーブル読み込み ==============================

// CSV: 1行につき "t_ns, L0_cm"。
// ヘッダ行（文字を含む行）、"#" 始まりの行は無視。
static void LoadTtoLTable(const char* fname,
                          TGraph& gr_tL,   // t -> L0(t)
                          TGraph& gr_Lt) { // L -> t
  std::ifstream fin(fname);
  if (!fin) {
    std::fprintf(stderr, "[ERR] cannot open %s\n", fname);
    return;
  }

  std::vector<double> vt, vL;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty()) continue;
    if (line[0] == '#') continue;

    bool hasAlpha = false;
    for (char c : line) {
      if (std::isalpha((unsigned char)c)) { hasAlpha = true; break; }
    }
    if (hasAlpha) continue;

    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line);
    double t, L;
    if (!(ss >> t >> L)) continue;
    vt.push_back(t);
    vL.push_back(L);
  }

  const int n = (int)vt.size();
  if (n == 0) {
    std::fprintf(stderr, "[ERR] no numeric data in %s\n", fname);
    return;
  }

  gr_tL = TGraph(n, vt.data(), vL.data());
  gr_tL.SetName("gr_tL");
  gr_tL.SetTitle("L_{0}(t);t [ns];L_{0} [cm]");

  gr_Lt = TGraph(n, vL.data(), vt.data());
  gr_Lt.SetName("gr_Lt");
  gr_Lt.SetTitle("t(L_{0});L_{0} [cm];t [ns]");

  std::printf("[info] loaded %d points from %s\n", n, fname);
}

// ==================== main ===================================================

int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);

  // ROOT style
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameFillColor(0);
  gROOT->ForceStyle();

  // --- L0(t) のテーブル読み込み -------------------------------------------
  TGraph gr_tL, gr_Lt;
  LoadTtoLTable("/home/fumiya/garfield/garfieldpp/Examples/LAS_V7/root/run_20251112_033311/analysis_L0/t_to_L0_table.csv", gr_tL, gr_Lt);
  if (gr_tL.GetN() == 0) {
    std::fprintf(stderr, "[ERR] empty L0(t) graph. abort.\n");
    return 1;
  }

  // --- Gas -------------------------------------------------------------------
  MediumMagboltz gas;
  gas.LoadGasFile("ic4H10_100_0.1Torr.gas");
  // gas.Initialise(true); // 必要なら

  // --- Geometry --------------------------------------------------------------
  std::puts("[build] geometry = PAP_WireCathode_Periodic");
  Detector::Geometry geo =
      Detector::PAP_WireCathode_Periodic(/*pitch=*/0.20, /*gap=*/1.00);
  const double gap = DetectGap(geo);

  // --- Field component & Sensor ---------------------------------------------
  auto comp = Detector::BuildField(geo);
  comp->SetMedium(&gas);

  Sensor sensor;
  sensor.AddComponent(comp.get());
  sensor.SetArea(geo.xmin, geo.ymin, -0.1,
                 geo.xmax, geo.ymax, +0.1);

  // --- AvalancheMC -----------------------------------------------------------
  AvalancheMC amc;
  amc.SetSensor(&sensor);
  // 必要なら細かいステップ指定など
  // amc.SetDistanceSteps(1.e-3); // 10 µm

  // --- 真のトラック定義 ----------------------------------------------------
  // トラック:  y = a_true * x + b_true
  const double a_true = 0.5;   // 傾き
  const double b_true = -0.5;  // 切片

  const double x_min_trk = geo.xmin;
  const double x_max_trk = geo.xmax;

  // トラック上のクラスター数（種の数）
  const int nClusters = 40;

  TRandom3 rnd(0);

  std::vector<double> x_true, y_true;
  std::vector<double> x_reco, y_reco;

  x_true.reserve(nClusters);
  y_true.reserve(nClusters);
  x_reco.reserve(nClusters);
  y_reco.reserve(nClusters);

  // --- クラスターごとのループ ---------------------------------------------
  for (int i = 0; i < nClusters; ++i) {

    // 真のクラスター位置 (x_true, y_true)
    const double x = rnd.Uniform(x_min_trk, x_max_trk);
    const double y = a_true * x + b_true;

    x_true.push_back(x);
    y_true.push_back(y);

    // 最近傍ワイヤを探す
    auto nw = NearestWireSurface(geo, x, y);
    if (!nw.hasWire) {
      std::fprintf(stderr, "[WARN] no wire found for cluster %d\n", i);
      continue;
    }

    // 1) 一次電子1個を出してアバランシェ
    const double t_start = 0.0; // トラック通過時刻を 0 ns とする
    amc.AvalancheElectron(x, y, 0.0, t_start);

    const size_t nEnd = amc.GetNumberOfElectronEndpoints();
    if (nEnd == 0) {
      std::fprintf(stderr, "[WARN] cluster %d: no electron endpoints\n", i);
      continue;
    }

    // 2) ワイヤに到達した電子のうち最も早い到着時間を探す
    double t_first = 1e9;  // 大きな値で初期化
    for (size_t k = 0; k < nEnd; ++k) {
      double x0,y0,z0,t0, x1,y1,z1,t1;
      int status = 0;
      amc.GetElectronEndpoint(k, x0,y0,z0,t0, x1,y1,z1,t1, status);

      // 終点がこのワイヤ表面付近かどうかを判定
      const double r_end = std::hypot(x1 - nw.xw_cm, y1 - nw.yw_cm);
      const double tol   = 1.0e-3; // 10 µm
      const bool atWire  = std::fabs(r_end - nw.r_cm) < tol;

      if (!atWire) continue;
      if (t1 < t_first) t_first = t1;
    }

    if (t_first > 1e8) {
      std::fprintf(stderr,
                   "[INFO] cluster %d: no electron reached the nearest wire\n", i);
      continue;
    }

    // 3) t_first → L0(t_first) でドリフト長を再構成
    const double L_meas = gr_tL.Eval(t_first);  // [cm]

    // 4) ワイヤ中心 → クラスター方向の単位ベクトル
    double ux = x - nw.xw_cm;
    double uy = y - nw.yw_cm;
    double norm = std::hypot(ux, uy);
    if (norm == 0.0) { ux = 0.0; uy = 1.0; norm = 1.0; }
    ux /= norm;
    uy /= norm;

    // 中心からの距離は r_wire + L_meas
    const double R_reco = nw.r_cm + L_meas;
    const double xr = nw.xw_cm + R_reco * ux;
    const double yr = nw.yw_cm + R_reco * uy;

    x_reco.push_back(xr);
    y_reco.push_back(yr);
  }

  std::printf("[info] clusters generated: %zu  reco hits: %zu\n",
              x_true.size(), x_reco.size());

  // --- 図: 検出器構造 + 真のトラック + 再構成ヒット + フィット直線 ----------
  TCanvas c("c", "track reco with avalanche", 1000, 800);
  c.SetGrid();

  TH2D frame("frame",
             "Detector geometry and reconstructed track;X [cm];Y [cm]",
             10, geo.xmin, geo.xmax,
             10, -gap, +gap);
  frame.SetStats(0);
  frame.Draw();

  // 電極を描画
  DrawElectrodes(geo);

  // 真のトラック
  TF1 f_true("f_true", "[0]*x + [1]", geo.xmin, geo.xmax);
  f_true.SetParameter(0, a_true);
  f_true.SetParameter(1, b_true);
  f_true.SetLineStyle(2);
  f_true.SetLineWidth(2);
  f_true.Draw("same");

  // 再構成ヒット
  if (!x_reco.empty()) {
    TGraph gr_hits(x_reco.size(), x_reco.data(), y_reco.data());
    gr_hits.SetName("gr_hits");
    gr_hits.SetMarkerStyle(20);
    gr_hits.SetMarkerSize(1.2);
    gr_hits.Draw("P SAME");

    // フィット用にコピー（pol1）
    TGraph gr_fit = gr_hits;
    gr_fit.Fit("pol1", "Q"); // Q: quiet

    TF1* f_fit = gr_fit.GetFunction("pol1");
    if (f_fit) {
      f_fit->SetLineColor(kRed);
      f_fit->SetLineWidth(2);
      f_fit->Draw("same");
      std::printf("[fit] y = a*x + b:  a_fit = %.4f,  b_fit = %.4f\n",
                  f_fit->GetParameter(1), f_fit->GetParameter(0));
    }
  } else {
    std::printf("[WARN] no reconstructed hits to fit.\n");
  }

  c.SaveAs("track_reco_geometry.png");
  std::puts("[done] saved figure as track_reco_geometry.png");

  return 0;
}
