// one_drift_wirecath.cc
#include <iostream>
#include <cmath>

#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"


#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"

#include "detector_geometry.hh"
#include "geometry/geom_pap_wirecath.hh"

using namespace Garfield;

// ======================================================
// WireNearest のローカル定義（fieldview_array.cc から移植）
// ======================================================
struct WireNearest {
  double d_surface_cm;   // 表面までの距離（負なら食い込み）
  double r_cm;           // ワイヤ半径
  double xw_cm, yw_cm;   // ワイヤ中心
  bool   hasWire;        // ワイヤがあるか
};

static WireNearest NearestWireSurface(const Detector::Geometry& g, double x, double y){
  WireNearest out;
  out.d_surface_cm = 1e300;
  out.r_cm = 0.0;
  out.xw_cm = 0.0;
  out.yw_cm = 0.0;
  out.hasWire = false;

  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    out.hasWire = true;

    // x 方向に周期配置
    const double k  = std::round((x - e.x0) / g.pitchX);
    const double xw = e.x0 + k * g.pitchX;

    const double d = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) {
      out.d_surface_cm = d;
      out.r_cm = e.radius;
      out.xw_cm = xw;
      out.yw_cm = e.y;
    }
  }
  return out;
}

// ======================================================
// main
// ======================================================
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);

  // Gas
  MediumMagboltz gas;
  gas.LoadGasFile("ic4H10_100_0.1Torr.gas");

  // Geometry (あなたの fieldview_array と同じ)
  Detector::Geometry geo =
      Detector::PAP_WireCathode_Periodic(/*pitch=*/0.20, /*gap=*/1.00);

  auto comp = Detector::BuildField(geo);
  comp->SetMedium(&gas);

  // Sensor
  Sensor sensor;
  sensor.AddComponent(comp.get());
  sensor.SetArea(geo.xmin, geo.ymin, -0.1,
                 geo.xmax, geo.ymax, +0.1);

  // Microscopic
  AvalancheMicroscopic mic;
  mic.SetSensor(&sensor);

  // --- テスト開始位置 ---
  double x0 = 0.00;
  double y0 = 0.50;
  double z0 = 0.00;
  double t0 = 0.00;

  // ガスチェック
  double ex,ey,ez,V; Medium* m=nullptr; int status=0;
  comp->ElectricField(x0,y0,z0, ex,ey,ez, V, m, status);

  if (m == nullptr || status != 0) {
    std::cout << "Start point not in gas region.\n";
    return 0;
  }

  std::cout << "=== Microscopic drift test ===\n";
  std::cout << "Start point: x=" << x0 << " y=" << y0 << "\n";
  std::cout << "|E| = " << std::sqrt(ex*ex+ey*ey+ez*ez) << " V/cm\n";

  // Microscopic drift
  mic.DriftElectron(
      x0, y0, z0,
      t0,
      0.0,      // initial energy
      0.0, 0.0, 0.0, // initial direction
      0         // unlimited collisions
  );

  size_t nEnd = mic.GetNumberOfElectronEndpoints();

  if (nEnd == 0) {
    std::cout << "No endpoints.\n";
    return 0;
  }

  double xa,ya,za,ta;
  double xb,yb,zb,tb;
  double e1,p1;
  int st;

  mic.GetElectronEndpoint(
      nEnd - 1,
      xa,ya,za,ta,
      xb,yb,zb,tb,
      e1,p1,
      st
  );

  std::cout << "\n=== Result ===\n";
  std::cout << "Final status = " << st << "\n";
  std::cout << "Start time   = " << ta << " ns\n";
  std::cout << "Final time   = " << tb << " ns\n";
  std::cout << "Drift time   = " << tb - ta << " ns\n";
  std::cout << "End point    = (" << xb << ", " << yb << ")\n";

  // ===============================
  //  追加：軌跡描画用 TGraph
  // ===============================
  TCanvas* c1 = new TCanvas("c1", "Electron drift trajectory", 600, 600);
  TGraph* gr = new TGraph();
  gr->SetTitle("Microscopic drift; x [cm]; y [cm]");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.8);
  gr->SetMarkerColor(kRed);

  // 全エンドポイントを読み出してグラフに追加
  size_t nSteps = mic.GetNumberOfElectronEndpoints();
  std::cout << "Number of steps = " << nSteps << "\n";

  for (size_t i = 0; i < nSteps; ++i) {
    double xa,ya,za,ta, xb,yb,zb,tb;
    double e1,p1;
    int st;

    mic.GetElectronEndpoint(i,
      xa,ya,za,ta,
      xb,yb,zb,tb,
      e1,p1,
      st
    );

    // 区間の終点 (xb, yb) をプロット
    gr->SetPoint(i, xb, yb);
  }

  gr->Draw("APL");

  c1->SaveAs("microscopic_track_xy.png");

  std::cout << "Saved track plot: microscopic_track_xy.png\n";


  return 0;
}
