// eval_maps.hh
#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "TH2D.h"
#include "TH1D.h"

#include "Garfield/ComponentAnalyticField.hh"
#include "detector_geometry.hh"

// ============================================================
// ヒートマップ/ヒスト用の δ / δy を作るユーティリティ
//   - δy(x,y) = | |Ey| - E0(y) | / E0(y)
//   - δ (x,y) = | |E|  - E0(y) | / E0(y)
// ここで E0(y) はアノード–上下カソードの差電位を縦距離で割った
// 理想平行平板の強さ（y のみの関数）
// ============================================================

// 周期境界対応の「ワイヤ近傍判定」(半径×multR 以内なら true)
inline bool NearAnyWire(double x, double y, const Detector::Geometry& g, double multR) {
  if (multR <= 0.0) return false;
  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    double xw = e.x0;
    if (g.periodicX) {
      const double n = std::round((x - e.x0) / g.pitchX);
      xw = e.x0 + n * g.pitchX;
    }
    const double r = std::max(1e-4, e.radius);
    const double d = std::hypot(x - xw, y - e.y);
    if (d < multR * r) return true;
  }
  return false;
}

// アノード/上下カソードを拾って基準 E0 を作るための情報
struct ACBaseline {
  double yAn = 0.0, yTop = 0.0, yBot = 0.0;
  double VAn = 0.0, VTop = 0.0, VBot = 0.0;
  double E0_top = 0.0, E0_bot = 0.0; // [V/cm]
};

// 幾何から基準を構築（anode名優先、見つからなければ y=0/V=0 をフォールバック）
inline ACBaseline MakeACBaseline(const Detector::Geometry& g) {
  ACBaseline B{};
  bool foundAn = false;
  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    std::string nm = e.name;
    std::transform(nm.begin(), nm.end(), nm.begin(), ::tolower);
    if (nm.find("anode") != std::string::npos) {
      B.yAn = e.y; B.VAn = e.V; foundAn = true; break;
    }
  }
  if (!foundAn) { B.yAn = 0.0; B.VAn = 0.0; }

  double yTop = -1e99, yBot = +1e99, vTop = 0.0, vBot = 0.0;
  for (const auto& e : g.electrodes) {
    if (e.y > yTop) { yTop = e.y; vTop = e.V; }
    if (e.y < yBot) { yBot = e.y; vBot = e.V; }
  }
  B.yTop = yTop; B.yBot = yBot; B.VTop = vTop; B.VBot = vBot;

  const double dyTop = std::abs(B.yTop - B.yAn);
  const double dyBot = std::abs(B.yAn  - B.yBot);
  B.E0_top = (dyTop > 0.0) ? std::abs(B.VAn - B.VTop) / dyTop : 0.0;
  B.E0_bot = (dyBot > 0.0) ? std::abs(B.VAn - B.VBot) / dyBot : 0.0;
  return B;
}

// 生成物（2D ヒートマップ + 1D ヒスト）
struct DeltaProducts {
  TH2D* dY_all = nullptr;   // δy（全点）
  TH2D* d_all  = nullptr;   // δ （全点）
  TH2D* dY_noh = nullptr;   // δy（No-Halo：ワイヤ近傍除外）
  TH2D* d_noh  = nullptr;   // δ （No-Halo）

  TH1D* h_dY_all = nullptr; // δy ヒスト（全点）
  TH1D* h_d_all  = nullptr; // δ  ヒスト（全点）
  TH1D* h_dY_noh = nullptr; // δy ヒスト（No-Halo）
  TH1D* h_d_noh  = nullptr; // δ  ヒスト（No-Halo）
};

// δ/δy の 2D マップと 1D ヒストを作成
inline DeltaProducts
MakeDeltaMapsAndHists(Garfield::ComponentAnalyticField& comp,
                      const Detector::Geometry& g,
                      int nx = 220, int ny = 160,
                      double halo_mult = 5.0) {
  DeltaProducts P{};

  // 2D マップ
  P.dY_all = new TH2D("dY_all", "delta_{y} (ALL);X [cm];Y [cm]",
                      nx, g.xmin, g.xmax, ny, g.ymin, g.ymax);
  P.d_all  = new TH2D("d_all",  "delta (ALL);X [cm];Y [cm]",
                      nx, g.xmin, g.xmax, ny, g.ymin, g.ymax);
  P.dY_noh = new TH2D("dY_noH", "delta_{y} (No-Halo);X [cm];Y [cm]",
                      nx, g.xmin, g.xmax, ny, g.ymin, g.ymax);
  P.d_noh  = new TH2D("d_noH",  "delta (No-Halo);X [cm];Y [cm]",
                      nx, g.xmin, g.xmax, ny, g.ymin, g.ymax);
  for (auto* h : {P.dY_all, P.d_all, P.dY_noh, P.d_noh}) h->SetDirectory(nullptr);

  // 1D ヒスト
  P.h_dY_all = new TH1D("h_dY_all", "delta_{y} ALL;relative deviation;counts", 200, 0.0, 2.0);
  P.h_d_all  = new TH1D("h_d_all",  "delta ALL;relative deviation;counts",      200, 0.0, 2.0);
  P.h_dY_noh = new TH1D("h_dY_noH", "delta_{y} No-Halo;relative deviation;counts", 200, 0.0, 2.0);
  P.h_d_noh  = new TH1D("h_d_noH",  "delta No-Halo;relative deviation;counts",      200, 0.0, 2.0);
  for (auto* h : {P.h_dY_all, P.h_d_all, P.h_dY_noh, P.h_d_noh}) h->SetDirectory(nullptr);

  // 基準 E0(y) を用意
  const auto B = MakeACBaseline(g);
  auto E0_of_y = [&](double y) { return (y >= B.yAn) ? B.E0_top : B.E0_bot; };

  // スキャンして δ/δy を格納
  for (int ix = 1; ix <= nx; ++ix) {
    const double x = P.dY_all->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= ny; ++iy) {
      const double y = P.dY_all->GetYaxis()->GetBinCenter(iy);

      double ex = 0, ey = 0, ez = 0, V = 0;
      Garfield::Medium* med = nullptr;
      int st = 0;
      comp.ElectricField(x, y, 0.0, ex, ey, ez, V, med, st);

      const double Emag = std::hypot(ex, std::hypot(ey, ez));
      const double E0   = E0_of_y(y);
      if (E0 <= 0.0) continue;

      const double dY = std::abs(std::abs(ey) - E0) / E0;
      const double dA = std::abs(Emag - E0) / E0;

      // ALL
      P.dY_all->SetBinContent(ix, iy, dY);
      P.d_all ->SetBinContent(ix, iy, dA);
      P.h_dY_all->Fill(dY);
      P.h_d_all ->Fill(dA);

      // No-Halo（ワイヤ半径×halo_mult 以内を除外）
      if (!NearAnyWire(x, y, g, halo_mult)) {
        P.dY_noh->SetBinContent(ix, iy, dY);
        P.d_noh ->SetBinContent(ix, iy, dA);
        P.h_dY_noh->Fill(dY);
        P.h_d_noh ->Fill(dA);
      } else {
        // 可視化で“除外領域”が分かるよう 0 を入れておく（後で Z 範囲を設定すると白抜けになる）
        P.dY_noh->SetBinContent(ix, iy, 0.0);
        P.d_noh ->SetBinContent(ix, iy, 0.0);
      }
    }
  }

  return P;
}
