#pragma once
#include <cmath>
#include <string>
#include <utility>
#include <stdexcept>
#include <cstdio>
#include "detector_geometry.hh"

namespace Detector {

inline std::pair<int,int> rational_approx(double r, int Nmax=24) {
  int best_m = 1, best_n = 1;
  double best_err = 1e300;
  for (int m = 1; m <= Nmax; ++m) {
    for (int n = 1; n <= Nmax; ++n) {
      const double err = std::abs(r - double(m)/double(n));
      if (err < best_err) { best_err = err; best_m = m; best_n = n; }
    }
  }
  return {best_m, best_n};
}

//============================================================
// AP（Anode 上 y=+gap、Potential 下 y=-gap）＋ワイヤカソード (Super-Cell 方式)
// Anode 側（+gap）と Potential 側（-gap）の非対称構成。
// * 周期境界(Periodic)をONにしつつ、全体を巨大な1つの周期(Super-Cell)とする。
// * 端の電場歪みを防ぎ、かつ全ワイヤーから独立した波形読み出しを可能にします。
//============================================================
inline Geometry PAP_WireCathode_Periodic( 
    double pitchSense = 0.40,   // [cm] A–A
    double gap        = 0.50,   // [cm]
    double pwOffset   = 0.0,    // [cm]
    double pitchCath  = 0.20,   // [cm] カソード希望ピッチ
    double cathPhase  = 0.0,    // [cm] カソード列の位相
    double rAn = 1.0e-3, double rPW = 2.5e-3, double rCat = 2.50e-3,
    double vAn = 0.0, double vPW = 0, double vCat = -700.0,
    bool pwSameAsAnode = false,
    bool printDebug    = true,
    int nWires         = 4     // 中心から左右に配置するセル数

) {
  (void)pwOffset; 
  const double L   = pitchSense;       
  const double off = L / 2.0;          

  if (off >= L) throw std::runtime_error("Invalid APA geometry: A–P must be < L.");
  // vPW = vCat * (rPW / gap - rAn / gap) / (1.0 - rAn / gap);

  Geometry g;
  g.periodicX = true; // ★ 魔法の復活：周期境界をONにする！
  
  // ★ スーパーセル全体の幅を計算し、それを新しい周期(pitchX)とする
  const int totalCells = 2 * nWires + 1;
  const double W = totalCells * L;
  g.pitchX = W; 

  // 描画用の表示範囲（スーパーセルの幅に合わせる）
  g.xmin = -W / 2.0; 
  g.xmax = +W / 2.0;
  const double yPad = std::max(5.0 * rCat, 1.0e-2); 
  g.ymin = -gap - yPad;
  g.ymax = +gap + yPad;

  //======================
  // Anode 側ワイヤ（y=+gap）と Potential 側ワイヤ（y=-gap）
  //======================
  int k = std::max(1, (int)std::round(L / pitchCath));
  const double pitchCathEff = L / k;  
  int numCathodes = totalCells * k;   // スーパーセルに入るカソードの総本数
  double startX = -W / 2.0;

  for (int j = 0; j < numCathodes; ++j) {
      // 境界での重複エラーを防ぐため、座標を必ず [-W/2, W/2) の間に収める
      double x0 = startX + j * pitchCathEff + cathPhase;
      while (x0 < -W / 2.0) x0 += W;
      while (x0 >= W / 2.0) x0 -= W;

      g.electrodes.push_back({ElectrodeKind::WireRow, x0, +gap, vCat, rCat, W, "AnodeCath_" + std::to_string(j)});  // Anode 側（上, +gap）
      g.electrodes.push_back({ElectrodeKind::WireRow, x0, -gap, vCat, rCat, W, "PotentialCath_" + std::to_string(j)});  // Potential 側（下, -gap）
  }

  //======================
  // 中央：Anode–PW ワイヤ列
  //======================
  const double vPW_eff = pwSameAsAnode ? vAn : vPW;
  for (int i = -nWires; i <= nWires; ++i) {
      double xa  = i * L;
      double xpr = i * L + off;
      
      std::string aName = "Anode_" + std::to_string(i);
      std::string pName = "PW_" + std::to_string(i);

      // ピッチ引数には、元のLではなくスーパーセル幅(W)を渡す
      g.electrodes.push_back({ElectrodeKind::WireRow, xa,  0.0, vAn,     rAn, W, aName});
      g.electrodes.push_back({ElectrodeKind::WireRow, xpr, 0.0, vPW_eff, rPW, W, pName});
  }

  if (printDebug) {
      std::fprintf(stderr,
        "[AP Super-Cell] L=%.3f, vPW=%.3fV, Wires: %d to %d (Total Width: %.3f cm)\n",
        L, vPW, -nWires, nWires, W);
  }
  return g;
}
} // namespace Detector