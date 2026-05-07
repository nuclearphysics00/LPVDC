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
// PAP（PW–Anode–PW）＋上下ワイヤカソード (Super-Cell 方式)
// * 周期境界(Periodic)をONにしつつ、全体を巨大な1つの周期(Super-Cell)とする。
//============================================================
inline Geometry PAP_WireCathode_Periodic(
    double pitchSense = 0.60,   // [cm] A–A
    double gap        = 0.50,   // [cm]
    double pwOffset   = 0.0,    // [cm]
    double pitchCath  = 0.20,   // [cm] カソード希望ピッチ
    double cathPhase  = 0.0,    // [cm] カソード列の位相
    double rAn = 1.0e-3, double rPW = 2.5e-3, double rCat = 2.50e-3,
    double vAn = 0.0, double vPW = -33.5, double vCat = -700.0,
    bool pwSameAsAnode = false,
    bool printDebug    = true,
    int nWires         = 4     // 中心から左右に配置するセル数
) {
  (void)pwOffset; 
  const double L   = pitchSense;       
  const double off = L / 3.0;          

  if (!(off > 0.0 && off < 0.5 * L))
    throw std::runtime_error("Invalid PAP geometry: A–P must be L/3.");

  //vPW = vCat * (rPW / gap - rAn / gap) / (1.0 - rAn / gap);

  Geometry g;
  g.periodicX = true; // 周期境界条件を有効化（Super-Cell 方式）
  
  // Super-Cell 全体幅を計算し、周期(pitchX)として設定する
  const int totalCells = 2 * nWires + 1;
  const double W = totalCells * L;
  g.pitchX = W; 

  g.xmin = -W / 2.0;
  g.xmax = +W / 2.0;
  const double yPad = std::max(5.0 * rCat, 1.0e-2); 
  g.ymin = -gap - yPad;
  g.ymax = +gap + yPad;

  //======================
  // 上下カソード（y=±gap）
  //======================
  int k = std::max(1, (int)std::round(L / pitchCath));
  const double pitchCathEff = L / k;  
  int numCathodes = totalCells * k;
  double startX = -W / 2.0;

  for (int j = 0; j < numCathodes; ++j) {
      double x0 = startX + j * pitchCathEff + cathPhase;
      while (x0 < -W / 2.0) x0 += W;
      while (x0 >= W / 2.0) x0 -= W;

      g.electrodes.push_back({ElectrodeKind::WireRow, x0, +gap, vCat, rCat, W, "CathTop_" + std::to_string(j)});
      g.electrodes.push_back({ElectrodeKind::WireRow, x0, -gap, vCat, rCat, W, "CathBot_" + std::to_string(j)});
  }

  //======================
  // 中央：P–A–P
  //======================
  const double vPW_eff = pwSameAsAnode ? vAn : vPW;
  for (int i = -nWires; i <= nWires; ++i) {
      double xa  = i * L;
      double xpl = i * L - off;
      double xpr = i * L + off;

      std::string aName   = "Anode_" + std::to_string(i);
      std::string pwLName = "PW_L_"  + std::to_string(i);
      std::string pwRName = "PW_R_"  + std::to_string(i);

      g.electrodes.push_back({ElectrodeKind::WireRow, xa,  0.0, vAn,     rAn, W, aName});
      g.electrodes.push_back({ElectrodeKind::WireRow, xpl, 0.0, vPW_eff, rPW, W, pwLName});
      g.electrodes.push_back({ElectrodeKind::WireRow, xpr, 0.0, vPW_eff, rPW, W, pwRName});
  }

  if (printDebug) {
      std::fprintf(stderr,
        "[PAP Super-Cell] L=%.3f, vPW=%.3fV, Wires: %d to %d (Total Width: %.3f cm)\n",
        L, vPW, -nWires, nWires, W);
  }
  return g;
}
} // namespace Detector