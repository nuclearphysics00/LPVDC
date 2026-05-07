#pragma once
#include "detector_geometry.hh"
#include <cmath>
#include <stdexcept>
#include <string>
#include <cstdio>

namespace Detector {

// 1セル = P–A–P。A–A = pitchSense、A–P = pitchSense/3 に固定 (Super-Cell 個別配置版)
inline Detector::Geometry PAP_PlaneCathode_Periodic(
    double pitchSense = 0.60,   // [cm]
    double gap        = 0.50,   // [cm]
    double phase      = 0.0,    // [cm]
    // 半径
    double rAn = 1.0e-3, double rPW = 2.5e-3,   // * 引数名を rPW に統一
    //double rAn = 1.0e-3, double rPW = 2.5e-3,   // * 川畑さんの論文値

    // 電位
     double vAn = 0.0, double vPW = -160, double vCat = -1100.0,
    //double vAn = 0.0, double vPW = -300, double vCat = -5000.0, // * 川畑さんの論文値

    bool pwSameAsAnode = false,
    int nWires = 4            // 中心から左右に配置するセル（アノード）の数
) {
  const double L   = pitchSense;      // 周期（=A–A）
  const double off = L / 3.0;         // A–P./
  if (off >= 0.5 * L) throw std::runtime_error("invalid PAP geometry");

   // vPW = vCat * (rPW / gap - rAn / gap) / (1.0 - rAn / gap);

  // Super-Cell 全体幅を計算する
  const int totalCells = 2 * nWires + 1;
  const double W = totalCells * L;

  std::fprintf(stderr,
      "[PAP Super-Cell] L(A–A)=%.3f cm, A–P=%.3f cm, vPW=%.3f V, vCat=%.1f V, Cells: %d to %d (Total Width: %.3f cm)\n",
      L, off, vPW, vCat, -nWires, nWires, W);

  Detector::Geometry g;
  g.periodicX = true; // 周期境界条件を有効化（Super-Cell 方式）
  g.pitchX    = W;    // Super-Cell 幅を周期として設定

  // 検出器の全体幅を設定 (スーパーセルの幅に合わせる)
  g.xmin = -W / 2.0;
  g.xmax = +W / 2.0;
  g.ymin = -gap; g.ymax = +gap;

  // カソード（上下）
  // 平面電極(PlaneY)はGarfield++内部でX方向に無限に広がるため、1つずつでOK
  // ピッチ引数には L ではなく Super-Cell 幅 W を渡す
  g.electrodes.push_back(Detector::Electrode{
      Detector::ElectrodeKind::PlaneY, 0.0, +gap, vCat, 0.0, W, "CathodeTop"});
  g.electrodes.push_back(Detector::Electrode{
      Detector::ElectrodeKind::PlaneY, 0.0, -gap, vCat, 0.0, W, "CathodeBottom"});

  const double vPW_eff = pwSameAsAnode ? vAn : vPW;
  
  // -nWires から +nWires まで各ワイヤを固有名で順に配置する
  for (int i = -nWires; i <= nWires; ++i) {
      // 1セル分（P-A-P）の各ワイヤーの座標を計算
      double xa   = phase + i * L;          // アノード
      double xpwL = phase + i * L - off;    // 左側のポテンシャルワイヤー
      double xpwR = phase + i * L + off;    // 右側のポテンシャルワイヤー
      
      // 一意なラベルを生成
      std::string aName   = "Anode_" + std::to_string(i);
      std::string pwLName = "PW_L_"  + std::to_string(i);
      std::string pwRName = "PW_R_"  + std::to_string(i);

      // アノード・PW配置 (ピッチ引数には L ではなく スーパーセル幅 W を渡す)
      g.electrodes.push_back(Detector::Electrode{
          Detector::ElectrodeKind::WireRow, xa, 0.0, vAn, rAn, W, aName});
      
      g.electrodes.push_back(Detector::Electrode{
          Detector::ElectrodeKind::WireRow, xpwL, 0.0, vPW_eff, rPW, W, pwLName});
          
      g.electrodes.push_back(Detector::Electrode{
          Detector::ElectrodeKind::WireRow, xpwR, 0.0, vPW_eff, rPW, W, pwRName});
  }

  return g;
}

} // namespace Detector