
#pragma once
#include "detector_geometry.hh"
#include <string>
#include <cmath>

namespace Detector {

// PAP（原点にAnode、±2mmにPW）。周期は Anode ピッチ=0.6 cm
inline Detector::Geometry PAP_PlaneCathode_Periodic(
    //double pitchSense = 0.60,    // [cm] = 6 mm（Anode→隣Anode）
    double pitchSense = 0.60,    // [cm] = 6 mm（Anode→隣Anode）
    double gap        = 1.00,    // [cm]
    double pwOffset   = 0.20,    // [cm] = 2 mm（原点からのPWの位置）
    // 半径
    double rAn = 6.0e-4, double rPW = 2.5e-3,
    // 電位
    //double vAn = 0.0, double vPW = -0.75, double vCat = -10.0,
    //double vAn = 0.0, double vPW = -0.00075, double vCat = -10.0,
    //double vAn = 0.0, double vPW = -0.555, double vCat = -10.0,
    double vAn = 0.0, double vPW = -10, double vCat = -10.0,
    // オプション：PWをアノードと等電位にする
    bool pwSameAsAnode = false
) {
  const double L = pitchSense; // 周期 = アノードピッチ

  Detector::Geometry g;
  g.periodicX = true;
  g.pitchX    = L;

  // 表示窓（必要に応じて調整）
  g.xmin = -0.6; g.xmax = +0.6;
  g.ymin = -gap; g.ymax = +gap;

  // 上下の平行平板
  g.electrodes.push_back({Detector::ElectrodeKind::PlaneY, 0.0, +gap, vCat, 0.0, L, "CathodeTop"});
  g.electrodes.push_back({Detector::ElectrodeKind::PlaneY, 0.0, -gap, vCat, 0.0, L, "CathodeBottom"});

  // 中央セル：原点にAnode、±2mmにPW（±4mmのPWと±6mmのAnodeは周期で生成される）
  const double vPW_eff = pwSameAsAnode ? vAn : vPW;
  g.electrodes.push_back({Detector::ElectrodeKind::WireRow,  0.0,        0.0, vAn,     rAn, L, "Anode"});
  g.electrodes.push_back({Detector::ElectrodeKind::WireRow, -pwOffset,   0.0, vPW_eff, rPW, L, "PW_L"});
  g.electrodes.push_back({Detector::ElectrodeKind::WireRow, +pwOffset,   0.0, vPW_eff, rPW, L, "PW_R"});

  return g;
}

} // namespace Detector
