#pragma once
#include "detector_geometry.hh"
#include <cmath>
#include <stdexcept>
#include <string>
#include <cstdio>

namespace Detector {

inline Detector::Geometry PAP_PlaneCathode_Periodic(
    double pitchSense = 0.40,   
    double gap        = 0.50,   
    double phase      = 0.0,    
    double rAn = 1.0e-3, double rPW = 2.5e-3,
    double vAn = 0.0, double vPW = 0, double vCat = -1100.0,
    bool pwSameAsAnode = false,
    int nWires = 4             // 中心から左右に配置する本数
) {
  const double L   = pitchSense;      
  const double off = L / 2.0;         
  
  if (off >= L) throw std::runtime_error("invalid APA geometry");
  // vPW = vCat * (rPW / gap - rAn / gap) / (1.0 - rAn / gap); // スイープ時は外部から vPW を渡すためコメントアウト
  
  Detector::Geometry g;
  // 周期境界条件を有効化（Super-Cell 方式）
  g.periodicX = true; 
  
  // Super-Cell 全体幅（2*nWires+1 セル分）を周期として設定する
  int totalCells = 2 * nWires + 1;
  const double W = totalCells * L; // Super-Cell 全幅
  g.pitchX = W; // ワイヤ重複防止のため Super-Cell 幅 W を周期として設定

  // 描画用の表示範囲
  g.xmin = -(nWires + 1) * L; 
  g.xmax = +(nWires + 1) * L;
  g.ymin = -gap; g.ymax = +gap;

  // カソード平面
  g.electrodes.push_back(Detector::Electrode{
      Detector::ElectrodeKind::PlaneY, 0.0, +gap, vCat, 0.0, W, "CathodeTop"}); // L を W に変更
  g.electrodes.push_back(Detector::Electrode{
      Detector::ElectrodeKind::PlaneY, 0.0, -gap, vCat, 0.0, W, "CathodeBottom"}); // L を W に変更

  const double vPW_eff = pwSameAsAnode ? vAn : vPW;
  
  // スーパーセルの中に全ワイヤーを個別の名前で登録
  for (int i = -nWires; i <= nWires; ++i) {
      double xa  = phase + i * L;
      double xpr = phase + i * L + off;
      
      std::string aName = "Anode_" + std::to_string(i);
      std::string pName = "PW_" + std::to_string(i);

      g.electrodes.push_back(Detector::Electrode{
          Detector::ElectrodeKind::WireRow, xa, 0.0, vAn, rAn, W, aName}); // L を W に変更
      //g.electrodes.push_back(Detector::Electrode{
      //    Detector::ElectrodeKind::WireRow, xpr, 0.0, vPW_eff, rPW, W, pName}); // L を W に変更
  }

  return g;
}
} // namespace Detector