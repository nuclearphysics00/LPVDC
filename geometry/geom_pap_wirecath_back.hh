//geom_pap_wirecath.hh
#pragma once
#include <cmath>
#include <string>
#include <utility>
#include "detector_geometry.hh"

namespace Detector {

// r ≈ m/n を小さな整数で近似（1..Nmax）
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

// ★ 追加：セル幅と本数を整合的に決める（PWが重ならない L を返す）
struct PitchPlan {
  double L;         // [cm] セル幅
  double margin;    // [cm] 境界→PW 距離 m = L/2 - d
  int nSense;       // セル内アノード本数
  int nCath;        // セル内カソード本数
  int a, b, t;      // 係数（デバッグ用）
};

inline PitchPlan compute_pitch_plan(double pitchSense, double pitchCath, double pwOffset,
                                    int Nmax = 24, double tol = 1e-12) {
  // 小整数比 a/b ≈ pC/pS をとる → 共通周期 L0 ≈ a*pS ≈ b*pC
  auto [a, b] = rational_approx(pitchCath / pitchSense, Nmax);
  double L0s = a * pitchSense;
  double L0c = b * pitchCath;
  // 誤差が十分小さければ平均にスナップ（数値誤差安定化）
  double L0 = (std::abs(L0s - L0c) <= tol * std::max(L0s, L0c)) ? 0.5 * (L0s + L0c) : L0s;

  // PWが重ならない最小 t を選ぶ（L = t*L0 > 2d）
  int t = int(std::floor((2.0 * pwOffset) / L0)) + 1;
  if (t < 1) t = 1;

  double L = t * L0;
  double m = 0.5 * L - pwOffset;

  PitchPlan plan;
  plan.L = L;
  plan.margin = m;
  plan.nSense = int(std::round(L / pitchSense)); // = t*a（近似スナップ時でも整数に落ちる）
  plan.nCath  = int(std::round(L / pitchCath));  // = t*b
  plan.a = a; plan.b = b; plan.t = t;
  return plan;
}

// PAP（PW–Anode–PW）＋上下ワイヤカソード（Periodic）
inline Geometry PAP_WireCathode_Periodic(
    // 幾何
    double pitchSense = 0.20,   // [cm] センス（アノード列）ピッチ
    double gap        = 1.00,   // [cm] 中央面→上下カソード距離
    double pwOffset   = 0.20,   // [cm] アノードから PW までの横距離（=2 mm）
    double pitchCath  = 0.10,   // [cm] カソードワイヤの目標ピッチ（=1 mm なら 0.10）
    double cathPhase  = 0.0,    // [cm] カソード列の位相（X軸方向オフセット）
    // 半径
    double rAn = 6.0e-4, double rPW = 2.5e-3, double rCat = 5.0e-4,
    //double rAn = 1.0e-1, double rPW = 2.5e-3, double rCat = 5.0e-3,
    // 電位
    //double vAn = 0.0, double vPW = -0.00075, double vCat = -10.0,
    double vAn = 0.0, double vPW = -0.555, double vCat = -10.0,
    // オプション: PW をアノードと等電位にするか
    bool pwSameAsAnode = false
) {
  // ★ 変更点：L を「共通周期×必要倍数」で決め、PW重なりを回避
  const PitchPlan plan = compute_pitch_plan(pitchSense, pitchCath, pwOffset);
  const double L = plan.L;

  Geometry g;
  g.periodicX = true;
  g.pitchX    = L;
  g.xmin = -0.6; g.xmax = +0.6;
  g.ymin = -gap; g.ymax = +gap;

  // 位相を [-L/2, L/2) へ折りたたむ
  auto wrap = [&](double x){
    while (x <  -0.5*L) x += L;
    while (x >=  0.5*L) x -= L;
    return x;
  };

  // 半開区間 [-L/2, L/2) に入る整数 i の範囲を返すヘルパ
  auto index_range = [&](double pitch, double anchor){
    const double eps = 1e-12; // 右端重複防止用
    int i_min = (int)std::ceil(( -0.5*L - anchor) / pitch);
    int i_max = (int)std::floor(( +0.5*L - anchor) / pitch - eps);
    return std::pair<int,int>(i_min, i_max);
  };

  // --- 上下カソード列：アンカー = cathPhase（★ pitchCath は指定通り・独立）
  {
    const auto [i_min, i_max] = index_range(pitchCath, cathPhase);
    for (int i = i_min; i <= i_max; ++i) {
      const double x = wrap(cathPhase + i * pitchCath);
      g.electrodes.push_back({ElectrodeKind::WireRow, x, +gap, vCat, rCat, L, "CathTop"});
      g.electrodes.push_back({ElectrodeKind::WireRow, x, -gap, vCat, rCat, L, "CathBot"});
    }
  }

  // --- 中央：アノード列と PW（アンカー = 0.0 → 原点に必ずアノード）
  {
    const double vPW_eff = pwSameAsAnode ? vAn : vPW;
    const double anchorAnode = 0.0;
    const auto [i_min, i_max] = index_range(pitchSense, anchorAnode);
    for (int i = i_min; i <= i_max; ++i) {
      const double xa = wrap(anchorAnode + i * pitchSense);
      g.electrodes.push_back({ElectrodeKind::WireRow, xa,             0.0, vAn,     rAn, L, "Anode"});
      g.electrodes.push_back({ElectrodeKind::WireRow, xa - pwOffset,  0.0, vPW_eff, rPW, L, "PW_L"});
      g.electrodes.push_back({ElectrodeKind::WireRow, xa + pwOffset,  0.0, vPW_eff, rPW, L, "PW_R"});
    }
  }

  return g;
}

} // namespace Detector
