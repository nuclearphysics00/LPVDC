/**
 * optimize_efield.cc
 *
 * [Plate / Wire 両対応 & ゲイン・飛跡評価 統合版 + L-t自動評価 + Ex境界評価 + Field score評価 + Summary出力]
 * - コンパイル時のマクロ (-DGEOM_PLATE / -DGEOM_WIRE) で処理を自動分岐。
 * - 6変数(Vc, Va, Vp, dc, da, dp)をスイープし、実行元のrunディレクトリ下に出力。
 * - 電場マップ画像へ設定パラメータを自動描画。
 * - 全コンフィギュレーションのemapを1つのディレクトリにまとめて出力。
 * - アバランシェシミュレーションによるゲイン評価、および飛跡再構成を実行。
 * - L-t (drift length vs drift time) の自動評価を追加し、ROOTへ保存。
 * - L-t residual RMS と R^2 を電場マップ画像へ追記。
 * - 電気力線描画に使った same seed から、増幅なしドリフトで
 *   吸い込まれたワイヤーIDごとの drift time / drift path を ROOT に保存。
 * - 各設定ごとに x=2 mm 断面で Ex vs Y 図を保存。
 * - 中心 A0 / A1 / PW_-1 / PW_0 から Field score を計算して PNG / CSV / ROOT へ保存。
 *   * S_A_trans   : A0 と A1 の並進対称性
 *   * S_A_align   : A0 周りの縦方向整列
 *   * S_PW_mirror : 左右 PW の鏡映対称性
 *   * S_PW_align  : 左右 PW 周りの縦方向整列
 *   * S_field_raw : Atrans + Aalign + PWmirror + PWalign (geometry-only raw diagnostic)
 *   * S_field_norm: median-normalized sum of Atrans, Aalign, PWmirror, PWalign, and <Ex^2>
 * - Summary として
 *   - ExProfile_x2mm_Summary
 *   - AllWire_Lt_Summary
 *   - FieldScore_Summary
 *   を sweepBaseDir 直下に作成し、全設定の図を保存。
 * - RMS/R^2/Ex^2/Field score のランキング md を出力。
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <fstream>
#include <sstream>
#include <numeric>
#include <filesystem>
#include <type_traits>
#include <utility>
#include <memory>

// ROOT includes
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TLine.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TGraph.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TPaveText.h"

// Garfield++ includes
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh"

// Local includes
#include "detector_geometry.hh"

#if defined(GEOM_PLATE)
  #include "geometry/geom_pap_plate.hh"
  //#include "geometry/geom_ap_plate.hh"
#elif defined(GEOM_WIRE)
  #include "geometry/geom_pap_wirecath.hh"
  //#include "geometry/geom_ap_wirecath.hh"
#else
  #error "Please define either GEOM_PLATE or GEOM_WIRE during compilation."
#endif

#define DEBUG_WAVEFORM
// #define DEBUG_LT_IO
// #define DEBUG_FIELDLINE_IO

using namespace Garfield;
namespace fs = std::filesystem;

// ==================== Preamplifier Response ====================
double TransferFunction(double t) {
  constexpr double tau = 20.0;
  constexpr double target_gain_mV_per_fC = 1.0;
  constexpr double gain = target_gain_mV_per_fC;
  if (t < 0.0) return 0.0;
  return gain * (t / tau) * std::exp(1.0 - t / tau);
}

// ===== Helpers ==========================
template <typename T, typename = void>
struct has_EnableSignalCalculation : std::false_type {};

template <typename T>
struct has_EnableSignalCalculation<T, std::void_t<decltype(std::declval<T&>().EnableSignalCalculation())>>
  : std::true_type {};

static void TryEnableSignalCalculation(Garfield::AvalancheMC& av) {
  if constexpr (has_EnableSignalCalculation<Garfield::AvalancheMC>::value) {
    av.EnableSignalCalculation();
  }
}

struct WireNearest {
  double d_surface_cm = 1.e300;
  double r_cm = 0.0;
  double xw_cm = 0.0, yw_cm = 0.0;
  bool hasWire = false;
  int wire_type = -1;  // 0: anode, 1: pw
  std::string wire_name;
  std::string wire_uid;
};

struct GainPlotInfo {
  int wType = -1;
  double xw = 0.0, yw = 0.0, gain = 0.0;
  int n_primary = 0, n_total = 0;
  std::string wire_name;
  std::string wire_uid;
};

struct TimeHist {
  std::vector<double> t_center_ns;
  std::vector<double> dt_ns;
  std::vector<double> count;
};

struct LtPoint {
  double t_ns = 0.0;
  double L_cm = 0.0;
  double x0_cm = 0.0;
  double y0_cm = 0.0;
  double x1_cm = 0.0;
  double y1_cm = 0.0;
  double path_cm = 0.0;
  double dx_end_cm = 0.0;
  double wire_x_cm = 0.0;
  double wire_y_cm = 0.0;
  std::string wire_name = "";
  std::string wire_uid = "";
  std::string target_wire_name = "";
  std::string target_wire_uid = "";
  int wire_type = 0;
  int is_target = 0;
  int status = 0;
};

struct LtFitSummary {
  int n = 0;
  double slope = 0.0;
  double intercept = 0.0;
  double r2 = 0.0;
  double residual_rms = 0.0;
  int monotonic_violations = 0;
  double mean_dx_end = 0.0;
  double mean_path_excess = 0.0;
  double collect_eff = 0.0;
};

struct ExBoundarySummary {
  double x_eval_cm = 0.20;
  double ex2_mean = 0.0;
  double ex2_int = 0.0;
  double ex_abs_max = 0.0;
  double ex_mean_abs = 0.0;
};

struct ExProfileForOverlay {
  std::string label;
  double vPW = 0.0;
  std::vector<double> y;
  std::vector<double> ex;
};

struct SignalWireRef {
  std::string name;
  std::string uid;
  double x = 0.0;
  double y = 0.0;
  double r = 0.0;
  int type = -1;  // 0: anode, 1: pw
};

struct AdjacentPwPair {
  SignalWireRef anode;
  SignalWireRef pwLeft;
  SignalWireRef pwRight;
  bool valid = false;
};

struct FieldScoreSummary {
  bool valid = false;
  double pitchSense = 0.0;
  double gap_cm = 0.0;
  double yMin_cm = 0.0;
  double yMax_cm = 0.0;
  double uMaxA_cm = 0.0;
  double uMaxPW_cm = 0.0;
  int nAtrans = 0;
  int nAalign = 0;
  int nPWmirror = 0;
  int nPWalign = 0;
  int nELine = 0;
  double sAtrans = 0.0;
  double sAalign = 0.0;
  double sPWmirror = 0.0;
  double sPWalign = 0.0;
  double sELine = 0.0;
  double sFieldRaw = 0.0;
  std::string a0Name, a1Name, pwLName, pwRName;
};

struct ConfigRankingEntry {
  std::string conf_name;
  std::string geom_tag;
  double vCat = 0.0;
  double vAnode = 0.0;
  double vPW = 0.0;
  double dCat_um = 0.0;
  double dAnode_um = 0.0;
  double dPW_um = 0.0;
  int n_points = 0;
  double residual_rms = 0.0;
  double r2 = 0.0;
  int monotonic_violations = 0;
  double mean_dx_end = 0.0;
  double mean_path_excess = 0.0;
  double collect_eff = 0.0;
  double ex_boundary_x_cm = 0.20;
  double ex2_mean = 0.0;
  double ex2_int = 0.0;
  double ex_abs_max = 0.0;
  double ex_mean_abs = 0.0;
  double sAtrans = -1.0;
  double sAalign = -1.0;
  double sPWmirror = -1.0;
  double sPWalign = -1.0;
  double sELine = -1.0;
  double sFieldRaw = -1.0;
  double sAtransNorm = -1.0;
  double sAalignNorm = -1.0;
  double sPWmirrorNorm = -1.0;
  double sPWalignNorm = -1.0;
  double sELineNorm = -1.0;
  double sEx2Norm = -1.0;
  double sFieldNorm = -1.0;
  int nAtrans = 0;
  int nAalign = 0;
  int nPWmirror = 0;
  int nPWalign = 0;
  int nELine = 0;
  std::string combo_png;
  std::string emap_png;
  std::string ex_profile_png;
  std::string allwire_lt_png;
  std::string field_score_png;
  std::string field_score_root;
  std::string field_score_csv;
};

static std::string SanitizeFileName(std::string s) {
  for (char& c : s) {
    const bool ok = std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-';
    if (!ok) c = '_';
  }
  return s;
}

static std::string ToLower(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

static double NormalizeSignedZero(double x) {
  return (std::abs(x) < 1.e-12) ? 0.0 : x;
}

// 座標ベースUID生成
static std::string BuildWireUid(const std::string& wire_name, double x_cm) {
  constexpr double quantum_cm = 1.e-6;
  const long long ix = std::llround(NormalizeSignedZero(x_cm) / quantum_cm);
  return wire_name + "@xq_" + std::to_string(ix);
}

static bool IsSignalWireName(const std::string& name) {
  const std::string nm = ToLower(name);
  return (nm.find("anode") != std::string::npos ||
          nm.find("pw")    != std::string::npos);
}

static bool IsCentralSignalWire(const Detector::Electrode& e) {
  if (e.kind != Detector::ElectrodeKind::WireRow) return false;
  if (!IsSignalWireName(e.name)) return false;
  return std::abs(e.y) < 1.e-9;
}

static void DrawElectrodes(const Detector::Geometry& g) {
  const double xmin = g.xmin, xmax = g.xmax;

  for (const auto& e : g.electrodes) {
    if (e.kind == Detector::ElectrodeKind::PlaneY) {
      auto ln = new TLine(xmin, e.y, xmax, e.y);
      ln->SetBit(kCanDelete);
      ln->SetLineColor(kBlack);
      ln->SetLineStyle(2);
      ln->SetLineWidth(2);
      ln->Draw("same");
      continue;
    }

    Color_t color = kBlack;
    auto draw_wire = [&](double x) {
      if (x < xmin || x > xmax) return;
      auto el = new TEllipse(x, e.y, e.radius, e.radius);
      el->SetBit(kCanDelete);
      el->SetLineColor(color);
      el->SetFillColor(color);
      el->SetFillStyle(1001);
      el->SetLineWidth(1);
      el->Draw("same");
    };

    if (!g.periodicX) {
      draw_wire(e.x0);
    } else {
      const double pitch = g.pitchX;
      int nL = int(std::floor((xmin - e.x0) / pitch)) - 1;
      int nR = int(std::ceil ((xmax - e.x0) / pitch)) + 1;
      for (int n = nL; n <= nR; ++n) draw_wire(e.x0 + n * pitch);
    }
  }
}

static double DetectGap(const Detector::Geometry& g) {
  double gap = 0.0;
  for (const auto& e : g.electrodes) gap = std::max(gap, std::abs(e.y));
  if (gap <= 0.0) gap = std::min(std::abs(g.ymax), std::abs(g.ymin));
  return gap;
}

// A-A ピッチを geometry から動的に推定
static double DetectSensePitchX(const Detector::Geometry& g) {
  std::vector<double> xs;
  xs.reserve(g.electrodes.size());

  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    const std::string nm = ToLower(e.name);
    if (nm.find("anode") == std::string::npos) continue;
    if (std::abs(e.y) > 1.e-9) continue;
    xs.push_back(e.x0);
  }

  if (xs.size() < 2) return 0.0;

  std::sort(xs.begin(), xs.end());

  double pitch = std::numeric_limits<double>::max();
  for (size_t i = 1; i < xs.size(); ++i) {
    const double dx = xs[i] - xs[i - 1];
    if (dx > 1.e-12) pitch = std::min(pitch, dx);
  }

  return (pitch < std::numeric_limits<double>::max()) ? pitch : 0.0;
}

// seed開始領域を動的に決定
static std::pair<double, double> ComputeSeedXRange(const Detector::Geometry& g,
                                                   int marginCells = 1) {
  const double L = DetectSensePitchX(g);
  if (L <= 0.0) return {g.xmin, g.xmax};

  double xmin = g.xmin + marginCells * L;
  double xmax = g.xmax - marginCells * L;

  if (xmin >= xmax) return {g.xmin, g.xmax};
  return {xmin, xmax};
}

// アノード中心で各セルを均等分割した x-seed を生成
static std::vector<double> BuildAnodeCenteredSeedX(const Detector::Geometry& g,
                                                   int marginCells,
                                                   int nPerCell) {
  std::vector<double> xs;
  const double L = DetectSensePitchX(g);
  if (L <= 0.0 || nPerCell <= 0) return xs;

  const double xmin = g.xmin + marginCells * L;
  const double xmax = g.xmax - marginCells * L;

  // アノード中心のセル番号範囲
  const int mMin = (int)std::ceil(xmin / L);
  const int mMax = (int)std::floor(xmax / L);

  for (int m = mMin; m <= mMax; ++m) {
    const double xa = m * L;  // anode center
    for (int j = 0; j < nPerCell; ++j) {
      const double u = -0.5 + (double(j) + 0.5) / double(nPerCell);
      const double x = xa + u * L;
      if (x >= xmin && x <= xmax) xs.push_back(x);
    }
  }

  std::sort(xs.begin(), xs.end());
  return xs;
}

static WireNearest NearestWireSurface(const Detector::Geometry& g, double x, double y) {
  WireNearest out;

  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;

    const std::string nm = ToLower(e.name);
    const bool isAnode = nm.find("anode") != std::string::npos;
    const bool isPW    = nm.find("pw")    != std::string::npos;
    if (!isAnode && !isPW) continue;

    double xw = e.x0;
    if (g.periodicX) {
      const double k = std::round((x - e.x0) / g.pitchX);
      xw = e.x0 + k * g.pitchX;
    }
    xw = NormalizeSignedZero(xw);

    const double d = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) {
      out.d_surface_cm = d;
      out.r_cm = e.radius;
      out.xw_cm = xw;
      out.yw_cm = e.y;
      out.hasWire = true;
      out.wire_type = isAnode ? 0 : 1;
      out.wire_name = e.name;
      out.wire_uid  = BuildWireUid(e.name, xw);
    }
  }
  return out;
}

static WireNearest FindTargetSignalWireByX(const Detector::Geometry& g, double x) {
  WireNearest out;
  for (const auto& e : g.electrodes) {
    if (!IsCentralSignalWire(e)) continue;

    double xw = e.x0;
    if (g.periodicX) {
      const double k = std::round((x - e.x0) / g.pitchX);
      xw = e.x0 + k * g.pitchX;
    }
    xw = NormalizeSignedZero(xw);

    const double d = std::abs(x - xw);
    if (d < out.d_surface_cm) {
      out.d_surface_cm = d;
      out.r_cm = e.radius;
      out.xw_cm = xw;
      out.yw_cm = e.y;
      out.hasWire = true;
      out.wire_name = e.name;
      out.wire_uid  = BuildWireUid(e.name, xw);

      const std::string nm = ToLower(e.name);
      out.wire_type = (nm.find("anode") != std::string::npos) ? 0 : 1;
    }
  }
  return out;
}

static bool LoadTimeHistogram(const char* fname, TimeHist& h) {
  std::ifstream fin(fname);
  if (!fin) return false;

  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;

    bool hasAlpha = false;
    for (char c : line) {
      if (std::isalpha((unsigned char)c)) {
        hasAlpha = true;
        break;
      }
    }
    if (hasAlpha) continue;

    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line);
    double tcen, dt, cnt;
    if (!(ss >> tcen >> dt >> cnt)) continue;

    h.t_center_ns.push_back(tcen);
    h.dt_ns.push_back(dt);
    h.count.push_back(cnt);
  }
  return !h.t_center_ns.empty();
}

static double CalculateL0Directly(const TimeHist& h, double gap_cm, double t_hit) {
  if (h.t_center_ns.empty()) return 0.0;

  double numerator = 0.0, denominator = 0.0;
  for (size_t i = 0; i < h.t_center_ns.size(); ++i) {
    denominator += h.count[i];
    if (h.t_center_ns[i] <= t_hit) numerator += h.count[i];
  }
  return denominator > 0.0 ? (numerator / denominator) * gap_cm : 0.0;
}

static fs::path GuessRunDirFromCsv(const fs::path& csvPath) {
  fs::path p = csvPath.parent_path();
  std::string leaf = p.filename().string();

  if (leaf == "analysis_L0_prim" || leaf == "analysis_L0" || leaf == "analysis" ||
      leaf == "analysis_L0_prim_exact") {
    p = p.parent_path();
  }

  fs::path q = p;
  while (!q.empty()) {
    if (q.filename().string().rfind("run_", 0) == 0) return q;
    if (!q.has_parent_path()) break;
    q = q.parent_path();
  }
  return p;
}

static LtFitSummary AnalyzeLt(const std::vector<LtPoint>& pts) {
  LtFitSummary out;
  out.n = (int)pts.size();
  if (pts.size() < 2) return out;

  double St = 0.0, SL = 0.0, Stt = 0.0, StL = 0.0;
  for (const auto& p : pts) {
    St  += p.t_ns;
    SL  += p.L_cm;
    Stt += p.t_ns * p.t_ns;
    StL += p.t_ns * p.L_cm;
    out.mean_dx_end += std::abs(p.dx_end_cm);
    out.mean_path_excess += std::max(0.0, p.path_cm - p.L_cm);
    out.collect_eff += p.is_target ? 1.0 : 0.0;
  }
  out.mean_dx_end /= pts.size();
  out.mean_path_excess /= pts.size();
  out.collect_eff /= pts.size();

  const double n = (double)pts.size();
  const double denom = n * Stt - St * St;
  if (std::abs(denom) < 1.e-20) return out;

  out.slope = (n * StL - St * SL) / denom;
  out.intercept = (Stt * SL - St * StL) / denom;

  double ss_res = 0.0;
  const double Lmean = SL / n;
  double ss_tot = 0.0;
  for (const auto& p : pts) {
    const double pred = out.slope * p.t_ns + out.intercept;
    const double res = p.L_cm - pred;
    ss_res += res * res;
    const double d = p.L_cm - Lmean;
    ss_tot += d * d;
  }
  out.residual_rms = std::sqrt(ss_res / n);
  out.r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;

  std::vector<LtPoint> sorted = pts;
  std::sort(sorted.begin(), sorted.end(),
            [](const LtPoint& a, const LtPoint& b) { return a.t_ns < b.t_ns; });

  for (size_t i = 1; i < sorted.size(); ++i) {
    if (sorted[i].L_cm + 1.e-12 < sorted[i - 1].L_cm) {
      out.monotonic_violations++;
    }
  }

  return out;
}

static ExBoundarySummary EvaluateBoundaryEx(
    Garfield::Component* comp,
    const double x_eval_cm,
    TH1D* hExY) {

  ExBoundarySummary out;
  out.x_eval_cm = x_eval_cm;
  if (!hExY || !comp) return out;

  const int n = hExY->GetNbinsX();
  double sum_ex2 = 0.0;
  double sum_abs = 0.0;
  double max_abs = 0.0;
  double dy_ref = 0.0;
  int nvalid = 0;

  for (int i = 1; i <= n; ++i) {
    const double y = hExY->GetXaxis()->GetBinCenter(i);
    const double dy = hExY->GetXaxis()->GetBinWidth(i);
    if (i == 1) dy_ref = dy;

    double ex = 0.0, ey = 0.0, ez = 0.0, v = 0.0;
    Garfield::Medium* med = nullptr;
    int status = 0;
    comp->ElectricField(x_eval_cm, y, 0.0, ex, ey, ez, v, med, status);

    if (status != 0 || med == nullptr) {
      hExY->SetBinContent(i, 0.0);
      continue;
    }

    hExY->SetBinContent(i, ex);
    sum_ex2 += ex * ex;
    sum_abs += std::abs(ex);
    max_abs = std::max(max_abs, std::abs(ex));
    ++nvalid;
  }

  if (nvalid > 0) {
    out.ex2_mean = sum_ex2 / nvalid;
    out.ex_mean_abs = sum_abs / nvalid;
    out.ex_abs_max = max_abs;
    out.ex2_int = sum_ex2 * dy_ref;
  }
  return out;
}

static std::vector<SignalWireRef> CollectCentralSignalWires(const Detector::Geometry& g) {
  std::vector<SignalWireRef> out;
  out.reserve(g.electrodes.size());

  for (const auto& e : g.electrodes) {
    if (!IsCentralSignalWire(e)) continue;

    const std::string nm = ToLower(e.name);
    const int type = (nm.find("anode") != std::string::npos) ? 0 : 1;

    SignalWireRef w;
    w.name = e.name;
    w.x = NormalizeSignedZero(e.x0);
    w.y = e.y;
    w.r = e.radius;
    w.type = type;
    w.uid = BuildWireUid(e.name, w.x);
    out.push_back(w);
  }

  std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) {
    return a.x < b.x;
  });
  return out;
}

static std::vector<AdjacentPwPair> BuildAdjacentPwPairs(const Detector::Geometry& g) {
  const auto wires = CollectCentralSignalWires(g);

  std::vector<SignalWireRef> anodes;
  std::vector<SignalWireRef> pws;
  for (const auto& w : wires) {
    if (w.type == 0) anodes.push_back(w);
    else if (w.type == 1) pws.push_back(w);
  }

  std::vector<AdjacentPwPair> pairs;
  pairs.reserve(anodes.size());

  for (const auto& a : anodes) {
    AdjacentPwPair p;
    p.anode = a;

    double bestLeftDx = 1.e300;
    double bestRightDx = 1.e300;
    bool hasLeft = false;
    bool hasRight = false;

    for (const auto& pw : pws) {
      const double dx = pw.x - a.x;
      if (dx < 0.0) {
        const double adx = std::abs(dx);
        if (adx < bestLeftDx) {
          bestLeftDx = adx;
          p.pwLeft = pw;
          hasLeft = true;
        }
      } else if (dx > 0.0) {
        const double adx = std::abs(dx);
        if (adx < bestRightDx) {
          bestRightDx = adx;
          p.pwRight = pw;
          hasRight = true;
        }
      }
    }

    p.valid = hasLeft && hasRight;
    if (p.valid) pairs.push_back(p);
  }

  return pairs;
}

static bool GetFieldAt(Garfield::Component* comp,
                       const double x, const double y, const double z,
                       double& ex, double& ey, double& ez, double& v,
                       int& status) {
  Garfield::Medium* med = nullptr;
  comp->ElectricField(x, y, z, ex, ey, ez, v, med, status);
  return (status == 0 && med != nullptr);
}

static double SafeMean(const double sum, const int n) {
  return (n > 0) ? (sum / double(n)) : 0.0;
}

static void DrawInfoPave(double x1, double y1, double x2, double y2,
                         std::initializer_list<std::string> lines);

static FieldScoreSummary EvaluateAndSaveFieldScores(
    Garfield::Component* comp,
    const Detector::Geometry& geo,
    const double sensePitchX,
    const fs::path& outDir,
    const TString& textVolt,
    const TString& textDiam,
    const std::string& confName) {

  FieldScoreSummary out;
  if (!comp || sensePitchX <= 0.0) return out;

  const auto pairs = BuildAdjacentPwPairs(geo);
  if (pairs.empty()) {
    std::cerr << "[warn] No valid anode-with-adjacent-PW pairs found for field-score evaluation.\n";
    return out;
  }

  fs::create_directories(outDir);

  constexpr double eps = 1.e-20;
  const double gap = DetectGap(geo);
  const double yMin = -gap;
  const double yMax = +gap;
  const int nY = 181;

  // Store per-anode line-comparison curves.
  std::vector<std::unique_ptr<TH1D>> hEAnodes;
  std::vector<std::unique_ptr<TH1D>> hEpwLefts;
  std::vector<std::unique_ptr<TH1D>> hEpwRights;
  std::vector<std::unique_ptr<TH1D>> hDiff2s;
  hEAnodes.reserve(pairs.size());
  hEpwLefts.reserve(pairs.size());
  hEpwRights.reserve(pairs.size());
  hDiff2s.reserve(pairs.size());

  std::ofstream lineCsv(outDir / "field_line_compare.csv");
  if (lineCsv) {
    lineCsv << "anode_name,pw_left_name,pw_right_name,y_cm,Eabs_A_Vpcm,Eabs_PWL_Vpcm,Eabs_PWR_Vpcm,rel_diff2\n";
  }

  double sumELine = 0.0;
  int nELine = 0;

  for (size_t ip = 0; ip < pairs.size(); ++ip) {
    const auto& pair = pairs[ip];
    const double rMask = std::max({pair.anode.r, pair.pwLeft.r, pair.pwRight.r});

    auto hEA = std::make_unique<TH1D>(Form("hElineA_%zu", ip),
        Form("|E|(x=%s,y);y [cm];|E| [V/cm]", pair.anode.name.c_str()),
        nY, yMin, yMax);
    auto hEL = std::make_unique<TH1D>(Form("hElinePWL_%zu", ip),
        Form("|E|(x=%s,y);y [cm];|E| [V/cm]", pair.pwLeft.name.c_str()),
        nY, yMin, yMax);
    auto hER = std::make_unique<TH1D>(Form("hElinePWR_%zu", ip),
        Form("|E|(x=%s,y);y [cm];|E| [V/cm]", pair.pwRight.name.c_str()),
        nY, yMin, yMax);
    auto hD2 = std::make_unique<TH1D>(Form("hElineDiff2_%zu", ip),
        Form("Relative field-balance score around %s;y [cm];0.5[(#Delta_{L})^{2}+(#Delta_{R})^{2}]", pair.anode.name.c_str()),
        nY, yMin, yMax);

    for (int iy = 1; iy <= nY; ++iy) {
      const double y = hEA->GetXaxis()->GetBinCenter(iy);
      if (std::abs(y - pair.anode.y) <= rMask) continue;

      double exA = 0.0, eyA = 0.0, ezA = 0.0, vA = 0.0;
      double exL = 0.0, eyL = 0.0, ezL = 0.0, vL = 0.0;
      double exR = 0.0, eyR = 0.0, ezR = 0.0, vR = 0.0;
      int stA = 0, stL = 0, stR = 0;

      const bool okA = GetFieldAt(comp, pair.anode.x,  y, 0.0, exA, eyA, ezA, vA, stA);
      const bool okL = GetFieldAt(comp, pair.pwLeft.x, y, 0.0, exL, eyL, ezL, vL, stL);
      const bool okR = GetFieldAt(comp, pair.pwRight.x,y, 0.0, exR, eyR, ezR, vR, stR);

      double eA = -1.0, eL = -1.0, eR = -1.0;
      if (okA) {
        eA = std::sqrt(exA * exA + eyA * eyA + ezA * ezA);
        hEA->SetBinContent(iy, eA);
      }
      if (okL) {
        eL = std::sqrt(exL * exL + eyL * eyL + ezL * ezL);
        hEL->SetBinContent(iy, eL);
      }
      if (okR) {
        eR = std::sqrt(exR * exR + eyR * eyR + ezR * ezR);
        hER->SetBinContent(iy, eR);
      }

      double diff2 = -1.0;
      if (eA >= 0.0 && eL >= 0.0 && eR >= 0.0) {
        const double dL = (eA - eL) / (eA + eL + eps);
        const double dR = (eA - eR) / (eA + eR + eps);
        diff2 = 0.5 * (dL * dL + dR * dR);
        hD2->SetBinContent(iy, diff2);
        sumELine += diff2;
        ++nELine;
      }

      if (lineCsv) {
        lineCsv << pair.anode.name << ","
                << pair.pwLeft.name << ","
                << pair.pwRight.name << ","
                << y << ","
                << eA << ","
                << eL << ","
                << eR << ","
                << diff2 << "\n";
      }
    }

    hEAnodes.push_back(std::move(hEA));
    hEpwLefts.push_back(std::move(hEL));
    hEpwRights.push_back(std::move(hER));
    hDiff2s.push_back(std::move(hD2));
  }

  out.valid = true;
  out.pitchSense = sensePitchX;
  out.gap_cm = gap;
  out.yMin_cm = yMin;
  out.yMax_cm = yMax;
  out.uMaxA_cm = 0.0;
  out.uMaxPW_cm = 0.0;
  out.nAtrans = 0;
  out.nAalign = 0;
  out.nPWmirror = 0;
  out.nPWalign = 0;
  out.nELine = nELine;
  out.sAtrans = -1.0;
  out.sAalign = -1.0;
  out.sPWmirror = -1.0;
  out.sPWalign = -1.0;
  out.sELine = SafeMean(sumELine, nELine);
  out.sFieldRaw = out.sELine;
  out.a0Name = pairs.front().anode.name;
  out.a1Name.clear();
  out.pwLName = pairs.front().pwLeft.name;
  out.pwRName = pairs.front().pwRight.name;

  const int nPads = std::max<int>(1, static_cast<int>(pairs.size()));
  const int nCols = (nPads <= 2) ? nPads : 2;
  const int nRows = (nPads + nCols - 1) / nCols;
  auto c = std::make_unique<TCanvas>("cFieldScores", "Field scores", 1500, std::max(700, 450 * nRows));
  c->Divide(nCols, nRows);

  for (size_t ip = 0; ip < pairs.size(); ++ip) {
    const auto& pair = pairs[ip];
    c->cd(static_cast<int>(ip) + 1);
    gPad->SetGrid();

    hEAnodes[ip]->SetLineColor(kBlue + 1);
    hEAnodes[ip]->SetLineWidth(2);
    hEpwLefts[ip]->SetLineColor(kRed + 1);
    hEpwLefts[ip]->SetLineWidth(2);
    hEpwRights[ip]->SetLineColor(kGreen + 2);
    hEpwRights[ip]->SetLineWidth(2);
    hDiff2s[ip]->SetLineColor(kMagenta + 2);
    hDiff2s[ip]->SetLineWidth(2);

    double ymax = std::max({hEAnodes[ip]->GetMaximum(), hEpwLefts[ip]->GetMaximum(), hEpwRights[ip]->GetMaximum()});
    if (ymax <= 0.0) ymax = 1.0;
    hEAnodes[ip]->SetTitle(Form("%s with adjacent %s / %s; y [cm]; |E| [V/cm]", 
                                pair.anode.name.c_str(), pair.pwLeft.name.c_str(), pair.pwRight.name.c_str()));
    hEAnodes[ip]->SetMaximum(1.10 * ymax);
    hEAnodes[ip]->Draw("HIST");
    hEpwLefts[ip]->Draw("HIST SAME");
    hEpwRights[ip]->Draw("HIST SAME");

    auto leg = new TLegend(0.50, 0.66, 0.88, 0.88);
    leg->SetBit(kCanDelete);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hEAnodes[ip].get(), Form("|E|(x=%s)", pair.anode.name.c_str()), "l");
    leg->AddEntry(hEpwLefts[ip].get(), Form("|E|(x=%s)", pair.pwLeft.name.c_str()), "l");
    leg->AddEntry(hEpwRights[ip].get(), Form("|E|(x=%s)", pair.pwRight.name.c_str()), "l");
    leg->Draw();

    TGaxis* axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(),
                              gPad->GetUxmax(), gPad->GetUymax(),
                              0.0, std::max(1.05 * hDiff2s[ip]->GetMaximum(), 1.e-9), 510, "+L");
    axis->SetLineColor(kMagenta + 2);
    axis->SetLabelColor(kMagenta + 2);
    axis->SetTitleColor(kMagenta + 2);
    axis->SetTitle("field-balance score");
    axis->Draw();

    const double x1 = gPad->GetUxmin();
    const double x2 = gPad->GetUxmax();
    const double y1 = gPad->GetUymin();
    const double y2 = gPad->GetUymax();
    const double scale = (y2 - y1) / std::max(1.05 * hDiff2s[ip]->GetMaximum(), 1.e-9);
    auto gDiff = new TGraph();
    gDiff->SetBit(kCanDelete);
    gDiff->SetLineColor(kMagenta + 2);
    gDiff->SetLineWidth(2);
    int np = 0;
    for (int ib = 1; ib <= hDiff2s[ip]->GetNbinsX(); ++ib) {
      const double xv = hDiff2s[ip]->GetXaxis()->GetBinCenter(ib);
      const double yv = hDiff2s[ip]->GetBinContent(ib);
      if (yv < 0.0) continue;
      gDiff->SetPoint(np++, xv, y1 + yv * scale);
    }
    gDiff->Draw("L SAME");

    DrawInfoPave(0.08, 0.62, 0.46, 0.88, {
      textVolt.Data(),
      textDiam.Data(),
      Form("pair %zu / %zu", ip + 1, pairs.size()),
      Form("exclude |y| < %.4g cm", std::max({pair.anode.r, pair.pwLeft.r, pair.pwRight.r}))
    });
  }

  const fs::path pngPath = outDir / "field_score_maps.png";
  const fs::path pdfPath = outDir / "field_score_maps.pdf";
  c->SaveAs(pngPath.string().c_str());
  c->SaveAs(pdfPath.string().c_str());

  {
    std::ofstream fout(outDir / "field_score_summary.csv");
    if (fout) {
      fout << "config_name,pitchSense_cm,gap_cm,yMin_cm,yMax_cm,S_E_line,S_field_raw,nELine,"
              "n_anode_pairs,example_anode_name,example_pwL_name,example_pwR_name\n";
      fout << confName << ","
           << out.pitchSense << ","
           << out.gap_cm << ","
           << out.yMin_cm << ","
           << out.yMax_cm << ","
           << out.sELine << ","
           << out.sFieldRaw << ","
           << out.nELine << ","
           << pairs.size() << ","
           << out.a0Name << ","
           << out.pwLName << ","
           << out.pwRName << "\n";
    }
  }

  {
    TFile f((outDir / "field_score.root").string().c_str(), "RECREATE");
    TTree t("fieldScoreSummary", "field score summary");
    double pitch_cm = out.pitchSense, gap_cm = out.gap_cm, yMin_cm = out.yMin_cm, yMax_cm = out.yMax_cm;
    double sELine = out.sELine, sFieldRaw = out.sFieldRaw;
    int nELine = out.nELine, nPairs = static_cast<int>(pairs.size());
    char a0Name[64], pwLName[64], pwRName[64];
    std::snprintf(a0Name, sizeof(a0Name), "%s", out.a0Name.c_str());
    std::snprintf(pwLName, sizeof(pwLName), "%s", out.pwLName.c_str());
    std::snprintf(pwRName, sizeof(pwRName), "%s", out.pwRName.c_str());
    t.Branch("pitchSense_cm", &pitch_cm);
    t.Branch("gap_cm", &gap_cm);
    t.Branch("yMin_cm", &yMin_cm);
    t.Branch("yMax_cm", &yMax_cm);
    t.Branch("S_E_line", &sELine);
    t.Branch("S_field_raw", &sFieldRaw);
    t.Branch("nELine", &nELine);
    t.Branch("n_anode_pairs", &nPairs);
    t.Branch("example_anode_name", a0Name, "example_anode_name/C");
    t.Branch("example_pwL_name", pwLName, "example_pwL_name/C");
    t.Branch("example_pwR_name", pwRName, "example_pwR_name/C");
    t.Fill();
    for (auto& h : hEAnodes) h->Write();
    for (auto& h : hEpwLefts) h->Write();
    for (auto& h : hEpwRights) h->Write();
    for (auto& h : hDiff2s) h->Write();
    t.Write();
    f.Close();
  }

  return out;
}


static void DrawInfoPave(double x1, double y1, double x2, double y2,
                         std::initializer_list<std::string> lines) {
  const double textSize = 0.038;
  auto pave = new TPaveText(x1, y1, x2, y2, "NDC");
  pave->SetBit(kCanDelete);
  pave->SetFillColor(kWhite);
  pave->SetFillStyle(1001);
  pave->SetFillColorAlpha(kWhite, 0.88);
  pave->SetLineColor(kGray + 2);
  pave->SetLineWidth(1);
  pave->SetTextColor(kBlack);
  pave->SetTextAlign(12);
  pave->SetTextFont(42);
  for (const auto& line : lines) {
    auto* txt = pave->AddText(line.c_str());
    if (txt) txt->SetTextSize(textSize);
  }
  pave->Draw();
}

static void SaveWireLtPlots(const fs::path& outDir,
                            const std::vector<LtPoint>& pts,
                            const TString& textVolt,
                            const TString& textDiam) {
  fs::create_directories(outDir);

  std::map<std::string, std::vector<LtPoint>> groups;
  for (const auto& p : pts) {
    if (p.wire_uid.empty() || p.wire_uid == "NO_WIRE") continue;
    if (p.L_cm < 0 || p.t_ns < 0) continue;
    groups[p.wire_uid].push_back(p);
  }

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.032);
  latex.SetTextColor(kBlack);
  latex.SetTextAlign(13);

  for (const auto& [wireUid, v] : groups) {
    if (v.empty()) continue;

    std::vector<double> xs(v.size()), ys(v.size());
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    int nTarget = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      xs[i] = v[i].t_ns;
      ys[i] = v[i].L_cm;
      xmin = std::min(xmin, xs[i]);
      xmax = std::max(xmax, xs[i]);
      ymin = std::min(ymin, ys[i]);
      ymax = std::max(ymax, ys[i]);
      if (v[i].is_target) ++nTarget;
    }

    if (!(xmin < xmax)) { xmin -= 1.0; xmax += 1.0; }
    if (!(ymin < ymax)) { ymin -= 0.01; ymax += 0.01; }

    const std::string dispName = v.front().wire_name;
    const double dispX = v.front().wire_x_cm;

    auto c = std::make_unique<TCanvas>(Form("cLt_%s", SanitizeFileName(wireUid).c_str()),
                                       wireUid.c_str(), 900, 700);
    c->SetGrid();

    auto g = std::make_unique<TGraph>((int)v.size(), xs.data(), ys.data());
    g->SetTitle(Form("L-t correlation (%s);Drift time t [ns];Drift length L [cm]", wireUid.c_str()));
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->Draw("AP");
    g->GetXaxis()->SetLimits(xmin - 0.02 * std::abs(xmax - xmin), xmax + 0.02 * std::abs(xmax - xmin));
    g->GetYaxis()->SetRangeUser(ymin - 0.05 * std::abs(ymax - ymin), ymax + 0.05 * std::abs(ymax - ymin));

    auto fit = std::make_unique<TF1>(Form("fLt_%s", SanitizeFileName(wireUid).c_str()), "[0]*x+[1]", xmin, xmax);
    fit->SetLineColor(kRed + 1);
    fit->SetLineWidth(2);
    g->Fit(fit.get(), "Q");

    double ss_res = 0.0, ss_tot = 0.0;
    const double ymean = std::accumulate(ys.begin(), ys.end(), 0.0) / std::max<size_t>(1, ys.size());
    for (size_t i = 0; i < ys.size(); ++i) {
      const double pred = fit->Eval(xs[i]);
      const double res = ys[i] - pred;
      ss_res += res * res;
      const double d = ys[i] - ymean;
      ss_tot += d * d;
    }
    const double rms = std::sqrt(ss_res / std::max<size_t>(1, ys.size()));
    const double r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
    const double eff = (v.empty() ? 0.0 : double(nTarget) / double(v.size()));

    latex.DrawLatex(0.15, 0.88, textVolt);
    latex.DrawLatex(0.15, 0.84, textDiam);
    latex.DrawLatex(0.15, 0.80, Form("Wire name: %s", dispName.c_str()));
    latex.DrawLatex(0.15, 0.76, Form("Wire x = %.6f cm", dispX));
    latex.DrawLatex(0.15, 0.72, Form("UID: %s", wireUid.c_str()));
    latex.DrawLatex(0.15, 0.68, Form("N = %zu", v.size()));
    latex.DrawLatex(0.15, 0.64, Form("RMS = %.5f cm", rms));
    latex.DrawLatex(0.15, 0.60, Form("R^{2} = %.5f", r2));
    latex.DrawLatex(0.15, 0.56, Form("same-wire frac = %.5f", eff));

    const fs::path outPng = outDir / (SanitizeFileName(wireUid) + "_Lt.png");
    c->SaveAs(outPng.string().c_str());
  }
}

static void SaveAllWireLtPlot(const fs::path& outPng,
                              const std::vector<LtPoint>& pts,
                              const TString& textVolt,
                              const TString& textDiam,
                              const std::string& title) {
  std::vector<double> xs, ys;
  xs.reserve(pts.size());
  ys.reserve(pts.size());

  int nTarget = 0;
  for (const auto& p : pts) {
    if (p.L_cm < 0 || p.t_ns < 0) continue;
    xs.push_back(p.t_ns);
    ys.push_back(p.L_cm);
    if (p.is_target) ++nTarget;
  }
  if (xs.empty()) return;

  double xmin = *std::min_element(xs.begin(), xs.end());
  double xmax = *std::max_element(xs.begin(), xs.end());
  double ymin = *std::min_element(ys.begin(), ys.end());
  double ymax = *std::max_element(ys.begin(), ys.end());

  if (!(xmin < xmax)) { xmin -= 1.0; xmax += 1.0; }
  if (!(ymin < ymax)) { ymin -= 0.01; ymax += 0.01; }

  auto c = std::make_unique<TCanvas>("cAllWireLt", "AllWire L-t", 900, 700);
  c->SetGrid();

  auto g = std::make_unique<TGraph>((int)xs.size(), xs.data(), ys.data());
  g->SetTitle(Form("%s;Drift time t [ns];Drift length L [cm]", title.c_str()));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.7);
  g->Draw("AP");
  g->GetXaxis()->SetLimits(xmin - 0.02 * std::abs(xmax - xmin), xmax + 0.02 * std::abs(xmax - xmin));
  g->GetYaxis()->SetRangeUser(ymin - 0.05 * std::abs(ymax - ymin), ymax + 0.05 * std::abs(ymax - ymin));

  auto fit = std::make_unique<TF1>("fAllWireLt", "[0]*x+[1]", xmin, xmax);
  fit->SetLineColor(kRed + 1);
  fit->SetLineWidth(2);
  g->Fit(fit.get(), "Q");

  double ss_res = 0.0, ss_tot = 0.0;
  const double ymean = std::accumulate(ys.begin(), ys.end(), 0.0) / std::max<size_t>(1, ys.size());
  for (size_t i = 0; i < ys.size(); ++i) {
    const double pred = fit->Eval(xs[i]);
    const double res = ys[i] - pred;
    ss_res += res * res;
    const double d = ys[i] - ymean;
    ss_tot += d * d;
  }
  const double rms = std::sqrt(ss_res / std::max<size_t>(1, ys.size()));
  const double r2 = (ss_tot > 0.0) ? (1.0 - ss_res / ss_tot) : 0.0;
  const double eff = double(nTarget) / double(xs.size());

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.032);
  latex.SetTextColor(kBlack);
  latex.SetTextAlign(13);
  latex.DrawLatex(0.15, 0.88, textVolt);
  latex.DrawLatex(0.15, 0.84, textDiam);
  latex.DrawLatex(0.15, 0.80, Form("N = %zu", xs.size()));
  latex.DrawLatex(0.15, 0.76, Form("RMS = %.5f cm", rms));
  latex.DrawLatex(0.15, 0.72, Form("R^{2} = %.5f", r2));
  latex.DrawLatex(0.15, 0.68, Form("same-wire frac = %.5f", eff));

  c->SaveAs(outPng.string().c_str());
}

static double Median(std::vector<double> v) {
  if (v.empty()) return 1.0;
  std::sort(v.begin(), v.end());
  const size_t n = v.size();
  if (n % 2 == 1) return v[n / 2];
  return 0.5 * (v[n / 2 - 1] + v[n / 2]);
}


static void ComputeNormalizedFieldScores(std::vector<ConfigRankingEntry>& entries) {
  std::vector<double> eLine, ex2;
  for (const auto& e : entries) {
    if (e.sELine >= 0.0) eLine.push_back(e.sELine);
    if (e.ex2_mean >= 0.0) ex2.push_back(e.ex2_mean);
  }

  const double medELine = std::max(Median(eLine), 1.e-30);
  const double medEx2   = std::max(Median(ex2),   1.e-30);

  for (auto& e : entries) {
    // Keep legacy fields available, but only E-line and Ex2 contribute to the total score.
    e.sAtransNorm = -1.0;
    e.sAalignNorm = -1.0;
    e.sPWmirrorNorm = -1.0;
    e.sPWalignNorm = -1.0;

    if (e.sELine < 0.0 || e.ex2_mean < 0.0) {
      e.sELineNorm = -1.0;
      e.sEx2Norm = -1.0;
      e.sFieldNorm = -1.0;
      continue;
    }
    e.sELineNorm = e.sELine / medELine;
    e.sEx2Norm = e.ex2_mean / medEx2;
    e.sFieldNorm = e.sELineNorm + e.sEx2Norm;
  }
}




static void RewritePerConfigFieldScoreSummaries(const std::vector<ConfigRankingEntry>& entries) {
  for (const auto& e : entries) {
    if (e.field_score_csv.empty()) continue;

    std::ofstream fout(e.field_score_csv);
    if (!fout) {
      std::cerr << "[warn] Failed to rewrite field-score summary: " << e.field_score_csv << "\n";
      continue;
    }

    fout << "config_name,Vp_V,"
            "S_E_line,Ex2_mean,S_field_raw,"
            "S_E_line_norm,Ex2_norm,S_field_norm,"
            "nELine,field_score_png,field_score_root\n";

    fout << e.conf_name << ","
         << e.vPW << ","
         << e.sELine << ","
         << e.ex2_mean << ","
         << e.sFieldRaw << ","
         << e.sELineNorm << ","
         << e.sEx2Norm << ","
         << e.sFieldNorm << ","
         << e.nELine << ","
         << e.field_score_png << ","
         << e.field_score_root
         << "\n";
  }
}



static void WriteRankingMarkdown(const fs::path& outPath,
                                 const std::vector<ConfigRankingEntry>& entries,
                                 const int topN = 10) {
  std::ofstream fout(outPath);
  if (!fout) {
    std::cerr << "[warn] Failed to write markdown: " << outPath << "\n";
    return;
  }

  fout << "# L-t, Ex boundary, and field-score ranking summary\n\n";
  fout << "- RMS: smaller is better\n";
  fout << "- R^2: larger is better\n";
  fout << "- <Ex^2> at boundary x_eval: smaller is better\n";
  fout << "- Field score: smaller is better\n";
  fout << "- Total field score = norm(S_E_line)+norm(<Ex^2>)\n\n";
  fout << "Total configurations: " << entries.size() << "\n\n";

  auto write_header = [&](const std::string& title, bool withEx = false, bool withField = false) {
    fout << "## " << title << "\n\n";
    if (withEx) {
      fout << "| Rank | Config | Vp [V] | <Ex^2> | max|Ex| | RMS [cm] | R^2 | Ex profile |\n";
      fout << "|---:|:---|---:|---:|---:|---:|---:|:---|\n";
    } else if (withField) {
      fout << "| Rank | Config | Vp [V] | S_field(norm) | S_E_line | <Ex^2> | Field map |\n";
      fout << "|---:|:---|---:|---:|---:|---:|:---|\n";
    } else {
      fout << "| Rank | Config | Vp [V] | RMS [cm] | R^2 | N | monoViol | collectEff | Combo map | AllWire Lt |\n";
      fout << "|---:|:---|---:|---:|---:|---:|---:|---:|:---|:---|\n";
    }
  };

  auto write_row_default = [&](int rank, const ConfigRankingEntry& e) {
    fout << "| " << rank
         << " | `" << e.conf_name << "`"
         << " | " << e.vPW
         << " | " << e.residual_rms
         << " | " << e.r2
         << " | " << e.n_points
         << " | " << e.monotonic_violations
         << " | " << e.collect_eff
         << " | " << e.combo_png
         << " | " << e.allwire_lt_png
         << " |\n";
  };

  auto write_row_ex = [&](int rank, const ConfigRankingEntry& e) {
    fout << "| " << rank
         << " | `" << e.conf_name << "`"
         << " | " << e.vPW
         << " | " << e.ex2_mean
         << " | " << e.ex_abs_max
         << " | " << e.residual_rms
         << " | " << e.r2
         << " | " << e.ex_profile_png
         << " |\n";
  };

  auto write_row_field = [&](int rank, const ConfigRankingEntry& e) {
    fout << "| " << rank
         << " | `" << e.conf_name << "`"
         << " | " << e.vPW
         << " | " << e.sFieldNorm
         << " | " << e.sELine
         << " | " << e.ex2_mean
         << " | " << e.field_score_png
         << " |\n";
  };

  {
    auto v = entries;
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
      return a.residual_rms < b.residual_rms;
    });
    write_header("Top 10 smallest RMS");
    for (int i = 0; i < std::min<int>(topN, (int)v.size()); ++i) write_row_default(i + 1, v[i]);
    fout << "\n";
  }

  {
    auto v = entries;
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
      return a.r2 > b.r2;
    });
    write_header("Top 10 largest R^2");
    for (int i = 0; i < std::min<int>(topN, (int)v.size()); ++i) write_row_default(i + 1, v[i]);
    fout << "\n";
  }

  {
    auto v = entries;
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
      return a.ex2_mean < b.ex2_mean;
    });
    write_header("Top 10 smallest <Ex^2> at x = 0.20 cm", true);
    for (int i = 0; i < std::min<int>(topN, (int)v.size()); ++i) write_row_ex(i + 1, v[i]);
    fout << "\n";
  }

  {
    auto v = entries;
    v.erase(std::remove_if(v.begin(), v.end(), [](const auto& e) { return e.sFieldNorm < 0.0; }), v.end());
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
      return a.sFieldNorm < b.sFieldNorm;
    });
    write_header("Top 10 smallest field score", false, true);
    for (int i = 0; i < std::min<int>(topN, (int)v.size()); ++i) write_row_field(i + 1, v[i]);
    fout << "\n";
  }
}



static void WriteFieldScoreTop10Markdown(const fs::path& outPath,
                                         const std::vector<ConfigRankingEntry>& entries,
                                         const int topN = 10) {
  std::vector<ConfigRankingEntry> v = entries;
  v.erase(std::remove_if(v.begin(), v.end(),
                         [](const auto& e) { return e.sFieldNorm < 0.0; }),
          v.end());
  std::sort(v.begin(), v.end(),
            [](const auto& a, const auto& b) { return a.sFieldNorm < b.sFieldNorm; });

  std::ofstream fout(outPath);
  if (!fout) {
    std::cerr << "[warn] Failed to write field-score markdown: " << outPath << "\n";
    return;
  }

  fout << "# Top field-score configurations\n\n";
  fout << "- Smaller is better.\n";
  fout << "- Total field score = norm(S_E_line)+norm(<Ex^2>)\n\n";
  fout << "| Rank | Config | Vp [V] | S_field(norm) | S_E_line | <Ex^2> | Field map |\n";
  fout << "|---:|:---|---:|---:|---:|---:|:---|\n";

  for (int i = 0; i < std::min<int>(topN, static_cast<int>(v.size())); ++i) {
    const auto& e = v[i];
    fout << "| " << (i + 1)
         << " | `" << e.conf_name << "`"
         << " | " << e.vPW
         << " | " << e.sFieldNorm
         << " | " << e.sELine
         << " | " << e.ex2_mean
         << " | " << e.field_score_png
         << " |\n";
  }
}



static void WriteFieldScoreTop10Csv(const fs::path& outPath,
                                    const std::vector<ConfigRankingEntry>& entries,
                                    const int topN = 10) {
  std::vector<ConfigRankingEntry> v = entries;
  v.erase(std::remove_if(v.begin(), v.end(),
                         [](const auto& e) { return e.sFieldNorm < 0.0; }),
          v.end());
  std::sort(v.begin(), v.end(),
            [](const auto& a, const auto& b) { return a.sFieldNorm < b.sFieldNorm; });

  std::ofstream fout(outPath);
  if (!fout) {
    std::cerr << "[warn] Failed to write field-score csv: " << outPath << "\n";
    return;
  }

  fout << "rank,config_name,Vp_V,S_field_norm,S_E_line,Ex2_mean,field_score_png\n";
  for (int i = 0; i < std::min<int>(topN, static_cast<int>(v.size())); ++i) {
    const auto& e = v[i];
    fout << (i + 1) << ","
         << e.conf_name << ","
         << e.vPW << ","
         << e.sFieldNorm << ","
         << e.sELine << ","
         << e.ex2_mean << ","
         << e.field_score_png
         << "\n";
  }
}


static void SaveExOverlayPlots(const fs::path& outDir,
                               const std::vector<ExProfileForOverlay>& profiles,
                               const std::string& titleBase) {
  fs::create_directories(outDir);
  if (profiles.empty()) return;

  auto save_one = [&](const std::vector<ExProfileForOverlay>& useProfiles,
                      const fs::path& outPng,
                      const std::string& title) {
    if (useProfiles.empty()) return;
    auto c = std::make_unique<TCanvas>(SanitizeFileName(title).c_str(), title.c_str(), 1000, 780);
    c->SetGrid();

    TLegend* leg = new TLegend(0.68, 0.62, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    std::vector<std::unique_ptr<TGraph>> graphs;
    graphs.reserve(useProfiles.size());

    int color = 1;
    bool first = true;
    double xmin = std::numeric_limits<double>::max();
    double xmax = -std::numeric_limits<double>::max();
    double ymin = std::numeric_limits<double>::max();
    double ymax = -std::numeric_limits<double>::max();

    for (const auto& p : useProfiles) {
      if (p.y.empty() || p.ex.empty()) continue;
      xmin = std::min(xmin, *std::min_element(p.y.begin(), p.y.end()));
      xmax = std::max(xmax, *std::max_element(p.y.begin(), p.y.end()));
      ymin = std::min(ymin, *std::min_element(p.ex.begin(), p.ex.end()));
      ymax = std::max(ymax, *std::max_element(p.ex.begin(), p.ex.end()));

      auto g = std::make_unique<TGraph>((int)p.y.size(), p.y.data(), p.ex.data());
      g->SetLineColor(color);
      g->SetLineWidth(2);
      g->SetTitle(Form("%s;Y [cm];E_{x} [V/cm]", title.c_str()));

      if (first) {
        g->Draw("AL");
        first = false;
      } else {
        g->Draw("L SAME");
      }
      leg->AddEntry(g.get(), p.label.c_str(), "l");
      graphs.push_back(std::move(g));

      ++color;
      if (color == 5 || color == 10) ++color;
    }

    if (!first) {
      if (!(xmin < xmax)) { xmin -= 0.1; xmax += 0.1; }
      if (!(ymin < ymax)) { ymin -= 1.; ymax += 1.; }
      graphs.front()->GetXaxis()->SetLimits(xmin, xmax);
      graphs.front()->GetYaxis()->SetRangeUser(ymin - 0.05 * std::abs(ymax - ymin),
                                               ymax + 0.05 * std::abs(ymax - ymin));
      auto l0 = new TLine(xmin, 0.0, xmax, 0.0);
      l0->SetBit(kCanDelete);
      l0->SetLineStyle(2);
      l0->SetLineColor(kRed + 1);
      l0->Draw("same");
      leg->Draw();
      c->SaveAs(outPng.string().c_str());
    }
  };

  save_one(profiles, outDir / "Ex_vs_Y_x2mm_overlay_all.png", titleBase + " (all)");

  auto sorted = profiles;
  std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) { return a.vPW < b.vPW; });
  if ((int)sorted.size() > 10) sorted.resize(10);
  save_one(sorted, outDir / "Ex_vs_Y_x2mm_overlay_first10.png", titleBase + " (first10)");
}

static void SaveAllWireLtSummaryPlots(const fs::path& outDir,
                                      const std::vector<ConfigRankingEntry>& rankingEntries,
                                      const int topN = 10) {
  fs::create_directories(outDir);
  if (rankingEntries.empty()) return;

  auto write_md = [&](const std::string& fname, const std::vector<ConfigRankingEntry>& v, const std::string& title) {
    std::ofstream fout(outDir / fname);
    if (!fout) return;
    fout << "# " << title << "\n\n";
    fout << "| Rank | Config | Vp [V] | RMS [cm] | R^2 | AllWire Lt |\n";
    fout << "|---:|:---|---:|---:|---:|:---|\n";
    for (size_t i = 0; i < v.size(); ++i) {
      fout << "| " << (i + 1)
           << " | `" << v[i].conf_name << "`"
           << " | " << v[i].vPW
           << " | " << v[i].residual_rms
           << " | " << v[i].r2
           << " | " << v[i].allwire_lt_png
           << " |\n";
    }
  };

  {
    auto v = rankingEntries;
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.residual_rms < b.residual_rms; });
    if ((int)v.size() > topN) v.resize(topN);
    write_md("AllWire_Lt_RMS_Top10.md", v, "AllWire Lt Summary by smallest RMS");
  }

  {
    auto v = rankingEntries;
    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.r2 > b.r2; });
    if ((int)v.size() > topN) v.resize(topN);
    write_md("AllWire_Lt_R2_Top10.md", v, "AllWire Lt Summary by largest R^2");
  }
}

static void SaveRmsVsPotentialPlot(const fs::path& outDir,
                                   const std::vector<ConfigRankingEntry>& entries,
                                   const std::string& title = "RMS vs Potential-wire voltage") {
  fs::create_directories(outDir);
  if (entries.empty()) return;

  std::vector<ConfigRankingEntry> v = entries;
  std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
    return a.vPW < b.vPW;
  });

  std::vector<double> xs, ys;
  xs.reserve(v.size());
  ys.reserve(v.size());

  for (const auto& e : v) {
    if (e.n_points <= 0) continue;
    if (e.residual_rms < 0) continue;
    xs.push_back(e.vPW);
    ys.push_back(e.residual_rms);
  }
  if (xs.empty()) return;

  double xmin = *std::min_element(xs.begin(), xs.end());
  double xmax = *std::max_element(xs.begin(), xs.end());
  double ymin = *std::min_element(ys.begin(), ys.end());
  double ymax = *std::max_element(ys.begin(), ys.end());

  if (!(xmin < xmax)) { xmin -= 1.0; xmax += 1.0; }
  if (!(ymin < ymax)) { ymin -= 0.001; ymax += 0.001; }

  auto it_best = std::min_element(v.begin(), v.end(), [](const auto& a, const auto& b) {
    return a.residual_rms < b.residual_rms;
  });

  auto c = std::make_unique<TCanvas>("cRmsVsPotential", "RMS vs Vp", 900, 700);
  c->SetGrid();

  auto g = std::make_unique<TGraph>((int)xs.size(), xs.data(), ys.data());
  g->SetTitle(Form("%s;Potential wire voltage V_{p} [V];L-t residual RMS [cm]", title.c_str()));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(2);
  g->Draw("APL");

  g->GetXaxis()->SetLimits(xmin, xmax);
  g->GetYaxis()->SetRangeUser(ymin - 0.05 * std::abs(ymax - ymin),
                              ymax + 0.08 * std::abs(ymax - ymin));

  std::unique_ptr<TGraph> gBest;
  std::unique_ptr<TLatex> labelBest;

  if (it_best != v.end() && it_best->n_points > 0 && it_best->residual_rms >= 0) {
    double xbest = it_best->vPW;
    double ybest = it_best->residual_rms;

    gBest = std::make_unique<TGraph>(1, &xbest, &ybest);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2.0);
    gBest->SetMarkerColor(kRed + 1);
    gBest->Draw("P SAME");

    labelBest = std::make_unique<TLatex>();
    labelBest->SetTextSize(0.028);
    labelBest->SetTextColor(kRed + 1);
    labelBest->DrawLatex(xbest, ybest,
                         Form("  Best (Vp=%.0f V, RMS=%.5f cm)", xbest, ybest));
  }

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.032);
  latex.SetTextAlign(13);
  latex.DrawLatex(0.15, 0.88, Form("N config = %zu", xs.size()));

  if (it_best != v.end() && it_best->n_points > 0 && it_best->residual_rms >= 0) {
    latex.DrawLatex(0.15, 0.84,
                    Form("Best: Vp = %.0f V, RMS = %.5f cm", it_best->vPW, it_best->residual_rms));
    latex.DrawLatex(0.15, 0.80,
                    Form("R^{2} = %.5f, collectEff = %.5f", it_best->r2, it_best->collect_eff));
  }

  c->SaveAs((outDir / "RMS_vs_Vp.png").string().c_str());

  std::ofstream fout(outDir / "RMS_vs_Vp.csv");
  if (fout) {
    fout << "Vp_V,residual_rms_cm,config_name,n_points,r2,collect_eff\n";
    for (const auto& e : v) {
      if (e.n_points <= 0 || e.residual_rms < 0) continue;
      fout << e.vPW << ","
           << e.residual_rms << ","
           << e.conf_name << ","
           << e.n_points << ","
           << e.r2 << ","
           << e.collect_eff << "\n";
    }
  }
}

static void SaveEx2VsPotentialPlot(const fs::path& outDir,
                                   const std::vector<ConfigRankingEntry>& entries,
                                   const std::string& title = "Mean Ex^{2} vs Potential-wire voltage") {
  fs::create_directories(outDir);
  if (entries.empty()) return;

  std::vector<ConfigRankingEntry> v = entries;
  std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
    return a.vPW < b.vPW;
  });

  std::vector<double> xs, ys;
  xs.reserve(v.size());
  ys.reserve(v.size());

  const ConfigRankingEntry* bestEntry = nullptr;
  for (const auto& e : v) {
    if (e.n_points <= 0) continue;
    if (e.ex2_mean < 0) continue;
    xs.push_back(e.vPW);
    ys.push_back(e.ex2_mean);

    if (!bestEntry || e.ex2_mean < bestEntry->ex2_mean) {
      bestEntry = &e;
    }
  }
  if (xs.empty()) return;

  double xmin = *std::min_element(xs.begin(), xs.end());
  double xmax = *std::max_element(xs.begin(), xs.end());
  double ymin = *std::min_element(ys.begin(), ys.end());
  double ymax = *std::max_element(ys.begin(), ys.end());

  if (!(xmin < xmax)) { xmin -= 1.0; xmax += 1.0; }
  if (!(ymin < ymax)) { ymin *= 0.95; ymax *= 1.05; }

  auto c = std::make_unique<TCanvas>("cEx2VsPotential", "Ex2 vs Vp", 900, 700);
  c->SetGrid();

  auto g = std::make_unique<TGraph>((int)xs.size(), xs.data(), ys.data());
  g->SetTitle(Form("%s;Potential wire voltage V_{p} [V];<E_{x}^{2}> [(V/cm)^{2}]",
                   title.c_str()));
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(2);
  g->Draw("APL");

  g->GetXaxis()->SetLimits(xmin, xmax);
  g->GetYaxis()->SetRangeUser(ymin - 0.05 * std::abs(ymax - ymin),
                              ymax + 0.08 * std::abs(ymax - ymin));
  // =========================================================
  // Quadratic fit for raw <Ex^2> vs Vp
  // f(Vp) = a Vp^2 + b Vp + c
  // Vp_min = -b / (2a)
  // =========================================================
  double fit_a = std::numeric_limits<double>::quiet_NaN();
  double fit_b = std::numeric_limits<double>::quiet_NaN();
  double fit_c = std::numeric_limits<double>::quiet_NaN();
  double vp_min_fit = std::numeric_limits<double>::quiet_NaN();
  double ex2_min_fit = std::numeric_limits<double>::quiet_NaN();

  std::unique_ptr<TF1> fQuad;
  std::unique_ptr<TGraph> gFitMin;

  if (xs.size() >= 3) {
    const double xfit_min = *std::min_element(xs.begin(), xs.end());
    const double xfit_max = *std::max_element(xs.begin(), xs.end());

    fQuad = std::make_unique<TF1>(
        "fQuadEx2",
        "[0]*x*x + [1]*x + [2]",
        xfit_min,
        xfit_max
    );

    fQuad->SetParameters(1.0, 1.0, ys.front());
    fQuad->SetLineColor(kRed + 1);
    fQuad->SetLineWidth(2);
    fQuad->SetLineStyle(1);

    g->Fit(fQuad.get(), "Q");

    fit_a = fQuad->GetParameter(0);
    fit_b = fQuad->GetParameter(1);
    fit_c = fQuad->GetParameter(2);

    if (std::abs(fit_a) > 1.0e-30) {
      vp_min_fit = -fit_b / (2.0 * fit_a);
      ex2_min_fit = fQuad->Eval(vp_min_fit);
    }

    // fit function
    fQuad->Draw("SAME");

    // draw Vmin only if it is inside the fit range
    if (std::isfinite(vp_min_fit) &&
        std::isfinite(ex2_min_fit) &&
        vp_min_fit >= xfit_min &&
        vp_min_fit <= xfit_max) {

      gFitMin = std::make_unique<TGraph>(1, &vp_min_fit, &ex2_min_fit);
      gFitMin->SetMarkerStyle(34);
      gFitMin->SetMarkerSize(2.2);
      gFitMin->SetMarkerColor(kBlue + 1);
      gFitMin->Draw("P SAME");

      auto lineVmin = new TLine(vp_min_fit,
                                ymin - 0.05 * std::abs(ymax - ymin),
                                vp_min_fit,
                                ex2_min_fit);
      lineVmin->SetBit(kCanDelete);
      lineVmin->SetLineColor(kBlue + 1);
      lineVmin->SetLineStyle(2);
      lineVmin->SetLineWidth(2);
      lineVmin->Draw("SAME");
    }

    // show fit equation and Vmin on the figure
    auto paveFit = new TPaveText(0.15, 0.56, 0.58, 0.76, "NDC");
    paveFit->SetBit(kCanDelete);
    paveFit->SetFillColor(kWhite);
    paveFit->SetFillStyle(1001);
    paveFit->SetFillColorAlpha(kWhite, 0.85);
    paveFit->SetLineColor(kGray + 2);
    paveFit->SetTextAlign(12);
    paveFit->SetTextFont(42);

    paveFit->AddText("Quadratic fit to raw <E_{x}^{2}>:");
    paveFit->AddText(Form("f(V_{p}) = %.3e V_{p}^{2} %+ .3e V_{p} %+ .3e",
                          fit_a, fit_b, fit_c));

    if (std::isfinite(vp_min_fit) && std::isfinite(ex2_min_fit)) {
      paveFit->AddText(Form("V_{min} = %.2f V", vp_min_fit));
      paveFit->AddText(Form("f(V_{min}) = %.3e", ex2_min_fit));
    } else {
      paveFit->AddText("V_{min}: not available");
    }

    paveFit->Draw("SAME");
  }
  std::unique_ptr<TGraph> gBest;
  if (bestEntry) {
    double xbest = bestEntry->vPW;
    double ybest = bestEntry->ex2_mean;
    gBest = std::make_unique<TGraph>(1, &xbest, &ybest);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2.0);
    gBest->SetMarkerColor(kRed + 1);
    gBest->Draw("P SAME");

    TLatex label;
    label.SetTextSize(0.028);
    label.SetTextColor(kRed + 1);
    label.DrawLatex(xbest, ybest,
                    Form("  Best (Vp=%.0f V, <Ex^{2}>=%.6e)", xbest, ybest));
  }

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.032);
  latex.SetTextAlign(13);
  latex.DrawLatex(0.15, 0.88, Form("N config = %zu", xs.size()));
  latex.DrawLatex(0.15, 0.84, Form("x = %.3f cm", v.front().ex_boundary_x_cm));

  if (bestEntry) {
    latex.DrawLatex(0.15, 0.80,
                    Form("Best: Vp = %.0f V, <Ex^{2}> = %.6e",
                         bestEntry->vPW, bestEntry->ex2_mean));
    latex.DrawLatex(0.15, 0.76,
                    Form("max|Ex| = %.6e", bestEntry->ex_abs_max));
  }

  c->SaveAs((outDir / "Ex2_vs_Vp.png").string().c_str());

  std::ofstream fout(outDir / "Ex2_vs_Vp.csv");
  if (fout) {
    fout << "Vp_V,ex2_mean,ex_abs_max,config_name,n_points,r2,collect_eff,x_eval_cm\n";
    for (const auto& e : v) {
      if (e.n_points <= 0 || e.ex2_mean < 0) continue;
      fout << e.vPW << ","
           << e.ex2_mean << ","
           << e.ex_abs_max << ","
           << e.conf_name << ","
           << e.n_points << ","
           << e.r2 << ","
           << e.collect_eff << ","
           << e.ex_boundary_x_cm << "\n";
    }
  }
}


static void SaveFieldScoreVsPotentialPlot(const fs::path& outDir,
                                          const std::vector<ConfigRankingEntry>& entries,
                                          const std::string& title = "Field scores vs Potential-wire voltage") {
  fs::create_directories(outDir);
  if (entries.empty()) return;

  std::vector<ConfigRankingEntry> v = entries;
  std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
    return a.vPW < b.vPW;
  });

  std::vector<double> x, yELine, yEx2, yFieldNorm, yELineNorm, yEx2Norm;
  for (const auto& e : v) {
    x.push_back(e.vPW);
    yELine.push_back(e.sELine);
    yEx2.push_back(e.ex2_mean);
    yFieldNorm.push_back(e.sFieldNorm);
    yELineNorm.push_back(e.sELineNorm);
    yEx2Norm.push_back(e.sEx2Norm);
  }

  auto makeGraph = [&](const std::vector<double>& ys, const char* name, const char* ttl) {
    auto g = std::make_unique<TGraph>((int)x.size(), x.data(), ys.data());
    g->SetName(name);
    g->SetTitle(ttl);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetLineWidth(2);
    return g;
  };

  auto g1 = makeGraph(yELine, "gELine",
                      "Anode-PW field-balance score (raw);Potential-wire voltage V_{p} [V];S_{E,line}");
  auto g2 = makeGraph(yEx2, "gEx2",
                      "Boundary <E_{x}^{2}> (raw);Potential-wire voltage V_{p} [V];<E_{x}^{2}> [(V/cm)^{2}]");
  auto g3 = makeGraph(yELineNorm, "gELineNorm",
                      "Anode-PW field-balance score (normalized);Potential-wire voltage V_{p} [V];Normalized S_{E,line}");
  auto g4 = makeGraph(yEx2Norm, "gEx2Norm",
                      "Boundary <E_{x}^{2}> (normalized);Potential-wire voltage V_{p} [V];Normalized <E_{x}^{2}>");
  auto g5 = makeGraph(yFieldNorm, "gFieldNorm",
                      "Total field score = norm(S_{E,line})+norm(<E_{x}^{2}>);Potential-wire voltage V_{p} [V];Normalized total field score");

  auto c = std::make_unique<TCanvas>("cFieldScoreVsVp", title.c_str(), 1400, 900);
  c->Divide(2, 2);
  std::vector<TGraph*> gs = {g1.get(), g2.get(), g3.get(), g4.get()};
  for (size_t i = 0; i < gs.size(); ++i) {
    c->cd((int)i + 1);
    gPad->SetGrid();
    gs[i]->Draw("APL");
  }
  c->SaveAs((outDir / "FieldScores_vs_Vp.png").string().c_str());

  auto cNorm = std::make_unique<TCanvas>("cFieldScoreVsVpNorm", "Field score normalized components", 1400, 900);
  cNorm->Divide(1, 3);
  std::vector<TGraph*> gsn = {g3.get(), g4.get(), g5.get()};
  for (size_t i = 0; i < gsn.size(); ++i) {
    cNorm->cd((int)i + 1);
    gPad->SetGrid();
    gsn[i]->Draw("APL");
  }
  cNorm->SaveAs((outDir / "FieldScores_normalized_vs_Vp.png").string().c_str());

  auto cTot = std::make_unique<TCanvas>("cTotalFieldScoreVsVp", "Total field score vs Vp", 1000, 760);
  cTot->SetGrid();
  g5->SetMarkerStyle(20);
  g5->SetMarkerSize(1.0);
  g5->SetLineWidth(2);
  g5->Draw("APL");
  TLatex latexTot;
  latexTot.SetNDC();
  latexTot.SetTextSize(0.032);
  latexTot.SetTextAlign(13);
  latexTot.DrawLatex(0.15, 0.88, "Smaller is better");
  latexTot.DrawLatex(0.15, 0.84,
    "S_{field,norm}=norm(S_{E,line})+norm(<E_{x}^{2}>)");
  auto itBest = std::min_element(v.begin(), v.end(), [](const auto& a, const auto& b) {
    return a.sFieldNorm < b.sFieldNorm;
  });
  if (itBest != v.end()) {
    const double xbest = itBest->vPW;
    const double ybest = itBest->sFieldNorm;
    auto gBest = std::make_unique<TGraph>(1, &xbest, &ybest);
    gBest->SetMarkerStyle(29);
    gBest->SetMarkerSize(2.0);
    gBest->SetMarkerColor(kRed + 1);
    gBest->Draw("P SAME");
    latexTot.SetTextColor(kRed + 1);
    latexTot.DrawLatex(0.15, 0.80,
      Form("Best: V_{p}=%.0f V, S_{field,norm}=%.6g", itBest->vPW, itBest->sFieldNorm));
    latexTot.SetTextColor(kBlack);
  }
  cTot->SaveAs((outDir / "TotalFieldScore_vs_Vp.png").string().c_str());
  cTot->SaveAs((outDir / "TotalFieldScore_vs_Vp.pdf").string().c_str());

  std::ofstream foutTot(outDir / "TotalFieldScore_vs_Vp.csv");
  if (foutTot) {
    foutTot << "Vp_V,S_field_norm,config_name\n";
    for (const auto& e : v) {
      foutTot << e.vPW << "," << e.sFieldNorm << "," << e.conf_name << "\n";
    }
  }

  std::ofstream fout(outDir / "FieldScores_vs_Vp.csv");
  if (fout) {
    fout << "Vp_V,S_E_line,Ex2_mean,S_E_line_norm,Ex2_norm,S_field_norm,config_name\n";
    for (const auto& e : v) {
      fout << e.vPW << ","
           << e.sELine << ","
           << e.ex2_mean << ","
           << e.sELineNorm << ","
           << e.sEx2Norm << ","
           << e.sFieldNorm << ","
           << e.conf_name << "\n";
    }
  }
}


// ---- main ----------------------------------------------------
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const char* timeHistFile =
      (argc > 1) ? argv[1] : "./root/run_20260302_161622/analysis_L0_prim_exact/t_hist_nt.csv";

  std::vector<double> vCat_list   = {-1500.0};
  std::vector<double> vAnode_list = {0};
  std::vector<double> vPW_list    = {
      0, -5, -10, -15, -20, -25, -30, -35, -40, -45, -50, -55, -60, -65, -70,
      -75, -80, -85, -90, -95, -100, -105, -110, -115, -120, -125, -130, -135,
      -140, -145, -150, -155, -160, -165, -170, -175, -180, -185, -190, -195,
      -200, -205, -210, -215, -220, -225, -230
  };/*, -235, -240, -245, -250, -255,
      -260, -265, -270, -275, -280, -285, -290, -295, -300*/

  std::vector<double> dAnode_list = {0.0020};
  std::vector<double> dPW_list    = {0.0050};

#if defined(GEOM_PLATE)
  std::vector<double> dCat_list = {0.0};
  std::string geom_tag = "plate";
#elif defined(GEOM_WIRE)
  std::vector<double> dCat_list = {0.0050};
  std::string geom_tag = "wirecath";
#endif

  fs::path csvPath(timeHistFile);
  // 入力元 run ディレクトリは確認用に残す
  fs::path inputRunDir = GuessRunDirFromCsv(csvPath);
  // 出力先をプログラム内で固定する
  // ディレクトリ名に ap / pap / a などを入れておくと後から分かりやすい
  fs::path outBaseDir = fs::path("opt_20260426_pap");
  fs::path sweepBaseDir = outBaseDir / geom_tag;
  fs::create_directories(sweepBaseDir);

  fs::path emapsSummaryDir = sweepBaseDir / "emaps_summary";
  fs::path exProfileSummaryDir = sweepBaseDir / "ExProfile_x2mm_Summary";
  fs::path allWireLtSummaryDir = sweepBaseDir / "AllWire_Lt_Summary";
  fs::path fieldScoreSummaryDir = sweepBaseDir / "FieldScore_Summary";
  fs::path totalScoreDir = sweepBaseDir / "TotalScore";
  fs::create_directories(emapsSummaryDir);
  fs::create_directories(exProfileSummaryDir);
  fs::create_directories(allWireLtSummaryDir);
  fs::create_directories(fieldScoreSummaryDir);
  fs::create_directories(totalScoreDir);

  TimeHist th;
  LoadTimeHistogram(timeHistFile, th);

  MediumMagboltz gas;
  gas.LoadGasFile("ic4H10_100_0.1atm.gas");

  std::vector<ConfigRankingEntry> rankingEntries;
  std::vector<ExProfileForOverlay> exProfiles;

  int config_count = 0;
  for (double vCat : vCat_list) {
  for (double vAnode : vAnode_list) {
  for (double vPW : vPW_list) {
  for (double dCat : dCat_list) {
  for (double dAnode : dAnode_list) {
  for (double dPW : dPW_list) {
    config_count++;

    const double rCat   = dCat   / 2.0;
    const double rAnode = dAnode / 2.0;
    const double rPW    = dPW    / 2.0;

    char conf_name[128];
    TString textVolt = TString::Format("Voltages: Vc=%.0fV, Va=%.0fV, Vp=%.0fV", vCat, vAnode, vPW);
    TString textDiam;

#if defined(GEOM_PLATE)
    std::snprintf(conf_name, sizeof(conf_name),
                  "vC_%.0f_vA_%.0f_vP_%.0f_dA_%.0f_dP_%.0f",
                  vCat, vAnode, vPW, dAnode * 1e4, dPW * 1e4);
    textDiam = TString::Format("Diameters: da=%.0fum, dp=%.0fum", dAnode * 1e4, dPW * 1e4);
#elif defined(GEOM_WIRE)
    std::snprintf(conf_name, sizeof(conf_name),
                  "vC_%.0f_vA_%.0f_vP_%.0f_dC_%.0f_dA_%.0f_dP_%.0f",
                  vCat, vAnode, vPW, dCat * 1e4, dAnode * 1e4, dPW * 1e4);
    textDiam = TString::Format("Diameters: dc=%.0fum, da=%.0fum, dp=%.0fum",
                               dCat * 1e4, dAnode * 1e4, dPW * 1e4);
#endif

    std::printf("\n======================================================\n");
    std::printf("[Sweep %d] Config: %s\n", config_count, conf_name);
    std::printf("======================================================\n");

    fs::path baseOutDir = sweepBaseDir / conf_name;
    std::string figDir = (baseOutDir / "png_maps").string();
    fs::create_directories(fs::path(figDir) / "waveform");
    fs::path trackDir = baseOutDir / "track";
    fs::path exDir = baseOutDir / "ex_boundary";
    fs::path ltPlotDir = baseOutDir / "lt_wire_plots";
    fs::path fieldScoreDir = baseOutDir / "field_scores";
    fs::create_directories(trackDir);
    fs::create_directories(exDir);
    fs::create_directories(ltPlotDir);
    fs::create_directories(fieldScoreDir);

    Detector::Geometry geo;
#if defined(GEOM_PLATE)
    geo = Detector::PAP_PlaneCathode_Periodic(0.60, 0.50, 0.0, rAnode, rPW, vAnode, vPW, vCat);
#elif defined(GEOM_WIRE)
    geo = Detector::PAP_WireCathode_Periodic(0.60, 0.50, 0.0, 0.20, 0.1, rAnode, rPW, rCat,
                                             vAnode, vPW, vCat);
#endif

    auto comp = Detector::BuildField(geo);
    comp->SetMedium(&gas);

    const double sensePitchX = DetectSensePitchX(geo);
    const int seedMarginCells = 1;
    const auto [xSeedMin, xSeedMax] = ComputeSeedXRange(geo, seedMarginCells);

    std::printf("[seed] sense pitch = %.6f cm, xSeedMin = %.6f cm, xSeedMax = %.6f cm, marginCells = %d\n",
                sensePitchX, xSeedMin, xSeedMax, seedMarginCells);

    std::string rootOutPath = (trackDir / "track_results.root").string();
    auto fOut = std::make_unique<TFile>(rootOutPath.c_str(), "RECREATE");
    TTree* tree = new TTree("tree", "Track Reconstruction Stats");
    double t_x_true = 0.0, t_x_reco = 0.0, t_diff = 0.0, t_abs_diff = 0.0;
    double t_a_true = 0.0, t_b_true = 0.0, t_a_reco = 0.0, t_b_reco = 0.0, t_chi2 = 0.0;
    int t_event_id = 0;
    tree->Branch("event_id", &t_event_id);
    tree->Branch("a_true", &t_a_true);
    tree->Branch("b_true", &t_b_true);
    tree->Branch("x_true", &t_x_true);
    tree->Branch("a_reco", &t_a_reco);
    tree->Branch("b_reco", &t_b_reco);
    tree->Branch("x_reco", &t_x_reco);
    tree->Branch("diff_x", &t_diff);
    tree->Branch("abs_diff", &t_abs_diff);
    tree->Branch("chi2", &t_chi2);

    std::string driftRootOutPath = (trackDir / "drift_info.root").string();
    auto fDriftOut = std::make_unique<TFile>(driftRootOutPath.c_str(), "RECREATE");
    TTree* tDrift = new TTree("driftTree", "Electron Drift Info");
    int d_event_id = 0;
    double d_x_start = 0.0, d_y_start = 0.0, d_x_end = 0.0, d_y_end = 0.0, d_t_drift = 0.0, d_dist = 0.0;
    char d_wire_name[64];
    char d_wire_uid[128];
    tDrift->Branch("event_id", &d_event_id);
    tDrift->Branch("x_start", &d_x_start);
    tDrift->Branch("y_start", &d_y_start);
    tDrift->Branch("x_end", &d_x_end);
    tDrift->Branch("y_end", &d_y_end);
    tDrift->Branch("t_drift", &d_t_drift);
    tDrift->Branch("dist", &d_dist);
    tDrift->Branch("wire_name", d_wire_name, "wire_name/C");
    tDrift->Branch("wire_uid", d_wire_uid, "wire_uid/C");

    std::string gainRootOutPath = (trackDir / "gain_info.root").string();
    auto fGainOut = std::make_unique<TFile>(gainRootOutPath.c_str(), "RECREATE");
    TTree* tGain = new TTree("gainTree", "Avalanche Gain Info");
    int g_event_id = 0, g_wire_type = 0, g_n_primary = 0, g_n_total = 0;
    double g_wire_tag = 0.0, g_x_wire = 0.0, g_y_wire = 0.0, g_gain = 0.0;
    char g_wire_name[64];
    char g_wire_uid[128];
    tGain->Branch("event_id", &g_event_id);
    tGain->Branch("wire_type", &g_wire_type);
    tGain->Branch("wire_tag", &g_wire_tag);
    tGain->Branch("x_wire", &g_x_wire);
    tGain->Branch("y_wire", &g_y_wire);
    tGain->Branch("n_primary", &g_n_primary);
    tGain->Branch("n_total", &g_n_total);
    tGain->Branch("gain", &g_gain);
    tGain->Branch("wire_name", g_wire_name, "wire_name/C");
    tGain->Branch("wire_uid", g_wire_uid, "wire_uid/C");

    std::string verifyRootOutPath = (trackDir / "verify_info.root").string();
    auto fVerifyOut = std::make_unique<TFile>(verifyRootOutPath.c_str(), "RECREATE");
    TTree* tVerify = new TTree("verifyTree", "Shockley-Ramo Verification");
    int v_event_id = 0;
    char v_wire_name[64];
    double v_x_start = 0.0, v_y_start = 0.0, v_x_end = 0.0, v_y_end = 0.0, v_wp_start = 0.0, v_wp_end = 0.0, v_wp_diff = 0.0;
    double v_raw_charge_e = 0.0, v_raw_charge_i = 0.0, v_raw_charge_total = 0.0;
    tVerify->Branch("event_id", &v_event_id);
    tVerify->Branch("wire_name", v_wire_name, "wire_name/C");
    tVerify->Branch("x_start", &v_x_start);
    tVerify->Branch("y_start", &v_y_start);
    tVerify->Branch("x_end", &v_x_end);
    tVerify->Branch("y_end", &v_y_end);
    tVerify->Branch("wp_start", &v_wp_start);
    tVerify->Branch("wp_end", &v_wp_end);
    tVerify->Branch("wp_diff", &v_wp_diff);
    tVerify->Branch("raw_charge_e", &v_raw_charge_e);
    tVerify->Branch("raw_charge_i", &v_raw_charge_i);
    tVerify->Branch("raw_charge_total", &v_raw_charge_total);

    std::string evalRootOutPath = (trackDir / "lt_eval.root").string();
    auto fEvalOut = std::make_unique<TFile>(evalRootOutPath.c_str(), "RECREATE");

    TTree* tLt = new TTree("ltTree", "L-t points");
    double lt_t_ns = 0.0, lt_L_cm = 0.0, lt_x0_cm = 0.0, lt_y0_cm = 0.0, lt_x1_cm = 0.0, lt_y1_cm = 0.0;
    double lt_path_cm = 0.0, lt_dx_end_cm = 0.0, lt_wire_x_cm = 0.0, lt_wire_y_cm = 0.0;
    int lt_wire_type = 0, lt_is_target = 0, lt_status = 0;
    char lt_wire_name[64];
    char lt_wire_uid[128];
    char lt_target_wire_name[64];
    char lt_target_wire_uid[128];
    tLt->Branch("t_ns", &lt_t_ns);
    tLt->Branch("L_cm", &lt_L_cm);
    tLt->Branch("x0_cm", &lt_x0_cm);
    tLt->Branch("y0_cm", &lt_y0_cm);
    tLt->Branch("x1_cm", &lt_x1_cm);
    tLt->Branch("y1_cm", &lt_y1_cm);
    tLt->Branch("path_cm", &lt_path_cm);
    tLt->Branch("dx_end_cm", &lt_dx_end_cm);
    tLt->Branch("wire_x_cm", &lt_wire_x_cm);
    tLt->Branch("wire_y_cm", &lt_wire_y_cm);
    tLt->Branch("wire_type", &lt_wire_type);
    tLt->Branch("is_target", &lt_is_target);
    tLt->Branch("status", &lt_status);
    tLt->Branch("wire_name", lt_wire_name, "wire_name/C");
    tLt->Branch("wire_uid", lt_wire_uid, "wire_uid/C");
    tLt->Branch("target_wire_name", lt_target_wire_name, "target_wire_name/C");
    tLt->Branch("target_wire_uid", lt_target_wire_uid, "target_wire_uid/C");

    TTree* tLtSummary = new TTree("ltSummary", "L-t fit summary");
    double s_vcat = 0.0, s_vanode = 0.0, s_vpw = 0.0;
    int s_npoints = 0, s_monotonic_violations = 0;
    double s_slope = 0.0, s_intercept = 0.0, s_r2 = 0.0, s_residual_rms = 0.0;
    double s_mean_dx_end = 0.0, s_mean_path_excess = 0.0, s_collect_eff = 0.0;
    tLtSummary->Branch("Vc", &s_vcat);
    tLtSummary->Branch("Va", &s_vanode);
    tLtSummary->Branch("Vp", &s_vpw);
    tLtSummary->Branch("n_points", &s_npoints);
    tLtSummary->Branch("slope", &s_slope);
    tLtSummary->Branch("intercept", &s_intercept);
    tLtSummary->Branch("r2", &s_r2);
    tLtSummary->Branch("residual_rms", &s_residual_rms);
    tLtSummary->Branch("monotonic_violations", &s_monotonic_violations);
    tLtSummary->Branch("mean_dx_end", &s_mean_dx_end);
    tLtSummary->Branch("mean_path_excess", &s_mean_path_excess);
    tLtSummary->Branch("collect_eff", &s_collect_eff);

    TTree* tFld = new TTree("fieldLineDriftTree",
                            "Drift info for the same seeds used in field-line plotting");
    int fld_seed_id = -1;
    double fld_x0_cm = 0.0, fld_y0_cm = 0.0, fld_z0_cm = 0.0;
    double fld_x1_cm = 0.0, fld_y1_cm = 0.0, fld_z1_cm = 0.0;
    double fld_t_drift_ns = 0.0, fld_drift_path_cm = 0.0;
    double fld_wire_x_cm = 0.0, fld_wire_y_cm = 0.0;
    int fld_wire_type = -1;
    char fld_wire_name[64];
    char fld_wire_uid[128];
    int fld_status = 0;
    tFld->Branch("seed_id", &fld_seed_id);
    tFld->Branch("x0_cm", &fld_x0_cm);
    tFld->Branch("y0_cm", &fld_y0_cm);
    tFld->Branch("z0_cm", &fld_z0_cm);
    tFld->Branch("x1_cm", &fld_x1_cm);
    tFld->Branch("y1_cm", &fld_y1_cm);
    tFld->Branch("z1_cm", &fld_z1_cm);
    tFld->Branch("t_drift_ns", &fld_t_drift_ns);
    tFld->Branch("drift_path_cm", &fld_drift_path_cm);
    tFld->Branch("wire_x_cm", &fld_wire_x_cm);
    tFld->Branch("wire_y_cm", &fld_wire_y_cm);
    tFld->Branch("wire_type", &fld_wire_type);
    tFld->Branch("wire_name", fld_wire_name, "wire_name/C");
    tFld->Branch("wire_uid", fld_wire_uid, "wire_uid/C");
    tFld->Branch("status", &fld_status);

    Sensor sensor;
    sensor.AddComponent(comp.get());
    sensor.SetArea(geo.xmin, geo.ymin, -5.0, geo.xmax, geo.ymax, +5.0);

    for (const auto& e : geo.electrodes) {
      if (e.kind == Detector::ElectrodeKind::WireRow &&
          ToLower(e.name).find("anode") != std::string::npos) {
        sensor.AddElectrode(comp.get(), e.name);
      }
    }

    sensor.SetTransferFunction(TransferFunction);

    const double gap_val   = DetectGap(geo);
    const double t0_sig_ns = 0.0;
    const double dt_sig_ns = 0.5;
    const double tmax_sig  = (th.t_center_ns.empty() ? 200.0 : th.t_center_ns.back()) + 150.0;
    const int nSigBins = std::max(200, (int)std::ceil((tmax_sig - t0_sig_ns) / dt_sig_ns));
    sensor.SetTimeWindow(t0_sig_ns, dt_sig_ns, nSigBins);

    TRandom3 gen(0);
    const int nEvents = 0;

    std::printf("[map] Generating high-resolution field maps & lines...\n");
    const int NX = 1000, NY = 1000;

    auto hE = std::make_unique<TH2D>("|E|", "Electric Field |E|;X [cm];Y [cm]",
                                     NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    hE->SetDirectory(nullptr);
    hE->SetStats(0);

    auto hV = std::make_unique<TH2D>("V", "Electric Potential V;X [cm];Y [cm]",
                                     NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    hV->SetDirectory(nullptr);
    hV->SetStats(0);

    double minPos = 1e98, maxVal_global = 0.0;
    for (int ix = 1; ix <= NX; ++ix) {
      const double x = hE->GetXaxis()->GetBinCenter(ix);
      for (int iy = 1; iy <= NY; ++iy) {
        const double y = hE->GetYaxis()->GetBinCenter(iy);
        double ex = 0, ey = 0, ez = 0, V = 0;
        Medium* m = nullptr;
        int status = 0;
        comp->ElectricField(x, y, 0, ex, ey, ez, V, m, status);
        if (status == 0) {
          const double E = std::hypot(ex, std::hypot(ey, ez));
          hE->SetBinContent(ix, iy, E);
          hV->SetBinContent(ix, iy, V);
          if (E > 0 && E < minPos) minPos = E;
          if (E > maxVal_global) maxVal_global = E;
        }
      }
    }

    const double plot_xmin = geo.xmin;
    const double plot_xmax = geo.xmax;
    const double plot_ymin = -gap_val * 1.05;
    const double plot_ymax =  gap_val * 1.05;

    TLatex latexConfig;
    latexConfig.SetNDC();
    latexConfig.SetTextSize(0.035);
    latexConfig.SetTextColor(kBlack);
    latexConfig.SetTextAlign(13);

    auto cE = std::make_unique<TCanvas>("cE", "|E| Map", 1000, 780);
    cE->SetRightMargin(0.15);
    cE->SetGrid();
    cE->SetLogz(1);

    if (minPos < 1e97) hE->SetMinimum(minPos * 0.9);
    hE->SetMaximum(maxVal_global * 1.05);

    hE->GetXaxis()->SetRangeUser(plot_xmin, plot_xmax);
    hE->GetYaxis()->SetRangeUser(plot_ymin, plot_ymax);
    hE->Draw("COLZ");

    gPad->SetFixedAspectRatio();
    DrawElectrodes(geo);

    latexConfig.DrawLatex(0.15, 0.85, textVolt);
    latexConfig.DrawLatex(0.15, 0.81, textDiam);

    const std::string emap_png = (fs::path(figDir) / "emap_absE_highres.png").string();
    cE->SaveAs(emap_png.c_str());
    cE->SaveAs((emapsSummaryDir / Form("emap_absE_%s.png", conf_name)).string().c_str());

    auto cCombo = std::make_unique<TCanvas>("cComboMap", "Combo Map", 1000, 780);
    cCombo->SetGrid();

    ViewField vfFL;
    vfFL.SetCanvas(cCombo.get());
    vfFL.SetComponent(comp.get());
    vfFL.SetArea(plot_xmin, plot_ymin, plot_xmax, plot_ymax);

    auto inMedium = [&](double x, double y) {
      double ex = 0, ey = 0, ez = 0, V = 0;
      Medium* m = nullptr;
      int status = 0;
      comp->ElectricField(x, y, 0.0, ex, ey, ez, V, m, status);
      return (m != nullptr && status == 0);
    };

    auto make_line_seeds_symmetric = [&](double yseed,
                                         std::vector<double>& xs,
                                         std::vector<double>& ys,
                                         std::vector<double>& zs) {
      const double pitch_wire = 0.20;
      const int nLinesPerPitch = 10;
      const double eps = 1.0e-4;
      const double x_center = 0.0;

      const int nCellL = int(std::ceil((x_center - plot_xmin) / pitch_wire)) + 1;
      const int nCellR = int(std::ceil((plot_xmax - x_center) / pitch_wire)) + 1;

      for (int n = -nCellL; n <= nCellR; ++n) {
        const double x_wire = x_center + n * pitch_wire;

        for (int j = 1; j <= nLinesPerPitch / 2; ++j) {
          const double u = double(j) / double(nLinesPerPitch / 2 + 1);
          const double offset = 0.5 * pitch_wire * u;

          const double xL = x_wire - offset + eps;
          const double xR = x_wire + offset + eps;

          if (xL >= plot_xmin && xL <= plot_xmax && inMedium(xL, yseed)) {
            xs.push_back(xL); ys.push_back(yseed); zs.push_back(0.0);
          }
          if (xR >= plot_xmin && xR <= plot_xmax && inMedium(xR, yseed)) {
            xs.push_back(xR); ys.push_back(yseed); zs.push_back(0.0);
          }
        }
      }
    };

    std::vector<double> xs2, ys2, zs2;
    make_line_seeds_symmetric(+0.25, xs2, ys2, zs2);
    make_line_seeds_symmetric(-0.25, xs2, ys2, zs2);

    if (!xs2.empty()) vfFL.PlotFieldLines(xs2, ys2, zs2, true, true);

    cCombo->SetRightMargin(0.15);
    cCombo->SetLogz(1);
    hE->Draw("COLZ SAME");
    gPad->SetFixedAspectRatio();

    hV->SetContour(40);
    hV->SetLineColor(kGray + 2);
    hV->SetLineWidth(1);
    hV->Draw("CONT3 SAME");

    DrawElectrodes(geo);

    TIter next(cCombo->GetListOfPrimitives());
    TObject* obj;
    std::vector<TObject*> linesToFront;
    while ((obj = next())) {
      if (obj->InheritsFrom("TGraph") || obj->InheritsFrom("TPolyLine")) {
        linesToFront.push_back(obj);
      }
    }
    for (auto* l : linesToFront) l->Draw("L");

    latexConfig.DrawLatex(0.15, 0.85, textVolt);
    latexConfig.DrawLatex(0.15, 0.81, textDiam);

    const std::string combo_png = (fs::path(figDir) / "fieldmap_combo_highres.png").string();
    cCombo->SaveAs(combo_png.c_str());
    cCombo->SaveAs((emapsSummaryDir / Form("fieldmap_combo_%s.png", conf_name)).string().c_str());

    // =========================================================
    // Ex boundary profile at x = 1 mm
    // =========================================================
    const double xBoundaryEval = 0.10;
    const int nExBins = 600;

    auto hExBoundary = std::make_unique<TH1D>(
        "hExBoundary",
        Form("E_{x} at x = %.3f cm;Y [cm];E_{x} [V/cm]", xBoundaryEval),
        nExBins, plot_ymin, plot_ymax);
    hExBoundary->SetDirectory(nullptr);
    hExBoundary->SetLineWidth(2);

    ExBoundarySummary exSumm = EvaluateBoundaryEx(comp.get(), xBoundaryEval, hExBoundary.get());

    auto cEx = std::make_unique<TCanvas>("cExBoundary", "Ex boundary", 900, 700);
    cEx->SetGrid();
    hExBoundary->Draw("HIST");

    auto l0 = new TLine(plot_ymin, 0.0, plot_ymax, 0.0);
    l0->SetBit(kCanDelete);
    l0->SetLineStyle(2);
    l0->SetLineColor(kRed + 1);
    l0->Draw("same");

    TLatex latexEx;
    latexEx.SetNDC();
    latexEx.SetTextSize(0.032);
    latexEx.SetTextAlign(13);
    latexEx.DrawLatex(0.15, 0.88, textVolt);
    latexEx.DrawLatex(0.15, 0.84, textDiam);
    latexEx.DrawLatex(0.15, 0.80, Form("x = %.3f cm", xBoundaryEval));
    latexEx.DrawLatex(0.15, 0.76, Form("<Ex^{2}> = %.6e", exSumm.ex2_mean));
    latexEx.DrawLatex(0.15, 0.72, Form("max|Ex| = %.6e", exSumm.ex_abs_max));

    const std::string ex_profile_png = (exDir / "Ex_vs_Y_x2mm.png").string();
    cEx->SaveAs(ex_profile_png.c_str());
    cEx->SaveAs((exProfileSummaryDir / Form("Ex_vs_Y_x2mm_%s.png", conf_name)).string().c_str());

    {
      ExProfileForOverlay p;
      p.label = conf_name;
      p.vPW = vPW;
      p.y.reserve(hExBoundary->GetNbinsX());
      p.ex.reserve(hExBoundary->GetNbinsX());
      for (int ib = 1; ib <= hExBoundary->GetNbinsX(); ++ib) {
        p.y.push_back(hExBoundary->GetXaxis()->GetBinCenter(ib));
        p.ex.push_back(hExBoundary->GetBinContent(ib));
      }
      exProfiles.push_back(std::move(p));
    }

    // =========================================================
    // Field score evaluation (A0/A1/PW_-1/PW_0)
    // =========================================================
    std::printf("[eval] Evaluating central field scores...\n");
    const FieldScoreSummary fieldSumm = EvaluateAndSaveFieldScores(
        comp.get(), geo, sensePitchX, fieldScoreDir, textVolt, textDiam, conf_name);

    const std::string field_score_png = (fieldScoreDir / "field_score_maps.png").string();
    const std::string field_score_root = (fieldScoreDir / "field_score.root").string();
    const std::string field_score_csv = (fieldScoreDir / "field_score_summary.csv").string();
    if (fs::exists(field_score_png)) {
      fs::copy_file(field_score_png,
                    fieldScoreSummaryDir / (std::string(conf_name) + "_field_score_maps.png"),
                    fs::copy_options::overwrite_existing);
    }

    // =========================================================
    // field-line same-seed drift-only save
    // =========================================================
    std::printf("[eval] Saving drift info for field-line seeds...\n");

    DriftLineRKF driftFieldLine;
    driftFieldLine.SetSensor(&sensor);

    for (size_t i = 0; i < xs2.size(); ++i) {
      const double x0 = xs2[i];
      const double y0 = ys2[i];
      const double z0 = zs2[i];

      double ex = 0.0, ey = 0.0, ez = 0.0, V = 0.0;
      Medium* med = nullptr;
      int estatus = 0;
      comp->ElectricField(x0, y0, z0, ex, ey, ez, V, med, estatus);

      fld_seed_id = static_cast<int>(i);
      fld_x0_cm = x0;
      fld_y0_cm = y0;
      fld_z0_cm = z0;
      fld_x1_cm = -999.;
      fld_y1_cm = -999.;
      fld_z1_cm = -999.;
      fld_t_drift_ns = -999.;
      fld_drift_path_cm = -999.;
      fld_wire_x_cm = -999.;
      fld_wire_y_cm = -999.;
      fld_wire_type = -999;
      std::snprintf(fld_wire_name, sizeof(fld_wire_name), "%s", "NO_WIRE");
      std::snprintf(fld_wire_uid,  sizeof(fld_wire_uid),  "%s", "NO_WIRE");
      fld_status = estatus;

      if (estatus == 0 && med != nullptr) {
        driftFieldLine.DriftElectron(x0, y0, z0, 0.0);

        double x1 = 0.0, y1 = 0.0, z1 = 0.0, t1 = 0.0;
        int status = 0;
        driftFieldLine.GetEndPoint(x1, y1, z1, t1, status);

        fld_x1_cm = x1;
        fld_y1_cm = y1;
        fld_z1_cm = z1;
        fld_t_drift_ns = t1;
        fld_drift_path_cm = driftFieldLine.GetPathLength();
        fld_status = status;

        auto end_wire = NearestWireSurface(geo, x1, y1);
        if (end_wire.hasWire) {
          fld_wire_x_cm = end_wire.xw_cm;
          fld_wire_y_cm = end_wire.yw_cm;
          fld_wire_type = end_wire.wire_type;
          std::snprintf(fld_wire_name, sizeof(fld_wire_name), "%s", end_wire.wire_name.c_str());
          std::snprintf(fld_wire_uid,  sizeof(fld_wire_uid),  "%s", end_wire.wire_uid.c_str());
        }
      }

      fEvalOut->cd();
      tFld->Fill();
    }

    // =========================================================
    // L-t evaluation (2D seed, drift-only)
    // =========================================================
    std::printf("[eval] Evaluating L-t relation (drift-only)...\n");

    std::vector<LtPoint> ltPoints;
    std::vector<LtPoint> ltPointsForPlots;
    DriftLineRKF driftLt;
    driftLt.SetSensor(&sensor);

    const int nSeedY = 41;
    const int nSeedPerCellX = 9;
    const double yMinSeed = 0.02;
    const double yMaxSeed = 0.48;

    const std::vector<double> seedXs = BuildAnodeCenteredSeedX(geo, seedMarginCells, nSeedPerCellX);
    std::printf("[seed] anode-centered x seeds = %zu columns (nPerCell=%d)\n",
                seedXs.size(), nSeedPerCellX);

    for (int iy = 0; iy < nSeedY; ++iy) {
      const double fy = (double(iy) + 0.5) / double(nSeedY);
      const double y0 = yMinSeed + fy * (yMaxSeed - yMinSeed);

      for (size_t ix = 0; ix < seedXs.size(); ++ix) {
        const double x0 = seedXs[ix];

        double ex = 0.0, ey = 0.0, ez = 0.0, V = 0.0;
        Medium* med = nullptr;
        int estatus = 0;
        comp->ElectricField(x0, y0, 0.0, ex, ey, ez, V, med, estatus);

        LtPoint p;
        p.x0_cm = x0;
        p.y0_cm = y0;
        p.x1_cm = -999.;
        p.y1_cm = -999.;
        p.t_ns = -999.;
        p.path_cm = -999.;
        p.dx_end_cm = -999.;
        p.L_cm = -999.;
        p.wire_x_cm = -999.;
        p.wire_y_cm = -999.;
        p.wire_name = "NO_WIRE";
        p.wire_uid = "NO_WIRE";
        p.target_wire_name = "NO_WIRE";
        p.target_wire_uid = "NO_WIRE";
        p.wire_type = -999;
        p.is_target = 0;
        p.status = estatus;

        const auto target_wire = FindTargetSignalWireByX(geo, x0);
        if (target_wire.hasWire) {
          p.target_wire_name = target_wire.wire_name;
          p.target_wire_uid  = target_wire.wire_uid;
        }

        if (estatus == 0 && med != nullptr) {
          driftLt.DriftElectron(x0, y0, 0.0, 0.0);

          double x1 = 0.0, y1 = 0.0, z1 = 0.0, t1 = 0.0;
          int status = 0;
          driftLt.GetEndPoint(x1, y1, z1, t1, status);

          p.status = status;
          p.x1_cm = x1;
          p.y1_cm = y1;
          p.t_ns = t1;
          p.path_cm = driftLt.GetPathLength();
          p.dx_end_cm = x1 - x0;

          auto end_wire = NearestWireSurface(geo, x1, y1);
          if (end_wire.hasWire) {
            p.L_cm = std::abs(y0 - end_wire.yw_cm);
            p.wire_x_cm = end_wire.xw_cm;
            p.wire_y_cm = end_wire.yw_cm;
            p.wire_name = end_wire.wire_name;
            p.wire_uid  = end_wire.wire_uid;
            p.wire_type = end_wire.wire_type;
            p.is_target = (p.wire_uid == p.target_wire_uid) ? 1 : 0;
          }
        }

        lt_t_ns = p.t_ns;
        lt_L_cm = p.L_cm;
        lt_x0_cm = p.x0_cm;
        lt_y0_cm = p.y0_cm;
        lt_x1_cm = p.x1_cm;
        lt_y1_cm = p.y1_cm;
        lt_path_cm = p.path_cm;
        lt_dx_end_cm = p.dx_end_cm;
        lt_wire_x_cm = p.wire_x_cm;
        lt_wire_y_cm = p.wire_y_cm;
        lt_wire_type = p.wire_type;
        lt_is_target = p.is_target;
        lt_status = p.status;
        std::snprintf(lt_wire_name, sizeof(lt_wire_name), "%s", p.wire_name.c_str());
        std::snprintf(lt_wire_uid, sizeof(lt_wire_uid), "%s", p.wire_uid.c_str());
        std::snprintf(lt_target_wire_name, sizeof(lt_target_wire_name), "%s", p.target_wire_name.c_str());
        std::snprintf(lt_target_wire_uid, sizeof(lt_target_wire_uid), "%s", p.target_wire_uid.c_str());

        fEvalOut->cd();
        tLt->Fill();

        if (p.L_cm >= 0.0 && p.t_ns >= 0.0 && p.wire_uid != "NO_WIRE") {
          ltPoints.push_back(p);
          ltPointsForPlots.push_back(p);
        }
      }
    }

    LtFitSummary summ = AnalyzeLt(ltPoints);

    s_vcat = vCat;
    s_vanode = vAnode;
    s_vpw = vPW;
    s_npoints = summ.n;
    s_slope = summ.slope;
    s_intercept = summ.intercept;
    s_r2 = summ.r2;
    s_residual_rms = summ.residual_rms;
    s_monotonic_violations = summ.monotonic_violations;
    s_mean_dx_end = summ.mean_dx_end;
    s_mean_path_excess = summ.mean_path_excess;
    s_collect_eff = summ.collect_eff;
    fEvalOut->cd();
    tLtSummary->Fill();

    std::printf("[eval] n=%d, R2=%.6f, residual RMS=%.6f cm, monoViol=%d, collectEff=%.4f\n",
                summ.n, summ.r2, summ.residual_rms, summ.monotonic_violations, summ.collect_eff);

    SaveWireLtPlots(ltPlotDir, ltPointsForPlots, textVolt, textDiam);

    const fs::path allWireLtPng = ltPlotDir / "AllWire_Lt.png";
    SaveAllWireLtPlot(allWireLtPng, ltPointsForPlots, textVolt, textDiam, "All-wire L-t correlation");
    SaveAllWireLtPlot(allWireLtSummaryDir / (std::string(conf_name) + "_AllWire_Lt.png"),
                      ltPointsForPlots, textVolt, textDiam,
                      Form("All-wire L-t correlation (%s)", conf_name));

    TLatex latexEval;
    latexEval.SetNDC();
    latexEval.SetTextSize(0.032);
    latexEval.SetTextColor(kBlack);
    latexEval.SetTextAlign(13);

    TString textRms = TString::Format("L-t residual RMS = %.5f cm", summ.residual_rms);
    TString textR2  = TString::Format("L-t R^{2} = %.5f", summ.r2);
    TString textEff = TString::Format("same-wire frac = %.5f", summ.collect_eff);
    TString textField = TString::Format("S_{field,raw} = %.6e", fieldSumm.sFieldRaw);

    cCombo->cd();
    latexEval.DrawLatex(0.15, 0.77, textRms);
    latexEval.DrawLatex(0.15, 0.73, textR2);
    latexEval.DrawLatex(0.15, 0.69, textEff);
    if (fieldSumm.valid) latexEval.DrawLatex(0.15, 0.65, textField);
    cCombo->SaveAs(combo_png.c_str());
    cCombo->SaveAs((emapsSummaryDir / Form("fieldmap_combo_%s.png", conf_name)).string().c_str());

    cE->cd();
    latexEval.DrawLatex(0.15, 0.77, textRms);
    latexEval.DrawLatex(0.15, 0.73, textR2);
    latexEval.DrawLatex(0.15, 0.69, textEff);
    if (fieldSumm.valid) latexEval.DrawLatex(0.15, 0.65, textField);
    cE->SaveAs(emap_png.c_str());
    cE->SaveAs((emapsSummaryDir / Form("emap_absE_%s.png", conf_name)).string().c_str());

    {
      ConfigRankingEntry e;
      e.conf_name = conf_name;
      e.geom_tag = geom_tag;
      e.vCat = vCat;
      e.vAnode = vAnode;
      e.vPW = vPW;
      e.dCat_um = dCat * 1e4;
      e.dAnode_um = dAnode * 1e4;
      e.dPW_um = dPW * 1e4;
      e.n_points = summ.n;
      e.residual_rms = summ.residual_rms;
      e.r2 = summ.r2;
      e.monotonic_violations = summ.monotonic_violations;
      e.mean_dx_end = summ.mean_dx_end;
      e.mean_path_excess = summ.mean_path_excess;
      e.collect_eff = summ.collect_eff;
      e.ex_boundary_x_cm = xBoundaryEval;
      e.ex2_mean = exSumm.ex2_mean;
      e.ex2_int = exSumm.ex2_int;
      e.ex_abs_max = exSumm.ex_abs_max;
      e.ex_mean_abs = exSumm.ex_mean_abs;
      e.sAtrans = fieldSumm.valid ? fieldSumm.sAtrans : -1.0;
      e.sAalign = fieldSumm.valid ? fieldSumm.sAalign : -1.0;
      e.sPWmirror = fieldSumm.valid ? fieldSumm.sPWmirror : -1.0;
      e.sPWalign = fieldSumm.valid ? fieldSumm.sPWalign : -1.0;
      e.sELine = fieldSumm.valid ? fieldSumm.sELine : -1.0;
      e.sFieldRaw = fieldSumm.valid ? fieldSumm.sFieldRaw : -1.0;
      e.nAtrans = fieldSumm.valid ? fieldSumm.nAtrans : 0;
      e.nAalign = fieldSumm.valid ? fieldSumm.nAalign : 0;
      e.nPWmirror = fieldSumm.valid ? fieldSumm.nPWmirror : 0;
      e.nPWalign = fieldSumm.valid ? fieldSumm.nPWalign : 0;
      e.combo_png = combo_png;
      e.emap_png = emap_png;
      e.ex_profile_png = ex_profile_png;
      e.allwire_lt_png = allWireLtPng.string();
      e.field_score_png = field_score_png;
      e.field_score_root = field_score_root;
      e.field_score_csv = field_score_csv;
      rankingEntries.push_back(e);
    }

    for (int ev = 0; ev < nEvents; ++ev) {
      std::printf("  >>> Processing Event %d / %d <<<\n", ev + 1, nEvents);

      const double true_angle_deg = -30.0;
      const double a_true = std::tan(true_angle_deg * M_PI / 180.0);
      const double b_true = 0.25 + gen.Uniform(-0.05, 0.05);

      const double x_min_trk = geo.xmin;
      const double x_max_trk = geo.xmax;
      const int nClusters = 879;

      AvalancheMC amcTrack;
      amcTrack.SetSensor(&sensor);
      TryEnableSignalCalculation(amcTrack);
      amcTrack.UseWeightingPotential(true);

      struct ClusterInfo { double x, y; };
      std::vector<ClusterInfo> clusters;

      for (int i = 0; i < nClusters; ++i) {
        const double frac = gen.Uniform(0.0, 1.0);
        const double x = x_min_trk + frac * (x_max_trk - x_min_trk);
        const double y = a_true * x + b_true;
        if (y >= geo.ymin && y <= geo.ymax) clusters.push_back({x, y});
      }

      std::map<std::string, GainPlotInfo> globalGainMap;
      sensor.ClearSignal();

      for (size_t i = 0; i < clusters.size(); ++i) {
        auto& cl = clusters[i];
        amcTrack.AvalancheElectron(cl.x, cl.y, 0, 0);

        const int n_e = amcTrack.GetNumberOfElectronEndpoints();
        if (n_e > 0) {
          double x0, y0, z0, t0, x1, y1, z1, t1;
          int status;
          amcTrack.GetElectronEndpoint(0, x0, y0, z0, t0, x1, y1, z1, t1, status);

          auto end_wire = NearestWireSurface(geo, x1, y1);
          if (end_wire.hasWire) {
            std::string wUid = end_wire.wire_uid;
            if (globalGainMap.find(wUid) == globalGainMap.end()) {
              globalGainMap[wUid] = {
                end_wire.wire_type, end_wire.xw_cm, end_wire.yw_cm, 0.0, 0, 0,
                end_wire.wire_name, end_wire.wire_uid
              };
            }
            globalGainMap[wUid].n_primary += 1;
            globalGainMap[wUid].n_total   += n_e;

            d_event_id = ev;
            d_x_start  = x0;
            d_y_start  = y0;
            d_x_end    = x1;
            d_y_end    = y1;
            d_t_drift  = t1 - t0;
            d_dist     = std::hypot(x1 - x0, y1 - y0);
            std::snprintf(d_wire_name, sizeof(d_wire_name), "%s", end_wire.wire_name.c_str());
            std::snprintf(d_wire_uid,  sizeof(d_wire_uid),  "%s", end_wire.wire_uid.c_str());
            fDriftOut->cd();
            tDrift->Fill();
          }
        }
      }

      sensor.ConvoluteSignals();

      for (auto& [wUid, data] : globalGainMap) {
        data.gain = (data.n_primary > 0) ? (double)data.n_total / data.n_primary : 0.0;

        g_event_id  = ev;
        g_wire_type = data.wType;
        g_wire_tag  = data.xw;
        g_x_wire    = data.xw;
        g_y_wire    = data.yw;
        g_n_primary = data.n_primary;
        g_n_total   = data.n_total;
        g_gain      = data.gain;
        std::snprintf(g_wire_name, sizeof(g_wire_name), "%s", data.wire_name.c_str());
        std::snprintf(g_wire_uid,  sizeof(g_wire_uid),  "%s", data.wire_uid.c_str());

        fGainOut->cd();
        tGain->Fill();
      }

      struct HitCand { double xwire, ywire, y_up, y_dn; };
      std::vector<HitCand> hits;
      const double signalThreshold = 4.0;

      for (const auto& e : geo.electrodes) {
        if (e.kind == Detector::ElectrodeKind::WireRow &&
            ToLower(e.name).find("anode") != std::string::npos) {
          std::string wName = e.name;
          double xw = e.x0;
          double yw = e.y;

#ifdef DEBUG_WAVEFORM
          if (ev == 0) {
            auto cWave = std::make_unique<TCanvas>("cWave", "", 600, 400);
            cWave->SetGrid();

            std::vector<double> tx(nSigBins), vy(nSigBins);
            for (int k = 0; k < nSigBins; ++k) {
              tx[k] = t0_sig_ns + k * dt_sig_ns;
              vy[k] = sensor.GetSignal(wName, k, true);
            }

            auto gW = std::make_unique<TGraph>(nSigBins, tx.data(), vy.data());
            gW->SetTitle((wName + " (Ev " + std::to_string(ev) + ");Time [ns];Signal [mV]").c_str());
            gW->SetLineWidth(1);
            gW->SetLineColor(kBlue);
            gW->Draw("AL");

            auto lth_p = new TLine(t0_sig_ns, signalThreshold, tmax_sig, signalThreshold);
            lth_p->SetBit(kCanDelete);
            lth_p->SetLineColor(kRed);
            lth_p->SetLineStyle(2);
            lth_p->Draw("same");

            auto lth_n = new TLine(t0_sig_ns, -signalThreshold, tmax_sig, -signalThreshold);
            lth_n->SetBit(kCanDelete);
            lth_n->SetLineColor(kRed);
            lth_n->SetLineStyle(2);
            lth_n->Draw("same");

            char wname_png[128];
            std::snprintf(wname_png, sizeof(wname_png),
                          "%s/waveform/wire_sig_ev%d_%s.png",
                          figDir.c_str(), ev, wName.c_str());
            cWave->SaveAs(wname_png);
          }
#endif

          double t_hit = -1.0;
          for (int k = 0; k < nSigBins; ++k) {
            if (std::abs(sensor.GetSignal(wName, k, true)) > signalThreshold) {
              t_hit = t0_sig_ns + k * dt_sig_ns;
              break;
            }
          }

          if (t_hit >= 0.0) {
            double L_meas = CalculateL0Directly(th, gap_val, t_hit);
            if (L_meas < 0) L_meas = 0;
            hits.push_back({xw, yw, yw + L_meas, yw - L_meas});
          }
        }
      }

      if (hits.size() < 3) continue;

      const int nHits = (int)hits.size();
      const unsigned long long nComb = 1ULL << nHits;
      double best_a = 0.0, best_b = 0.0, min_chi2 = 1.0e300;
      unsigned long long best_mask = 0;

      for (unsigned long long mask = 0; mask < nComb; ++mask) {
        double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
        for (int i = 0; i < nHits; ++i) {
          const double x = hits[i].xwire;
          const double y = ((mask >> i) & 1) ? hits[i].y_dn : hits[i].y_up;
          Sx += x;
          Sy += y;
          Sxx += x * x;
          Sxy += x * y;
        }

        const double delta = nHits * Sxx - Sx * Sx;
        if (std::abs(delta) < 1e-12) continue;

        const double a_tmp = (nHits * Sxy - Sx * Sy) / delta;
        const double b_tmp = (Sxx * Sy - Sx * Sxy) / delta;

        double chi2 = 0.0;
        for (int i = 0; i < nHits; ++i) {
          const double x = hits[i].xwire;
          const double y = ((mask >> i) & 1) ? hits[i].y_dn : hits[i].y_up;
          const double diff = y - (a_tmp * x + b_tmp);
          chi2 += diff * diff;
        }

        if (chi2 < min_chi2) {
          min_chi2 = chi2;
          best_a = a_tmp;
          best_b = b_tmp;
          best_mask = mask;
        }
      }

      (void)best_mask;

      auto cTrk = std::make_unique<TCanvas>("cTrk", "track reco", 1000, 780);
      cTrk->SetGrid();

      auto frame = std::make_unique<TH2D>("frame", "Reconstructed Track;X;Y",
                                          10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
      frame->SetDirectory(nullptr);
      frame->SetStats(0);
      frame->Draw();

      gPad->SetFixedAspectRatio();
      DrawElectrodes(geo);

      auto fTrue = new TF1("fTrue", "[0]*x+[1]", geo.xmin, geo.xmax);
      fTrue->SetBit(kCanDelete);
      fTrue->SetParameters(a_true, b_true);
      fTrue->SetLineColor(kBlue);
      fTrue->SetLineStyle(2);
      fTrue->Draw("same");

      auto fReco = new TF1("fReco", "[0]*x+[1]", geo.xmin, geo.xmax);
      fReco->SetBit(kCanDelete);
      fReco->SetParameters(best_a, best_b);
      fReco->SetLineColor(kRed);
      fReco->SetLineWidth(2);
      fReco->Draw("same");

      for (const auto& h : hits) {
        const double L = std::abs(h.y_up - h.ywire);
        auto el = new TEllipse(h.xwire, h.ywire, L, L);
        el->SetBit(kCanDelete);
        el->SetFillStyle(0);
        el->SetLineColor(kGreen + 2);
        el->Draw("same");
      }

      latexConfig.DrawLatex(0.15, 0.85, textVolt);
      latexConfig.DrawLatex(0.15, 0.81, textDiam);

      cTrk->SaveAs(Form("%s/track_reco_ev%d.png", figDir.c_str(), ev));

      t_event_id = ev;
      t_a_true   = a_true;
      t_b_true   = b_true;
      t_a_reco   = best_a;
      t_b_reco   = best_b;
      t_chi2     = min_chi2;
      fOut->cd();
      tree->Fill();
    }

    std::printf("[debug] ltTree entries = %lld\n", tLt->GetEntries());
    std::printf("[debug] fieldLineDriftTree entries = %lld\n", tFld->GetEntries());
    std::printf("[debug] ltSummary entries = %lld\n", tLtSummary->GetEntries());

    fOut->cd();
    tree->Write();
    fOut->Close();

    fDriftOut->cd();
    tDrift->Write();
    fDriftOut->Close();

    fGainOut->cd();
    tGain->Write();
    fGainOut->Close();

    fVerifyOut->cd();
    tVerify->Write();
    fVerifyOut->Close();

    fEvalOut->cd();
    tLt->Write();
    tLtSummary->Write();
    tFld->Write();
    fEvalOut->Write();
    fEvalOut->Close();

    std::printf("[debug] lt_eval.root written and closed\n");

  }}}}}}

  ComputeNormalizedFieldScores(rankingEntries);
  RewritePerConfigFieldScoreSummaries(rankingEntries);
  SaveExOverlayPlots(exProfileSummaryDir, exProfiles, "Ex vs Y at x = 2 mm");
  SaveAllWireLtSummaryPlots(allWireLtSummaryDir, rankingEntries, 10);

  fs::path rankingMdPath = sweepBaseDir / "lt_ranking_summary.md";
  WriteRankingMarkdown(rankingMdPath, rankingEntries, 10);
  WriteFieldScoreTop10Markdown(totalScoreDir / "field_score_top10.md", rankingEntries, 10);
  WriteFieldScoreTop10Csv(totalScoreDir / "field_score_top10.csv", rankingEntries, 10);

  SaveRmsVsPotentialPlot(sweepBaseDir, rankingEntries, "L-t residual RMS vs potential-wire voltage");
  SaveEx2VsPotentialPlot(sweepBaseDir, rankingEntries,
                         "Mean Ex^{2} vs potential-wire voltage");
  SaveFieldScoreVsPotentialPlot(totalScoreDir, rankingEntries,
                                "Central field scores vs potential-wire voltage");

  std::printf("\n[Info] All sweep configurations completed. (%d configurations for %s)\n",
              config_count, geom_tag.c_str());
  std::printf("[Info] All summary maps are saved in: %s\n", emapsSummaryDir.string().c_str());
  std::printf("[Info] Ex profile summary dir: %s\n", exProfileSummaryDir.string().c_str());
  std::printf("[Info] All-wire Lt summary dir: %s\n", allWireLtSummaryDir.string().c_str());
  std::printf("[Info] Field-score summary dir: %s\n", fieldScoreSummaryDir.string().c_str());
  std::printf("[Info] Total-score dir: %s\n", totalScoreDir.string().c_str());
  std::printf("[Info] Field-score Top 10 saved to: %s\n", (totalScoreDir / "field_score_top10.md").string().c_str());
  std::printf("[Info] Ranking markdown saved to: %s\n", rankingMdPath.string().c_str());
  return 0;
}
