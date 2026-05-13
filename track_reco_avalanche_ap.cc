/**
 * track_reco_combinatorial.cc
 *
 * [完全版: Finiteジオメトリ / 全アノード波形 / 高解像度ズームマップ / PW対応 / 画面キープ]
 * - 周期境界を廃止した「個別配置ジオメトリ」に対応。
 * - 現実のDAQと同様に全チャンネルの波形を独立して取得し閾値(4mV)判定。
 * - 1 fC = 1 mV のコンバージョンゲインを物理定数から計算し適用。
 * - アノード(赤)とポテンシャルワイヤー(青)の電場プロファイルを独立して計算・描画。
 * (電場が10,000 V/cmを超えるY座標の範囲を自動検出しグラフに描画)
 * - Shockley-Ramo定理の検算用ツリー(verifyTree)を追加。生電流(電子・イオン)の積分値とWP差分を記録。
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

// ROOT includes
#include "TApplication.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
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
#include "TMultiGraph.h" 

// Garfield++ includes
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewDrift.hh" 

// Local includes
#include "detector_geometry.hh"
#include "geometry/eval_uniform.hh"
#include "eval_maps.hh"

// --- ★ Debug Macro Switch ★ ------------------------------------------------
#define DEBUG_WAVEFORM 

// --- ROOT output switch ----------------------------------------------------
#define SAVE_DRIFT_ROOT   0
#define SAVE_GAIN_ROOT    0
#define SAVE_VERIFY_ROOT  0
// メインROOT (track_results_job*.root) は残す
// ----------------------------------------------------------------------------

// --- geometry choice --------------------------------------------------------
#if defined(GEOM_PLATE)
  #include "geometry/geom_ap_plate.hh"
#elif defined(GEOM_WIRE)
  #include "geometry/geom_ap_wirecath.hh"
#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

using namespace Garfield;
namespace fs = std::filesystem;

// ==================== Preamplifier Response ====================
double TransferFunction(double t) {
  constexpr double tau = 20.0; // プリアンプの時定数 [ns]
  constexpr double e_charge = 1.602176634e-19; // 素電荷 [C]
  constexpr double fC_in_C  = 1.0e-15;         // 1 fC のクーロン単位 [C]
  constexpr double electrons_per_fC = fC_in_C / e_charge; 
  constexpr double target_gain_mV_per_fC = 1.0; // コンバージョンゲイン [mV/fC]
  constexpr double gain = target_gain_mV_per_fC; // そのまま 1.0 にする
  
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
  double d_surface_cm;
  double r_cm;          
  double xw_cm, yw_cm;  
  bool   hasWire;       
  int    wire_type; // 0: Anode, 1: PW
  std::string wire_name; 
};

static void DrawElectrodes(const Detector::Geometry& g) {
  const double xmin = g.xmin, xmax = g.xmax;
  for (const auto& e : g.electrodes) {
    if (e.kind == Detector::ElectrodeKind::PlaneY) {
      auto ln = new TLine(xmin, e.y, xmax, e.y);
      ln->SetLineColor(kGray + 3); ln->SetLineStyle(2); ln->SetLineWidth(3);
      ln->Draw("same");
      continue;
    }

    auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
    const std::string nm = lower(e.name);
    const bool isAnode = nm.find("anode") != std::string::npos;
    const bool isPW    = nm.find("pw")    != std::string::npos;

    if (!g.periodicX) {
      const double x = e.x0; if (x < xmin || x > xmax) continue;
      auto cir = new TEllipse(x, e.y, e.radius, e.radius); 
      cir->SetFillStyle(1001);
      if (isAnode) cir->SetFillColor(kOrange + 7);
      else if (isPW) cir->SetFillColor(kAzure  + 2);
      else cir->SetFillColor(kBlack);
      cir->SetLineColor(kBlack); cir->SetLineWidth(1); cir->Draw("same");
      continue;
    }

    const double pitch = g.pitchX;
    int nL = int(std::floor((xmin - e.x0) / pitch)) - 1;
    int nR = int(std::ceil ((xmax - e.x0) / pitch)) + 1;
    for (int n = nL; n <= nR; ++n) {
      const double x = e.x0 + n * pitch; if (x < xmin || x > xmax) continue;
      auto cir = new TEllipse(x, e.y, e.radius, e.radius);
      cir->SetFillStyle(1001);
      if (isAnode) cir->SetFillColor(kOrange + 7);
      else if (isPW) cir->SetFillColor(kAzure  + 2);
      else cir->SetFillColor(kBlack);
      cir->SetLineColor(kBlack); cir->SetLineWidth(1); cir->Draw("same");
    }
  }
}

static double DetectGap(const Detector::Geometry& g) {
  double gap = 0.0;
  for (const auto& e : g.electrodes) gap = std::max(gap, std::abs(e.y));
  if (gap <= 0.0) gap = std::min(std::abs(g.ymax), std::abs(g.ymin));
  return gap;
}

static WireNearest NearestWireSurface(const Detector::Geometry& g, double x, double y){
  WireNearest out; out.d_surface_cm = 1e300; out.r_cm = 0.0; out.xw_cm = 0.0; out.yw_cm = 0.0; out.hasWire=false; out.wire_type = -1;
  auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };

  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    std::string nm = lower(e.name);
    bool isAnode = nm.find("anode") != std::string::npos;
    bool isPW    = nm.find("pw")    != std::string::npos;
    if (!isAnode && !isPW) continue; 

    double xw = e.x0;
    if (g.periodicX) {
        const double k  = std::round((x - e.x0) / g.pitchX);
        xw = e.x0 + k * g.pitchX;
    }
    const double d  = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) { 
        out.d_surface_cm = d; out.r_cm = e.radius; out.xw_cm = xw; out.yw_cm = e.y; 
        out.hasWire = true; out.wire_type = isAnode ? 0 : 1; out.wire_name = e.name; 
    }
  }
  return out;
}

// ==================== L0 Calculation Logic ====================
struct TimeHist {
  std::vector<double> t_center_ns; std::vector<double> dt_ns; std::vector<double> count;
};

static bool LoadTimeHistogram(const char* fname, TimeHist& h) {
  std::ifstream fin(fname);
  if (!fin) return false;
  std::string line;
  while (std::getline(fin, line)) {
    if (line.empty() || line[0] == '#') continue;
    bool hasAlpha = false; for (char c : line) if (std::isalpha((unsigned char)c)) { hasAlpha = true; break; }
    if (hasAlpha) continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    std::stringstream ss(line);
    double tcen, dt, cnt;
    if (!(ss >> tcen >> dt >> cnt)) continue;
    h.t_center_ns.push_back(tcen); h.dt_ns.push_back(dt); h.count.push_back(cnt);
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
  if (leaf == "analysis_L0_prim" || leaf == "analysis_L0" || leaf == "analysis") p = p.parent_path();
  fs::path q = p;
  while (!q.empty()) {
    if (q.filename().string().rfind("run_", 0) == 0) return q;
    if (!q.has_parent_path()) break; q = q.parent_path();
  }
  return p;
}

struct GainPlotInfo {
    int wType; double xw, yw, gain; int n_primary, n_total;
};

static std::string SanitizePathComponent(std::string s) {
  if (s.empty()) return "manual";
  for (char& c : s) {
    const bool ok =
        std::isalnum(static_cast<unsigned char>(c)) ||
        c == '_' || c == '-' || c == '.';
    if (!ok) c = '_';
  }
  return s;
}

static std::string ResolveArrayId() {
  if (const char* p = std::getenv("PBS_ARRAY_ID")) {
    return SanitizePathComponent(p);
  }
  return "manual";
}

static bool TryParsePositiveInt(const char* s, int& out) {
  if (s == nullptr || *s == '\0') return false;
  char* end = nullptr;
  long v = std::strtol(s, &end, 10);
  if (end == s || *end != '\0') return false;
  if (v <= 0 || v > std::numeric_limits<int>::max()) return false;
  out = static_cast<int>(v);
  return true;
}


struct TrackSegment {
  double x0 = 0.0;
  double y0 = 0.0;
  double x1 = 0.0;
  double y1 = 0.0;
  bool valid = false;
};

static TrackSegment ClipTrackToActiveArea(const Detector::Geometry& geo,
                                          double slope,
                                          double intercept) {
  constexpr double kEps = 1.0e-9;
  std::vector<std::pair<double, double>> points;

  auto inside = [&](double x, double y) {
    return x >= geo.xmin - kEps && x <= geo.xmax + kEps &&
           y >= geo.ymin - kEps && y <= geo.ymax + kEps;
  };
  auto pushPoint = [&](double x, double y) {
    if (!std::isfinite(x) || !std::isfinite(y) || !inside(x, y)) return;
    for (const auto& [px, py] : points) {
      if (std::hypot(px - x, py - y) < 1.0e-7) return;
    }
    points.push_back({x, y});
  };

  pushPoint(geo.xmin, slope * geo.xmin + intercept);
  pushPoint(geo.xmax, slope * geo.xmax + intercept);
  if (std::abs(slope) < kEps) {
    if (intercept >= geo.ymin - kEps && intercept <= geo.ymax + kEps) {
      pushPoint(geo.xmin, intercept);
      pushPoint(geo.xmax, intercept);
    }
  } else {
    pushPoint((geo.ymin - intercept) / slope, geo.ymin);
    pushPoint((geo.ymax - intercept) / slope, geo.ymax);
  }

  TrackSegment segment;
  if (points.size() < 2) return segment;

  size_t iBest = 0;
  size_t jBest = 1;
  double bestDistance = -1.0;
  for (size_t i = 0; i < points.size(); ++i) {
    for (size_t j = i + 1; j < points.size(); ++j) {
      const double distance = std::hypot(points[j].first - points[i].first,
                                         points[j].second - points[i].second);
      if (distance > bestDistance) {
        bestDistance = distance;
        iBest = i;
        jBest = j;
      }
    }
  }

  segment.x0 = points[iBest].first;
  segment.y0 = points[iBest].second;
  segment.x1 = points[jBest].first;
  segment.y1 = points[jBest].second;
  segment.valid = true;
  if (segment.x0 > segment.x1 ||
      (std::abs(segment.x0 - segment.x1) < kEps && segment.y0 > segment.y1)) {
    std::swap(segment.x0, segment.x1);
    std::swap(segment.y0, segment.y1);
  }
  return segment;
}

// ---- main ----------------------------------------------------
int main(int argc, char** argv) {
  //TApplication app("app", &argc, argv);
  
  // ★ バッチジョブ向けに画面の乱立を防ぐためバッチモードをON(非表示)にする
  gROOT->SetBatch(kTRUE);  gStyle->SetOptStat(0);

  // ★ 引数処理:
  //   argv[1] = worker/job ID
  //   argv[2] = eventsPerJob または time histogram file
  //   argv[3] = time histogram file または eventsPerJob
  //
  // 例:
  //   ./track_reco_avalanche_ap_wirecath 1
  //   ./track_reco_avalanche_ap_wirecath 1 10
  //   ./track_reco_avalanche_ap_wirecath 1 root/.../t_hist_nt.csv 10
  //   EVENTS_PER_JOB=10 ./track_reco_avalanche_ap_wirecath 1
  int job_seed = 1;
  int eventsPerJob = 1;
  std::string timeHistFile;

  if (argc > 1) {
      int parsedJob = 0;
      if (TryParsePositiveInt(argv[1], parsedJob)) {
          job_seed = parsedJob;
      } else {
          std::cerr << "[WARN] invalid job id: " << argv[1]
                    << " -> use job_seed = 1" << std::endl;
      }
  }

  if (const char* envEvents = std::getenv("EVENTS_PER_JOB")) {
      int parsedEvents = 0;
      if (TryParsePositiveInt(envEvents, parsedEvents)) {
          eventsPerJob = parsedEvents;
      } else {
          std::cerr << "[WARN] invalid EVENTS_PER_JOB=" << envEvents
                    << " -> use eventsPerJob = 1" << std::endl;
      }
  }
  
  #if defined(GEOM_WIRE)
    const char* geomTag = "wirecath";
    const char* defaultTimeHistFile = "root/run_20260322_205404_ap_wire/analysis_L0_prim/t_hist_nt.csv";
  #elif defined(GEOM_PLATE)
    const char* geomTag = "plate";
    const char* defaultTimeHistFile = "root/run_20260322_202630_ap_plate/analysis_L0_prim/t_hist_nt.csv";
  #else
    const char* geomTag = "unknown";
    const char* defaultTimeHistFile = "";
  #endif

  timeHistFile = defaultTimeHistFile;

  // Flexible optional arguments:
  //   argv[2] can be either eventsPerJob or timeHistFile.
  //   argv[3] can be either eventsPerJob or timeHistFile.
  if (argc > 2) {
      int parsedEvents = 0;
      if (TryParsePositiveInt(argv[2], parsedEvents)) {
          eventsPerJob = parsedEvents;
      } else {
          timeHistFile = argv[2];
      }
  }

  if (argc > 3) {
      int parsedEvents = 0;
      if (TryParsePositiveInt(argv[3], parsedEvents)) {
          eventsPerJob = parsedEvents;
      } else {
          timeHistFile = argv[3];
      }
  }

  const int evStart = (job_seed - 1) * eventsPerJob + 1;
  const int evEnd   = evStart + eventsPerJob;

  std::printf("[env] Job Seed (ID) = %d\n", job_seed);
  std::printf("[env] eventsPerJob = %d\n", eventsPerJob);
  std::printf("[env] event range  = %d to %d\n", evStart, evEnd - 1);
  std::printf("[env] timeHistFile=%s\n", timeHistFile.c_str());

  // 1. ジオメトリとカソード電圧 (vCat) 取得
  MediumMagboltz gas; gas.LoadGasFile("ic4H10_100_0.1atm.gas");
  Detector::Geometry geo;
  double vCat_val = 0.0;

#if defined(GEOM_PLATE)
  // ===== Plane cathode 用 =====
  const double pitchSense    = 0.40;   // [cm]
  const double gap           = 0.50;   // [cm]
  const double vAn           = 0.0;    // [V]
  const double vCat          = -1500.0; // [V]
  const bool   pwSameAsAnode = false;
  const int    nWires        = 4;

  const double phase = 0.0;      // [cm]
  const double rAn   = 1.0e-3;   // [cm]
  const double rPW   = 2.5e-3;   // [cm]
  const double vPW   = -121.48; // [V]
  //const double vPW   = vCat * (55.0 / 700.0);

  geo = Detector::PAP_PlaneCathode_Periodic(
      pitchSense,
      gap,
      phase,
      rAn,
      rPW,
      vAn,
      vPW,
      vCat,
      pwSameAsAnode,
      nWires
  );
  vCat_val = vCat;

#elif defined(GEOM_WIRE)
  // ===== Wire cathode 用 =====
  const double pitchSense    = 0.40;   // [cm]
  const double gap           = 0.50;   // [cm]
  const double vAn           = 0.0;    // [V]
  const double vCat          = -1500.0; // [V]
  const bool   pwSameAsAnode = false;
  const int    nWires        = 4;

  const double pwOffset   = 0.0;     // [cm]
  const double pitchCath  = 0.20;    // [cm]
  const double cathPhase  = 0;     // [cm]
  const double rAn        = 1.0e-3;  // [cm]
  const double rPW        = 2.5e-3;  // [cm]
  const double rCat       = 2.5e-3;  // [cm]
  const bool   printDebug = true;
  const double vPW   = -109.20; // [V]
  //const double vPW        = vCat * (50.0 / 700.0);

  geo = Detector::PAP_WireCathode_Periodic(
      pitchSense,
      gap,
      pwOffset,
      pitchCath,
      cathPhase,
      rAn,
      rPW,
      rCat,
      vAn,
      vPW,
      vCat,
      pwSameAsAnode,
      printDebug,
      nWires
  );
  vCat_val = vCat;

#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

  auto comp = Detector::BuildField(geo); comp->SetMedium(&gas);

  // 2. 出力ディレクトリ作成
  fs::path csvPath(timeHistFile);
  fs::path runDir  = GuessRunDirFromCsv(csvPath);
  char vol_dir_name[64]; std::snprintf(vol_dir_name, sizeof(vol_dir_name), "vCat_%.0fV", vCat_val);

  const std::string arrayId = ResolveArrayId();
  std::printf("[env] PBS_ARRAY_ID = %s\n", arrayId.c_str());

  fs::path vCatDir = runDir / vol_dir_name;
  fs::path arrayOutDir = vCatDir / ("array_" + arrayId);
  std::string figDir = (arrayOutDir / ("png_" + std::string(geomTag))).string();
  gSystem->mkdir(figDir.c_str(), kTRUE);
  fs::create_directories(fs::path(figDir) / "waveform");
  fs::create_directories(fs::path(figDir) / "gain");

  // track も array_<PBS_ARRAY_ID> 以下に保存する
  fs::path trackDir = arrayOutDir / "track";
  fs::create_directories(trackDir);

  // 3. データ記録用ROOTファイルセットアップ (★ファイル名にjob_seedを付与)
  std::string rootOutPath = (trackDir / Form("track_results_job%d.root", job_seed)).string();
  TFile* fOut = new TFile(rootOutPath.c_str(), "RECREATE"); fOut->cd(); 
  TTree* tree = new TTree("tree", "Track Reconstruction Stats"); tree->SetDirectory(fOut); 
  double t_x_true, t_x_reco, t_diff, t_a_true, t_b_true, t_a_reco, t_b_reco, t_chi2;
  int t_event_id;
  int t_nHits;
  tree->Branch("event_id", &t_event_id); tree->Branch("a_true", &t_a_true); tree->Branch("b_true", &t_b_true);
  tree->Branch("x_true", &t_x_true); tree->Branch("a_reco", &t_a_reco); tree->Branch("b_reco", &t_b_reco);
  tree->Branch("x_reco", &t_x_reco); tree->Branch("diff_x", &t_diff); tree->Branch("chi2",&t_chi2);
  tree->Branch("nHits", &t_nHits);

  // ===== Additional event-level diagnostic branches =====
  // 1 event = 1 entry.
  double t_theta_true, t_theta_reco;
  double t_chi2_ndf;
  double t_mean_t_hit, t_mean_L_meas, t_mean_gain;
  double t_rms_hit_residual;
  double t_max_abs_hit_residual;
  int    t_n_side_wrong;

  // Active-track segment used for primary-electron injection.
  double t_active_x0, t_active_y0, t_active_x1, t_active_y1;
  int    t_nClusters;

  tree->Branch("theta_true", &t_theta_true);
  tree->Branch("theta_reco", &t_theta_reco);
  tree->Branch("chi2_ndf", &t_chi2_ndf);
  tree->Branch("mean_t_hit", &t_mean_t_hit);
  tree->Branch("mean_L_meas", &t_mean_L_meas);
  tree->Branch("mean_gain", &t_mean_gain);
  tree->Branch("rms_hit_residual", &t_rms_hit_residual);
  tree->Branch("max_abs_hit_residual", &t_max_abs_hit_residual);
  tree->Branch("n_side_wrong", &t_n_side_wrong);
  tree->Branch("active_x0", &t_active_x0);
  tree->Branch("active_y0", &t_active_y0);
  tree->Branch("active_x1", &t_active_x1);
  tree->Branch("active_y1", &t_active_y1);
  tree->Branch("nClusters", &t_nClusters);

  // ===== Hit-level diagnostic tree =====
  // 1 hit = 1 entry. Used for correlations:
  // diff_x vs t_hit, diff_x vs wire_x, hit_residual vs gain, etc.
  TTree* hitTree = new TTree("hitTree", "Hit-level reconstruction information");
  hitTree->SetDirectory(fOut);

  int h_event_id;
  int h_nHits;
  int h_hit_id;
  int h_wire_index;
  int h_selected_side;  // +1: up, -1: down
  int h_true_side;      // +1: up, -1: down
  int h_side_ok;        // 1: correct, 0: wrong

  double h_diff_x;
  double h_chi2;
  double h_wire_x;
  double h_wire_y;
  double h_t_hit;
  double h_L_meas;
  double h_gain;
  int    h_n_primary;
  int    h_n_total;
  double h_y_selected;
  double h_y_reco;
  double h_y_true;
  double h_hit_residual;
  double h_truth_residual;

  hitTree->Branch("event_id", &h_event_id);
  hitTree->Branch("nHits", &h_nHits);
  hitTree->Branch("hit_id", &h_hit_id);
  hitTree->Branch("wire_index", &h_wire_index);
  hitTree->Branch("selected_side", &h_selected_side);
  hitTree->Branch("true_side", &h_true_side);
  hitTree->Branch("side_ok", &h_side_ok);

  hitTree->Branch("diff_x", &h_diff_x);
  hitTree->Branch("chi2", &h_chi2);
  hitTree->Branch("wire_x", &h_wire_x);
  hitTree->Branch("wire_y", &h_wire_y);
  hitTree->Branch("t_hit", &h_t_hit);
  hitTree->Branch("L_meas", &h_L_meas);
  hitTree->Branch("gain", &h_gain);
  hitTree->Branch("n_primary", &h_n_primary);
  hitTree->Branch("n_total", &h_n_total);
  hitTree->Branch("y_selected", &h_y_selected);
  hitTree->Branch("y_reco", &h_y_reco);
  hitTree->Branch("y_true", &h_y_true);
  hitTree->Branch("hit_residual", &h_hit_residual);
  hitTree->Branch("truth_residual", &h_truth_residual);

  Sensor sensor; sensor.AddComponent(comp.get()); sensor.SetArea(geo.xmin, geo.ymin, -5.0, geo.xmax, geo.ymax, +5.0);
  int numAnodes = 0;
  for (const auto& e : geo.electrodes) {
      auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
      if (e.kind == Detector::ElectrodeKind::WireRow && lower(e.name).find("anode") != std::string::npos) {
          sensor.AddElectrode(comp.get(), e.name); numAnodes++;
      }
  }
  sensor.SetTransferFunction(TransferFunction);

  const double gap_val = DetectGap(geo);
  TimeHist th;
  if (!LoadTimeHistogram(timeHistFile.c_str(), th)) {
      std::cerr << "[ERROR] failed to read time histogram: " << timeHistFile << std::endl;
      return 1;
  }
  const double t0_sig_ns = 0.0, dt_sig_ns = 0.5; 
  const double tmax_sig  = (th.t_center_ns.empty() ? 200.0 : th.t_center_ns.back()) + 150.0; 
  const int nSigBins = std::max(200, (int)std::ceil((tmax_sig - t0_sig_ns) / dt_sig_ns));
  sensor.SetTimeWindow(t0_sig_ns, dt_sig_ns, nSigBins);

  // ★ workerごとに乱数列を変える。
  // 同じworker内の複数イベントは、この乱数列を連続して使う。
  TRandom3 gen(job_seed); 

  TFile* fDriftOut = nullptr;
  TTree* tDrift = nullptr;
  int d_event_id, d_wire_idx; double d_x_start, d_y_start, d_x_end, d_y_end, d_t_drift, d_dist;
#if SAVE_DRIFT_ROOT
  {
    std::string driftRootOutPath = (trackDir / Form("drift_info_job%d.root", job_seed)).string();
    fDriftOut = new TFile(driftRootOutPath.c_str(), "RECREATE"); fDriftOut->cd();
    tDrift = new TTree("driftTree", "Electron Drift Info"); tDrift->SetDirectory(fDriftOut); 
    tDrift->Branch("event_id", &d_event_id); tDrift->Branch("wire_idx", &d_wire_idx);
    tDrift->Branch("x_start", &d_x_start); tDrift->Branch("y_start", &d_y_start); tDrift->Branch("x_end", &d_x_end);
    tDrift->Branch("y_end", &d_y_end); tDrift->Branch("t_drift", &d_t_drift); tDrift->Branch("dist", &d_dist);
  }
#endif

  TFile* fGainOut = nullptr;
  TTree* tGain = nullptr;
  int g_event_id, g_wire_type, g_n_primary, g_n_total; double g_wire_tag, g_x_wire, g_y_wire, g_gain;
#if SAVE_GAIN_ROOT
  {
    char gain_filename[128]; std::snprintf(gain_filename, sizeof(gain_filename), "gain_info_vCat%.0fV_job%d.root", vCat_val, job_seed);
    std::string gainRootOutPath = (trackDir / gain_filename).string();
    fGainOut = new TFile(gainRootOutPath.c_str(), "RECREATE"); fGainOut->cd();
    tGain = new TTree("gainTree", "Avalanche Gain Info"); tGain->SetDirectory(fGainOut);
    tGain->Branch("event_id", &g_event_id); tGain->Branch("wire_type", &g_wire_type); tGain->Branch("wire_tag", &g_wire_tag);   
    tGain->Branch("x_wire", &g_x_wire); tGain->Branch("y_wire", &g_y_wire); tGain->Branch("n_primary", &g_n_primary); 
    tGain->Branch("n_total", &g_n_total); tGain->Branch("gain", &g_gain);
  }
#endif

  TFile* fVerifyOut = nullptr;
  TTree* tVerify = nullptr;
  int v_event_id; 
  char v_wire_name[64];
  double v_x_start, v_y_start, v_x_end, v_y_end;
  double v_wp_start, v_wp_end, v_wp_diff;
  double v_raw_charge_e, v_raw_charge_i, v_raw_charge_total;
#if SAVE_VERIFY_ROOT
  {
    std::string verifyRootOutPath = (trackDir / Form("verify_info_job%d.root", job_seed)).string();
    fVerifyOut = new TFile(verifyRootOutPath.c_str(), "RECREATE"); fVerifyOut->cd();
    tVerify = new TTree("verifyTree", "Shockley-Ramo Verification"); tVerify->SetDirectory(fVerifyOut);

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
  }
#endif

  // =========================================================
  // (A) Map Gen (Field Lines & Maps)
  // =========================================================
  // ★ シード1のジョブでのみ重いマップ生成を行う
  bool make_maps = (job_seed == 1); 
  if (make_maps) {
    std::printf("[map] Generating high-resolution field maps & lines...\n");
    const int NX = 1500, NY = 1500;
    TH2D hE("|E|","Electric Field |E|;X [cm];Y [cm]", NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    TH2D hV("V","Electric Potential V;X [cm];Y [cm]", NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    hE.SetStats(0); hV.SetStats(0);

    double minPos = 1e98, maxVal_global = 0.0;
    for (int ix = 1; ix <= NX; ++ix) {
      const double x = hE.GetXaxis()->GetBinCenter(ix);
      for (int iy = 1; iy <= NY; ++iy) {
        const double y = hE.GetYaxis()->GetBinCenter(iy);
        double ex=0, ey=0, ez=0, V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(x, y, 0, ex, ey, ez, V, m, status);
        if (status == 0) {
          const double E = std::hypot(ex, std::hypot(ey, ez));
          hE.SetBinContent(ix, iy, E); hV.SetBinContent(ix, iy, V);
          if (E > 0 && E < minPos) minPos = E;
          if (E > maxVal_global) maxVal_global = E;
        }
      }
    }

    std::vector<double> anode_xs, pw_xs;
    for (const auto& e : geo.electrodes) {
      if (e.kind != Detector::ElectrodeKind::WireRow) continue;
      auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
      std::string ename = lower(e.name);
      bool isAnode = ename.find("anode") != std::string::npos;
      bool isPW    = ename.find("pw")    != std::string::npos;
      if (!isAnode && !isPW) continue;

      if (!geo.periodicX) {
        if (e.x0 >= geo.xmin && e.x0 <= geo.xmax) {
            if (isAnode) anode_xs.push_back(e.x0);
            if (isPW) pw_xs.push_back(e.x0);
        }
      } else {
        const double pitch = geo.pitchX;
        int nL = int(std::floor((geo.xmin - e.x0) / pitch));
        int nR = int(std::ceil ((geo.xmax - e.x0) / pitch));
        for (int n = nL; n <= nR; ++n) {
          const double x = e.x0 + n * pitch;
          if (x >= geo.xmin && x <= geo.xmax) {
              if (isAnode) anode_xs.push_back(x);
              if (isPW) pw_xs.push_back(x);
          }
        }
      }
    }
    auto remove_duplicates = [](std::vector<double>& vec) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end(), [](double a, double b){ return std::abs(a-b)<1e-5; }), vec.end());
    };
    remove_duplicates(anode_xs); remove_duplicates(pw_xs);

    struct MaxEInfo { double x, y, E; int type; };
    std::vector<MaxEInfo> maxE_list;

// --- 2.1 アノードの1Dプロファイル ---
    TCanvas* cProfileAnode = new TCanvas("cProfileAnode", "Anode E Field Profiles", 800, 600);
    cProfileAnode->SetGrid(); cProfileAnode->SetLogy(0); 
    TMultiGraph* mg_anode = new TMultiGraph();
    mg_anode->SetTitle("Electric Field Profile at Anode X positions;Y [cm];|E| [V/cm]");
    
    int color_idx = 1;
    // ★ 描画用のテキスト情報を一時保存する構造体
    struct LabelInfo { std::string text; int color; };
    std::vector<LabelInfo> labels10k;

    for (double xw : anode_xs) {
      std::vector<double> y_vals, e_vals; double local_maxE = 0.0, local_maxY = 0.0;
      
      double y_10k_min = 1e9;
      double y_10k_max = -1e9;
      bool found_10k = false;

      for (int i = 0; i <= 2000; ++i) {
        double y = geo.ymin + i * (geo.ymax - geo.ymin) / 2000.0;
        double ex=0, ey=0, ez=0, V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(xw, y, 0, ex, ey, ez, V, m, status);
        
        if (status == 0) {
          double E = std::hypot(ex, std::hypot(ey, ez));
          y_vals.push_back(y); e_vals.push_back(E);
          if (E > local_maxE) { local_maxE = E; local_maxY = y; }
          
          if (E >= 10000.0) {
              if (y < y_10k_min) y_10k_min = y;
              if (y > y_10k_max) y_10k_max = y;
              found_10k = true;
          }
        }
      }
      maxE_list.push_back({xw, local_maxY, local_maxE, 0});
      
      if (!y_vals.empty()) {
        TGraph* gr = new TGraph(y_vals.size(), y_vals.data(), e_vals.data());
        gr->SetTitle(Form("Anode at x = %.2f cm", xw));
        gr->SetLineColor(color_idx); gr->SetLineWidth(2); mg_anode->Add(gr, "L");
        
        // ★ テキスト情報をベクターに保存しておく
        if (found_10k) {
            labels10k.push_back({Form("x=%.2f: Y [%.4f, %.4f] cm", xw, y_10k_min, y_10k_max), color_idx});
        }

        color_idx++; 
        if(color_idx == 5) color_idx++; if(color_idx > 9) color_idx = 1;
      }
    }

    cProfileAnode->cd();
    if (mg_anode->GetListOfGraphs() != nullptr) {
        // ★ 先にグラフ本体と軸を描画（ここで座標系が確定する）
        mg_anode->Draw("A"); 
        cProfileAnode->BuildLegend(0.65, 0.75, 0.90, 0.90);

        // ★ グラフ描画後にテキストを NDC (0~1の相対座標) で画面に配置する
        TLatex latex_10k; 
        latex_10k.SetNDC(); // 絶対座標系を有効化
        latex_10k.SetTextSize(0.025); 
        latex_10k.SetTextAlign(12); // 左寄せ
        
        double text_y_pos = 0.70; // 凡例の少し下から書き始める
        for (const auto& lbl : labels10k) {
            latex_10k.SetTextColor(lbl.color);
            latex_10k.DrawLatex(0.55, text_y_pos, lbl.text.c_str());
            text_y_pos -= 0.04; // 次の行が重ならないように下へずらす
        }
    }
    cProfileAnode->SaveAs((figDir + "/field_profiles_1D_anode.png").c_str());
    
    // --- 2.2 PWの1Dプロファイル ---
    TCanvas* cProfilePW = new TCanvas("cProfilePW", "PW E Field Profiles", 800, 600);
    cProfilePW->SetGrid(); cProfilePW->SetLogy(0); 
    TMultiGraph* mg_pw = new TMultiGraph();
    mg_pw->SetTitle("Electric Field Profile at Potential Wire X positions;Y [cm];|E| [V/cm]");
    
    color_idx = 1;
    for (double xw : pw_xs) {
      std::vector<double> y_vals, e_vals; double local_maxE = 0.0, local_maxY = 0.0;
      for (int i = 0; i <= 2000; ++i) {
        double y = geo.ymin + i * (geo.ymax - geo.ymin) / 2000.0;
        double ex=0, ey=0, ez=0, V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(xw, y, 0, ex, ey, ez, V, m, status);
        if (status == 0) {
          double E = std::hypot(ex, std::hypot(ey, ez));
          y_vals.push_back(y); e_vals.push_back(E);
          if (E > local_maxE) { local_maxE = E; local_maxY = y; }
        }
      }
      maxE_list.push_back({xw, local_maxY, local_maxE, 1}); 
      if (!y_vals.empty()) {
        TGraph* gr = new TGraph(y_vals.size(), y_vals.data(), e_vals.data());
        gr->SetTitle(Form("PW at x = %.2f cm", xw));
        gr->SetLineColor(color_idx++); gr->SetLineWidth(2); mg_pw->Add(gr, "L");
        if(color_idx == 5) color_idx++; if(color_idx > 9) color_idx = 1;
      }
    }
    cProfilePW->cd();
    if (mg_pw->GetListOfGraphs() != nullptr) {
        mg_pw->Draw("A"); cProfilePW->BuildLegend(0.70, 0.70, 0.90, 0.90);
    }
    cProfilePW->SaveAs((figDir + "/field_profiles_1D_pw.png").c_str());
    
    // --- 2.3 Townsend係数(alpha)の1Dプロファイル ---
    TCanvas* cProfileAlpha = new TCanvas("cProfileAlpha", "Townsend Coefficient Profile", 800, 600);
    cProfileAlpha->SetGrid(); 
    cProfileAlpha->SetLogy(1); // alphaは急激に変化するためLogスケールが見やすい

    TMultiGraph* mg_alpha = new TMultiGraph();
    mg_alpha->SetTitle("Townsend Coefficient (#alpha) at Anode X;Y [cm];#alpha [1/cm]");

    int color_idx_alpha = 1;
    struct AlphaLabelInfo { std::string text; int color; };
    std::vector<AlphaLabelInfo> labelsAlpha;

    for (double xw : anode_xs) {
      std::vector<double> y_vals_a, alpha_vals; 
      double y_alpha_min = 1e9, y_alpha_max = -1e9;
      bool found_alpha = false;

      for (int i = 0; i <= 2000; ++i) {
        double y = geo.ymin + i * (geo.ymax - geo.ymin) / 2000.0;
        double ex=0, ey=0, ez=0, V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(xw, y, 0, ex, ey, ez, V, m, status);
        
        if (status == 0 && m != nullptr) {
          double alpha = 0.0;
          // 磁場(Bx, By, Bz)はゼロとしてTownsend係数を取得
          m->ElectronTownsend(ex, ey, ez, 0., 0., 0., alpha);
          
          // Logスケール描画のため、ゼロや極端に小さい値は弾く
          if (alpha > 1e-3) {
              y_vals_a.push_back(y);
              alpha_vals.push_back(alpha);
              
              // 電離が「実質的に意味を持つレベル」の目安として alpha > 1.0 cm^-1 を記録
              if (alpha > 1.0) {
                  if (y < y_alpha_min) y_alpha_min = y;
                  if (y > y_alpha_max) y_alpha_max = y;
                  found_alpha = true;
              }
          }
        }
      }
      
      if (!y_vals_a.empty()) {
        TGraph* gr = new TGraph(y_vals_a.size(), y_vals_a.data(), alpha_vals.data());
        gr->SetTitle(Form("Anode at x = %.2f cm", xw));
        gr->SetLineColor(color_idx_alpha); gr->SetLineWidth(2); mg_alpha->Add(gr, "L");
        
        if (found_alpha) {
            labelsAlpha.push_back({Form("x=%.2f: #alpha>1.0 at Y [%.4f, %.4f] cm", xw, y_alpha_min, y_alpha_max), color_idx_alpha});
        }

        color_idx_alpha++; 
        if(color_idx_alpha == 5) color_idx_alpha++; if(color_idx_alpha > 9) color_idx_alpha = 1;
      }
    }

    cProfileAlpha->cd();
    if (mg_alpha->GetListOfGraphs() != nullptr) {
        mg_alpha->Draw("A"); 
        cProfileAlpha->BuildLegend(0.65, 0.75, 0.90, 0.90);

        TLatex latex_alpha; 
        latex_alpha.SetNDC(); 
        latex_alpha.SetTextSize(0.025); 
        latex_alpha.SetTextAlign(12); 
        
        double text_y_pos = 0.70; 
        for (const auto& lbl : labelsAlpha) {
            latex_alpha.SetTextColor(lbl.color);
            latex_alpha.DrawLatex(0.55, text_y_pos, lbl.text.c_str());
            text_y_pos -= 0.04; 
        }
    }
    cProfileAlpha->SaveAs((figDir + "/field_profiles_1D_alpha.png").c_str());

    // --- 2.4 Python積分計算用の空間プロファイル出力 (中心のワイヤーを対象) ---
    std::string csv_out_path = figDir + "/avalanche_profile.csv";
    std::ofstream fout_csv(csv_out_path);
    fout_csv << "r_cm,E_Vcm,alpha_invcm,phi_w\n";
    
    // 中心のワイヤー(x=0)の名前を探す
    std::string center_wire_name = "";
    for (const auto& e : geo.electrodes) {
        auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
        if (e.kind == Detector::ElectrodeKind::WireRow && lower(e.name).find("anode") != std::string::npos) {
            if (std::abs(e.x0) < 1e-3) { center_wire_name = e.name; break; }
        }
    }
    if(center_wire_name.empty()) center_wire_name = "anode_0"; // 念のためのフォールバック

    const double r_anode_cm = 0.0006;  // アノード半径 (6 um)
    const double r_start_cm = 0.0400;  // 積分開始位置 (400 um)
    const int n_steps_csv = 4000;      // 0.1 um刻みの超高解像度でスキャン
    const double dr_csv = (r_start_cm - r_anode_cm) / n_steps_csv;

    for (int i = 0; i <= n_steps_csv; ++i) {
        double r = r_anode_cm + i * dr_csv;
        double y = r; // 中心(0,0)から+Y方向へスキャン
        
        double ex=0, ey=0, ez=0, V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(0.0, y, 0.0, ex, ey, ez, V, m, status);
        
        double E = std::hypot(ex, std::hypot(ey, ez));
        double alpha = 0.0;
        if (status == 0 && m != nullptr) {
            m->ElectronTownsend(ex, ey, ez, 0., 0., 0., alpha);
        }
        
        double phi_w = sensor.WeightingPotential(0.0, y, 0.0, center_wire_name);
        fout_csv << r << "," << E << "," << alpha << "," << phi_w << "\n";
    }
    fout_csv.close();
    std::printf("[map] Exported spatial profile to %s\n", csv_out_path.c_str());

    // --- 3. 2D電場マップ（単独） ---
    TCanvas* cE = new TCanvas("cE", "|E| Map", 1000, 780);
    cE->SetRightMargin(0.15); cE->SetGrid(); cE->SetLogz(1);
    if (minPos < 1e97) hE.SetMinimum(minPos * 0.9); hE.SetMaximum(maxVal_global * 1.05);
    
    // ★ ここで描画レンジを絞る（ズーム）
    hE.GetXaxis()->SetRangeUser(-0.2, 0.2);
    hE.GetYaxis()->SetRangeUser(-0.25, 0.25);
    
    hE.Draw("COLZ"); DrawElectrodes(geo);

    TLatex latexMaxE; latexMaxE.SetTextSize(0.025); latexMaxE.SetTextAlign(21);
    for (const auto& maxE : maxE_list) {
      TMarker* mMax = new TMarker(maxE.x, maxE.y, 29);
      mMax->SetMarkerColor(maxE.type == 0 ? kRed : kBlue); mMax->SetMarkerSize(2.0); mMax->Draw("same");
      latexMaxE.SetTextColor(maxE.type == 0 ? kRed : kBlue);
      latexMaxE.DrawLatex(maxE.x, maxE.y + 0.04, Form("%.1e", maxE.E));
    }
    cE->SaveAs((figDir + "/emap_absE_highres.png").c_str());

    // --- 4. コンボマップ ---
    TCanvas* cCombo = new TCanvas("cComboMap","Combo Map", 1000, 780); 
    cCombo->SetRightMargin(0.15); cCombo->SetGrid(); cCombo->SetLogz(1);
    hE.Draw("COLZ"); // RangeUserの指定はhEに反映済み
    hV.SetContour(40); hV.SetLineColor(kGray+2); hV.SetLineWidth(1); hV.Draw("CONT3 SAME");
    DrawElectrodes(geo);

    for (const auto& maxE : maxE_list) {
      TMarker* mMax = new TMarker(maxE.x, maxE.y, 29);
      mMax->SetMarkerColor(maxE.type == 0 ? kRed : kBlue); mMax->SetMarkerSize(2.0); mMax->Draw("same");
      latexMaxE.SetTextColor(maxE.type == 0 ? kRed : kBlue);
      latexMaxE.DrawLatex(maxE.x, maxE.y + 0.04, Form("%.1e", maxE.E));
    }

    ViewField vfFL; vfFL.SetCanvas(cCombo); vfFL.SetComponent(comp.get());
    vfFL.SetArea(geo.xmin, geo.ymin, geo.xmax, geo.ymax);
    const double gapv = DetectGap(geo);
    auto inMedium = [&](double x, double y){
      double ex=0, ey=0, ez=0, V=0; Medium* m=nullptr; int status=0;
      comp->ElectricField(x,y,0.0,ex,ey,ez,V,m,status); return (m!=nullptr && status==0);
    };
    auto make_line_seeds = [&](double yseed, std::vector<double>& xs2, std::vector<double>& ys2, std::vector<double>& zs2){
      const double Lx = (geo.xmax - 1e-3) - (geo.xmin + 1e-3);
      for (int i = 0; i < 64; ++i) {
        const double x = geo.xmin + 1e-3 + (i + 0.5) * Lx / 64;
        if (!inMedium(x, yseed)) continue;
        xs2.push_back(x); ys2.push_back(yseed); zs2.push_back(0.0);
      }
    };
    std::vector<double> xs2, ys2, zs2;
    make_line_seeds(+0.98 * gapv, xs2, ys2, zs2);
    make_line_seeds(-0.98 * gapv, xs2, ys2, zs2);
    if (!xs2.empty()) vfFL.PlotFieldLines(xs2, ys2, zs2, true, true);

    cCombo->SaveAs((figDir + "/fieldmap_combo_highres.png").c_str());
  }

  // =========================================================
  // イベントループ
  //   worker/job ID = 1 -> ev = 1 ... eventsPerJob
  //   worker/job ID = 2 -> ev = eventsPerJob+1 ... 2*eventsPerJob
  // =========================================================
  for (int ev = evStart; ev < evEnd; ++ev) {
        std::printf("\n>>> Processing Event/Job ID %d <<<\n", ev);
        
        // --- True track generation with angular straggling ---
        // 物理的には 57 度側に δθ が付く
        const double theta_las0_deg = 57.0;
        const double theta_las0_rad = theta_las0_deg * M_PI / 180.0;

        // 合計で 5 mrad なので、半幅は 2.5 mrad
        const double theta_total_rad = 5.0e-3;
        const double theta_half_rad  = 0.5 * theta_total_rad;

        // δθ を -2.5 mrad から +2.5 mrad の範囲で振る
        const double dtheta_rad = gen.Uniform(-theta_half_rad, theta_half_rad);

        // 57度側: 57° + δθ
        const double theta_las_rad = theta_las0_rad + dtheta_rad;

        // コード内の直線角度は負の 33 度側
        // 33° -> 33° - δθ なので、
        // theta_code = -(33° - δθ) = -33° + δθ
        const double true_angle_rad = -(0.5 * M_PI - theta_las_rad);
        const double true_angle_deg = true_angle_rad * 180.0 / M_PI;

        // 真トラック y = a_true x + b_true
        const double a_true = std::tan(true_angle_rad);
        const double b_true = 0.25 + gen.Uniform(-0.05, 0.05);

        const int nClusters = 100; // 879, 2707
        const auto activeTrack = ClipTrackToActiveArea(geo, a_true, b_true);
        if (!activeTrack.valid) {
          std::fprintf(stderr, "[WARN] true track is outside the active area. Skip event %d.\n", ev);
          continue;
        }

        AvalancheMC amcTrack; amcTrack.SetSensor(&sensor);
        TryEnableSignalCalculation(amcTrack); amcTrack.UseWeightingPotential(true);

        ViewDrift driftView; driftView.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, 0.1);

        struct ClusterInfo { double x, y; };
        std::vector<ClusterInfo> clusters;
        for (int i = 0; i < nClusters; ++i) {
          const double frac = gen.Uniform(0.0, 1.0); 
          const double x = activeTrack.x0 + frac * (activeTrack.x1 - activeTrack.x0);
          const double y = activeTrack.y0 + frac * (activeTrack.y1 - activeTrack.y0);
          clusters.push_back({x, y});
        }

        std::map<std::string, GainPlotInfo> globalGainMap;
        sensor.ClearSignal(); 

        for (size_t i = 0; i < clusters.size(); ++i) {
            auto& cl = clusters[i];

            #ifdef DEBUG_WAVEFORM
            if (i % 1 == 0) { 
                amcTrack.EnablePlotting(&driftView);
                amcTrack.AvalancheElectron(cl.x, cl.y, 0, 0); 
                amcTrack.DisablePlotting();
            } else { amcTrack.AvalancheElectron(cl.x, cl.y, 0, 0); }
            #else
            amcTrack.AvalancheElectron(cl.x, cl.y, 0, 0);
            #endif
            
            int n_e = amcTrack.GetNumberOfElectronEndpoints();
            if (n_e > 0) {
                double x0, y0, z0, t0, x1, y1, z1, t1; int status;
                amcTrack.GetElectronEndpoint(0, x0, y0, z0, t0, x1, y1, z1, t1, status);
                auto end_wire = NearestWireSurface(geo, x1, y1);
                if (end_wire.hasWire) {
                    std::string wName = end_wire.wire_name;
                    if (globalGainMap.find(wName) == globalGainMap.end()) {
                        globalGainMap[wName] = {end_wire.wire_type, end_wire.xw_cm, end_wire.yw_cm, 0.0, 0, 0};
                    }
                    globalGainMap[wName].n_primary += 1; globalGainMap[wName].n_total += n_e; 
                    
                    d_event_id = ev; d_wire_idx = std::round(end_wire.xw_cm / geo.pitchX); 
                    d_x_start = x0; d_y_start = y0; d_x_end = x1; d_y_end = y1;
                    d_t_drift = t1 - t0; d_dist = std::hypot(x1 - x0, y1 - y0); 
#if SAVE_DRIFT_ROOT
                    fDriftOut->cd(); tDrift->Fill();
#endif
                }
            }
        }

        // =========================================================
        // 検算データの取得 (ConvoluteSignals の直前に行う)
        // =========================================================
        for (auto& [wName, data] : globalGainMap) {
            std::strncpy(v_wire_name, wName.c_str(), sizeof(v_wire_name) - 1);
            v_event_id = ev;

            // 代表点として、最後のアバランシェの飛跡(d_x_start 等)を使用
            v_x_start = d_x_start; v_y_start = d_y_start;
            v_x_end = d_x_end;     v_y_end = d_y_end;

            // 重み付きポテンシャル (Weighting Potential) の取得 [0.0 ~ 1.0]
            v_wp_start = sensor.WeightingPotential(v_x_start, v_y_start, 0.0, wName);
            v_wp_end   = sensor.WeightingPotential(v_x_end, v_y_end, 0.0, wName);
            v_wp_diff  = v_wp_end - v_wp_start;

            // 生電流の積分 (Raw Signal Integration)
            v_raw_charge_e = 0.0;
            v_raw_charge_i = 0.0;
            for (int k = 0; k < nSigBins; ++k) {
                // GetElectronSignal / GetIonSignal は畳み込み前の純粋な誘導電流
                v_raw_charge_e += sensor.GetElectronSignal(wName, k) * dt_sig_ns;
                v_raw_charge_i += sensor.GetIonSignal(wName, k) * dt_sig_ns;
            }
            v_raw_charge_total = v_raw_charge_e + v_raw_charge_i;
            
#if SAVE_VERIFY_ROOT
            fVerifyOut->cd();
            tVerify->Fill();
#endif
        }
        
        // --- プリアンプ応答の畳み込み ---
        sensor.ConvoluteSignals();

        std::vector<GainPlotInfo> evGainData;
        for (auto& [wName, data] : globalGainMap) {
            data.gain = (data.n_primary > 0) ? (double)data.n_total / data.n_primary : 0.0;
            g_event_id = ev; g_wire_type = data.wType;
            g_wire_tag = data.xw;
            g_x_wire = data.xw; g_y_wire = data.yw;
            g_n_primary = data.n_primary; g_n_total = data.n_total; g_gain = data.gain;
#if SAVE_GAIN_ROOT
            fGainOut->cd(); tGain->Fill();
#endif
            evGainData.push_back(data);
        }

        struct HitCand {
          double xwire;
          double ywire;
          double y_up;
          double y_dn;
          double t_hit;
          double L_meas;
          int wire_index;

          double gain;
          int n_primary;
          int n_total;
        };
        std::vector<HitCand> hits;
        
        // 閾値判定 (4mV)
        const double signalThreshold = 0.05;

        for (const auto& e : geo.electrodes) {
            auto lower = [](std::string s){ std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s; };
            if (e.kind == Detector::ElectrodeKind::WireRow && lower(e.name).find("anode") != std::string::npos) {
                std::string wName = e.name; double xw = e.x0; double yw = e.y;

                #ifdef DEBUG_WAVEFORM
                {
                  TCanvas cWave("cWave", "", 600, 400); cWave.SetGrid();
                  std::vector<double> tx(nSigBins), vy(nSigBins);
                  for(int k=0; k<nSigBins; ++k) {
                    tx[k] = t0_sig_ns + k*dt_sig_ns; vy[k] = sensor.GetSignal(wName, k, true);
                  }
                  TGraph gW(nSigBins, tx.data(), vy.data());
                  gW.SetTitle((wName + " (Ev " + std::to_string(ev) + ");Time [ns];Signal [mV]").c_str());
                  gW.SetLineWidth(1); gW.SetLineColor(kBlue); gW.Draw("AL");

                  TLine lth_n(t0_sig_ns, -signalThreshold, tmax_sig, -signalThreshold);
                  lth_n.SetLineColor(kRed); lth_n.SetLineStyle(2); lth_n.Draw("same");

                  char wname_png[128]; 
                  std::snprintf(wname_png, sizeof(wname_png), "%s/waveform/wire_sig_ev%d_%s.png", figDir.c_str(), ev, wName.c_str());
                  cWave.SaveAs(wname_png);
                }
                #endif

                double t_hit = -1.0;
                for (int k = 0; k < nSigBins; ++k) {
                  //if (std::abs(sensor.GetSignal(wName, k, true)) > signalThreshold) {
                  if (sensor.GetSignal(wName, k, true) < -signalThreshold) {
                    t_hit = t0_sig_ns + k * dt_sig_ns; break; 
                  }
                }

                if (t_hit >= 0.0) {
                    double L_meas = CalculateL0Directly(th, gap_val, t_hit);
                    if (L_meas < 0) L_meas = 0;

                    int wire_index = 0;
                    if (geo.pitchX > 0.0) {
                      wire_index = (int)std::round(xw / geo.pitchX);
                    }

                    double hit_gain = 0.0;
                    int hit_n_primary = 0;
                    int hit_n_total = 0;

                    auto itGain = globalGainMap.find(wName);
                    if (itGain != globalGainMap.end()) {
                      hit_gain = itGain->second.gain;
                      hit_n_primary = itGain->second.n_primary;
                      hit_n_total = itGain->second.n_total;
                    }

                    hits.push_back({
                      xw,
                      yw,
                      yw + L_meas,
                      yw - L_meas,
                      t_hit,
                      L_meas,
                      wire_index,
                      hit_gain,
                      hit_n_primary,
                      hit_n_total
                    });
                }
            }
        }
        
        #ifdef DEBUG_WAVEFORM
        {
            TCanvas cGain("cGain", "Gain Map", 800, 500); cGain.SetGrid();
            TH2D frGain("frGain", Form("Avalanche Gain Map (Event %d);X [cm];Y [cm]", ev), 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
            frGain.SetStats(0); frGain.Draw(); DrawElectrodes(geo);

            TLatex latex; latex.SetTextSize(0.035); latex.SetTextAlign(21); 
            for (const auto& gdata : evGainData) {
                if (gdata.n_primary <= 0) continue;
                TMarker mt(gdata.xw, gdata.yw, 20);
                mt.SetMarkerColor(gdata.wType == 0 ? kRed : kBlue);   
                mt.SetMarkerSize(gdata.wType == 0 ? 1.5 : 1.0); mt.DrawClone("same");
                latex.SetTextColor(gdata.wType == 0 ? kRed : kBlue);
                latex.DrawLatex(gdata.xw, gdata.yw + 0.04, Form("%.0f", gdata.gain)); 
            }
            cGain.SaveAs(Form("%s/gain/gain_map_ev%d.png", figDir.c_str(), ev));
        }

        {
            TCanvas cClEnd("cClEnd", "cluster endpoint tracks", 1100, 850); cClEnd.SetGrid();
            TH2D frClEnd("frClEnd", "Electron Drift Lines;X [cm];Y [cm]", 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
            frClEnd.SetStats(0); frClEnd.Draw(); DrawElectrodes(geo);
            driftView.SetCanvas(&cClEnd); driftView.Plot(true, false);
            
            TF1 fTrueDrift("fTrueDrift","[0]*x+[1]", geo.xmin, geo.xmax);
            fTrueDrift.SetParameters(a_true, b_true);
            fTrueDrift.SetLineColor(kBlue); fTrueDrift.SetLineStyle(7); fTrueDrift.SetLineWidth(2);
            fTrueDrift.Draw("same");

            cClEnd.SaveAs(Form("%s/cluster_endtracks_ev%d.png", figDir.c_str(), ev));
        }
        #endif

        if (hits.size() < 3) continue; 

        int nHits = hits.size(); unsigned long long nComb = 1ULL << nHits;
        double best_a = 0.0, best_b = 0.0, min_chi2 = 1.0e300; unsigned long long best_mask = 0;
        
        for (unsigned long long mask = 0; mask < nComb; ++mask) {
            double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
            for (int i = 0; i < nHits; ++i) {
                double x = hits[i].xwire; double y = ((mask >> i) & 1) ? hits[i].y_dn : hits[i].y_up;
                Sx += x; Sy += y; Sxx += x * x; Sxy += x * y;
            }
            double delta = nHits * Sxx - Sx * Sx; if (std::abs(delta) < 1e-12) continue;
            double a_tmp = (nHits * Sxy - Sx * Sy) / delta; double b_tmp = (Sxx * Sy - Sx * Sxy) / delta;

            double chi2 = 0.0;
            for (int i = 0; i < nHits; ++i) {
                double x = hits[i].xwire; double y = ((mask >> i) & 1) ? hits[i].y_dn : hits[i].y_up;
                double diff = y - (a_tmp * x + b_tmp); chi2 += diff * diff;
            }
            if (chi2 < min_chi2) { min_chi2 = chi2; best_a = a_tmp; best_b = b_tmp; best_mask = mask; }
        }

        TCanvas cTrk("cTrk","track reco",1000,780); cTrk.SetGrid();
        TH2D frame("frame", "Reconstructed Track;X;Y", 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
        frame.SetStats(0); frame.Draw(); DrawElectrodes(geo);

        TF1 fTrue("fTrue","[0]*x+[1]", geo.xmin, geo.xmax);
        fTrue.SetParameters(a_true, b_true); fTrue.SetLineColor(kBlue); fTrue.SetLineStyle(2); fTrue.Draw("same");

        //for (const auto& h : hits) {
        //    double L = std::abs(h.y_up - h.ywire); TEllipse *el = new TEllipse(h.xwire, h.ywire, L, L);
        //    el->SetFillStyle(0); el->SetLineColor(kGreen+2); el->SetLineWidth(1); el->Draw("same");
        //}

        TF1 fReco("fReco","[0]*x+[1]", geo.xmin, geo.xmax);
        fReco.SetParameters(best_a, best_b); fReco.SetLineColor(kRed); fReco.SetLineWidth(2); fReco.Draw("same");

        for (int i = 0; i < nHits; ++i) {
            double x = hits[i].xwire; double y = ((best_mask >> i) & 1) ? hits[i].y_dn : hits[i].y_up;
            TMarker *mt = new TMarker(x, y, 20); mt->SetMarkerColor(kRed); mt->SetMarkerSize(1.0); mt->Draw("same");
        }
        
        TLatex latex; latex.SetNDC(); latex.SetTextSize(0.035); latex.SetTextAlign(13);
        double x_intercept_true = (std::abs(a_true) > 1e-9) ? -b_true / a_true : 0.0;
        double x_intercept_reco = (std::abs(best_a) > 1e-9) ? -best_b / best_a : 0.0;
        double diff_val = x_intercept_reco - x_intercept_true;

        latex.SetTextColor(kBlue); latex.DrawLatex(0.15, 0.85, Form("True: y = %.4fx + %.4f", a_true, b_true));
        latex.DrawLatex(0.15, 0.81, Form("      x(y=0) = %.4f", x_intercept_true)); 
        latex.SetTextColor(kRed); latex.DrawLatex(0.15, 0.75, Form("Reco: y = %.4fx + %.4f", best_a, best_b));
        latex.DrawLatex(0.15, 0.71, Form("      x(y=0) = %.4f", x_intercept_reco)); 
        latex.SetTextColor(kBlack); latex.DrawLatex(0.15, 0.65, Form("Diff: x_reco - x_true = %.4f cm", diff_val));

        cTrk.SaveAs(Form("%s/track_reco_combinatorial_ev%d.png", figDir.c_str(), ev));

        t_event_id = ev;

        t_a_true = a_true;
        t_b_true = b_true;
        t_x_true = x_intercept_true;

        t_a_reco = best_a;
        t_b_reco = best_b;
        t_x_reco = x_intercept_reco;

        t_diff = diff_val;
        t_chi2 = min_chi2;
        t_nHits = nHits;

        t_theta_true = true_angle_rad;
        t_theta_reco = std::atan(best_a);
        t_chi2_ndf = min_chi2 / std::max(1, nHits - 2);

        t_active_x0 = activeTrack.x0;
        t_active_y0 = activeTrack.y0;
        t_active_x1 = activeTrack.x1;
        t_active_y1 = activeTrack.y1;
        t_nClusters = nClusters;

        double sum_t_hit = 0.0;
        double sum_L_meas = 0.0;
        double sum_gain = 0.0;
        double sum_res2 = 0.0;
        double max_abs_res = 0.0;
        int n_side_wrong = 0;

        for (int i = 0; i < nHits; ++i) {
          const double x = hits[i].xwire;

          // Candidate selected by the best mask.
          const bool use_down = ((best_mask >> i) & 1);
          const double y_selected = use_down ? hits[i].y_dn : hits[i].y_up;

          const double y_reco = best_a * x + best_b;
          const double y_true = a_true * x + b_true;

          const double hit_residual = y_selected - y_reco;
          const double truth_residual = y_selected - y_true;

          // Side definition: +1 = up, -1 = down.
          const int selected_side = use_down ? -1 : +1;
          const int true_side = (y_true >= hits[i].ywire) ? +1 : -1;
          const int side_ok = (selected_side == true_side) ? 1 : 0;

          if (!side_ok) n_side_wrong++;

          sum_t_hit += hits[i].t_hit;
          sum_L_meas += hits[i].L_meas;
          sum_gain += hits[i].gain;
          sum_res2 += hit_residual * hit_residual;
          max_abs_res = std::max(max_abs_res, std::abs(hit_residual));

          // Fill hit-level tree.
          h_event_id = ev;
          h_nHits = nHits;
          h_hit_id = i;
          h_wire_index = hits[i].wire_index;
          h_selected_side = selected_side;
          h_true_side = true_side;
          h_side_ok = side_ok;

          h_diff_x = diff_val;
          h_chi2 = min_chi2;
          h_wire_x = hits[i].xwire;
          h_wire_y = hits[i].ywire;
          h_t_hit = hits[i].t_hit;
          h_L_meas = hits[i].L_meas;
          h_gain = hits[i].gain;
          h_n_primary = hits[i].n_primary;
          h_n_total = hits[i].n_total;
          h_y_selected = y_selected;
          h_y_reco = y_reco;
          h_y_true = y_true;
          h_hit_residual = hit_residual;
          h_truth_residual = truth_residual;

          hitTree->Fill();
        }

        t_mean_t_hit = sum_t_hit / nHits;
        t_mean_L_meas = sum_L_meas / nHits;
        t_mean_gain = sum_gain / nHits;
        t_rms_hit_residual = std::sqrt(sum_res2 / nHits);
        t_max_abs_hit_residual = max_abs_res;
        t_n_side_wrong = n_side_wrong;

        tree->Fill();
    } 

    fOut->cd();
    tree->Write();
    hitTree->Write();
    fOut->Close();
#if SAVE_DRIFT_ROOT
    fDriftOut->cd(); tDrift->Write(); fDriftOut->Close();
#endif
#if SAVE_GAIN_ROOT
    fGainOut->cd(); tGain->Write(); fGainOut->Close();
#endif
#if SAVE_VERIFY_ROOT
    fVerifyOut->cd(); tVerify->Write(); fVerifyOut->Close();
#endif

  return 0;
}