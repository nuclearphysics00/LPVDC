// fieldview_array.cc
// [修正版: Finite(個別配置)ジオメトリ対応]
// - 周期境界条件(Periodic)をOFFにしたジオメトリに対応
// - DrawElectrodes と NearestWireSurface を個別ワイヤー対応に改修
// - ViewField::SetLineColor() を削除
// - drift_all / drift_primary 
// - アノード信号のみ読み出し

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <limits>
#include <random>
#include <cstdint>
#include <iostream>

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
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/FundamentalConstants.hh"  
#include "Garfield/ViewDrift.hh" // ★これを追加


#include "detector_geometry.hh"
#include "geometry/eval_uniform.hh"
#include "eval_maps.hh"

#if defined(GEOM_PLATE)
  #include "geometry/geom_pap_plate.hh"
  //#include "geometry/geom_ap_plate.hh"
#elif defined(GEOM_WIRE)
  #include "geometry/geom_pap_wirecath.hh"
  //#include "geometry/geom_ap_wirecath.hh"

#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

#define DEBUG_WAVEFORM

using namespace Garfield;

// ==================== Parameters ====================
const double SIGNAL_THRESHOLD = 0.005; 

double TransferFunction(double t) {
  constexpr double tau = 20.0; 
  constexpr double gain = 100.0; 
  if (t < 0.0) return 0.0;
  return gain * (t / tau) * std::exp(1.0 - t / tau);
}

// ==================== Helpers ====================
static void DrawElectrodes(const Detector::Geometry& g) {
  const double xmin = g.xmin, xmax = g.xmax;
  for (const auto& e : g.electrodes) {
    if (e.kind == Detector::ElectrodeKind::PlaneY) {
      auto ln = new TLine(xmin, e.y, xmax, e.y);
      ln->SetLineColor(kGray + 3); ln->SetLineStyle(2); ln->SetLineWidth(3);
      ln->Draw("same");
      continue;
    }
    
    auto lower = [](std::string s){
      std::transform(s.begin(), s.end(), s.begin(), ::tolower); return s;
    };
    std::string nm = lower(e.name);
    const bool isAnode = nm.find("anode") != std::string::npos;
    const bool isPW    = nm.find("pw")    != std::string::npos;

    // ★ 周期設定がOFFならコピーを展開せず、自分自身だけを描画
    if (!g.periodicX) {
      const double x = e.x0; 
      if (x < xmin || x > xmax) continue;
      auto cir = new TEllipse(x, e.y, e.radius, e.radius);
      cir->SetFillStyle(1001);
      if      (isAnode) cir->SetFillColor(kOrange + 7);
      else if (isPW)    cir->SetFillColor(kAzure + 2);
      else              cir->SetFillColor(kBlack);
      
      cir->SetLineColor(kBlack);
      cir->SetLineWidth(1);
      cir->Draw("same");
      continue;
    }

    // (互換用) 周期ON時の描画
    const double pitch = g.pitchX;
    int nL = int(std::floor((xmin - e.x0) / pitch)) - 1;
    int nR = int(std::ceil ((xmax - e.x0) / pitch)) + 1;
    for (int n = nL; n <= nR; ++n) {
      const double x = e.x0 + n * pitch; if (x < xmin || x > xmax) continue;
      auto cir = new TEllipse(x, e.y, e.radius, e.radius);
      cir->SetFillStyle(1001);
      
      if      (isAnode) cir->SetFillColor(kOrange + 7);
      else if (isPW)    cir->SetFillColor(kAzure + 2);
      else              cir->SetFillColor(kBlack);
      
      cir->SetLineColor(kBlack);
      cir->SetLineWidth(1);
      cir->Draw("same");
    }
  }
}

static double DetectGap(const Detector::Geometry& g) {
  double gap = 0.0;
  for (const auto& e : g.electrodes) gap = std::max(gap, std::abs(e.y));
  if (gap <= 0.0) gap = std::min(std::abs(g.ymax), std::abs(g.ymin));
  return gap;
}

static void AddSeedsUniformOnLine(const Detector::Geometry& g, double x0, double y0, double ux, double uy, double spacing, std::vector<double>& xs, std::vector<double>& ys, std::vector<double>& zs, double margin=1.0e-3) {
  const double L = std::hypot(ux, uy); if (L == 0) return;
  ux /= L; uy /= L;
  const double xmin = g.xmin + margin, xmax = g.xmax - margin;
  const double ymin = g.ymin + margin, ymax = g.ymax - margin;
  for(double s=-1.0; s<1.0; s+=spacing) {
      double px = x0+s*ux, py=y0+s*uy;
      if(px>xmin && px<xmax && py>ymin && py<ymax) { 
          xs.push_back(px); ys.push_back(py); zs.push_back(0); 
      }
  }
}

struct WireNearest {
  double d_surface_cm;
  double r_cm;
  double xw_cm, yw_cm;
  std::string name;
  bool   hasWire;
};

// ★ 周期設定OFFに対応した最寄りワイヤー探索
static WireNearest NearestWireSurface(const Detector::Geometry& g, double x, double y){
  WireNearest out;
  out.d_surface_cm = 1e300; out.r_cm = 0.0; out.xw_cm = 0.0; out.yw_cm = 0.0; out.name = ""; out.hasWire = false;
  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    
    double xw = e.x0;
    if (g.periodicX) {
        const double k  = std::round((x - e.x0) / g.pitchX);
        xw = e.x0 + k * g.pitchX;
    }
    const double d  = std::hypot(x - xw, y - e.y) - e.radius;
    if (d < out.d_surface_cm) {
      out.d_surface_cm = d; out.r_cm = e.radius; out.xw_cm = xw; out.yw_cm = e.y; out.name = e.name; out.hasWire = true;
    }
  }
  return out;
}

// ---- main ----------------------------------------------------
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();

  // --- env ---
  const char* p_shards = std::getenv("N_SHARDS");
  const char* p_shard  = std::getenv("PBS_ARRAY_INDEX");
  const char* p_shard2 = std::getenv("SHARD_ID");
  const char* p_out    = std::getenv("OUT_ROOT");
  const char* p_maps   = std::getenv("MAKE_MAPS");
  const char* p_halo   = std::getenv("HALO_MULT");   
  const char* p_nevt   = std::getenv("N_EVENTS");    
  const char* p_seed   = std::getenv("SEED_BASE");   

  const int N_SHARDS = p_shards ? std::max(1, std::atoi(p_shards)) : 1;
  int SHARD_ID = -1;
  if (p_shard && std::strlen(p_shard)) SHARD_ID = std::atoi(p_shard);
  else if (p_shard2 && std::strlen(p_shard2)) SHARD_ID = std::atoi(p_shard2);
  else SHARD_ID = 0; 

  const double HALO_MULT = (p_halo && std::strlen(p_halo)) ? std::atof(p_halo) : 3.0;
  const long long N_EVENTS_TOTAL = p_nevt ? std::atoll(p_nevt) : 10000LL;

  unsigned long long baseSeed = p_seed ? std::strtoull(p_seed, nullptr, 10) : 0x9e3779b97f4a7c15ULL;
  unsigned long long seed = baseSeed + (unsigned long long)SHARD_ID * 1000003ULL;
  std::mt19937_64 rng(seed);

  char outname[512];
  if (p_out && std::strlen(p_out)) std::snprintf(outname, sizeof(outname), "%s", p_out);
  else std::snprintf(outname, sizeof(outname), "grid_times.shard%04d.root", SHARD_ID);

  std::string outpath = outname;
  auto slashPos = outpath.find_last_of('/');
  std::string baseDir = (slashPos == std::string::npos) ? std::string(".") : outpath.substr(0, slashPos);
  
#if defined(GEOM_WIRE)
  const char* geomTag = "wirecath";
#elif defined(GEOM_PLATE)
  const char* geomTag = "plate";
#endif
  std::string figDir = baseDir + "/png_" + std::string(geomTag);
  const char* p_fig = std::getenv("FIG_DIR");
  if (p_fig && std::strlen(p_fig)) figDir = p_fig;
  gSystem->mkdir(figDir.c_str(), kTRUE);

  auto saveFig = [&](TCanvas& c, const char* fname){ c.SaveAs((figDir + "/" + fname).c_str()); };
  const bool make_maps = (p_maps && std::atoi(p_maps) != 0) && (SHARD_ID == 0);

  std::printf("[env] SHARD_ID=%d  OUT=%s  MAKE_MAPS=%d\n", SHARD_ID, outname, (int)make_maps);

  // --- Gas & Geometry ---
  MediumMagboltz gas;
  gas.LoadGasFile("ic4H10_100_0.1atm.gas");
  // gas.LoadGasFile("ic4H10_100_1atm.gas");

  std::cout << "[Gas Info] Pressure: " << gas.GetPressure() << " Torr ("
            << gas.GetPressure() / 760.0 << " atm)" << std::endl;

  Detector::Geometry geo;

  // ===== 共通設定 =====
  const double gap           = 0.50;     // [cm]
  const double vAn           = 0.0;      // [V]
  const double vCat          = -1100.0;  // [V]
  const bool   pwSameAsAnode = false;
  const int    nWires        = 4;

#if defined(GEOM_PLATE)

  // ===== Plane cathode 用 =====
  const double pitchSense = 0.60;      // [cm]
  const double phase      = 0.0;       // [cm]
  const double rAn        = 1.0e-3;    // [cm]
  const double rPW        = 2.5e-3;    // [cm]
  const double vPW        = vCat * (55.0 / 700.0);

  std::puts("[build] geometry = PAP_PlaneCathode_Periodic");
  std::printf("[build] pitchSense=%.3f cm gap=%.3f cm vAn=%.1f V vPW=%.3f V vCat=%.1f V\n",
              pitchSense, gap, vAn, vPW, vCat);

  geo = Detector::PAP_PlaneCathode_Periodic(
      pitchSense,   // pitchSense
      gap,          // gap
      phase,        // phase
      rAn,          // rAn
      rPW,          // rPW
      vAn,          // vAn
      vPW,          // vPW
      vCat,         // vCat
      pwSameAsAnode,
      nWires
  );

#elif defined(GEOM_WIRE)

  // ===== Wire cathode 用 =====
  const double pitchSense = 0.60;      // [cm]
  const double pwOffset   = 0.0;       // [cm]
  const double pitchCath  = 0.20;      // [cm]
  const double cathPhase  = 0.10;      // [cm]
  const double rAn        = 1.0e-3;    // [cm]
  const double rPW        = 2.5e-3;    // [cm]
  const double rCat       = 2.5e-3;    // [cm]
  const double vPW        = vCat * (50.0 / 700.0);

  std::puts("[build] geometry = PAP_WireCathode_Periodic");
  std::printf("[build] pitchSense=%.3f cm pitchCath=%.3f cm gap=%.3f cm vAn=%.1f V vPW=%.3f V vCat=%.1f V\n",
              pitchSense, pitchCath, gap, vAn, vPW, vCat);

  geo = Detector::PAP_WireCathode_Periodic(
      pitchSense,   // pitchSense
      gap,          // gap
      pwOffset,     // pwOffset
      pitchCath,    // pitchCath
      cathPhase,    // cathPhase
      rAn,          // rAn
      rPW,          // rPW
      rCat,         // rCat
      vAn,          // vAn
      vPW,          // vPW
      vCat,         // vCat
      pwSameAsAnode,
      true,         // printDebug
      nWires
  );

#else
  #error "Define one of: -DGEOM_PLATE or -DGEOM_WIRE"
#endif

  auto comp = Detector::BuildField(geo);
  comp->SetMedium(&gas);
  
  // --- Sensor ---
  Sensor sensorRnd;
  sensorRnd.AddComponent(comp.get());
  sensorRnd.SetArea(geo.xmin, geo.ymin, -0.1, geo.xmax, geo.ymax, +0.1);
  std::vector<std::string> anode_names; // アノード名を保存するリスト
  std::cout << "[Info] Adding electrodes to Sensor (Filtering 'Anode' only):" << std::endl;
  bool anodeFound = false;
  for(const auto& el : geo.electrodes) {
      if (el.kind == Detector::ElectrodeKind::WireRow) {
          std::string nm = el.name;
          std::string nm_lower = nm;
          std::transform(nm_lower.begin(), nm_lower.end(), nm_lower.begin(), ::tolower);
          
          if (nm_lower.find("anode") != std::string::npos) {
              std::cout << "  -> Added: " << el.name << std::endl;
              sensorRnd.AddElectrode(comp.get(), el.name); 
              anode_names.push_back(el.name); // 見つけたアノード名をリストに記録
              anodeFound = true;
          }
      }
  }
  if(!anodeFound) std::cerr << "[Warning] No 'Anode' found in geometry!" << std::endl;

  const double tMin = 0., tMax = 500., tStep = 0.05;
  const int nTimeBins = (int)((tMax - tMin) / tStep);
  sensorRnd.SetTimeWindow(tMin, tStep, nTimeBins);
  sensorRnd.SetTransferFunction(TransferFunction);

  // ===== can_drift() =====
  struct DriftCheck { bool driftable; bool in_halo; double d_wire_cm; double r_wire_cm; };
  auto can_drift = [&](double x, double y)->DriftCheck{
    double Ex=0, Ey=0, Ez=0, V=0; Medium* m=nullptr; int status=0;
    comp->ElectricField(x, y, 0.0, Ex, Ey, Ez, V, m, status);
    const bool in_area = sensorRnd.IsInArea(x, y, 0.0);
    WireNearest wn = NearestWireSurface(geo, x, y);
    double dpos = wn.d_surface_cm; if (!std::isfinite(dpos)) dpos = 1e300;
    const bool in_halo = (wn.hasWire && HALO_MULT > 0.0 && dpos >= 0.0 && dpos < HALO_MULT * wn.r_cm);
    const bool driftable = (m != nullptr) && (status == 0) && in_area;
    return {driftable, in_halo, dpos, wn.r_cm};
  };

  // =========================================================
  // (A) Map Gen (Field Lines & Maps)
  // =========================================================
  if (make_maps) {
    std::printf("[map] Generating field maps & lines...\n");

    const int NX=200, NY=160;
    TH2D hE("|E|","|E|;X [cm];Y [cm]", NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    TH2D hV("V","V;X [cm];Y [cm]",     NX, geo.xmin, geo.xmax, NY, geo.ymin, geo.ymax);
    double minPos=1e98, maxVal=0.0;
    for (int ix=1; ix<=NX; ++ix) {
      const double x = hE.GetXaxis()->GetBinCenter(ix);
      for (int iy=1; iy<=NY; ++iy) {
        const double y = hE.GetYaxis()->GetBinCenter(iy);
        double ex=0,ey=0,ez=0,V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(x,y,0,ex,ey,ez,V,m,status);
        const double E = std::hypot(ex,std::hypot(ey,ez));
        hE.SetBinContent(ix,iy,E);
        hV.SetBinContent(ix,iy,V);
        if (E>0 && E<minPos) minPos=E;
        if (E>maxVal) maxVal=E;
      }
    }

    {
      TCanvas cE("cE","|E| map",900,720); cE.SetRightMargin(0.14); cE.SetGrid(); cE.SetLogz();
      if (minPos < 1e97) hE.SetMinimum(minPos * 0.9);
      hE.SetMaximum(maxVal * 1.05); hE.Draw("COLZ");
      saveFig(cE, "emap_absE.png");
    }

    TCanvas c("cField","Combo Map",1000,780); c.SetRightMargin(0.14); c.SetGrid(); c.SetLogz();
    if (minPos<1e97) hE.SetMinimum(minPos*0.9); hE.SetMaximum(maxVal*1.05);
    hE.Draw("COLZ");
    hV.SetContour(40); hV.SetLineColor(kGray+2); hV.SetLineWidth(2); hV.Draw("CONT3 SAME");
    ViewField vf; vf.SetCanvas(&c); vf.SetComponent(comp.get()); vf.SetArea(geo.xmin, geo.ymin, geo.xmax, geo.ymax);
    std::vector<double> xs,ys,zs;
    const double th = +45.0 * M_PI / 180.0;
    AddSeedsUniformOnLine(geo, 0.0, 0.0, std::cos(th), std::sin(th), 0.05, xs,ys,zs);
    if (!xs.empty()) vf.PlotFieldLines(xs,ys,zs,true,true);
    DrawElectrodes(geo);
    saveFig(c, "fieldmap_combo.png");

    {
      std::printf("[map] Drawing pure field lines (no track)...\n");
      const int    FL_N    = 64;   
      const double yfrac   = 0.98; 
      const bool   both    = true; 
      const double edge_eps2 = 1e-3;
      const double gapv     = DetectGap(geo);
      const double yTop     = +yfrac * gapv;
      const double yBot     = -yfrac * gapv;

      auto inMedium = [&](double x, double y){
        double ex=0,ey=0,ez=0,V=0; Medium* m=nullptr; int status=0;
        comp->ElectricField(x,y,0.0,ex,ey,ez,V,m,status);
        return (m!=nullptr && status==0);
      };

      auto make_line_seeds = [&](double yseed, std::vector<double>& xs2,
                                 std::vector<double>& ys2, std::vector<double>& zs2){
        const double Lx = (geo.xmax - edge_eps2) - (geo.xmin + edge_eps2);
        for (int i = 0; i < FL_N; ++i) {
          const double x = geo.xmin + edge_eps2 + (i + 0.5) * Lx / FL_N;
          if (!inMedium(x, yseed)) continue;
          xs2.push_back(x);
          ys2.push_back(yseed);
          zs2.push_back(0.0);
        }
      };

      TCanvas cFL("cFieldLines","Field lines from y = +/- yfrac*gap", 1000, 780);
      cFL.SetGrid();
      TH2D frame("fl_frame","Field lines; x [cm]; y [cm]", 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
      frame.SetStats(0); frame.Draw();

      ViewField vfFL; vfFL.SetCanvas(&cFL); vfFL.SetComponent(comp.get());
      vfFL.SetArea(geo.xmin, geo.ymin, geo.xmax, geo.ymax);

      std::vector<double> xs2, ys2, zs2;
      make_line_seeds(yTop, xs2, ys2, zs2);
      if (both) make_line_seeds(yBot, xs2, ys2, zs2);

      if (!xs2.empty()) {
          vfFL.PlotFieldLines(xs2, ys2, zs2, true, true);
      }
      DrawElectrodes(geo);

      char fname[128];
      std::snprintf(fname, sizeof(fname), "fieldlines_yfrac_%03d%s.png", int(std::round(yfrac*100)), both? "_both":"");
      saveFig(cFL, fname);
      std::printf("[map] Saved %s\n", fname);
    }
  }

  // ===== (B) Loop (Drift Simulation) =====
  const double edge_eps = 1e-3;
  const double gap_val  = DetectGap(geo);
  const double xminU = geo.xmin + edge_eps;
  const double xmaxU = geo.xmax - edge_eps;
  const double yminU = -gap_val + edge_eps;
  const double ymaxU = +gap_val - edge_eps;

  auto inGas = [&](double x, double y){
    double ex,ey,ez,V; Medium* m=nullptr; int status=0;
    comp->ElectricField(x,y,0.0,ex,ey,ez,V,m,status);
    return (m!=nullptr && status==0);
  };
  auto nudge_to_gas = [&](double& x, double& y)->bool{
    if (inGas(x,y)) return true;
    const double dx = (xmaxU - xminU), dy = (ymaxU - yminU);
    const double ox[9] = {0, +0.25*dx, -0.25*dx, 0, 0, +0.25*dx, +0.25*dx, -0.25*dx, -0.25*dx};
    const double oy[9] = {0, 0, 0, +0.25*dy, -0.25*dy, +0.25*dy, -0.25*dy, +0.25*dy, -0.25*dy};
    for (int k=1;k<9;++k){
      const double xx = x + ox[k], yy = y + oy[k];
      if (xx<=xminU || xx>=xmaxU || yy<=yminU || yy>=ymaxU) continue;
      if (inGas(xx,yy)) { x=xx; y=yy; return true; }
    }
    return false;
  };

  std::uniform_real_distribution<double> distX(xminU, xmaxU);
  std::uniform_real_distribution<double> distY(yminU, ymaxU);

  AvalancheMC amcRnd;
  amcRnd.SetSensor(&sensorRnd);
  amcRnd.UseWeightingPotential(true); 
  amcRnd.EnableSignalCalculation(); 

  struct RowData {
    double xwire;
    double t0_ns; 
    double t_th;  
    int    in_halo;
    double dwire_cm;
  };
  struct RowAll { 
    double xwire;
    double t_ns;
    int    in_halo;
    double dwire_cm;
  };

  std::vector<RowData> rows_prim;
  std::vector<RowAll>  rows_all;
  rows_prim.reserve(N_EVENTS_TOTAL / N_SHARDS + 1000);
  rows_all.reserve(N_EVENTS_TOTAL / N_SHARDS * 2);

  long long visited_evt = 0;

  for (long long iev = 0; iev < N_EVENTS_TOTAL; ++iev) {
    if ((iev % N_SHARDS) != SHARD_ID) continue;
    ++visited_evt;

    double x = distX(rng);
    double y = distY(rng);
    if (!nudge_to_gas(x, y)) continue;

    auto chk = can_drift(x, y);
    if (!chk.driftable) continue;

    sensorRnd.ClearSignal();

    // --- ★ Debug Plot: ViewDrift setup ★ ---
    #ifdef DEBUG_WAVEFORM
        ViewDrift driftView;
        if (visited_evt == 1 && SHARD_ID == 0) {
            amcRnd.EnablePlotting(&driftView);
        }
    #endif

    amcRnd.DriftElectron(x, y, 0.0, 0.0);
    // --- ★ Debug Plot: Draw Avalanche Tracks ★ ---
    #ifdef DEBUG_WAVEFORM
        if (visited_evt == 1  && SHARD_ID == 0) {
            TCanvas cClEnd("cClEnd", "cluster endpoint tracks", 1100, 850);
            cClEnd.SetGrid();
            TH2D frClEnd("frClEnd", "Electron Drift Lines;X [cm];Y [cm]", 10, geo.xmin, geo.xmax, 10, geo.ymin, geo.ymax);
            frClEnd.SetStats(0); 
            frClEnd.Draw();
            
            // 電極を描画
            DrawElectrodes(geo);
            
            // ドリフト軌跡を描画
            driftView.SetCanvas(&cClEnd);
            driftView.Plot(true, false);
            
            TMarker mStart(x, y, 20); 
            mStart.SetMarkerColor(kBlue);
            mStart.SetMarkerSize(1.5);
            mStart.Draw("same");
            
            char fname[128];
            std::snprintf(fname, sizeof(fname), "cluster_endtracks_ev%lld.png", iev); 
            saveFig(cClEnd, fname);
        }
    #endif
    
    double t_exact_first = 1e300;
    bool hit_anode = false;
    double end_x = 0, end_y = 0;
    std::string hit_wire_name = "";

    WireNearest wn_start = NearestWireSurface(geo, x, y);

    const size_t nEnd = amcRnd.GetNumberOfElectronEndpoints();
    for(size_t k=0; k<nEnd; ++k){
        double xa,ya,za,ta, xb,yb,zb,tb; int st;
        amcRnd.GetElectronEndpoint(k, xa,ya,za,ta, xb,yb,zb,tb, st);
        WireNearest wn_end = NearestWireSurface(geo, xb, yb);
        
        std::string nm = wn_end.name;
        std::transform(nm.begin(), nm.end(), nm.begin(), ::tolower);

        if (wn_end.hasWire && nm.find("anode") != std::string::npos) {
            hit_anode = true;
            if (tb < t_exact_first) t_exact_first = tb;
            end_x = xb; end_y = yb;
            hit_wire_name = wn_end.name;

            rows_all.push_back({
                wn_start.xw_cm,
                tb,
                chk.in_halo ? 1 : 0,
                wn_start.d_surface_cm
            });
        }
    }

    if (!hit_anode) continue; 

    sensorRnd.ConvoluteSignals();

    std::string readoutLabel = hit_wire_name; 
    double t_threshold_found = -1.0;
    
    for (int k = 0; k < nTimeBins; ++k) {
        double sig = sensorRnd.GetSignal(readoutLabel, k, true); 
        if (std::abs(sig) > SIGNAL_THRESHOLD) {
            t_threshold_found = tMin + k * tStep;
            break;
        }
    }

    // --- ★ Debug Waveform Plot ★ ---
    #ifdef DEBUG_WAVEFORM
        if (visited_evt == 1 && SHARD_ID == 0) {
            for(const auto& a_name : anode_names) {
                std::string cName = "cWave_" + a_name;
                TCanvas cWave(cName.c_str(), "", 600, 400);
                cWave.SetGrid();

                std::vector<double> tx(nTimeBins), vy(nTimeBins);
                for(int k=0; k<nTimeBins; ++k) {
                    tx[k] = tMin + k*tStep;
                    vy[k] = sensorRnd.GetSignal(a_name, k, true);
                }

                TGraph gW(nTimeBins, tx.data(), vy.data());
                std::string title = a_name + " (Ev " + std::to_string(iev) + ");Time [ns];Signal";
                gW.SetTitle(title.c_str());
                gW.SetLineWidth(1); 
                gW.SetLineColor(kBlue);
                gW.Draw("AL");

                TLine lth_p(tMin, SIGNAL_THRESHOLD, tMax, SIGNAL_THRESHOLD);
                lth_p.SetLineColor(kRed); lth_p.SetLineStyle(2); lth_p.Draw("same");

                TLine lth_n(tMin, -SIGNAL_THRESHOLD, tMax, -SIGNAL_THRESHOLD);
                lth_n.SetLineColor(kRed); lth_n.SetLineStyle(2); lth_n.Draw("same");

                char wname[128]; 
                std::snprintf(wname, sizeof(wname), "wire_sig_ev%lld_%s.png", iev, a_name.c_str());
                
                saveFig(cWave, wname);
            }
        }
    #endif    

    rows_prim.push_back({
        wn_start.xw_cm,
        t_exact_first,
        t_threshold_found,
        chk.in_halo ? 1 : 0,
        wn_start.d_surface_cm
    });
  }

  // ===== ROOT output =====
  TFile fout(outname, "RECREATE");
  fout.SetCompressionSettings(207);
  
  TTree t_prim("drift_primary","primary-like first-arrival times");
  int    br_shard = SHARD_ID, br_nshards = N_SHARDS;
  double br_xwire = 0.0, br_t0ns = 0.0, br_t_th = 0.0, br_dwire = 0.0;
  int    br_inhalo = 0;

  t_prim.Branch("shard_id", &br_shard,    "shard_id/I");
  t_prim.Branch("n_shards", &br_nshards,  "n_shards/I");
  t_prim.Branch("xwire",    &br_xwire,    "xwire/D");
  t_prim.Branch("t0_ns",    &br_t0ns,     "t0_ns/D");
  t_prim.Branch("t_th",     &br_t_th,     "t_th/D");
  t_prim.Branch("in_halo",  &br_inhalo,   "in_halo/I");
  t_prim.Branch("dwire_cm", &br_dwire,    "dwire_cm/D");

  for (const auto& r : rows_prim) {
    br_xwire  = r.xwire;
    br_t0ns   = r.t0_ns; 
    br_t_th   = r.t_th;  
    br_inhalo = r.in_halo;
    br_dwire  = r.dwire_cm;
    t_prim.Fill();
  }

  TTree t_all("drift_all","all electron arrival times");
  double br_tns_all = 0.0;
  
  t_all.Branch("shard_id", &br_shard,   "shard_id/I");
  t_all.Branch("xwire",    &br_xwire,   "xwire/D");
  t_all.Branch("t_ns",     &br_tns_all, "t_ns/D");
  t_all.Branch("in_halo",  &br_inhalo,  "in_halo/I");
  t_all.Branch("dwire_cm", &br_dwire,   "dwire_cm/D");

  for (const auto& r : rows_all) {
      br_xwire   = r.xwire;
      br_tns_all = r.t_ns;
      br_inhalo  = r.in_halo;
      br_dwire   = r.dwire_cm;
      t_all.Fill();
  }
  
  TTree meta("meta","run metadata");
  double meta_xmin=xminU, meta_xmax=xmaxU, meta_ymin=yminU, meta_ymax=ymaxU;
  meta.Branch("xmin",&meta_xmin,"xmin/D");
  meta.Branch("xmax",&meta_xmax,"xmax/D");
  meta.Branch("ymin",&meta_ymin,"ymin/D");
  meta.Branch("ymax",&meta_ymax,"ymax/D");
  meta.Fill();

  t_prim.Write();
  t_all.Write();
  meta.Write(); 
  fout.Close();
  
  std::printf("[done] Processed %lld events.\n", visited_evt);
  std::printf("       drift_primary: %zu rows\n", rows_prim.size());
  std::printf("       drift_all    : %zu rows\n", rows_all.size());
  std::printf("       Written to %s\n", outname);

  return 0;
}