// geometry/eval_uniform.hh
#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "Garfield/ComponentAnalyticField.hh"
#include "detector_geometry.hh"

namespace Eval {

inline void prepare_logz(TH2* h, double min_floor = 1e-12) {
  if (!h) return;
  const int nx = h->GetNbinsX();
  const int ny = h->GetNbinsY();
  double minPos = std::numeric_limits<double>::infinity();
  double maxVal = 0.0;

  for (int ix = 1; ix <= nx; ++ix) {
    for (int iy = 1; iy <= ny; ++iy) {
      const double z = h->GetBinContent(ix, iy);
      if (z > 0.0 && z < minPos) minPos = z;
      if (z > maxVal) maxVal = z;
    }
  }
  if (!std::isfinite(minPos)) minPos = min_floor;       // 全部0等の場合の保険
  if (minPos <= 0.0) minPos = min_floor;                // logZでは >0 が必須
  if (maxVal > 0.0 && maxVal <= minPos) maxVal = minPos * 1.01;

  h->SetMinimum(minPos);
  if (maxVal > 0.0) h->SetMaximum(maxVal);
  // 表示面の微調整（任意）
  h->GetZaxis()->SetMoreLogLabels();
  h->GetZaxis()->SetNoExponent();
}

// ---------------- 共通：ワイヤ近傍判定（周期境界対応） ----------------
inline bool NearAnyWire(double x, double y, const Detector::Geometry& g, double multR=5.0){
  for (const auto& e : g.electrodes) {
    if (e.kind != Detector::ElectrodeKind::WireRow) continue;
    double xw = e.x0;
    if (g.periodicX) {
      const double k = std::round((x - e.x0)/g.pitchX);
      xw = e.x0 + k * g.pitchX;
    }
    const double r = std::max(1e-4, e.radius);
    if (std::hypot(x - xw, y - e.y) < multR * r) return true;
  }
  return false;
}

// ---------------- 基準：アノード–カソード（上/下） -------------------
struct AC_Baseline {
  double yAn=0, vAn=0;
  double yTop=0, vTop=0, yBot=0, vBot=0;
  double E0_top=0, E0_bot=0;  // [V/cm]
  bool   ok=false;
};

// アノード優先規則：名前に "anode" を含む電極。無ければ最高電位を採用。
inline AC_Baseline MakeACBaseline(const Detector::Geometry& g) {
  AC_Baseline b{};
  // anode 探索
  bool foundAn = false; double vMax=-1e300; size_t iMax=0;
  for (size_t i=0;i<g.electrodes.size();++i){
    const auto& e = g.electrodes[i];
    std::string name=e.name; std::transform(name.begin(),name.end(),name.begin(),::tolower);
    if (name.find("anode")!=std::string::npos){ b.yAn=e.y; b.vAn=e.V; foundAn=true; break; }
    if (e.V>vMax){ vMax=e.V; iMax=i; }
  }
  if (!foundAn && !g.electrodes.empty()){ const auto& e=g.electrodes[iMax]; b.yAn=e.y; b.vAn=e.V; }

  // 上下端
  b.yTop=-1e99; b.yBot=+1e99;
  for (const auto& e: g.electrodes){
    if (e.y > b.yTop){ b.yTop=e.y; b.vTop=e.V; }
    if (e.y < b.yBot){ b.yBot=e.y; b.vBot=e.V; }
  }
  const double dTop = std::abs(b.yTop - b.yAn);
  const double dBot = std::abs(b.yAn - b.yBot);
  if (dTop>0) b.E0_top = std::abs(b.vAn - b.vTop)/dTop;
  if (dBot>0) b.E0_bot = std::abs(b.vAn - b.vBot)/dBot;
  b.ok = (b.E0_top>0 || b.E0_bot>0);
  return b;
}

// ---------------- 集計構造体 -------------------
struct UniformStats {
  double E0{};  // 代表値（上/下平均の目安）
  // Ey 基準
  double meanEy{}, rmsRelEy{}, p95RelEy{}, maxRelEy{};
  double frac_lt1p{}, frac_lt2p{}, frac_lt5p{};
  // |E| 基準（参考）
  double meanE{},  rmsRelE{},  p95RelE{},  maxRelE{};
  // 1ピッチ横方向リップル（y=y0）
  double ripple_pp_Ey_over_E0{}, ripple_fundA_over_E0{};
};

// ---------------- 全域評価（端は切らない想定OK） -------------------
// 旧API名を維持：内部を AC 基準に変更
inline UniformStats EvaluateUniformityMasked(
    Garfield::ComponentAnalyticField& comp,
    const Detector::Geometry& g,
    double roi_margin_frac,            // 端を切らないなら 0.0
    double wire_exclude_mult,          // 0.0=除外なし、例:5.0=半径の5倍除外
    int nx=180, int ny=140,
    double y0_for_ripple=0.0) {

  const AC_Baseline B = MakeACBaseline(g);
  UniformStats s{};
  if (!B.ok) return s;

  // 代表 E0（情報用）：上/下の平均
  s.E0 = 0.5 * ( (B.E0_top>0?B.E0_top:0) + (B.E0_bot>0?B.E0_bot:0) );

  // ROI
  const double xm = g.xmin + roi_margin_frac*(g.xmax-g.xmin);
  const double xM = g.xmax - roi_margin_frac*(g.xmax-g.xmin);
  const double ym = g.ymin + roi_margin_frac*(g.ymax-g.ymin);
  const double yM = g.ymax - roi_margin_frac*(g.ymax-g.ymin);

  auto E0loc = [&](double y){ return (y >= B.yAn) ? B.E0_top : B.E0_bot; };

  std::vector<double> dEy, dE; dEy.reserve(nx*ny); dE.reserve(nx*ny);
  double sumEy=0, sumE=0; long nKeep=0;

  for (int ix=0; ix<nx; ++ix) {
    const double x = xm + (xM-xm) * (ix+0.5)/nx;
    for (int iy=0; iy<ny; ++iy) {
      const double y = ym + (yM-ym) * (iy+0.5)/ny;
      if (wire_exclude_mult>0.0 && NearAnyWire(x,y,g,wire_exclude_mult)) continue;

      const double E0y = E0loc(y);
      if (E0y<=0) continue;

      double ex=0,ey=0,ez=0,V=0; Garfield::Medium* m=nullptr; int st=0;
      comp.ElectricField(x,y,0,ex,ey,ez,V,m,st);
      const double E  = std::hypot(ex,std::hypot(ey,ez));
      // |Ey| と局所 E0 を比較
      const double d1 = std::abs(std::abs(ey) - E0y)/E0y;
      const double d2 = std::abs(E - E0y)/E0y;

      dEy.push_back(d1); dE.push_back(d2);
      sumEy += std::abs(ey); sumE += E; ++nKeep;
    }
  }

  auto percentile = [](std::vector<double>& v, double p){
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const double t = std::clamp(p,0.0,1.0)*(v.size()-1);
    size_t i=(size_t)t; double f=t-i;
    return (1.0-f)*v[i] + f*((i+1<v.size())? v[i+1] : v[i]);
  };
  auto rms = [](const std::vector<double>& v){
    if (v.empty()) return 0.0; double s2=0; for(double x:v) s2+=x*x; return std::sqrt(s2/v.size());
  };
  auto frac_lt = [](const std::vector<double>& v, double thr){
    if (v.empty()) return 0.0; long c=0; for(double x:v) if (x<thr) ++c; return double(c)/v.size();
  };

  s.meanEy=(nKeep?sumEy/nKeep:0.0); s.meanE=(nKeep?sumE/nKeep:0.0);
  s.rmsRelEy=rms(dEy); s.rmsRelE=rms(dE);
  s.p95RelEy=percentile(dEy,0.95); s.p95RelE=percentile(dE,0.95);
  s.maxRelEy=dEy.empty()?0:*std::max_element(dEy.begin(),dEy.end());
  s.maxRelE =dE .empty()?0:*std::max_element(dE .begin(),dE .end());
  s.frac_lt1p=frac_lt(dEy,0.01); s.frac_lt2p=frac_lt(dEy,0.02); s.frac_lt5p=frac_lt(dEy,0.05);

  // 1ピッチの横方向リップル（y=y0）：局所 E0(y0) で正規化
  const int NXl=600;
  std::vector<double> Ey; Ey.reserve(NXl);
  const double x0 = -0.5*g.pitchX, x1=+0.5*g.pitchX;
  const double E0y = E0loc(y0_for_ripple);
  for (int i=0;i<NXl;i++){
    const double x = x0 + (x1-x0)*i/(NXl-1);
    if (wire_exclude_mult>0.0 && NearAnyWire(x,y0_for_ripple,g,wire_exclude_mult)) continue;
    double ex=0,ey=0,ez=0,V=0; Garfield::Medium* m=nullptr; int st=0;
    comp.ElectricField(x,y0_for_ripple,0,ex,ey,ez,V,m,st);
    Ey.push_back(std::abs(ey));
  }
  if (Ey.size()>2 && E0y>0){
    const auto [mn,mx] = std::minmax_element(Ey.begin(),Ey.end());
    s.ripple_pp_Ey_over_E0 = (*mx - *mn)/E0y;
    double Scc=0, Scy=0;
    for (int i=0;i<(int)Ey.size();++i){
      const double x = x0 + (x1-x0)*i/(NXl-1);
      const double c = std::cos(2*M_PI*x/g.pitchX);
      Scc += c*c; Scy += c*(Ey[i]-E0y);
    }
    s.ripple_fundA_over_E0 = (Scc>0)? std::abs(Scy/Scc)/E0y : 0.0;
  }
  return s;
}

// ---------------- 端/中央を抽出してレポート（AC基準） --------------
struct RegionStats {
  const char* name{};
  double areaFrac{};
  double rmsRel{}, p95Rel{}, maxRel{};
  double frac_lt1p{}, frac_lt2p{}, frac_lt5p{};
};
struct EdgeReport {
  double E0{};  // 代表値（上/下平均）
  RegionStats bulk, top, bottom, left, right, wireHalo;
};

// 互換のため引数はそのまま（内部は AC 基準）。wireHaloRmult 等は同様に解釈。
inline RegionStats AccumulateRegion(Garfield::ComponentAnalyticField& comp,
                                    const Detector::Geometry& g, double /*E0_unused*/,
                                    double xmin,double xmax,double ymin,double ymax,
                                    int nx,int ny,
                                    double wireHaloRmult=0.0, bool invertHalo=false){
  const AC_Baseline B = MakeACBaseline(g);
  auto E0loc = [&](double y){ return (y >= B.yAn) ? B.E0_top : B.E0_bot; };
  auto inHalo=[&](double x,double y){
    const bool near=NearAnyWire(x,y,g,wireHaloRmult);
    return invertHalo? !near : near;
  };

  long n=0,nKeep=0,n1=0,n2=0,n5=0; double s2=0,p95=0,mx=0;
  std::vector<double> vals; vals.reserve(nx*ny);
  for(int ix=0;ix<nx;++ix){
    const double x=xmin + (xmax-xmin)*(ix+0.5)/nx;
    for(int iy=0;iy<ny;++iy){
      const double y=ymin + (ymax-ymin)*(iy+0.5)/ny;
      ++n;
      if (wireHaloRmult>0.0 && !inHalo(x,y)) continue;

      const double E0y = E0loc(y);
      if (E0y<=0) continue;

      double ex=0,ey=0,ez=0,V=0; Garfield::Medium* m=nullptr; int st=0;
      comp.ElectricField(x,y,0,ex,ey,ez,V,m,st);
      const double d = std::abs(std::abs(ey)-E0y)/E0y;

      vals.push_back(d); s2+=d*d; if (d>mx) mx=d;
      if (d<0.01) ++n1; if (d<0.02) ++n2; if (d<0.05) ++n5; ++nKeep;
    }
  }
  if (!vals.empty()){
    std::sort(vals.begin(),vals.end());
    const double t=0.95*(vals.size()-1);
    size_t i=(size_t)t; double f=t-i;
    p95=(1.0-f)*vals[i] + f*((i+1<vals.size())? vals[i+1]:vals[i]);
  }
  RegionStats r{};
  r.areaFrac=(n>0)? double(nKeep)/n : 0.0;
  r.rmsRel=(nKeep>0)? std::sqrt(s2/nKeep) : 0.0;
  r.p95Rel=p95; r.maxRel=mx;
  r.frac_lt1p=(nKeep>0)? double(n1)/nKeep : 0.0;
  r.frac_lt2p=(nKeep>0)? double(n2)/nKeep : 0.0;
  r.frac_lt5p=(nKeep>0)? double(n5)/nKeep : 0.0;
  return r;
}

inline EdgeReport EvaluateEdges(Garfield::ComponentAnalyticField& comp,
                                const Detector::Geometry& g,
                                double edgeBandFrac=0.10, int nx=220, int ny=160){
  const AC_Baseline B = MakeACBaseline(g);
  EdgeReport R{}; R.E0 = 0.5 * ( (B.E0_top>0?B.E0_top:0) + (B.E0_bot>0?B.E0_bot:0) );

  const double dx = edgeBandFrac*(g.xmax-g.xmin);
  const double dy = edgeBandFrac*(g.ymax-g.ymin);

  R.bulk   = AccumulateRegion(comp,g, /*E0_unused*/0.0, g.xmin+dx,g.xmax-dx,g.ymin+dy,g.ymax-dy, nx,ny, 5.0,true);  R.bulk.name="bulk";
  R.top    = AccumulateRegion(comp,g, /*E0_unused*/0.0, g.xmin,g.xmax,g.ymax-dy,g.ymax,           nx,ny, 0.0,false); R.top.name="top";
  R.bottom = AccumulateRegion(comp,g, /*E0_unused*/0.0, g.xmin,g.xmax,g.ymin,     g.ymin+dy,      nx,ny, 0.0,false); R.bottom.name="bottom";
  R.left   = AccumulateRegion(comp,g, /*E0_unused*/0.0, g.xmin,g.xmin+dx,g.ymin,  g.ymax,         nx,ny, 0.0,false); R.left.name="left";
  R.right  = AccumulateRegion(comp,g, /*E0_unused*/0.0, g.xmax-dx,g.xmax,g.ymin,  g.ymax,         nx,ny, 0.0,false); R.right.name="right";
  R.wireHalo=AccumulateRegion(comp,g, /*E0_unused*/0.0, g.xmin,g.xmax,g.ymin,     g.ymax,         nx,ny, 3.0,false); R.wireHalo.name="wire-halo(3r)";
  return R;
}

} // namespace Eval
