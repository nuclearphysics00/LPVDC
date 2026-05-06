// alpha_simple.cc
// Magboltz gas ファイルから α(E), v(E) と α(E/p), v(E/p)
// をプロットする実行ファイル（軸レンジの細工はしない簡略版）。
//  - ガス圧は 0.1 atm (= 76 Torr) に上書き
//  - E/p は atm ベース（V/(cm*atm)）
//  - ROOT の GUI で軸レンジなどを自由に調整できる
//  - ★α≈0 から立ち上がる E のしきい値もターミナルに表示する★

#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/FundamentalConstants.hh"

using namespace Garfield;

//--------------------------------------------------
// パッドの基本設定（グリッド + 余白 + log/lin）
//--------------------------------------------------
static void SetupPad(bool logx = false, bool logy = false) {
  if (logx) gPad->SetLogx();
  if (logy) gPad->SetLogy();
  gPad->SetGrid();
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.16);
  gPad->SetTopMargin(0.08);
}

//--------------------------------------------------
// 指定した E 配列から α, v, E/p(Torr) を計算
//--------------------------------------------------
static void SampleGas(MediumMagboltz& gas,
                      const std::vector<double>& E,
                      double P_torr,
                      std::vector<double>& Ep_torr,
                      std::vector<double>& alpha,
                      std::vector<double>& vdrift) {
  const int N = static_cast<int>(E.size());
  Ep_torr.resize(N);
  alpha.resize(N);
  vdrift.resize(N);

  for (int i = 0; i < N; ++i) {
    const double Ei = E[i];  // [V/cm]

    double vx = 0., vy = 0., vz = 0.;
    double a  = 0.;

    const bool okV = gas.ElectronVelocity(Ei, 0., 0., 0., 0., 0., vx, vy, vz);
    const bool okA = gas.ElectronTownsend(Ei, 0., 0., 0., 0., 0., a);

    const double vabs     = std::sqrt(vx*vx + vy*vy + vz*vz);      // [cm/ns]
    const double Ei_overP = (P_torr > 0.) ? (Ei / P_torr) : 0.;    // [V/(cm*Torr)]

    Ep_torr[i] = Ei_overP;
    vdrift[i]  = okV ? vabs : 0.0;
    alpha[i]   = okA ? a    : 0.0;
  }
}

//--------------------------------------------------
// α≈0 の最大 E を求めるヘルパー
//  - α が単調増加だと仮定
//  - 最初に α > tol となる手前の点の E を返す
//--------------------------------------------------
static double FindEthreshold(const std::vector<double>& E,
                             const std::vector<double>& alpha,
                             double tol = 1e-3) {
  const int N = static_cast<int>(E.size());
  if (N == 0 || alpha.size() != E.size()) return -1.0;

  int idxFirstPositive = -1;
  for (int i = 0; i < N; ++i) {
    if (alpha[i] > tol) {
      idxFirstPositive = i;
      break;
    }
  }
  if (idxFirstPositive <= 0) {
    // ぜんぶ 0 か、最初から >tol の場合
    return -1.0;
  }
  // alpha[idxFirstPositive-1] ≲ tol, alpha[idxFirstPositive] > tol
  return E[idxFirstPositive - 1];
}

//--------------------------------------------------
// メイン処理
//--------------------------------------------------
void gas_alpha_v_vsE(const char* gasfile) {
  // 1) ガス読み込み
  MediumMagboltz gas;
  if (!gas.LoadGasFile(gasfile)) {
    std::cerr << "Error: cannot load gas file " << gasfile << std::endl;
    return;
  }

  // 圧力を 0.1 atm (= 76 Torr) に上書き
  const double Patm_target = 0.1;                                // [atm]
  const double P_torr      = Patm_target * AtmosphericPressure;  // [Torr]
  gas.SetPressure(P_torr);
  const double temperature = gas.GetTemperature();
  const double P_bar       = P_torr / 750.061683;

  std::printf("Loaded gas: %s\n", gasfile);
  std::printf("  P = %.4g Torr (%.3g bar, %.2f atm), T = %.2f K\n",
              P_torr, P_bar, Patm_target, temperature);

  // 2) gas ファイルの E グリッドを取得（レンジをそのまま使う）
  std::vector<double> efields, bfields, angles;
  gas.GetFieldGrid(efields, bfields, angles);
  if (efields.empty()) {
    std::cerr << "Error: field grid is empty in gas file.\n";
    return;
  }

  const size_t Ntab  = efields.size();
  const double Emin  = efields.front();          // [V/cm]
  const double Emax  = efields.back();           // [V/cm]
  const int    N     = static_cast<int>(Ntab);

  std::printf("  E grid from gas file: Emin = %.3g V/cm, Emax = %.3g V/cm, N = %zu\n",
              Emin, Emax, Ntab);

  // --- (A) log 等間隔の E ---
  std::vector<double> E_log(N);
  {
    const double logEmin = std::log10(Emin);
    const double logEmax = std::log10(Emax);
    for (int i = 0; i < N; ++i) {
      const double logE = logEmin + (logEmax - logEmin) * double(i) / double(N - 1);
      E_log[i] = std::pow(10.0, logE);
    }
  }

  // --- (B) linear 等間隔の E ---
  std::vector<double> E_lin(N);
  {
    for (int i = 0; i < N; ++i) {
      E_lin[i] = Emin + (Emax - Emin) * double(i) / double(N - 1);
    }
  }

  // 3) ガス物性サンプリング（まず Torr ベースの E/p）
  std::vector<double> Ep_log_torr, alpha_log, vdrift_log;
  SampleGas(gas, E_log, P_torr, Ep_log_torr, alpha_log, vdrift_log);

  std::vector<double> Ep_lin_torr, alpha_lin, vdrift_lin;
  SampleGas(gas, E_lin, P_torr, Ep_lin_torr, alpha_lin, vdrift_lin);

  // 4) atm ベースの E/p = E / p[atm] に変換（描画用）
  std::vector<double> Ep_log_atm(N), Ep_lin_atm(N);
  for (int i = 0; i < N; ++i) {
    if (Patm_target > 0.) {
      Ep_log_atm[i] = E_log[i] / Patm_target; // [V/(cm*atm)]
      Ep_lin_atm[i] = E_lin[i] / Patm_target;
    } else {
      Ep_log_atm[i] = 0.;
      Ep_lin_atm[i] = 0.;
    }
  }

  // ★ 5) α ≈ 0 から立ち上がる E のしきい値を推定
  const double Eth = FindEthreshold(E_log, alpha_log, 1e-3); // tol=1e-3 [1/cm] くらい
  if (Eth > 0.0) {
    const double Eth_overP_torr = Eth / P_torr;
    const double Eth_overP_atm  = Eth / Patm_target;
    std::printf("\n[Threshold estimate]\n");
    std::printf("  alpha ~ 0 up to about E_th ≈ %.3g V/cm\n", Eth);
    std::printf("    E_th/p ≈ %.3g V/(cm*Torr) ≈ %.3g V/(cm*atm)\n",
                Eth_overP_torr, Eth_overP_atm);
  } else {
    std::printf("\n[Threshold estimate]\n  Could not determine E_th in this E range.\n");
  }

  // ----------------------------------------------------
  // 6-A) x = E, log-x 版
  // ----------------------------------------------------
  {
    TCanvas* cElog = new TCanvas("c_gas_E_log", "alpha & v vs E (log x)", 1200, 500);
    cElog->Divide(2, 1);

    // 左: α(E)（log-x, log-y）
    cElog->cd(1);
    SetupPad(true, true);
    auto* gA_Elog = new TGraph(N, E_log.data(), alpha_log.data());
    gA_Elog->SetTitle(";electric field [V/cm];#alpha [1/cm]");
    gA_Elog->SetMarkerStyle(20);
    gA_Elog->SetLineWidth(2);
    gA_Elog->Draw("APL");

    // 右: v(E)（log-x, linear-y）
    cElog->cd(2);
    SetupPad(true, false);
    auto* gV_Elog = new TGraph(N, E_log.data(), vdrift_log.data());
    gV_Elog->SetTitle(";electric field [V/cm];drift velocity [cm/ns]");
    gV_Elog->SetMarkerStyle(20);
    gV_Elog->SetLineWidth(2);
    gV_Elog->Draw("APL");

    cElog->SaveAs("gas_alpha_v_vsE_logx.png");
  }

  // ----------------------------------------------------
  // 6-B) x = E, linear-x 版
  // ----------------------------------------------------
  {
    TCanvas* cElin = new TCanvas("c_gas_E_lin", "alpha & v vs E (linear x)", 1200, 500);
    cElin->Divide(2, 1);

    // 左: α(E)（linear-x, log-y）
    cElin->cd(1);
    SetupPad(false, true);
    auto* gA_Elin = new TGraph(N, E_lin.data(), alpha_lin.data());
    gA_Elin->SetTitle(";electric field [V/cm];#alpha [1/cm]");
    gA_Elin->SetMarkerStyle(20);
    gA_Elin->SetLineWidth(2);
    gA_Elin->Draw("APL");

    // 右: v(E)（linear-x, linear-y）
    cElin->cd(2);
    SetupPad(false, false);
    auto* gV_Elin = new TGraph(N, E_lin.data(), vdrift_lin.data());
    gV_Elin->SetTitle(";electric field [V/cm];drift velocity [cm/ns]");
    gV_Elin->SetMarkerStyle(20);
    gV_Elin->SetLineWidth(2);
    gV_Elin->Draw("APL");

    cElin->SaveAs("gas_alpha_v_vsE_linearx.png");
  }

  // ----------------------------------------------------
  // 6-C) x = E/p(atm), log-x 版
  // ----------------------------------------------------
  {
    TCanvas* cEp_log = new TCanvas("c_gas_Ep_log", "alpha & v vs E/p (log x, atm)", 1200, 500);
    cEp_log->Divide(2, 1);

    // 左: α(E/p)（log-x, log-y）
    cEp_log->cd(1);
    SetupPad(true, true);
    auto* gA_Ep_log = new TGraph(N, Ep_log_atm.data(), alpha_log.data());
    gA_Ep_log->SetTitle(";E/p [V/(cm#timesatm)];#alpha [1/cm]");
    gA_Ep_log->SetMarkerStyle(20);
    gA_Ep_log->SetLineWidth(2);
    gA_Ep_log->Draw("APL");

    // 右: v(E/p)（log-x, linear-y）
    cEp_log->cd(2);
    SetupPad(true, false);
    auto* gV_Ep_log = new TGraph(N, Ep_log_atm.data(), vdrift_log.data());
    gV_Ep_log->SetTitle(";E/p [V/(cm#timesatm)];drift velocity [cm/ns]");
    gV_Ep_log->SetMarkerStyle(20);
    gV_Ep_log->SetLineWidth(2);
    gV_Ep_log->Draw("APL");

    cEp_log->SaveAs("gas_alpha_v_vsEoverP_logx_atm.png");
  }

  // ----------------------------------------------------
  // 6-D) x = E/p(atm), linear-x 版
  // ----------------------------------------------------
  {
    TCanvas* cEp_lin = new TCanvas("c_gas_Ep_lin", "alpha & v vs E/p (linear x, atm)", 1200, 500);
    cEp_lin->Divide(2, 1);

    // 左: α(E/p)（linear-x, log-y）
    cEp_lin->cd(1);
    SetupPad(false, true);
    auto* gA_Ep_lin = new TGraph(N, Ep_lin_atm.data(), alpha_lin.data());
    gA_Ep_lin->SetTitle(";E/p [V/(cm#timesatm)];#alpha [1/cm]");
    gA_Ep_lin->SetMarkerStyle(20);
    gA_Ep_lin->SetLineWidth(2);
    gA_Ep_lin->Draw("APL");

    // 右: v(E/p)（linear-x, linear-y）
    cEp_lin->cd(2);
    SetupPad(false, false);
    auto* gV_Ep_lin = new TGraph(N, Ep_lin_atm.data(), vdrift_lin.data());
    gV_Ep_lin->SetTitle(";E/p [V/(cm#timesatm)];drift velocity [cm/ns]");
    gV_Ep_lin->SetMarkerStyle(20);
    gV_Ep_lin->SetLineWidth(2);
    gV_Ep_lin->Draw("APL");

    cEp_lin->SaveAs("gas_alpha_v_vsEoverP_linearx_atm.png");
  }
}

//--------------------------------------------------
// main（GUI 有り）
//--------------------------------------------------
int main(int argc, char** argv) {
  TApplication app("app", &argc, argv);
  gStyle->SetOptStat(0);

  const char* gasfile = "ic4H10_100_0.1Torr.gas";
  if (argc > 1) {
    gasfile = argv[1];
  }

  gas_alpha_v_vsE(gasfile);

  app.Run();   // GUI で軸・レンジは自由にいじる
  return 0;
}
