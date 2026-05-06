#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TLatex.h>

#include <fstream>
#include <iostream>
#include <iomanip>

void draw_diff_x(const char* input_root = "merged_track_results.root",
                 const char* output_png = "diff_x_hist.png",
                 int nbins = 100,
                 double xmin = -0.2,
                 double xmax = 0.2,
                 const char* output_csv = "diff_x_hist.csv") {
  TFile* f = TFile::Open(input_root, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "[ERROR] cannot open file: " << input_root << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(f->Get("tree"));
  if (!tree) {
    std::cerr << "[ERROR] tree 'tree' not found." << std::endl;
    f->Close();
    return;
  }

  gStyle->SetOptStat(1110);

  TH1D* h = new TH1D("h_diff_x",
                     "Position residual;diff_x = x_{reco} - x_{true} [cm];Entries",
                     nbins, xmin, xmax);

  tree->Draw("diff_x>>h_diff_x", "", "goff");

  TCanvas* c = new TCanvas("c_diff_x", "diff_x", 800, 600);
  c->SetGrid();
  h->SetLineWidth(2);
  h->Draw();

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.03);
  latex.DrawLatex(0.60, 0.85, Form("Entries = %.0f", h->GetEntries()));
  latex.DrawLatex(0.60, 0.80, Form("Mean = %.6f cm", h->GetMean()));
  latex.DrawLatex(0.60, 0.75, Form("RMS = %.6f cm", h->GetRMS()));

  c->SaveAs(output_png);

  std::ofstream fout(output_csv);
  if (fout) {
    fout << "x_center_cm,bin_width_cm,count\n";
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
      fout << std::setprecision(10)
           << h->GetBinCenter(i) << ","
           << h->GetBinWidth(i) << ","
           << h->GetBinContent(i) << "\n";
    }
  }

  std::cout << "[DONE] input   = " << input_root << std::endl;
  std::cout << "[DONE] output  = " << output_png << std::endl;
  std::cout << "[DONE] csv     = " << output_csv << std::endl;
  std::cout << "[DONE] entries = " << h->GetEntries() << std::endl;
  std::cout << "[DONE] mean    = " << h->GetMean() << " cm" << std::endl;
  std::cout << "[DONE] rms     = " << h->GetRMS() << " cm" << std::endl;

  f->Close();
}