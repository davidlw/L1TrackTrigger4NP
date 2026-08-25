// ---------------------------------------------------------------------------
// plot_offline_perf.C -- draw the offline efficiency, fake rate and duplicate
// rate produced by offline_perf.C, and dump the binned numbers as text.
//
// One curve per panel: offline reconstruction is a single algorithm, so there is
// nothing to compare against. Axis, binning and styling match the L1 macros in
// ../../L1_tracking_performance/scripts, so an offline panel and the corresponding
// L1 panel can be put side by side, or overlaid, without rescaling.
//
// The y axis is fixed to 0-1.35 with a guide line at 1 for all three quantities,
// as in the L1 plots, so efficiency and fake rate are read on the same scale.
//
// Usage: root -l -b -q 'plot_offline_perf.C("../output/offperf_qed_mumu.root")'
// ---------------------------------------------------------------------------

#include <cstdio>

namespace {

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

TGraphAsymmErrors* MakeRate(TFile* f, const char* num, const char* den, int color, int marker) {
  TH1D* hn = (TH1D*)f->Get(num);
  TH1D* hd = (TH1D*)f->Get(den);
  if (!hn || !hd) { printf("[error] missing %s or %s\n", num, den); return nullptr; }
  auto* g = new TGraphAsymmErrors(hn, hd, "cl=0.683 b(1,1) mode");
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  return g;
}

// numPfx/denPfx select the quantity:
//   ("num","den")   efficiency        ("fnum","fden") fake rate
//   ("numD","denM") duplicate rate
void DrawPanel(TFile* f, const char* var, const char* xtitle, const char* title,
               const char* outfile, double xlo, double xhi,
               const char* numPfx, const char* denPfx, const char* ytitle,
               const char* legend) {
  auto* c = new TCanvas(Form("c_%s_%s", numPfx, var), "", 800, 650);
  c->SetGridy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);

  auto* frame = c->DrawFrame(xlo, 0.0, xhi, 1.35);
  auto* one = new TLine(xlo, 1.0, xhi, 1.0);
  one->SetLineStyle(3);
  one->SetLineColor(kGray + 1);
  one->Draw();
  frame->GetXaxis()->SetTitle(xtitle);
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.35);

  auto* g = MakeRate(f, Form("%s_%s_off", numPfx, var), Form("%s_%s_off", denPfx, var),
                     kBlack, 20);
  if (g) g->Draw("P same");

  auto* leg = new TLegend(0.17, 0.78, 0.48, 0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.036);
  if (g) leg->AddEntry(g, legend, "lp");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.17, 0.94, title);

  SaveBoth(c, outfile);
}

void DumpTable(TFile* f, const char* var, const char* numPfx, const char* denPfx,
               const char* label) {
  TH1D* hd = (TH1D*)f->Get(Form("%s_%s_off", denPfx, var));
  TH1D* hn = (TH1D*)f->Get(Form("%s_%s_off", numPfx, var));
  if (!hd || !hn) return;
  printf("\n=== %s vs %s ===\n", label, var);
  printf("%14s %12s %12s %10s\n", "bin", "numerator", "denominator", "rate");
  for (int i = 1; i <= hd->GetNbinsX(); ++i) {
    double a = hn->GetBinContent(i), b = hd->GetBinContent(i);
    if (b < 1) continue;
    printf("%6.2f-%6.2f %12.0f %12.0f %10.4f\n",
           hd->GetBinLowEdge(i), hd->GetBinLowEdge(i + 1), a, b, b > 0 ? a / b : 0);
  }
}

}  // namespace

void plot_offline_perf(const char* fname = "../output/offperf_qed_mumu.root",
                       const char* tag = "",
                       const char* sample = "STARlight QED #mu#mu, offline tracks") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("[error] cannot open %s\n", fname); return; }

  const double P = TMath::Pi();
  TString t0 = TString(sample);
  TString t1 = TString(sample) + ", p_{T} > 1 GeV, |#eta| < 2.4";
  TString t2 = TString(sample) + ", p_{T} > 2 GeV, |#eta| < 2.4";

  // ---------------- efficiency ----------------
  DrawPanel(f, "pt", "tracking particle p_{T} [GeV]", t0.Data(),
            Form("../figures/off_eff_vs_pt%s.pdf", tag), 0.0, 10.0,
            "num", "den", "offline tracking efficiency", "Offline");
  DrawPanel(f, "eta2", "tracking particle #eta", t2.Data(),
            Form("../figures/off_eff_vs_eta%s.pdf", tag), -2.4, 2.4,
            "num", "den", "offline tracking efficiency", "Offline");
  DrawPanel(f, "phi2", "tracking particle #phi", t2.Data(),
            Form("../figures/off_eff_vs_phi%s.pdf", tag), -P, P,
            "num", "den", "offline tracking efficiency", "Offline");
  DrawPanel(f, "z02", "tracking particle z_{0} [cm]", t2.Data(),
            Form("../figures/off_eff_vs_z0%s.pdf", tag), -20, 20,
            "num", "den", "offline tracking efficiency", "Offline");
  DrawPanel(f, "eta1", "tracking particle #eta", t1.Data(),
            Form("../figures/off_eff_vs_eta_pt1%s.pdf", tag), -2.4, 2.4,
            "num", "den", "offline tracking efficiency", "Offline");
  DrawPanel(f, "phi1", "tracking particle #phi", t1.Data(),
            Form("../figures/off_eff_vs_phi_pt1%s.pdf", tag), -P, P,
            "num", "den", "offline tracking efficiency", "Offline");

  // ---------------- fake rate ----------------
  DrawPanel(f, "pt", "track p_{T} [GeV]", t0.Data(),
            Form("../figures/off_fake_vs_pt%s.pdf", tag), 0.0, 10.0,
            "fnum", "fden", "offline fake rate", "Offline");
  DrawPanel(f, "eta2", "track #eta", t2.Data(),
            Form("../figures/off_fake_vs_eta%s.pdf", tag), -2.4, 2.4,
            "fnum", "fden", "offline fake rate", "Offline");
  DrawPanel(f, "phi2", "track #phi", t2.Data(),
            Form("../figures/off_fake_vs_phi%s.pdf", tag), -P, P,
            "fnum", "fden", "offline fake rate", "Offline");
  DrawPanel(f, "z02", "track z_{0} [cm]", t2.Data(),
            Form("../figures/off_fake_vs_z0%s.pdf", tag), -20, 20,
            "fnum", "fden", "offline fake rate", "Offline");
  DrawPanel(f, "eta1", "track #eta", t1.Data(),
            Form("../figures/off_fake_vs_eta_pt1%s.pdf", tag), -2.4, 2.4,
            "fnum", "fden", "offline fake rate", "Offline");
  DrawPanel(f, "phi1", "track #phi", t1.Data(),
            Form("../figures/off_fake_vs_phi_pt1%s.pdf", tag), -P, P,
            "fnum", "fden", "offline fake rate", "Offline");

  // ---------------- duplicate rate ----------------
  TString d0 = TString(sample) + " -- duplicate rate";
  TString d1 = TString(sample) + " -- duplicate rate, p_{T} > 1 GeV, |#eta| < 2.4";
  TString d2 = TString(sample) + " -- duplicate rate, p_{T} > 2 GeV, |#eta| < 2.4";
  DrawPanel(f, "pt", "tracking particle p_{T} [GeV]", d0.Data(),
            Form("../figures/off_dup_vs_pt%s.pdf", tag), 0.0, 10.0,
            "numD", "denM", "duplicate rate", "Offline");
  DrawPanel(f, "eta2", "tracking particle #eta", d2.Data(),
            Form("../figures/off_dup_vs_eta%s.pdf", tag), -2.4, 2.4,
            "numD", "denM", "duplicate rate", "Offline");
  DrawPanel(f, "phi2", "tracking particle #phi", d2.Data(),
            Form("../figures/off_dup_vs_phi%s.pdf", tag), -P, P,
            "numD", "denM", "duplicate rate", "Offline");
  DrawPanel(f, "eta1", "tracking particle #eta", d1.Data(),
            Form("../figures/off_dup_vs_eta_pt1%s.pdf", tag), -2.4, 2.4,
            "numD", "denM", "duplicate rate", "Offline");
  DrawPanel(f, "phi1", "tracking particle #phi", d1.Data(),
            Form("../figures/off_dup_vs_phi_pt1%s.pdf", tag), -P, P,
            "numD", "denM", "duplicate rate", "Offline");

  DumpTable(f, "pt", "num", "den", "efficiency");
  DumpTable(f, "pt", "fnum", "fden", "fake rate");
  DumpTable(f, "pt", "numD", "denM", "duplicate rate");
  DumpTable(f, "eta2", "num", "den", "efficiency (pT > 2, |eta| < 2.4)");

  printf("\nfigures written to ../figures/\n");
}
