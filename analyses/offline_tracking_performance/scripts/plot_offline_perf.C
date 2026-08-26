// ---------------------------------------------------------------------------
// plot_offline_perf.C -- draw the offline efficiency, fake rate and duplicate
// rate produced by offline_perf.C, and dump the binned numbers as text.
//
// pT panels use a LOG x axis from 0.3 to 10 GeV, matching the fine low-pT binning
// where the turn-on lives.
//
// eta, phi and Nch panels OVERLAY the four pT thresholds (0.3, 0.6, 1, 2 GeV)
// that offline_perf.C filled. Everything is inside |eta| < 2.4.
//
// The Nch x range is read from the histograms rather than hard-coded, because the
// multiplicity axis is sample dependent (set by nchMax when they were filled).
//
// Offline reconstruction is a single algorithm, so within a panel the curves
// differ only by the pT threshold, never by algorithm.
//
// The y axis is fixed to 0-1.35 with a guide line at 1 for all three quantities,
// as in the L1 plots, so efficiency and fake rate are read on the same scale.
// The legend sits low on efficiency panels and high on fake/duplicate panels,
// since those curves live at opposite ends of that range.
//
// Usage: root -l -b -q 'plot_offline_perf.C("../output/offperf_qedmumu_hp.root","_hp","STARlight QED #mu#mu, offline highPurity")'
// ---------------------------------------------------------------------------

#include <cstdio>

namespace {

const int kNThr = 4;
const char* kThrLabel[kNThr] = {"p_{T} > 0.3 GeV", "p_{T} > 0.6 GeV",
                                "p_{T} > 1 GeV", "p_{T} > 2 GeV"};
const int kCol[kNThr] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2};
const int kMrk[kNThr] = {20, 21, 22, 33};

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

TLegend* MakeLegend(bool low) {
  // efficiency sits near 1 -> put the legend low; fake/duplicate sit near 0 -> high
  auto* leg = low ? new TLegend(0.17, 0.16, 0.52, 0.40)
                  : new TLegend(0.17, 0.64, 0.52, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.034);
  return leg;
}

// numPfx/denPfx select the quantity:
//   ("num","den")   efficiency   ("fnum","fden") fake rate   ("numD","denM") duplicate
void DrawPt(TFile* f, const char* xtitle, const char* title, const char* outfile,
            const char* numPfx, const char* denPfx, const char* ytitle,
            const char* legend, bool legLow) {
  auto* c = new TCanvas(Form("cpt_%s", numPfx), "", 800, 650);
  c->SetGridy();
  c->SetLogx();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);

  auto* frame = c->DrawFrame(0.3, 0.0, 10.0, 1.35);
  auto* one = new TLine(0.3, 1.0, 10.0, 1.0);
  one->SetLineStyle(3); one->SetLineColor(kGray + 1); one->Draw();
  frame->GetXaxis()->SetTitle(xtitle);
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->GetXaxis()->SetMoreLogLabels();
  frame->GetXaxis()->SetNoExponent();

  auto* g = MakeRate(f, Form("%s_pt_off", numPfx), Form("%s_pt_off", denPfx), kBlack, 20);
  if (g) g->Draw("P same");

  auto* leg = MakeLegend(legLow);
  if (g) leg->AddEntry(g, legend, "lp");
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
  tx.DrawLatex(0.17, 0.94, title);
  SaveBoth(c, outfile);
}

// eta / phi with the four pT thresholds overlaid
void DrawOverlay(TFile* f, const char* var, const char* xtitle, const char* title,
                 const char* outfile, double xlo, double xhi,
                 const char* numPfx, const char* denPfx, const char* ytitle, bool legLow) {
  auto* c = new TCanvas(Form("c_%s_%s", numPfx, var), "", 800, 650);
  c->SetGridy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);

  // xhi <= xlo means "take the range from the histogram" -- used for Nch, whose
  // axis is sample dependent and set when the histograms were filled.
  if (xhi <= xlo) {
    TH1D* h0 = (TH1D*)f->Get(Form("%s_%s0_off", denPfx, var));
    if (h0) { xlo = h0->GetXaxis()->GetXmin(); xhi = h0->GetXaxis()->GetXmax(); }
    else { xlo = 0; xhi = 1; }
  }

  auto* frame = c->DrawFrame(xlo, 0.0, xhi, 1.35);
  auto* one = new TLine(xlo, 1.0, xhi, 1.0);
  one->SetLineStyle(3); one->SetLineColor(kGray + 1); one->Draw();
  frame->GetXaxis()->SetTitle(xtitle);
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.35);

  auto* leg = MakeLegend(legLow);
  for (int i = 0; i < kNThr; ++i) {
    auto* g = MakeRate(f, Form("%s_%s%d_off", numPfx, var, i),
                       Form("%s_%s%d_off", denPfx, var, i), kCol[i], kMrk[i]);
    if (!g) continue;
    g->Draw("P same");
    leg->AddEntry(g, kThrLabel[i], "lp");
  }
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
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
  TString t0 = TString(sample) + ", |#eta| < 2.4";

  // ---------------- vs pT (log x) ----------------
  DrawPt(f, "tracking particle p_{T} [GeV]", t0.Data(),
         Form("../figures/off_eff_vs_pt%s.pdf", tag),
         "num", "den", "offline tracking efficiency", "Offline", true);
  DrawPt(f, "track p_{T} [GeV]", t0.Data(),
         Form("../figures/off_fake_vs_pt%s.pdf", tag),
         "fnum", "fden", "offline fake rate", "Offline", false);
  DrawPt(f, "tracking particle p_{T} [GeV]", t0.Data(),
         Form("../figures/off_dup_vs_pt%s.pdf", tag),
         "numD", "denM", "duplicate rate", "Offline", false);

  // ---------------- vs eta / phi, four thresholds overlaid ----------------
  DrawOverlay(f, "eta", "tracking particle #eta", t0.Data(),
              Form("../figures/off_eff_vs_eta%s.pdf", tag), -2.4, 2.4,
              "num", "den", "offline tracking efficiency", true);
  DrawOverlay(f, "phi", "tracking particle #phi", t0.Data(),
              Form("../figures/off_eff_vs_phi%s.pdf", tag), -P, P,
              "num", "den", "offline tracking efficiency", true);

  DrawOverlay(f, "eta", "track #eta", t0.Data(),
              Form("../figures/off_fake_vs_eta%s.pdf", tag), -2.4, 2.4,
              "fnum", "fden", "offline fake rate", false);
  DrawOverlay(f, "phi", "track #phi", t0.Data(),
              Form("../figures/off_fake_vs_phi%s.pdf", tag), -P, P,
              "fnum", "fden", "offline fake rate", false);

  DrawOverlay(f, "eta", "tracking particle #eta", t0.Data(),
              Form("../figures/off_dup_vs_eta%s.pdf", tag), -2.4, 2.4,
              "numD", "denM", "duplicate rate", false);
  DrawOverlay(f, "phi", "tracking particle #phi", t0.Data(),
              Form("../figures/off_dup_vs_phi%s.pdf", tag), -P, P,
              "numD", "denM", "duplicate rate", false);

  // ---------------- vs Nch, four thresholds overlaid ----------------
  // x range is read from the histograms, since Nch is sample dependent.
  const char* nchTitle = "N_{ch} (truth, p_{T} > 0.4, |#eta| < 2.4)";
  DrawOverlay(f, "nch", nchTitle, t0.Data(),
              Form("../figures/off_eff_vs_nch%s.pdf", tag), 0, 0,
              "num", "den", "offline tracking efficiency", true);
  DrawOverlay(f, "nch", nchTitle, t0.Data(),
              Form("../figures/off_fake_vs_nch%s.pdf", tag), 0, 0,
              "fnum", "fden", "offline fake rate", false);
  DrawOverlay(f, "nch", nchTitle, t0.Data(),
              Form("../figures/off_dup_vs_nch%s.pdf", tag), 0, 0,
              "numD", "denM", "duplicate rate", false);

  DumpTable(f, "pt", "num", "den", "efficiency");
  DumpTable(f, "pt", "fnum", "fden", "fake rate");
  DumpTable(f, "pt", "numD", "denM", "duplicate rate");

  printf("\nfigures written to ../figures/\n");
}
