// ---------------------------------------------------------------------------
// plot_eff.C -- draw the Default vs Dummy stub efficiency comparison produced
// by compare_eff.C, and dump the binned numbers as text.
//
// Usage: root -l -b -q 'plot_eff.C("eff_qed_mumu.root")'
// ---------------------------------------------------------------------------

#include <cstdio>

namespace {

// Write every figure as both PDF and PNG.
void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

TGraphAsymmErrors* MakeEff(TFile* f, const char* num, const char* den, int color, int marker) {
  TH1D* hn = (TH1D*)f->Get(num);
  TH1D* hd = (TH1D*)f->Get(den);
  if (!hn || !hd) {
    printf("[error] missing %s or %s\n", num, den);
    return nullptr;
  }
  auto* g = new TGraphAsymmErrors(hn, hd, "cl=0.683 b(1,1) mode");
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  return g;
}

// numPfx/denPfx select the quantity: ("num","den") = efficiency,
// ("numD","denM") = duplicate rate. Axis and styling are identical either way.
void DrawPanelGen(TFile* f, const char* var, const char* xtitle, const char* title,
                  const char* outfile, double xlo, double xhi, const char* numPfx,
                  const char* denPfx, const char* ytitle) {
  auto* c = new TCanvas(Form("c_%s_%s", numPfx, var), "", 800, 650);
  c->SetGridy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);

  // Ceiling at 1.35 rather than 1.05: efficiency plateaus near 0.93, which with a
  // 1.05 ceiling sat directly under the legend. The extra headroom keeps the legend
  // clear of both the curves and the title.
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

  auto* gDef = MakeEff(f, Form("%s_%s_def", numPfx, var), Form("%s_%s_def", denPfx, var),
                       kBlue + 1, 20);
  auto* gDum = MakeEff(f, Form("%s_%s_dum", numPfx, var), Form("%s_%s_dum", denPfx, var),
                       kRed + 1, 21);
  if (gDef) gDef->Draw("P same");
  if (gDum) gDum->Draw("P same");

  auto* leg = new TLegend(0.17, 0.74, 0.48, 0.87);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.036);
  if (gDef) leg->AddEntry(gDef, "Default Stub", "lp");
  if (gDum) leg->AddEntry(gDum, "Dummy Stub", "lp");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.17, 0.94, title);

  SaveBoth(c, outfile);
}

void DrawPanel(TFile* f, const char* var, const char* xtitle, const char* title,
               const char* outfile, double xlo, double xhi) {
  DrawPanelGen(f, var, xtitle, title, outfile, xlo, xhi, "num", "den",
               "L1 tracking efficiency");
}

void DumpTable(TFile* f, const char* var, const char* label) {
  TH1D* nd = (TH1D*)f->Get(Form("num_%s_def", var));
  TH1D* dd = (TH1D*)f->Get(Form("den_%s_def", var));
  TH1D* nu = (TH1D*)f->Get(Form("num_%s_dum", var));
  TH1D* du = (TH1D*)f->Get(Form("den_%s_dum", var));
  if (!nd || !dd || !nu || !du) return;

  printf("\n=== efficiency vs %s (%s) ===\n", var, label);
  printf("%14s %10s %10s %9s %10s %10s %9s   %s\n", "bin", "num(def)", "den(def)", "eff(def)",
         "num(dum)", "den(dum)", "eff(dum)", "dummy-default");
  for (int i = 1; i <= dd->GetNbinsX(); ++i) {
    double a = nd->GetBinContent(i), b = dd->GetBinContent(i);
    double x = nu->GetBinContent(i), y = du->GetBinContent(i);
    if (b < 1 && y < 1) continue;
    double ea = b > 0 ? a / b : 0, ex = y > 0 ? x / y : 0;
    printf("%6.2f-%6.2f %10.0f %10.0f %9.4f %10.0f %10.0f %9.4f   %+9.4f\n",
           dd->GetBinLowEdge(i), dd->GetBinLowEdge(i + 1), a, b, ea, x, y, ex, ex - ea);
  }
}

// Fold the fine phi efficiency onto an N-fold symmetry and measure how
// non-flat the result is. If the underlying modulation really has period
// 2pi/N, folding on N stacks the structure coherently and chi2/ndf against a
// flat line blows up; folding on the wrong N averages it away.
//
// The 144-bin histogram folds cleanly onto N in {8,9,12,16,18,24}.
void FoldScan(TFile* f, const char* smp, const char* label) {
  auto* hn = (TH1D*)f->Get(Form("num_phiF_%s", smp));
  auto* hd = (TH1D*)f->Get(Form("den_phiF_%s", smp));
  if (!hn || !hd) {
    printf("[warn] no fine phi histograms for %s\n", smp);
    return;
  }
  const int nb = hd->GetNbinsX();  // 144
  printf("\n  %s\n", label);
  printf("  %5s %8s %10s %10s %12s\n", "N", "sub-bins", "chi2/ndf", "p-value", "amplitude");
  for (int N : {8, 9, 12, 16, 18, 24}) {
    if (nb % N != 0) continue;
    int sub = nb / N;
    std::vector<double> num(sub, 0), den(sub, 0);
    for (int i = 1; i <= nb; ++i) {
      int k = (i - 1) % sub;
      num[k] += hn->GetBinContent(i);
      den[k] += hd->GetBinContent(i);
    }
    double tn = 0, td = 0;
    for (int k = 0; k < sub; ++k) { tn += num[k]; td += den[k]; }
    if (td <= 0) continue;
    double e0 = tn / td;
    double chi2 = 0;
    int ndf = 0;
    double lo = 1e9, hi = -1e9;
    for (int k = 0; k < sub; ++k) {
      if (den[k] < 50) continue;
      double e = num[k] / den[k];
      double s = std::sqrt(e0 * (1 - e0) / den[k]);
      chi2 += (e - e0) * (e - e0) / (s * s);
      ++ndf;
      lo = std::min(lo, e);
      hi = std::max(hi, e);
    }
    if (ndf < 2) continue;
    ndf -= 1;  // one parameter (the mean) taken from the data
    printf("  %5d %8d %10.2f %10.3g %12.4f\n", N, sub, chi2 / ndf,
           TMath::Prob(chi2, ndf), hi - lo);
  }
}

}  // namespace

void plot_eff(const char* fname = "eff_qed_mumu.root",
              const char* tag = "",
              const char* sample = "STARlight QED #mu#mu") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("[error] cannot open %s\n", fname);
    return;
  }

  TString common = TString(sample);
  DrawPanel(f, "pt", "tracking particle p_{T} [GeV]", common.Data(),
            Form("../figures/eff_vs_pt%s.pdf", tag), 0.0, 10.0);

  // eta / phi / z0 are only meaningful above the track finder's threshold, so
  // these use the pT > 2 GeV histograms.
  TString common2 = TString(sample) + ", p_{T} > 2 GeV, |#eta| < 2.4";
  DrawPanel(f, "eta2", "tracking particle #eta", common2.Data(),
            Form("../figures/eff_vs_eta%s.pdf", tag), -2.4, 2.4);
  DrawPanel(f, "phi2", "tracking particle #phi", common2.Data(),
            Form("../figures/eff_vs_phi%s.pdf", tag), -TMath::Pi(), TMath::Pi());
  DrawPanel(f, "z02", "tracking particle z_{0} [cm]", common2.Data(),
            Form("../figures/eff_vs_z0%s.pdf", tag), -20, 20);

  // pT > 1 GeV, |eta| < 2.4 -- the comparison for the next dummy-tracking iteration
  TString common1 = TString(sample) + ", p_{T} > 1 GeV, |#eta| < 2.4";
  DrawPanel(f, "eta1", "tracking particle #eta", common1.Data(),
            Form("../figures/eff_vs_eta_pt1%s.pdf", tag), -2.4, 2.4);
  DrawPanel(f, "phi1", "tracking particle #phi", common1.Data(),
            Form("../figures/eff_vs_phi_pt1%s.pdf", tag), -TMath::Pi(), TMath::Pi());

  // Duplicate rate: fraction of matched TPs picked up by more than one track.
  TString d0 = TString(sample) + " -- duplicate rate";
  DrawPanelGen(f, "pt", "tracking particle p_{T} [GeV]", d0.Data(),
               Form("../figures/dup_vs_pt%s.pdf", tag), 0.0, 10.0, "numD", "denM",
               "duplicate rate");
  TString d2 = TString(sample) + " -- duplicate rate, p_{T} > 2 GeV, |#eta| < 2.4";
  DrawPanelGen(f, "eta2", "tracking particle #eta", d2.Data(),
               Form("../figures/dup_vs_eta%s.pdf", tag), -2.4, 2.4, "numD", "denM",
               "duplicate rate");
  DrawPanelGen(f, "phi2", "tracking particle #phi", d2.Data(),
               Form("../figures/dup_vs_phi%s.pdf", tag), -TMath::Pi(), TMath::Pi(),
               "numD", "denM", "duplicate rate");

  DumpTable(f, "pt", "all pT");
  DumpTable(f, "eta2", "pT > 2 GeV, |eta| < 2.4");
  DumpTable(f, "phi2", "pT > 2 GeV, |eta| < 2.4");
  DumpTable(f, "eta1", "pT > 1 GeV, |eta| < 2.4");
  DumpTable(f, "phi1", "pT > 1 GeV, |eta| < 2.4");

  printf("\n=== phi periodicity scan (pT > 2 GeV, |eta| < 2.4) ===\n");
  printf("A large chi2/ndf at one N means the efficiency really is modulated with\n");
  printf("period 2pi/N. A flat scan across all N means no phi structure.\n");
  FoldScan(f, "def", "Default Stub");
  FoldScan(f, "dum", "Dummy Stub");

  printf("\nfigures written to ../figures/\n");
}
