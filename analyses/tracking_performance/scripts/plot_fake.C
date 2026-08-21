// ---------------------------------------------------------------------------
// plot_fake.C -- draw the Default vs Dummy fake-rate comparison produced by
// compare_fake.C, and dump the binned numbers as text.
//
// Same style as plot_eff.C so fake rate and efficiency can be shown side by side.
// The y axis is fixed to 0-1.35 with a guide line at 1, identical to plot_eff.C,
// so fake rate and efficiency can be shown side by side on the same scale.
//
// Usage: root -l -b -q 'plot_fake.C("fake_qed_mumu.root")'
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

void DrawPanel(TFile* f, const char* var, const char* xtitle, const char* title,
               const char* outfile, double xlo, double xhi) {
  auto* c = new TCanvas(Form("cf_%s", var), "", 800, 650);
  c->SetGridy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);

  // Same axis as plot_eff.C: fixed 0-1.35 with a guide line at 1, so fake rate and
  // efficiency can be read side by side on identical scales.
  auto* frame = c->DrawFrame(xlo, 0.0, xhi, 1.35);
  auto* one = new TLine(xlo, 1.0, xhi, 1.0);
  one->SetLineStyle(3);
  one->SetLineColor(kGray + 1);
  one->Draw();
  frame->GetXaxis()->SetTitle(xtitle);
  frame->GetYaxis()->SetTitle("L1 track fake rate");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.35);

  auto* gDef = MakeRate(f, Form("numF_%s_def", var), Form("den_%s_def", var), kBlue + 1, 20);
  auto* gDum = MakeRate(f, Form("numF_%s_dum", var), Form("den_%s_dum", var), kRed + 1, 21);
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

void DumpTable(TFile* f, const char* var, const char* label) {
  TH1D* dd = (TH1D*)f->Get(Form("den_%s_def", var));
  TH1D* fd = (TH1D*)f->Get(Form("numF_%s_def", var));
  TH1D* du = (TH1D*)f->Get(Form("den_%s_dum", var));
  TH1D* fu = (TH1D*)f->Get(Form("numF_%s_dum", var));
  if (!dd || !fd || !du || !fu) return;
  printf("\n=== fake rate vs %s (%s) ===\n", var, label);
  printf("%14s %10s %10s %10s %10s %10s %10s   %s\n", "bin", "fake(def)", "trk(def)",
         "rate(def)", "fake(dum)", "trk(dum)", "rate(dum)", "dummy-default");
  for (int i = 1; i <= dd->GetNbinsX(); ++i) {
    double a = fd->GetBinContent(i), b = dd->GetBinContent(i);
    double x = fu->GetBinContent(i), y = du->GetBinContent(i);
    if (b < 1 && y < 1) continue;
    double ra = b > 0 ? a / b : 0, rx = y > 0 ? x / y : 0;
    printf("%6.2f-%6.2f %10.0f %10.0f %10.4f %10.0f %10.0f %10.4f   %+9.4f\n",
           dd->GetBinLowEdge(i), dd->GetBinLowEdge(i + 1), a, b, ra, x, y, rx, rx - ra);
  }
}

}  // namespace

void plot_fake(const char* fname = "fake_qed_mumu.root", const char* tag = "",
               const char* sample = "STARlight QED #mu#mu") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("[error] cannot open %s\n", fname); return; }

  TString t0 = TString(sample);
  DrawPanel(f, "pt", "track p_{T} [GeV]", t0.Data(),
            Form("../figures/fake_vs_pt%s.pdf", tag), 0.0, 10.0);

  TString t2 = TString(sample) + ", p_{T} > 2 GeV, |#eta| < 2.4";
  DrawPanel(f, "eta2", "track #eta", t2.Data(),
            Form("../figures/fake_vs_eta%s.pdf", tag), -2.4, 2.4);
  DrawPanel(f, "phi2", "track #phi", t2.Data(),
            Form("../figures/fake_vs_phi%s.pdf", tag), -TMath::Pi(), TMath::Pi());
  DrawPanel(f, "z02", "track z_{0} [cm]", t2.Data(),
            Form("../figures/fake_vs_z0%s.pdf", tag), -20, 20);

  TString t1 = TString(sample) + ", p_{T} > 1 GeV, |#eta| < 2.4";
  DrawPanel(f, "eta1", "track #eta", t1.Data(),
            Form("../figures/fake_vs_eta_pt1%s.pdf", tag), -2.4, 2.4);
  DrawPanel(f, "phi1", "track #phi", t1.Data(),
            Form("../figures/fake_vs_phi_pt1%s.pdf", tag), -TMath::Pi(), TMath::Pi());

  DumpTable(f, "pt", "all pT");
  DumpTable(f, "eta2", "pT > 2 GeV, |eta| < 2.4");
  DumpTable(f, "eta1", "pT > 1 GeV, |eta| < 2.4");

  printf("\nfigures written to ../figures/\n");
}
