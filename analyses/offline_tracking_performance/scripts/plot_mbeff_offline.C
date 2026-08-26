// ---------------------------------------------------------------------------
// plot_mbeff_offline.C -- draw the event trigger efficiency produced by
// mbeff_offline.C: truth vs offline highPurity reco, against the pT threshold.
//
// The gap between the two curves is what offline reconstruction costs you
// relative to a perfect tracker; the truth curve itself is the ceiling.
//
// Usage:
//   root -l -b -q 'plot_mbeff_offline.C("../output/mbeff_hydjet.root","_hydjet","HYDJET PbPb")'
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

}  // namespace

void plot_mbeff_offline(const char* fname = "../output/mbeff.root",
                        const char* tag = "",
                        const char* sample = "HYDJET PbPb",
                        bool logx = false) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("[error] cannot open %s\n", fname); return; }

  auto* c = new TCanvas("cmb", "", 800, 650);
  c->SetGridy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
  if (logx) c->SetLogx();

  auto* frame = c->DrawFrame(0.3, 0.0, 5.0, 1.15);
  auto* one = new TLine(0.3, 1.0, 5.0, 1.0);
  one->SetLineStyle(3); one->SetLineColor(kGray + 1); one->Draw();
  frame->GetXaxis()->SetTitle("p_{T} threshold [GeV]");
  frame->GetYaxis()->SetTitle("event trigger efficiency");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.35);
  if (logx) { frame->GetXaxis()->SetMoreLogLabels(); frame->GetXaxis()->SetNoExponent(); }

  auto* gT = MakeRate(f, "numT_mb", "den_mb", kBlack, 20);
  auto* gR = MakeRate(f, "numR_mb", "den_mb", kRed + 1, 21);
  if (gT) gT->Draw("P same");
  if (gR) gR->Draw("P same");

  auto* leg = new TLegend(0.17, 0.18, 0.58, 0.36);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.034);
  if (gT) leg->AddEntry(gT, "truth (TrackingParticles)", "lp");
  if (gR) leg->AddEntry(gR, "offline reco (highPurity)", "lp");
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
  tx.DrawLatex(0.17, 0.94, Form("%s -- at least one track, |#eta| < 2.4", sample));

  SaveBoth(c, Form("../figures/mbeff_vs_ptmin%s.pdf", tag));

  // text table
  TH1D* hd = (TH1D*)f->Get("den_mb");
  TH1D* hT = (TH1D*)f->Get("numT_mb");
  TH1D* hR = (TH1D*)f->Get("numR_mb");
  if (hd && hT && hR) {
    printf("\n%10s %14s %14s %10s\n", "pTmin", "truth", "reco(hp)", "reco/truth");
    for (int b = 1; b <= hd->GetNbinsX(); ++b) {
      double d = hd->GetBinContent(b);
      if (d < 1) continue;
      double t = hd->GetBinCenter(b);
      double a = hT->GetBinContent(b) / d, r = hR->GetBinContent(b) / d;
      printf("%10.2f %14.4f %14.4f %10.4f\n", t, a, r, a > 0 ? r / a : 0);
    }
  }
  printf("\nfigure written to ../figures/mbeff_vs_ptmin%s.pdf\n", tag);
}
