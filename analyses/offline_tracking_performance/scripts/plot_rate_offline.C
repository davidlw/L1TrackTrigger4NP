// ---------------------------------------------------------------------------
// plot_rate_offline.C -- trigger rate from the histograms made by mbeff_offline.C.
//
//   rate = (total collision rate) x (fraction of events passing the requirement)
//
// Two plots per sample:
//   rate_vs_ptmin   require >= 1 object with |eta| < 2.4 above a pT threshold
//   rate_vs_nchmin  require Nch >= N, counting objects with pT > 0.4, |eta| < 2.4
//
// Each shows two curves:
//   truth   TrackingParticles (primaries only for the Nch plot) -- the ceiling
//   reco    offline highPurity tracks -- what a trigger could actually select on
//
// The total rate is an ASSUMPTION supplied by the caller, not something measured
// from the sample; the shape comes from the MC, the normalisation from you.
// Both axes are logarithmic in rate, since a threshold scan covers decades.
//
// Usage:
//   root -l -b -q 'plot_rate_offline.C("../output/mbeff_hydjet.root","_hydjet","HYDJET PbPb",50e3)'
//   root -l -b -q 'plot_rate_offline.C("../output/mbeff_epospb.root","_epospb","EPOS pPb",3e6)'
// ---------------------------------------------------------------------------

#include <cstdio>
#include <vector>

namespace {

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

// rate in kHz at each bin, from a numerator/denominator pair
TGraph* RateFromRatio(TH1D* hn, TH1D* hd, double totHz, int color, int marker) {
  auto* g = new TGraph();
  int k = 0;
  for (int b = 1; b <= hd->GetNbinsX(); ++b) {
    double d = hd->GetBinContent(b);
    if (d < 1) continue;
    double r = hn->GetBinContent(b) / d * totHz / 1000.0;   // kHz
    if (r <= 0) continue;                                   // log axis
    g->SetPoint(k++, hd->GetBinCenter(b), r);
  }
  g->SetLineColor(color); g->SetMarkerColor(color);
  g->SetMarkerStyle(marker); g->SetMarkerSize(1.0); g->SetLineWidth(2);
  return g;
}

// rate in kHz for "Nch >= bin low edge", i.e. the survival function of a
// per-event multiplicity distribution
TGraph* RateFromCumulative(TH1D* h, double nev, double totHz, int color, int marker) {
  auto* g = new TGraph();
  int k = 0;
  const int nb = h->GetNbinsX();
  for (int b = 1; b <= nb; ++b) {
    double surv = h->Integral(b, nb + 1);      // include overflow
    if (surv <= 0) continue;
    g->SetPoint(k++, h->GetBinLowEdge(b), surv / nev * totHz / 1000.0);
  }
  g->SetLineColor(color); g->SetMarkerColor(color);
  g->SetMarkerStyle(marker); g->SetMarkerSize(1.0); g->SetLineWidth(2);
  return g;
}

void Draw(TGraph* gT, TGraph* gR, const char* xtitle, const char* title,
          const char* outfile, double xlo, double xhi, double totHz,
          const char* legT, const char* legR) {
  auto* c = new TCanvas(Form("c%s", outfile), "", 800, 650);
  c->SetLogy();
  c->SetGridy();
  c->SetLeftMargin(0.14);
  c->SetBottomMargin(0.13);

  const double top = totHz / 1000.0 * 3.0;
  const double bot = totHz / 1000.0 * 1e-5;
  auto* frame = c->DrawFrame(xlo, bot, xhi, top);
  frame->GetXaxis()->SetTitle(xtitle);
  frame->GetYaxis()->SetTitle("trigger rate [kHz]");
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.45);

  // the total collision rate: nothing can sit above this
  auto* tot = new TLine(xlo, totHz / 1000.0, xhi, totHz / 1000.0);
  tot->SetLineStyle(2); tot->SetLineColor(kGray + 2); tot->Draw();

  if (gT) gT->Draw("PL same");
  if (gR) gR->Draw("PL same");

  auto* leg = new TLegend(0.50, 0.72, 0.90, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
  if (gT) leg->AddEntry(gT, legT, "lp");
  if (gR) leg->AddEntry(gR, legR, "lp");
  leg->AddEntry(tot, Form("total rate = %.3g kHz", totHz / 1000.0), "l");
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
  tx.DrawLatex(0.16, 0.94, title);
  SaveBoth(c, outfile);
}

}  // namespace

void plot_rate_offline(const char* fname, const char* tag, const char* sample,
                       double totalRateHz, double nchAxisMax = 0) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("[error] cannot open %s\n", fname); return; }

  TH1D* den  = (TH1D*)f->Get("den_mb");
  TH1D* numT = (TH1D*)f->Get("numT_mb");
  TH1D* numR = (TH1D*)f->Get("numR_mb");
  TH1D* nchT = (TH1D*)f->Get("nch_truth");
  TH1D* nchR = (TH1D*)f->Get("nch_reco");
  if (!den || !numT || !numR) { printf("[error] missing mb histograms\n"); return; }

  // ---------------- rate vs pT threshold ----------------
  Draw(RateFromRatio(numT, den, totalRateHz, kBlack, 20),
       RateFromRatio(numR, den, totalRateHz, kRed + 1, 21),
       "p_{T} threshold [GeV]",
       Form("%s -- at least one track, |#eta| < 2.4", sample),
       Form("../figures/rate_vs_ptmin%s.pdf", tag),
       den->GetXaxis()->GetXmin(), den->GetXaxis()->GetXmax(), totalRateHz,
       "truth (primary particles)", "offline reco (highPurity)");

  // ---------------- rate vs Nch threshold ----------------
  if (nchT && nchR) {
    const double nev = den->GetBinContent(1);   // every bin holds all events
    double xhi = nchAxisMax > 0 ? nchAxisMax : nchT->GetXaxis()->GetXmax();
    Draw(RateFromCumulative(nchT, nev, totalRateHz, kBlack, 20),
         RateFromCumulative(nchR, nev, totalRateHz, kRed + 1, 21),
         "N_{ch} threshold (p_{T} > 0.4, |#eta| < 2.4)",
         Form("%s -- N_{ch} #geq threshold", sample),
         Form("../figures/rate_vs_nchmin%s.pdf", tag),
         0, xhi, totalRateHz,
         "truth (primary particles)", "offline reco (highPurity)");

    printf("\n%s: %.0f events, total rate %.3g kHz\n", sample, nev, totalRateHz / 1000.);
    printf("%10s %14s %14s\n", "Nch >=", "truth [kHz]", "reco [kHz]");
    for (int b = 1; b <= nchT->GetNbinsX(); ++b) {
      double sT = nchT->Integral(b, nchT->GetNbinsX() + 1);
      double sR = nchR->Integral(b, nchR->GetNbinsX() + 1);
      if (sT < 1 && sR < 1) continue;
      if (b % 5 && b != 1) continue;
      printf("%10.0f %14.4g %14.4g\n", nchT->GetBinLowEdge(b),
             sT / nev * totalRateHz / 1000., sR / nev * totalRateHz / 1000.);
    }
  }

  printf("\n%10s %14s %14s\n", "pTmin", "truth [kHz]", "reco [kHz]");
  for (int b = 1; b <= den->GetNbinsX(); ++b) {
    double d = den->GetBinContent(b);
    if (d < 1) continue;
    if (b % 5 && b != 1) continue;
    printf("%10.2f %14.4g %14.4g\n", den->GetBinCenter(b),
           numT->GetBinContent(b) / d * totalRateHz / 1000.,
           numR->GetBinContent(b) / d * totalRateHz / 1000.);
  }
  printf("\nfigures written to ../figures/\n");
}
