// ---------------------------------------------------------------------------
// plot_res.C -- momentum resolution plots from compare_res.C output.
//
//   1. 2D: pT(reco)-pT(truth) vs truth pT, one canvas per sample
//   2. profile mean  vs truth pT, both samples overlaid  (the bias)
//   3. profile RMS   vs truth pT, both samples overlaid  (the resolution)
//   4. relative resolution vs truth pT
//   5. resolution vs |eta| and vs phi, pT > 2 GeV
//   6. residual shape in a representative bin
//
// Usage: root -l -b -q 'plot_res.C("res_qed_mumu.root")'
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

const double kPtBins[] = {1.8, 2.0, 2.25, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
const int kNPt = sizeof(kPtBins) / sizeof(double) - 1;
const double kEtaBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int kNEta = sizeof(kEtaBins) / sizeof(double) - 1;
const int kNPhi = 12;

void Style(TGraphErrors* g, int color, int marker) {
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.15);
  g->SetLineWidth(2);
}

// Mean (or RMS) of each x-slice of the 2D distribution.
TGraphErrors* ProfileGraph(TFile* f, const char* smp, const char* set, bool wantRMS,
                           int color, int marker) {
  auto* h2 = (TH2D*)f->Get(Form("h2_dpt_vs_pt_%s_%s", smp, set));
  if (!h2) return nullptr;
  auto* g = new TGraphErrors();
  int p = 0;
  for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
    std::unique_ptr<TH1D> slice(h2->ProjectionY(Form("_py%d_%s_%s_%d", ix, smp, set, wantRMS),
                                                ix, ix));
    double n = slice->GetEntries();
    if (n < 50) continue;
    double v = wantRMS ? slice->GetRMS() : slice->GetMean();
    double e = wantRMS ? slice->GetRMSError() : slice->GetMeanError();
    double c = h2->GetXaxis()->GetBinCenter(ix);
    double w = 0.5 * h2->GetXaxis()->GetBinWidth(ix);
    g->SetPoint(p, c, v);
    g->SetPointError(p, w, e);
    ++p;
  }
  Style(g, color, marker);
  return g;
}

// 68% half-width of a residual histogram -- robust against tails.
double Q68(TH1D* h, double& err) {
  err = 0;
  if (!h || h->GetEntries() < 50) return 0;
  double q[2] = {0.1587, 0.8413}, v[2];
  h->GetQuantiles(2, v, q);
  double s = 0.5 * (v[1] - v[0]);
  err = s / std::sqrt(2.0 * h->GetEntries());
  return s;
}

TGraphErrors* BinnedGraph(TFile* f, const char* kind, const char* axis, const char* smp,
                          const char* set, const double* edges, int n, int color, int marker) {
  auto* g = new TGraphErrors();
  int p = 0;
  for (int i = 0; i < n; ++i) {
    auto* h = (TH1D*)f->Get(Form("%s_%s%d_%s_%s", kind, axis, i, smp, set));
    double err = 0, s = Q68(h, err);
    if (s <= 0) continue;
    g->SetPoint(p, 0.5 * (edges[i] + edges[i + 1]), s);
    g->SetPointError(p, 0.5 * (edges[i + 1] - edges[i]), err);
    ++p;
  }
  Style(g, color, marker);
  return g;
}

void Overlay(TGraphErrors* gd, TGraphErrors* gu, const char* xtitle, const char* ytitle,
             const char* title, const char* out, double xlo, double xhi, double ylo,
             double yhi, bool zeroLine) {
  auto* c = new TCanvas(Form("c_%s", out), "", 820, 620);
  c->SetGridy();
  c->SetLeftMargin(0.145);
  c->SetBottomMargin(0.13);
  auto* fr = c->DrawFrame(xlo, ylo, xhi, yhi);
  fr->GetXaxis()->SetTitle(xtitle);
  fr->GetYaxis()->SetTitle(ytitle);
  fr->GetXaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleOffset(1.5);
  if (zeroLine) {
    auto* l = new TLine(xlo, 0, xhi, 0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
  }
  if (gd) gd->Draw("P same");
  if (gu) gu->Draw("P same");
  auto* leg = new TLegend(0.18, 0.74, 0.58, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  if (gd) leg->AddEntry(gd, "Default Stub", "lp");
  if (gu) leg->AddEntry(gu, "Dummy Stub", "lp");
  leg->Draw();
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.18, 0.92, title);
  SaveBoth(c, out);
}

void Draw2D(TFile* f, const char* smp, const char* set, const char* label, const char* out) {
  auto* h2 = (TH2D*)f->Get(Form("h2_dpt_vs_pt_%s_%s", smp, set));
  if (!h2) return;
  auto* c = new TCanvas(Form("c2_%s_%s", smp, set), "", 820, 620);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.16);
  c->SetBottomMargin(0.13);
  c->SetLogz();
  h2->GetYaxis()->SetRangeUser(-0.35, 0.35);
  h2->GetXaxis()->SetTitle("truth p_{T} [GeV]");
  h2->GetYaxis()->SetTitle("p_{T}^{reco} - p_{T}^{truth} [GeV]");
  h2->GetYaxis()->SetTitleOffset(1.35);
  h2->Draw("colz");
  auto* l = new TLine(kPtBins[0], 0, 10.0, 0);
  l->SetLineStyle(2);
  l->SetLineColor(kGray + 3);
  l->Draw();
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.13, 0.92, label);
  SaveBoth(c, out);
}

}  // namespace

void plot_res(const char* fname = "res_qed_mumu.root", const char* tag = "") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    printf("[error] cannot open %s\n", fname);
    return;
  }

  const char* set = "common";  // identical particles on both sides
  const char* title1 = "QED #mu#mu, particles reconstructed in both samples";

  // 1. the 2D distributions
  Draw2D(f, "def", set, "Default Stub -- QED #mu#mu, common particles",
         Form("../figures/res2d_dpt_vs_pt_default%s.pdf", tag));
  Draw2D(f, "dum", set, "Dummy Stub -- QED #mu#mu, common particles",
         Form("../figures/res2d_dpt_vs_pt_dummy%s.pdf", tag));

  // 2. profile mean -- the bias
  Overlay(ProfileGraph(f, "def", set, false, kBlue + 1, 20),
          ProfileGraph(f, "dum", set, false, kRed + 1, 21), "truth p_{T} [GeV]",
          "#LTp_{T}^{reco} - p_{T}^{truth}#GT [GeV]", title1,
          Form("../figures/res_profile_mean_vs_pt%s.pdf", tag), 1.8, 10.0, -0.05, 0.12, true);

  // 3. profile RMS -- the resolution
  Overlay(ProfileGraph(f, "def", set, true, kBlue + 1, 20),
          ProfileGraph(f, "dum", set, true, kRed + 1, 21), "truth p_{T} [GeV]",
          "RMS(p_{T}^{reco} - p_{T}^{truth}) [GeV]", title1,
          Form("../figures/res_profile_rms_vs_pt%s.pdf", tag), 1.8, 10.0, 0.0, 0.30, false);

  // 4. relative resolution
  Overlay(BinnedGraph(f, "rel", "pt", "def", set, kPtBins, kNPt, kBlue + 1, 20),
          BinnedGraph(f, "rel", "pt", "dum", set, kPtBins, kNPt, kRed + 1, 21),
          "truth p_{T} [GeV]", "#sigma(p_{T})/p_{T}  (68% half-width)", title1,
          Form("../figures/res_rel_vs_pt%s.pdf", tag), 1.8, 10.0, 0.0, 0.030, false);

  // 5. angular dependence, above threshold
  const char* title2 = "QED #mu#mu, p_{T} > 2 GeV, particles reconstructed in both samples";
  Overlay(BinnedGraph(f, "rel", "eta", "def", set, kEtaBins, kNEta, kBlue + 1, 20),
          BinnedGraph(f, "rel", "eta", "dum", set, kEtaBins, kNEta, kRed + 1, 21),
          "truth |#eta|", "#sigma(p_{T})/p_{T}  (68% half-width)", title2,
          Form("../figures/res_rel_vs_eta%s.pdf", tag), 0.0, 2.4, 0.0, 0.030, false);

  double phiEdges[kNPhi + 1];
  for (int i = 0; i <= kNPhi; ++i) phiEdges[i] = -TMath::Pi() + i * (2 * TMath::Pi() / kNPhi);
  Overlay(BinnedGraph(f, "rel", "phi", "def", set, phiEdges, kNPhi, kBlue + 1, 20),
          BinnedGraph(f, "rel", "phi", "dum", set, phiEdges, kNPhi, kRed + 1, 21),
          "truth #phi", "#sigma(p_{T})/p_{T}  (68% half-width)", title2,
          Form("../figures/res_rel_vs_phi%s.pdf", tag), -TMath::Pi(), TMath::Pi(), 0.0, 0.030, false);

  // 6. residual shape in the 2.0-2.25 GeV bin
  auto* hd = (TH1D*)f->Get(Form("abs_pt1_def_%s", set));
  auto* hu = (TH1D*)f->Get(Form("abs_pt1_dum_%s", set));
  if (hd && hu) {
    auto* c = new TCanvas("c_shape", "", 820, 620);
    c->SetLeftMargin(0.145);
    c->SetBottomMargin(0.13);
    hd->Rebin(4); hu->Rebin(4);
    hd->SetLineColor(kBlue + 1); hd->SetLineWidth(2);
    hu->SetLineColor(kRed + 1);  hu->SetLineWidth(2);
    hd->GetXaxis()->SetRangeUser(-0.25, 0.25);
    hd->GetXaxis()->SetTitle("p_{T}^{reco} - p_{T}^{truth} [GeV]");
    hd->GetYaxis()->SetTitle("tracks");
    hd->GetYaxis()->SetTitleOffset(1.5);
    hd->Draw("hist");
    hu->Draw("hist same");
    auto* leg = new TLegend(0.62, 0.74, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hd, "Default Stub", "l");
    leg->AddEntry(hu, "Dummy Stub", "l");
    leg->Draw();
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.DrawLatex(0.18, 0.92, "QED #mu#mu, 2.00 < p_{T}^{truth} < 2.25 GeV, common particles");
    SaveBoth(c, Form("../figures/res_shape_2p0_2p25%s.pdf", tag));
  }

  // Numeric dump of the profile, since that is what gets quoted.
  printf("\n=== profile of pT(reco)-pT(truth) vs truth pT (common particles) ===\n");
  printf("%16s %10s %10s %10s %10s\n", "truth pT [GeV]", "mean(def)", "RMS(def)",
         "mean(dum)", "RMS(dum)");
  auto* h2d = (TH2D*)f->Get(Form("h2_dpt_vs_pt_def_%s", set));
  auto* h2u = (TH2D*)f->Get(Form("h2_dpt_vs_pt_dum_%s", set));
  if (h2d && h2u) {
    for (int ix = 1; ix <= h2d->GetNbinsX(); ++ix) {
      std::unique_ptr<TH1D> sd(h2d->ProjectionY(Form("_a%d", ix), ix, ix));
      std::unique_ptr<TH1D> su(h2u->ProjectionY(Form("_b%d", ix), ix, ix));
      if (sd->GetEntries() < 50) continue;
      printf("%6.2f-%9.2f %10.4f %10.4f %10.4f %10.4f\n", h2d->GetXaxis()->GetBinLowEdge(ix),
             h2d->GetXaxis()->GetBinUpEdge(ix), sd->GetMean(), sd->GetRMS(), su->GetMean(),
             su->GetRMS());
    }
  }
  printf("\nfigures written to ../figures/\n");
}
