// ---------------------------------------------------------------------------
// plot_projection_error.C
//
// Why the tracklet projection fails at low pT. Three panels, no data needed:
// everything here is the helix and the arithmetic the emulation actually does.
//
// The emulation predicts the azimuth of a track at a layer radius r as
//
//     phi(r) = phi0 - asin(x),            x = r * rinv / 2
//
// but an FPGA has no arcsine, so IMATH_TrackletCalculator.h:250-258 evaluates
// it as a third-order series,
//
//     asin(x)  ->  x + x^3/6                                   [PANEL 1, 2]
//
// and the stub's own radius differs from the layer's nominal radius by up to
// drmax = 3.75 cm, which MatchProcessor.cc:532 corrects with one linear step
// using the small-angle slope,
//
//     dphi/dr  ->  -rinv/2       instead of    -rinv/2 / sqrt(1-x^2)   [PANEL 3]
//
// Both are excellent at 2 GeV, where x < 0.31 even at L6, and both fail as
// x -> 1, where the arcsine has a square-root singularity that no polynomial
// can follow.
//
// Usage:  root -l -b -q 'plot_projection_error.C()'
// ---------------------------------------------------------------------------

#include <cmath>

namespace {

// 0.01 * c * B = 0.0114 GeV*cm^-1; rinv [1/cm] = kCurv / pT [GeV]
const double kCurv = 0.01139536;
const double kRmean[6] = {24.9, 37.2, 52.3, 68.7, 86.0, 108.3};
// rphimatchcut_[layer][seed 0 = L1L2], cm -> mm
const double kWindow[6] = {0.0, 0.0, 1.0, 1.9, 4.0, 5.0};
const double kDrMax = 3.75;  // Settings::drmax()

double asin3(double x) { return x + x * x * x / 6.; }                       // what the code does
double asin5(double x) { return asin3(x) + 3. * std::pow(x, 5) / 40.; }     // one more term
double asin7(double x) { return asin5(x) + 15. * std::pow(x, 7) / 336.; }   // two more

// Position error on the layer, in mm, from truncating the series.
double ErrMM(double r, double pt, int order) {
  double x = 0.5 * r * (kCurv / pt);
  if (x >= 1.)
    return -1.;
  double approx = (order == 3) ? asin3(x) : (order == 5) ? asin5(x) : asin7(x);
  return std::abs(std::asin(x) - approx) * r * 10.;
}

void Frame(TVirtualPad* p, double l, double r, double b, double t) {
  p->SetLeftMargin(l);
  p->SetRightMargin(r);
  p->SetBottomMargin(b);
  p->SetTopMargin(t);
}

}  // namespace

void plot_projection_error(const char* outBase = "../figures/projection_error") {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  auto* c = new TCanvas("c", "", 1500, 500);
  c->Divide(3, 1);

  // =======================================================================
  // PANEL 1 -- the transverse plane. The dashed curve is where the algorithm
  // thinks the track is: the locus of (r, x + x^3/6) instead of (r, asin(x)).
  // =======================================================================
  c->cd(1);
  Frame(gPad, 0.13, 0.04, 0.13, 0.09);
  auto* fr = gPad->DrawFrame(0, -8, 118, 118);
  fr->GetXaxis()->SetTitle("x [cm]");
  fr->GetYaxis()->SetTitle("y [cm]");
  fr->GetYaxis()->SetTitleOffset(1.35);

  for (int i = 0; i < 6; ++i) {  // the six barrel layers
    auto* arc = new TArc(0, 0, kRmean[i], 0, 92);
    arc->SetFillStyle(0);
    arc->SetLineColor(kGray + 1);
    arc->SetLineStyle(3);
    arc->Draw("only");
  }

  // For a helix from the origin the position azimuth at radius r is exactly
  // asin(x) -- the same quantity the projection has to compute. So the true
  // trajectory and the algorithm's belief can be drawn on the same axes.
  const int kNpt = 2;
  const double kPt[kNpt] = {2.0, 0.5};
  const int kCol[kNpt] = {kAzure + 2, kRed + 1};
  auto* leg = new TLegend(0.42, 0.63, 0.95, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.040);

  for (int k = 0; k < kNpt; ++k) {
    double rinv = kCurv / kPt[k];
    auto* gT = new TGraph();
    auto* gA = new TGraph();
    for (double r = 0; r <= 118; r += 0.25) {
      double x = 0.5 * r * rinv;
      if (x < 1.0) {
        double a = std::asin(x);
        gT->SetPoint(gT->GetN(), r * std::cos(a), r * std::sin(a));
      }
      double b = asin3(x);  // the series keeps returning a value past x = 1
      if (x < 1.0)
        gA->SetPoint(gA->GetN(), r * std::cos(b), r * std::sin(b));
    }
    gT->SetLineColor(kCol[k]);
    gT->SetLineWidth(3);
    gT->Draw("L same");
    gA->SetLineColor(kCol[k]);
    gA->SetLineWidth(2);
    gA->SetLineStyle(7);
    gA->Draw("L same");
    leg->AddEntry(gT, Form("p_{T} = %.1f GeV, true helix", kPt[k]), "l");
    leg->AddEntry(gA, Form("p_{T} = %.1f GeV, as projected", kPt[k]), "l");
  }
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.047);
  tx.DrawLatex(0.13, 0.945, "Where the algorithm thinks the track is");
  tx.SetTextSize(0.036);
  tx.SetTextColor(kGray + 3);
  tx.DrawLatex(0.17, 0.30, "dashed: #phi_{0} #minus (x + x^{3}/6)");
  tx.DrawLatex(0.17, 0.24, "solid:  #phi_{0} #minus asin(x)");
  tx.SetTextColor(kBlack);

  // =======================================================================
  // PANEL 2 -- the same error in millimetres at the layer, against the
  // actual match windows for the L1L2 seed.
  // =======================================================================
  c->cd(2);
  Frame(gPad, 0.13, 0.04, 0.13, 0.09);
  gPad->SetLogy();
  gPad->SetLogx();
  auto* fr2 = gPad->DrawFrame(0.35, 2e-3, 3.0, 3e3);
  fr2->GetXaxis()->SetTitle("track p_{T} [GeV]");
  fr2->GetYaxis()->SetTitle("projection error at the layer [mm]");
  fr2->GetYaxis()->SetTitleOffset(1.35);
  fr2->GetXaxis()->SetMoreLogLabels();
  fr2->GetXaxis()->SetNoExponent();

  const int kLcol[6] = {0, 0, kGreen + 2, kOrange + 7, kMagenta + 2, kRed + 1};
  auto* leg2 = new TLegend(0.17, 0.62, 0.58, 0.88);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  leg2->SetTextSize(0.038);
  for (int L = 2; L < 6; ++L) {
    auto* g = new TGraph();
    for (double pt = 0.35; pt <= 3.0; pt *= 1.01) {
      double e = ErrMM(kRmean[L], pt, 3);
      if (e > 0)
        g->SetPoint(g->GetN(), pt, e);
    }
    g->SetLineColor(kLcol[L]);
    g->SetLineWidth(3);
    g->Draw("L same");
    leg2->AddEntry(g, Form("L%d  (window %.1f mm)", L + 1, kWindow[L]), "l");
    auto* w = new TLine(0.35, kWindow[L], 3.0, kWindow[L]);
    w->SetLineColor(kLcol[L]);
    w->SetLineStyle(3);
    w->Draw();
  }
  leg2->Draw();
  tx.SetTextSize(0.047);
  tx.DrawLatex(0.13, 0.945, "Error vs the match window");
  tx.SetTextSize(0.034);
  tx.SetTextColor(kGray + 3);
  tx.DrawLatex(0.58, 0.30, "dotted: window");
  tx.DrawLatex(0.58, 0.245, "above it #Rightarrow no match");
  tx.SetTextColor(kBlack);

  // =======================================================================
  // PANEL 3 -- inside one layer. Two separate errors: the offset at the
  // nominal radius (the series) and the tilt across the layer (the slope).
  // =======================================================================
  c->cd(3);
  Frame(gPad, 0.13, 0.04, 0.13, 0.09);
  const double pt3 = 0.7, L3r = kRmean[3];  // L4 at 0.7 GeV: both terms matter
  const double rinv3 = kCurv / pt3;
  auto* fr3 = gPad->DrawFrame(L3r - kDrMax - 0.4, -7, L3r + kDrMax + 0.4, 11);
  fr3->GetXaxis()->SetTitle("stub radius [cm]");
  fr3->GetYaxis()->SetTitle("predicted #minus true position [mm]");
  fr3->GetYaxis()->SetTitleOffset(1.35);

  const double x0 = 0.5 * L3r * rinv3;
  const double derApprox = -0.5 * rinv3;                            // what the code uses
  const double derExact = -0.5 * rinv3 / std::sqrt(1 - x0 * x0);    // the true slope
  auto* gBoth = new TGraph();
  auto* gSlopeOnly = new TGraph();
  for (double r = L3r - kDrMax; r <= L3r + kDrMax; r += 0.05) {
    double truePhi = -std::asin(0.5 * r * rinv3);
    double dr = r - L3r;
    double both = -asin3(x0) + dr * derApprox;       // series offset + wrong slope
    double slope = -std::asin(x0) + dr * derApprox;  // exact at rmean, wrong slope only
    gBoth->SetPoint(gBoth->GetN(), r, (both - truePhi) * L3r * 10.);
    gSlopeOnly->SetPoint(gSlopeOnly->GetN(), r, (slope - truePhi) * L3r * 10.);
  }
  auto* zero = new TLine(L3r - kDrMax - 0.4, 0, L3r + kDrMax + 0.4, 0);
  zero->SetLineColor(kGray + 2);
  zero->Draw();
  for (int s : {-1, 1}) {  // the match window
    auto* w = new TLine(L3r - kDrMax - 0.4, s * kWindow[3], L3r + kDrMax + 0.4, s * kWindow[3]);
    w->SetLineStyle(2);
    w->SetLineColor(kGray + 3);
    w->Draw();
  }
  gSlopeOnly->SetLineColor(kAzure + 2);
  gSlopeOnly->SetLineWidth(3);
  gSlopeOnly->SetLineStyle(7);
  gSlopeOnly->Draw("L same");
  gBoth->SetLineColor(kRed + 1);
  gBoth->SetLineWidth(3);
  gBoth->Draw("L same");

  auto* leg3 = new TLegend(0.16, 0.16, 0.95, 0.33);
  leg3->SetBorderSize(0);
  leg3->SetFillStyle(0);
  leg3->SetTextSize(0.033);
  leg3->AddEntry(gBoth, "series offset + small-angle slope (the code)", "l");
  leg3->AddEntry(gSlopeOnly, "small-angle slope only (series made exact)", "l");
  leg3->Draw();
  tx.SetTextSize(0.047);
  tx.DrawLatex(0.13, 0.945, Form("Inside L4, p_{T} = %.1f GeV", pt3));
  tx.SetTextSize(0.034);
  tx.SetTextColor(kGray + 3);
  tx.DrawLatex(0.56, 0.44, "dashed: #pm window");
  tx.DrawLatex(0.56, 0.385, "offset at r_{mean} = series");
  tx.DrawLatex(0.56, 0.33, "tilt = wrong slope");

  c->SaveAs(Form("%s.pdf", outBase));
  c->SaveAs(Form("%s.png", outBase));

  printf("\nprojection error [mm], 3rd order, against the window\n");
  printf("%8s", "pT");
  for (int L = 2; L < 6; ++L)
    printf("%14s", Form("L%d (%.1f mm)", L + 1, kWindow[L]));
  printf("\n");
  for (double pt : {2.0, 1.5, 1.0, 0.8, 0.6, 0.5, 0.4}) {
    printf("%8.1f", pt);
    for (int L = 2; L < 6; ++L) {
      double e = ErrMM(kRmean[L], pt, 3);
      printf("%14s", e < 0 ? "--" : Form("%.2f", e));
    }
    printf("\n");
  }
}
