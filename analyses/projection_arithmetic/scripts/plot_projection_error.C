// ---------------------------------------------------------------------------
// plot_projection_error.C
//
// Why the tracklet projection fails at low pT. Three standalone figures, no
// data needed: everything is the helix and the arithmetic the emulation does.
//
// The emulation predicts the azimuth of a track at a layer radius r as
//
//     phi(r) = phi0 - asin(x),            x = r * rinv / 2 = r / 2R
//
// but an FPGA has no arcsine, so IMATH_TrackletCalculator.h:250-258 evaluates
// a third-order series, asin(x) -> x + x^3/6; and the stub's own radius
// differs from the layer's nominal radius by up to drmax = 3.75 cm, which
// MatchProcessor.cc:532 corrects with one linear step using the small-angle
// slope, dphi/dr -> -rinv/2 instead of -rinv/2 / sqrt(1-x^2).
//
// Both are excellent at 2 GeV, where x < 0.31 even at L6, and both fail as
// x -> 1, where the arcsine has a square-root singularity no polynomial can
// follow. See ../README.md for the equations.
//
// Usage:
//   root -l -b -q 'plot_projection_error.C()'            // all three
//   root -l -b -q 'plot_projection_error.C(1)'           // just figure 1
// ---------------------------------------------------------------------------

#include <cmath>

namespace {

// 0.01 * c * B = 0.0114 GeV cm^-1;  rinv [1/cm] = kCurv / pT [GeV]
const double kCurv = 0.01139536;
const double kRmean[6] = {24.9, 37.2, 52.3, 68.7, 86.0, 108.3};
// rphimatchcut_[layer][seed 0 = L1L2], converted cm -> mm
const double kWindow[6] = {0.0, 0.0, 1.0, 1.9, 4.0, 5.0};
const double kDrMax = 3.75;  // Settings::drmax()

double asin3(double x) { return x + x * x * x / 6.; }  // what the code computes

double ErrMM(double r, double pt) {  // position error on the layer, mm
  double x = 0.5 * r * (kCurv / pt);
  if (x >= 1.)
    return -1.;
  return std::abs(std::asin(x) - asin3(x)) * r * 10.;
}

void Save(TCanvas* c, const char* base) {
  c->SaveAs(Form("%s.pdf", base));
  c->SaveAs(Form("%s.png", base));
}

// ---------------------------------------------------------------------------
// Figure 1 -- the transverse plane.
//
// x has a simple geometric meaning: the turning circle has radius R = 1/rinv,
// so the track never gets further than 2R from the beam line, and
// x = r/(2R) is how far along toward that turning point a layer sits.
// ---------------------------------------------------------------------------
void FigGeometry(const char* out) {
  auto* c = new TCanvas("cgeo", "", 720, 700);
  c->SetLeftMargin(0.12);
  c->SetRightMargin(0.03);
  c->SetBottomMargin(0.11);
  c->SetTopMargin(0.08);

  auto* fr = gPad->DrawFrame(0, -12, 120, 116);
  fr->GetXaxis()->SetTitle("x [cm]");
  fr->GetYaxis()->SetTitle("y [cm]");
  fr->GetYaxis()->SetTitleOffset(1.25);

  for (int i = 0; i < 6; ++i) {
    auto* arc = new TArc(0, 0, kRmean[i], 0, 94);
    arc->SetFillStyle(0);
    arc->SetLineColor(kGray + 1);
    arc->SetLineStyle(3);
    arc->Draw("only");
    auto* lt = new TLatex(2.0, kRmean[i] + 1.8, Form("L%d", i + 1));
    lt->SetTextSize(0.028);
    lt->SetTextColor(kGray + 2);
    lt->Draw();
  }

  const double ptA = 0.5, rinvA = kCurv / ptA;
  const double Rturn = 1. / rinvA, reach = 2. * Rturn;

  auto* circ = new TArc(0, Rturn, Rturn);  // the full turning circle
  circ->SetFillStyle(0);
  circ->SetLineColor(kRed - 9);
  circ->SetLineStyle(2);
  circ->Draw("only");

  auto *gT = new TGraph(), *gA = new TGraph();
  for (double r = 0; r <= reach; r += 2.0) {
    double x = 0.5 * r * rinvA;
    if (x > 1.)
      break;
    gT->SetPoint(gT->GetN(), r * std::cos(std::asin(x)), r * std::sin(std::asin(x)));
    gA->SetPoint(gA->GetN(), r * std::cos(asin3(x)), r * std::sin(asin3(x)));
  }
  gT->SetLineColor(kRed + 1);
  gT->SetLineWidth(4);
  gT->Draw("L same");
  gA->SetLineColor(kRed + 1);
  gA->SetLineWidth(2);
  gA->SetLineStyle(11);
  gA->Draw("L same");

  const double rinvB = kCurv / 2.0;  // 2 GeV, for contrast
  auto *gT2 = new TGraph(), *gA2 = new TGraph();
  for (double r = 0; r <= 118; r += 3.0) {
    double x = 0.5 * r * rinvB;
    gT2->SetPoint(gT2->GetN(), r * std::cos(std::asin(x)), r * std::sin(std::asin(x)));
    gA2->SetPoint(gA2->GetN(), r * std::cos(asin3(x)), r * std::sin(asin3(x)));
  }
  gT2->SetLineColor(kAzure + 2);
  gT2->SetLineWidth(4);
  gT2->Draw("L same");
  gA2->SetLineColor(kAzure - 9);
  gA2->SetLineWidth(2);
  gA2->SetLineStyle(11);
  gA2->Draw("L same");

  auto *mT = new TGraph(), *mA = new TGraph();  // the layer crossings
  for (int L = 2; L < 6; ++L) {
    double x = 0.5 * kRmean[L] * rinvA;
    if (x >= 1.)
      continue;
    mT->SetPoint(mT->GetN(), kRmean[L] * std::cos(std::asin(x)), kRmean[L] * std::sin(std::asin(x)));
    mA->SetPoint(mA->GetN(), kRmean[L] * std::cos(asin3(x)), kRmean[L] * std::sin(asin3(x)));
  }
  mT->SetMarkerStyle(20);
  mT->SetMarkerSize(1.6);
  mT->SetMarkerColor(kRed + 1);
  mT->Draw("P same");
  mA->SetMarkerStyle(24);
  mA->SetMarkerSize(2.3);
  mA->SetMarkerColor(kRed + 1);
  mA->SetLineWidth(2);
  mA->Draw("P same");

  // phi0 is the track direction at the beam line; the projection swings away
  // from it by asin(x). Draw both, using L4 as the illustration.
  const double xL4 = 0.5 * kRmean[3] * rinvA, aL4 = std::asin(xL4);
  auto* ray = new TLine(0, 0, kRmean[3] * std::cos(aL4), kRmean[3] * std::sin(aL4));
  ray->SetLineColor(kGray + 2);
  ray->SetLineStyle(3);
  ray->Draw();
  auto* ar = new TArrow(0, 0, 30, 0, 0.018, "|>");
  ar->SetLineColor(kGray + 3);
  ar->SetLineWidth(2);
  ar->SetFillColor(kGray + 3);
  ar->Draw();
  auto* sweep = new TArc(0, 0, 21, 0, aL4 * TMath::RadToDeg());
  sweep->SetFillStyle(0);
  sweep->SetLineColor(kGray + 3);
  sweep->Draw("only");
  TLatex an;
  an.SetTextColor(kGray + 3);
  an.SetTextSize(0.034);
  an.DrawLatex(31, -1.5, "#phi_{0}");
  an.SetTextSize(0.030);
  an.DrawLatex(23.5, 11.0, "asin(x)");

  auto* leg = new TLegend(0.40, 0.715, 0.985, 0.915);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.09);
  leg->SetTextSize(0.0262);
  leg->SetHeader("x = r #upoint r_{inv}/2 = r / 2R");
  leg->AddEntry(gT, "0.5 GeV, true helix   #phi_{0} #minus asin(x)", "l");
  leg->AddEntry(gA, "0.5 GeV, as projected   #phi_{0} #minus (x + x^{3}/6)", "l");
  leg->AddEntry(gT2, "2 GeV, true helix", "l");
  leg->AddEntry(gA2, "2 GeV, as projected (coincides, 0.2 mm at L6)", "l");
  leg->Draw();

  Save(c, out);
}

// ---------------------------------------------------------------------------
// Figure 2 -- the same error in millimetres at the layer, against the real
// match windows for the L1L2 seed.
// ---------------------------------------------------------------------------
void FigVsPt(const char* out) {
  auto* c = new TCanvas("cpt", "", 760, 650);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);
  gPad->SetLogy();
  gPad->SetLogx();

  auto* fr = gPad->DrawFrame(0.35, 2e-3, 3.0, 3e3);
  fr->GetXaxis()->SetTitle("track p_{T} [GeV]");
  fr->GetYaxis()->SetTitle("projection error at the layer [mm]");
  fr->GetYaxis()->SetTitleOffset(1.30);
  fr->GetXaxis()->SetMoreLogLabels();
  fr->GetXaxis()->SetNoExponent();

  const int col[6] = {0, 0, kGreen + 2, kOrange + 7, kMagenta + 2, kRed + 1};
  auto* leg = new TLegend(0.45, 0.655, 0.97, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.09);
  leg->SetTextSize(0.029);
  leg->SetHeader("3rd-order asin(x) truncation, L1L2 seed");
  for (int L = 2; L < 6; ++L) {
    auto* g = new TGraph();
    for (double pt = 0.35; pt <= 3.0; pt *= 1.01) {
      double e = ErrMM(kRmean[L], pt);
      if (e > 0)
        g->SetPoint(g->GetN(), pt, e);
    }
    g->SetLineColor(col[L]);
    g->SetLineWidth(3);
    g->Draw("L same");

    // where this curve meets its own window: below that pT the layer stops
    // matching. One marker per layer is easier to read than four look-alike
    // horizontal lines bunched between 1 and 5 mm.
    double lo = 0.35, hi = 3.0;
    for (int it = 0; it < 60; ++it) {
      double mid = 0.5 * (lo + hi);
      double e = ErrMM(kRmean[L], mid);
      ((e < 0 || e > kWindow[L]) ? lo : hi) = mid;
    }
    auto* mk = new TGraph(1);
    mk->SetPoint(0, hi, kWindow[L]);
    mk->SetMarkerStyle(20);
    mk->SetMarkerSize(1.7);
    mk->SetMarkerColor(col[L]);
    mk->Draw("P same");

    leg->AddEntry(g, Form("L%d   window %.1f mm,  fails below %.2f GeV", L + 1, kWindow[L], hi), "l");
  }
  leg->Draw();

  Save(c, out);
}

// ---------------------------------------------------------------------------
// Figure 3 -- inside one layer. Two separate errors: the offset at the nominal
// radius (the series) and the tilt across the layer (the wrong slope).
// ---------------------------------------------------------------------------
void FigInLayer(const char* out, double pt = 0.7, int layer = 3) {
  auto* c = new TCanvas("clay", "", 760, 650);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.08);

  const double rm = kRmean[layer], rinv = kCurv / pt, x0 = 0.5 * rm * rinv;
  auto* fr = gPad->DrawFrame(rm - kDrMax - 0.4, -7, rm + kDrMax + 0.4, 11);
  fr->GetXaxis()->SetTitle("stub radius [cm]");
  fr->GetYaxis()->SetTitle("predicted #minus true position [mm]");
  fr->GetYaxis()->SetTitleOffset(1.30);

  const double derApprox = -0.5 * rinv;  // what the code uses

  auto *gBoth = new TGraph(), *gSlope = new TGraph();
  for (double r = rm - kDrMax; r <= rm + kDrMax; r += 0.75) {
    double truePhi = -std::asin(0.5 * r * rinv), dr = r - rm;
    gBoth->SetPoint(gBoth->GetN(), r, (-asin3(x0) + dr * derApprox - truePhi) * rm * 10.);
    gSlope->SetPoint(gSlope->GetN(), r, (-std::asin(x0) + dr * derApprox - truePhi) * rm * 10.);
  }

  auto* zero = new TLine(rm - kDrMax - 0.4, 0, rm + kDrMax + 0.4, 0);
  zero->SetLineColor(kGray + 2);
  zero->Draw();
  for (int s : {-1, 1}) {
    auto* w = new TLine(rm - kDrMax - 0.4, s * kWindow[layer], rm + kDrMax + 0.4, s * kWindow[layer]);
    w->SetLineStyle(2);
    w->SetLineColor(kGray + 3);
    w->Draw();
  }
  auto* vm = new TLine(rm, -7, rm, 11);
  vm->SetLineStyle(3);
  vm->SetLineColor(kGray + 2);
  vm->Draw();

  gSlope->SetLineColor(kAzure + 2);
  gSlope->SetLineWidth(3);
  gSlope->SetLineStyle(11);
  gSlope->Draw("L same");
  gBoth->SetLineColor(kRed + 1);
  gBoth->SetLineWidth(4);
  gBoth->Draw("L same");

  auto* leg = new TLegend(0.165, 0.775, 0.97, 0.895);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.09);
  leg->SetTextSize(0.029);
  leg->AddEntry(gBoth, "series offset + small-angle slope  (what the code does)", "l");
  leg->AddEntry(gSlope, "small-angle slope only  (series made exact)", "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.029);
  tx.SetTextColor(kGray + 2);
  tx.DrawLatex(0.545, 0.715, "r_{mean}");
  Save(c, out);
}

}  // namespace

void plot_projection_error(int which = 0) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  // ROOT's built-in dashed styles are too fine to read at these line widths
  gStyle->SetLineStyleString(11, "70 34");
  if (which == 0 || which == 1)
    FigGeometry("../figures/projection_geometry");
  if (which == 0 || which == 2)
    FigVsPt("../figures/projection_error_vs_pt");
  if (which == 0 || which == 3)
    FigInLayer("../figures/projection_error_in_layer");

  if (which == 0) {
    printf("\nprojection error [mm], 3rd order, against the window\n%8s", "pT");
    for (int L = 2; L < 6; ++L)
      printf("%14s", Form("L%d (%.1f mm)", L + 1, kWindow[L]));
    printf("\n");
    for (double pt : {2.0, 1.5, 1.0, 0.8, 0.6, 0.5, 0.4}) {
      printf("%8.1f", pt);
      for (int L = 2; L < 6; ++L) {
        double e = ErrMM(kRmean[L], pt);
        printf("%14s", e < 0 ? "--" : Form("%.2f", e));
      }
      printf("\n");
    }
  }
}
