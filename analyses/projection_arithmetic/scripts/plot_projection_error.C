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
// Figure 3 -- one layer, with its real geometry, in the frame of the true track.
//
// L4 is built from 2S modules on ladders staggered between two radii so that
// they overlap in phi without gaps, and each module is a pair of sensors a few
// millimetres apart. The radii below are measured from the stub and cluster
// positions in an EPOS pPb ntuple (Extended2026D110):
//
//     inner ladder   sensors at 66.80 and 67.25 cm
//     outer ladder   sensors at 69.97 and 70.41 cm
//     rmean = 68.7 cm lies in the gap BETWEEN them
//
// so no stub is ever at rmean, and every stub needs the linear step -- by
// about +-1.7 cm. The track's azimuth sweeps tens of millimetres across that,
// far more than the 1.9 mm window, so the plot is drawn relative to the true
// track: the track is the line y = 0 with its stubs on it, and the question is
// whether the window around the PREDICTED position still contains them.
// ---------------------------------------------------------------------------
void FigInLayer(const char* out, double pt = 0.80, int layer = 3) {
  auto* c = new TCanvas("clay", "", 820, 650);
  c->SetLeftMargin(0.125);
  c->SetRightMargin(0.035);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.06);

  // measured sensor radii of L4; stubs sit at the inner sensor of each module
  const double kSensor[4] = {66.80, 67.25, 69.97, 70.41};
  const double kStub[2] = {66.80, 69.97};

  const double rm = kRmean[layer], rinv = kCurv / pt, x0 = 0.5 * rm * rinv;
  const double win = kWindow[layer], lo = 65.9, hi = 71.7;
  const double derApprox = -0.5 * rinv;  // the small-angle slope the code uses

  auto* fr = gPad->DrawFrame(lo, -4.3, hi, 5.3);
  fr->GetXaxis()->SetTitle("r  [cm]   (radially outward #rightarrow)");
  fr->GetYaxis()->SetTitle("azimuthal position, relative to the true track [mm]");
  fr->GetYaxis()->SetTitleOffset(1.25);

  // predicted position measured from the true track: series offset + slope*dr
  auto Pred = [&](double r) {
    return (-asin3(x0) + (r - rm) * derApprox + std::asin(0.5 * r * rinv)) * rm * 10.;
  };

  // Two ladders, staggered in radius so that neighbouring ladders overlap in
  // phi without gaps. Each is drawn as one module block with its two sensor
  // planes inside. A track crosses ONE ladder, or both where they overlap.
  for (int m = 0; m < 2; ++m) {
    auto* mod = new TBox(kSensor[2 * m] - 0.06, -2.45, kSensor[2 * m + 1] + 0.06, 5.3);
    mod->SetFillColorAlpha(kGray + 1, 0.22);
    mod->Draw();
    for (int i = 0; i < 2; ++i) {
      auto* sb = new TBox(kSensor[2 * m + i] - 0.022, -2.45, kSensor[2 * m + i] + 0.022, 5.3);
      sb->SetFillColorAlpha(kGray + 3, 0.65);
      sb->Draw();
    }
    auto* mb = new TBox(kSensor[2 * m] - 0.06, -2.78, kSensor[2 * m + 1] + 0.06, -2.58);
    mb->SetFillColorAlpha(kGray + 2, 0.6);
    mb->Draw();
  }

  auto* gWin = new TGraph();  // the match window, around where the code looks
  for (double r = lo; r <= hi; r += 0.25)
    gWin->SetPoint(gWin->GetN(), r, Pred(r) + win);
  for (double r = hi; r >= lo; r -= 0.25)
    gWin->SetPoint(gWin->GetN(), r, Pred(r) - win);
  gWin->SetFillColorAlpha(kOrange + 7, 0.25);
  gWin->Draw("F same");

  auto* gP = new TGraph();
  for (double r = lo; r <= hi; r += 0.58)
    gP->SetPoint(gP->GetN(), r, Pred(r));
  gP->SetLineColor(kOrange + 8);
  gP->SetLineWidth(3);
  gP->SetLineStyle(11);
  gP->Draw("L same");

  auto* trk = new TLine(lo, 0, hi, 0);  // the true track, by construction
  trk->SetLineColor(kRed + 1);
  trk->SetLineWidth(4);
  trk->Draw();

  auto* vm = new TLine(rm, -3.5, rm, 5.3);
  vm->SetLineStyle(3);
  vm->SetLineColor(kGray + 2);
  vm->Draw();

  auto *gOk = new TGraph(), *gBad = new TGraph();
  for (double r : kStub)
    ((std::abs(Pred(r)) < win) ? gOk : gBad)->SetPoint(((std::abs(Pred(r)) < win) ? gOk : gBad)->GetN(), r, 0.);
  gOk->SetMarkerStyle(20);
  gOk->SetMarkerSize(2.2);
  gOk->SetMarkerColor(kGreen + 3);
  gBad->SetMarkerStyle(20);
  gBad->SetMarkerSize(2.2);
  gBad->SetMarkerColor(kGray + 3);
  gOk->Draw("P same");
  gBad->Draw("P same");
  for (double r : kStub) {
    if (std::abs(Pred(r)) < win)
      continue;
    auto* xx = new TMarker(r, 0, 5);
    xx->SetMarkerColor(kRed + 2);
    xx->SetMarkerSize(2.6);
    xx->Draw();
  }

  TLatex an;
  an.SetTextSize(0.026);
  an.SetTextAlign(22);
  an.SetTextColor(kGray + 3);
  for (int m = 0; m < 2; ++m)
    an.DrawLatex(0.5 * (kSensor[2 * m] + kSensor[2 * m + 1]), -3.25,
                 m == 0 ? "inner ladder" : "outer ladder");
  an.SetTextSize(0.0235);
  an.DrawLatex(0.5 * (kSensor[0] + kSensor[1]), -3.72, "one module = 2 sensors");
  an.SetTextSize(0.026);
  an.SetTextColor(kGray + 2);
  an.DrawLatex(rm, -4.05, "r_{mean} : no stub here");
  an.SetTextSize(0.0245);
  an.SetTextColor(kGreen + 3);
  an.SetTextAlign(12);
  an.DrawLatex(kStub[0] + 0.12, -0.75, "inside the window");
  an.SetTextColor(kRed + 2);
  an.SetTextAlign(32);
  an.DrawLatex(kStub[1] - 0.12, -0.75, "outside: lost");

  auto* leg = new TLegend(0.125, 0.745, 0.86, 0.935);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.06);
  leg->SetTextSize(0.0255);
  leg->SetHeader(Form("L%d, p_{T} = %.2f GeV,  match window %.1f mm", layer + 1, pt, win));
  leg->AddEntry(trk, "true track; markers are its stub on each ladder", "l");
  leg->AddEntry(gP, "where the algorithm looks", "l");
  leg->AddEntry(gWin, "match window around it", "f");
  leg->Draw();

  Save(c, out);
}

// ---------------------------------------------------------------------------
// Figure 4 -- where 1/sqrt(1-x^2) comes from, and why it only matters at low pT.
//
// Differentiating phi(r) = phi0 - asin(x) with x = r*rinv/2 gives
//
//     dphi/dr = -(rinv/2) / sqrt(1-x^2)
//
// and the factor has a geometric meaning: asin(x) is the angle between the
// track and the radial direction where it crosses the layer, so
//
//     1/sqrt(1-x^2) = 1 / cos(crossing angle) = sec(crossing angle)
//
// The code keeps only -(rinv/2), i.e. it assumes the track crosses the layer
// radially. At 2 GeV the crossing angle at L4 is 11 deg and the secant is
// 1.02 -- a 2% error nobody would notice. At 0.5 GeV the angle is 52 deg and
// the secant is 1.61, so the assumed direction is 60% too shallow and the
// prediction drifts away across the thickness of the layer.
// ---------------------------------------------------------------------------
void FigSlope(const char* out, int layer = 3) {
  // equal x and y scales, so the crossing angle on the page is the real one
  auto* c = new TCanvas("cslp", "", 700, 700);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.04);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.06);

  const double rm = kRmean[layer], win = kWindow[layer];
  const double H = 5.6;
  auto* fr = gPad->DrawFrame(-H, -H, H, H);
  fr->GetXaxis()->SetTitle("r_{stub} #minus r_{mean}  [cm]   (radially outward #rightarrow)");
  fr->GetYaxis()->SetTitle("azimuthal displacement across the layer, r#Delta#phi [cm]");
  fr->GetYaxis()->SetTitleOffset(1.35);

  auto* band = new TBox(-kDrMax, -H, kDrMax, H);  // the thickness of one layer
  band->SetFillColorAlpha(kAzure + 1, 0.05);
  band->Draw();
  auto* rad = new TLine(-H, 0, 2.4, 0);  // the radial direction
  rad->SetLineColor(kGray + 2);
  rad->SetLineStyle(3);
  rad->Draw();
  auto* vz = new TLine(0, -H, 0, H);  // r = rmean
  vz->SetLineColor(kGray + 2);
  vz->Draw();

  const int kN = 2;
  const double pts[kN] = {2.0, 0.5};
  const int cols[kN] = {kAzure + 2, kRed + 1};
  auto* leg = new TLegend(0.145, 0.775, 0.97, 0.945);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetMargin(0.06);
  leg->SetTextSize(0.0255);
  leg->SetHeader(Form("L%d:  d#phi/dr = #minus(r_{inv}/2) / #sqrt{1#minus x^{2}} = #minus(r_{inv}/2) #upoint sec#theta", layer + 1));

  TLatex an;
  an.SetTextSize(0.026);

  for (int k = 0; k < kN; ++k) {
    const double rinv = kCurv / pts[k], x0 = 0.5 * rm * rinv;
    const double theta = std::asin(x0), secT = 1. / std::cos(theta);

    auto *gT = new TGraph(), *gA = new TGraph();
    for (double dr = -kDrMax; dr <= kDrMax + 1e-9; dr += 0.75) {
      double r = rm + dr;
      gT->SetPoint(gT->GetN(), dr, (std::asin(0.5 * r * rinv) - std::asin(x0)) * rm);
      gA->SetPoint(gA->GetN(), dr, (0.5 * rinv * dr) * rm);
    }
    gT->SetLineColor(cols[k]);
    gT->SetLineWidth(4);
    gT->Draw("L same");
    gA->SetLineColor(cols[k]);
    gA->SetLineWidth(3);
    gA->SetLineStyle(11);
    gA->Draw("L same");

    double gap = (std::asin(0.5 * (rm + kDrMax) * rinv) - std::asin(x0) - 0.5 * rinv * kDrMax) * rm * 10.;
    leg->AddEntry(gT, Form("%.1f GeV:  #theta = %.0f#circ,  sec#theta = %.2f,  misses by %.1f mm at the layer edge",
                           pts[k], theta * TMath::RadToDeg(), secT, gap), "l");

    // the crossing angle, measured from the radial direction
    double arad = (k == 0) ? 2.2 : 3.6;
    auto* arc = new TArc(0, 0, arad, 0, theta * TMath::RadToDeg());
    arc->SetFillStyle(0);
    arc->SetLineColor(cols[k]);
    arc->Draw("only");
    an.SetTextColor(cols[k]);
    an.SetTextAlign(12);
    double am = 0.5 * theta;
    an.DrawLatex((arad + 0.3) * std::cos(am), (arad + 0.3) * std::sin(am), Form("#theta = %.0f#circ", theta * TMath::RadToDeg()));
  }
  leg->Draw();

  // the match window, for scale
  auto* wb = new TBox(4.75, -win / 10. - 2.6, 5.05, win / 10. - 2.6);
  wb->SetFillColorAlpha(kOrange + 7, 0.7);
  wb->Draw();
  an.SetTextColor(kOrange + 8);
  an.SetTextAlign(32);
  an.DrawLatex(4.60, -2.0, Form("#pm%.1f mm", win));
  an.DrawLatex(4.60, -2.6, "window,");
  an.DrawLatex(4.60, -3.2, "to scale");

  an.SetTextColor(kGray + 3);
  an.SetTextAlign(12);
  an.DrawLatex(-5.25, -3.6, "solid: true track");
  an.DrawLatex(-5.25, -4.25, "dashed: the direction the code assumes");
  an.SetTextAlign(12);
  an.DrawLatex(-5.25, 0.42, "radial");
  an.SetTextAlign(22);
  an.SetTextColor(kAzure + 3);
  an.DrawLatex(0, -5.15, Form("one layer:  |r_{stub} #minus r_{mean}| #leq %.2f cm", kDrMax));

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
  if (which == 0 || which == 4)
    FigSlope("../figures/projection_slope");

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
