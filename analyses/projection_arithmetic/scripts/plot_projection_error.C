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
  auto* c = new TCanvas("clay", "", 1250, 660);
  auto* p1 = new TPad("p1", "", 0.000, 0.155, 0.505, 1.0);
  auto* p2 = new TPad("p2", "", 0.505, 0.155, 1.000, 1.0);
  auto* p3 = new TPad("p3", "", 0.000, 0.000, 1.000, 0.155);
  p1->Draw();
  p2->Draw();
  p3->Draw();

  // measured sensor radii of L4; a stub sits at the inner sensor of its module
  const double kSensor[4] = {66.80, 67.25, 69.97, 70.41};
  const double kStub[2] = {66.80, 69.97};

  const double rm = kRmean[layer], rinv = kCurv / pt, x0 = 0.5 * rm * rinv;
  const double win = kWindow[layer], lo = 65.9, hi = 71.7;
  const double derApprox = -0.5 * rinv;  // the small-angle slope the code uses

  // predicted position measured from the true track: series offset + slope*dr
  auto Pred = [&](double r) {
    return (-asin3(x0) + (r - rm) * derApprox + std::asin(0.5 * r * rinv)) * rm * 10.;
  };

  // the error grows steeply as pT falls, so scale the frame with it and place
  // the labels as fractions of the range rather than at fixed millimetres
  const double yTop = std::max(5.3, Pred(hi) + win + 0.6);
  const double yBot = -0.45 * yTop;

  // Two tracks at different azimuth, each crossing one ladder. Each panel is
  // drawn in the frame of its own track, so that track is the line y = 0.
  for (int which = 0; which < 2; ++which) {
    (which == 0 ? p1 : p2)->cd();
    gPad->SetLeftMargin(which == 0 ? 0.155 : 0.090);
    gPad->SetRightMargin(0.035);
    gPad->SetBottomMargin(0.125);
    gPad->SetTopMargin(0.085);

    auto* fr = gPad->DrawFrame(lo, yBot, hi, yTop);
    fr->GetXaxis()->SetTitle("r  [cm]   (radially outward #rightarrow)");
    fr->GetXaxis()->SetTitleSize(0.047);
    fr->GetXaxis()->SetLabelSize(0.042);
    fr->GetYaxis()->SetTitleSize(0.047);
    fr->GetYaxis()->SetLabelSize(0.042);
    if (which == 0) {
      fr->GetYaxis()->SetTitle("azimuthal position, relative to this track [mm]");
      fr->GetYaxis()->SetTitleOffset(1.25);
    }

    for (int m = 0; m < 2; ++m) {  // both ladders; the crossed one is solid
      const bool hit = (m == which);
      auto* mod = new TBox(kSensor[2 * m] - 0.06, 0.62 * yBot, kSensor[2 * m + 1] + 0.06, yTop);
      mod->SetFillColorAlpha(kGray + 1, hit ? 0.30 : 0.10);
      mod->Draw();
      for (int i = 0; i < 2; ++i) {
        auto* sb = new TBox(kSensor[2 * m + i] - 0.022, 0.62 * yBot, kSensor[2 * m + i] + 0.022, yTop);
        sb->SetFillColorAlpha(kGray + 3, hit ? 0.75 : 0.22);
        sb->Draw();
      }
      auto* mb = new TBox(kSensor[2 * m] - 0.06, 0.705 * yBot, kSensor[2 * m + 1] + 0.06, 0.655 * yBot);
      mb->SetFillColorAlpha(kGray + 2, hit ? 0.7 : 0.22);
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

    auto* trk = new TLine(lo, 0, hi, 0);
    trk->SetLineColor(kRed + 1);
    trk->SetLineWidth(4);
    trk->Draw();

    auto* vm = new TLine(rm, 0.88 * yBot, rm, yTop);
    vm->SetLineStyle(3);
    vm->SetLineColor(kGray + 2);
    vm->Draw();

    const double rs = kStub[which];
    const bool ok = std::abs(Pred(rs)) < win;
    auto* mk = new TMarker(rs, 0, 20);
    mk->SetMarkerSize(2.3);
    mk->SetMarkerColor(ok ? kGreen + 3 : kGray + 3);
    mk->Draw();
    if (!ok) {
      auto* xx = new TMarker(rs, 0, 5);
      xx->SetMarkerColor(kRed + 2);
      xx->SetMarkerSize(2.7);
      xx->Draw();
    }
    // how far the stub is from where the algorithm looked
    auto* amiss = new TArrow(rs, 0, rs, Pred(rs), 0.013, "<|>");
    amiss->SetLineColor(ok ? kGreen + 3 : kRed + 2);
    amiss->SetFillColor(ok ? kGreen + 3 : kRed + 2);
    amiss->SetLineWidth(2);
    amiss->Draw();

    TLatex an;
    an.SetTextSize(0.040);
    an.SetTextAlign(22);
    an.SetTextColor(kGray + 3);
    for (int m = 0; m < 2; ++m)
      an.DrawLatex(0.5 * (kSensor[2 * m] + kSensor[2 * m + 1]), 0.80 * yBot, m == 0 ? "inner ladder" : "outer ladder");
    an.SetTextColor(kGray + 2);
    an.DrawLatex(rm, 0.95 * yBot, "r_{mean}");
    an.SetTextSize(0.043);
    an.SetTextColor(ok ? kGreen + 3 : kRed + 2);
    an.SetTextAlign(which == 0 ? 12 : 32);
    an.DrawLatex(rs + (which == 0 ? 0.18 : -0.18), 0.5 * Pred(rs),
                 ok ? Form("%.1f mm: matched", std::abs(Pred(rs))) : Form("%.1f mm: lost", std::abs(Pred(rs))));
    an.SetTextSize(0.044);
    an.SetTextColor(kBlack);
    an.SetTextAlign(12);
    an.DrawLatex(lo + 0.15, 0.90 * yTop,
                 which == 0 ? "track A, crossing the inner ladder" : "track B, crossing the outer ladder");
    an.SetTextSize(0.036);
    an.SetTextColor(kGray + 3);
    an.DrawLatex(lo + 0.15, 0.785 * yTop, Form("dr = %+.2f cm", rs - rm));

  }

  // one shared legend, in its own strip so nothing overlaps the panels
  p3->cd();
  auto* lt = new TLine();
  lt->SetLineColor(kRed + 1);
  lt->SetLineWidth(4);
  auto* lp = new TLine();
  lp->SetLineColor(kOrange + 8);
  lp->SetLineWidth(3);
  lp->SetLineStyle(11);
  auto* lw = new TBox();
  lw->SetFillColorAlpha(kOrange + 7, 0.25);
  auto* leg = new TLegend(0.03, 0.08, 0.99, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->SetMargin(0.12);
  leg->SetTextSize(0.26);
  leg->AddEntry(lt, "the track, with its stub", "l");
  leg->AddEntry(lp, "where the algorithm looks", "l");
  leg->AddEntry(lw, Form("#pm%.1f mm match window", kWindow[layer]), "f");
  leg->Draw();

  c->cd();
  TLatex hd;
  hd.SetNDC();
  hd.SetTextSize(0.030);
  hd.SetTextAlign(32);
  hd.DrawLatex(0.965, 0.965, Form("L%d, p_{T} = %.2f GeV,  layer geometry measured from the ntuple", layer + 1, pt));

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

// ---------------------------------------------------------------------------
// Figure 5 -- the canonical transverse view, zoomed onto one patch of L4.
//
// Radius runs upward, azimuth across. The ladders alternate between two radii
// around phi, so at any given azimuth there is EITHER an inner-ladder module
// OR an outer-ladder one; a track crosses one of them. Both are drawn here,
// side by side in phi, with a track through each.
//
// For each track: the true helix, the straight line the algorithm predicts
// (series offset at rmean, then the small-angle slope), the stub where the
// true track meets the sensor, and the match window around the prediction at
// that radius. The gap between them is what decides whether the stub is found.
//
// At these momenta the errors are centimetres, so for once everything can be
// drawn to scale in the plane.
// ---------------------------------------------------------------------------
void FigTransverse(const char* out, double pt = 0.60, int layer = 3) {
  auto* c = new TCanvas("ctr", "", 1200, 700);
  c->SetLeftMargin(0.075);
  c->SetRightMargin(0.02);
  c->SetBottomMargin(0.115);
  c->SetTopMargin(0.165);

  const double kSensor[4] = {66.80, 67.25, 69.97, 70.41};
  const double rm = kRmean[layer], rinv = kCurv / pt, x0 = 0.5 * rm * rinv;
  const double win = kWindow[layer] / 10.;  // cm

  // azimuthal arc of the TRUE track, measured from its own crossing of rmean
  auto Vtrue = [&](double r) { return rm * (std::asin(x0) - std::asin(0.5 * r * rinv)); };
  // ... and of the straight line the algorithm predicts
  auto Vpred = [&](double r) { return rm * (std::asin(x0) - asin3(x0) - (r - rm) * 0.5 * rinv); };

  const double vLo = -4.2, vHi = 9.3, rLo = 65.5, rHi = 71.7;
  auto* fr = gPad->DrawFrame(vLo, rLo, vHi, rHi);
  fr->GetXaxis()->SetTitle("azimuthal arc  r#Delta#phi  [cm]");
  fr->GetYaxis()->SetTitle("r  [cm]");
  fr->GetXaxis()->SetTitleSize(0.046);
  fr->GetYaxis()->SetTitleSize(0.046);
  fr->GetXaxis()->SetLabelSize(0.040);
  fr->GetYaxis()->SetLabelSize(0.040);
  fr->GetYaxis()->SetTitleOffset(0.72);

  // the two ladders, adjacent in phi, staggered in radius
  const double vSplit = 2.0;
  auto drawModule = [&](double r0, double r1, double va, double vb) {
    auto* box = new TBox(va, r0, vb, r1);
    box->SetFillColorAlpha(kGray + 1, 0.30);
    box->Draw();
    for (double rr : {r0, r1}) {
      auto* sl = new TLine(va, rr, vb, rr);
      sl->SetLineColor(kGray + 3);
      sl->SetLineWidth(3);
      sl->Draw();
    }
  };
  drawModule(kSensor[0], kSensor[1], vLo, vSplit + 0.25);   // inner ladder
  drawModule(kSensor[2], kSensor[3], vSplit - 0.25, vHi);   // outer ladder

  auto* vmean = new TLine(vLo, rm, vHi, rm);
  vmean->SetLineStyle(3);
  vmean->SetLineColor(kGray + 2);
  vmean->Draw();

  TLatex an;
  an.SetTextSize(0.034);

  // two tracks: A crosses the inner ladder, B the outer one
  const double vHitA = -0.6, vHitB = 5.6;
  const double rHitA = kSensor[0], rHitB = kSensor[2];
  for (int k = 0; k < 2; ++k) {
    const double rHit = (k == 0) ? rHitA : rHitB;
    const double shift = ((k == 0) ? vHitA : vHitB) - Vtrue(rHit);

    auto* gT = new TGraph();
    for (double r = rLo; r <= rHi; r += 0.25)
      gT->SetPoint(gT->GetN(), Vtrue(r) + shift, r);
    gT->SetLineColor(kRed + 1);
    gT->SetLineWidth(4);
    gT->Draw("L same");

    auto* gP = new TGraph();
    for (double r = rLo; r <= rHi; r += 0.62)
      gP->SetPoint(gP->GetN(), Vpred(r) + shift, r);
    gP->SetLineColor(kOrange + 8);
    gP->SetLineWidth(3);
    gP->SetLineStyle(11);
    gP->Draw("L same");

    // the match window, at the stub's radius, centred on the prediction
    const double vP = Vpred(rHit) + shift;
    auto* wb = new TBox(vP - win, rHit - 0.07, vP + win, rHit + 0.07);
    wb->SetFillColorAlpha(kOrange + 7, 0.75);
    wb->Draw();

    const double miss = std::abs(vP - (Vtrue(rHit) + shift));
    auto* mk = new TMarker(Vtrue(rHit) + shift, rHit, 20);
    mk->SetMarkerSize(2.0);
    mk->SetMarkerColor(miss < win ? kGreen + 3 : kGray + 3);
    mk->Draw();
    if (miss >= win) {
      auto* xx = new TMarker(Vtrue(rHit) + shift, rHit, 5);
      xx->SetMarkerColor(kRed + 2);
      xx->SetMarkerSize(2.4);
      xx->Draw();
    }
    auto* ar = new TArrow(Vtrue(rHit) + shift, rHit, vP, rHit, 0.012, "<|>");
    ar->SetLineColor(miss < win ? kGreen + 3 : kRed + 2);
    ar->SetFillColor(miss < win ? kGreen + 3 : kRed + 2);
    ar->SetLineWidth(2);
    ar->Draw();

    an.SetTextColor(miss < win ? kGreen + 3 : kRed + 2);
    an.SetTextAlign(23);
    an.DrawLatex(0.5 * (Vtrue(rHit) + shift + vP), rHit - 0.42,
                 Form("%.1f mm %s", miss * 10., miss < win ? "- matched" : "- lost"));
    // the line is anchored at rmean -- a radius with no module on it -- and is
    // already displaced there by the series truncation; from that point it is
    // extended with the small-angle slope.
    auto* anch = new TMarker(Vpred(rm) + shift, rm, 21);
    anch->SetMarkerColor(kOrange + 8);
    anch->SetMarkerSize(1.6);
    anch->Draw();
    if (k == 0) {
      auto* aser = new TArrow(Vtrue(rm) + shift, rm, Vpred(rm) + shift, rm, 0.011, "<|>");
      aser->SetLineColor(kOrange + 9);
      aser->SetFillColor(kOrange + 9);
      aser->SetLineWidth(2);
      aser->Draw();
      an.SetTextColor(kOrange + 9);
      an.SetTextAlign(23);
      an.DrawLatex(0.5 * (Vtrue(rm) + Vpred(rm)) + shift, rm - 0.18,
                   Form("series error, %.1f mm", (Vpred(rm) - Vtrue(rm)) * 10.));
      an.SetTextAlign(12);
      an.DrawLatex(Vpred(rm) + shift + 0.22, rm + 0.42, "line starts here, at r_{mean}");
    }

    an.SetTextColor(kRed + 1);
    an.SetTextAlign(k == 0 ? 32 : 32);
    an.DrawLatex(Vtrue(rLo + 0.25) + shift - 0.22, rLo + 0.30, k == 0 ? "track A" : "track B");
  }

  an.SetTextColor(kGray + 3);
  an.SetTextAlign(12);
  an.DrawLatex(vLo + 0.25, 0.5 * (kSensor[0] + kSensor[1]) - 0.75, "inner ladder");
  an.SetTextAlign(32);
  an.DrawLatex(vHi - 0.25, 0.5 * (kSensor[2] + kSensor[3]) + 0.75, "outer ladder");
  an.SetTextColor(kGray + 2);
  an.SetTextAlign(12);
  an.DrawLatex(vLo + 0.25, rm + 0.25, "r_{mean}");

  auto* leg = new TLegend(0.075, 0.845, 0.98, 0.925);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->SetMargin(0.10);
  leg->SetTextSize(0.0295);
  auto* lt = new TLine();
  lt->SetLineColor(kRed + 1);
  lt->SetLineWidth(4);
  auto* lp = new TLine();
  lp->SetLineColor(kOrange + 8);
  lp->SetLineWidth(3);
  lp->SetLineStyle(11);
  auto* lw = new TBox();
  lw->SetFillColorAlpha(kOrange + 7, 0.75);
  leg->AddEntry(lt, "true helix, and its stub", "l");
  leg->AddEntry(lp, "the straight line the algorithm predicts", "l");
  leg->AddEntry(lw, Form("#pm%.1f mm match window", kWindow[layer]), "f");
  leg->Draw();

  TLatex hd;
  hd.SetNDC();
  hd.SetTextSize(0.034);
  hd.SetTextAlign(12);
  hd.DrawLatex(0.075, 0.955, Form("L%d, p_{T} = %.2f GeV #minus transverse plane, drawn to scale", layer + 1, pt));

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
  if (which == 0 || which == 3) {
    FigInLayer("../figures/projection_error_in_layer", 0.80);
    FigInLayer("../figures/projection_error_in_layer_0p6", 0.60);
  }
  if (which == 0 || which == 4)
    FigSlope("../figures/projection_slope");
  if (which == 0 || which == 5) {
    FigTransverse("../figures/projection_transverse", 0.60);
    FigTransverse("../figures/projection_transverse_0p8", 0.80);
  }

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
