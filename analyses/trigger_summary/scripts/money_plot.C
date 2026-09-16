// ---------------------------------------------------------------------------
// money_plot.C
//
// One-glance summary: event trigger efficiency per physics process, for
//
//   L1 default        the stock tracklet algorithm (pT floor ~2 GeV)
//   L1 low-pT         the new low-pT algorithm (dummy stubs + widened windows)
//   offline           highPurity offline tracks -- the ceiling
//
// with a lower panel showing either
//   ratioMode = "gain"     low-pT L1 / default L1  -- the signal gain of the new
//                          trigger, on a log scale (the headline number)
//   ratioMode = "offline"  L1 / offline for both algorithms -- how close each
//                          is to the ceiling
//
// The efficiency for a process is
//
//   eff = N(events with >= k tracks, pT > pT_algo, |eta| < 2.4)
//         ------------------------------------------------------
//         N(events with >= k PRIMARY charged particles, pT > pT_fid, |eta| < 2.4)
//
// The DENOMINATOR is the same for all three series: the fiducial truth at the
// low-pT threshold (pT_fid = 0.4 GeV), NOT all generated events. Otherwise a
// soft process (rho -> pi pi) looks inefficient because of acceptance rather
// than tracking, and processes are not comparable.
//
// The NUMERATOR threshold is each trigger's own: the default L1 tracking only
// exists above ~2 GeV, so its trigger is ">= k tracks with pT > 2 GeV"; the
// low-pT L1 and offline trigger at pT > 0.4. This compares the triggers as
// they would be deployed, so the default column shows the acceptance loss of
// the high threshold AND its tracking inefficiency together -- which is the
// point of the plot.
//
// Input is a text table (money_plot_input.txt) so that the counts can be
// filled by whichever macro produced them without touching this file. Errors
// are Clopper-Pearson via TEfficiency.
//
// Usage:
//   root -l -b -q 'money_plot.C("money_plot_input.txt", "../figures/money_plot", "gain", 1, 0.4, 2.0, 0.4, 0.4)'
//
//   arg 3 = ratioMode; arg 4 = k (tracks required); arg 5 = pT_fid
//   (denominator); args 6-8 = numerator thresholds for default L1, low-pT L1,
//   offline [GeV]. The thresholds are only used for the labels on the canvas --
//   the counts in the input file must have been made with the same selection.
// ---------------------------------------------------------------------------

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Process {
  TString label;
  long nfid = 0, ndef = 0, nlow = 0, noff = 0;
};

// A line "--- title" in the table starts a new group: a heavier divider is
// drawn before the next process and the title is written above the group.
struct Group {
  int first;  // index of the first process in the group
  TString title;
};

bool ReadTable(const char* path, std::vector<Process>& procs, std::vector<Group>& groups, bool& isTemplate) {
  std::ifstream in(path);
  if (!in) {
    printf("cannot open %s\n", path);
    return false;
  }
  std::string line;
  while (std::getline(in, line)) {
    if (line.find("TEMPLATE") != std::string::npos)
      isTemplate = true;
    if (line.rfind("---", 0) == 0) {
      TString t = line.substr(3);
      t = t.Strip(TString::kBoth);
      t.ReplaceAll("_", " ");
      groups.push_back({(int)procs.size(), t});
      continue;
    }
    // a leading '#' is a comment UNLESS it starts a TLatex label like #gamma
    if (line.empty() || (line[0] == '#' && line.size() > 1 && line[1] == ' '))
      continue;
    if (line[0] == '#' && !isalpha(line[1]))
      continue;
    std::istringstream ss(line);
    Process p;
    std::string lab;
    if (!(ss >> lab >> p.nfid >> p.ndef >> p.nlow >> p.noff))
      continue;
    p.label = lab;
    p.label.ReplaceAll("_", " ");
    procs.push_back(p);
  }
  return !procs.empty();
}

// Category axis: process i occupies [i, i+1); the three series are offset
// inside the bin so the markers do not sit on top of each other.
TGraphAsymmErrors* MakeSeries(const std::vector<Process>& procs, int which, double xoff) {
  TH1D num("num", "", procs.size(), 0, procs.size());
  TH1D den("den", "", procs.size(), 0, procs.size());
  for (size_t i = 0; i < procs.size(); ++i) {
    long n = which == 0 ? procs[i].ndef : which == 1 ? procs[i].nlow : procs[i].noff;
    num.SetBinContent(i + 1, n);
    den.SetBinContent(i + 1, procs[i].nfid);
  }
  TEfficiency eff(num, den);
  eff.SetStatisticOption(TEfficiency::kFCP);
  auto* g = new TGraphAsymmErrors(procs.size());
  for (size_t i = 0; i < procs.size(); ++i) {
    g->SetPoint(i, i + 0.5 + xoff, eff.GetEfficiency(i + 1));
    g->SetPointError(i, 0, 0, eff.GetEfficiencyErrorLow(i + 1), eff.GetEfficiencyErrorUp(i + 1));
  }
  return g;
}

// L1 / offline, per process. Errors propagated as uncorrelated ratio of two
// binomials, which is conservative (the numerators share events).
TGraphAsymmErrors* MakeRatio(TGraphAsymmErrors* l1, TGraphAsymmErrors* off) {
  auto* g = new TGraphAsymmErrors(l1->GetN());
  for (int i = 0; i < l1->GetN(); ++i) {
    double x, a, b;
    l1->GetPoint(i, x, a);
    off->GetPoint(i, x, b);
    b = off->GetY()[i];
    double r = b > 0 ? a / b : 0;
    auto rel = [&](double ea, double eb) {
      return (a > 0 && b > 0) ? r * std::sqrt(ea * ea / (a * a) + eb * eb / (b * b)) : 0;
    };
    g->SetPoint(i, x, r);
    g->SetPointError(i, 0, 0, rel(l1->GetEYlow()[i], off->GetEYhigh()[i]),
                     rel(l1->GetEYhigh()[i], off->GetEYlow()[i]));
  }
  return g;
}

void SetStyle(TGraphAsymmErrors* g, int color, int marker) {
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.6);
  g->SetLineWidth(2);
}

void LabelAxis(TH1* frame, const std::vector<Process>& procs, double size) {
  for (size_t i = 0; i < procs.size(); ++i)
    frame->GetXaxis()->SetBinLabel(i + 1, procs[i].label);
  frame->GetXaxis()->SetLabelSize(size);
  frame->GetXaxis()->LabelsOption("h");
}

}  // namespace

void money_plot(const char* input = "money_plot_input.txt",
                const char* outBase = "../figures/money_plot",
                const char* ratioMode = "gain",
                int kTracks = 1,
                double ptFid = 0.4,
                double ptDef = 2.0,
                double ptLow = 0.4,
                double ptOff = 0.4) {
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(4);

  std::vector<Process> procs;
  std::vector<Group> groups;
  bool isTemplate = false;
  if (!ReadTable(input, procs, groups, isTemplate))
    return;
  const int n = procs.size();

  const int cDef = kGray + 2, cLow = kRed + 1, cOff = kAzure + 2;
  auto* gDef = MakeSeries(procs, 0, -0.2);
  auto* gLow = MakeSeries(procs, 1, 0.0);
  auto* gOff = MakeSeries(procs, 2, +0.2);
  SetStyle(gDef, cDef, 25);  // open square: the thing being replaced
  SetStyle(gLow, cLow, 21);  // filled square: the new result
  SetStyle(gOff, cOff, 34);  // star: the reference
  const bool gain = TString(ratioMode) == "gain";
  auto* rDef = gain ? nullptr : MakeRatio(gDef, gOff);
  auto* rLow = gain ? MakeRatio(gLow, gDef) : MakeRatio(gLow, gOff);
  if (rDef)
    SetStyle(rDef, cDef, 25);
  SetStyle(rLow, cLow, 21);

  auto* c = new TCanvas("c", "", std::max(900, 110 * n), 750);
  auto* pTop = new TPad("pTop", "", 0, 0.33, 1, 1);
  auto* pBot = new TPad("pBot", "", 0, 0, 1, 0.33);
  pTop->SetBottomMargin(0.02);
  pTop->SetLeftMargin(0.12);
  pTop->SetRightMargin(0.04);
  pBot->SetTopMargin(0.03);
  pBot->SetBottomMargin(0.36);
  pBot->SetLeftMargin(0.12);
  pBot->SetRightMargin(0.04);
  pTop->Draw();
  pBot->Draw();

  // ---- top: efficiency
  pTop->cd();
  auto* fTop = new TH1D("fTop", "", n, 0, n);
  fTop->SetMinimum(0);
  fTop->SetMaximum(1.30);  // room for group titles and the legend above the points
  fTop->GetYaxis()->SetTitle("Event trigger efficiency");
  fTop->GetYaxis()->SetTitleSize(0.045);
  fTop->GetYaxis()->SetTitleOffset(1.2);
  fTop->GetYaxis()->SetLabelSize(0.045);
  fTop->GetXaxis()->SetLabelSize(0);
  fTop->GetYaxis()->SetNdivisions(506);
  fTop->Draw();
  auto isBoundary = [&](int i) {
    for (auto& g : groups)
      if (g.first == i)
        return true;
    return false;
  };
  for (int i = 1; i < n; ++i) {
    auto* l = new TLine(i, 0, i, isBoundary(i) ? 1.25 : 1.0);
    l->SetLineStyle(isBoundary(i) ? 1 : 3);
    l->SetLineColor(isBoundary(i) ? kGray + 3 : kGray + 1);
    l->Draw();
  }
  TLatex gt;
  gt.SetTextSize(0.038);
  gt.SetTextAlign(21);
  gt.SetTextColor(kGray + 3);
  for (size_t k = 0; k < groups.size(); ++k) {
    int last = k + 1 < groups.size() ? groups[k + 1].first : n;
    if (last > groups[k].first)
      gt.DrawLatex(0.5 * (groups[k].first + last), 1.10, groups[k].title);
  }
  auto* one = new TLine(0, 1, n, 1);
  one->SetLineStyle(2);
  one->Draw();
  gDef->Draw("PZ same");
  gOff->Draw("PZ same");
  gLow->Draw("PZ same");

  auto* leg = new TLegend(0.14, 0.85, 0.95, 0.96);
  leg->SetNColumns(3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.038);
  leg->AddEntry(gDef, Form("L1 default, p_{T} > %g GeV", ptDef), "p");
  leg->AddEntry(gLow, Form("L1 low-p_{T}, p_{T} > %g GeV", ptLow), "p");
  leg->AddEntry(gOff, Form("Offline, p_{T} > %g GeV", ptOff), "p");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.05);
  tx.SetTextFont(62);
  tx.DrawLatex(0.12, 0.965, "CMS Phase-2");
  tx.SetTextFont(52);
  tx.SetTextSize(0.042);
  tx.DrawLatex(0.30, 0.965, "Simulation");
  tx.SetTextFont(42);
  tx.SetTextAlign(31);
  tx.DrawLatex(0.96, 0.965,
               Form("#geq%d track%s, |#eta| < 2.4;  denominator: #geq%d primary charged, p_{T} > %g GeV", kTracks,
                    kTracks > 1 ? "s" : "", kTracks, ptFid));

  if (isTemplate) {
    TLatex wm;
    wm.SetNDC();
    wm.SetTextColor(kRed - 4);
    wm.SetTextSize(0.09);
    wm.SetTextAngle(20);
    wm.SetTextAlign(22);
    wm.DrawLatex(0.5, 0.5, "TEMPLATE -- placeholder numbers");
  }

  // ---- bottom: L1 / offline
  pBot->cd();
  auto* fBot = new TH1D("fBot", "", n, 0, n);
  const double botMin = gain ? 0.5 : 0, botMax = gain ? 60 : 1.19;
  fBot->SetMinimum(botMin);
  fBot->SetMaximum(botMax);
  fBot->GetYaxis()->SetTitle(gain ? "low-p_{T} / default" : "L1 / offline");
  if (gain)
    pBot->SetLogy();
  fBot->GetYaxis()->SetTitleSize(0.09);
  fBot->GetYaxis()->SetTitleOffset(0.6);
  fBot->GetYaxis()->SetLabelSize(0.085);
  fBot->GetYaxis()->SetNdivisions(504);
  LabelAxis(fBot, procs, n > 8 ? 0.10 : 0.13);
  fBot->Draw();
  for (int i = 1; i < n; ++i) {
    auto* l = new TLine(i, botMin, i, botMax);
    l->SetLineStyle(isBoundary(i) ? 1 : 3);
    l->SetLineColor(isBoundary(i) ? kGray + 3 : kGray + 1);
    l->Draw();
  }
  auto* oneB = new TLine(0, 1, n, 1);
  oneB->SetLineStyle(2);
  oneB->Draw();
  if (rDef)
    rDef->Draw("PZ same");
  rLow->Draw("PZ same");
  if (gain) {
    // print the gain factor next to each point: it is the number people quote
    TLatex gl;
    gl.SetTextSize(0.075);
    gl.SetTextAlign(21);
    gl.SetTextColor(cLow);
    for (int i = 0; i < n; ++i) {
      double r = rLow->GetY()[i];
      if (r > 0)
        gl.DrawLatex(rLow->GetX()[i], r * 1.6, r < 10 ? Form("#times%.1f", r) : Form("#times%.0f", r));
    }
  }

  c->SaveAs(Form("%s.pdf", outBase));
  c->SaveAs(Form("%s.png", outBase));

  printf("%-40s %8s %8s %8s   %7s %6s\n", "process", "default", "low-pT", "offline", "low/def", "l/off");
  for (int i = 0; i < n; ++i) {
    double d = gDef->GetY()[i], l = gLow->GetY()[i], o = gOff->GetY()[i];
    printf("%-40s %8.3f %8.3f %8.3f   %7.2f %6.3f\n", procs[i].label.Data(), d, l, o, d > 0 ? l / d : 0,
           o > 0 ? l / o : 0);
  }
}
