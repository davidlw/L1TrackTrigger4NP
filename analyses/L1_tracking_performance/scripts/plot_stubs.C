// ---------------------------------------------------------------------------
// plot_stubs.C -- draw the stub distributions filled by compare_stubs.C.
//
// Reads the histogram file, so axis ranges, binning and styling can be changed
// and redrawn in seconds without touching the samples again.
//
// The multiplicity is stored with one bin per stub. This macro rebins it onto a
// VARIABLE axis: one bin per stub up to 30, where the rejected spike lives and
// integers must stay resolved, then ~25 geometric bins per decade above. Counts
// are divided by bin width, so the y axis is events per stub -- without that the
// curve steps upward wherever the bin width jumps.
//
// N = 0 was clamped into the first bin when filling (it cannot sit on a log
// axis); the exact zero fraction is read from the stats histogram and printed in
// the legend, so nothing is hidden.
//
// Usage:
//   root -l -b -q 'plot_stubs.C("../output/stub_epos.root","_epos","EPOS pPb","epos")'
//   root -l -b -q 'plot_stubs.C("../output/stub_hydjet.root","_hydjet","HYDJET PbPb","hydjet",3e5)'
//
//   xmax  optional override of the multiplicity axis; 0 = 25% above the largest
//         count seen, which is what keeps the dummy tail inside the frame
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>

namespace {

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

// rebin the one-bin-per-stub histogram onto the variable axis, then divide by width
TH1D* Rebin(TH1D* fine, const char* name, const std::vector<double>& edges) {
  auto* h = new TH1D(name, "", (int)edges.size() - 1, edges.data());
  h->Sumw2();
  for (int b = 1; b <= fine->GetNbinsX(); ++b) {
    const double c = fine->GetBinContent(b);
    if (c > 0) h->Fill(fine->GetBinCenter(b), c);
  }
  h->Scale(1.0, "width");
  return h;
}

void Style(TH1D* h, int col, int ls) {
  h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(ls);
}

void Draw4(TH1D* da, TH1D* dr, TH1D* ua, TH1D* ur, const char* title,
           const char* outfile, const char* xtitle, const char* ytitle,
           bool logx, const char* legDA, const char* legDR,
           const char* legUA, const char* legUR,
           const char* legDR2 = "", const char* legUR2 = "") {
  auto* c = new TCanvas(Form("c%s", outfile), "", 950, 700);
  if (logx) c->SetLogx();
  c->SetLogy();
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetRightMargin(0.30);   // the legend sits outside the frame
  c->SetTopMargin(0.10);

  Style(da, kBlue + 2, 1); Style(dr, kCyan + 1, 2);
  Style(ua, kRed + 1, 1);  Style(ur, kGreen + 2, 2);

  // y range from the contents: these are densities, so with wide bins the values
  // run far below 1 and a fixed minimum would clip whole curves away
  double mx = 0, mn = 1e300;
  for (auto* h : {da, dr, ua, ur})
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
      const double v = h->GetBinContent(b);
      if (v > 0) { mx = std::max(mx, v); mn = std::min(mn, v); }
    }
  if (mx <= 0) { mx = 1; mn = 0.1; }

  da->SetTitle("");
  da->GetXaxis()->SetTitle(xtitle);
  da->GetYaxis()->SetTitle(ytitle);
  da->GetXaxis()->SetTitleSize(0.045);
  da->GetYaxis()->SetTitleSize(0.045);
  da->GetYaxis()->SetTitleOffset(1.25);
  const double xlo = da->GetXaxis()->GetXmin(), xhi = da->GetXaxis()->GetXmax();
  if (logx && std::log10(xhi / xlo) < 3.0) {   // extra labels collide on a long axis
    da->GetXaxis()->SetMoreLogLabels();
    da->GetXaxis()->SetNoExponent();
  }
  da->SetMaximum(mx * 30);
  da->SetMinimum(std::max(mn * 0.5, mx * 1e-9));
  da->Draw("HIST");
  dr->Draw("HIST SAME");
  ua->Draw("HIST SAME");
  ur->Draw("HIST SAME");

  auto* leg = new TLegend(0.71, 0.30, 0.995, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.027);
  leg->AddEntry(da, "Default accepted", "l");
  leg->AddEntry((TObject*)nullptr, legDA, "");
  leg->AddEntry((TObject*)nullptr, " ", "");
  leg->AddEntry(dr, "Default rejected", "l");
  leg->AddEntry((TObject*)nullptr, legDR, "");
  if (strlen(legDR2)) leg->AddEntry((TObject*)nullptr, legDR2, "");
  leg->AddEntry((TObject*)nullptr, " ", "");
  leg->AddEntry(ua, "Dummy accepted", "l");
  leg->AddEntry((TObject*)nullptr, legUA, "");
  leg->AddEntry((TObject*)nullptr, " ", "");
  leg->AddEntry(ur, "Dummy rejected", "l");
  leg->AddEntry((TObject*)nullptr, legUR, "");
  if (strlen(legUR2)) leg->AddEntry((TObject*)nullptr, legUR2, "");
  leg->Draw();

  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042); tx.SetTextAlign(22);
  tx.DrawLatex(0.42, 0.95, title);
  SaveBoth(c, outfile);
}

TH1D* Get(TFile* f, const char* n) {
  auto* h = (TH1D*)f->Get(n);
  if (!h) printf("[error] missing %s\n", n);
  return h;
}

}  // namespace

void plot_stubs(const char* fname = "../output/stub_epos.root", const char* tag = "_epos",
                const char* label = "EPOS pPb", const char* pfx = "epos", double xmax = 0) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) { printf("[error] cannot open %s\n", fname); return; }

  TH1D* st = Get(f, "stats");
  if (!st) return;
  const double nevD = st->GetBinContent(1), naD = st->GetBinContent(2),
               nrD = st->GetBinContent(3), zD = st->GetBinContent(4),
               nevU = st->GetBinContent(5), naU = st->GetBinContent(6),
               nrU = st->GetBinContent(7), zU = st->GetBinContent(8);

  TH1D* fDA = Get(f, "multfine_defA"); TH1D* fDR = Get(f, "multfine_defR");
  TH1D* fUA = Get(f, "multfine_dumA"); TH1D* fUR = Get(f, "multfine_dumR");
  if (!fDA || !fDR || !fUA || !fUR) return;

  // largest populated bin across all four, then 25% headroom
  int nmax = 1;
  for (auto* h : {fDA, fDR, fUA, fUR})
    for (int b = h->GetNbinsX(); b >= 1; --b)
      if (h->GetBinContent(b) > 0) { nmax = std::max(nmax, (int)h->GetBinLowEdge(b + 1)); break; }
  const double hi = xmax > 0 ? xmax : std::max(60.0, nmax * 1.25);
  printf("largest stub count = %d  ->  multiplicity axis to %.0f\n", nmax, hi);

  std::vector<double> me;
  for (int i = 0; i <= 30; ++i) me.push_back(i + 0.5);
  const double lo = 30.5;
  const int ng = std::max(45, (int)std::lround(25.0 * std::log10(hi / lo)));
  for (int i = 1; i <= ng; ++i) me.push_back(lo * std::pow(hi / lo, (double)i / ng));

  Draw4(Rebin(fDA, "mDA", me), Rebin(fDR, "mDR", me),
        Rebin(fUA, "mUA", me), Rebin(fUR, "mUR", me),
        Form("Stubs per event, %s", label),
        Form("../figures/%s_nstub%s.pdf", pfx, tag),
        "Stubs per event", "Events / stub", true,
        Form("#LTN#GT = %.1f/evt", nevD ? naD / nevD : 0),
        Form("#LTN#GT = %.1f/evt, %.1f%% rej.", nevD ? nrD / nevD : 0,
             (naD + nrD) ? 100 * nrD / (naD + nrD) : 0),
        Form("#LTN#GT = %.1f/evt", nevU ? naU / nevU : 0),
        Form("#LTN#GT = %.1f/evt, %.1f%% rej.", nevU ? nrU / nevU : 0,
             (naU + nrU) ? 100 * nrU / (naU + nrU) : 0),
        Form("0 rej. in %.1f%% of evts", nevD ? 100 * zD / nevD : 0),
        Form("0 rej. in %.1f%% of evts", nevU ? 100 * zU / nevU : 0));

  // the per-event stub-position distributions
  auto perEvt = [&](const char* base, const char* xt, const char* out, bool lx) {
    TH1D* a = (TH1D*)Get(f, Form("%s_defA", base))->Clone(Form("c%s1", base));
    TH1D* b = (TH1D*)Get(f, Form("%s_defR", base))->Clone(Form("c%s2", base));
    TH1D* cc = (TH1D*)Get(f, Form("%s_dumA", base))->Clone(Form("c%s3", base));
    TH1D* d = (TH1D*)Get(f, Form("%s_dumR", base))->Clone(Form("c%s4", base));
    if (nevD) { a->Scale(1.0 / nevD); b->Scale(1.0 / nevD); }
    if (nevU) { cc->Scale(1.0 / nevU); d->Scale(1.0 / nevU); }
    Draw4(a, b, cc, d, label, out, xt, "stubs / event / bin", lx,
          "accepted", "rejected", "accepted", "rejected");
  };
  perEvt("r", "stub r [cm]", Form("../figures/%s_stub_r%s.pdf", pfx, tag), false);
  perEvt("z", "stub z [cm]", Form("../figures/%s_stub_z%s.pdf", pfx, tag), false);
  perEvt("ld", "0-5 = barrel L1-L6, 6-10 = disk D1-D5",
         Form("../figures/%s_stub_layerdisk%s.pdf", pfx, tag), false);
  perEvt("bend", "stub trigBend (dummy is the 999999 sentinel, off scale)",
         Form("../figures/%s_stub_bend%s.pdf", pfx, tag), false);

  printf("\nfigures written to ../figures/\n");
}
