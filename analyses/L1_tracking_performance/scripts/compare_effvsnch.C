// ---------------------------------------------------------------------------
// compare_effvsnch.C
//
// Event-level trigger efficiency vs charged multiplicity.
//
//   x : Nch = number of truth charged particles with pT > 0.4 GeV, |eta| < 2.4
//   y : fraction of events at that Nch with >= 1 RECONSTRUCTED track with
//       |eta| < 2.4 and pT > threshold, for thresholds 0.4 / 1.0 / 2.0 GeV
//
// Denominator is every event in the Nch bin, so events the trigger cannot fire on
// stay in the denominator -- the same convention as compare_mbeff.C.
//
// TWO LIMITATIONS OF THE CURRENT SAMPLES, both lifted by regenerating:
//
//  1. ptMinTP = 1.0 GeV in the mixing module means the truth collection has
//     essentially nothing below 1 GeV (~1% of TPs). So "Nch with pT > 0.4" is in
//     practice "Nch with pT > 1.0" -- the 0.4 threshold on the x axis does not
//     yet do anything. Regenerate with a lower ptMinTP for the intended quantity.
//
//  2. The track finder in these productions emits nothing below 1.95 GeV, so the
//     0.4 and 1.0 GeV trigger thresholds are DEGENERATE -- both reduce to ">= 1
//     reconstructed track at all" and their curves lie exactly on top of each
//     other. Only the 2.0 GeV curve differs. On a sample with the lowered pT
//     threshold the three curves separate, which is the point of the plot.
//
// Nch is taken from each production's OWN truth collection, so the two stub
// configurations use slightly different x-axis definitions (their TP collections
// differ through the >=3 stub / >=3 layer preselection). Noted rather than
// corrected, because it also disappears once the preselection is relaxed.
//
// Usage:
//   root -l -b -q 'compare_effvsnch.C("<def>","<dum>","HYDJet PbPb","hydjet",300,30,100)'
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";
const double kEtaMax = 2.4;
const double kNchPtMin = 0.4;                       // pT cut defining Nch
const int kNThr = 3;
const double kThr[kNThr] = {0.4, 1.0, 2.0};         // trigger thresholds

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

struct Set {
  TH1D* den;                 // events per Nch bin
  TH1D* num[kNThr];          // events firing, per threshold
  long nEvt = 0;
  double sumNch = 0;
  void Book(const char* tag, int nbins, double nchMax) {
    den = new TH1D(Form("den_nch_%s", tag), "", nbins, 0, nchMax);
    den->Sumw2();
    for (int j = 0; j < kNThr; ++j) {
      num[j] = new TH1D(Form("num_nch_%s_t%d", tag, j), "", nbins, 0, nchMax);
      num[j]->Sumw2();
    }
  }
};

void Fill(const char* dir, int nfiles, Set& s) {
  for (const auto& p : ListRootFiles(dir, nfiles)) {
    TFile* f = TFile::Open(p);
    if (!f || f->IsZombie()) continue;
    TTree* t = (TTree*)f->Get(kTreePath);
    if (!t || !t->GetBranch("tp_pt") || !t->GetBranch("trk_pt")) {
      if (f) f->Close();
      continue;
    }
    std::vector<float> *tpPt = nullptr, *tpEta = nullptr, *trkPt = nullptr,
                       *trkEta = nullptr;
    t->SetBranchStatus("*", 0);
    for (auto b : {"tp_pt", "tp_eta", "trk_pt", "trk_eta"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("tp_pt", &tpPt);
    t->SetBranchAddress("tp_eta", &tpEta);
    t->SetBranchAddress("trk_pt", &trkPt);
    t->SetBranchAddress("trk_eta", &trkEta);

    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      // Nch from truth
      int nch = 0;
      size_t nt = std::min(tpPt->size(), tpEta->size());
      for (size_t k = 0; k < nt; ++k)
        if (std::abs(tpEta->at(k)) < kEtaMax && tpPt->at(k) > kNchPtMin) ++nch;
      // leading reco track in acceptance
      double lead = -1;
      size_t nr = std::min(trkPt->size(), trkEta->size());
      for (size_t k = 0; k < nr; ++k)
        if (std::abs(trkEta->at(k)) < kEtaMax && trkPt->at(k) > lead) lead = trkPt->at(k);

      ++s.nEvt;
      s.sumNch += nch;
      s.den->Fill(nch);
      for (int j = 0; j < kNThr; ++j)
        if (lead > kThr[j]) s.num[j]->Fill(nch);
    }
    f->Close();
  }
  printf("[info] %s: %ld events, <Nch(pT>%.1f,|eta|<%.1f)> = %.1f\n", dir, s.nEvt,
         kNchPtMin, kEtaMax, s.nEvt ? s.sumNch / s.nEvt : 0);
  fflush(stdout);
}

TGraphAsymmErrors* Eff(Set& s, int j, int color, int style, int marker) {
  auto* g = new TGraphAsymmErrors(s.num[j], s.den, "cl=0.683 b(1,1) mode");
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(2);
  return g;
}

}  // namespace

void compare_effvsnch(const char* dirDef, const char* dirDum, const char* sample,
                      const char* tag, double nchMax = 300, int nbins = 30,
                      int nfiles = 100) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Set D, U;
  D.Book(Form("def_%s", tag), nbins, nchMax);
  U.Book(Form("dum_%s", tag), nbins, nchMax);
  Fill(dirDef, nfiles, D);
  Fill(dirDum, nfiles, U);

  printf("\n===== %s: trigger efficiency vs Nch =====\n", sample);
  printf("Nch: truth charged, pT > %.1f GeV, |eta| < %.1f\n", kNchPtMin, kEtaMax);
  printf("fire: >= 1 reconstructed track, |eta| < %.1f, pT > threshold\n\n", kEtaMax);
  printf("%12s %10s", "Nch bin", "N(evt)");
  for (int j = 0; j < kNThr; ++j) printf("%12s", Form("def>%.1f", kThr[j]));
  for (int j = 0; j < kNThr; ++j) printf("%12s", Form("dum>%.1f", kThr[j]));
  printf("\n");
  for (int b = 1; b <= D.den->GetNbinsX(); ++b) {
    double nd = D.den->GetBinContent(b), nu = U.den->GetBinContent(b);
    if (nd < 1 && nu < 1) continue;
    printf("%5.0f-%6.0f %10.0f", D.den->GetBinLowEdge(b), D.den->GetBinLowEdge(b + 1), nd);
    for (int j = 0; j < kNThr; ++j)
      printf("%12s", nd > 0 ? Form("%.4f", D.num[j]->GetBinContent(b) / nd) : "-");
    for (int j = 0; j < kNThr; ++j)
      printf("%12s", nu > 0 ? Form("%.4f", U.num[j]->GetBinContent(b) / nu) : "-");
    printf("\n");
  }

  const int col[kNThr] = {kGreen + 3, kBlue + 1, kRed + 1};
  auto* c = new TCanvas(Form("c_effnch_%s", tag), "", 950, 680);
  c->SetGridy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
  auto* fr = c->DrawFrame(0, 0, nchMax, 1.38);
  fr->GetXaxis()->SetTitle(Form("N_{ch}  (truth, p_{T} > %.1f GeV, |#eta| < %.1f)", kNchPtMin, kEtaMax));
  fr->GetYaxis()->SetTitle("trigger efficiency");
  fr->GetXaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleOffset(1.3);
  auto* one = new TLine(0, 1, nchMax, 1);
  one->SetLineStyle(3);
  one->SetLineColor(kGray + 1);
  one->Draw();

  auto* leg = new TLegend(0.40, 0.66, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.028);
  leg->SetNColumns(2);
  for (int j = 0; j < kNThr; ++j) {
    auto* gd = Eff(D, j, col[j], 1, 20);
    auto* gu = Eff(U, j, col[j], 2, 24);
    gd->Draw("PL same");
    gu->Draw("PL same");
    leg->AddEntry(gd, Form("Default, p_{T} > %.1f", kThr[j]), "lp");
    leg->AddEntry(gu, Form("Dummy, p_{T} > %.1f", kThr[j]), "lp");
  }
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.033);
  tx.DrawLatex(0.13, 0.94, Form("%s -- #geq 1 reconstructed track vs N_{ch}", sample));

  SaveBoth(c, Form("../figures/eff_vs_nch_%s.pdf", tag));

  TFile* fout = TFile::Open(Form("effvsnch_%s.root", tag), "RECREATE");
  D.den->Write();
  U.den->Write();
  for (int j = 0; j < kNThr; ++j) { D.num[j]->Write(); U.num[j]->Write(); }
  fout->Close();
  printf("\nwrote ../figures/eff_vs_nch_%s.pdf\n", tag);
}
