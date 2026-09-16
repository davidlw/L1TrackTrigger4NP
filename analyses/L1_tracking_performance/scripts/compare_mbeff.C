// ---------------------------------------------------------------------------
// compare_mbeff.C
//
// Minimum-bias trigger efficiency vs the pT threshold, for HYDJet PbPb.
//
//   MB efficiency(pTmin) = fraction of events containing at least one charged
//                          track with |eta| < 2.4 and pT > pTmin
//
// Four curves:
//   truth (default TPs) / truth (dummy TPs)  -- what fraction of events HAVE
//       such a particle, i.e. the ceiling a perfect trigger could reach
//   reco (default) / reco (dummy)            -- what the L1 track finding
//       actually delivers, using trk_pt / trk_eta
//
// IMPORTANT CAVEAT on "truth": the ntuple stores a tracking particle only if it
// has >=3 stubs in >=3 layers (TP_minNStub / TP_minNStubLayer). That is a
// stub-based cut, so the stored TP collection is NOT the full generator record,
// and it differs between the two productions -- the dummy sample keeps 1,084,113
// TPs where the default keeps 83,452. Both truth curves are therefore shown: they
// bracket the true generator-level answer rather than pinning it. For a genuine
// generator-level ceiling the ntuples would need regenerating with the stub
// preselection removed.
//
// All TPs stored are charged (the producer applies charge != 0; verified
// tp_charge != 0 for 100% of entries in both samples).
//
// Usage: root -l -b -q 'compare_mbeff.C(100)'
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

// Set on the command line (arguments 2 and 3); a directory of *.root or one file.
const char* kDirDefault = "";
const char* kDirDummy = "";
const char* kSample = "HYDJet PbPb";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

const double kEtaMax = 2.4;
const int kNBins = 240;      // 0 to 6 GeV in 25 MeV steps
const double kPtMaxScan = 6.0;

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

// Requiring >= N tracks above threshold is equivalent to asking whether the Nth
// HIGHEST pT in the event exceeds the threshold. So storing the Nth-highest per
// event and taking the reverse-cumulative gives the whole efficiency curve for
// that N in one pass.
//
// N = 1 saturates immediately in PbPb: no HYDJet event has a leading charged
// particle below 1 GeV, so that curve is flat below ~1.2 GeV by construction and
// says nothing about low-pT reach. The higher N are where the reach shows up.
const int kNMult = 5;
const int kMultReq[kNMult] = {1, 2, 5, 10, 20};

struct Scan {
  TH1D* lead;               // N = 1, kept for the truth-vs-reco comparison
  TH1D* nth[kNMult];        // Nth-highest pT, for the >= N curves
  long nEvt = 0, nEmpty = 0;
  void Book(const char* name) {
    lead = new TH1D(name, "", kNBins, 0, kPtMaxScan);
    lead->Sumw2();
    for (int j = 0; j < kNMult; ++j) {
      nth[j] = new TH1D(Form("%s_n%d", name, kMultReq[j]), "", kNBins, 0, kPtMaxScan);
      nth[j]->Sumw2();
    }
  }
};

// which = 0 -> tracking particles (tp_*), 1 -> reconstructed tracks (trk_*)
void Fill(const char* dir, int nfiles, int which, Scan& s) {
  for (const auto& path : ListRootFiles(dir, nfiles)) {
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) continue;
    TTree* t = (TTree*)f->Get(kTreePath);
    if (!t) { f->Close(); continue; }

    const char* bPt = which == 0 ? "tp_pt" : "trk_pt";
    const char* bEta = which == 0 ? "tp_eta" : "trk_eta";
    if (!t->GetBranch(bPt) || !t->GetBranch(bEta)) { f->Close(); continue; }

    std::vector<float> *pt = nullptr, *eta = nullptr;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus(bPt, 1);
    t->SetBranchStatus(bEta, 1);
    t->SetBranchAddress(bPt, &pt);
    t->SetBranchAddress(bEta, &eta);

    std::vector<float> inAcc;
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      size_t n = std::min(pt->size(), eta->size());
      inAcc.clear();
      for (size_t k = 0; k < n; ++k)
        if (std::abs(eta->at(k)) < kEtaMax) inAcc.push_back(pt->at(k));
      ++s.nEvt;
      if (inAcc.empty()) { ++s.nEmpty; continue; }
      std::sort(inAcc.begin(), inAcc.end(), std::greater<float>());
      auto clip = [](double v) { return std::min(v, kPtMaxScan - 1e-6); };
      s.lead->Fill(clip(inAcc[0]));
      // Nth-highest: an event with fewer than N tracks in acceptance can never
      // satisfy ">= N above threshold", so it simply does not fill that histogram.
      for (int j = 0; j < kNMult; ++j) {
        int need = kMultReq[j];
        if ((int)inAcc.size() >= need) s.nth[j]->Fill(clip(inAcc[need - 1]));
      }
    }
    f->Close();
  }
  printf("[info] %s (%s): %ld events, %ld with nothing in |eta|<2.4\n", dir,
         which == 0 ? "truth" : "reco", s.nEvt, s.nEmpty);
  fflush(stdout);
}

// Reverse-cumulative of the leading-pT histogram, with binomial errors.
TGraphErrors* EffCurve(const Scan& s, int color, int style) {
  auto* g = new TGraphErrors();
  int p = 0;
  double N = (double)s.nEvt;
  if (N <= 0) return g;
  for (int i = 1; i <= s.lead->GetNbinsX(); ++i) {
    double thr = s.lead->GetBinLowEdge(i);
    double pass = s.lead->Integral(i, s.lead->GetNbinsX() + 1);
    double eff = pass / N;
    double err = std::sqrt(std::max(0.0, eff * (1 - eff) / N));
    g->SetPoint(p, thr, eff);
    g->SetPointError(p, 0, err);
    ++p;
  }
  g->SetLineColor(color);
  g->SetLineStyle(style);
  g->SetLineWidth(3);
  g->SetMarkerColor(color);
  return g;
}

double EffAt(const Scan& s, double thr) {
  if (s.nEvt <= 0) return 0;
  int b = s.lead->FindBin(thr + 1e-9);
  return s.lead->Integral(b, s.lead->GetNbinsX() + 1) / (double)s.nEvt;
}

}  // namespace

void compare_mbeff(int nfiles = 100, const char* dirDef = "", const char* dirDum = "",
                   const char* sample = "", const char* tag = "hydjet") {
  if (strlen(dirDef)) kDirDefault = dirDef;
  if (strlen(dirDum)) kDirDummy = dirDum;
  if (!strlen(kDirDefault) || !strlen(kDirDummy)) {
    printf("give the Default-stub and Dummy-stub ntuple locations (a directory of *.root, or one file)\n");
    return;
  }
  if (strlen(sample)) kSample = sample;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Scan tD, tU, rD, rU;
  tD.Book("lead_truth_def");
  tU.Book("lead_truth_dum");
  rD.Book("lead_reco_def");
  rU.Book("lead_reco_dum");

  Fill(kDirDefault, nfiles, 0, tD);
  Fill(kDirDummy, nfiles, 0, tU);
  Fill(kDirDefault, nfiles, 1, rD);
  Fill(kDirDummy, nfiles, 1, rU);

  printf("\n===== %s minimum-bias trigger efficiency =====\n", kSample);
  printf("at least one charged track with |eta| < %.1f and pT > threshold\n", kEtaMax);
  printf("truth = stored tracking particles (see caveat in header); reco = L1 tracks\n\n");
  printf("%10s %14s %14s %14s %14s\n", "pTmin", "truth(def)", "truth(dum)", "reco(def)",
         "reco(dum)");
  for (double thr : {0.2, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0})
    printf("%10.2f %14.4f %14.4f %14.4f %14.4f\n", thr, EffAt(tD, thr), EffAt(tU, thr),
           EffAt(rD, thr), EffAt(rU, thr));

  auto* gTD = EffCurve(tD, kGray + 2, 2);
  auto* gTU = EffCurve(tU, kBlack, 2);
  auto* gRD = EffCurve(rD, kBlue + 1, 1);
  auto* gRU = EffCurve(rU, kRed + 1, 1);

  auto* c = new TCanvas("c_mbeff", "", 900, 660);
  c->SetGridy();
  c->SetGridx();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
  auto* fr = c->DrawFrame(0, 0, 5.0, 1.35);
  fr->GetXaxis()->SetTitle("p_{T} threshold [GeV]");
  fr->GetYaxis()->SetTitle("minimum-bias trigger efficiency");
  fr->GetXaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleSize(0.045);
  fr->GetYaxis()->SetTitleOffset(1.35);

  auto* one = new TLine(0, 1, 5.0, 1);
  one->SetLineStyle(3);
  one->SetLineColor(kGray + 1);
  one->Draw();

  gTU->Draw("L same");
  gTD->Draw("L same");
  gRD->Draw("L same");
  gRU->Draw("L same");

  auto* leg = new TLegend(0.40, 0.68, 0.90, 0.89);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->AddEntry(gTU, "truth ceiling  (dummy TP collection)", "l");
  leg->AddEntry(gTD, "truth ceiling  (default TP collection)", "l");
  leg->AddEntry(gRD, "reconstructed  (Default Stub)", "l");
  leg->AddEntry(gRU, "reconstructed  (Dummy Stub)", "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.13, 0.94,
               Form("%s -- MB efficiency: #geq 1 charged track, |#eta| < 2.4", kSample));

  SaveBoth(c, Form("../figures/mbeff_vs_ptmin_%s.pdf", tag));

  // ---- >= N tracks, truth only -------------------------------------------
  // Uses the dummy TP collection: the default's TPs are a 99.7% subset of it, so
  // it is effectively the union of the two and the closest available stand-in for
  // generator truth.
  auto* c2 = new TCanvas("c_mbeff_n", "", 900, 660);
  c2->SetGridy();
  c2->SetGridx();
  c2->SetLeftMargin(0.13);
  c2->SetBottomMargin(0.13);
  auto* fr2 = c2->DrawFrame(0, 0, 5.0, 1.35);
  fr2->GetXaxis()->SetTitle("p_{T} threshold [GeV]");
  fr2->GetYaxis()->SetTitle("fraction of events with #geq N such tracks");
  fr2->GetXaxis()->SetTitleSize(0.045);
  fr2->GetYaxis()->SetTitleSize(0.045);
  fr2->GetYaxis()->SetTitleOffset(1.35);
  auto* one2 = new TLine(0, 1, 5.0, 1);
  one2->SetLineStyle(3);
  one2->SetLineColor(kGray + 1);
  one2->Draw();

  const int colN[kNMult] = {kBlack, kBlue + 1, kGreen + 3, kOrange + 7, kRed + 1};
  auto* legN = new TLegend(0.62, 0.62, 0.90, 0.89);
  legN->SetBorderSize(0);
  legN->SetFillStyle(0);
  legN->SetTextSize(0.032);
  for (int j = 0; j < kNMult; ++j) {
    auto* g = new TGraphErrors();
    int p = 0;
    double N = (double)tU.nEvt;
    for (int i = 1; i <= tU.nth[j]->GetNbinsX(); ++i) {
      double pass = tU.nth[j]->Integral(i, tU.nth[j]->GetNbinsX() + 1);
      double eff = pass / N;
      g->SetPoint(p, tU.nth[j]->GetBinLowEdge(i), eff);
      g->SetPointError(p, 0, std::sqrt(std::max(0.0, eff * (1 - eff) / N)));
      ++p;
    }
    g->SetLineColor(colN[j]);
    g->SetLineWidth(3);
    g->Draw("L same");
    legN->AddEntry(g, Form("#geq %d tracks", kMultReq[j]), "l");
  }
  legN->Draw();
  TLatex tx2;
  tx2.SetNDC();
  tx2.SetTextSize(0.034);
  tx2.DrawLatex(0.13, 0.94,
                Form("%s -- truth: #geq N charged tracks, |#eta| < 2.4", kSample));
  SaveBoth(c2, Form("../figures/mbeff_vs_ptmin_%s_multN.pdf", tag));

  printf("\n--- truth, >= N charged tracks (dummy TP collection) ---\n");
  printf("%10s", "pTmin");
  for (int j = 0; j < kNMult; ++j) printf("%12s", Form(">=%d", kMultReq[j]));
  printf("\n");
  for (double thr : {0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0}) {
    printf("%10.2f", thr);
    for (int j = 0; j < kNMult; ++j) {
      int b = tU.nth[j]->FindBin(thr + 1e-9);
      printf("%12.4f", tU.nth[j]->Integral(b, tU.nth[j]->GetNbinsX() + 1) / (double)tU.nEvt);
    }
    printf("\n");
  }

  TFile* fout = TFile::Open(Form("mbeff_%s.root", tag), "RECREATE");
  for (auto* h : {tD.lead, tU.lead, rD.lead, rU.lead}) h->Write();
  for (int j = 0; j < kNMult; ++j) { tD.nth[j]->Write(); tU.nth[j]->Write(); }
  gTD->Write("eff_truth_def");
  gTU->Write("eff_truth_dum");
  gRD->Write("eff_reco_def");
  gRU->Write("eff_reco_dum");
  fout->Close();
  printf("\nwrote ../figures/mbeff_vs_ptmin_%s.pdf and mbeff_%s.root\n", tag, tag);
}
