// ---------------------------------------------------------------------------
// compare_stubmult.C
//
// Stub multiplicity per event, split by the allstub_isRejected flag (0 = accepted,
// 1 = rejected), overlaid for the Default Stub and Dummy Stub productions.
//
// All four distributions go on one canvas:
//   Default accepted / Default rejected / Dummy accepted / Dummy rejected
//
// Normalised to unit area, so each curve reads as "fraction of events with N
// stubs of this class". Log y, because the rejected classes are rare and would
// otherwise be invisible next to the accepted ones.
//
// Defaults reproduce the STARlight QED mumu plot. Pass different directories to
// run the same comparison on another sample; pythia pp needs log-x binning
// because its four classes span ~70 to ~75000 stubs per event.
//
// Usage:
//   root -l -b -q 'compare_stubmult.C(10)'                        // QED mumu
//   root -l -b -q 'compare_stubmult.C(2,"","",1,1000,2e5,"pp")'   // pythia pp
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

// Set on the command line (arguments 2 and 3); a directory of *.root or one file.
const char* kQEDDefault = "";
const char* kQEDDummy = "";
const char* kPPDefault = "";
const char* kPPDummy = "";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

// Log-spaced bins, with the first bin holding zero (log(0) is undefined but
// "no rejected stubs this event" is a real and common outcome).
TH1D* BookHist(const char* name, bool logx, double xmax) {
  if (!logx) return new TH1D(name, "", 80, 0, xmax);
  const int nb = 100;
  std::vector<double> edges(nb + 1);
  edges[0] = 0.0;
  double lo = std::log10(0.5), hi = std::log10(xmax);
  for (int i = 1; i <= nb; ++i)
    edges[i] = std::pow(10, lo + (hi - lo) * (i - 1) / (nb - 1));
  return new TH1D(name, "", nb, edges.data());
}

// refDir restricts the loop to file names that exist in BOTH samples, so the two
// sides cover the same underlying events. Needed for HYDJet, where the default
// production has only 14 of the 100 files.
void Fill(const char* dir, const char* refDir, int nfiles, int maxEvents, TH1D* hAcc,
          TH1D* hRej, long& nEvt, long& nAcc, long& nRej) {
  std::vector<TString> files = ListRootFiles(dir, nfiles);
  if (strlen(refDir)) {
    std::vector<TString> ref = ListRootFiles(refDir, nfiles);
    KeepCommonFiles(files, ref);
  }
  for (const auto& path : files) {
    if (maxEvents > 0 && nEvt >= maxEvents) break;
    TFile* f = TFile::Open(path);
    if (!f || f->IsZombie()) continue;
    TTree* t = (TTree*)f->Get(kTreePath);
    if (!t) { f->Close(); continue; }
    std::vector<int>* rej = nullptr;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("allstub_isRejected", 1);
    t->SetBranchAddress("allstub_isRejected", &rej);
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      if (maxEvents > 0 && nEvt >= maxEvents) break;
      t->GetEntry(e);
      long a = 0, r = 0;
      for (size_t k = 0; k < rej->size(); ++k) (rej->at(k) ? r : a)++;
      hAcc->Fill(a);
      hRej->Fill(r);
      ++nEvt;
      nAcc += a;
      nRej += r;
    }
    f->Close();
    printf("[info] %s done (%ld events)\n", gSystem->BaseName(path), nEvt);
    fflush(stdout);
  }
}

// Write every figure as both PDF (for notes/talks) and PNG (for quick viewing).
void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

void SetLineStyle(TH1D* h, int color, int style) {
  h->SetLineColor(color);
  h->SetLineStyle(style);
  h->SetLineWidth(2);
}

void Report(const char* lab, TH1D* h, long nEvt, long ntot) {
  printf("  %-22s mean=%10.3f  RMS=%10.3f  total=%10ld  per event=%10.3f\n", lab,
         h->GetMean(), h->GetRMS(), ntot, nEvt ? (double)ntot / nEvt : 0);
}

}  // namespace

// commonIndicesOnly restricts both samples to file indices present in both. Off by
// default: the curves are unit-normalised, so unequal event counts are fine and
// using every available file gives the better-populated sample smoother statistics.
// Turn it on only if you specifically want both sides over identical events.
void compare_stubmult(int nfiles = 10, const char* dirDef = "", const char* dirDum = "",
                      bool logx = false, int maxEvents = -1, double xmax = 80.0,
                      const char* tag = "qedmumu", bool commonIndicesOnly = false) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  // Empty strings select the QED mumu sample; "pp" selects pythia pp.
  const char* dDef = strlen(dirDef) ? dirDef : (strcmp(tag, "pp") == 0 ? kPPDefault : kQEDDefault);
  const char* dDum = strlen(dirDum) ? dirDum : (strcmp(tag, "pp") == 0 ? kPPDummy : kQEDDummy);
  if (!strlen(dDef) || !strlen(dDum)) {
    printf("give the Default-stub and Dummy-stub ntuple locations (a directory of *.root, or one file)\n");
    return;
  }

  auto* hAccD = BookHist(Form("nstub_acc_def_%s", tag), logx, xmax);
  auto* hRejD = BookHist(Form("nstub_rej_def_%s", tag), logx, xmax);
  auto* hAccU = BookHist(Form("nstub_acc_dum_%s", tag), logx, xmax);
  auto* hRejU = BookHist(Form("nstub_rej_dum_%s", tag), logx, xmax);

  long nEvtD = 0, nAccD = 0, nRejD = 0, nEvtU = 0, nAccU = 0, nRejU = 0;
  const char* refD = commonIndicesOnly ? dDum : "";
  const char* refU = commonIndicesOnly ? dDef : "";
  printf("=== Default Stub ===\n");
  Fill(dDef, refD, nfiles, maxEvents, hAccD, hRejD, nEvtD, nAccD, nRejD);
  printf("=== Dummy Stub ===\n");
  Fill(dDum, refU, nfiles, maxEvents, hAccU, hRejU, nEvtU, nAccU, nRejU);

  printf("\n=========== stub multiplicity per event (%s) ===========\n", tag);
  printf("Default Stub  (%ld events)\n", nEvtD);
  Report("accepted", hAccD, nEvtD, nAccD);
  Report("rejected", hRejD, nEvtD, nRejD);
  printf("  rejected fraction: %.4f%%\n", 100.0 * nRejD / (nAccD + nRejD));
  printf("Dummy Stub  (%ld events)\n", nEvtU);
  Report("accepted", hAccU, nEvtU, nAccU);
  Report("rejected", hRejU, nEvtU, nRejU);
  printf("  rejected fraction: %.4f%%\n", 100.0 * nRejU / (nAccU + nRejU));

  // Normalise to unit area so all four are directly comparable as shapes.
  for (auto* h : {hAccD, hRejD, hAccU, hRejU}) {
    double n = h->Integral(0, h->GetNbinsX() + 1);
    if (n > 0) h->Scale(1.0 / n);
  }

  SetLineStyle(hAccD, kBlue + 1, 1);
  SetLineStyle(hRejD, kBlue + 1, 2);
  SetLineStyle(hAccU, kRed + 1, 1);
  SetLineStyle(hRejU, kRed + 1, 2);

  auto* c = new TCanvas(Form("c_stubmult_%s", tag), "", 900, 660);
  c->SetLogy();
  if (logx) c->SetLogx();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);

  hAccD->GetXaxis()->SetTitle("number of stubs per event");
  hAccD->GetYaxis()->SetTitle("fraction of events");
  hAccD->GetXaxis()->SetTitleSize(0.045);
  hAccD->GetYaxis()->SetTitleSize(0.045);
  hAccD->GetYaxis()->SetTitleOffset(1.35);

  // Push the ceiling well above the tallest peak so the legend sits in clear
  // space rather than on top of a distribution. With log y this costs nothing
  // but empty decades at the top.
  double peak = 0;
  for (auto* h : {hAccD, hRejD, hAccU, hRejU}) peak = std::max(peak, h->GetMaximum());
  hAccD->SetMaximum(peak * 3000);
  hAccD->SetMinimum(2e-6);
  hAccD->Draw("hist");
  hRejD->Draw("hist same");
  hAccU->Draw("hist same");
  hRejU->Draw("hist same");

  auto* leg = new TLegend(0.40, 0.73, 0.90, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.030);
  leg->AddEntry(hAccD, Form("Default  accepted  (#LTN#GT = %.1f)", (double)nAccD / nEvtD), "l");
  leg->AddEntry(hRejD, Form("Default  rejected  (#LTN#GT = %.1f)", (double)nRejD / nEvtD), "l");
  leg->AddEntry(hAccU, Form("Dummy   accepted  (#LTN#GT = %.1f)", (double)nAccU / nEvtU), "l");
  leg->AddEntry(hRejU, Form("Dummy   rejected  (#LTN#GT = %.1f)", (double)nRejU / nEvtU), "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.034);
  const char* sampleName = strcmp(tag, "pp") == 0        ? "pythia pp"
                           : strcmp(tag, "hydjet") == 0  ? "HYDJet PbPb"
                           : strcmp(tag, "epos") == 0    ? "EPOS pPb"
                                                         : "STARlight QED #mu#mu";
  tx.DrawLatex(0.13, 0.94,
               Form("%s -- stub multiplicity, accepted vs rejected", sampleName));

  // Rejected fraction of all stubs, quoted on the plot itself.
  double fracD = (nAccD + nRejD) > 0 ? 100.0 * nRejD / (nAccD + nRejD) : 0;
  double fracU = (nAccU + nRejU) > 0 ? 100.0 * nRejU / (nAccU + nRejU) : 0;
  TLatex tf;
  tf.SetNDC();
  tf.SetTextSize(0.030);
  tf.SetTextColor(kBlue + 1);
  tf.DrawLatex(0.40, 0.675, Form("Default: %.3f%% of stubs rejected", fracD));
  tf.SetTextColor(kRed + 1);
  tf.DrawLatex(0.40, 0.635, Form("Dummy:  %.3f%% of stubs rejected", fracU));

  SaveBoth(c, Form("../figures/stub_multiplicity_overlay_%s.pdf", tag));

  TFile* fout = TFile::Open(Form("stubmult_%s.root", tag), "RECREATE");
  for (auto* h : {hAccD, hRejD, hAccU, hRejU}) h->Write();
  fout->Close();
  printf("\nwrote ../figures/stub_multiplicity_overlay_%s.pdf\n", tag);
}
