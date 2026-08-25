// ---------------------------------------------------------------------------
// plot_tpspectrum.C
//
// Truth (tracking particle) pT spectrum within |eta| < 2.4, Default vs Dummy.
//
// Purpose: find out whether the sub-GeV depletion in the stored TP collection is
// a hard cut (a step edge at some pT) or a gradual turn-on (a smooth efficiency
// roll-off from the >=3 stub / >=3 layer preselection).
//
// A real EPOS pPb or HYDJet charged spectrum peaks near 0.3-0.5 GeV and falls
// steeply, so the overwhelming majority of charged particles are below 1 GeV.
// The stored collections contain ~1% below 1 GeV, which is a large depletion --
// this plot shows its shape.
//
// Note on geometry: at B = 3.8 T the helix radius is R[cm] ~ 87.7 x pT[GeV], so
// a 0.5 GeV particle reaches 2R ~ 88 cm, i.e. barrel layers 1-5 (r up to 86 cm).
// Reaching 3 layers (r = 51 cm) only needs pT > ~0.29 GeV. So the 3-layer
// requirement alone does NOT explain a depletion setting in around 1 GeV.
//
// Usage: root -l -b -q 'plot_tpspectrum.C("/path/def","/path/dum","EPOS pPb","epospb",50)'
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>

namespace {

const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";
const double kEtaMax = 2.4;

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

// Fine binning at low pT, where the question is.
TH1D* Book(const char* name) { return new TH1D(name, "", 500, 0, 5.0); }

long Fill(const char* dir, int nfiles, TH1D* h, TH1D* hns) {
  long nEvt = 0;
  for (int i = 1; i <= nfiles; ++i) {
    TString p = Form("%s/L1TrackHitNtuple_UPC_v4_%d.root", dir, i);
    if (gSystem->AccessPathName(p)) continue;
    TFile* f = TFile::Open(p);
    if (!f || f->IsZombie()) continue;
    TTree* t = (TTree*)f->Get(kTreePath);
    if (!t || !t->GetBranch("tp_pt")) { if (f) f->Close(); continue; }
    std::vector<float> *pt = nullptr, *eta = nullptr;
    std::vector<int>* ns = nullptr;
    t->SetBranchStatus("*", 0);
    for (auto b : {"tp_pt", "tp_eta", "tp_nstub"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("tp_pt", &pt);
    t->SetBranchAddress("tp_eta", &eta);
    t->SetBranchAddress("tp_nstub", &ns);
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      ++nEvt;
      size_t n = std::min(pt->size(), eta->size());
      for (size_t k = 0; k < n; ++k) {
        if (std::abs(eta->at(k)) >= kEtaMax) continue;
        h->Fill(pt->at(k));
        // <nstub> vs pT: if the depletion is the stub requirement biting, the
        // mean stub count should fall towards 3 as pT drops.
        if (ns && k < ns->size()) hns->Fill(pt->at(k), ns->at(k));
      }
    }
    f->Close();
  }
  return nEvt;
}

}  // namespace

void plot_tpspectrum(const char* dirDef, const char* dirDum, const char* sample,
                     const char* tag, int nfiles = 50) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  TH1D* hD = Book(Form("tpspec_def_%s", tag));
  TH1D* hU = Book(Form("tpspec_dum_%s", tag));
  auto* pD = new TProfile(Form("nstub_vs_pt_def_%s", tag), "", 500, 0, 5.0);
  auto* pU = new TProfile(Form("nstub_vs_pt_dum_%s", tag), "", 500, 0, 5.0);

  long nD = Fill(dirDef, nfiles, hD, pD);
  long nU = Fill(dirDum, nfiles, hU, pU);
  printf("[info] %s: default %ld events, dummy %ld events\n", sample, nD, nU);

  // Per-event normalisation so the two are directly comparable.
  if (nD > 0) hD->Scale(1.0 / nD);
  if (nU > 0) hU->Scale(1.0 / nU);

  printf("\n=== %s: truth pT spectrum, |eta| < %.1f (TPs per event per 10 MeV) ===\n",
         sample, kEtaMax);
  printf("%14s %14s %14s %10s\n", "pT [GeV]", "default", "dummy", "dum/def");
  for (double thr = 0.05; thr < 2.0; thr += (thr < 0.5 ? 0.05 : 0.1)) {
    int b = hD->FindBin(thr);
    double a = hD->GetBinContent(b), c2 = hU->GetBinContent(b);
    printf("%14.2f %14.6f %14.6f %10s\n", thr, a, c2,
           a > 0 ? Form("%.2f", c2 / a) : "-");
  }
  printf("\nfirst non-empty bin: default %.3f GeV, dummy %.3f GeV\n",
         hD->GetBinLowEdge(hD->FindFirstBinAbove(0)),
         hU->GetBinLowEdge(hU->FindFirstBinAbove(0)));
  printf("peak of spectrum:    default %.3f GeV, dummy %.3f GeV\n",
         hD->GetBinCenter(hD->GetMaximumBin()), hU->GetBinCenter(hU->GetMaximumBin()));

  hD->SetLineColor(kBlue + 1);
  hD->SetLineWidth(2);
  hU->SetLineColor(kRed + 1);
  hU->SetLineWidth(2);

  auto* c = new TCanvas(Form("c_tpspec_%s", tag), "", 900, 660);
  c->SetLogy();
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.13);
  hU->GetXaxis()->SetTitle("truth p_{T} [GeV]");
  hU->GetYaxis()->SetTitle("tracking particles / event / 10 MeV");
  hU->GetXaxis()->SetTitleSize(0.045);
  hU->GetYaxis()->SetTitleSize(0.045);
  hU->GetYaxis()->SetTitleOffset(1.35);
  hU->SetMaximum(hU->GetMaximum() * 60);
  hU->Draw("hist");
  hD->Draw("hist same");

  auto* leg = new TLegend(0.52, 0.74, 0.90, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.032);
  leg->AddEntry(hU, "Dummy Stub", "l");
  leg->AddEntry(hD, "Default Stub", "l");
  leg->Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.13, 0.94,
               Form("%s -- stored truth p_{T} spectrum, |#eta| < 2.4", sample));
  SaveBoth(c, Form("../figures/tp_spectrum_%s.pdf", tag));

  // <tp_nstub> vs pT -- shows whether the stub requirement is what bites.
  auto* c2 = new TCanvas(Form("c_nstub_%s", tag), "", 900, 660);
  c2->SetGridy();
  c2->SetLeftMargin(0.13);
  c2->SetBottomMargin(0.13);
  pU->SetLineColor(kRed + 1);  pU->SetLineWidth(2); pU->SetMarkerColor(kRed + 1);
  pD->SetLineColor(kBlue + 1); pD->SetLineWidth(2); pD->SetMarkerColor(kBlue + 1);
  pU->GetXaxis()->SetRangeUser(0, 5);
  pU->SetMinimum(0);
  pU->GetXaxis()->SetTitle("truth p_{T} [GeV]");
  pU->GetYaxis()->SetTitle("#LTtp_nstub#GT");
  pU->GetXaxis()->SetTitleSize(0.045);
  pU->GetYaxis()->SetTitleSize(0.045);
  pU->GetYaxis()->SetTitleOffset(1.35);
  pU->Draw("hist");
  pD->Draw("hist same");
  auto* l3 = new TLine(0, 3, 5, 3);
  l3->SetLineStyle(2); l3->SetLineColor(kGray + 2); l3->Draw();
  auto* leg2 = new TLegend(0.52, 0.20, 0.90, 0.34);
  leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.032);
  leg2->AddEntry(pU, "Dummy Stub", "l");
  leg2->AddEntry(pD, "Default Stub", "l");
  leg2->AddEntry(l3, "TP_minNStub = 3", "l");
  leg2->Draw();
  TLatex tx2;
  tx2.SetNDC(); tx2.SetTextSize(0.034);
  tx2.DrawLatex(0.13, 0.94, Form("%s -- #LTstubs per truth particle#GT vs p_{T}", sample));
  SaveBoth(c2, Form("../figures/tp_nstub_vs_pt_%s.pdf", tag));

  TFile* fout = TFile::Open(Form("tpspectrum_%s.root", tag), "RECREATE");
  hD->Write(); hU->Write(); pD->Write(); pU->Write();
  fout->Close();
  printf("\nwrote ../figures/tp_spectrum_%s.pdf and tp_nstub_vs_pt_%s.pdf\n", tag, tag);
}
