// ---------------------------------------------------------------------------
// compare_allstub_pp.C
//
// Stub-level comparison between the Default Stub and Dummy Stub productions for
// the pythia pp sample.
//
// The dummy pp ntuples contain ONLY the allstub_* branches (17 of them) -- no
// tp_*, trk_* or matchtrk_*. So no efficiency/resolution/chi2 comparison is
// possible for pp; this compares every allstub_* branch that is usable.
//
// TWO BRANCHES ARE NOT USABLE IN pp, verified before writing this:
//
//   allstub_genuine, allstub_matchTP_{pt,eta,phi,pdgid}
//     * Default pp : vectors are the right length but every value is -999 / 0,
//       i.e. the stub MC-truth association was never run.
//     * Dummy pp   : the vectors are SHORTER than allstub_x (e.g. 1,395,962 vs
//       1,494,751 over 20 events), so they cannot be safely indexed against the
//       position branches.
//   Both are therefore skipped. This removes the "what pT of particle made this
//   stub" study, which is unfortunately the most interesting one.
//
// Layer numbering: allstub_layer runs 1-6 for barrel AND 1-5 for disks, so it is
// only meaningful combined with allstub_isBarrel. Handled here as a combined
// layer/disk index: 0-5 = barrel L1-L6, 6-10 = disk D1-D5.
//
// Usage: root -l -b -q 'compare_allstub_pp.C(1000)'   // events per sample
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>

namespace {

// Write every figure as both PDF and PNG.
void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

const char* kDirDefault =
    "/Users/wl33/Documents/DefaultStub/pythia_pp";
const char* kDirDummy =
    "/Users/wl33/Documents/DummyStub/pythia_pp";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

const char* kLayerName[11] = {"L1", "L2", "L3", "L4", "L5", "L6",
                              "D1", "D2", "D3", "D4", "D5"};

struct Set {
  TH1D *mult, *layerdisk, *r, *z, *absz, *rz_ps;
  long nEvt = 0, nStub = 0, nRej = 0, nPS = 0, nBarrel = 0, nTilt = 0;
  long perLD[11] = {0};
  long rejPerLD[11] = {0};
  long psPerLD[11] = {0};

  void Book(const char* tag) {
    // Common binning for both samples so ratios are well defined.
    mult = new TH1D(Form("mult_%s", tag), "", 300, 0, 150000);
    layerdisk = new TH1D(Form("ld_%s", tag), "", 11, -0.5, 10.5);
    r = new TH1D(Form("r_%s", tag), "", 120, 0, 120);
    z = new TH1D(Form("z_%s", tag), "", 140, -300, 300);
    absz = new TH1D(Form("absz_%s", tag), "", 120, 0, 300);
    for (auto h : {mult, layerdisk, r, z, absz}) h->Sumw2();
  }
};

void Fill(const char* dir, int maxEvents, Set& s) {
  for (int i = 1; i <= 2; ++i) {
    if (s.nEvt >= maxEvents) break;
    TFile* f = TFile::Open(Form("%s/L1TrackHitNtuple_UPC_v4_%d.root", dir, i));
    if (!f || f->IsZombie()) continue;
    TTree* t = (TTree*)f->Get(kTreePath);
    if (!t) { f->Close(); continue; }

    std::vector<float> *x = nullptr, *y = nullptr, *zz = nullptr;
    std::vector<int> *lay = nullptr, *isPS = nullptr, *isBar = nullptr, *isTilt = nullptr,
                     *rej = nullptr;
    t->SetBranchStatus("*", 0);
    for (auto b : {"allstub_x", "allstub_y", "allstub_z", "allstub_layer",
                   "allstub_isPSmodule", "allstub_isBarrel", "allstub_isTiltedBarrel",
                   "allstub_isRejected"})
      t->SetBranchStatus(b, 1);
    t->SetBranchAddress("allstub_x", &x);
    t->SetBranchAddress("allstub_y", &y);
    t->SetBranchAddress("allstub_z", &zz);
    t->SetBranchAddress("allstub_layer", &lay);
    t->SetBranchAddress("allstub_isPSmodule", &isPS);
    t->SetBranchAddress("allstub_isBarrel", &isBar);
    t->SetBranchAddress("allstub_isTiltedBarrel", &isTilt);
    t->SetBranchAddress("allstub_isRejected", &rej);

    for (Long64_t e = 0; e < t->GetEntries() && s.nEvt < maxEvents; ++e) {
      t->GetEntry(e);
      size_t n = x ? x->size() : 0;
      s.mult->Fill(n);
      ++s.nEvt;
      for (size_t k = 0; k < n; ++k) {
        ++s.nStub;
        s.r->Fill(std::hypot(x->at(k), y->at(k)));
        s.z->Fill(zz->at(k));
        s.absz->Fill(std::abs(zz->at(k)));

        int barrel = (isBar && k < isBar->size()) ? isBar->at(k) : -1;
        int L = (lay && k < lay->size()) ? lay->at(k) : -1;
        int ld = -1;
        if (barrel == 1 && L >= 1 && L <= 6) ld = L - 1;        // barrel L1-L6 -> 0-5
        else if (barrel == 0 && L >= 1 && L <= 5) ld = 5 + L;   // disk D1-D5   -> 6-10
        if (ld >= 0) {
          s.layerdisk->Fill(ld);
          s.perLD[ld]++;
          if (rej && k < rej->size() && rej->at(k)) s.rejPerLD[ld]++;
          if (isPS && k < isPS->size() && isPS->at(k)) s.psPerLD[ld]++;
        }
        if (isPS && k < isPS->size() && isPS->at(k)) ++s.nPS;
        if (barrel == 1) ++s.nBarrel;
        if (isTilt && k < isTilt->size() && isTilt->at(k)) ++s.nTilt;
        if (rej && k < rej->size() && rej->at(k)) ++s.nRej;
      }
    }
    f->Close();
    printf("[info] file %d done (%ld events)\n", i, s.nEvt);
    fflush(stdout);
  }
}

void Report(const char* lab, Set& s) {
  printf("\n--- %s ---\n", lab);
  printf("  events          : %ld\n", s.nEvt);
  printf("  stubs           : %ld   (%.1f per event)\n", s.nStub,
         s.nEvt ? (double)s.nStub / s.nEvt : 0);
  printf("  rejected        : %.3f%%\n", s.nStub ? 100.0 * s.nRej / s.nStub : 0);
  printf("  PS modules      : %.2f%%\n", s.nStub ? 100.0 * s.nPS / s.nStub : 0);
  printf("  barrel          : %.2f%%   (tilted barrel: %.2f%%)\n",
         s.nStub ? 100.0 * s.nBarrel / s.nStub : 0, s.nStub ? 100.0 * s.nTilt / s.nStub : 0);
  printf("  <r>             : %.2f cm    <|z|> : %.2f cm\n", s.r->GetMean(),
         s.absz->GetMean());
}

void Overlay(TH1D* hd, TH1D* hu, const char* xtitle, const char* title, const char* out,
             bool logy, bool ratio, double xlo = 0, double xhi = 0) {
  auto* d = (TH1D*)hd->Clone(Form("%s_c", hd->GetName()));
  auto* u = (TH1D*)hu->Clone(Form("%s_c", hu->GetName()));
  double nd = d->Integral(0, d->GetNbinsX() + 1), nu = u->Integral(0, u->GetNbinsX() + 1);
  if (nd > 0) d->Scale(1.0 / nd);
  if (nu > 0) u->Scale(1.0 / nu);
  d->SetLineColor(kBlue + 1); d->SetLineWidth(2);
  u->SetLineColor(kRed + 1);  u->SetLineWidth(2);

  auto* c = new TCanvas(Form("c_%s", out), "", 860, ratio ? 760 : 640);
  TPad *p1 = nullptr, *p2 = nullptr;
  if (ratio) {
    p1 = new TPad(Form("p1%s", out), "", 0, 0.30, 1, 1);
    p2 = new TPad(Form("p2%s", out), "", 0, 0.0, 1, 0.30);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.13);
    p2->SetTopMargin(0.02); p2->SetBottomMargin(0.34);
    p2->SetLeftMargin(0.13); p2->SetGridy();
    if (logy) p1->SetLogy();
    p1->Draw(); p2->Draw(); p1->cd();
  } else {
    c->SetLeftMargin(0.13); c->SetBottomMargin(0.13);
    if (logy) c->SetLogy();
  }

  if (xhi > xlo) d->GetXaxis()->SetRangeUser(xlo, xhi);
  if (ratio) d->GetXaxis()->SetLabelSize(0);
  else { d->GetXaxis()->SetTitle(xtitle); d->GetXaxis()->SetTitleSize(0.045); }
  d->GetYaxis()->SetTitle("fraction of stubs");
  d->GetYaxis()->SetTitleSize(ratio ? 0.05 : 0.045);
  d->GetYaxis()->SetTitleOffset(1.3);
  d->SetMaximum(std::max(d->GetMaximum(), u->GetMaximum()) * (logy ? 8 : 1.45));
  if (logy) d->SetMinimum(1e-7);
  d->Draw("hist"); u->Draw("hist same");

  auto* leg = new TLegend(0.56, 0.74, 0.88, 0.88);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.036);
  leg->AddEntry(d, "Default Stub", "l");
  leg->AddEntry(u, "Dummy Stub", "l");
  leg->Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.037);
  tx.DrawLatex(0.13, 0.94, title);

  if (ratio) {
    p2->cd();
    auto* rr = (TH1D*)u->Clone(Form("r%s", out));
    rr->Divide(d);
    rr->SetLineColor(kBlack); rr->SetMarkerColor(kBlack);
    rr->SetMarkerStyle(20); rr->SetMarkerSize(0.8);
    if (xhi > xlo) rr->GetXaxis()->SetRangeUser(xlo, xhi);
    rr->GetYaxis()->SetRangeUser(0, 3.0);
    rr->GetYaxis()->SetTitle("dummy / default");
    rr->GetYaxis()->SetNdivisions(505);
    rr->GetYaxis()->SetTitleSize(0.115); rr->GetYaxis()->SetTitleOffset(0.52);
    rr->GetYaxis()->SetLabelSize(0.10);
    rr->GetXaxis()->SetTitle(xtitle);
    rr->GetXaxis()->SetTitleSize(0.125); rr->GetXaxis()->SetTitleOffset(1.15);
    rr->GetXaxis()->SetLabelSize(0.10);
    rr->Draw("ep");
  }
  SaveBoth(c, out);
}

}  // namespace

void compare_allstub_pp(int maxEvents = 1000) {
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Set D, U;
  D.Book("def");
  U.Book("dum");

  printf("=== Default Stub pp ===\n");
  Fill(kDirDefault, maxEvents, D);
  printf("=== Dummy Stub pp ===\n");
  Fill(kDirDummy, maxEvents, U);

  printf("\n============ pythia pp, allstub_* comparison ============\n");
  printf("NOTE: allstub_genuine and allstub_matchTP_* are unusable in pp\n");
  printf("      (default: all -999/0; dummy: vectors misaligned with allstub_x)\n");
  Report("Default Stub", D);
  Report("Dummy Stub", U);
  printf("\n  stub multiplicity ratio (dummy/default): %.2fx\n",
         (D.nStub && D.nEvt && U.nEvt) ? ((double)U.nStub / U.nEvt) / ((double)D.nStub / D.nEvt) : 0);

  printf("\n--- stubs per event by layer/disk ---\n");
  printf("  %6s %12s %12s %8s | %9s %9s | %9s %9s\n", "", "def/evt", "dum/evt", "ratio",
         "PS(def)", "PS(dum)", "rej(def)", "rej(dum)");
  for (int i = 0; i < 11; ++i) {
    double a = (double)D.perLD[i] / D.nEvt, b = (double)U.perLD[i] / U.nEvt;
    if (D.perLD[i] == 0 && U.perLD[i] == 0) continue;
    printf("  %6s %12.1f %12.1f %8s | %8.2f%% %8.2f%% | %8.3f%% %8.3f%%\n", kLayerName[i], a,
           b, a > 0 ? Form("%.2f", b / a) : "-",
           D.perLD[i] ? 100.0 * D.psPerLD[i] / D.perLD[i] : 0,
           U.perLD[i] ? 100.0 * U.psPerLD[i] / U.perLD[i] : 0,
           D.perLD[i] ? 100.0 * D.rejPerLD[i] / D.perLD[i] : 0,
           U.perLD[i] ? 100.0 * U.rejPerLD[i] / U.perLD[i] : 0);
  }

  for (int i = 0; i < 11; ++i) {
    D.layerdisk->GetXaxis()->SetBinLabel(i + 1, kLayerName[i]);
    U.layerdisk->GetXaxis()->SetBinLabel(i + 1, kLayerName[i]);
  }

  Overlay(D.mult, U.mult, "stubs per event", "pythia pp -- stub multiplicity",
          "../figures/pp_allstub_mult.pdf", true, false, 0, 130000);
  Overlay(D.layerdisk, U.layerdisk, "layer / disk",
          "pythia pp -- stub fraction by layer and disk",
          "../figures/pp_allstub_layerdisk.pdf", false, true);
  Overlay(D.r, U.r, "stub r [cm]", "pythia pp -- stub radius",
          "../figures/pp_allstub_r.pdf", false, true);
  Overlay(D.z, U.z, "stub z [cm]", "pythia pp -- stub z",
          "../figures/pp_allstub_z.pdf", false, true);

  TFile* fout = TFile::Open("allstub_pp.root", "RECREATE");
  for (auto* h : {D.mult, D.layerdisk, D.r, D.z, D.absz, U.mult, U.layerdisk, U.r, U.z,
                  U.absz})
    h->Write();
  fout->Close();
  printf("\nfigures written to ../figures/, histograms to allstub_pp.root\n");
}
