// ---------------------------------------------------------------------------
// compare_trkquality.C
//
// Compare L1 track quality variables between the Default Stub and Dummy Stub
// productions. Unlike the tp_/matchtrk_ side, the trk_* collection is fully
// filled in both productions, including trk_bendchi2 and the chi2 family.
//
// The point of this script is to test whether the tracklet bend cuts were
// actually active in the dummy-stub production:
//
//   * If the bend cut was ON in both, the dummy trk_bendchi2 distribution is
//     bounded the same way the default one is.
//   * If the bend cut had been OFF for dummy, dummy would keep tracks the
//     default rejects, producing a visible high-bendchi2 tail.
//
// Usage: root -l -b -q 'compare_trkquality.C(10)'
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
    "/Users/wl33/Documents/DefaultStub/STARlight_QED_mumu";
const char* kDirDummy =
    "/Users/wl33/Documents/DummyStub/STARligt_QED_mumu";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

struct Q {
  TH1D *bendchi2, *chi2dof, *chi2rphi, *chi2rz, *nstub;
  TH1D *chi2rphidof, *chi2rzdof;
  TH1D* seed;
  long n = 0, ngen = 0, nfake = 0;
  double maxbend = -1;

  void Book(const char* tag) {
    bendchi2 = new TH1D(Form("bendchi2_%s", tag), ";trk_bendchi2;tracks", 400, 0, 20);
    chi2dof = new TH1D(Form("chi2dof_%s", tag), ";#chi^{2}/ndf;tracks", 200, 0, 20);
    chi2rphi = new TH1D(Form("chi2rphi_%s", tag), ";trk_chi2rphi;tracks", 400, 0, 100);
    chi2rz = new TH1D(Form("chi2rz_%s", tag), ";trk_chi2rz;tracks", 400, 0, 100);
    // trk_chi2rphi_dof / trk_chi2rz_dof are unfilled in these productions, so they
    // are recomputed here with the standard L1TrackNtupleMaker convention:
    //   ndofrphi = nstub - nPar + 2 = nstub - 2  (nPar = 4)
    //   ndofrz   = nstub - 2
    chi2rphidof = new TH1D(Form("chi2rphidof_%s", tag), ";#chi^{2}_{r#phi}/ndf;tracks", 200, 0, 20);
    chi2rzdof = new TH1D(Form("chi2rzdof_%s", tag), ";#chi^{2}_{rz}/ndf;tracks", 200, 0, 20);
    nstub = new TH1D(Form("nstub_%s", tag), ";trk_nstub;tracks", 12, -0.5, 11.5);
    seed = new TH1D(Form("seed_%s", tag), ";trk_seed;tracks", 24, -0.5, 23.5);
  }
};

const int kNPar = 4;  // L1Tk_nPar in L1TrackHitNtupleMaker_cfi.py

void Report(const char* label, Q& q) {
  printf("\n--- %s ---\n", label);
  printf("  tracks: %ld   genuine: %ld (%.1f%%)   fake: %ld (%.1f%%)\n", q.n, q.ngen,
         q.n ? 100.0 * q.ngen / q.n : 0, q.nfake, q.n ? 100.0 * q.nfake / q.n : 0);
  double p[5] = {0.5, 0.9, 0.99, 0.999, 1.0}, v[5];
  q.bendchi2->GetQuantiles(4, v, p);
  printf("  trk_bendchi2: mean %.3f  median %.3f  p90 %.3f  p99 %.3f  p99.9 %.3f  max %.3f\n",
         q.bendchi2->GetMean(), v[0], v[1], v[2], v[3], q.maxbend);
  printf("    overflow (>20): %.0f tracks (%.4f%%)\n",
         q.bendchi2->GetBinContent(q.bendchi2->GetNbinsX() + 1),
         q.n ? 100.0 * q.bendchi2->GetBinContent(q.bendchi2->GetNbinsX() + 1) / q.n : 0);
  q.chi2dof->GetQuantiles(4, v, p);
  printf("  trk_chi2_dof: mean %.3f  median %.3f  p90 %.3f  p99 %.3f\n", q.chi2dof->GetMean(),
         v[0], v[1], v[2]);
  printf("  trk_chi2rphi mean %.2f   trk_chi2rz mean %.2f   <nstub> %.3f\n",
         q.chi2rphi->GetMean(), q.chi2rz->GetMean(), q.nstub->GetMean());
  printf("  seed composition:");
  for (int i = 1; i <= q.seed->GetNbinsX(); ++i)
    if (q.seed->GetBinContent(i) > 0)
      printf(" %d:%.1f%%", (int)q.seed->GetBinLowEdge(i) + 0,
             100.0 * q.seed->GetBinContent(i) / q.n);
  printf("\n");
}

void Fill(const char* dir, int nfiles, Q& q) {
  for (int i = 1; i <= nfiles; ++i) {
    TFile* f = TFile::Open(Form("%s/L1TrackHitNtuple_UPC_v4_%d.root", dir, i));
    if (!f || f->IsZombie()) continue;
    TTree* t = (TTree*)f->Get(kTreePath);
    if (!t) { f->Close(); continue; }
    std::vector<float> *pt = nullptr, *bend = nullptr, *c2d = nullptr, *c2rphi = nullptr,
                       *c2rz = nullptr;
    std::vector<int> *ns = nullptr, *sd = nullptr, *gen = nullptr, *fake = nullptr;
    t->SetBranchStatus("*", 0);
    for (auto b : {"trk_pt", "trk_bendchi2", "trk_chi2_dof", "trk_chi2rphi", "trk_chi2rz",
                   "trk_nstub", "trk_seed", "trk_genuine", "trk_fake"})
      t->SetBranchStatus(b, 1);
    t->SetBranchAddress("trk_pt", &pt);
    t->SetBranchAddress("trk_bendchi2", &bend);
    t->SetBranchAddress("trk_chi2_dof", &c2d);
    t->SetBranchAddress("trk_chi2rphi", &c2rphi);
    t->SetBranchAddress("trk_chi2rz", &c2rz);
    t->SetBranchAddress("trk_nstub", &ns);
    t->SetBranchAddress("trk_seed", &sd);
    t->SetBranchAddress("trk_genuine", &gen);
    t->SetBranchAddress("trk_fake", &fake);
    for (Long64_t e = 0; e < t->GetEntries(); ++e) {
      t->GetEntry(e);
      for (size_t k = 0; k < pt->size(); ++k) {
        q.n++;
        // trk_genuine is 0/1. trk_fake follows the ntuple maker's convention:
        // 0 = fake (no TP), 1 = genuine primary, 2 = combinatoric/pileup.
        if (gen->size() == pt->size() && gen->at(k) == 1) q.ngen++;
        if (fake->size() == pt->size() && fake->at(k) == 0) q.nfake++;
        if (bend->size() == pt->size()) {
          q.bendchi2->Fill(bend->at(k));
          if (bend->at(k) > q.maxbend) q.maxbend = bend->at(k);
        }
        if (c2d->size() == pt->size()) q.chi2dof->Fill(c2d->at(k));
        if (c2rphi->size() == pt->size()) q.chi2rphi->Fill(c2rphi->at(k));
        if (c2rz->size() == pt->size()) q.chi2rz->Fill(c2rz->at(k));
        if (ns->size() == pt->size()) {
          q.nstub->Fill(ns->at(k));
          int ndofrphi = ns->at(k) - kNPar + 2;
          int ndofrz = ns->at(k) - 2;
          if (ndofrphi > 0 && c2rphi->size() == pt->size())
            q.chi2rphidof->Fill(c2rphi->at(k) / ndofrphi);
          if (ndofrz > 0 && c2rz->size() == pt->size())
            q.chi2rzdof->Fill(c2rz->at(k) / ndofrz);
        }
        if (sd->size() == pt->size()) q.seed->Fill(sd->at(k));
      }
    }
    f->Close();
    printf("[info] file %d done\n", i);
    fflush(stdout);
  }
}

}  // namespace

void compare_trkquality(int nfiles = 10, const char* dirDef = "", const char* dirDum = "",
                        const char* outname = "trkquality_qed_mumu.root",
                        const char* tag = "") {
  if (strlen(dirDef)) kDirDefault = dirDef;
  if (strlen(dirDum)) kDirDummy = dirDum;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  Q qd, qu;
  qd.Book("def");
  qu.Book("dum");
  printf("=== reading Default Stub ===\n");
  Fill(kDirDefault, nfiles, qd);
  printf("=== reading Dummy Stub ===\n");
  Fill(kDirDummy, nfiles, qu);

  printf("\n============ L1 track quality, %d files/sample ============\n", nfiles);
  Report("Default Stub", qd);
  Report("Dummy Stub", qu);

  printf("\n--- interpretation of trk_bendchi2 ---\n");
  printf("  If the tracklet bend cuts were disabled for the dummy sample, dummy would\n");
  printf("  retain tracks the default rejects and show a heavier high-bendchi2 tail.\n");
  printf("  Similar distributions => the bend cuts were active in BOTH productions.\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (auto* h : {qd.bendchi2, qd.chi2dof, qd.chi2rphi, qd.chi2rz, qd.chi2rphidof,
                  qd.chi2rzdof, qd.nstub, qd.seed, qu.bendchi2, qu.chi2dof, qu.chi2rphi,
                  qu.chi2rz, qu.chi2rphidof, qu.chi2rzdof, qu.nstub, qu.seed})
    h->Write();
  fout->Close();

  auto* c = new TCanvas("c_bend", "", 820, 620);
  c->SetLogy();
  c->SetLeftMargin(0.13);
  qd.bendchi2->SetLineColor(kBlue + 1);
  qd.bendchi2->SetLineWidth(2);
  qu.bendchi2->SetLineColor(kRed + 1);
  qu.bendchi2->SetLineWidth(2);
  qd.bendchi2->GetXaxis()->SetRangeUser(0, 10);
  qd.bendchi2->Draw("hist");
  qu.bendchi2->Draw("hist same");
  auto* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(qd.bendchi2, "Default Stub", "l");
  leg->AddEntry(qu.bendchi2, "Dummy Stub", "l");
  leg->Draw();
  SaveBoth(c, Form("../figures/trk_bendchi2%s.pdf", tag));
  printf("\nwrote trkquality_qed_mumu.root and ../figures/trk_bendchi2.pdf\n");
}
