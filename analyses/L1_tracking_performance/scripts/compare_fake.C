// ---------------------------------------------------------------------------
// compare_fake.C
//
// L1 track fake rate, Default Stub vs Dummy Stub. The companion to compare_eff.C:
// efficiency asks "were the real particles found", this asks "are the found
// tracks real". Both are needed -- a configuration can buy efficiency purely by
// producing more combinatorial tracks.
//
//   fake rate    = N(trk_fake == 0)    / N(tracks)   -- no associated TP at all
//   non-genuine  = N(trk_genuine == 0) / N(tracks)   -- looser: includes tracks
//                  that have a TP but whose stubs do not all come from it
//
// trk_fake convention from the ntuplizer: 0 = fake (no TP), 1 = genuine primary,
// 2 = combinatoric/pileup. With no pileup overlay in these samples only 0 and 1
// occur, so "fake" and "non-genuine" differ solely through the isGenuine stub
// requirement -- both are reported because that gap is itself informative.
//
// Denominator is every reconstructed track. Same binning, acceptance and
// conventions as compare_eff.C, so the two sets of plots overlay directly.
//
// Usage:
//   root -l -b -q 'compare_fake.C(50)'                                // QED mumu
//   root -l -b -q 'compare_fake.C(100,"fake_hydjet.root","<def>","<dum>")'
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

// Set on the command line (arguments 3 and 4); a directory of *.root or one file.
const char* kDirDefault = "";
const char* kDirDummy = "";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

const double kPtCut = 2.0;    // matches compare_eff.C
const double kPtCut1 = 1.0;
const double kEtaAcc = 2.4;

struct Trk {
  float pt, eta, phi, z0;
  int fake, genuine;
  bool isFake() const { return fake == 0; }
  bool isNonGenuine() const { return genuine == 0; }
};

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *z0 = nullptr;
  std::vector<int> *fake = nullptr, *gen = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree || !tree->GetBranch("trk_fake")) return false;
    tree->SetBranchStatus("*", 0);
    for (auto b : {"trk_pt", "trk_eta", "trk_phi", "trk_z0", "trk_fake", "trk_genuine"})
      tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("trk_pt", &pt);
    tree->SetBranchAddress("trk_eta", &eta);
    tree->SetBranchAddress("trk_phi", &phi);
    tree->SetBranchAddress("trk_z0", &z0);
    tree->SetBranchAddress("trk_fake", &fake);
    tree->SetBranchAddress("trk_genuine", &gen);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr; tree = nullptr;
  }
  void Flatten(std::vector<Trk>& out) const {
    out.clear();
    size_t n = pt->size();
    for (size_t i = 0; i < n; ++i) {
      Trk t;
      t.pt = pt->at(i);  t.eta = eta->at(i);  t.phi = phi->at(i);
      t.z0 = (z0->size() == n) ? z0->at(i) : -999;
      t.fake = (fake->size() == n) ? fake->at(i) : 1;
      t.genuine = (gen->size() == n) ? gen->at(i) : 1;
      out.push_back(t);
    }
  }
};

// den = all tracks; numF = fake (no TP); numG = non-genuine.
struct Hists {
  TH1D *den_pt, *numF_pt, *numG_pt;
  TH1D *den_eta1, *numF_eta1, *numG_eta1, *den_phi1, *numF_phi1, *numG_phi1;
  TH1D *den_eta2, *numF_eta2, *numG_eta2, *den_phi2, *numF_phi2, *numG_phi2;
  TH1D *den_z02, *numF_z02, *numG_z02;

  TH1D* B(const char* w, const char* v, const char* tag, int nb, double lo, double hi) {
    return new TH1D(Form("%s_%s_%s", w, v, tag), "", nb, lo, hi);
  }
  TH1D* Bv(const char* w, const char* v, const char* tag, int nb, const double* e) {
    return new TH1D(Form("%s_%s_%s", w, v, tag), "", nb, e);
  }

  void Book(const char* tag, int npt, const double* ptbins) {
    den_pt  = Bv("den", "pt", tag, npt, ptbins);
    numF_pt = Bv("numF", "pt", tag, npt, ptbins);
    numG_pt = Bv("numG", "pt", tag, npt, ptbins);
    const double P = TMath::Pi();
    den_eta1  = B("den", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    numF_eta1 = B("numF", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    numG_eta1 = B("numG", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    den_phi1  = B("den", "phi1", tag, 32, -P, P);
    numF_phi1 = B("numF", "phi1", tag, 32, -P, P);
    numG_phi1 = B("numG", "phi1", tag, 32, -P, P);
    den_eta2  = B("den", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    numF_eta2 = B("numF", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    numG_eta2 = B("numG", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    den_phi2  = B("den", "phi2", tag, 32, -P, P);
    numF_phi2 = B("numF", "phi2", tag, 32, -P, P);
    numG_phi2 = B("numG", "phi2", tag, 32, -P, P);
    den_z02   = B("den", "z02", tag, 40, -20, 20);
    numF_z02  = B("numF", "z02", tag, 40, -20, 20);
    numG_z02  = B("numG", "z02", tag, 40, -20, 20);
    for (auto h : All()) h->Sumw2();
  }

  std::vector<TH1D*> All() {
    return {den_pt, numF_pt, numG_pt,
            den_eta1, numF_eta1, numG_eta1, den_phi1, numF_phi1, numG_phi1,
            den_eta2, numF_eta2, numG_eta2, den_phi2, numF_phi2, numG_phi2,
            den_z02, numF_z02, numG_z02};
  }

  void Fill(const Trk& t) {
    bool F = t.isFake(), G = t.isNonGenuine();
    den_pt->Fill(t.pt);
    if (F) numF_pt->Fill(t.pt);
    if (G) numG_pt->Fill(t.pt);
    if (t.pt > kPtCut1) {
      den_eta1->Fill(t.eta);
      if (F) numF_eta1->Fill(t.eta);
      if (G) numG_eta1->Fill(t.eta);
      if (std::abs(t.eta) < kEtaAcc) {
        den_phi1->Fill(t.phi);
        if (F) numF_phi1->Fill(t.phi);
        if (G) numG_phi1->Fill(t.phi);
      }
    }
    if (t.pt > kPtCut) {
      den_eta2->Fill(t.eta);
      if (F) numF_eta2->Fill(t.eta);
      if (G) numG_eta2->Fill(t.eta);
      if (std::abs(t.eta) < kEtaAcc) {
        den_phi2->Fill(t.phi);  den_z02->Fill(t.z0);
        if (F) { numF_phi2->Fill(t.phi); numF_z02->Fill(t.z0); }
        if (G) { numG_phi2->Fill(t.phi); numG_z02->Fill(t.z0); }
      }
    }
  }
};

struct Totals {
  long nFiles = 0, nEvents = 0, nTrk = 0, nFake = 0, nNonGen = 0;
  long nTrkHi = 0, nFakeHi = 0, nNonGenHi = 0;
};

void Run(const char* dir, int nfiles, Hists& h, Totals& tot) {
  std::vector<Trk> trks;
  for (const auto& path : ListRootFiles(dir, nfiles)) {
    Reader R;
    if (!R.Open(path)) { R.Close(); continue; }
    ++tot.nFiles;
    for (Long64_t e = 0; e < R.tree->GetEntries(); ++e) {
      R.tree->GetEntry(e);
      R.Flatten(trks);
      ++tot.nEvents;
      for (const auto& t : trks) {
        h.Fill(t);
        ++tot.nTrk;
        tot.nFake += t.isFake();
        tot.nNonGen += t.isNonGenuine();
        if (t.pt > kPtCut && std::abs(t.eta) < kEtaAcc) {
          ++tot.nTrkHi;
          tot.nFakeHi += t.isFake();
          tot.nNonGenHi += t.isNonGenuine();
        }
      }
    }
    R.Close();
    printf("[info] %s done\n", gSystem->BaseName(path));
    fflush(stdout);
  }
}

void PrintRate(const char* label, double num, double den) {
  double r = den > 0 ? num / den : 0;
  double err = den > 0 ? std::sqrt(r * (1 - r) / den) : 0;
  printf("  %-22s %10.0f / %10.0f = %7.4f +/- %.4f\n", label, num, den, r, err);
}

}  // namespace

void compare_fake(int nfiles = 50, const char* outname = "fake_qed_mumu.root",
                  const char* dirDef = "", const char* dirDum = "") {
  if (strlen(dirDef)) kDirDefault = dirDef;
  if (strlen(dirDum)) kDirDummy = dirDum;
  if (!strlen(kDirDefault) || !strlen(kDirDummy)) {
    printf("give the Default-stub and Dummy-stub ntuple locations (a directory of *.root, or one file)\n");
    return;
  }
  gROOT->SetBatch(true);

  const double ptbins[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
                           2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
  const int npt = sizeof(ptbins) / sizeof(double) - 1;

  Hists hDef, hDum;
  hDef.Book("def", npt, ptbins);
  hDum.Book("dum", npt, ptbins);
  Totals tDef, tDum;

  printf("=== Default Stub ===\n");
  Run(kDirDefault, nfiles, hDef, tDef);
  printf("=== Dummy Stub ===\n");
  Run(kDirDummy, nfiles, hDum, tDum);

  printf("\n================ L1 track fake rate ================\n");
  printf("fake        = N(trk_fake == 0)    / N(tracks)   -- no associated TP\n");
  printf("non-genuine = N(trk_genuine == 0) / N(tracks)   -- looser\n");
  printf("default: %ld files, %ld events    dummy: %ld files, %ld events\n",
         tDef.nFiles, tDef.nEvents, tDum.nFiles, tDum.nEvents);

  printf("\n-- fake rate, all pT --\n");
  PrintRate("Default Stub", tDef.nFake, tDef.nTrk);
  PrintRate("Dummy Stub", tDum.nFake, tDum.nTrk);

  printf("\n-- fake rate, pT > %.1f GeV, |eta| < %.1f --\n", kPtCut, kEtaAcc);
  PrintRate("Default Stub", tDef.nFakeHi, tDef.nTrkHi);
  PrintRate("Dummy Stub", tDum.nFakeHi, tDum.nTrkHi);

  printf("\n-- non-genuine rate, pT > %.1f GeV, |eta| < %.1f --\n", kPtCut, kEtaAcc);
  PrintRate("Default Stub", tDef.nNonGenHi, tDef.nTrkHi);
  PrintRate("Dummy Stub", tDum.nNonGenHi, tDum.nTrkHi);

  printf("\n-- tracks reconstructed --\n");
  printf("  default: %ld  (%.2f per event)\n", tDef.nTrk,
         tDef.nEvents ? (double)tDef.nTrk / tDef.nEvents : 0);
  printf("  dummy  : %ld  (%.2f per event)\n", tDum.nTrk,
         tDum.nEvents ? (double)tDum.nTrk / tDum.nEvents : 0);
  printf("====================================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (auto* hs : {&hDef, &hDum})
    for (auto* h : hs->All()) h->Write();
  fout->Close();
  printf("wrote %s\n", outname);
}
