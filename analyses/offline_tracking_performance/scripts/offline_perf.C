// ---------------------------------------------------------------------------
// offline_perf.C
//
// Offline tracking performance from the OfflineTrackNtupleMaker tree: efficiency,
// fake rate and duplicate rate, differentially in pT, eta, phi and z0.
//
//   efficiency     = N(tp_nmatch > 0)  / N(tp)              -- were real particles found
//   fake rate      = N(trk_isTrue == 0)/ N(trk)             -- are found tracks real
//   duplicate rate = N(tp_nmatch > 1)  / N(tp_nmatch > 0)   -- found more than once
//
// All three are needed together: a configuration can buy efficiency by producing
// more tracks, which shows up as fakes or duplicates rather than in the efficiency.
//
// Binning, acceptance (|eta| < 2.4) and the pT > 1 / pT > 2 variants are identical
// to the L1 macros in ../../L1_tracking_performance/scripts, so the offline curves
// overlay the L1 ones directly.
//
// Note the denominator: the TP selection in the ntuplizer carries NO hit or stub
// requirement, so efficiency here is not biased by the reconstruction. The low-pT
// reach is set upstream by ptMinTP in the step-2 mixing module.
//
// Usage:
//   root -l -b -q 'offline_perf.C("../output/OfflineTrackNtuple.root")'
//   root -l -b -q 'offline_perf.C("/path/to/dir","../output/offperf.root",50)'
//
//   input   single .root file when nfiles = 0, otherwise the DIRECTORY holding
//           files named by `pattern` (1 .. nfiles)
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>

namespace {

const char* kTreePath = "OfflineTrackNtupleMaker/eventTree";

const double kPtCut = 2.0;   // same as the L1 macros
const double kPtCut1 = 1.0;
const double kEtaAcc = 2.4;

struct TP {
  float pt, eta, phi, z0;
  int pdgid, nmatch;
};
struct Trk {
  float pt, eta, phi, z0;
  int isTrue;
  bool isFake() const { return isTrue == 0; }
};

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *tp_pt = nullptr, *tp_eta = nullptr, *tp_phi = nullptr, *tp_z0 = nullptr;
  std::vector<int> *tp_pdgid = nullptr, *tp_nmatch = nullptr;
  std::vector<float> *t_pt = nullptr, *t_eta = nullptr, *t_phi = nullptr, *t_z0 = nullptr;
  std::vector<int> *t_isTrue = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree || !tree->GetBranch("tp_nmatch")) return false;
    tree->SetBranchStatus("*", 0);
    for (auto b : {"tp_pt", "tp_eta", "tp_phi", "tp_z0", "tp_pdgid", "tp_nmatch",
                   "trk_pt", "trk_eta", "trk_phi", "trk_z0", "trk_isTrue"})
      tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("tp_pt", &tp_pt);
    tree->SetBranchAddress("tp_eta", &tp_eta);
    tree->SetBranchAddress("tp_phi", &tp_phi);
    tree->SetBranchAddress("tp_z0", &tp_z0);
    tree->SetBranchAddress("tp_pdgid", &tp_pdgid);
    tree->SetBranchAddress("tp_nmatch", &tp_nmatch);
    tree->SetBranchAddress("trk_pt", &t_pt);
    tree->SetBranchAddress("trk_eta", &t_eta);
    tree->SetBranchAddress("trk_phi", &t_phi);
    tree->SetBranchAddress("trk_z0", &t_z0);
    tree->SetBranchAddress("trk_isTrue", &t_isTrue);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr; tree = nullptr;
  }
  void FlattenTP(std::vector<TP>& out) const {
    out.clear();
    size_t n = tp_pt->size();
    for (size_t i = 0; i < n; ++i) {
      TP t;
      t.pt = tp_pt->at(i); t.eta = tp_eta->at(i); t.phi = tp_phi->at(i);
      t.z0 = tp_z0->at(i); t.pdgid = tp_pdgid->at(i); t.nmatch = tp_nmatch->at(i);
      out.push_back(t);
    }
  }
  void FlattenTrk(std::vector<Trk>& out) const {
    out.clear();
    size_t n = t_pt->size();
    for (size_t i = 0; i < n; ++i) {
      Trk t;
      t.pt = t_pt->at(i); t.eta = t_eta->at(i); t.phi = t_phi->at(i);
      t.z0 = t_z0->at(i); t.isTrue = t_isTrue->at(i);
      out.push_back(t);
    }
  }
};

// den/num   = efficiency (per TP)
// denM/numD = duplicate rate (per matched TP)
// fden/fnum = fake rate (per track)
struct Hists {
  TH1D *den_pt, *num_pt, *den_eta, *num_eta, *den_phi, *num_phi, *den_z0, *num_z0;
  TH1D *den_eta1, *num_eta1, *den_phi1, *num_phi1;
  TH1D *den_eta2, *num_eta2, *den_phi2, *num_phi2, *den_z02, *num_z02;
  TH1D *denM_pt, *numD_pt, *denM_eta1, *numD_eta1, *denM_phi1, *numD_phi1;
  TH1D *denM_eta2, *numD_eta2, *denM_phi2, *numD_phi2;
  TH1D *fden_pt, *fnum_pt, *fden_eta1, *fnum_eta1, *fden_phi1, *fnum_phi1;
  TH1D *fden_eta2, *fnum_eta2, *fden_phi2, *fnum_phi2, *fden_z02, *fnum_z02;

  TH1D* B(const char* w, const char* v, const char* tag, int nb, double lo, double hi) {
    return new TH1D(Form("%s_%s_%s", w, v, tag), "", nb, lo, hi);
  }
  TH1D* Bv(const char* w, const char* v, const char* tag, int nb, const double* e) {
    return new TH1D(Form("%s_%s_%s", w, v, tag), "", nb, e);
  }

  void Book(const char* tag, int npt, const double* ptbins) {
    const double P = TMath::Pi();
    den_pt = Bv("den", "pt", tag, npt, ptbins);  num_pt = Bv("num", "pt", tag, npt, ptbins);
    den_eta = B("den", "eta", tag, 48, -kEtaAcc, kEtaAcc);
    num_eta = B("num", "eta", tag, 48, -kEtaAcc, kEtaAcc);
    den_phi = B("den", "phi", tag, 32, -P, P);   num_phi = B("num", "phi", tag, 32, -P, P);
    den_z0 = B("den", "z0", tag, 40, -20, 20);   num_z0 = B("num", "z0", tag, 40, -20, 20);
    den_eta1 = B("den", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    num_eta1 = B("num", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    den_phi1 = B("den", "phi1", tag, 32, -P, P); num_phi1 = B("num", "phi1", tag, 32, -P, P);
    den_eta2 = B("den", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    num_eta2 = B("num", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    den_phi2 = B("den", "phi2", tag, 32, -P, P); num_phi2 = B("num", "phi2", tag, 32, -P, P);
    den_z02 = B("den", "z02", tag, 40, -20, 20); num_z02 = B("num", "z02", tag, 40, -20, 20);

    denM_pt = Bv("denM", "pt", tag, npt, ptbins); numD_pt = Bv("numD", "pt", tag, npt, ptbins);
    denM_eta1 = B("denM", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    numD_eta1 = B("numD", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    denM_phi1 = B("denM", "phi1", tag, 32, -P, P); numD_phi1 = B("numD", "phi1", tag, 32, -P, P);
    denM_eta2 = B("denM", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    numD_eta2 = B("numD", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    denM_phi2 = B("denM", "phi2", tag, 32, -P, P); numD_phi2 = B("numD", "phi2", tag, 32, -P, P);

    fden_pt = Bv("fden", "pt", tag, npt, ptbins); fnum_pt = Bv("fnum", "pt", tag, npt, ptbins);
    fden_eta1 = B("fden", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    fnum_eta1 = B("fnum", "eta1", tag, 48, -kEtaAcc, kEtaAcc);
    fden_phi1 = B("fden", "phi1", tag, 32, -P, P); fnum_phi1 = B("fnum", "phi1", tag, 32, -P, P);
    fden_eta2 = B("fden", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    fnum_eta2 = B("fnum", "eta2", tag, 48, -kEtaAcc, kEtaAcc);
    fden_phi2 = B("fden", "phi2", tag, 32, -P, P); fnum_phi2 = B("fnum", "phi2", tag, 32, -P, P);
    fden_z02 = B("fden", "z02", tag, 40, -20, 20); fnum_z02 = B("fnum", "z02", tag, 40, -20, 20);

    for (auto h : All()) h->Sumw2();
  }

  std::vector<TH1D*> All() {
    return {den_pt, num_pt, den_eta, num_eta, den_phi, num_phi, den_z0, num_z0,
            den_eta1, num_eta1, den_phi1, num_phi1,
            den_eta2, num_eta2, den_phi2, num_phi2, den_z02, num_z02,
            denM_pt, numD_pt, denM_eta1, numD_eta1, denM_phi1, numD_phi1,
            denM_eta2, numD_eta2, denM_phi2, numD_phi2,
            fden_pt, fnum_pt, fden_eta1, fnum_eta1, fden_phi1, fnum_phi1,
            fden_eta2, fnum_eta2, fden_phi2, fnum_phi2, fden_z02, fnum_z02};
  }

  void FillTP(const TP& t) {
    const bool m = t.nmatch > 0;
    const bool dup = t.nmatch > 1;
    den_pt->Fill(t.pt);  den_eta->Fill(t.eta);  den_phi->Fill(t.phi);  den_z0->Fill(t.z0);
    if (m) { num_pt->Fill(t.pt); num_eta->Fill(t.eta); num_phi->Fill(t.phi); num_z0->Fill(t.z0); }
    if (m) {
      denM_pt->Fill(t.pt);
      if (dup) numD_pt->Fill(t.pt);
    }
    if (t.pt > kPtCut1) {
      den_eta1->Fill(t.eta);
      if (m) num_eta1->Fill(t.eta);
      if (m) { denM_eta1->Fill(t.eta); if (dup) numD_eta1->Fill(t.eta); }
      if (std::abs(t.eta) < kEtaAcc) {
        den_phi1->Fill(t.phi);
        if (m) num_phi1->Fill(t.phi);
        if (m) { denM_phi1->Fill(t.phi); if (dup) numD_phi1->Fill(t.phi); }
      }
    }
    if (t.pt > kPtCut) {
      den_eta2->Fill(t.eta);
      if (m) num_eta2->Fill(t.eta);
      if (m) { denM_eta2->Fill(t.eta); if (dup) numD_eta2->Fill(t.eta); }
      if (std::abs(t.eta) < kEtaAcc) {
        den_phi2->Fill(t.phi);  den_z02->Fill(t.z0);
        if (m) { num_phi2->Fill(t.phi); num_z02->Fill(t.z0); }
        if (m) { denM_phi2->Fill(t.phi); if (dup) numD_phi2->Fill(t.phi); }
      }
    }
  }

  void FillTrk(const Trk& t) {
    const bool F = t.isFake();
    fden_pt->Fill(t.pt);
    if (F) fnum_pt->Fill(t.pt);
    if (t.pt > kPtCut1) {
      fden_eta1->Fill(t.eta);
      if (F) fnum_eta1->Fill(t.eta);
      if (std::abs(t.eta) < kEtaAcc) {
        fden_phi1->Fill(t.phi);
        if (F) fnum_phi1->Fill(t.phi);
      }
    }
    if (t.pt > kPtCut) {
      fden_eta2->Fill(t.eta);
      if (F) fnum_eta2->Fill(t.eta);
      if (std::abs(t.eta) < kEtaAcc) {
        fden_phi2->Fill(t.phi);  fden_z02->Fill(t.z0);
        if (F) { fnum_phi2->Fill(t.phi); fnum_z02->Fill(t.z0); }
      }
    }
  }
};

struct Totals {
  long nFiles = 0, nEvents = 0;
  long nTP = 0, nMatched = 0, nDup = 0;
  long nTPhi = 0, nMatchedHi = 0, nDupHi = 0;
  long nTrk = 0, nFake = 0, nTrkHi = 0, nFakeHi = 0;
};

void RunOne(const char* path, bool muonsOnly, Hists& h, Totals& tot) {
  Reader R;
  if (!R.Open(path)) { printf("[warn] cannot read %s\n", path); R.Close(); return; }
  ++tot.nFiles;
  std::vector<TP> tps;
  std::vector<Trk> trks;
  for (Long64_t e = 0; e < R.tree->GetEntries(); ++e) {
    R.tree->GetEntry(e);
    ++tot.nEvents;
    R.FlattenTP(tps);
    for (const auto& t : tps) {
      if (muonsOnly && std::abs(t.pdgid) != 13) continue;
      h.FillTP(t);
      ++tot.nTP;
      tot.nMatched += (t.nmatch > 0);
      tot.nDup += (t.nmatch > 1);
      if (t.pt > kPtCut && std::abs(t.eta) < kEtaAcc) {
        ++tot.nTPhi;
        tot.nMatchedHi += (t.nmatch > 0);
        tot.nDupHi += (t.nmatch > 1);
      }
    }
    R.FlattenTrk(trks);
    for (const auto& t : trks) {
      h.FillTrk(t);
      ++tot.nTrk;
      tot.nFake += t.isFake();
      if (t.pt > kPtCut && std::abs(t.eta) < kEtaAcc) {
        ++tot.nTrkHi;
        tot.nFakeHi += t.isFake();
      }
    }
  }
  R.Close();
  printf("[info] %s done\n", path);
  fflush(stdout);
}

void PrintRate(const char* label, double num, double den) {
  double r = den > 0 ? num / den : 0;
  double err = den > 0 ? std::sqrt(r * (1 - r) / den) : 0;
  printf("  %-30s %9.0f / %9.0f = %7.4f +/- %.4f\n", label, num, den, r, err);
}

}  // namespace

void offline_perf(const char* input = "../output/OfflineTrackNtuple.root",
                  const char* outname = "../output/offperf_qed_mumu.root",
                  int nfiles = 0, bool muonsOnly = false,
                  const char* pattern = "OfflineTrackNtuple_%d.root") {
  gROOT->SetBatch(true);

  const double ptbins[] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
                           2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
  const int npt = sizeof(ptbins) / sizeof(double) - 1;

  Hists h;
  h.Book("off", npt, ptbins);
  Totals tot;

  if (nfiles <= 0) {
    RunOne(input, muonsOnly, h, tot);            // input is one file
  } else {
    for (int i = 1; i <= nfiles; ++i) {          // input is a directory
      TString path = TString(input) + "/" + Form(pattern, i);
      if (gSystem->AccessPathName(path)) continue;
      RunOne(path, muonsOnly, h, tot);
    }
  }

  printf("\n============== offline tracking performance ==============\n");
  printf("efficiency     = N(tp_nmatch > 0)   / N(tp)\n");
  printf("fake rate      = N(trk_isTrue == 0) / N(trk)\n");
  printf("duplicate rate = N(tp_nmatch > 1)   / N(tp_nmatch > 0)\n");
  printf("%ld files, %ld events%s\n", tot.nFiles, tot.nEvents,
         muonsOnly ? "  (muons only)" : "");

  printf("\n-- all pT --\n");
  PrintRate("efficiency", tot.nMatched, tot.nTP);
  PrintRate("fake rate", tot.nFake, tot.nTrk);
  PrintRate("duplicate rate", tot.nDup, tot.nMatched);

  printf("\n-- pT > %.1f GeV, |eta| < %.1f --\n", kPtCut, kEtaAcc);
  PrintRate("efficiency", tot.nMatchedHi, tot.nTPhi);
  PrintRate("fake rate", tot.nFakeHi, tot.nTrkHi);
  PrintRate("duplicate rate", tot.nDupHi, tot.nMatchedHi);

  printf("\n-- yields --\n");
  printf("  tracking particles: %ld  (%.2f per event)\n", tot.nTP,
         tot.nEvents ? (double)tot.nTP / tot.nEvents : 0);
  printf("  offline tracks    : %ld  (%.2f per event)\n", tot.nTrk,
         tot.nEvents ? (double)tot.nTrk / tot.nEvents : 0);
  printf("=========================================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (auto* x : h.All()) x->Write();
  fout->Close();
  printf("wrote %s\n", outname);
}
