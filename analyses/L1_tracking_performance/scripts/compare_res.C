// ---------------------------------------------------------------------------
// compare_res.C
//
// Momentum resolution vs truth pT, Default Stub vs Dummy Stub, STARlight QED mumu.
//
// Two residual definitions are filled:
//   absolute : pT(reco) - pT(truth)          [GeV]   -- the primary one
//   relative : (pT(reco) - pT(truth)) / pT(truth)
//
// evaluated for tracking particles with tp_nmatch > 0. The matchtrk_* vectors are
// index-aligned with the tp_* vectors (verified: same total length, and every
// tp_nmatch>0 entry has matchtrk_pt > 0).
//
// The absolute residual is also stored as a 2D distribution vs truth pT, from
// which plot_res.C builds the profile (mean and RMS vs truth pT).
//
// Two populations are reported:
//   "all"    : every matched TP in each sample. Mixes a resolution difference
//              with a population difference.
//   "common" : TPs both samples stored AND both reconstructed, cross-matched on
//              generator kinematics. Identical particles on both sides, so any
//              difference is purely resolution. This is the fair comparison.
//
// The track finder emits nothing below ~1.95 GeV, so the eta and phi breakdowns
// apply pT > 2 GeV.
//
// Usage: root -l -b -q 'compare_res.C(50,"res_qed_mumu.root")'
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

// Set on the command line (arguments 3 and 4); a directory of *.root or one file.
// The two productions must hold the SAME events, file for file.
const char* kDirDefault = "";
const char* kDirDummy = "";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

// Coarse bins for the quoted tables. Nothing is reconstructed below ~1.8 GeV.
const double kPtBins[] = {1.8, 2.0, 2.25, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
const int kNPt = sizeof(kPtBins) / sizeof(double) - 1;

// Finer binning for the 2D distribution and its profile.
const double kPt2D[] = {1.8, 1.9, 2.0,  2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                        3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 10.0};
const int kNPt2D = sizeof(kPt2D) / sizeof(double) - 1;

// Angular breakdowns, all with pT > 2 GeV.
const double kPtCut = 2.0;
const double kEtaBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
const int kNEta = sizeof(kEtaBins) / sizeof(double) - 1;
const int kNPhi = 12;

int FindBin(double x, const double* edges, int n) {
  for (int i = 0; i < n; ++i)
    if (x >= edges[i] && x < edges[i + 1]) return i;
  return -1;
}
int PhiBin(double phi) {
  int i = (int)std::floor((phi + TMath::Pi()) / (2 * TMath::Pi() / kNPhi));
  return (i >= 0 && i < kNPhi) ? i : -1;
}
double PhiLow(int i) { return -TMath::Pi() + i * (2 * TMath::Pi() / kNPhi); }

struct TP {
  float pt, eta, phi;
  float mpt, meta, mz0;
  int pdgid, nmatch, nstub, mnstub;
  bool matched() const { return nmatch > 0 && mpt > 0; }
  float dabs() const { return mpt - pt; }         // GeV
  float drel() const { return (mpt - pt) / pt; }
};

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr;
  std::vector<float> *mpt = nullptr, *meta = nullptr, *mz0 = nullptr;
  std::vector<int> *pdgid = nullptr, *nmatch = nullptr, *nstub = nullptr, *mnstub = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree) return false;
    tree->SetBranchStatus("*", 0);
    const char* on[] = {"tp_pt",       "tp_eta",        "tp_phi",      "tp_pdgid",
                        "tp_nmatch",   "tp_nstub",      "matchtrk_pt", "matchtrk_eta",
                        "matchtrk_z0", "matchtrk_nstub"};
    for (auto b : on) tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("tp_pt", &pt);
    tree->SetBranchAddress("tp_eta", &eta);
    tree->SetBranchAddress("tp_phi", &phi);
    tree->SetBranchAddress("tp_pdgid", &pdgid);
    tree->SetBranchAddress("tp_nmatch", &nmatch);
    tree->SetBranchAddress("tp_nstub", &nstub);
    tree->SetBranchAddress("matchtrk_pt", &mpt);
    tree->SetBranchAddress("matchtrk_eta", &meta);
    tree->SetBranchAddress("matchtrk_z0", &mz0);
    tree->SetBranchAddress("matchtrk_nstub", &mnstub);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr;
    tree = nullptr;
  }
  void Flatten(std::vector<TP>& out) const {
    out.clear();
    size_t n = pt->size();
    if (mpt->size() != n) return;  // guard: some matchtrk_* branches are unfilled
    for (size_t i = 0; i < n; ++i) {
      TP t;
      t.pt = pt->at(i);
      t.eta = eta->at(i);
      t.phi = phi->at(i);
      t.pdgid = pdgid->at(i);
      t.nmatch = nmatch->at(i);
      t.nstub = nstub->at(i);
      t.mpt = mpt->at(i);
      t.meta = (meta->size() == n) ? meta->at(i) : -999;
      t.mz0 = (mz0->size() == n) ? mz0->at(i) : -999;
      t.mnstub = (mnstub->size() == n) ? mnstub->at(i) : -1;
      out.push_back(t);
    }
  }
};

float dphi(float a, float b) {
  float d = a - b;
  while (d > TMath::Pi()) d -= 2 * TMath::Pi();
  while (d < -TMath::Pi()) d += 2 * TMath::Pi();
  return d;
}
bool SameTP(const TP& a, const TP& b) {
  if (a.pdgid != b.pdgid) return false;
  if (std::abs(a.pt - b.pt) > 1e-3 * std::max(1.f, a.pt)) return false;
  if (std::abs(a.eta - b.eta) > 1e-3) return false;
  if (std::abs(dphi(a.phi, b.phi)) > 1e-3) return false;
  return true;
}

struct Res {
  double n = 0, mean = 0, meanErr = 0, rms = 0, gsig = 0, q68 = 0, tail = 0;
};

Res Extract(TH1D* h) {
  Res r;
  if (!h) return r;
  r.n = h->GetEntries();
  if (r.n < 1) return r;
  r.tail = (h->GetBinContent(0) + h->GetBinContent(h->GetNbinsX() + 1)) / r.n;
  r.mean = h->GetMean();
  r.meanErr = h->GetMeanError();
  r.rms = h->GetRMS();
  if (r.n >= 50) {
    double q[2] = {0.1587, 0.8413}, v[2];
    h->GetQuantiles(2, v, q);
    r.q68 = 0.5 * (v[1] - v[0]);
    TF1 f("f", "gaus");
    double lo = r.mean - 2 * r.rms, hi = r.mean + 2 * r.rms;
    for (int it = 0; it < 3; ++it) {
      if (h->Fit(&f, "QN0", "", lo, hi) != 0) break;
      double m = f.GetParameter(1), s = std::abs(f.GetParameter(2));
      if (s <= 0 || !std::isfinite(s)) break;
      r.gsig = s;
      lo = m - 2 * s;
      hi = m + 2 * s;
    }
  }
  return r;
}

const char* kSetName[2] = {"all", "common"};
const char* kSmpName[2] = {"def", "dum"};

}  // namespace

void compare_res(int nfiles = 50, const char* outname = "res_qed_mumu.root",
                 const char* dirDef = "", const char* dirDum = "") {
  if (strlen(dirDef)) kDirDefault = dirDef;
  if (strlen(dirDum)) kDirDummy = dirDum;
  if (!strlen(kDirDefault) || !strlen(kDirDummy)) {
    printf("give the Default-stub and Dummy-stub ntuple locations (a directory of *.root, or one file)\n");
    return;
  }
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  // Absolute residual (GeV) -- primary. Range covers the observed spread with room
  // for tails; overflow is counted and reported.
  TH1D* hAbsPt[2][2][16];
  TH1D* hAbsEta[2][2][8];
  TH1D* hAbsPhi[2][2][16];
  TH1D* hRelPt[2][2][16];
  TH1D* hRelEta[2][2][8];
  TH1D* hRelPhi[2][2][16];
  TH2D* h2D[2][2];

  for (int s = 0; s < 2; ++s)
    for (int m = 0; m < 2; ++m) {
      for (int i = 0; i < kNPt; ++i) {
        hAbsPt[s][m][i] = new TH1D(Form("abs_pt%d_%s_%s", i, kSmpName[m], kSetName[s]),
                                   ";p_{T}^{reco}-p_{T}^{truth} [GeV];tracks", 400, -1.0, 1.0);
        hRelPt[s][m][i] = new TH1D(Form("rel_pt%d_%s_%s", i, kSmpName[m], kSetName[s]),
                                   ";(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};tracks", 400, -1, 1);
      }
      for (int i = 0; i < kNEta; ++i) {
        hAbsEta[s][m][i] = new TH1D(Form("abs_eta%d_%s_%s", i, kSmpName[m], kSetName[s]),
                                    ";p_{T}^{reco}-p_{T}^{truth} [GeV];tracks", 400, -1.0, 1.0);
        hRelEta[s][m][i] = new TH1D(Form("rel_eta%d_%s_%s", i, kSmpName[m], kSetName[s]),
                                    ";(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};tracks", 400, -1, 1);
      }
      for (int i = 0; i < kNPhi; ++i) {
        hAbsPhi[s][m][i] = new TH1D(Form("abs_phi%d_%s_%s", i, kSmpName[m], kSetName[s]),
                                    ";p_{T}^{reco}-p_{T}^{truth} [GeV];tracks", 400, -1.0, 1.0);
        hRelPhi[s][m][i] = new TH1D(Form("rel_phi%d_%s_%s", i, kSmpName[m], kSetName[s]),
                                    ";(p_{T}^{reco}-p_{T}^{truth})/p_{T}^{truth};tracks", 400, -1, 1);
      }
      h2D[s][m] = new TH2D(Form("h2_dpt_vs_pt_%s_%s", kSmpName[m], kSetName[s]),
                           ";truth p_{T} [GeV];p_{T}^{reco}-p_{T}^{truth} [GeV]",
                           kNPt2D, kPt2D, 400, -1.0, 1.0);
    }

  long nEvents = 0, nFilesOK = 0, nMatchDef = 0, nMatchDum = 0, nCommon = 0;

  // same basename on both sides = same events; anything unpaired is dropped
  std::vector<TString> filesA = ListRootFiles(kDirDefault, nfiles), filesB = ListRootFiles(kDirDummy, nfiles);
  KeepCommonFiles(filesA, filesB);
  for (size_t i = 0; i < filesA.size(); ++i) {
    Reader A, B;
    bool okA = A.Open(filesA[i]);
    bool okB = B.Open(filesB[i]);
    if (!okA || !okB) {
      printf("[warn] skipping %s (def=%d dum=%d)\n", gSystem->BaseName(filesA[i]), okA, okB);
      A.Close(); B.Close();
      continue;
    }
    Long64_t nA = A.tree->GetEntries(), nB = B.tree->GetEntries();
    if (nA != nB) {
      printf("[warn] %s entry mismatch %lld vs %lld -- skipping\n", gSystem->BaseName(filesA[i]), nA, nB);
      A.Close(); B.Close();
      continue;
    }
    ++nFilesOK;

    std::vector<TP> a, b;
    auto fillOne = [&](const TP& t, int set, int smp) {
      int ip = FindBin(t.pt, kPtBins, kNPt);
      if (ip >= 0) {
        hAbsPt[set][smp][ip]->Fill(t.dabs());
        hRelPt[set][smp][ip]->Fill(t.drel());
      }
      h2D[set][smp]->Fill(t.pt, t.dabs());
      if (t.pt > kPtCut) {
        int ie = FindBin(std::abs(t.eta), kEtaBins, kNEta);
        if (ie >= 0) {
          hAbsEta[set][smp][ie]->Fill(t.dabs());
          hRelEta[set][smp][ie]->Fill(t.drel());
        }
        int ih = PhiBin(t.phi);
        if (ih >= 0) {
          hAbsPhi[set][smp][ih]->Fill(t.dabs());
          hRelPhi[set][smp][ih]->Fill(t.drel());
        }
      }
    };

    for (Long64_t e = 0; e < nA; ++e) {
      A.tree->GetEntry(e);
      B.tree->GetEntry(e);
      A.Flatten(a);
      B.Flatten(b);
      ++nEvents;

      for (const auto& t : a)
        if (t.matched()) { ++nMatchDef; fillOne(t, 0, 0); }
      for (const auto& t : b)
        if (t.matched()) { ++nMatchDum; fillOne(t, 0, 1); }

      std::vector<bool> usedB(b.size(), false);
      for (size_t ia = 0; ia < a.size(); ++ia) {
        if (!a[ia].matched()) continue;
        for (size_t ib = 0; ib < b.size(); ++ib) {
          if (usedB[ib] || !b[ib].matched()) continue;
          if (!SameTP(a[ia], b[ib])) continue;
          usedB[ib] = true;
          ++nCommon;
          fillOne(a[ia], 1, 0);
          fillOne(b[ib], 1, 1);
          break;
        }
      }
    }
    A.Close();
    B.Close();
    printf("[info] %s done (%lld events)\n", gSystem->BaseName(filesA[i]), nA);
    fflush(stdout);
  }

  printf("\n=========== QED mumu momentum resolution ===========\n");
  printf("files: %ld   events: %ld\n", nFilesOK, nEvents);
  printf("matched TPs: default=%ld  dummy=%ld  common (both matched)=%ld\n",
         nMatchDef, nMatchDum, nCommon);

  for (int s = 0; s < 2; ++s) {
    printf("\n############ population: %s ############\n",
           s == 0 ? "all matched TPs in each sample" : "TPs reconstructed in BOTH samples");

    printf("\n-- ABSOLUTE residual pT(reco)-pT(truth) [GeV] vs truth pT --\n");
    printf("%14s %8s | %9s %8s %8s %8s | %9s %8s %8s %8s | %8s\n", "pT bin [GeV]", "N",
           "mean(def)", "RMS", "gaus", "q68", "mean(dum)", "RMS", "gaus", "q68", "dRMS");
    for (int i = 0; i < kNPt; ++i) {
      Res rd = Extract(hAbsPt[s][0][i]), ru = Extract(hAbsPt[s][1][i]);
      if (rd.n < 1 && ru.n < 1) continue;
      printf("%6.2f-%7.2f %8.0f | %+9.4f %8.4f %8.4f %8.4f | %+9.4f %8.4f %8.4f %8.4f | %+7.1f%%\n",
             kPtBins[i], kPtBins[i + 1], rd.n, rd.mean, rd.rms, rd.gsig, rd.q68, ru.mean,
             ru.rms, ru.gsig, ru.q68, rd.rms > 0 ? 100 * (ru.rms - rd.rms) / rd.rms : 0);
    }

    printf("\n-- RELATIVE residual (pT(reco)-pT(truth))/pT(truth) vs truth pT --\n");
    printf("%14s %8s | %9s %8s %8s | %9s %8s %8s | %8s\n", "pT bin [GeV]", "N",
           "mean(def)", "RMS", "q68", "mean(dum)", "RMS", "q68", "dq68");
    for (int i = 0; i < kNPt; ++i) {
      Res rd = Extract(hRelPt[s][0][i]), ru = Extract(hRelPt[s][1][i]);
      if (rd.n < 1 && ru.n < 1) continue;
      printf("%6.2f-%7.2f %8.0f | %+9.4f %8.4f %8.4f | %+9.4f %8.4f %8.4f | %+7.1f%%\n",
             kPtBins[i], kPtBins[i + 1], rd.n, rd.mean, rd.rms, rd.q68, ru.mean, ru.rms,
             ru.q68, rd.q68 > 0 ? 100 * (ru.q68 - rd.q68) / rd.q68 : 0);
    }

    printf("\n-- vs |eta|, pT > 2 GeV (relative residual) --\n");
    printf("%14s %8s | %8s %8s | %8s %8s | %8s\n", "|eta| bin", "N", "RMS(def)", "q68(def)",
           "RMS(dum)", "q68(dum)", "dq68");
    for (int i = 0; i < kNEta; ++i) {
      Res rd = Extract(hRelEta[s][0][i]), ru = Extract(hRelEta[s][1][i]);
      if (rd.n < 1 && ru.n < 1) continue;
      printf("%6.2f-%7.2f %8.0f | %8.4f %8.4f | %8.4f %8.4f | %+7.1f%%\n", kEtaBins[i],
             kEtaBins[i + 1], rd.n, rd.rms, rd.q68, ru.rms, ru.q68,
             rd.q68 > 0 ? 100 * (ru.q68 - rd.q68) / rd.q68 : 0);
    }

    printf("\n-- vs phi, pT > 2 GeV (relative residual) --\n");
    printf("%14s %8s | %8s %8s | %8s %8s | %8s\n", "phi bin", "N", "RMS(def)", "q68(def)",
           "RMS(dum)", "q68(dum)", "dq68");
    for (int i = 0; i < kNPhi; ++i) {
      Res rd = Extract(hRelPhi[s][0][i]), ru = Extract(hRelPhi[s][1][i]);
      if (rd.n < 1 && ru.n < 1) continue;
      printf("%6.2f-%7.2f %8.0f | %8.4f %8.4f | %8.4f %8.4f | %+7.1f%%\n", PhiLow(i),
             PhiLow(i + 1), rd.n, rd.rms, rd.q68, ru.rms, ru.q68,
             rd.q68 > 0 ? 100 * (ru.q68 - rd.q68) / rd.q68 : 0);
    }

    printf("\n-- overflow fraction, |pT(reco)-pT(truth)| > 1 GeV --\n");
    for (int i = 0; i < kNPt; ++i) {
      Res rd = Extract(hAbsPt[s][0][i]), ru = Extract(hAbsPt[s][1][i]);
      if (rd.n < 1 && ru.n < 1) continue;
      printf("  %5.2f-%5.2f GeV   default %.5f   dummy %.5f\n", kPtBins[i], kPtBins[i + 1],
             rd.tail, ru.tail);
    }
  }
  printf("\n====================================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (int s = 0; s < 2; ++s)
    for (int m = 0; m < 2; ++m) {
      for (int i = 0; i < kNPt; ++i) { hAbsPt[s][m][i]->Write(); hRelPt[s][m][i]->Write(); }
      for (int i = 0; i < kNEta; ++i) { hAbsEta[s][m][i]->Write(); hRelEta[s][m][i]->Write(); }
      for (int i = 0; i < kNPhi; ++i) { hAbsPhi[s][m][i]->Write(); hRelPhi[s][m][i]->Write(); }
      h2D[s][m]->Write();
    }
  fout->Close();
  printf("wrote %s\n", outname);
}
