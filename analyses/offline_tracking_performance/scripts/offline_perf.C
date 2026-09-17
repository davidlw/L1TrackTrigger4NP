// ---------------------------------------------------------------------------
// offline_perf.C
//
// Offline tracking performance from the OfflineTrackNtupleMaker tree: efficiency,
// fake rate and duplicate rate, differentially in pT, eta and phi.
//
//   efficiency     = N(TP matched by >=1 selected track) / N(tp)
//   fake rate      = N(trk_isTrue == 0)                  / N(selected trk)
//   duplicate rate = N(TP matched by  >1 selected track) / N(matched TP)
//
// All three are needed together: a configuration can buy efficiency by producing
// more tracks, which then shows up as fakes or duplicates instead.
//
// EVERYTHING here is inside |eta| < 2.4, the pT spectra included.
//
// pT runs from 0.3 GeV (the ptMinTP of the samples) to 10, finely binned at the
// low end where the turn-on is, and is meant to be drawn on a log x axis.
//
// eta, phi and Nch are filled for FOUR pT thresholds -- 0.3, 0.6, 1, 2 GeV -- so
// the four can be overlaid on one panel. Histogram names carry the threshold
// index: den_eta0_off is pT > 0.3, den_eta3_off is pT > 2.
//
// Nch is the event's TRUTH charged multiplicity (pT > 0.4, |eta| < 2.4), so the
// multiplicity axis is fixed by the event and does not shift when the
// reconstruction changes. Its range is the nchMax argument and is sample
// dependent -- roughly 10 for QED mumu, 200 for EPOS pPb, 2000+ for HYDJET PbPb.
// Entries above nchMax pile into the overflow and are dropped from the plot, so
// set it generously.
//
// hpOnly restricts to highPurity tracks. That selection applies to the efficiency
// numerator too, so the per-TP match count is rebuilt by joining tracks back to
// their TP rather than read from tp_nmatch -- that branch counts every match
// regardless of quality and would inflate both efficiency and duplicates. The
// join key is the TP's own truth values, written to trk_matchtp_* and tp_* from
// the same TrackingParticle, so they are identical floats and the match is exact.
// Validated with no selection: identical to tp_nmatch on 31301 TPs, duplicates too.
//
// Usage:
//   root -l -b -q 'offline_perf.C("../output/OfflineTrackNtuple.root")'
//   root -l -b -q 'offline_perf.C("/eos/.../0000","../output/offperf_hp.root",100,false,true)'
//
//   input   single .root file when nfiles = 0, otherwise the DIRECTORY holding
//           files named by `pattern` (1 .. nfiles)
//   nchMax  upper edge of the Nch axis; see above, it must suit the sample
//   tpClass 0 = all truth particles, 1 = PRIMARY only (tp_ngenpart > 0, i.e. from
//           the generator), 2 = SECONDARY only (tp_ngenpart == 0, made by GEANT
//           during detector simulation).
//
//           This matters a lot in the heavy-ion samples: roughly a quarter of the
//           stored TPs are GEANT secondaries and they are reconstructed far less
//           often than primaries, so the all-particle efficiency sits well below
//           the primary one. Run tpClass = 1 and 2 into separate output files to
//           see the two populations side by side.
//
//           The fake rate follows the same selection: with tpClass = 1 a track
//           matching only a secondary counts as fake, which is the usual
//           convention when the efficiency denominator is primaries.
// ---------------------------------------------------------------------------

#include <vector>
#include <map>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

const char* kTreePath = "OfflineTrackNtupleMaker/eventTree";

const double kEtaAcc = 2.4;

// pT thresholds for the eta/phi/Nch panels, overlaid by the plotting macro.
const int kNThr = 4;
const double kThr[kNThr] = {0.3, 0.6, 1.0, 2.0};

// Nch = PRIMARY truth charged particles with pT > 0.4, |eta| < 2.4, per event.
// Primary means tp_ngenpart > 0, i.e. from the event generator; GEANT secondaries
// are excluded, so the multiplicity axis reflects the collision rather than the
// detector material. Counting truth rather than tracks also means the axis does
// not move when the reconstruction changes.
const double kNchPtMin = 0.4;
const int kNchBins = 40;

struct TP {
  float pt, eta, phi, z0;
  int pdgid, nmatch, ngenpart;
  // tp_ngenpart is the number of associated GenParticles: 0 means the particle was
  // made by GEANT during detector simulation (a secondary), >0 means it came from
  // the event generator (a primary).
  int cls() const { return ngenpart > 0 ? 1 : 2; }   // 1 = primary, 2 = secondary
};
struct Trk {
  float pt, eta, phi, z0;
  int isTrue, highPurity;
  float mtp_pt, mtp_eta, mtp_phi;   // truth values of the TP this track matched
  bool isFake() const { return isTrue == 0; }
};

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *tp_pt = nullptr, *tp_eta = nullptr, *tp_phi = nullptr, *tp_z0 = nullptr;
  std::vector<int> *tp_pdgid = nullptr, *tp_nmatch = nullptr, *tp_ngen = nullptr;
  std::vector<float> *t_pt = nullptr, *t_eta = nullptr, *t_phi = nullptr, *t_z0 = nullptr;
  std::vector<int> *t_isTrue = nullptr, *t_hp = nullptr;
  std::vector<float> *t_mpt = nullptr, *t_meta = nullptr, *t_mphi = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree || !tree->GetBranch("tp_nmatch")) return false;
    tree->SetBranchStatus("*", 0);
    for (auto b : {"tp_pt", "tp_eta", "tp_phi", "tp_z0", "tp_pdgid", "tp_nmatch", "tp_ngenpart",
                   "trk_pt", "trk_eta", "trk_phi", "trk_z0", "trk_isTrue",
                   "trk_highPurity", "trk_matchtp_pt", "trk_matchtp_eta", "trk_matchtp_phi"})
      tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("tp_pt", &tp_pt);
    tree->SetBranchAddress("tp_eta", &tp_eta);
    tree->SetBranchAddress("tp_phi", &tp_phi);
    tree->SetBranchAddress("tp_z0", &tp_z0);
    tree->SetBranchAddress("tp_pdgid", &tp_pdgid);
    tree->SetBranchAddress("tp_nmatch", &tp_nmatch);
    tree->SetBranchAddress("tp_ngenpart", &tp_ngen);
    tree->SetBranchAddress("trk_pt", &t_pt);
    tree->SetBranchAddress("trk_eta", &t_eta);
    tree->SetBranchAddress("trk_phi", &t_phi);
    tree->SetBranchAddress("trk_z0", &t_z0);
    tree->SetBranchAddress("trk_isTrue", &t_isTrue);
    tree->SetBranchAddress("trk_highPurity", &t_hp);
    tree->SetBranchAddress("trk_matchtp_pt", &t_mpt);
    tree->SetBranchAddress("trk_matchtp_eta", &t_meta);
    tree->SetBranchAddress("trk_matchtp_phi", &t_mphi);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr; tree = nullptr;
  }
  void FlattenTP(std::vector<TP>& out) const {
    out.clear();
    for (size_t i = 0; i < tp_pt->size(); ++i) {
      TP t;
      t.pt = tp_pt->at(i); t.eta = tp_eta->at(i); t.phi = tp_phi->at(i);
      t.z0 = tp_z0->at(i); t.pdgid = tp_pdgid->at(i); t.nmatch = tp_nmatch->at(i);
      t.ngenpart = tp_ngen->at(i);
      out.push_back(t);
    }
  }
  void FlattenTrk(std::vector<Trk>& out) const {
    out.clear();
    for (size_t i = 0; i < t_pt->size(); ++i) {
      Trk t;
      t.pt = t_pt->at(i); t.eta = t_eta->at(i); t.phi = t_phi->at(i);
      t.z0 = t_z0->at(i); t.isTrue = t_isTrue->at(i); t.highPurity = t_hp->at(i);
      t.mtp_pt = t_mpt->at(i); t.mtp_eta = t_meta->at(i); t.mtp_phi = t_mphi->at(i);
      out.push_back(t);
    }
  }
};

// den/num   = efficiency (per TP)      denM/numD = duplicate rate (per matched TP)
// fden/fnum = fake rate (per track)
struct Hists {
  TH1D *den_pt, *num_pt, *denM_pt, *numD_pt, *fden_pt, *fnum_pt;
  TH1D *den_z0, *num_z0, *fden_z0, *fnum_z0;
  TH1D *den_eta[kNThr], *num_eta[kNThr], *denM_eta[kNThr], *numD_eta[kNThr];
  TH1D *fden_eta[kNThr], *fnum_eta[kNThr];
  TH1D *den_phi[kNThr], *num_phi[kNThr], *denM_phi[kNThr], *numD_phi[kNThr];
  TH1D *fden_phi[kNThr], *fnum_phi[kNThr];
  TH1D *den_nch[kNThr], *num_nch[kNThr], *denM_nch[kNThr], *numD_nch[kNThr];
  TH1D *fden_nch[kNThr], *fnum_nch[kNThr];
  std::vector<TH1D*> all;

  TH1D* B(const char* w, const char* v, const char* tag, int nb, double lo, double hi) {
    auto* h = new TH1D(Form("%s_%s_%s", w, v, tag), "", nb, lo, hi);
    h->Sumw2(); all.push_back(h); return h;
  }
  TH1D* Bv(const char* w, const char* v, const char* tag, int nb, const double* e) {
    auto* h = new TH1D(Form("%s_%s_%s", w, v, tag), "", nb, e);
    h->Sumw2(); all.push_back(h); return h;
  }

  void Book(const char* tag, int npt, const double* ptbins, double nchMax) {
    const double P = TMath::Pi();
    // Nch is an integer. For a narrow range give it one bin per unit rather than
    // 40 bins of 0.1, which would leave most of them empty (QED mumu has Nch <= 2).
    const int nchBins = (nchMax <= kNchBins) ? (int)nchMax : kNchBins;
    den_pt  = Bv("den",  "pt", tag, npt, ptbins);
    num_pt  = Bv("num",  "pt", tag, npt, ptbins);
    denM_pt = Bv("denM", "pt", tag, npt, ptbins);
    numD_pt = Bv("numD", "pt", tag, npt, ptbins);
    fden_pt = Bv("fden", "pt", tag, npt, ptbins);
    fnum_pt = Bv("fnum", "pt", tag, npt, ptbins);

    den_z0  = B("den",  "z0", tag, 40, -20, 20);
    num_z0  = B("num",  "z0", tag, 40, -20, 20);
    fden_z0 = B("fden", "z0", tag, 40, -20, 20);
    fnum_z0 = B("fnum", "z0", tag, 40, -20, 20);

    for (int i = 0; i < kNThr; ++i) {
      den_eta[i]  = B("den",  Form("eta%d", i), tag, 48, -kEtaAcc, kEtaAcc);
      num_eta[i]  = B("num",  Form("eta%d", i), tag, 48, -kEtaAcc, kEtaAcc);
      denM_eta[i] = B("denM", Form("eta%d", i), tag, 48, -kEtaAcc, kEtaAcc);
      numD_eta[i] = B("numD", Form("eta%d", i), tag, 48, -kEtaAcc, kEtaAcc);
      fden_eta[i] = B("fden", Form("eta%d", i), tag, 48, -kEtaAcc, kEtaAcc);
      fnum_eta[i] = B("fnum", Form("eta%d", i), tag, 48, -kEtaAcc, kEtaAcc);
      den_phi[i]  = B("den",  Form("phi%d", i), tag, 32, -P, P);
      num_phi[i]  = B("num",  Form("phi%d", i), tag, 32, -P, P);
      denM_phi[i] = B("denM", Form("phi%d", i), tag, 32, -P, P);
      numD_phi[i] = B("numD", Form("phi%d", i), tag, 32, -P, P);
      fden_phi[i] = B("fden", Form("phi%d", i), tag, 32, -P, P);
      fnum_phi[i] = B("fnum", Form("phi%d", i), tag, 32, -P, P);
      den_nch[i]  = B("den",  Form("nch%d", i), tag, nchBins, 0, nchMax);
      num_nch[i]  = B("num",  Form("nch%d", i), tag, nchBins, 0, nchMax);
      denM_nch[i] = B("denM", Form("nch%d", i), tag, nchBins, 0, nchMax);
      numD_nch[i] = B("numD", Form("nch%d", i), tag, nchBins, 0, nchMax);
      fden_nch[i] = B("fden", Form("nch%d", i), tag, nchBins, 0, nchMax);
      fnum_nch[i] = B("fnum", Form("nch%d", i), tag, nchBins, 0, nchMax);
    }
  }

  // Everything is inside |eta| < kEtaAcc. Entries below the first pT bin edge
  // (0.3) land in underflow and are excluded automatically.
  void FillTP(const TP& t, int nch) {
    if (std::abs(t.eta) >= kEtaAcc) return;
    const bool m = t.nmatch > 0, dup = t.nmatch > 1;
    den_pt->Fill(t.pt);
    if (m) {
      num_pt->Fill(t.pt);
      denM_pt->Fill(t.pt);
      if (dup) numD_pt->Fill(t.pt);
    }
    if (t.pt > kThr[0]) {
      den_z0->Fill(t.z0);
      if (m) num_z0->Fill(t.z0);
    }
    for (int i = 0; i < kNThr; ++i) {
      if (t.pt <= kThr[i]) continue;
      den_eta[i]->Fill(t.eta);
      den_phi[i]->Fill(t.phi);
      if (m) {
        num_eta[i]->Fill(t.eta);   num_phi[i]->Fill(t.phi);
        denM_eta[i]->Fill(t.eta);  denM_phi[i]->Fill(t.phi);
        if (dup) { numD_eta[i]->Fill(t.eta); numD_phi[i]->Fill(t.phi); }
      }
      den_nch[i]->Fill(nch);
      if (m) {
        num_nch[i]->Fill(nch);
        denM_nch[i]->Fill(nch);
        if (dup) numD_nch[i]->Fill(nch);
      }
    }
  }

  void FillTrk(const Trk& t, int nch) {
    if (std::abs(t.eta) >= kEtaAcc) return;
    const bool F = t.isFake();
    fden_pt->Fill(t.pt);
    if (F) fnum_pt->Fill(t.pt);
    if (t.pt > kThr[0]) {
      fden_z0->Fill(t.z0);
      if (F) fnum_z0->Fill(t.z0);
    }
    for (int i = 0; i < kNThr; ++i) {
      if (t.pt <= kThr[i]) continue;
      fden_eta[i]->Fill(t.eta);
      fden_phi[i]->Fill(t.phi);
      if (F) { fnum_eta[i]->Fill(t.eta); fnum_phi[i]->Fill(t.phi); }
      fden_nch[i]->Fill(nch);
      if (F) fnum_nch[i]->Fill(nch);
    }
  }
};

struct Totals {
  long nFiles = 0, nEvents = 0;
  long nTP = 0, nMatched = 0, nDup = 0;
  long nTPhi = 0, nMatchedHi = 0, nDupHi = 0;
  long nTrk = 0, nFake = 0, nTrkHi = 0, nFakeHi = 0;
};

void RunOne(const char* path, int pdgSel, bool hpOnly, int tpClass,
            Hists& h, Totals& tot) {
  Reader R;
  if (!R.Open(path)) { printf("[warn] cannot read %s\n", path); R.Close(); return; }
  ++tot.nFiles;
  std::vector<TP> tps;
  std::vector<Trk> trks;
  for (Long64_t e = 0; e < R.tree->GetEntries(); ++e) {
    R.tree->GetEntry(e);
    ++tot.nEvents;
    R.FlattenTP(tps);
    R.FlattenTrk(trks);

    // Event multiplicity: PRIMARY truth particles only, independent of both the
    // reconstruction and of GEANT secondary production.
    int nch = 0;
    for (const auto& t : tps)
      if (t.cls() == 1 && t.pt > kNchPtMin && std::abs(t.eta) < kEtaAcc) ++nch;

    // Per-TP count of SELECTED matching tracks (see the header note on the join).
    std::map<float, std::vector<int>> byPt;   // pt first, to keep this near-linear
    for (size_t i = 0; i < tps.size(); ++i) byPt[tps[i].pt].push_back(i);
    std::vector<int> nsel(tps.size(), 0);
    std::vector<int> mIdx(trks.size(), -1);   // TP each track matched, -1 if none
    for (size_t k = 0; k < trks.size(); ++k) {
      const auto& q = trks[k];
      if (q.isTrue != 1) continue;
      if (hpOnly && q.highPurity != 1) continue;
      auto it = byPt.find(q.mtp_pt);
      if (it == byPt.end()) continue;
      for (int i : it->second)
        if (tps[i].eta == q.mtp_eta && tps[i].phi == q.mtp_phi) {
          ++nsel[i]; mIdx[k] = i; break;
        }
    }

    for (size_t i = 0; i < tps.size(); ++i) {
      TP t = tps[i];
      if (pdgSel && std::abs(t.pdgid) != pdgSel) continue;
      if (tpClass && t.cls() != tpClass) continue;
      t.nmatch = nsel[i];               // selection-aware, replaces tp_nmatch
      h.FillTP(t, nch);
      if (t.pt <= kThr[0] || std::abs(t.eta) >= kEtaAcc) continue;   // totals: same acceptance
      ++tot.nTP;
      tot.nMatched += (t.nmatch > 0);
      tot.nDup += (t.nmatch > 1);
      if (t.pt > 2.0) {
        ++tot.nTPhi;
        tot.nMatchedHi += (t.nmatch > 0);
        tot.nDupHi += (t.nmatch > 1);
      }
    }
    for (size_t k = 0; k < trks.size(); ++k) {
      Trk t = trks[k];
      if (hpOnly && t.highPurity != 1) continue;
      // A track is "fake" for this class if it is not matched to a TP of it.
      // With tpClass = 0 that is the usual "matched to nothing". With tpClass = 1
      // a track matching only a GEANT secondary also counts as fake, which is the
      // standard convention when the efficiency denominator is primaries.
      const int cls = (mIdx[k] >= 0) ? tps[mIdx[k]].cls() : 0;
      t.isTrue = tpClass ? (cls == tpClass) : (cls != 0);
      h.FillTrk(t, nch);
      if (t.pt <= kThr[0] || std::abs(t.eta) >= kEtaAcc) continue;
      ++tot.nTrk;
      tot.nFake += t.isFake();
      if (t.pt > 2.0) {
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
                  int nfiles = 0, bool muonsOnly = false, bool hpOnly = false,
                  double nchMax = 200, int tpClass = 0, int pdgSel = 0) {
  // pdgSel: keep only truth particles with |pdgId| == pdgSel (11 electrons,
  // 13 muons, 211 pions ...); 0 = all. muonsOnly is the old spelling of 13.
  if (muonsOnly && !pdgSel) pdgSel = 13;
  gROOT->SetBatch(true);

  // Fine below 1 GeV, where the turn-on is; drawn on a log x axis.
  const double ptbins[] = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9,
                           1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0};
  const int npt = sizeof(ptbins) / sizeof(double) - 1;

  Hists h;
  h.Book("off", npt, ptbins, nchMax);
  Totals tot;

  // one file, or a directory of *.root (first nfiles of them; 0 = all)
  for (const auto& path : ListRootFiles(input, nfiles))
    RunOne(path, pdgSel, hpOnly, tpClass, h, tot);

  printf("\n============== offline tracking performance ==============\n");
  printf("efficiency     = N(TP matched by >=1 selected track) / N(tp)\n");
  printf("fake rate      = N(trk_isTrue == 0) / N(selected trk)\n");
  printf("duplicate rate = N(TP matched by  >1 selected track) / N(matched TP)\n");
  printf("%ld files, %ld events%s%s\n", tot.nFiles, tot.nEvents,
         pdgSel ? Form("  (|pdgId| == %d only)", pdgSel) : "",
         hpOnly ? "  [highPurity tracks only]" : "  [all tracks]");
  printf("truth particles: %s\n", tpClass == 1 ? "PRIMARY only (tp_ngenpart > 0)"
                                 : tpClass == 2 ? "SECONDARY only (tp_ngenpart == 0, made by GEANT)"
                                                : "all");

  printf("\n-- pT > %.1f GeV, |eta| < %.1f --\n", kThr[0], kEtaAcc);
  PrintRate("efficiency", tot.nMatched, tot.nTP);
  PrintRate("fake rate", tot.nFake, tot.nTrk);
  PrintRate("duplicate rate", tot.nDup, tot.nMatched);

  printf("\n-- pT > 2.0 GeV, |eta| < %.1f --\n", kEtaAcc);
  PrintRate("efficiency", tot.nMatchedHi, tot.nTPhi);
  PrintRate("fake rate", tot.nFakeHi, tot.nTrkHi);
  PrintRate("duplicate rate", tot.nDupHi, tot.nMatchedHi);

  printf("\n-- yields (pT > %.1f, |eta| < %.1f) --\n", kThr[0], kEtaAcc);
  printf("  tracking particles: %ld  (%.2f per event)\n", tot.nTP,
         tot.nEvents ? (double)tot.nTP / tot.nEvents : 0);
  printf("  offline tracks    : %ld  (%.2f per event)\n", tot.nTrk,
         tot.nEvents ? (double)tot.nTrk / tot.nEvents : 0);
  printf("=========================================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (auto* x : h.all) x->Write();
  fout->Close();
  printf("wrote %s\n", outname);
}
