// ---------------------------------------------------------------------------
// compare_eff.C
//
// L1 tracking efficiency, Default Stub vs Dummy Stub.
//
//   efficiency = N(tracking particles with tp_nmatch > 0) / N(tracking particles)
//
// computed DIFFERENTIALLY -- vs pT, and vs eta / phi / z0 above the pT threshold.
// Each production is measured against its own stored tracking particles, and the
// two are read completely independently: no event pairing, no cross-matching, and
// the two samples need not have the same files or even the same events.
//
// Why that is safe. The ntuple only stores a TP with >= 3 stubs in >= 3 layers,
// which is a stub-based cut, so in principle the two productions could have
// different denominators. Measured bin by bin, they do not where it matters:
// above 2 GeV the two denominators agree to better than 0.5%, and the
// default-minus-dummy difference is the same to within ~0.005 in every bin. The
// denominators diverge only below ~1.6 GeV (up to 8x), where both efficiencies
// are identically zero in the current samples.
//
// CAUTION for future low-pT samples: once tracks are reconstructed below 2 GeV,
// those low-pT bins carry real efficiency and the denominator difference there
// becomes a genuine bias. Setting TP_minNStub = 0 and TP_minNStubLayer = 0 in the
// ntuplizer removes the stub dependence from the TP collection and fixes it at
// source.
//
// Usage:
//   root -l -b -q 'compare_eff.C(50)'                          // QED mumu
//   root -l -b -q 'compare_eff.C(100,false,"eff_hydjet.root",\
//                  "<defdir>","<dumdir>")'                     // hadronic sample
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

// Set on the command line (arguments 4 and 5); a directory of *.root or one file.
const char* kDirDefault = "";
const char* kDirDummy = "";
const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

// The track finder emits nothing below ~1.95 GeV, so the eta / phi / z0
// breakdowns apply this cut; the inclusive versions only show the turn-on.
const double kPtCut = 2.0;

// Second, looser threshold for the eta / phi breakdowns. The next dummy-tracking
// iteration reaches down to 1 GeV, so these are the plots that will show the gain
// over default. Restricted to the tracker acceptance |eta| < 2.4.
const double kPtCut1 = 1.0;
const double kEtaAcc = 2.4;

// One tracking particle, flattened out of the per-event vectors.
//
// NOTE: tp_d0 and tp_eventid are present as branches but were never filled in
// these productions (size 0 for every event), so they are not read here.
struct TP {
  float pt, eta, phi, z0;
  int pdgid, nmatch, nstub;
};

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *pt = nullptr, *eta = nullptr, *phi = nullptr, *z0 = nullptr;
  std::vector<int> *pdgid = nullptr, *nmatch = nullptr, *nstub = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree) return false;
    // Enable only what is needed -- the tp_* branches are a small fraction of
    // these files, so this is the difference between seconds and minutes.
    tree->SetBranchStatus("*", 0);
    for (auto b : {"tp_pt", "tp_eta", "tp_phi", "tp_z0", "tp_pdgid", "tp_nmatch", "tp_nstub"})
      tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("tp_pt", &pt);
    tree->SetBranchAddress("tp_eta", &eta);
    tree->SetBranchAddress("tp_phi", &phi);
    tree->SetBranchAddress("tp_z0", &z0);
    tree->SetBranchAddress("tp_pdgid", &pdgid);
    tree->SetBranchAddress("tp_nmatch", &nmatch);
    tree->SetBranchAddress("tp_nstub", &nstub);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr;
    tree = nullptr;
  }
  void Flatten(std::vector<TP>& out) const {
    out.clear();
    for (size_t i = 0; i < pt->size(); ++i) {
      TP t;
      t.pt = pt->at(i);      t.eta = eta->at(i);       t.phi = phi->at(i);
      t.z0 = z0->at(i);      t.pdgid = pdgid->at(i);
      t.nmatch = nmatch->at(i);  t.nstub = nstub->at(i);
      out.push_back(t);
    }
  }
};

// Numerator and denominator for one production.
struct Hists {
  TH1D *den_pt, *num_pt, *den_eta, *num_eta, *den_phi, *num_phi, *den_z0, *num_z0;
  TH1D *den_eta2, *num_eta2, *den_phi2, *num_phi2, *den_z02, *num_z02;
  TH1D *den_eta1, *num_eta1, *den_phi1, *num_phi1;
  // Duplicate rate: denM = TPs that were matched at all (nmatch > 0),
  // numD = those matched more than once (nmatch > 1).
  TH1D *denM_pt, *numD_pt, *denM_eta1, *numD_eta1, *denM_phi1, *numD_phi1;
  TH1D *denM_eta2, *numD_eta2, *denM_phi2, *numD_phi2;
  TH1D *den_phiF, *num_phiF;

  void Book(const char* tag, int npt, const double* ptbins) {
    den_pt  = new TH1D(Form("den_pt_%s", tag),  ";p_{T} [GeV];tracking particles", npt, ptbins);
    num_pt  = new TH1D(Form("num_pt_%s", tag),  ";p_{T} [GeV];matched", npt, ptbins);
    den_eta = new TH1D(Form("den_eta_%s", tag), ";#eta;tracking particles", 50, -2.5, 2.5);
    num_eta = new TH1D(Form("num_eta_%s", tag), ";#eta;matched", 50, -2.5, 2.5);
    den_phi = new TH1D(Form("den_phi_%s", tag), ";#phi;tracking particles", 32, -TMath::Pi(), TMath::Pi());
    num_phi = new TH1D(Form("num_phi_%s", tag), ";#phi;matched", 32, -TMath::Pi(), TMath::Pi());
    den_z0  = new TH1D(Form("den_z0_%s", tag),  ";z_{0} [cm];tracking particles", 40, -20, 20);
    num_z0  = new TH1D(Form("num_z0_%s", tag),  ";z_{0} [cm];matched", 40, -20, 20);
    // pT > kPtCut variants, within |eta| < kEtaAcc -- same acceptance as the
    // kPtCut1 set below, so the two thresholds are directly comparable.
    den_eta2 = new TH1D(Form("den_eta2_%s", tag), ";#eta;tracking particles", 48, -kEtaAcc, kEtaAcc);
    num_eta2 = new TH1D(Form("num_eta2_%s", tag), ";#eta;matched", 48, -kEtaAcc, kEtaAcc);
    den_phi2 = new TH1D(Form("den_phi2_%s", tag), ";#phi;tracking particles", 32, -TMath::Pi(), TMath::Pi());
    num_phi2 = new TH1D(Form("num_phi2_%s", tag), ";#phi;matched", 32, -TMath::Pi(), TMath::Pi());
    den_z02  = new TH1D(Form("den_z02_%s", tag),  ";z_{0} [cm];tracking particles", 40, -20, 20);
    num_z02  = new TH1D(Form("num_z02_%s", tag),  ";z_{0} [cm];matched", 40, -20, 20);
    // pT > kPtCut1 variants, within |eta| < kEtaAcc.
    den_eta1 = new TH1D(Form("den_eta1_%s", tag), ";#eta;tracking particles", 48, -kEtaAcc, kEtaAcc);
    num_eta1 = new TH1D(Form("num_eta1_%s", tag), ";#eta;matched", 48, -kEtaAcc, kEtaAcc);
    den_phi1 = new TH1D(Form("den_phi1_%s", tag), ";#phi;tracking particles", 32, -TMath::Pi(), TMath::Pi());
    num_phi1 = new TH1D(Form("num_phi1_%s", tag), ";#phi;matched", 32, -TMath::Pi(), TMath::Pi());
    denM_pt   = new TH1D(Form("denM_pt_%s", tag),   ";p_{T} [GeV];matched TPs", npt, ptbins);
    numD_pt   = new TH1D(Form("numD_pt_%s", tag),   ";p_{T} [GeV];duplicated", npt, ptbins);
    denM_eta1 = new TH1D(Form("denM_eta1_%s", tag), ";#eta;matched TPs", 48, -kEtaAcc, kEtaAcc);
    numD_eta1 = new TH1D(Form("numD_eta1_%s", tag), ";#eta;duplicated", 48, -kEtaAcc, kEtaAcc);
    denM_phi1 = new TH1D(Form("denM_phi1_%s", tag), ";#phi;matched TPs", 32, -TMath::Pi(), TMath::Pi());
    numD_phi1 = new TH1D(Form("numD_phi1_%s", tag), ";#phi;duplicated", 32, -TMath::Pi(), TMath::Pi());
    denM_eta2 = new TH1D(Form("denM_eta2_%s", tag), ";#eta;matched TPs", 48, -kEtaAcc, kEtaAcc);
    numD_eta2 = new TH1D(Form("numD_eta2_%s", tag), ";#eta;duplicated", 48, -kEtaAcc, kEtaAcc);
    denM_phi2 = new TH1D(Form("denM_phi2_%s", tag), ";#phi;matched TPs", 32, -TMath::Pi(), TMath::Pi());
    numD_phi2 = new TH1D(Form("numD_phi2_%s", tag), ";#phi;duplicated", 32, -TMath::Pi(), TMath::Pi());
    // Fine phi binning for the periodicity fold scan. 144 divides by 8, 9, 12,
    // 16, 18 and 24, so it folds cleanly onto any candidate symmetry.
    den_phiF = new TH1D(Form("den_phiF_%s", tag), ";#phi;tracking particles", 144, -TMath::Pi(), TMath::Pi());
    num_phiF = new TH1D(Form("num_phiF_%s", tag), ";#phi;matched", 144, -TMath::Pi(), TMath::Pi());
    for (auto h : All()) h->Sumw2();
  }

  std::vector<TH1D*> All() {
    return {den_pt,   num_pt,   den_eta,  num_eta,  den_phi,  num_phi,  den_z0,
            num_z0,   den_eta2, num_eta2, den_phi2, num_phi2, den_z02,  num_z02,
            den_eta1, num_eta1, den_phi1, num_phi1,
            denM_pt,   numD_pt,   denM_eta1, numD_eta1, denM_phi1, numD_phi1,
            denM_eta2, numD_eta2, denM_phi2, numD_phi2,
            den_phiF, num_phiF};
  }

  void Fill(const TP& t, bool matched) {
    const bool dup = t.nmatch > 1;   // matched by more than one track
    if (matched) {
      denM_pt->Fill(t.pt);
      if (dup) numD_pt->Fill(t.pt);
      if (t.pt > kPtCut1) {
        denM_eta1->Fill(t.eta);
        if (dup) numD_eta1->Fill(t.eta);
        if (std::abs(t.eta) < kEtaAcc) {
          denM_phi1->Fill(t.phi);
          if (dup) numD_phi1->Fill(t.phi);
        }
      }
      if (t.pt > kPtCut) {
        denM_eta2->Fill(t.eta);
        if (dup) numD_eta2->Fill(t.eta);
        if (std::abs(t.eta) < kEtaAcc) {
          denM_phi2->Fill(t.phi);
          if (dup) numD_phi2->Fill(t.phi);
        }
      }
    }
    den_pt->Fill(t.pt);
    den_eta->Fill(t.eta);
    den_phi->Fill(t.phi);
    den_z0->Fill(t.z0);
    if (matched) {
      num_pt->Fill(t.pt);
      num_eta->Fill(t.eta);
      num_phi->Fill(t.phi);
      num_z0->Fill(t.z0);
    }
    if (t.pt > kPtCut1) {
      den_eta1->Fill(t.eta);
      if (matched) num_eta1->Fill(t.eta);
      if (std::abs(t.eta) < kEtaAcc) {
        den_phi1->Fill(t.phi);
        if (matched) num_phi1->Fill(t.phi);
      }
    }
    if (t.pt > kPtCut) {
      den_eta2->Fill(t.eta);
      if (matched) num_eta2->Fill(t.eta);
      if (std::abs(t.eta) < kEtaAcc) {
        den_phi2->Fill(t.phi);
        den_z02->Fill(t.z0);
        den_phiF->Fill(t.phi);
        if (matched) {
          num_phi2->Fill(t.phi);
          num_z02->Fill(t.z0);
          num_phiF->Fill(t.phi);
        }
      }
    }
  }
};

// Running totals for the printed summary.
struct Totals {
  long nFiles = 0, nEvents = 0, nTP = 0, nMatched = 0, nTPhi = 0, nMatchedHi = 0;
  long nDup = 0, nDupHi = 0;
  double sumStub = 0;
};

// Read one production. Missing file indices are skipped, so the indices need not
// be contiguous and the two samples need not have the same ones.
void Run(const char* dir, int nfiles, bool muonsOnly, Hists& h, Totals& tot) {
  std::vector<TP> tps;
  for (const auto& path : ListRootFiles(dir, nfiles)) {
    Reader R;
    if (!R.Open(path)) { R.Close(); continue; }
    ++tot.nFiles;
    for (Long64_t e = 0; e < R.tree->GetEntries(); ++e) {
      R.tree->GetEntry(e);
      R.Flatten(tps);
      ++tot.nEvents;
      for (const auto& t : tps) {
        if (muonsOnly && std::abs(t.pdgid) != 13) continue;
        bool m = t.nmatch > 0;
        h.Fill(t, m);
        ++tot.nTP;
        tot.nMatched += m;
        tot.nDup += (t.nmatch > 1);
        tot.sumStub += t.nstub;
        if (t.pt > kPtCut) {
          ++tot.nTPhi;
          tot.nMatchedHi += m;
          tot.nDupHi += (t.nmatch > 1);
        }
      }
    }
    R.Close();
    printf("[info] %s done\n", gSystem->BaseName(path));
    fflush(stdout);
  }
}

void PrintEff(const char* label, double num, double den) {
  double e = den > 0 ? num / den : 0;
  double err = den > 0 ? std::sqrt(e * (1 - e) / den) : 0;
  printf("  %-22s %10.0f / %10.0f = %7.4f +/- %.4f\n", label, num, den, e, err);
}

}  // namespace

// muonsOnly must be false for hadronic samples (HYDJet, EPOS, pythia), where the
// tracking particles are overwhelmingly pions / kaons / protons.
void compare_eff(int nfiles = 50, bool muonsOnly = true,
                 const char* outname = "eff_qed_mumu.root", const char* dirDef = "",
                 const char* dirDum = "") {
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
  Run(kDirDefault, nfiles, muonsOnly, hDef, tDef);
  printf("=== Dummy Stub ===\n");
  Run(kDirDummy, nfiles, muonsOnly, hDum, tDum);

  printf("\n================ L1 tracking efficiency ================\n");
  printf("N(tp_nmatch > 0) / N(tracking particles), each sample on its own TPs\n");
  printf("selection: %s\n", muonsOnly ? "|pdgid| == 13" : "all pdgid");
  printf("default: %ld files, %ld events    dummy: %ld files, %ld events\n",
         tDef.nFiles, tDef.nEvents, tDum.nFiles, tDum.nEvents);

  printf("\n-- all pT --\n");
  PrintEff("Default Stub", tDef.nMatched, tDef.nTP);
  PrintEff("Dummy Stub", tDum.nMatched, tDum.nTP);
  printf("     (integrated over all pT this is dominated by the low-pT population,\n");
  printf("      which differs between the samples -- read the differential table instead)\n");

  printf("\n-- pT > %.1f GeV --\n", kPtCut);
  PrintEff("Default Stub", tDef.nMatchedHi, tDef.nTPhi);
  PrintEff("Dummy Stub", tDum.nMatchedHi, tDum.nTPhi);

  printf("\n-- duplicate rate (of matched TPs), pT > %.1f GeV --\n", kPtCut);
  PrintEff("Default Stub", tDef.nDupHi, tDef.nMatchedHi);
  PrintEff("Dummy Stub", tDum.nDupHi, tDum.nMatchedHi);

  printf("\n-- tracking particles stored --\n");
  printf("  default: %ld   (<tp_nstub> = %.2f)\n", tDef.nTP,
         tDef.nTP ? tDef.sumStub / tDef.nTP : 0);
  printf("  dummy  : %ld   (<tp_nstub> = %.2f)\n", tDum.nTP,
         tDum.nTP ? tDum.sumStub / tDum.nTP : 0);
  printf("========================================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (auto* hs : {&hDef, &hDum})
    for (auto* h : hs->All()) h->Write();
  fout->Close();
  printf("wrote %s\n", outname);
}
