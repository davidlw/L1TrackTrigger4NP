// ---------------------------------------------------------------------------
// mbeff_offline.C
//
// Event (minimum-bias) trigger efficiency vs the pT threshold:
//
//   eff(pTmin) = N(events with >=1 object, |eta| < 2.4, pT > pTmin) / N(all events)
//
// The denominator is EVERY event, the same for all thresholds -- this is an event
// trigger efficiency, not a per-track quantity. Events with nothing above the
// threshold count in the denominator and not the numerator, which is the point.
//
// Two curves:
//   truth   TrackingParticles -- the ceiling a perfect tracker could reach
//   reco    offline highPurity tracks -- what offline reconstruction delivers
//
// The truth curve here is meaningfully better than the one in the older L1 macro
// (../../L1_tracking_performance/scripts/compare_mbeff.C). That one had to warn
// that its TPs were stub-preselected (TP_minNStub / TP_minNStubLayer >= 3), so it
// bracketed rather than pinned the generator-level answer. The offline ntuplizer
// applies no hit or stub requirement, and these samples were produced with
// ptMinTP = 0.3, so the truth curve is a genuine ceiling down to 0.3 GeV.
//
// Thresholds are the bin centres of the output histograms: 0.3 to 5.0 GeV in
// 0.1 steps, so TGraphAsymmErrors(num, den) puts each point at its own threshold.
//
// Usage:
//   root -l -b -q 'mbeff_offline.C("/eos/.../0000","../output/mbeff_hydjet.root",900)'
// ---------------------------------------------------------------------------

#include <vector>
#include <cmath>
#include <cstdio>

namespace {

const char* kTreePath = "OfflineTrackNtupleMaker/eventTree";
const double kEtaAcc = 2.4;

// bin centres = the scanned thresholds: 0.30, 0.40, ... 5.00
const int    kNBin = 48;
const double kLo = 0.25, kHi = 5.05;

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *tp_pt = nullptr, *tp_eta = nullptr;
  std::vector<float> *t_pt = nullptr, *t_eta = nullptr;
  std::vector<int> *t_hp = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree || !tree->GetBranch("tp_pt")) return false;
    tree->SetBranchStatus("*", 0);
    for (auto b : {"tp_pt", "tp_eta", "trk_pt", "trk_eta", "trk_highPurity"})
      tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("tp_pt", &tp_pt);
    tree->SetBranchAddress("tp_eta", &tp_eta);
    tree->SetBranchAddress("trk_pt", &t_pt);
    tree->SetBranchAddress("trk_eta", &t_eta);
    tree->SetBranchAddress("trk_highPurity", &t_hp);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr; tree = nullptr;
  }
};

}  // namespace

void mbeff_offline(const char* input, const char* outname = "../output/mbeff.root",
                   int nfiles = 0,
                   const char* pattern = "OfflineTrackNtuple_%d.root") {
  gROOT->SetBatch(true);

  auto* den    = new TH1D("den_mb",  "", kNBin, kLo, kHi);
  auto* numT   = new TH1D("numT_mb", "", kNBin, kLo, kHi);   // truth
  auto* numR   = new TH1D("numR_mb", "", kNBin, kLo, kHi);   // offline highPurity
  for (auto* h : {den, numT, numR}) h->Sumw2();

  long nev = 0, nfile = 0;
  std::vector<TString> paths;
  if (nfiles <= 0) paths.push_back(input);
  else for (int i = 1; i <= nfiles; ++i) {
    TString p = TString(input) + "/" + Form(pattern, i);
    if (!gSystem->AccessPathName(p)) paths.push_back(p);
  }

  for (const auto& p : paths) {
    Reader R;
    if (!R.Open(p)) { R.Close(); continue; }
    ++nfile;
    for (Long64_t e = 0; e < R.tree->GetEntries(); ++e) {
      R.tree->GetEntry(e);
      ++nev;
      // highest-pT object in acceptance; 0 if the event has none
      double maxT = 0, maxR = 0;
      for (size_t i = 0; i < R.tp_pt->size(); ++i)
        if (std::abs(R.tp_eta->at(i)) < kEtaAcc && R.tp_pt->at(i) > maxT)
          maxT = R.tp_pt->at(i);
      for (size_t i = 0; i < R.t_pt->size(); ++i)
        if (R.t_hp->at(i) == 1 && std::abs(R.t_eta->at(i)) < kEtaAcc && R.t_pt->at(i) > maxR)
          maxR = R.t_pt->at(i);
      // an event fires threshold t iff its highest-pT object is above t
      for (int b = 1; b <= kNBin; ++b) {
        const double t = den->GetBinCenter(b);
        den->Fill(t);
        if (maxT > t) numT->Fill(t);
        if (maxR > t) numR->Fill(t);
      }
    }
    R.Close();
    if (nfile % 50 == 0) { printf("[info] %ld files, %ld events\n", nfile, nev); fflush(stdout); }
  }

  printf("\n========== event trigger efficiency ==========\n");
  printf("eff(pTmin) = N(events with >=1 object, |eta| < %.1f, pT > pTmin) / N(events)\n", kEtaAcc);
  printf("%ld files, %ld events\n\n", nfile, nev);
  printf("%10s %14s %14s\n", "pTmin", "truth (TP)", "reco (highPur)");
  for (int b = 1; b <= kNBin; ++b) {
    double d = den->GetBinContent(b);
    if (d < 1) continue;
    double t = den->GetBinCenter(b);
    const int t10 = (int)std::lround(t * 10);
    if (t > 1.0 && t10 % 5) continue;   // every 0.1 up to 1 GeV, then every 0.5
    printf("%10.2f %14.4f %14.4f\n", t, numT->GetBinContent(b) / d, numR->GetBinContent(b) / d);
  }
  printf("==============================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  den->Write(); numT->Write(); numR->Write();
  fout->Close();
  printf("wrote %s\n", outname);
}
