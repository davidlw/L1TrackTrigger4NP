// ---------------------------------------------------------------------------
// compare_stubs.C
//
// Stub-level comparison between the Default and Dummy stub productions, split
// into ACCEPTED and REJECTED stubs. Sample-agnostic: pass the two directories,
// a label for the plot title and a prefix for the file names.
//
// Both live in the same allstub_* collection and are separated by the flag
// allstub_isRejected (0 = accepted, non-zero = rejected).
//
// PRODUCTIONS DIFFER IN WHETHER DUMMY REJECTS ANYTHING -- check before reading a
// zero. Fingerprinted over nalewis's EOS area (a dummy production always has
// allstub_trigBend at the 999999 sentinel; the rejected fraction is what varies):
//
//   EPOS   dummy 260825_213152  rej 0.0000     dummy 260828_193234  rej 0.0062
//   HYDJET dummy 260825_213428  rej 0.0000     dummy 260828_193408  rej 0.0265
//   pp     dummy 260720_222437  rej 0.0335
//   QEDee  dummy 260825_212958  rej 0.0000     dummy 260825_163451  rej 0.0060
//
// So an empty rejected curve means "this production applied no stub selection",
// not "dummy stubs cannot be rejected". Prefer the newer productions unless you
// specifically want the no-selection variant.
//
// WHY THE 260825 DUMMY PRODUCTIONS HAVE ZERO REJECTED -- measured, not assumed.
//
// Measured on EPOS pPb, 25000 events, the same events in both productions:
//
//   default : 157.8 accepted + 1.04 rejected stubs/evt   (0.7% rejected)
//             clusters 1408.0/evt on sensor 0, 1368.4/evt on sensor 1
//   dummy   : 1408.0 accepted + ZERO rejected stubs/evt
//
// dummy stubs / default sensor-0 clusters = 1.0000. In THAT production every
// inner-sensor cluster became exactly one dummy stub: nothing was paired and no
// window was applied, so no selection existed for a stub to fail. The zero is
// real, not a missing collection or a broken branch -- allstub_isRejected is
// present, its length always matches allstub_x, and it is 0 for all 7.1M stubs.
// The 260828 productions of the same samples do reject stubs, so this is a
// statement about that production and not about dummy stubs in general.
//
// The same fact explains the other two signatures: the ~8.9x larger accepted
// yield (default needs a matched pair passing the window) and allstub_trigBend
// sitting at the 999999 sentinel (bend is undefined without a pair).
//
// Note the dummy ntuples also lack the eight ttclus_* branches the default ones
// carry, so the two were produced with different ntuplizer configurations.
//
// This macro only FILLS and SAVES. Drawing is plot_stubs.C, so axis ranges,
// binning and styling can be changed without re-reading the samples -- the same
// split the other compare_*/plot_* pairs in this directory use.
//
// The multiplicity is saved with ONE BIN PER STUB (no rebinning yet), so the
// plot step is free to choose any binning or upper limit afterwards.
//
// Distributions, all stored raw so the plot step can normalise as it likes
// rather than being normalised away:
//   multiplicity, r = sqrt(x^2+y^2), z, layer/disk index, trigBend
//
// Layer numbering: allstub_layer runs 1-6 for barrel AND 1-5 for disks, so it is
// only meaningful combined with allstub_isBarrel. Encoded here as a combined
// index: 0-5 = barrel L1-L6, 6-10 = disk D1-D5.
//
// trigBend is kept even though dummy has no usable value: an empty dummy curve
// on that panel is the point, showing the bend information really is absent.
//
// Usage:
//   root -l -b -q 'compare_stub_epos.C(2)'                 // files per sample
//   root -l -b -q 'compare_stub_epos.C(2,"<defDir>","<dumDir>","_epos")'
// ---------------------------------------------------------------------------

#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include "../../common/InputFiles.h"

namespace {

const char* kTreePath = "L1TrackHitNtupleMaker/eventTree";

// Set on the command line (arguments 2 and 3); a directory of *.root or one file.
// The EPOS pPb productions used in the README are nalewis 260825_201453 (default)
// and 260825_213152 (dummy) under /eos/cms/store/group/phys_heavyions/nalewis/.
const char* kDirDefault = "";
const char* kDirDummy = "";

void SaveBoth(TCanvas* c, const char* pdfPath) {
  c->SaveAs(pdfPath);
  TString png(pdfPath);
  png.ReplaceAll(".pdf", ".png");
  c->SaveAs(png);
}

struct Hists {
  TH1D *r, *z, *ld, *bend;
  std::vector<int> counts;      // per-event stub count, kept at full resolution
  long nev = 0;

  void Book(const char* tag) {
    r    = new TH1D(Form("r_%s", tag), "", 60, 0, 120);
    z    = new TH1D(Form("z_%s", tag), "", 140, -280, 280);
    ld   = new TH1D(Form("ld_%s", tag), "", 11, -0.5, 10.5);
    bend = new TH1D(Form("bend_%s", tag), "", 65, -8, 8);
    for (auto* h : {r, z, ld, bend}) h->Sumw2();
  }
  std::vector<TH1D*> All() { return {r, z, ld, bend}; }

  // one bin per stub; the plot step rebins from this
  TH1D* MakeFine(const char* tag, int nmax) {
    auto* h = new TH1D(Form("multfine_%s", tag), "", nmax, 0.5, nmax + 0.5);
    h->Sumw2();
    // N = 0 cannot sit on a log axis; it is clamped into the first bin and the
    // exact zero fraction is carried in the stats histogram instead.
    for (int n : counts) h->Fill(n < 1 ? 1 : n);
    return h;
  }
};

struct Reader {
  TFile* file = nullptr;
  TTree* tree = nullptr;
  std::vector<float> *x = nullptr, *y = nullptr, *zz = nullptr, *bend = nullptr;
  std::vector<int> *isRej = nullptr, *layer = nullptr, *isBarrel = nullptr;

  bool Open(const char* path) {
    file = TFile::Open(path);
    if (!file || file->IsZombie()) return false;
    tree = (TTree*)file->Get(kTreePath);
    if (!tree || !tree->GetBranch("allstub_isRejected")) return false;
    tree->SetBranchStatus("*", 0);
    for (auto b : {"allstub_x", "allstub_y", "allstub_z", "allstub_isRejected",
                   "allstub_layer", "allstub_isBarrel", "allstub_trigBend"})
      tree->SetBranchStatus(b, 1);
    tree->SetBranchAddress("allstub_x", &x);
    tree->SetBranchAddress("allstub_y", &y);
    tree->SetBranchAddress("allstub_z", &zz);
    tree->SetBranchAddress("allstub_isRejected", &isRej);
    tree->SetBranchAddress("allstub_layer", &layer);
    tree->SetBranchAddress("allstub_isBarrel", &isBarrel);
    tree->SetBranchAddress("allstub_trigBend", &bend);
    return true;
  }
  void Close() {
    if (file) file->Close();
    file = nullptr; tree = nullptr;
  }
};

void Run(const char* dir, int nfiles, Hists& acc, Hists& rej, long& nstubAcc,
         long& nstubRej, long& nZeroRej) {
  for (const auto& path : ListRootFiles(dir, nfiles)) {
    Reader R;
    if (!R.Open(path)) { R.Close(); continue; }
    for (Long64_t e = 0; e < R.tree->GetEntries(); ++e) {
      R.tree->GetEntry(e);
      ++acc.nev; ++rej.nev;
      int na = 0, nr = 0;
      const size_t n = R.x->size();
      for (size_t j = 0; j < n; ++j) {
        const bool rejected = R.isRej->at(j) != 0;
        Hists& h = rejected ? rej : acc;
        (rejected ? nr : na)++;
        h.r->Fill(std::sqrt(R.x->at(j) * R.x->at(j) + R.y->at(j) * R.y->at(j)));
        h.z->Fill(R.zz->at(j));
        const int lay = R.layer->at(j);           // barrel L1-L6 -> 0-5, disk D1-D5 -> 6-10
        h.ld->Fill(R.isBarrel->at(j) ? lay - 1 : 5 + lay);
        h.bend->Fill(R.bend->at(j));
      }
      acc.counts.push_back(na);
      rej.counts.push_back(nr);
      nstubAcc += na; nstubRej += nr;
      if (nr == 0) ++nZeroRej;
    }
    R.Close();
    printf("[info] %s done\n", path.Data());
    fflush(stdout);
  }
}

}  // namespace

void compare_stubs(int nfiles = 5, const char* dirDef = "", const char* dirDum = "",
                   const char* outname = "../output/stub_epos.root") {
  if (strlen(dirDef)) kDirDefault = dirDef;
  if (strlen(dirDum)) kDirDummy = dirDum;
  if (!strlen(kDirDefault) || !strlen(kDirDummy)) {
    printf("give the Default-stub and Dummy-stub ntuple locations (a directory of *.root, or one file)\n");
    return;
  }
  gROOT->SetBatch(true);

  Hists defA, defR, dumA, dumR;
  defA.Book("defA"); defR.Book("defR"); dumA.Book("dumA"); dumR.Book("dumR");

  long naDef = 0, nrDef = 0, naDum = 0, nrDum = 0, zDef = 0, zDum = 0;
  printf("=== Default Stub ===\n");
  Run(kDirDefault, nfiles, defA, defR, naDef, nrDef, zDef);
  printf("=== Dummy Stub ===\n");
  Run(kDirDummy, nfiles, dumA, dumR, naDum, nrDum, zDum);

  int nmax = 1;
  for (auto* h : {&defA, &defR, &dumA, &dumR})
    for (int n : h->counts) nmax = std::max(nmax, n);

  printf("\n============ stub yields ============\n");
  printf("%-10s %8s %12s %12s %10s\n", "sample", "events", "accepted/ev", "rejected/ev", "rej frac");
  auto line = [](const char* s, long nev, long na, long nr) {
    printf("%-10s %8ld %12.1f %12.2f %10.4f\n", s, nev,
           nev ? (double)na / nev : 0, nev ? (double)nr / nev : 0,
           (na + nr) ? (double)nr / (na + nr) : 0);
  };
  line("default", defA.nev, naDef, nrDef);
  line("dummy", dumA.nev, naDum, nrDum);
  printf("largest stub count in any event: %d\n", nmax);
  printf("====================================\n");

  TFile* fout = TFile::Open(outname, "RECREATE");
  for (auto* h : {&defA, &defR, &dumA, &dumR}) for (auto* x : h->All()) x->Write();
  defA.MakeFine("defA", nmax)->Write();
  defR.MakeFine("defR", nmax)->Write();
  dumA.MakeFine("dumA", nmax)->Write();
  dumR.MakeFine("dumR", nmax)->Write();

  // scalars the plot step needs, one per bin
  auto* st = new TH1D("stats", "", 8, 0.5, 8.5);
  st->SetBinContent(1, defA.nev); st->SetBinContent(2, naDef);
  st->SetBinContent(3, nrDef);    st->SetBinContent(4, zDef);
  st->SetBinContent(5, dumA.nev); st->SetBinContent(6, naDum);
  st->SetBinContent(7, nrDum);    st->SetBinContent(8, zDum);
  st->Write();
  fout->Close();
  printf("wrote %s\n", outname);
}
