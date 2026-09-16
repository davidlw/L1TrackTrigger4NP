// ---------------------------------------------------------------------------
// InputFiles.h -- shared by every fill macro under analyses/.
//
// Resolve the `input` argument of a macro to the list of ntuple files to read:
//
//   a single .root file   -> that file
//   a directory           -> every *.root inside it, sorted by name, and the
//                            first `nfiles` of them (nfiles <= 0 means all)
//   a comma-separated mix -> the union, in order
//
// Nothing here assumes lxplus or a naming pattern: three files called
// _4, _17 and _22 copied to a laptop are read just as well as /eos/.../0000.
// What was found is printed, so a wrong path shows up as a message instead of
// as an empty histogram.
//
// Include from a macro run inside its scripts/ directory:
//   #include "../../common/InputFiles.h"
// ---------------------------------------------------------------------------
#ifndef L1TT4NP_INPUTFILES_H
#define L1TT4NP_INPUTFILES_H

#include <algorithm>
#include <cstdio>
#include <vector>

#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

inline std::vector<TString> ListRootFiles(const char* input, int nfiles = 0) {
  std::vector<TString> out;
  TString in(input);
  in = in.Strip(TString::kBoth);
  if (in.IsNull()) {
    printf("[InputFiles] no input given\n");
    return out;
  }
  TObjArray* parts = in.Tokenize(",");
  for (int ip = 0; ip < parts->GetEntries(); ++ip) {
    TString item = ((TObjString*)parts->At(ip))->GetString().Strip(TString::kBoth);
    if (item.IsNull())
      continue;
    while (item.Length() > 1 && item.EndsWith("/"))
      item.Remove(item.Length() - 1);

    if (gSystem->AccessPathName(item)) {  // true when it does NOT exist
      printf("[InputFiles] not found: %s\n", item.Data());
      continue;
    }
    FileStat_t st;
    gSystem->GetPathInfo(item, st);
    if (!R_ISDIR(st.fMode)) {
      out.push_back(item);
      continue;
    }

    std::vector<TString> names;
    TSystemDirectory d("d", item);
    TList* fl = d.GetListOfFiles();
    if (fl) {
      TIter next(fl);
      while (auto* o = (TSystemFile*)next()) {
        TString n = o->GetName();
        if (!o->IsDirectory() && n.EndsWith(".root"))
          names.push_back(n);
      }
    }
    std::sort(names.begin(), names.end());
    for (auto& n : names) {
      if (nfiles > 0 && (int)out.size() >= nfiles)
        break;
      out.push_back(item + "/" + n);
    }
  }
  delete parts;
  if (out.empty())
    printf("[InputFiles] no .root files found for '%s'\n", input);
  else
    printf("[InputFiles] %zu file(s) from '%s' (first: %s)\n", out.size(), input,
           gSystem->BaseName(out.front()));
  return out;
}

// For paired productions of the SAME events (default vs dummy stubs): keep in
// both lists only the basenames present in both, so that file k of one is
// file k of the other.
inline void KeepCommonFiles(std::vector<TString>& a, std::vector<TString>& b) {
  auto base = [](const TString& p) { return TString(gSystem->BaseName(p)); };
  std::vector<TString> ba, bb;
  for (auto& p : a) ba.push_back(base(p));
  for (auto& p : b) bb.push_back(base(p));
  std::vector<TString> a2, b2;
  for (size_t i = 0; i < a.size(); ++i)
    if (std::find(bb.begin(), bb.end(), ba[i]) != bb.end()) a2.push_back(a[i]);
  for (size_t i = 0; i < b.size(); ++i)
    if (std::find(ba.begin(), ba.end(), bb[i]) != ba.end()) b2.push_back(b[i]);
  if (a2.size() != a.size() || b2.size() != b.size())
    printf("[InputFiles] paired samples: keeping %zu common file(s) (had %zu and %zu)\n", a2.size(), a.size(),
           b.size());
  a.swap(a2);
  b.swap(b2);
}

#endif
