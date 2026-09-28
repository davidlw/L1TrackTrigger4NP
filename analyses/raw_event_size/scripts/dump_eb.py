#!/usr/bin/env python3
"""ECAL barrel occupancy for the Phase-2 zero-suppression estimate.

The step-2 files carry the legacy 10-sample EB digis; the Phase-2 barrel has no
packer. For each event, count crystals whose pedestal-subtracted peak exceeds a
scan of ADC thresholds, and accumulate the noise RMS of one in-time sample so
the threshold can be expressed in units of sigma.

Usage (inside cmsenv):  MAXEV=n python3 dump_eb.py out.jsonl file1.root [...]
"""
import os, sys, json
import ROOT
ROOT.gROOT.SetBatch(True)
from DataFormats.FWLite import Events, Handle

EB_THR = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50]
ROOT.gInterpreter.Declare('''
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include <vector>
#include <algorithm>
std::vector<double> ebScan(const EBDigiCollection& c, const std::vector<int>& thr) {
  // out = [n above thr_0, ..., n above thr_k, sum(d), sum(d^2), n] with d = sample(5) - ped
  std::vector<double> out(thr.size() + 3, 0.);
  for (unsigned i = 0; i < c.size(); ++i) {
    EBDataFrame f(c[i]);
    if (f.size() < 10) continue;
    bool g12 = true;
    for (int s = 0; s < f.size(); ++s) if (f.sample(s).gainId() != 1) g12 = false;
    double ped = (f.sample(0).adc() + f.sample(1).adc() + f.sample(2).adc()) / 3.;
    double peak = 1e9;
    if (g12) {
      peak = -1e9;
      for (int s = 3; s < f.size(); ++s) peak = std::max(peak, f.sample(s).adc() - ped);
      double d = f.sample(5).adc() - ped;
      out[thr.size()] += d; out[thr.size() + 1] += d * d; out[thr.size() + 2] += 1;
    }
    for (unsigned t = 0; t < thr.size(); ++t) if (peak >= thr[t]) out[t] += 1;
  }
  return out;
}
''')
thr = ROOT.std.vector("int")()
for t in EB_THR: thr.push_back(t)

h = Handle("EBDigiCollection")
MAXEV = int(os.environ.get("MAXEV", "-1"))
out = open(sys.argv[1], "w")
for fname in sys.argv[2:]:
    for iev, ev in enumerate(Events("file:" + fname)):
        if 0 <= MAXEV <= iev: break
        if not ev.getByLabel(("simEcalUnsuppressedDigis", "", "HLT"), h): continue
        r = list(ROOT.ebScan(h.product(), thr))
        n = len(EB_THR)
        out.write(json.dumps({"file": fname, "evt": ev.eventAuxiliary().event(),
                              "thr": EB_THR, "above": r[:n],
                              "noiseSum": r[n], "noiseSum2": r[n + 1], "noiseN": r[n + 2]}) + "\n")
out.close()
