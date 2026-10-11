#!/usr/bin/env python3
"""Per-event content of the re-slimmed MINIAODSIM files: packed PF candidates (charged,
with track details), lost tracks, and generated charged particles, to judge how
complete the standard (pp-flavour) miniAOD is for heavy-ion events.
Usage: python3 count_miniaod.py file_sim.root [maxEvents]"""
import sys, statistics as st
import ROOT
ROOT.gROOT.SetBatch(True)
from DataFormats.FWLite import Events, Handle
pf, lost, gen = Handle("std::vector<pat::PackedCandidate>"), Handle("std::vector<pat::PackedCandidate>"), Handle("std::vector<pat::PackedGenParticle>")
maxev = int(sys.argv[2]) if len(sys.argv) > 2 else -1
rows = []
for i, ev in enumerate(Events(sys.argv[1])):
    if 0 <= maxev <= i: break
    ev.getByLabel("packedPFCandidates", pf); ev.getByLabel("lostTracks", lost); ev.getByLabel("packedGenParticles", gen)
    c = pf.product()
    nch = sum(1 for p in c if p.charge() != 0)
    ntd = sum(1 for p in c if p.charge() != 0 and p.hasTrackDetails())
    nlost = lost.product().size()
    g = gen.product()
    ngch = sum(1 for p in g if p.charge() != 0 and p.status() == 1 and abs(p.eta()) < 2.5 and p.pt() > 0.3)
    nrch = sum(1 for p in c if p.charge() != 0 and abs(p.eta()) < 2.5 and p.pt() > 0.3)
    rows.append((c.size(), nch, ntd, nlost, ngch, nrch))
names = ["packedPF", "charged", "withTrackDetails", "lostTracks", "gen ch |eta|<2.5 pT>0.3", "PF ch |eta|<2.5 pT>0.3"]
print(sys.argv[1], len(rows), "events")
for k, n in enumerate(names):
    v = [r[k] for r in rows]
    print(f"  {n:28s} mean {st.mean(v):8.0f}  max {max(v):8.0f}")
