#!/usr/bin/env python3
"""Counts of locally reconstructed HGCAL objects (rechits, layer clusters) per event,
to size a 'local reco instead of raw' format.  Usage: MAXEV=n python3 dump_localreco.py out.jsonl files..."""
import os, sys, json
import ROOT
ROOT.gROOT.SetBatch(True)
from DataFormats.FWLite import Events, Handle
RH = "edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >"
prods = {k: (Handle(RH), ("HGCalRecHit", k, "HLT")) for k in ("HGCEERecHits", "HGCHEFRecHits", "HGCHEBRecHits")}
prods["layerClusters"] = (Handle("std::vector<reco::CaloCluster>"), ("hgcalMergeLayerClusters", "", "HLT"))
MAXEV = int(os.environ.get("MAXEV", "-1"))
out = open(sys.argv[1], "w")
for f in sys.argv[2:]:
    for i, ev in enumerate(Events("file:" + f)):
        if 0 <= MAXEV <= i: break
        r = {"file": f, "evt": ev.eventAuxiliary().event()}
        for k, (h, lab) in prods.items():
            r[k] = h.product().size() if ev.getByLabel(lab, h) and h.isValid() else -1
        out.write(json.dumps(r) + "\n")
out.close()
