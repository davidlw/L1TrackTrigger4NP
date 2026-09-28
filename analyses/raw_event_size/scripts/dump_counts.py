#!/usr/bin/env python3
"""Per-event inventory for the Phase-2 raw event size estimate.

Reads step-2 (GEN-SIM-DIGI-RAW) files and writes one JSON line per event with
  - fed:    {fedId: bytes} for every non-empty FED in rawDataCollector
  - counts: hit / digi / cluster multiplicities for the subdetectors that have
            no Phase-2 packer in CMSSW (IT, OT, HGCAL, MTD), plus the ones that
            do, as a cross-check.

Usage (inside cmsenv):  python3 dump_counts.py out.jsonl file1.root [file2.root ...]
Set MAXEV=<n> to stop after n events per file.
"""
import os, sys, json
import ROOT
ROOT.gROOT.SetBatch(True)
from DataFormats.FWLite import Events, Handle

# MuonDigiCollection ranges are not usable from Python; count them in C++.
ROOT.gInterpreter.Declare('''
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include <iterator>
template <class C> size_t countMuonDigis(const C& c) {
  size_t n = 0;
  for (auto it = c.begin(); it != c.end(); ++it) {
    auto r = (*it).second;
    n += std::distance(r.first, r.second);
  }
  return n;
}
''')

H = lambda t: Handle(t)
products = {
    "raw":   (H("FEDRawDataCollection"), ("rawDataCollector", "", "HLT")),
    "itDigi":(H("edm::DetSetVector<PixelDigi>"), ("simSiPixelDigis", "Pixel", "HLT")),
    "itClu": (H("edmNew::DetSetVector<SiPixelCluster>"), ("siPixelClusters", "", "HLT")),
    "otDigi":(H("edm::DetSetVector<Phase2TrackerDigi>"), ("mix", "Tracker", "HLT")),
    "otClu": (H("edmNew::DetSetVector<Phase2TrackerCluster1D>"), ("siPhase2Clusters", "", "HLT")),
    "hgcEE": (H("edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >"), ("simHGCalUnsuppressedDigis", "EE", "HLT")),
    "hgcHEF":(H("edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >"), ("simHGCalUnsuppressedDigis", "HEfront", "HLT")),
    "hgcHEB":(H("edm::SortedCollection<HGCDataFrame<DetId,HGCSample>,edm::StrictWeakOrdering<HGCDataFrame<DetId,HGCSample> > >"), ("simHGCalUnsuppressedDigis", "HEback", "HLT")),
    "btl":   (H("edm::SortedCollection<FTLDataFrameT<BTLDetId,BTLSample,mtdhelpers::BTLRowColDecode>,edm::StrictWeakOrdering<FTLDataFrameT<BTLDetId,BTLSample,mtdhelpers::BTLRowColDecode> > >"), ("mix", "FTLBarrel", "HLT")),
    "etl":   (H("edm::SortedCollection<FTLDataFrameT<ETLDetId,ETLSample,mtdhelpers::ETLRowColDecode>,edm::StrictWeakOrdering<FTLDataFrameT<ETLDetId,ETLSample,mtdhelpers::ETLRowColDecode> > >"), ("mix", "FTLEndcap", "HLT")),
    "eb":    (H("EBDigiCollection"), ("simEcalDigis", "ebDigis", "HLT")),
    "ebU":   (H("EBDigiCollection"), ("simEcalUnsuppressedDigis", "", "HLT")),
    "hbhe":  (H("HcalDataFrameContainer<QIE11DataFrame>"), ("simHcalDigis", "HBHEQIE11DigiCollection", "HLT")),
    "hf":    (H("HcalDataFrameContainer<QIE10DataFrame>"), ("simHcalDigis", "HFQIE10DigiCollection", "HLT")),
    "zdc":   (H("edm::SortedCollection<ZDCDataFrame,edm::StrictWeakOrdering<ZDCDataFrame> >"), ("simHcalUnsuppressedDigis", "", "HLT")),
    "gem":   (H("MuonDigiCollection<GEMDetId,GEMDigi>"), ("simMuonGEMDigis", "", "HLT")),
    "dt":    (H("MuonDigiCollection<DTLayerId,DTDigi>"), ("simMuonDTDigis", "", "HLT")),
    "rpc":   (H("MuonDigiCollection<RPCDetId,RPCDigi>"), ("simMuonRPCDigis", "", "HLT")),
    "cscS":  (H("MuonDigiCollection<CSCDetId,CSCStripDigi>"), ("simMuonCSCDigis", "MuonCSCStripDigi", "HLT")),
    "cscW":  (H("MuonDigiCollection<CSCDetId,CSCWireDigi>"), ("simMuonCSCDigis", "MuonCSCWireDigi", "HLT")),
}

def get(ev, key):
    h, lab = products[key]
    ok = ev.getByLabel(lab, h)
    return h.product() if ok and h.isValid() else None

def dsv_counts(dsv):
    """(#modules with hits, #entries) for an edm::DetSetVector."""
    nmod = nent = 0
    for ds in dsv:
        n = ds.size()
        if n: nmod += 1; nent += n
    return nmod, nent

def new_dsv_counts(dsv):
    nmod = nent = 0
    for ds in dsv:
        n = ds.size()
        if n: nmod += 1; nent += n
    return nmod, nent

def ot_split(dsv):
    """Split OT digis into macro-pixel (PS-p, column>1 seen) vs strip sensors."""
    out = {"stripMod": 0, "stripDigi": 0, "mpxMod": 0, "mpxDigi": 0}
    for ds in dsv:
        n = ds.size()
        if not n: continue
        mpx = any(ds.data[i].column() > 1 for i in range(n))
        k = "mpx" if mpx else "strip"
        out[k + "Mod"] += 1; out[k + "Digi"] += n
    return out

def ot_clu_split(dsv):
    """Clusters per module, with a per-cluster width sum (strip count)."""
    ncl = nmod = width = 0
    for ds in dsv:
        n = ds.size()
        if not n: continue
        nmod += 1; ncl += n
        for c in ds: width += c.size()
    return nmod, ncl, width

# In-time ADC thresholds scanned for HGCAL zero suppression. ADC LSB is
# 68.75 fC / 1024 = 0.067 fC; 1 MIP in 300 um Si is ~3.6 fC ~ 54 ADC, less
# in 200/120 um sensors, so these span roughly 0.15 - 2 MIP.
HGC_THR = [8, 12, 16, 20, 27, 40, 54, 80, 108]

def hgc_counts(coll):
    """Frames stored, and frames whose in-time sample (or TOT) passes each threshold."""
    nframe = ntot = 0
    above = [0] * len(HGC_THR)
    for f in coll:
        nframe += 1
        # in-time sample: index 2 of the 5 stored by the 14_0_X digitizer
        s = f.sample(2 if f.size() >= 5 else f.size() // 2)
        if s.mode():             # TOT mode: large signal, always kept
            ntot += 1
            a = 1 << 30
        else:
            a = s.data()
        for i, t in enumerate(HGC_THR):
            if a >= t: above[i] += 1
    return {"frames": nframe, "tot": ntot, "above": above}

MAXEV = int(os.environ.get("MAXEV", "-1"))
out = open(sys.argv[1], "w")
for fname in sys.argv[2:]:
    events = Events(fname if fname.startswith("root:") or fname.startswith("file:") else "file:" + fname)
    for iev, ev in enumerate(events):
        if 0 <= MAXEV <= iev: break
        rec = {"file": fname, "run": ev.eventAuxiliary().run(), "evt": ev.eventAuxiliary().event()}
        raw = get(ev, "raw")
        fed = {}
        if raw is not None:
            for i in range(0, 4096):
                s = raw.FEDData(i).size()
                if s: fed[i] = s
        rec["fed"] = fed
        c = {}
        p = get(ev, "itDigi");  c["itMod"], c["itDigi"] = dsv_counts(p) if p is not None else (-1, -1)
        p = get(ev, "itClu");   c["itCluMod"], c["itClu"] = new_dsv_counts(p) if p is not None else (-1, -1)
        p = get(ev, "otDigi");  c.update({"ot_" + k: v for k, v in (ot_split(p) if p is not None else {}).items()})
        p = get(ev, "otClu")
        if p is not None: c["otCluMod"], c["otClu"], c["otCluWidth"] = ot_clu_split(p)  # absent in some RelVals
        for k in ("hgcEE", "hgcHEF", "hgcHEB"):
            p = get(ev, k)
            if p is not None: c[k] = hgc_counts(p)
        for k in ("btl", "etl", "eb", "ebU", "hbhe", "hf", "zdc"):
            p = get(ev, k); c[k] = p.size() if p is not None else -1
        for k in ("gem", "dt", "rpc", "cscS", "cscW"):
            p = get(ev, k)
            if p is None:
                c[k] = -1
                continue
            c[k] = int(ROOT.countMuonDigis(p))
        rec["counts"] = c
        out.write(json.dumps(rec) + "\n")
        out.flush()
out.close()
