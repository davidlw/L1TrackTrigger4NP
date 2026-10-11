#!/usr/bin/env python3
"""Hybrid efficiency, fakes and duplicates per cut variant, from the diag_ntuple.root files.

  python3 compare_variants.py <runs dir> job_ me6_ win4_ both_ ...

Each tag selects <runs dir>/<tag><i>/diag_ntuple.root (the baseline jobs are 'job').
Efficiency: TPs with |eta| < 2.4, found = tp_nmatch > 0 (compare_eff.C definition).
Fake: tracks that are not loosely genuine. Duplicate: TPs matched by more than one track.
"""
import glob, os, sys
import numpy as np
import uproot

BR = ["tp_pt", "tp_eta", "tp_nmatch", "trk_pt", "trk_eta", "trk_loose", "trk_nstub", "trk_matchtp_pt"]
PT = [(0.3, 0.5), (0.5, 0.75), (0.75, 1.0), (1.0, 2.0), (2.0, 100.0)]


def load(files):
    out = {b: [] for b in BR}
    nev = 0
    for f in files:
        a = uproot.open(f)["L1TrackHitNtupleMaker/eventTree"].arrays(BR, library="np")
        nev += len(a["tp_pt"])
        for b in BR:
            out[b].append(np.concatenate(a[b]) if len(a[b]) else np.zeros(0))
    return {b: np.concatenate(v) for b, v in out.items()}, nev


def main():
    rundir, tags = sys.argv[1], sys.argv[2:]
    print(f"{'variant':8s} {'events':>6s}  " + "  ".join(f"{r} {lo:g}-{hi:g}" for r in ["bar", "fwd"] for lo, hi in PT[1:4])
          + "   <nstub> bar 0.5-1  fake  dup  trk/evt")
    for tag in tags:
        files = sorted(glob.glob(os.path.join(rundir, f"{tag}[0-9]*", "diag_ntuple.root")))
        if not files:
            print(f"{tag:8s} no files")
            continue
        a, nev = load(files)
        acc = np.abs(a["tp_eta"]) < 2.4
        line = f"{tag:8s} {nev:6d}  "
        for reg in [np.abs(a["tp_eta"]) < 0.8, np.abs(a["tp_eta"]) > 1.8]:
            for lo, hi in PT[1:4]:
                s = acc & reg & (a["tp_pt"] >= lo) & (a["tp_pt"] < hi)
                line += f"   {np.mean(a['tp_nmatch'][s] > 0):.2f}({s.sum()})"
        gb = (a["trk_loose"] == 1) & (np.abs(a["trk_eta"]) < 0.8) & (a["trk_matchtp_pt"] >= 0.5) & (a["trk_matchtp_pt"] < 1.0)
        found = acc & (a["tp_nmatch"] > 0)
        line += f"   {a['trk_nstub'][gb].mean() if gb.any() else float('nan'):6.2f}" \
                f"  {np.mean(a['trk_loose'] == 0) if len(a['trk_loose']) else float('nan'):.3f}" \
                f"  {np.mean(a['tp_nmatch'][found] > 1) if found.any() else float('nan'):.3f}" \
                f"  {len(a['trk_pt']) / max(nev, 1):.2f}"
        print(line)


if __name__ == "__main__":
    main()
