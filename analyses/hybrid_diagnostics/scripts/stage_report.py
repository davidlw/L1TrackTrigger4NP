#!/usr/bin/env python3
"""Stage-by-stage report from the truth trace (run analyze_trace.py first).

  python3 stage_report.py <dir with trace_rows.pkl> [--eta-max 0.8]

Prints, for true muons in |eta| < eta-max per pT bin: the funnel (seedable ->
seed pair -> calculator -> any match -> >= 2 matched layers -> found), and per
seed -> layer: the fraction of routed projections whose true stub becomes a
match candidate (MatchEngine) and is matched (MatchProcessor), with the r-phi /
z pass fractions and the r-phi residual in units of the match window.
"""
import argparse, os, pickle
from collections import defaultdict
import numpy as np

SEEDN = {0: "L1L2", 1: "L2L3", 2: "L3L4", 3: "L5L6"}
BINS = [(0.5, 0.75), (0.75, 1.0), (1.0, 2.0), (2.0, 100.0)]


def pbin(pt):
    for lo, hi in BINS:
        if lo <= pt < hi:
            return f"{lo:g}-{hi:g}"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--eta-max", type=float, default=0.8)
    a = ap.parse_args()
    rows = [r for r in pickle.load(open(os.path.join(a.dir, "trace_rows.pkl"), "rb")) if abs(r["eta"]) < a.eta_max]

    print(f"== funnel, true muons |eta| < {a.eta_max} (barrel seeds only are traced)")
    print("   pT bin        n  seedable  seed-pair  calc-pass  any-match  >=2 match  found")
    for lo, hi in BINS:
        R = [r for r in rows if lo <= r["pt"] < hi]
        if not R:
            continue
        f = lambda k: np.mean([bool(r[k]) for r in R])
        fd = np.mean([r["found"] for r in R if r["found"] is not None])
        print(f"   {lo:g}-{hi:<6g} {len(R):5d}  {f('seedable'):8.2f}  {f('tried'):9.2f}  {f('passed'):9.2f}"
              f"  {f('reached_builder'):9.2f}  {f('builder_ok'):9.2f}  {fd:5.2f}")

    cnt = defaultdict(lambda: [0, 0, 0])
    mt = defaultdict(list)
    for r in rows:
        pb = pbin(r["pt"])
        if pb is None:
            continue
        m, tp = r["m"], r["key"][2]
        for (s, sec, lay) in m["A"]:
            if (s, sec) in m["F"] or lay not in m["stubs"]:
                continue
            c = cnt[(s, lay, pb)]
            c[0] += 1
            cands = [x for x in m["M"] if x[0] == s and x[2] == lay and x[3] == tp]
            c[1] += bool(cands)
            c[2] += any(x[8] for x in cands)
        for (s, pm, lay, stp, dI, cut, dzI, zc, im, dEx, dAp) in m["M"]:
            if stp == tp:
                mt[(s, lay, pb)].append((abs(dI), cut, abs(dzI), zc, im))

    print("\n== routed projections with a true stub on the layer")
    print("   seed->layer pT          routed  candidate  matched | pass r-phi  pass z  |dphi|/window q50 q90  window r-phi")
    order = [f"{lo:g}-{hi:g}" for lo, hi in BINS]
    for k in sorted(cnt, key=lambda k: (k[0], k[1], order.index(k[2]))):
        n, c, m_ = cnt[k]
        if n < 20:
            continue
        s, lay, pb = k
        line = f"   {SEEDN[s]}->L{lay}    {pb:10s} {n:6d}  {c / n:9.2f}  {m_ / n:7.2f} |"
        if len(mt[k]) >= 10:
            x = np.array(mt[k])
            q = np.quantile(x[:, 0] / x[:, 1], [0.5, 0.9])
            line += f"  {np.mean(x[:, 0] <= x[:, 1]):9.2f}  {np.mean(x[:, 2] < x[:, 3]):6.2f}  {q[0]:12.1f} {q[1]:4.1f}  {x[0, 1]:8.2f} cm"
        print(line)


if __name__ == "__main__":
    main()
