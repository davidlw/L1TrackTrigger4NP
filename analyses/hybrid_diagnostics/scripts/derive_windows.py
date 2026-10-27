#!/usr/bin/env python3
"""Derive |rinv|-binned barrel r-phi match windows from the truth trace.

  python3 derive_windows.py <runs dir with */truthtrace.txt> [--q 0.95]

For each true particle, seed, barrel layer and |rinv| bin (from the particle's
pT, |rinv| = 0.003 B / pT), the smallest |r dphi| of its true-stub match
candidates (the matcher's own fixed-point residual) is the one that decides the
match. Prints per (layer, seed, bin): candidates, coverage of today's window,
the q and 99% quantiles, and the proposed window max(today, q quantile) --
today's value is never lowered, so high-pT behaviour is unchanged. Windows above
the 12-bit residual word (2047 kphi1 rmean) are flagged.

Use traces from a run with an open MatchEngine cut (TT_IDPHICUT=8): with the
default cut the low-pT tails never reach the matcher and the quantiles are biased low.
"""
import argparse, glob, os
from collections import defaultdict
import numpy as np

C = 0.003 * 3.8                                     # |rinv| [1/cm] = C / pT
RINV_EDGES = np.array([0.0, 0.0057, 0.0114, 0.0171, 0.0228])   # pT > 2, 1-2, 0.67-1, 0.5-0.67 GeV
SEEDS = ["L1L2", "L2L3", "L3L4", "L5L6", "D1D2", "D3D4", "L1D1", "L2D1", "L2L3L4", "L4L5L6", "L2L3D1", "D1D2L2"]
# Settings.h (repo, v2): rphimatchcut_[layer][seed] in cm
RPHIMATCHCUT = [
    [0.0, 0.1, 0.07, 0.08, 0.07, 0.05, 0.0, 0.05, 0.08, 0.15, 0.125, 0.15],
    [0.0, 0.0, 0.06, 0.08, 0.05, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0],
    [0.1, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0],
    [0.19, 0.19, 0.0, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.4, 0.4, 0.08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.08, 0.0, 0.0, 0.0],
    [0.5, 0.0, 0.19, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.0]]
IRMEAN = [851, 1269, 1784, 2347, 2936, 3697]
RMEAN = [i * 120.0 / 4096 for i in IRMEAN]
KPHI1 = 2 * (2 * np.pi / 9) / (1 << 17)            # dphisectorHG / 2^nphibitsstub(L6)
WORD_MAX = [2047 * KPHI1 * r for r in RMEAN]         # 12-bit signed residual, in cm


def collect(rundir):
    best = {}                                        # (file, ev, tp, seed, phimin, layer) -> (pt, min |dphi|)
    for path in glob.glob(os.path.join(rundir, "*", "truthtrace.txt")):
        ev, pt = None, {}
        for line in open(path):
            f = line.split()
            if f[0] == "E":
                ev, pt = int(f[1]), {}
            elif f[0] == "T":
                pt[int(f[1])] = float(f[3])
            elif f[0] == "M" and int(f[5]) == int(f[1]) and int(f[1]) in pt:
                key = (path, ev, int(f[1]), int(f[2]), f[3], int(f[4]))
                d = abs(float(f[6]))
                if key not in best or d < best[key][1]:
                    best[key] = (pt[int(f[1])], d)
    res = defaultdict(list)
    for (_, _, _, seed, _, layer), (p, d) in best.items():
        b = np.searchsorted(RINV_EDGES, C / p, side="right") - 1
        if 0 <= b < 4:
            res[(layer, seed, b)].append(d)
    return res


def collect_raw(rundir):
    """(layer, seed) -> list of (|rinv| from the particle pT, smallest true-stub |r dphi|)."""
    best = {}
    for path in glob.glob(os.path.join(rundir, "*", "truthtrace.txt")):
        ev, pt = None, {}
        for line in open(path):
            f = line.split()
            if f[0] == "E":
                ev, pt = int(f[1]), {}
            elif f[0] == "T":
                pt[int(f[1])] = float(f[3])
            elif f[0] == "M" and int(f[5]) == int(f[1]) and int(f[1]) in pt:
                key = (path, ev, int(f[1]), int(f[2]), f[3], int(f[4]))
                d = abs(float(f[6]))
                if key not in best or d < best[key][1]:
                    best[key] = (C / pt[int(f[1])], d)
    out = defaultdict(list)
    for (_, _, _, seed, _, layer), v in best.items():
        out[(layer, seed)].append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--q", type=float, default=0.95)
    ap.add_argument("--min-entries", type=int, default=30)
    a = ap.parse_args()
    res = collect(a.rundir)
    prop = [[[RPHIMATCHCUT[l][s]] * 4 for s in range(12)] for l in range(6)]
    print(f"layer seed     bin  pT[GeV]   n   today  cover  q{int(a.q*100)}   q99   proposed  (word max)")
    for l in range(6):
        for s in range(12):
            today = RPHIMATCHCUT[l][s]
            if today == 0:
                continue
            for b, ptr in enumerate([">2", "1-2", "0.67-1", "0.5-0.67"]):
                v = np.array(res.get((l + 1, s, b), []))
                if len(v) < a.min_entries:
                    flag = " (too few entries: kept)" if b else ""
                    print(f"L{l+1}  {SEEDS[s]:7s}  {b}  {ptr:8s} {len(v):4d}  {today:5.2f}{flag}")
                    if b:  # inherit the previous bin, never narrower than today
                        prop[l][s][b] = prop[l][s][b - 1]
                    continue
                q, q99 = np.quantile(v, [a.q, 0.99])
                w = max(today, round(q, 2)) if b else today
                prop[l][s][b] = max(w, prop[l][s][b - 1]) if b else w
                over = " OVER WORD" if prop[l][s][b] > WORD_MAX[l] else ""
                print(f"L{l+1}  {SEEDS[s]:7s}  {b}  {ptr:8s} {len(v):4d}  {today:5.2f}  {np.mean(v <= today):5.2f}  {q:5.2f}  {q99:5.2f}   {prop[l][s][b]:5.2f}   ({WORD_MAX[l]:.2f}){over}")
    print("\n// proposed: rphimatchcut_[layer][seed][rinvbin], rinvbin = |rinv| in 4 equal bins up to 0.0228/cm")
    print("// (pT > 2, 1-2, 0.67-1, 0.5-0.67 GeV); seeds without data keep today's value in every bin")
    print("std::array<std::array<std::array<double, N_RINVBIN>, N_SEED>, N_LAYER> rphimatchcut_{{")
    for l in range(6):
        cells = ", ".join("{{" + ", ".join(f"{x:g}" for x in prop[l][s]) + "}}" for s in range(12))
        print(f"    {{{{{cells}}}}},  // Layer {l+1}")
    print("}};")


if __name__ == "__main__":
    main()
