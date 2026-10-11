#!/usr/bin/env python3
"""Fit the barrel r-phi match window w(rinv) = sqrt(a^2 + (b (|rinv|/rinvmax)^n)^2) per (layer, seed).

  python3 fit_window_slope.py <runs dir with */truthtrace.txt> [--q 0.95] [--n 3]

a = today's rphimatchcut_ (kept). b [cm] = extra window at |rinv| = rinvmax (pT = ptmin = 0.5 GeV).
The true-stub residual is the matcher's fixed-point one for inner layers and the
exact-helix one where the fixed-point projection is known to be off (--exact-for,
default: L4-L6 from L1L2, L5-L6 from L3L4), i.e. assuming that projection is fixed.
Prints the fitted b, measured vs fitted windows per |rinv| bin, and the windows
at the LUT bin centres with the 10-bit table / 12-bit residual limits.
"""
import argparse, glob, os
from collections import defaultdict
import numpy as np
from derive_windows import RPHIMATCHCUT, RMEAN, KPHI1, SEEDS, C

RINVMAX = C / 0.5
KRINV = KPHI1 / (120.0 / 4096) * 2 ** -8
SHIFT = 12
CENTRES = (np.arange(4) + 0.5) * (1 << SHIFT) * KRINV          # LUT bin centres
FITBINS = np.linspace(0.0, 0.0228, 9)
EXACT = {(4, 0), (5, 0), (6, 0), (5, 2), (6, 2)}                # (layer, seed): use exact-helix residual


def collect(rundir):
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
                fixed, exact = abs(float(f[6])), abs(float(f[11]))
                if key not in best or fixed < best[key][1]:
                    best[key] = (C / pt[int(f[1])], fixed, exact)
    out = defaultdict(list)
    for (_, _, _, seed, _, layer), v in best.items():
        out[(layer, seed)].append(v)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--q", type=float, default=0.95)
    ap.add_argument("--n", type=float, default=3.0)
    ap.add_argument("--min-entries", type=int, default=25)
    a_ = ap.parse_args()
    raw = collect(a_.rundir)
    print(f"LUT |rinv| bin centres {np.round(CENTRES, 4)} /cm = pT {np.round(C / CENTRES, 2)} GeV;  n = {a_.n}")
    print("layer seed     resid   a     b[cm] | measured q95 vs fit in fit bins (pT)            | window at LUT centres [cm]  table/resid max")
    B = np.zeros((6, 12))
    for l in range(6):
        for s in range(12):
            a = RPHIMATCHCUT[l][s]
            v = raw.get((l + 1, s))
            if a == 0 or not v:
                continue
            v = np.array(v)
            col = 2 if (l + 1, s) in EXACT else 1
            r, q = [], []
            for lo, hi in zip(FITBINS[:-1], FITBINS[1:]):
                m = (v[:, 0] >= lo) & (v[:, 0] < hi)
                if m.sum() >= a_.min_entries:
                    r.append(0.5 * (lo + hi)); q.append(np.quantile(v[m, col], a_.q))
            r, q = np.array(r), np.array(q)
            x = (r / RINVMAX) ** a_.n
            use = q > a
            b = np.sqrt(np.sum((q[use] ** 2 - a ** 2) * x[use] ** 2) / np.sum(x[use] ** 4)) if use.any() else 0.0
            B[l, s] = b
            fit = np.sqrt(a ** 2 + (b * x) ** 2)
            w = np.sqrt(a ** 2 + (b * (CENTRES / RINVMAX) ** a_.n) ** 2); w[0] = a
            tmax, rmax = 1023 * KPHI1 * RMEAN[l], 2047 * KPHI1 * RMEAN[l]
            cmp = " ".join(f"{qq:.2f}/{ff:.2f}" for qq, ff in zip(q[-4:], fit[-4:]))
            flag = (" >table" if w.max() > tmax else "") + (" >resid" if w.max() > rmax else "")
            print(f"L{l+1}  {SEEDS[s]:7s} {'exact' if col == 2 else 'fixed'}  {a:4.2f}  {b:5.2f} | {cmp:46s} | {np.round(w, 2)}  {tmax:.2f}/{rmax:.2f}{flag}")
    print(f"\n// extra barrel r-phi match window at |rinv| = rinvmax [cm]; w = sqrt(a^2 + (b (|rinv|/rinvmax)^{a_.n:g})^2)")
    print("std::array<std::array<double, N_SEED>, N_LAYER> rphimatchcutlowpt_{")
    for l in range(6):
        print("    {{" + ", ".join(f"{x:.2f}" for x in B[l]) + f"}}}},  //Layer {l+1}")
    print("};")


if __name__ == "__main__":
    main()
