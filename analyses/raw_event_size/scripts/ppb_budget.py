#!/usr/bin/env python3
"""pPb (EPOS MB) stored event size and throughput at a given L1 rate, step by step.

Steps 0-1 come from the model. Steps 2-5 add assumptions, marked below.
Usage: python3 ppb_budget.py [datadir] [rate_Hz] [budget_GBps]   (defaults ../data 750e3 51)
"""
import math, os, statistics as st, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
datadir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, "..", "data")
R = float(sys.argv[2]) if len(sys.argv) > 2 else 750e3
B = float(sys.argv[3]) * 1e9 if len(sys.argv) > 3 else 51e9
sys.argv = [sys.argv[0], datadir]
import size_model as M

m = M.Model(M.load("pu140.jsonl"), M.load("eb_pu140.jsonl"), 0)
ev, eb = M.load("ppb.jsonl"), M.load("eb_ppb.jsonl")
C = lambda k: st.mean(e["counts"][k] for e in ev)
neb3, neb5 = M.eb_mean_above(eb, 12), M.eb_mean_above(eb, 20)
base = st.mean(sum(M.scenario_sizes(m, e["counts"], neb3, False, False).values()) for e in ev)
rows = [M.scenario_sizes(m, e["counts"], neb3, True, False) for e in ev]
mc = {d: st.mean(r[d] for r in rows) for d in rows[0]}
hgc = st.mean(M.hgc_hits(e["counts"], m.hgc_thr) for e in ev)

HCAL_ZS, L1_FINAL = 0.05e6, 0.02e6          # assumptions
IT_B_PER_CLU, OT_B_PER_CLU, HGC_B_PER_HIT = 6, 3, 5   # compact local-reco formats (assumption)
COMPRESSION = 1.5                            # assumption, must be measured


def line(name, b):
    print(f"  {name:58s} {b / 1e3:7.0f} kB {b * R / 1e9:8.0f} GB/s {b * R / B:6.1f}x budget")


print(f"EPOS pPb MB, {len(ev)} events, L1 rate {R / 1e3:.0f} kHz, storage budget {B / 1e9:.0f} GB/s")
line("0. Baseline raw", base)
s1 = sum(mc.values()); line("1. Empty headers dropped (lossless)", s1)
s2 = s1 - mc["HCAL (HB+HO+HF)"] + HCAL_ZS - mc["L1 trigger + TPGs"] + L1_FINAL
line("2. + HCAL ZS, L1 final objects only", s2)
lr_it, lr_ot, lr_hg = C("itClu") * IT_B_PER_CLU, C("otClu") * OT_B_PER_CLU, hgc * HGC_B_PER_HIT
s3 = s2 - mc["Inner Tracker"] - mc["Outer Tracker"] - mc["HGCAL"] + lr_it + lr_ot + lr_hg
line("3. + local reco IT/OT/HGCAL", s3)
s4 = s3 - mc["ECAL Barrel"] + neb5 * m.k["ECAL Barrel"]
line("4. + ECAL barrel ZS at 5 sigma instead of 3", s4)
line(f"5. + lossless compression x{COMPRESSION}", s4 / COMPRESSION)
mtd = mc["MTD BTL"] + mc["MTD ETL"]
print("  tracker (+MTD) stream:")
line("T1. IT+OT compact raw + MTD + L1 final", mc["Inner Tracker"] + mc["Outer Tracker"] + mtd + L1_FINAL)
line("T2. IT+OT clusters + MTD + L1 final", lr_it + lr_ot + mtd + L1_FINAL)
hits = s1 - mc["HCAL (HB+HO+HF)"] - mc["L1 trigger + TPGs"] - mc["ECAL Barrel"]
for mu in (0.06, 0.28):
    nu = mu / (1 - math.exp(-mu))
    print(f"  pileup mu={mu}: {nu:.3f} collisions per triggered crossing -> +{(nu - 1) * hits / 1e3:.0f} kB of hits")
