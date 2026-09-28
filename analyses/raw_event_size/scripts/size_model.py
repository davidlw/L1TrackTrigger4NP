#!/usr/bin/env python3
"""Phase-2 raw event size model, per subdetector.

Input: the per-event JSON lines written by dump_counts.py and dump_eb.py.
Each subdetector's size is modelled as

    size = fixed + per_hit * N_hits

`fixed` is the readout overhead sent for every L1 accept regardless of
occupancy: front-end chip headers, empty sub-packets, and so on. It is
built bottom-up from the front-end data formats. For HCAL and L1 the whole
sub-event is fixed. `per_hit` (or, for HGCAL, the zero-suppression
threshold) is calibrated so that the model reproduces the Phase-2 DAQ TDR
(CMS-TDR-022, Table 3.2) sub-event size. The TDR sizes were simulated at
PU200 and scaled linearly to PU140 (i.e. assuming no fixed part), so we
calibrate against the PU200 column, using the official PU140 ttbar hit
counts scaled by 200/140.

This matters because only the muon, HCAL, legacy-ECAL and L1 FEDs have
packers in CMSSW; the IT, OT, HGCAL, MTD and Phase-2 EB have none, so their
sizes cannot be read from rawDataCollector.

Usage: python3 size_model.py [datadir]   (defaults to ../data)
"""
import json, math, os, statistics as st, sys

DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "data")
MB = 1e6

# ---------------------------------------------------------------------------
# Reference: DAQ TDR Table 3.2, sub-event sizes (MB).
TDR_PU140 = {
    "Inner Tracker": 1.01, "Outer Tracker": 0.80, "MTD BTL": 0.17, "MTD ETL": 0.31,
    "ECAL Barrel": 0.42, "HCAL (HB+HO+HF)": 0.33, "HGCAL": 2.10,
    "Muons (DT+CSC+GEM+RPC)": 0.53, "L1 trigger + TPGs": 0.42,
}
TDR_PU200 = {
    "Inner Tracker": 1.44, "Outer Tracker": 1.15, "MTD BTL": 0.24, "MTD ETL": 0.44,
    "ECAL Barrel": 0.60, "HCAL (HB+HO+HF)": 0.33, "HGCAL": 3.00,
    "Muons (DT+CSC+GEM+RPC)": 0.755, "L1 trigger + TPGs": 0.47,
}
PU_SCALE = 200 / 140      # PU140 sample occupancy -> PU200

# ---------------------------------------------------------------------------
# Fixed-overhead ingredients (per L1 accept). Each has a low/high variant
# used for the systematic band; see README for the reasoning.
P = {
    # IT: ~10k CROC chips (TDR); each sends at least an event header word.
    "it_chips": 10000, "it_bytes_per_chip": (8, 4, 16),
    # OT: 7608 2S + 5592 PS modules, 2 CICs each, header per CIC per L1A.
    "ot_modules": 13200, "ot_bytes_per_module": (8, 4, 33),
    # HGCAL (ECON-D format, EventFilter/HGCalRawToDigi):
    #   ECON-D: 2x32b header + 32b CRC; eRx: 32b if empty, 64b if not;
    #   channel data 16-32 b.
    #   Counts from the CMSSW electronics map (data-Geometry-HGCalMapping
    #   V00-07-00, modulelocator_CEminus_V15p5 x 2 endcaps): 29.8k modules =
    #   ECON-Ds; 209k eRx (half-ROCs): 176k Si + 33k SiPM-on-tile.
    #   Low variant: the back-end drops empty eRx sub-packets.
    "hgc_econd": 29796, "hgc_econd_bytes": 12,
    "hgc_erx": 209116, "hgc_empty_erx_bytes": (4, 0, 4),
    "hgc_bytes_per_hit": (3.5, 3.0, 4.0),
    # The zero-suppression threshold (in-time ADC, LSB 0.067 fC) is solved for
    # so that PU140 reproduces the TDR sub-event size.
    # ETL: ~33k ETROCs (8.5M pads / 256), 40b header + 40b trailer per L1A.
    "etl_etrocs": 33000, "etl_bytes_per_etroc": (10, 0, 10), "etl_bytes_per_hit": 5,
    # BTL: TOFHIR readout; small fixed part assumed.
    "btl_fixed": (0.01 * MB, 0.0, 0.03 * MB),
    # EB: zero suppression at 3 sigma; sigma = 4.0 ADC measured in UPC (noise only).
    "eb_thr_adc": 12, "eb_fixed": (0.01 * MB, 0.0, 0.05 * MB),
    # HCAL: no zero suppression in the Phase-2 baseline -> fixed (TDR).
    "hcal_fixed": 0.33 * MB,
    # Muons: small fixed headers.
    "mu_fixed": (0.01 * MB, 0.0, 0.03 * MB),
    # L1 trigger readout is pileup independent (TDR 0.26 MB); Track-finder TPG
    # 0.01 MB; HGCAL TPG 0.15 MB at PU140 scaled with HGCAL occupancy.
    "l1_fixed": 0.27 * MB, "hgctpg_pu200": 0.20 * MB,
}

THR_HGC = [8, 12, 16, 20, 27, 40, 54, 80, 108]


def interp_log(thr_list, counts, t):
    """Counts above threshold t, log-linear interpolation in the scan."""
    if t <= thr_list[0]:
        return counts[0]
    for i in range(1, len(thr_list)):
        if t <= thr_list[i]:
            a, b = counts[i - 1], counts[i]
            f = (t - thr_list[i - 1]) / (thr_list[i] - thr_list[i - 1])
            if a > 0 and b > 0:
                return a * (b / a) ** f
            return a + (b - a) * f
    return counts[-1]


def load(name):
    p = os.path.join(DATA, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def pick(v, var):
    """Parameter value for variant var in {0: nominal, 1: low, 2: high}."""
    return v[var] if isinstance(v, tuple) else v


def hgc_hits(c, thr):
    return sum(interp_log(THR_HGC, c[k]["above"], thr) for k in ("hgcEE", "hgcHEF", "hgcHEB") if k in c)


def eb_mean_above(ebrecs, thr):
    if not ebrecs:
        return None
    t = ebrecs[0]["thr"]
    return st.mean(interp_log(t, e["above"], thr) for e in ebrecs)


class Model:
    """Per-hit coefficients are calibrated once on (scaled) PU140, then frozen."""

    def __init__(self, pu140, eb_pu140, var=0):
        self.var = var
        v = lambda k: pick(P[k], var)
        m = lambda f: PU_SCALE * st.mean(f(e["counts"]) for e in pu140)   # PU200-equivalent
        self.fix = {
            "Inner Tracker": P["it_chips"] * v("it_bytes_per_chip"),
            "Outer Tracker": P["ot_modules"] * v("ot_bytes_per_module"),
            "MTD BTL": v("btl_fixed"),
            "MTD ETL": P["etl_etrocs"] * v("etl_bytes_per_etroc"),
            "ECAL Barrel": v("eb_fixed"),
            "HCAL (HB+HO+HF)": P["hcal_fixed"],
            "HGCAL": P["hgc_econd"] * P["hgc_econd_bytes"] + P["hgc_erx"] * v("hgc_empty_erx_bytes"),
            "Muons (DT+CSC+GEM+RPC)": v("mu_fixed"),
            "L1 trigger + TPGs": P["l1_fixed"],
        }
        # calibration to the TDR PU200 sizes
        self.k = {}
        for det, nfun in (("Inner Tracker", lambda c: c["itDigi"]),
                          ("Outer Tracker", lambda c: c["otClu"]),
                          ("MTD BTL", lambda c: c["btl"])):
            self.k[det] = (TDR_PU200[det] * MB - self.fix[det]) / m(nfun)
        n_eb = PU_SCALE * eb_mean_above(eb_pu140, P["eb_thr_adc"])
        self.k["ECAL Barrel"] = (TDR_PU200["ECAL Barrel"] * MB - self.fix["ECAL Barrel"]) / n_eb
        self.k["MTD ETL"] = (TDR_PU200["MTD ETL"] * MB - self.fix["MTD ETL"]) / m(lambda c: c["etl"])
        self.k["HGCAL"] = v("hgc_bytes_per_hit")
        self.hgc_thr = self.solve_hgc_threshold(pu140, TDR_PU200["HGCAL"] * MB)
        # muons: TDR PU200 sub-sizes scaled by digi counts, per system
        self.mu_ref = {
            "dt": (0.15 * MB, m(lambda c: c["dt"])),
            "csc": (0.47 * MB, m(lambda c: c["cscS"] + c["cscW"])),
            "gem": (0.125 * MB, m(lambda c: c["gem"])),
            "rpc": (0.011 * MB, m(lambda c: c["rpc"])),
        }
        self.hgc_pu200 = m(lambda c: hgc_hits(c, self.hgc_thr))
        self.ot_clu_per_digi = m(lambda c: c["otClu"]) / m(lambda c: c["ot_stripDigi"] + c["ot_mpxDigi"])

    def hgcal(self, c, thr, scale=1.0):
        nh = scale * hgc_hits(c, thr)
        n_erx = P["hgc_erx"]
        nonempty = n_erx * (1 - math.exp(-nh / n_erx))   # Poisson occupancy of eRx
        extra = 8 - pick(P["hgc_empty_erx_bytes"], self.var)  # a non-empty eRx header is 64 bits
        return self.fix["HGCAL"] + self.k["HGCAL"] * nh + extra * nonempty, nh

    def solve_hgc_threshold(self, pu140, target):
        mean_at = lambda t: st.mean(self.hgcal(e["counts"], t, PU_SCALE)[0] for e in pu140)
        lo, hi = THR_HGC[0], THR_HGC[-1]
        if mean_at(lo) <= target:
            return lo
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if mean_at(mid) > target else (lo, mid)
        return 0.5 * (lo + hi)

    def sizes(self, c, n_eb):
        """Bytes per subdetector for one event's counts."""
        s = {}
        s["Inner Tracker"] = self.fix["Inner Tracker"] + self.k["Inner Tracker"] * c["itDigi"]
        # some RelVals do not keep siPhase2Clusters: use digis x the PU140 cluster/digi ratio
        n_ot = c["otClu"] if "otClu" in c else (c["ot_stripDigi"] + c["ot_mpxDigi"]) * self.ot_clu_per_digi
        s["Outer Tracker"] = self.fix["Outer Tracker"] + self.k["Outer Tracker"] * n_ot
        s["MTD BTL"] = self.fix["MTD BTL"] + self.k["MTD BTL"] * c["btl"]
        s["MTD ETL"] = self.fix["MTD ETL"] + self.k["MTD ETL"] * c["etl"]
        s["ECAL Barrel"] = self.fix["ECAL Barrel"] + self.k["ECAL Barrel"] * n_eb
        s["HCAL (HB+HO+HF)"] = self.fix["HCAL (HB+HO+HF)"]
        s["HGCAL"], nh = self.hgcal(c, self.hgc_thr)
        mu = self.fix["Muons (DT+CSC+GEM+RPC)"]
        for key, n in (("dt", c["dt"]), ("csc", c["cscS"] + c["cscW"]), ("gem", c["gem"]), ("rpc", c["rpc"])):
            ref, nref = self.mu_ref[key]
            mu += ref * n / nref if nref > 0 else 0
        s["Muons (DT+CSC+GEM+RPC)"] = mu
        s["L1 trigger + TPGs"] = self.fix["L1 trigger + TPGs"] + P["hgctpg_pu200"] * nh / self.hgc_pu200
        return s


# ---------------------------------------------------------------------------
# Readout scenarios. Channel-level zero suppression is assumed in all of them
# (except HCAL, which has none in the Phase-2 baseline). They differ in two
# choices that no published design fixes:
#   drop_headers: back-ends send headers only for readout units with hits
#                 (IT chips, OT modules, ETL ETROCs, HGCAL ECON-D and eRx),
#                 instead of for every unit on every L1A;
#   eb_full:      ECAL barrel without the planned sparse readout
#                 (fixed 2.1 MB, TDR section 3.2.4).
# Per-hit costs keep the nominal calibration.
SCENARIOS = [
    ("Baseline: channel ZS, all headers kept, EB sparse readout", dict(drop_headers=False, eb_full=False)),
    ("Headers only for units with hits", dict(drop_headers=True, eb_full=False)),
    ("EB full readout (no sparsification), all headers kept", dict(drop_headers=False, eb_full=True)),
    ("EB full readout, headers only for units with hits", dict(drop_headers=True, eb_full=True)),
]
EB_FULL = 2.1 * MB
IT_CHIPS_PER_MODULE = 10000 / 3892


def occupied(n_units, n_hits):
    """Expected number of units with at least one hit (Poisson)."""
    return n_units * (1 - math.exp(-n_hits / n_units)) if n_units else 0.0


def scenario_sizes(model, c, n_eb, drop_headers, eb_full):
    s = model.sizes(c, n_eb)
    if drop_headers:
        # IT: chip headers only for chips on modules with hits (upper bound)
        it_chips = min(P["it_chips"], c["itMod"] * IT_CHIPS_PER_MODULE)
        s["Inner Tracker"] += (it_chips - P["it_chips"]) * pick(P["it_bytes_per_chip"], 0)
        # OT: module headers only for modules with hits (sensors with clusters is an upper bound)
        ot_mod = min(P["ot_modules"], c.get("otCluMod", c["ot_stripMod"] + c["ot_mpxMod"]))
        s["Outer Tracker"] += (ot_mod - P["ot_modules"]) * pick(P["ot_bytes_per_module"], 0)
        # ETL: header + trailer only for ETROCs with hits
        etl = occupied(P["etl_etrocs"], c["etl"])
        s["MTD ETL"] += (etl - P["etl_etrocs"]) * pick(P["etl_bytes_per_etroc"], 0)
        # HGCAL: ECON-D packet only for modules with hits, eRx sub-packet only if non-empty
        nh = hgc_hits(c, model.hgc_thr)
        s["HGCAL"] = (occupied(P["hgc_econd"], nh) * P["hgc_econd_bytes"]
                      + occupied(P["hgc_erx"], nh) * 8 + model.k["HGCAL"] * nh)
        # small fixed parts
        for d in ("MTD BTL", "ECAL Barrel", "Muons (DT+CSC+GEM+RPC)"):
            s[d] -= model.fix[d]
    if eb_full:
        s["ECAL Barrel"] = EB_FULL
    return s


def evaluate_scenarios(model, events, eb_recs):
    n_eb = eb_mean_above(eb_recs, P["eb_thr_adc"])
    by_key = {(e["file"], e["evt"]): interp_log(e["thr"], e["above"], P["eb_thr_adc"]) for e in eb_recs}
    out = {}
    for name, kw in SCENARIOS:
        rows = [scenario_sizes(model, e["counts"], by_key.get((e["file"], e["evt"]), n_eb), **kw) for e in events]
        out[name] = {d: st.mean(r[d] for r in rows) for d in rows[0]}
        out[name]["TOTAL"] = st.mean(sum(r.values()) for r in rows)
    return out


# ---------------------------------------------------------------------------
# Options for the size of the *stored* event (HLT output), per event class.
# Each option says, per class ("pp", "hadronic", "upc"), which readout it uses:
#   "full"    : baseline readout (all headers kept)
#   "compact" : headers only for readout units with hits (lossless for hits)
#   "nohgcal" : baseline readout with the HGCAL FEDs dropped (lossy)
#   "aggr"    : compact + HGCAL ZS at ~1 MIP (54 ADC) + 16-bit HGCAL words (lossy)
# plus an optional prescaled fraction of UPC events kept in full format.
OPTIONS = [
    ("1. Baseline (TDR-like readout)",
     dict(pp="full", hadronic="full", upc="full")),
    ("2. Back-ends drop empty headers, all events",
     dict(pp="compact", hadronic="compact", upc="compact")),
    ("3. HLT repacks UPC stream, 1% kept full [recommended]",
     dict(pp="full", hadronic="full", upc="compact", upc_full_fraction=0.01)),
    ("4. HLT drops HGCAL FEDs for UPC (lossy)",
     dict(pp="full", hadronic="full", upc="nohgcal")),
    ("5. Aggressive HGCAL, all events (lossy)",
     dict(pp="aggr", hadronic="aggr", upc="aggr")),
]


def mode_sizes(model, c, n_eb, mode):
    if mode == "full":
        return scenario_sizes(model, c, n_eb, drop_headers=False, eb_full=False)
    s = scenario_sizes(model, c, n_eb, drop_headers=True, eb_full=False)
    if mode == "nohgcal":
        s = scenario_sizes(model, c, n_eb, drop_headers=False, eb_full=False)
        s["HGCAL"] = 0.0
    elif mode == "aggr":
        nh = hgc_hits(c, 54)
        s["HGCAL"] = (occupied(P["hgc_econd"], nh) * P["hgc_econd_bytes"]
                      + occupied(P["hgc_erx"], nh) * 8 + 2.0 * nh)
    return s


def evaluate_options(model, events, eb_recs, cls):
    n_eb = eb_mean_above(eb_recs, P["eb_thr_adc"])
    by_key = {(e["file"], e["evt"]): interp_log(e["thr"], e["above"], P["eb_thr_adc"]) for e in eb_recs}
    out = {}
    for name, opt in OPTIONS:
        rows = [mode_sizes(model, e["counts"], by_key.get((e["file"], e["evt"]), n_eb), opt[cls]) for e in events]
        mean = {d: st.mean(r[d] for r in rows) for d in rows[0]}
        f = opt.get("upc_full_fraction", 0.0) if cls == "upc" else 0.0
        if f:
            full = [mode_sizes(model, e["counts"], by_key.get((e["file"], e["evt"]), n_eb), "full") for e in events]
            for d in mean:
                mean[d] = (1 - f) * mean[d] + f * st.mean(r[d] for r in full)
        mean["TOTAL"] = sum(mean.values())
        out[name] = mean
    return out


def evaluate(model, events, eb_recs):
    n_eb = eb_mean_above(eb_recs, P["eb_thr_adc"])
    by_key = {(e["file"], e["evt"]): interp_log(e["thr"], e["above"], P["eb_thr_adc"]) for e in eb_recs}
    out = []
    for e in events:
        neb = by_key.get((e["file"], e["evt"]), n_eb)
        out.append(model.sizes(e["counts"], neb))
    return out


def summarize(rows):
    dets = list(rows[0].keys())
    mean = {d: st.mean(r[d] for r in rows) for d in dets}
    tot = [sum(r.values()) for r in rows]
    return mean, tot


def scale_events(events, f):
    """Copy of events with every hit multiplicity scaled by f (occupancy ~ pileup)."""
    out = []
    for e in events:
        c = {}
        for k, v in e["counts"].items():
            if isinstance(v, dict):
                c[k] = {kk: ([f * x for x in vv] if isinstance(vv, list) else f * vv) for kk, vv in v.items()}
            elif k in ("ebU", "hbhe", "hf", "zdc"):   # channel counts, not occupancies
                c[k] = v
            else:
                c[k] = f * v
        out.append({**e, "counts": c})
    return out


def scale_eb(recs, f):
    return [{**e, "above": [f * x for x in e["above"]]} for e in recs]


def main():
    pu200 = load("pu200.jsonl")
    if pu200:   # a real PU200 sample, if one has been dumped
        pu200_entry = ("PU200 ttbar", (pu200, load("eb_pu200.jsonl")))
    else:       # otherwise the TDR's own assumption: occupancy linear in pileup
        pu200_entry = ("PU200 ttbar (PU140 x 200/140)",
                       (scale_events(load("pu140.jsonl"), PU_SCALE), scale_eb(load("eb_pu140.jsonl"), PU_SCALE)))
    samples = {
        "PU140 ttbar (validation)": (load("pu140.jsonl"), load("eb_pu140.jsonl")),
        pu200_entry[0]: pu200_entry[1],
        "Hydjet PbPb MB": (load("hydjet.jsonl"), load("eb_hydjet.jsonl")),
        "UPC (STARlight QED mumu)": (load("upc.jsonl"), load("eb_upc.jsonl")),
    }
    pu140, eb140 = samples["PU140 ttbar (validation)"]
    models = [Model(pu140, eb140, v) for v in (0, 1, 2)]
    nom = models[0]

    print("Calibrated per-hit coefficients (bytes):")
    for d, k in nom.k.items():
        print(f"  {d:24s} {k:7.2f}")
    for m in models:
        print(f"  variant {m.var}: HGCAL ZS threshold {m.hgc_thr:.1f} ADC = {m.hgc_thr * 0.0671:.2f} fC"
              f" ({m.hgc_pu200:,.0f} cells at PU200); ETL {m.k['MTD ETL']:.1f} B/hit;"
              f" OT {m.k['Outer Tracker']:.2f} B/cluster; IT {m.k['Inner Tracker']:.2f} B/digi")
    print("Fixed overheads (MB):", {d: round(f / MB, 3) for d, f in nom.fix.items()})

    result = {}
    for name, (ev, eb) in samples.items():
        rows = [evaluate(m, ev, eb) for m in models]
        (mean, tot), (mlo, tlo), (mhi, thi) = (summarize(r) for r in rows)
        result[name] = {"n": len(ev), "mean": mean, "low": mlo, "high": mhi,
                        "total": st.mean(tot), "total_low": st.mean(tlo), "total_high": st.mean(thi),
                        "median": st.median(tot), "per_event_total": tot,
                        "per_event": [{d: r[d] for d in r} for r in rows[0]],
                        "per_event_itDigi": [e["counts"]["itDigi"] for e in ev]}
        print(f"\n=== {name}  ({len(ev)} events)")
        print(f"  {'subdetector':24s} {'mean MB':>8s} {'low':>7s} {'high':>7s} {'frac':>6s}" +
              ("   TDR" if "PU" in name else ""))
        T = st.mean(tot)
        for d in mean:
            tdr = TDR_PU200 if "PU200" in name else TDR_PU140
            extra = f"   {tdr[d]:.2f}" if "PU" in name else ""
            print(f"  {d:24s} {mean[d] / MB:8.3f} {mlo[d] / MB:7.3f} {mhi[d] / MB:7.3f} {mean[d] / T:6.1%}{extra}")
        print(f"  {'TOTAL':24s} {T / MB:8.3f} {st.mean(tlo) / MB:7.3f} {st.mean(thi) / MB:7.3f}" +
              (f"          {sum((TDR_PU200 if 'PU200' in name else TDR_PU140).values()):.2f}" if "PU" in name else ""))
        print(f"  median {st.median(tot) / MB:.3f} MB, max {max(tot) / MB:.3f} MB")

    print("\n=== Readout scenarios: mean total MB (channel-level ZS in all; HCAL never zero-suppressed)")
    scen = {name: evaluate_scenarios(nom, ev, eb) for name, (ev, eb) in samples.items()}
    names = list(samples)
    print(f"  {'scenario':58s}" + "".join(f"{n[:14]:>16s}" for n in names))
    for sc, _ in SCENARIOS:
        print(f"  {sc:58s}" + "".join(f"{scen[n][sc]['TOTAL'] / MB:16.2f}" for n in names))

    classes = {"pp": next(n for n in samples if n.startswith("PU200")),
               "hadronic": "Hydjet PbPb MB", "upc": "UPC (STARlight QED mumu)"}
    opts = {cls: evaluate_options(nom, *samples[name], cls) for cls, name in classes.items()}
    print("\n=== Options: mean stored event size, MB")
    print(f"  {'option':56s}{'pp PU200':>10s}{'PbPb had':>10s}{'PbPb UPC':>10s}")
    for o, _ in OPTIONS:
        print(f"  {o:56s}" + "".join(f"{opts[c][o]['TOTAL'] / MB:10.2f}" for c in ("pp", "hadronic", "upc")))

    json.dump({"params": {k: v for k, v in P.items()}, "tdr_pu140": TDR_PU140, "tdr_pu200": TDR_PU200,
               "scenarios": scen, "options": opts, "option_samples": classes,
               "coeff": nom.k, "fixed": nom.fix, "hgc_thr_adc": [m.hgc_thr for m in models],
               "results": result},
              open(os.path.join(DATA, "size_model_results.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
