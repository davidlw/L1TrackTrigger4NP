#!/usr/bin/env python3
"""GNN track building vs the hybrid L1 tracks already in the ntuple, same events.

Efficiency, on the GNN test events, for
  GNN          a candidate with >= N stubs (or distinct layers, --count layers)
               and >= 2/3 of its stubs from the particle, for N = 3 and 4;
  no classifier  the same track building with every graph edge kept, N = 4:
               what the graph construction alone gives;
  + 1 crossing/layer  both of the above with the constrained union-find of
               evaluate.py (--crossing-dist, --max-layers): edges added by
               decreasing GNN weight (edges above --split-threshold), or by
               increasing Delta R without the classifier; a union that would
               give two crossings on a layer or > 8 layers is refused;
  hybrid       the trk_* tracks of the same ntuple entries (>= 4 stubs, one
               per layer; matched when loosely genuine, i.e. at most one stub
               from elsewhere).
Denominator (--denominator):
  stubs        particles with >= 3 stubs (or layers) in the graph; hybrid found
               = a loosely genuine track whose matched TP is the particle;
  tp           every TP in tp_* with |eta| < 2.4, the definition of
               compare_eff.C; hybrid found = tp_nmatch > 0. A TP without stubs
               in the graph cannot be found by the GNN.
Panels: (a) vs pT, (b) vs eta and (c) vs phi above --eta-pt-min, and (d)
either the hybrid efficiency of a --reference sample (standard definition, own
TPs) or, without one, the fraction of particles with stubs on >= 3 / >= 4 layers
and the GNN efficiency with an ideal edge classifier.

--offline-ntuple adds 2D efficiency maps (eta vs pT) from the OfflineTrackNtuple
files themselves, with the offline_perf.C definition: primary TPs (tp_ngenpart
> 0) with |eta| < 2.4 and |pdgId| = --offline-pdg, found if matched by at least
one highPurity track.

--offline overlays, as a line, the offline (generalTracks) efficiency from an
output of ../../offline_tracking_performance/scripts/offline_perf.C. Its eta and
phi histograms exist for pT > 0.3, 0.6, 1 and 2 GeV; --eta-pt-min must be one of
those for them to be drawn. The offline ntuples are not the same events.

Usage:
  python3 compare_hybrid.py --checkpoint ../output/model_ppb/model.pt --threshold 0.9 \
      --reference /my/ntuples/DefaultStub/EPOS_pPb
  python3 compare_hybrid.py --checkpoint ../output/model_mumu/model.pt --count layers --denominator tp
"""
import argparse
import json
import os
from collections import defaultdict

import numpy as np
import torch
import uproot

from gnn_common import TREE, list_root_files, load_shards, read_json, split_events
from evaluate import binned_eff, build_and_match, crossing_ids, particle_nlayers
from gnn_common import EDGE_FEATURES
from model import InteractionNetwork
from train import GraphSet, predict

TRK = ["trk_matchtp_pdgid", "trk_matchtp_pt", "trk_matchtp_eta", "trk_matchtp_phi", "trk_loose"]
TP = ["tp_pdgid", "tp_pt", "tp_eta", "tp_phi", "tp_nmatch"]
PT_BINS = np.array([0.3, 0.4, 0.5, 0.6, 0.75, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 20])
ETA_BINS = np.linspace(-2.5, 2.5, 21)
PHI_BINS = np.linspace(-np.pi, np.pi, 17)
OFFLINE_THRESHOLDS = [0.3, 0.6, 1.0, 2.0]  # index of den_eta<k>_off in offline_perf.C
MAX_ETA = 2.4


def read_entries(files, entries_by_file, branches):
    """{(file, entry): {branch: array}} for the requested entries."""
    out = {}
    for ifile, entries in entries_by_file.items():
        with uproot.open(files[ifile]) as rf:
            a = rf[TREE].arrays(branches, library="np")
        for e in entries:
            out[(ifile, e)] = {b: a[b][e] for b in branches}
    return out


def keys_of(pdg, pt, eta, phi):
    return list(zip(np.asarray(pdg, dtype=float), pt, eta, phi))


def offline_curves(path, pt_min):
    """{variable: (x, efficiency)} from an offline_perf.C output; eta/phi only if pt_min is a threshold there."""
    out = {}
    with uproot.open(path) as rf:
        names = {"pt": "pt"}
        if pt_min in OFFLINE_THRESHOLDS:
            k = OFFLINE_THRESHOLDS.index(pt_min)
            names.update({"eta": f"eta{k}", "phi": f"phi{k}"})
        else:
            print(f"[offline] no eta/phi histograms for pT > {pt_min}; choose one of {OFFLINE_THRESHOLDS}")
        for var, h in names.items():
            num, den = rf[f"num_{h}_off"], rf[f"den_{h}_off"]
            edges = den.axis().edges()
            n, d = num.values(), den.values()
            x = np.sqrt(edges[1:] * edges[:-1]) if var == "pt" else 0.5 * (edges[1:] + edges[:-1])
            ok = d > 0
            out[var] = (x[ok], n[ok] / d[ok])
    return out


def standard_tp_eff(files):
    """(pt, found) of every TP with |eta| < MAX_ETA, found = tp_nmatch > 0 (compare_eff.C)."""
    pts, found = [], []
    for f in files:
        with uproot.open(f) as rf:
            t = rf[TREE]
            if "tp_pt" not in t.keys():
                print(f"[warn] {f} has no tp_* branches")
                continue
            a = t.arrays(["tp_pt", "tp_eta", "tp_nmatch"], library="np")
        for pt, eta, nm in zip(a["tp_pt"], a["tp_eta"], a["tp_nmatch"]):
            m = np.abs(eta) < MAX_ETA
            pts.append(pt[m])
            found.append(nm[m] > 0)
    return np.concatenate(pts), np.concatenate(found).astype(float)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--graphs", help="graph directory to evaluate (default: the training graphs, test split; "
                                     "another directory is evaluated on all its events)")
    ap.add_argument("--threshold", type=float, help="GNN edge-weight cut (default: delta from training)")
    ap.add_argument("--min-purity", type=float, default=2.0 / 3.0)
    ap.add_argument("--count", choices=["stubs", "layers"], default="stubs",
                    help="GNN candidate size counted in stubs or in distinct layers")
    ap.add_argument("--denominator", choices=["stubs", "tp"], default="stubs")
    ap.add_argument("--eta-pt-min", type=float, default=2.0, help="pT cut for the efficiency vs eta panel")
    ap.add_argument("--reference", help="ntuples for panel (c), e.g. the DefaultStub sample")
    ap.add_argument("--reference-label", default="DefaultStub")
    ap.add_argument("--crossing-dist", type=float, default=5.0,
                    help="stubs on a layer within this distance [cm] are one crossing (0 = no split curves)")
    ap.add_argument("--max-layers", type=int, default=8)
    ap.add_argument("--split-threshold", type=float, default=0.5,
                    help="GNN edge cut for the one-crossing-per-layer track building")
    ap.add_argument("--offline", help="offline_perf.C output to overlay as a line")
    ap.add_argument("--offline-label", default="offline generalTracks (highPurity)")
    ap.add_argument("--offline-ntuple", help="OfflineTrackNtuple file or directory, for the 2D maps")
    ap.add_argument("--offline-pdg", type=int, default=13, help="|pdgId| of the offline truth particles (0 = all)")
    ap.add_argument("--nfiles", type=int, default=2, help="files per sample in panel (c)")
    ap.add_argument("--out", help="output directory (default: next to the checkpoint)")
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    ck = torch.load(args.checkpoint, map_location=args.device, weights_only=False)
    graphs = os.path.abspath(args.graphs) if args.graphs else ck["graphs"]
    cfg = read_json(os.path.join(graphs, "graph_config.json"))
    if "files" not in cfg:
        raise SystemExit("graphs were built before event bookkeeping was added; rerun build_graphs.py")
    files = cfg["files"]
    events = split_events(load_shards(graphs))
    idx = ck["test_idx"] if graphs == ck["graphs"] else list(range(len(events)))
    out = args.out or os.path.join(os.path.dirname(args.checkpoint), f"compare_hybrid_{args.denominator}_{args.count}")
    os.makedirs(out, exist_ok=True)
    thr = args.threshold if args.threshold is not None else ck["delta"]
    unit = args.count

    norm = {k: np.array(v, dtype=np.float32) for k, v in ck["norm"].items()}
    data = GraphSet(events, ck["node_cols"], ck["edge_cols"], norm)
    model = InteractionNetwork(len(ck["node_cols"]), len(ck["edge_cols"]), hidden=ck["args"]["hidden"],
                               iterations=ck["args"]["iterations"]).to(args.device)
    model.load_state_dict(ck["model"])
    s, _ = predict(model, data, idx, args.batch_size, torch.device(args.device), ck["signed_cols"])
    scores = np.split(s, np.cumsum([events[i]["edge_index"].shape[1] for i in idx])[:-1])

    by_file = defaultdict(list)
    for i in idx:
        by_file[events[i]["file"]].append(events[i]["entry"])
    ntuple = read_entries(files, by_file, TRK + (TP if args.denominator == "tp" else []))

    # columns: pt, eta, GNN >= 3, GNN >= 4, ideal >= 3, hybrid, n (stubs or layers), no classifier >= 4, phi,
    #          GNN + split >= 4, no classifier + split >= 4, GNN + split >= 3
    dR = EDGE_FEATURES.index("dR")
    do_split = args.crossing_dist > 0
    rows, n_cand, n_true, n_trk, n_loose = [], 0, 0, 0, 0
    fk = defaultdict(lambda: [0, 0, 0])  # variant -> [candidates, true tracks, duplicates]

    def tally(name, res):
        fk[name][0] += res[0]
        fk[name][1] += res[1]
        fk[name][2] += res[2]
        return res[3]

    for k, i in enumerate(idx):
        ev = events[i]
        layers = ev["layer"] if unit == "layers" else None
        common = (len(ev["pid"]), ev["edge_index"])
        rest = (ev["pid"], ev["part_nstub"])
        nc, nt, _, f3 = build_and_match(*common, scores[k] > thr, *rest, 3, args.min_purity, layers)
        f4 = tally("GNN >= 4", build_and_match(*common, scores[k] > thr, *rest, 4, args.min_purity, layers))
        _, _, _, fi = build_and_match(*common, ev["y"] > 0.5, *rest, 3, args.min_purity, layers)
        allk = np.ones(len(ev["y"]), bool)
        fa = tally("no classifier >= 4", build_and_match(*common, allk, *rest, 4, args.min_purity, layers))
        fs = fas = fs3 = np.zeros(len(ev["part_nstub"]), bool)
        if do_split and len(ev["pid"]):
            cr = crossing_ids(ev["x"], ev["layer"], args.crossing_dist)
            fs = tally("GNN + 1 crossing/layer >= 4",
                       build_and_match(*common, scores[k] > args.split_threshold, *rest, 4, args.min_purity,
                                       layers, (scores[k], ev["layer"], cr, args.max_layers)))
            fs3 = tally("GNN + 1 crossing/layer >= 3",
                        build_and_match(*common, scores[k] > args.split_threshold, *rest, 3, args.min_purity,
                                        layers, (scores[k], ev["layer"], cr, args.max_layers)))
            fas = tally("no classifier + 1 crossing/layer >= 4",
                        build_and_match(*common, allk, *rest, 4, args.min_purity, layers,
                                        (-ev["edge_attr"][:, dR], ev["layer"], cr, args.max_layers)))
        n_cand += nc
        n_true += nt
        npart = particle_nlayers(ev) if unit == "layers" else ev["part_nstub"]
        nt_ = ntuple[(ev["file"], ev["entry"])]
        loose = nt_["trk_loose"] == 1
        n_trk += len(loose)
        n_loose += int(loose.sum())
        pkeys = keys_of(ev["part_pdg"], ev["part_pt"], ev["part_eta"], ev["part_phi"])

        if args.denominator == "stubs":
            hk = set(keys_of(*(nt_[b][loose] for b in TRK[:4])))
            fh = np.array([key in hk for key in pkeys], dtype=bool)
            sel = npart >= 3
            rows.append(np.stack([ev["part_pt"][sel], ev["part_eta"][sel], f3[sel], f4[sel], fi[sel],
                                  fh[sel], npart[sel], fa[sel], ev["part_phi"][sel], fs[sel], fas[sel], fs3[sel]], axis=1))
        else:
            where = {key: j for j, key in enumerate(pkeys)}
            m = np.abs(nt_["tp_eta"]) < MAX_ETA
            for key, pt, eta, phi, nm in zip(keys_of(*(nt_[b][m] for b in TP[:4])), nt_["tp_pt"][m],
                                             nt_["tp_eta"][m], nt_["tp_phi"][m], nt_["tp_nmatch"][m]):
                j = where.get(key)
                g = ((f3[j], f4[j], fi[j], npart[j], fa[j], fs[j], fas[j], fs3[j]) if j is not None
                     else (False, False, False, 0, False, False, False, False))
                rows.append(np.array([[pt, eta, g[0], g[1], g[2], nm > 0, g[3], g[4], phi, g[5], g[6], g[7]]],
                                     dtype=float))
    rows = np.concatenate(rows).astype(float)
    pt, eta, phi = rows[:, 0], rows[:, 1], rows[:, 8]

    summary = {"threshold": thr, "count": unit, "denominator": args.denominator, "events": len(idx),
               "particles": len(rows), "gnn_fake_rate": (n_cand - n_true) / max(n_cand, 1),
               "hybrid_fake_rate": (n_trk - n_loose) / max(n_trk, 1), "efficiency": {}}
    edges = [0.3, 0.5, 0.75, 1, 2, np.inf]
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (pt >= lo) & (pt < hi)
        if not m.any():
            continue
        name = f"pT {lo}-{hi}" if np.isfinite(hi) else f"pT > {lo}"
        e = {"particles": int(m.sum()), "gnn_3": float(rows[m, 2].mean()), "gnn_4": float(rows[m, 3].mean()),
             "ideal_3": float(rows[m, 4].mean()), "no_classifier_4": float(rows[m, 7].mean()),
             "gnn_split_4": float(rows[m, 9].mean()), "no_classifier_split_4": float(rows[m, 10].mean()),
             "hybrid": float(rows[m, 5].mean()),
             f"frac_ge3_{unit}": float((rows[m, 6] >= 3).mean()), f"frac_ge4_{unit}": float((rows[m, 6] >= 4).mean())}
        summary["efficiency"][name] = e
        print(f"[{name:10s}] {e['particles']:6d}  GNN>=3 {e['gnn_3']:.3f}  GNN>=4 {e['gnn_4']:.3f}  "
              f"ideal>=3 {e['ideal_3']:.3f}  no-clf>=4 {e['no_classifier_4']:.3f}  "
              f"GNN+split>=4 {e['gnn_split_4']:.3f}  no-clf+split>=4 {e['no_classifier_split_4']:.3f}  "
              f"hybrid {e['hybrid']:.3f}   "
              f"(>=3 {unit}: {e[f'frac_ge3_{unit}']:.3f}, >=4: {e[f'frac_ge4_{unit}']:.3f})")
    summary["variants"] = {}
    for name, (c, t, d) in fk.items():
        summary["variants"][name] = {"fake_rate": (c - t) / max(c, 1), "duplicate_rate": d / max(t, 1)}
        print(f"[fakes] {name:40s} fake {(c - t) / max(c, 1):.4f}  duplicate {d / max(t, 1):.4f}  ({c} candidates)")
    print(f"[fakes] {'hybrid':40s} fake {summary['hybrid_fake_rate']:.4f}")

    ref_curves = []
    if args.reference:
        own = files[: args.nfiles]
        ref = list_root_files(args.reference, args.nfiles)
        for lab, fl in [("hybrid, graph sample", own), (f"hybrid, {args.reference_label}", ref)]:
            p, f = standard_tp_eff(fl)
            ref_curves.append((lab, p, f))
            hi = p >= 2
            summary[f"standard_eff_pT>2: {lab}"] = float(f[hi].mean()) if hi.any() else None
            print(f"[standard] {lab}: TP efficiency for pT > 2 GeV = {f[hi].mean():.3f} ({hi.sum()} TPs)")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    ax = [axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]]
    xpt = np.sqrt(PT_BINS[1:] * PT_BINS[:-1])
    xeta = 0.5 * (ETA_BINS[1:] + ETA_BINS[:-1])
    xphi = 0.5 * (PHI_BINS[1:] + PHI_BINS[:-1])
    curves = [(f"GNN, >= 4 {unit} (cut {thr:.2f})", 3, "^", "C2"),
              (f"no classifier (all edges), >= 4 {unit}", 7, "x", "C7")]
    if do_split:
        curves += [(f"GNN + 1 crossing/layer, >= 4 (cut {args.split_threshold:.2f})", 9, "o", "C0"),
                   ("no classifier + 1 crossing/layer, >= 4", 10, "D", "C1")]
    curves += [("hybrid L1 tracks", 5, "s", "C3")]
    hi = pt >= args.eta_pt_min
    for lab, col, mk, c in curves:
        e, err = binned_eff(pt, rows[:, col], PT_BINS)
        ax[0].errorbar(xpt, e, yerr=err, fmt=mk, ms=4, color=c, label=lab)
        e, err = binned_eff(eta[hi], rows[hi, col], ETA_BINS)
        ax[1].errorbar(xeta, e, yerr=err, fmt=mk, ms=4, color=c, label=lab)
        e, err = binned_eff(phi[hi], rows[hi, col], PHI_BINS)
        ax[2].errorbar(xphi, e, yerr=err, fmt=mk, ms=4, color=c, label=lab)
    if args.offline:
        off = offline_curves(args.offline, args.eta_pt_min)
        for a, var in zip(ax[:3], ["pt", "eta", "phi"]):
            if var in off:
                a.plot(*off[var], "-", color="k", lw=1.5, label=args.offline_label, zorder=0)
        summary["offline"] = {v: {"x": x.tolist(), "eff": e.tolist()} for v, (x, e) in off.items()}
    den = "all TPs, |eta| < 2.4" if args.denominator == "tp" else f"particles with >= 3 {unit}"
    ax[0].set_xscale("log"); ax[0].set_xlabel(r"particle $p_T$ [GeV]")
    ax[0].set_title(f"(a) same events; denominator: {den}", fontsize=9)
    ax[1].set_xlabel(r"particle $\eta$")
    ax[1].set_title(rf"(b) same, $p_T > {args.eta_pt_min}$ GeV", fontsize=9)
    ax[2].set_xlabel(r"particle $\phi$")
    ax[2].set_title(rf"(c) same, $p_T > {args.eta_pt_min}$ GeV", fontsize=9)
    for a in ax[:3]:
        a.set_ylabel("tracking efficiency"); a.set_ylim(0, 1.05); a.legend(fontsize=8, loc="lower right")
    if ref_curves:
        for (lab, p, f), mk, c in zip(ref_curves, ["s", "D"], ["C3", "C1"]):
            e, err = binned_eff(p, f, PT_BINS)
            ax[3].errorbar(xpt, e, yerr=err, fmt=mk, ms=4, color=c, label=lab)
        ax[3].set_ylabel(r"TP efficiency ($n_{match}>0$ / TPs)")
        ax[3].set_title("(d) hybrid, standard definition, own TPs", fontsize=9)
    else:
        for lab, vals, mk, c in [(f"stubs on >= 3 layers", rows[:, 6] >= 3, "v", "gray"),
                                 (f"stubs on >= 4 layers", rows[:, 6] >= 4, "^", "k"),
                                 ("GNN, ideal edge classifier, >= 3", rows[:, 4], "o", "C0")]:
            if unit == "stubs" and "layers" in lab:
                lab = lab.replace("layers", "stubs").replace("stubs on", "")
            e, err = binned_eff(pt, vals.astype(float), PT_BINS)
            ax[3].errorbar(xpt, e, yerr=err, fmt=mk, ms=4, color=c, mfc="none" if c == "C0" else c, label=lab)
        ax[3].set_ylabel("fraction of denominator")
        ax[3].set_title("(d) acceptance and graph-construction limit", fontsize=9)
    ax[3].set_xscale("log"); ax[3].set_ylim(0, 1.05); ax[3].set_xlabel(r"particle $p_T$ [GeV]")
    ax[3].legend(fontsize=8, loc="lower right")
    for a in (ax[0], ax[3]):
        a.set_xticks([0.3, 0.5, 1, 2, 5, 10]); a.set_xticklabels(["0.3", "0.5", "1", "2", "5", "10"])
    with open(os.path.join(out, "summary.json"), "w") as fj:
        json.dump(summary, fj, indent=1)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "gnn_vs_hybrid.png"), dpi=130)
    fig.savefig(os.path.join(out, "gnn_vs_hybrid.pdf"))
    np.savez(os.path.join(out, "rows.npz"), rows=rows,
             columns=np.array(["pt", "eta", "gnn_3", "gnn_4", "ideal_3", "hybrid", "n", "no_classifier_4", "phi",
                               "gnn_split_4", "no_classifier_split_4", "gnn_split_3"]))
    gnn_map = ("GNN + 1 crossing/layer, >= 4", 9) if do_split else (f"GNN, >= 4 {unit}", 3)
    maps = [("hybrid L1 tracks", rows[:, 0], rows[:, 1], rows[:, 5]),
            (gnn_map[0], rows[:, 0], rows[:, 1], rows[:, gnn_map[1]])]
    if args.offline_ntuple:
        o = offline_rows(list_root_files(args.offline_ntuple), args.offline_pdg)
        maps.insert(0, (args.offline_label, o[:, 0], o[:, 1], o[:, 2]))
        print(f"[offline] {len(o)} truth particles, efficiency {o[:, 2].mean():.3f}")
    plot_2d(maps, out)
    print(f"[done] {out}/gnn_vs_hybrid.png, eff2d.png, summary.json, rows.npz")


def offline_rows(files, pdg):
    """(pt, eta, found, phi) per primary TP in OfflineTrackNtuple files (offline_perf.C definition)."""
    out = []
    br = ["tp_pt", "tp_eta", "tp_phi", "tp_pdgid", "tp_ngenpart", "trk_highPurity", "trk_matchtp_pt",
          "trk_matchtp_eta", "trk_matchtp_phi"]
    for f in files:
        with uproot.open(f) as rf:
            a = rf["OfflineTrackNtupleMaker/eventTree"].arrays(br, library="np")
        for i in range(len(a["tp_pt"])):
            hp = a["trk_highPurity"][i] == 1
            matched = set(zip(a["trk_matchtp_pt"][i][hp], a["trk_matchtp_eta"][i][hp], a["trk_matchtp_phi"][i][hp]))
            sel = (np.abs(a["tp_eta"][i]) < MAX_ETA) & (a["tp_ngenpart"][i] > 0)
            if pdg:
                sel &= np.abs(a["tp_pdgid"][i]) == pdg
            for pt, eta, phi in zip(a["tp_pt"][i][sel], a["tp_eta"][i][sel], a["tp_phi"][i][sel]):
                out.append((pt, eta, (pt, eta, phi) in matched, phi))
    return np.array(out, dtype=float)


PT2 = np.array([0.3, 0.4, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0])
ETA2 = np.linspace(-2.4, 2.4, 13)


def eff_map(pt, eta, found, min_count=10):
    den, _, _ = np.histogram2d(eta, pt, [ETA2, PT2])
    num, _, _ = np.histogram2d(eta, pt, [ETA2, PT2], weights=found)
    with np.errstate(invalid="ignore", divide="ignore"):
        e = num / den
    e[den < min_count] = np.nan
    return e, den


def plot_2d(maps, out):
    """Efficiency vs (eta, pT) for each (label, pt, eta, found); second row: differences.

    With two maps (hybrid, GNN): one difference, GNN - hybrid. With three
    (offline, hybrid, GNN): hybrid - offline, GNN - offline, GNN - hybrid.
    """
    import matplotlib.pyplot as plt

    effs = [(lab, eff_map(pt, eta, f)[0]) for lab, pt, eta, f in maps]
    pairs = [(1, 0)] if len(effs) == 2 else [(1, 0), (2, 0), (2, 1)]
    ncol = len(effs)
    fig, axs = plt.subplots(2, ncol, figsize=(5.6 * ncol, 9.2))
    ax = list(axs[0]) + list(axs[1])
    X, Y = np.meshgrid(np.arange(len(PT2)), ETA2)
    xt = [f"{v:g}" for v in PT2]

    def draw(a, m, title, cmap, vmin, vmax, fmt):
        pc = a.pcolormesh(X, Y, m, cmap=cmap, vmin=vmin, vmax=vmax)
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                if np.isfinite(m[i, j]):
                    a.text(j + 0.5, 0.5 * (ETA2[i] + ETA2[i + 1]), fmt.format(m[i, j]), ha="center", va="center",
                           fontsize=6.5, color="k" if (m[i, j] - vmin) / (vmax - vmin) > 0.35 else "w")
        a.set_xticks(np.arange(len(PT2))); a.set_xticklabels(xt, fontsize=7)
        a.set_xlabel(r"particle $p_T$ [GeV]"); a.set_ylabel(r"particle $\eta$"); a.set_title(title, fontsize=9)
        fig.colorbar(pc, ax=a, fraction=0.046)

    short = lambda lab: lab.split(",")[0].split(" (")[0]
    for a, (lab, m) in zip(ax, effs):
        draw(a, m, f"efficiency: {lab}", "viridis", 0, 1, "{:.2f}")
    for a, (i, j) in zip(ax[ncol:], pairs):
        draw(a, effs[i][1] - effs[j][1], f"{short(effs[i][0])} - {short(effs[j][0])}", "RdBu", -1, 1, "{:+.2f}")
    for a in ax[ncol + len(pairs):]:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "eff2d.png"), dpi=130)
    fig.savefig(os.path.join(out, "eff2d.pdf"))


if __name__ == "__main__":
    main()
