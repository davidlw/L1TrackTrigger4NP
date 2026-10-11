#!/usr/bin/env python3
"""Edge-classification and track-building performance of a trained model.

Track building: keep the edges with weight > delta, take connected components
of the kept edges (union-find, as in DeZoort et al. Sec. 4.4), and call each
component with at least --min-nodes stubs a track candidate.

With --crossing-dist D > 0 the union-find is constrained: a track that does not
loop crosses each layer once and spans at most --max-layers (8) layers. Stubs on
one layer within D cm of each other (single linkage) form one crossing (DummyStub
gives several stubs per crossing). Kept edges are added in order of decreasing
weight, and an edge is refused if the merged candidate would then have two
different crossings on a layer, or more than --max-layers layers. Merged
tracks are split this way instead of lost.

Matching: a candidate is a true track if at least --min-purity of its stubs come
from one particle; that particle is then found. Efficiency is the fraction of
particles with >= --min-nodes stubs in the graph that are found; the fake rate
is the fraction of candidates that are not true tracks; the duplicate rate is the
fraction of true tracks whose particle was already found by another candidate.

The same track building is also run with the truth labels as edge weights
("ideal classifier"), which isolates the losses due to graph construction.

Usage:
  python3 evaluate.py --checkpoint ../output/model_ppb/model.pt              # its test split
  python3 evaluate.py --checkpoint ../output/model_ppb/model.pt --graphs ../output/graphs_hydjet
"""
import argparse
import json
import os

import numpy as np
import torch
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from gnn_common import load_shards, split_events
from model import InteractionNetwork
from train import GraphSet, auc, predict, roc_curve

SPECIES = {"electron": [11], "muon": [13], "hadron": None}


def species_of(pdg):
    a = np.abs(pdg)
    return np.where(a == 11, "electron", np.where(a == 13, "muon", "hadron"))


def particle_nlayers(ev):
    """Number of distinct layers with a stub of each particle of the event."""
    m = ev["pid"] >= 0
    pairs = np.unique(np.stack([ev["pid"][m], ev["layer"][m]]), axis=1)
    return np.bincount(pairs[0], minlength=len(ev["part_nstub"]))


def crossing_ids(x, layer, dist):
    """Label the stubs of each layer that belong to one crossing: single linkage within dist [cm]."""
    r, phi, z = x[:, 0], x[:, 1], x[:, 2]
    pos = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)
    ids = np.full(len(layer), -1, dtype=np.int64)
    nxt = 0
    for l in np.unique(layer):
        m = np.where(layer == l)[0]
        pairs = cKDTree(pos[m]).query_pairs(dist, output_type="ndarray")
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(m), len(m)))
        nc, lab = connected_components(g, directed=False)
        ids[m] = lab + nxt
        nxt += nc
    return ids


def constrained_components(n_nodes, edge_index, order, layer, crossing, max_layers):
    """Union-find adding edges in the given order, refusing any union that would put two
    crossings on one layer or more than max_layers layers into a component."""
    parent = np.arange(n_nodes)
    comp = {i: {int(layer[i]): int(crossing[i])} for i in range(n_nodes)}  # root -> {layer: crossing}

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for e in order:
        ra, rb = find(edge_index[0, e]), find(edge_index[1, e])
        if ra == rb:
            continue
        da, db = comp[ra], comp[rb]
        if len(da) < len(db):
            ra, rb, da, db = rb, ra, db, da
        if any(da.get(l, c) != c for l, c in db.items()):
            continue
        if len(da) + sum(l not in da for l in db) > max_layers:
            continue
        parent[rb] = ra
        da.update(db)
        del comp[rb]
    return np.array([find(i) for i in range(n_nodes)])


def build_and_match(n_nodes, edge_index, keep, pid, part_nstub, min_nodes, min_purity, layers=None,
                    split=None):
    """Connected components of the kept edges, matched to particles.

    A component is a candidate if it has >= min_nodes stubs, or, when `layers`
    (the layer of each node) is given, stubs on >= min_nodes distinct layers.
    split = (order_score, node_layer, crossing, max_layers) switches to the
    constrained union-find, adding kept edges by decreasing order_score.
    Returns (n_candidates, n_true_tracks, n_duplicates, found[particle]).
    """
    found = np.zeros(len(part_nstub), dtype=bool)
    if n_nodes == 0:
        return 0, 0, 0, found
    if split is None:
        ei = edge_index[:, keep]
        g = coo_matrix((np.ones(ei.shape[1]), (ei[0], ei[1])), shape=(n_nodes, n_nodes))
        _, label = connected_components(g, directed=False)
    else:
        score, node_layer, crossing, max_layers = split
        kept = np.where(keep)[0]
        order = kept[np.argsort(-score[kept], kind="stable")]
        label = constrained_components(n_nodes, edge_index, order, node_layer, crossing, max_layers)
        label = np.unique(label, return_inverse=True)[1]
    size = np.bincount(label)
    count = size
    if layers is not None:
        count = np.bincount(np.unique(np.stack([label, layers]), axis=1)[0], minlength=len(size))
    ncand = ntrue = ndup = 0
    for c in np.where(count >= min_nodes)[0]:
        ncand += 1
        p = pid[label == c]
        p = p[p >= 0]
        if len(p) == 0:
            continue
        counts = np.bincount(p)
        major = counts.argmax()
        if counts[major] >= min_purity * size[c] - 1e-9:
            ntrue += 1
            if found[major]:
                ndup += 1
            found[major] = True
    return ncand, ntrue, ndup, found


def split_args(ev, score, crossing_dist, max_layers):
    """The split= argument of build_and_match for one event (None if the constraint is off)."""
    if crossing_dist <= 0 or len(ev["pid"]) == 0:
        return None
    return score, ev["layer"], crossing_ids(ev["x"], ev["layer"], crossing_dist), max_layers


def track_performance(events, idx, scores_per_event, threshold, min_nodes, min_purity, count="stubs",
                      crossing_dist=0.0, max_layers=8):
    tot = {"candidates": 0, "true_tracks": 0, "duplicates": 0}
    rows = []  # (pt, eta, species, found) for every reconstructable particle
    for k, i in enumerate(idx):
        ev = events[i]
        keep = scores_per_event[k] > threshold
        layers = ev["layer"] if count == "layers" else None
        split = split_args(ev, scores_per_event[k], crossing_dist, max_layers)
        nc, nt, nd, found = build_and_match(len(ev["pid"]), ev["edge_index"], keep, ev["pid"],
                                            ev["part_nstub"], min_nodes, min_purity, layers, split)
        tot["candidates"] += nc
        tot["true_tracks"] += nt
        tot["duplicates"] += nd
        sel = (particle_nlayers(ev) if count == "layers" else ev["part_nstub"]) >= min_nodes
        rows.append(np.stack([ev["part_pt"][sel], ev["part_eta"][sel],
                              np.abs(ev["part_pdg"][sel]), found[sel]], axis=1))
    rows = np.concatenate(rows) if rows else np.zeros((0, 4))
    res = {
        "threshold": float(threshold),
        "particles": int(len(rows)),
        "efficiency": float(rows[:, 3].mean()) if len(rows) else float("nan"),
        "fake_rate": (tot["candidates"] - tot["true_tracks"]) / max(tot["candidates"], 1),
        "duplicate_rate": tot["duplicates"] / max(tot["true_tracks"], 1),
        **tot,
    }
    sp = species_of(rows[:, 2].astype(int))
    for name in SPECIES:
        m = sp == name
        res[f"efficiency_{name}"] = float(rows[m, 3].mean()) if m.any() else None
        res[f"particles_{name}"] = int(m.sum())
    return res, rows


def binned_eff(val, found, bins):
    idx = np.digitize(val, bins) - 1
    eff, err = np.full(len(bins) - 1, np.nan), np.zeros(len(bins) - 1)
    for b in range(len(bins) - 1):
        m = idx == b
        if m.sum() > 0:
            e = found[m].mean()
            eff[b], err[b] = e, np.sqrt(max(e * (1 - e), 1.0 / m.sum()) / m.sum())
    return eff, err


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--graphs", help="graph directory (default: the training graphs, test split)")
    ap.add_argument("--out", help="output directory (default: next to the checkpoint)")
    ap.add_argument("--threshold", type=float, help="edge-weight cut (default: delta from training)")
    ap.add_argument("--min-nodes", type=int, default=3, help="minimum stubs (or layers) in a track candidate")
    ap.add_argument("--count", choices=["stubs", "layers"], default="stubs",
                    help="apply --min-nodes to stubs, or to distinct layers (several stubs per layer "
                         "per particle are common with DummyStub)")
    ap.add_argument("--min-purity", type=float, default=2.0 / 3.0,
                    help="minimum fraction of a candidate's stubs from one particle (default 2/3)")
    ap.add_argument("--crossing-dist", type=float, default=0.0,
                    help="> 0: one crossing per layer, stubs within this distance [cm] form a crossing")
    ap.add_argument("--max-layers", type=int, default=8, help="with --crossing-dist: most layers per track")
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--batch-size", type=int, default=64)
    ap.add_argument("--device", default="cpu")
    args = ap.parse_args()

    ck = torch.load(args.checkpoint, map_location=args.device, weights_only=False)
    targs = ck["args"]
    graphs = args.graphs or ck["graphs"]
    same_sample = os.path.abspath(graphs) == ck["graphs"]
    events = split_events(load_shards(graphs))
    idx = ck["test_idx"] if same_sample else list(range(len(events)))
    if args.max_events:
        idx = idx[: args.max_events]
    out = args.out or os.path.join(os.path.dirname(args.checkpoint),
                                   "eval_" + os.path.basename(os.path.normpath(graphs)))
    os.makedirs(out, exist_ok=True)
    print(f"[data] {graphs}: {len(idx)} events ({'test split' if same_sample else 'all events'})")

    norm = {k: np.array(v, dtype=np.float32) for k, v in ck["norm"].items()}
    data = GraphSet(events, ck["node_cols"], ck["edge_cols"], norm)
    model = InteractionNetwork(len(ck["node_cols"]), len(ck["edge_cols"]),
                               hidden=targs["hidden"], iterations=targs["iterations"]).to(args.device)
    model.load_state_dict(ck["model"])
    delta = args.threshold if args.threshold is not None else ck["delta"]

    s, y = predict(model, data, idx, args.batch_size, torch.device(args.device), ck["signed_cols"])
    split_at = np.cumsum([events[i]["edge_index"].shape[1] for i in idx])[:-1]
    scores = np.split(s, split_at)
    truth = [events[i]["y"].astype(float) for i in idx]

    thr, tpr, tnr = roc_curve(s, y)
    pred = s > delta
    edge = {"edges": int(len(y)), "true_fraction": float(y.mean()), "auc": auc(tpr, tnr), "delta": delta,
            "accuracy": float((pred == y.astype(bool)).mean()),
            "tpr": float(pred[y == 1].mean()), "tnr": float(1 - pred[y == 0].mean())}
    print(f"[edges] {edge['edges']} edges, true fraction {edge['true_fraction']:.4f}, AUC {edge['auc']:.5f}")
    print(f"[edges] at delta = {delta:.4f}: accuracy {edge['accuracy']:.4f}, TPR {edge['tpr']:.4f}, "
          f"TNR {edge['tnr']:.4f}")

    tp = lambda sc, t: track_performance(events, idx, sc, t, args.min_nodes, args.min_purity, args.count,
                                         args.crossing_dist, args.max_layers)
    gnn, rows = tp(scores, delta)
    ideal, rows_ideal = tp(truth, 0.5)
    fmt = lambda r: (f"eff {r['efficiency']:.4f} (e {r['efficiency_electron'] or float('nan'):.3f}, "
                     f"mu {r['efficiency_muon'] or float('nan'):.3f}, "
                     f"had {r['efficiency_hadron'] or float('nan'):.3f})  fake {r['fake_rate']:.4f}  "
                     f"dup {r['duplicate_rate']:.4f}  [{r['particles']} particles, {r['candidates']} candidates]")
    print(f"[tracks] >= {args.min_nodes} {args.count}, purity >= {args.min_purity:.3f}")
    print(f"[tracks] GNN   {fmt(gnn)}")
    print(f"[tracks] ideal {fmt(ideal)}")

    scan = []
    for t in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
        r, _ = tp(scores, t)
        scan.append(r)
        print(f"[scan]  threshold {t:.2f}: {fmt(r)}")

    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump({"edge": edge, "gnn": gnn, "ideal": ideal, "threshold_scan": scan,
                   "min_nodes": args.min_nodes, "count": args.count, "min_purity": args.min_purity,
                   "graphs": graphs,
                   "checkpoint": os.path.abspath(args.checkpoint)}, f, indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    bins = np.linspace(0, 1, 51)
    ax[0, 0].hist(s[y == 1], bins, histtype="step", label="true edges", density=True)
    ax[0, 0].hist(s[y == 0], bins, histtype="step", label="false edges", density=True)
    ax[0, 0].axvline(delta, color="k", ls="--", lw=1)
    ax[0, 0].set_yscale("log"); ax[0, 0].set_xlabel("edge weight"); ax[0, 0].legend()

    t_s = [r["threshold"] for r in scan]
    ax[0, 1].plot(t_s, [r["efficiency"] for r in scan], "o-", label="efficiency")
    ax[0, 1].plot(t_s, [r["fake_rate"] for r in scan], "s-", label="fake rate")
    ax[0, 1].plot(t_s, [r["duplicate_rate"] for r in scan], "^-", label="duplicate rate")
    ax[0, 1].axhline(ideal["efficiency"], color="gray", ls=":", label="ideal classifier eff.")
    ax[0, 1].set_xlabel("edge-weight threshold"); ax[0, 1].set_ylim(0, 1.02); ax[0, 1].legend(fontsize=8)

    ptbins = np.array([0.5, 1, 1.25, 1.5, 2, 2.5, 3, 4, 5, 7, 10, 20])
    etabins = np.linspace(-2.5, 2.5, 21)
    for col, (var, b, lab) in enumerate([(0, ptbins, r"particle $p_T$ [GeV]"), (1, etabins, r"particle $\eta$")]):
        a = ax[1, col]
        for name, r, sty in [("GNN", rows, "o"), ("ideal classifier", rows_ideal, "s")]:
            if not len(r):
                continue
            e, err = binned_eff(r[:, var], r[:, 3], b)
            a.errorbar(0.5 * (b[1:] + b[:-1]), e, yerr=err, fmt=sty, ms=4, label=name)
        a.set_xlabel(lab); a.set_ylabel("tracking efficiency"); a.set_ylim(0, 1.05); a.legend(fontsize=8)
        if var == 0:
            a.set_xscale("log")
    fig.suptitle(f">= {args.min_nodes} {args.count}, purity >= {args.min_purity:.2f}; delta = {delta:.3f}")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "performance.png"), dpi=120)
    print(f"[done] {out}/summary.json, performance.png")


if __name__ == "__main__":
    main()
