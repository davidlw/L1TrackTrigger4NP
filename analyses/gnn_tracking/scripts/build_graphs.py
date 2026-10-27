#!/usr/bin/env python3
"""Build stub graphs from L1TrackHitNtuple files (step 1 of the GNN pipeline).

Nodes are stubs; candidate edges connect stubs on an allowed pair of layers
within a Delta-phi and Delta-eta window. Both the layer pairs and the windows
are derived from truth rather than tuned by hand ("data-driven" graphs in
DeZoort et al., arXiv:2103.16701, Sec. 4.2):

  * for every particle, order its stubs by distance from the origin and take
    stub pairs on consecutive layers and on layers one apart (one missing stub),
    using only stubs within --helix-tol of the particle's production helix in phi
    (electrons: ~20% of stubs are post-bremsstrahlung or looping segments);
  * a layer pair is allowed if it occurs at least --min-pair-count times;
  * per layer pair, the |dphi| and |deta| windows are the --quantile quantile of
    the closest such pair of each particle, times --safety.

Two physics preconditions then remove pairs that no prompt track above pT_min
from the luminous region can produce:

  * pT:  the circle through the beamline (d0 = 0) and the two stubs has
         curvature |kappa| <= 0.3 B / pT_min (+ tolerance)          --pt-min
  * z0:  stub b lies inside the cone of lines from |z0| <= Z through stub a,
         in the r-z plane with the arc length of that circle      --z0-max
         (+ tolerance, in cm of z at stub b)

Their tolerances are per layer pair, the --quantile quantile of the measured
excess over truth for prompt particles (|z0| < Z), e.g. 2S strips give z only
to 5 cm in the barrel and r only to 5 cm in the disks. Set either to 0 to
disable it.

An edge is labelled true (y=1) when both stubs come from the same particle, so
connected components of the true edges are exactly the particles.

The windows are written to <out>/graph_config.json. Pass --config to reuse the
windows of another sample (e.g. build the test sample with the training cuts).

Usage:
  python3 build_graphs.py --input ~/ntuples/DummyStub/STARlight_QED_ee --out ../output/graphs_qedee
  python3 build_graphs.py --input ~/ntuples/DummyStub/EPOS_pPb --out ../output/graphs_ppb \
      --config ../output/graphs_qedee/graph_config.json
"""
import argparse
import json
import os
from collections import defaultdict

import numpy as np

from gnn_common import (CURV_PER_INV_GEV, EDGE_FEATURES, NODE_FEATURES, eta_of, event_stubs,
                        helix_phi_residual, iter_events, list_root_files, pair_geometry,
                        particle_z0, vertex_cone_residual, wrap_phi)

R, PHI, Z = 0, 1, 2


def true_neighbour_pairs(nodes, max_gap=2):
    """Same-particle stub pairs whose layers are <= max_gap apart along the track.

    Returns an array (M, 3): inner stub, outer stub, gap (1 = consecutive layers).
    """
    feat, pid, lay = nodes["feat"], nodes["pid"], nodes["layer"]
    rho = np.hypot(feat[:, R], feat[:, Z])
    out = []
    for p in np.unique(pid[pid >= 0]):
        idx = np.where(pid == p)[0]
        idx = idx[np.argsort(rho[idx])]
        # layers in the order the particle crosses them (first appearance)
        seq = list(dict.fromkeys(lay[idx].tolist()))
        for k, la in enumerate(seq):
            for gap in range(1, max_gap + 1):
                if k + gap >= len(seq):
                    break
                lb = seq[k + gap]
                for a in idx[lay[idx] == la]:
                    for b in idx[lay[idx] == lb]:
                        out.append((a, b, gap))
    return np.array(out, dtype=np.int64).reshape(-1, 3)


def derive_config(files, args):
    dphi, deta, dcurv, cone = defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list)
    ntrue = 0
    for _, _, ev in iter_events(files, args.derive_events):
        nodes, parts = event_stubs(ev, drop_unmatched=True)
        pairs = true_neighbour_pairs(nodes)
        if len(pairs) == 0:
            continue
        pt = parts["pt"][nodes["pid"][pairs[:, 0]]]
        on_helix = np.abs(helix_phi_residual(nodes, parts)) < args.helix_tol
        pairs = pairs[(pt >= args.derive_pt_min) & on_helix[pairs[:, 0]] & on_helix[pairs[:, 1]]]
        if len(pairs) == 0:
            continue
        z0_true = particle_z0(nodes, parts)
        f, lay, pid = nodes["feat"], nodes["layer"], nodes["pid"]
        a, b = pairs[:, 0], pairs[:, 1]
        dp = np.abs(wrap_phi(f[b, PHI] - f[a, PHI]))
        de = np.abs(eta_of(f[b, R], f[b, Z]) - eta_of(f[a, R], f[a, Z]))
        # A particle stays connected if ONE of its stub pairs on a layer pair survives,
        # so take each particle's closest pair. Using every pair would let looping and
        # bremsstrahlung electrons (tens of stubs on one TP) set the windows.
        best = {}
        for i, (la, lb) in enumerate(zip(lay[a], lay[b])):
            key = (pid[a[i]], min(la, lb), max(la, lb))
            d = np.hypot(dp[i], de[i])
            if key not in best or d < best[key][0]:
                best[key] = (d, dp[i], de[i], i)
        # precondition tolerances, from the same closest pairs of prompt particles
        for (p, la, lb), (_, x, y, i) in best.items():
            dphi[(la, lb)].append(x)
            deta[(la, lb)].append(y)
            if abs(z0_true[p]) < args.z0_ref:
                ia, ibb = a[i], b[i]
                k, _ = pair_geometry(f[ia, R], f[ia, PHI], f[ia, Z], f[ibb, R], f[ibb, PHI], f[ibb, Z])
                dcurv[(la, lb)].append(abs(k) - CURV_PER_INV_GEV / parts["pt"][p])
                cone[(la, lb)].append(vertex_cone_residual(f[ia, R], f[ia, Z], f[ibb, R], f[ibb, Z],
                                                           k, args.z0_ref))
        ntrue += len(best)

    config = {"quantile": args.quantile, "safety": args.safety, "derive_pt_min": args.derive_pt_min,
              "helix_tol": args.helix_tol, "z0_ref": args.z0_ref,
              "min_pair_count": args.min_pair_count, "n_true_pairs_used": ntrue,
              "node_features": NODE_FEATURES, "edge_features": EDGE_FEATURES,
              "drop_unmatched": not args.keep_unmatched, "layer_pairs": []}
    for key in sorted(dphi):
        n = len(dphi[key])
        if n < args.min_pair_count:
            continue
        q = lambda v: float(np.quantile(v, args.quantile) * args.safety) if len(v) else float("inf")
        config["layer_pairs"].append({
            "layers": [int(key[0]), int(key[1])], "count": n,
            "dphi_max": float(min(np.quantile(dphi[key], args.quantile) * args.safety, 0.99 * np.pi)),
            "deta_max": q(deta[key]),
            # excess |curvature| [1/cm] and vertex-cone residual [cm] tolerated by the preconditions
            "curv_tol": q(np.maximum(dcurv[key], 0)), "cone_tol": q(cone[key]),
        })
    kept = sum(p["count"] for p in config["layer_pairs"])
    print(f"[config] {len(config['layer_pairs'])} layer pairs kept, covering "
          f"{kept}/{ntrue} (particle, layer pair) links ({kept / max(ntrue, 1):.4f})")
    return config


def pairs_in_window(ia, ib, phi, cut):
    """All (a, b) with a in ia, b in ib and |dphi| < cut, via a phi-sorted window search."""
    order = ib[np.argsort(phi[ib])]
    pb = phi[order]
    pext = np.concatenate([pb - 2 * np.pi, pb, pb + 2 * np.pi])
    iext = np.concatenate([order, order, order])
    lo = np.searchsorted(pext, phi[ia] - cut, "left")
    hi = np.searchsorted(pext, phi[ia] + cut, "right")
    cnt = hi - lo
    tot = int(cnt.sum())
    if tot == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    src = np.repeat(ia, cnt)
    start = np.repeat(lo - np.concatenate([[0], np.cumsum(cnt)[:-1]]), cnt)
    return src, iext[np.arange(tot) + start]


def build_event_graph(nodes, pair_cuts, on_helix, pt_min=0.0, z0_max=0.0):
    f, lay, pid = nodes["feat"], nodes["layer"], nodes["pid"]
    r, phi, z = f[:, R], f[:, PHI], f[:, Z]
    eta = eta_of(r, z)
    rho = np.hypot(r, z)
    by_layer = {l: np.where(lay == l)[0] for l in np.unique(lay)}

    def same_particle_pairs(ia, ib):
        pa, pb = pid[ia], pid[ib]
        pa, pb = pa[pa >= 0], pb[pb >= 0]
        if len(pa) == 0 or len(pb) == 0:
            return 0
        m = max(pa.max(), pb.max()) + 1
        return int(np.dot(np.bincount(pa, minlength=m), np.bincount(pb, minlength=m)))

    srcs, dsts = [], []
    n_true_pairs = n_true_pairs_helix = n_window = n_window_true = 0
    curv_max = CURV_PER_INV_GEV / pt_min if pt_min > 0 else np.inf
    for (la, lb), (cphi, ceta, ctol_curv, ctol_cone) in pair_cuts.items():
        ia, ib = by_layer.get(la), by_layer.get(lb)
        if ia is None or ib is None:
            continue
        # truth pairs available on this layer pair, before any window cut
        n_true_pairs += same_particle_pairs(ia, ib)
        n_true_pairs_helix += same_particle_pairs(ia[on_helix[ia]], ib[on_helix[ib]])

        a, b = pairs_in_window(ia, ib, phi, cphi)
        ok = np.abs(eta[b] - eta[a]) < ceta
        a, b = a[ok], b[ok]
        n_window += len(a)
        n_window_true += int(((pid[a] == pid[b]) & (pid[a] >= 0)).sum())
        if pt_min > 0 or z0_max > 0:
            inner = rho[a] <= rho[b]
            i, o = np.where(inner, a, b), np.where(inner, b, a)
            curv, _ = pair_geometry(r[i], phi[i], z[i], r[o], phi[o], z[o])
            ok = np.abs(curv) <= curv_max + ctol_curv
            if z0_max > 0:
                ok &= vertex_cone_residual(r[i], z[i], r[o], z[o], curv, z0_max) <= ctol_cone
            a, b = a[ok], b[ok]
        srcs.append(a)
        dsts.append(b)

    if srcs:
        a, b = np.concatenate(srcs), np.concatenate(dsts)
    else:
        a = b = np.zeros(0, np.int64)
    # orient every edge inner -> outer
    flip = rho[a] > rho[b]
    a, b = np.where(flip, b, a), np.where(flip, a, b)

    dphi = wrap_phi(phi[b] - phi[a])
    deta = eta[b] - eta[a]
    curv, z0 = pair_geometry(r[a], phi[a], z[a], r[b], phi[b], z[b])
    edge_feat = np.stack([r[b] - r[a], dphi, z[b] - z[a], deta, np.hypot(deta, dphi),
                          np.clip(z0, -100, 100), curv / CURV_PER_INV_GEV], axis=1)
    y = ((pid[a] == pid[b]) & (pid[a] >= 0)).astype(np.uint8)
    counts = {"true_pairs": n_true_pairs, "true_pairs_helix": n_true_pairs_helix,
              "window_edges": n_window, "window_true_edges": n_window_true,
              "true_edges_helix": int((y.astype(bool) & on_helix[a] & on_helix[b]).sum())}
    return np.stack([a, b]).astype(np.int32), edge_feat.astype(np.float32), y, counts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, help="ntuple file, directory, or comma-separated list")
    ap.add_argument("--out", required=True, help="output directory for graph shards")
    ap.add_argument("--nfiles", type=int, default=0, help="use the first N files of a directory (0 = all)")
    ap.add_argument("--max-events", type=int, default=0, help="stop after N events in total (0 = all)")
    ap.add_argument("--config", help="reuse layer pairs and windows from this graph_config.json")
    ap.add_argument("--derive-events", type=int, default=50000, help="events used to derive the windows")
    ap.add_argument("--quantile", type=float, default=0.995, help="window = this quantile of true pairs")
    ap.add_argument("--safety", type=float, default=1.1, help="window scale factor on top of the quantile")
    ap.add_argument("--pt-min", type=float, default=None,
                    help="precondition: pair compatible with a prompt track above this pT [GeV] "
                         "(default 0.5; 0 = off; with --config, the config's value)")
    ap.add_argument("--z0-max", type=float, default=None,
                    help="precondition: pair compatible with a track from |z0| < this [cm] "
                         "(default 15; 0 = off; with --config, the config's value)")
    ap.add_argument("--z0-ref", type=float, default=15.0,
                    help="prompt particles (|z0| < this [cm]) used to measure the precondition tolerances")
    ap.add_argument("--derive-pt-min", type=float, default=0.0,
                    help="derive windows from particles above this pT [GeV]")
    ap.add_argument("--helix-tol", type=float, default=0.05,
                    help="derive windows from stubs within this |dphi| [rad] of the production helix")
    ap.add_argument("--min-pair-count", type=int, default=20, help="minimum particles linking a layer pair to allow it")
    ap.add_argument("--keep-unmatched", action="store_true",
                    help="keep stubs with no TrackingParticle (labelled as noise). By default they are "
                         "dropped: with ptMinTP = 1 GeV they include real low-pT particles")
    ap.add_argument("--events-per-shard", type=int, default=20000)
    args = ap.parse_args()

    files = list_root_files(args.input, args.nfiles)
    if not files:
        raise SystemExit("no input files")
    os.makedirs(args.out, exist_ok=True)

    if args.config:
        with open(args.config) as fcfg:
            config = json.load(fcfg)
        print(f"[config] reusing {args.config}")
        if config.get("drop_unmatched", True) == args.keep_unmatched:
            print("[config] WARNING: --keep-unmatched differs from the config's sample")
    else:
        config = derive_config(files, args)
    for key, default in (("pt_min", 0.5), ("z0_max", 15.0)):
        val = getattr(args, key)
        config[key] = val if val is not None else config.get(key, default)
    print(f"[config] preconditions: pT > {config['pt_min']} GeV, |z0| < {config['z0_max']} cm (0 = off)")
    config["drop_unmatched"] = not args.keep_unmatched
    config["input"] = args.input
    config["files"] = [os.path.abspath(f) for f in files]
    with open(os.path.join(args.out, "graph_config.json"), "w") as fcfg:
        json.dump(config, fcfg, indent=1)
    pair_cuts = {tuple(p["layers"]): (p["dphi_max"], p["deta_max"], p["curv_tol"], p["cone_tol"])
                 for p in config["layer_pairs"]}

    buf = defaultdict(list)
    nshard = 0
    tot = defaultdict(int)

    def flush():
        nonlocal nshard
        if not buf["node_feat"]:
            return
        cat = lambda k, axis=0: np.concatenate(buf[k], axis=axis)
        ptr = lambda k: np.concatenate([[0], np.cumsum(buf[k])]).astype(np.int64)
        np.savez_compressed(
            os.path.join(args.out, f"graphs_{nshard:03d}.npz"),
            node_feat=cat("node_feat"), node_pid=cat("node_pid"), node_layer=cat("node_layer"),
            node_helix_res=cat("node_helix_res"),
            node_ptr=ptr("n_nodes"), edge_index=cat("edge_index", 1), edge_feat=cat("edge_feat"),
            edge_y=cat("edge_y"), edge_ptr=ptr("n_edges"), part_ptr=ptr("n_parts"),
            part_pt=cat("part_pt"), part_eta=cat("part_eta"), part_phi=cat("part_phi"),
            part_pdg=cat("part_pdg"), part_nstub=cat("part_nstub"),
            n_true_pairs=np.array(buf["n_true_pairs"], dtype=np.int64),
            event_file=np.array(buf["event_file"], dtype=np.int32),
            event_entry=np.array(buf["event_entry"], dtype=np.int64))
        nshard += 1
        buf.clear()

    for ifile, ientry, ev in iter_events(files, args.max_events):
        nodes, parts = event_stubs(ev, drop_unmatched=not args.keep_unmatched)
        helix_res = helix_phi_residual(nodes, parts)
        ei, ef, y, counts = build_event_graph(nodes, pair_cuts, np.abs(helix_res) < config["helix_tol"],
                                              config["pt_min"], config["z0_max"])
        n_true_pairs = counts["true_pairs"]
        buf["node_feat"].append(nodes["feat"].astype(np.float32))
        buf["node_pid"].append(nodes["pid"])
        buf["node_layer"].append(nodes["layer"])
        buf["node_helix_res"].append(helix_res.astype(np.float32))
        buf["edge_index"].append(ei)
        buf["edge_feat"].append(ef)
        buf["edge_y"].append(y)
        for k in ("pt", "eta", "phi", "pdg", "nstub"):
            buf["part_" + k].append(parts[k])
        buf["n_nodes"].append(len(nodes["pid"]))
        buf["n_edges"].append(len(y))
        buf["n_parts"].append(len(parts["pt"]))
        buf["n_true_pairs"].append(n_true_pairs)
        buf["event_file"].append(ifile)
        buf["event_entry"].append(ientry)
        tot["events"] += 1
        tot["nodes"] += len(nodes["pid"])
        tot["edges"] += len(y)
        tot["true_edges"] += int(y.sum())
        for k, v in counts.items():
            tot[k] += v
        if len(buf["n_nodes"]) >= args.events_per_shard:
            flush()
    flush()

    ne = max(tot["events"], 1)
    print(f"[graphs] {tot['events']} events in {nshard} shard(s): "
          f"{tot['nodes'] / ne:.1f} nodes, {tot['edges'] / ne:.1f} edges per event")
    print(f"[graphs] graph efficiency = true edges kept / same-particle pairs on allowed layer pairs:\n"
          f"         all stubs {tot['true_edges'] / max(tot['true_pairs'], 1):.4f}, "
          f"on-helix stubs {tot['true_edges_helix'] / max(tot['true_pairs_helix'], 1):.4f}")
    print(f"[graphs] graph purity (true / all edges) = {tot['true_edges'] / max(tot['edges'], 1):.4f}")
    print(f"[graphs] preconditions keep {tot['edges'] / max(tot['window_edges'], 1):.4f} of the "
          f"window edges and {tot['true_edges'] / max(tot['window_true_edges'], 1):.4f} of the true ones")


if __name__ == "__main__":
    main()
