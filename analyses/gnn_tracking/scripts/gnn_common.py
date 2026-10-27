"""Shared helpers for the GNN edge-classification study.

Reading: the `input` argument follows analyses/common/InputFiles.h -- a single
.root file, a directory (every *.root inside, sorted), or a comma-separated mix.

Graph shards (written by build_graphs.py, read by train.py / evaluate.py) are
.npz files holding many events as concatenated arrays plus *_ptr offsets:

  node_feat  (N, len(NODE_FEATURES))   node_pid (N)  node_layer (N)   node_ptr (nev+1)
  node_helix_res (N): phi residual to the particle's production helix (NaN if no particle)
  event_file, event_entry (nev): source file (index into graph_config.json "files") and entry
  edge_index (2, E), local to event    edge_feat (E, len(EDGE_FEATURES))
  edge_y (E)                           edge_ptr (nev+1)
  part_pt/eta/phi/pdg/nstub (P)        part_ptr (nev+1)   n_true_pairs (nev)

node_pid indexes the event's particle table; -1 = stub with no TrackingParticle.
"""
import glob
import json
import os

import numpy as np

TREE = "L1TrackHitNtupleMaker/eventTree"
BEND_SENTINEL = 999999.0  # DummyStub production: no bend information

NODE_FEATURES = ["r", "phi", "z", "bend", "isPS", "isBarrel"]
# z0: beamline intercept of the pair [cm], clipped to +-100; curv: signed q/pT of the
# prompt helix through the beamline and the pair [1/GeV] (see pair_geometry)
EDGE_FEATURES = ["dr", "dphi", "dz", "deta", "dR", "z0", "curv"]
SIGNED_EDGE_FEATURES = ["dr", "dphi", "dz", "deta", "curv"]  # flip sign when the edge is reversed
BFIELD = 3.8  # T
CURV_PER_INV_GEV = 0.003 * BFIELD  # curvature [1/cm] of a 1 GeV track: 0.3 B / (100 pT)

STUB_BRANCHES = [
    "allstub_x", "allstub_y", "allstub_z", "allstub_layer", "allstub_isBarrel",
    "allstub_isPSmodule", "allstub_trigBend", "allstub_matchTP_pdgid",
    "allstub_matchTP_pt", "allstub_matchTP_eta", "allstub_matchTP_phi",
]


def list_root_files(inp, nfiles=0):
    out = []
    for item in (s.strip().rstrip("/") for s in inp.split(",")):
        if not item:
            continue
        item = os.path.expanduser(item)
        if os.path.isdir(item):
            files = sorted(glob.glob(os.path.join(item, "*.root")))
            out += files[: nfiles] if nfiles > 0 else files
        elif os.path.exists(item):
            out.append(item)
        else:
            print(f"[input] not found: {item}")
    print(f"[input] {len(out)} file(s) from {inp}")
    return out


def iter_events(files, max_events=0):
    """Yield (file_index, entry, dict of per-stub numpy arrays) for every event."""
    import uproot

    n = 0
    for ifile, f in enumerate(files):
        with uproot.open(f) as rf:
            arr = rf[TREE].arrays(STUB_BRANCHES, library="np")
        for iev in range(len(arr["allstub_x"])):
            yield ifile, iev, {b: arr[b][iev] for b in STUB_BRANCHES}
            n += 1
            if max_events and n >= max_events:
                return


def layer_id(is_barrel, layer, z):
    """Barrel layers 1-6 -> 1..6; endcap disks 1-5 -> 11..15 (+z) or 21..25 (-z)."""
    lid = np.where(is_barrel == 1, layer, np.where(z > 0, 10 + layer, 20 + layer))
    return lid.astype(np.int16)


def event_stubs(ev, drop_unmatched=True):
    """Convert one event's branches into node arrays and a particle table."""
    x, y, z = (ev[k].astype(np.float64) for k in ("allstub_x", "allstub_y", "allstub_z"))
    r = np.hypot(x, y)
    phi = np.arctan2(y, x)
    bend = ev["allstub_trigBend"].astype(np.float64)
    bend = np.where(np.abs(bend) >= BEND_SENTINEL, 0.0, bend)
    is_barrel = ev["allstub_isBarrel"].astype(np.int16)
    lid = layer_id(is_barrel, ev["allstub_layer"], z)

    pdg = ev["allstub_matchTP_pdgid"]
    matched = pdg != -999
    # No TP index is stored per stub; (pdgid, pt, eta, phi) identifies the particle.
    keys = np.stack([pdg.astype(np.float64), ev["allstub_matchTP_pt"],
                     ev["allstub_matchTP_eta"], ev["allstub_matchTP_phi"]], axis=1)
    pid = np.full(len(x), -1, dtype=np.int32)
    if matched.any():
        uniq, inv = np.unique(keys[matched], axis=0, return_inverse=True)
        pid[matched] = inv.reshape(-1)
    else:
        uniq = np.zeros((0, 4))

    keep = matched if drop_unmatched else np.ones(len(x), dtype=bool)
    nodes = {
        "feat": np.stack([r, phi, z, bend, ev["allstub_isPSmodule"], is_barrel], axis=1)[keep],
        "pid": pid[keep],
        "layer": lid[keep],
    }
    nstub = np.bincount(nodes["pid"][nodes["pid"] >= 0], minlength=len(uniq))
    parts = {"pdg": uniq[:, 0].astype(np.int32), "pt": uniq[:, 1], "eta": uniq[:, 2],
             "phi": uniq[:, 3], "nstub": nstub.astype(np.int32)}
    return nodes, parts


def helix_phi_residual(nodes, parts, bfield=BFIELD):
    """phi(stub) minus phi of the particle's production helix at the stub radius.

    Prompt helix: phi(r) = phi0 - q * asin(r / 2R), R [cm] = 100 pT / (0.3 B).
    Large residuals are bremsstrahlung, looping or displaced segments of the TP.
    NaN for stubs without a particle.
    """
    res = np.full(len(nodes["pid"]), np.nan)
    m = nodes["pid"] >= 0
    if m.any():
        p = nodes["pid"][m]
        pdg = parts["pdg"][p]
        q = np.where(np.isin(np.abs(pdg), [11, 13, 15]), -np.sign(pdg), np.sign(pdg))
        radius = 100.0 * parts["pt"][p] / (0.3 * bfield)
        r, phi = nodes["feat"][m, 0], nodes["feat"][m, 1]
        res[m] = wrap_phi(phi - parts["phi"][p] + q * np.arcsin(np.clip(r / (2 * radius), -1, 1)))
    return res


def arc_length(r, curv):
    """Transverse path length from the beamline to radius r on a circle of curvature curv [1/cm]."""
    k = np.abs(curv)
    x = np.clip(r * k / 2, 0, 1)
    return np.where(k * r < 1e-6, r, 2 * np.arcsin(x) / np.maximum(k, 1e-12))


def pair_geometry(ra, phia, za, rb, phib, zb):
    """Curvature and z0 of the prompt helix through the beamline and two stubs a -> b.

    curv: signed curvature [1/cm] of the circle through (0, 0), a and b in x-y
          (d0 = 0 assumed); |curv| = 0.3 B / (100 pT). The sign follows the charge
          for an edge oriented inner -> outer.
    z0:   intercept at the beamline of the straight line through a and b in the
          (s, z) plane, with s the path length along that circle (not r, which
          matters at low pT).
    """
    xa, ya, xb, yb = ra * np.cos(phia), ra * np.sin(phia), rb * np.cos(phib), rb * np.sin(phib)
    cross = xa * yb - ya * xb
    chord = np.maximum(np.hypot(xb - xa, yb - ya), 1e-6)
    curv = 2 * cross / (ra * rb * chord)
    sa, sb = arc_length(ra, curv), arc_length(rb, curv)
    ds = sb - sa
    ds = np.where(np.abs(ds) < 1e-3, np.where(ds < 0, -1e-3, 1e-3), ds)
    z0 = za - sa * (zb - za) / ds
    return curv, z0


def vertex_cone_residual(ra, za, rb, zb, curv, z0_max):
    """How far stub b lies outside the cone of lines from the luminous region through stub a.

    The lines from (s=0, z0), |z0| <= z0_max, through a reach z_b in
    z_a t -+ z0_max (t - 1), with t = s_b / s_a. Returns the distance of z_b outside
    that interval [cm]; 0 means compatible with a prompt track from |z0| <= z0_max.
    Unlike extrapolating z0, this has no small denominator for disk-to-disk pairs.
    """
    t = arc_length(rb, curv) / arc_length(ra, curv)
    return np.maximum(0.0, np.abs(zb - za * t) - z0_max * np.abs(t - 1))


def particle_z0(nodes, parts):
    """Per-particle production z from its own stubs: median of z - s sinh(eta_TP).

    Uses the particle's pT (for s) and eta from the TrackingParticle, and only PS
    stubs when there are any (z in barrel 2S modules is known only to 5 cm).
    NaN for particles without stubs.
    """
    f, pid = nodes["feat"], nodes["pid"]
    z0 = np.full(len(parts["pt"]), np.nan)
    m = pid >= 0
    if not m.any():
        return z0
    p = pid[m]
    s = arc_length(f[m, 0], CURV_PER_INV_GEV / parts["pt"][p])
    est = f[m, 2] - s * np.sinh(parts["eta"][p])
    ps = f[m, 4] > 0
    for q in np.unique(p):
        sel = p == q
        use = sel & ps if (sel & ps).any() else sel
        z0[q] = np.median(est[use])
    return z0


def wrap_phi(dphi):
    return (dphi + np.pi) % (2 * np.pi) - np.pi


def eta_of(r, z):
    return np.arcsinh(z / r)


def load_shards(graph_dir):
    files = sorted(glob.glob(os.path.join(graph_dir, "graphs_*.npz")))
    if not files:
        raise SystemExit(f"no graphs_*.npz in {graph_dir}")
    return [dict(np.load(f)) for f in files]


def split_events(shards):
    """Flatten shards into a list of per-event dicts of numpy views."""
    events = []
    for s in shards:
        for i in range(len(s["node_ptr"]) - 1):
            n0, n1 = s["node_ptr"][i], s["node_ptr"][i + 1]
            e0, e1 = s["edge_ptr"][i], s["edge_ptr"][i + 1]
            p0, p1 = s["part_ptr"][i], s["part_ptr"][i + 1]
            events.append({
                "x": s["node_feat"][n0:n1], "pid": s["node_pid"][n0:n1],
                "layer": s["node_layer"][n0:n1], "helix_res": s["node_helix_res"][n0:n1],
                "edge_index": s["edge_index"][:, e0:e1], "edge_attr": s["edge_feat"][e0:e1],
                "y": s["edge_y"][e0:e1],
                "part_pt": s["part_pt"][p0:p1], "part_eta": s["part_eta"][p0:p1],
                "part_phi": s["part_phi"][p0:p1],
                "part_pdg": s["part_pdg"][p0:p1], "part_nstub": s["part_nstub"][p0:p1],
                "n_true_pairs": int(s["n_true_pairs"][i]),
                "file": int(s["event_file"][i]) if "event_file" in s else -1,
                "entry": int(s["event_entry"][i]) if "event_entry" in s else -1,
            })
    return events


def read_json(path):
    with open(path) as f:
        return json.load(f)
