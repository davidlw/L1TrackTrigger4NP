"""Stage-by-stage loss of true muons in the hybrid tracklet chain, from truthtrace.txt.

Records (one muon = (job, event, tp)):
  E ev nsim nstub | T tp type pt eta phi vz | S tp layerdisk region r phi z
  C tp iseed sector phimin rinv phi0 z0 t          seed pair reached the calculator
  F tp iseed sector rz|sec ...                     calculator rejected it
  P tp iseed sector layer status phiEx phiApprox phiInt zEx zApprox zInt rinvInt
  A tp iseed sector layer added                    projection routed to matching
  M tp iseed phimin layer stubtp dphiInt cut dzInt zcut imatch dphiEx dphiApprox dz
  B tp iseed sector nMatches nMatchesUniq layermask
"""
import glob, sys, os, numpy as np
from collections import defaultdict

BARREL_SEEDS = {0: (1, 2), 1: (2, 3), 2: (3, 4), 3: (5, 6)}
RMEAN = {1: 25.18, 2: 37.36, 3: 52.70, 4: 68.71, 5: 86.24, 6: 108.2}  # nominal, only for printing

def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi

def read(path, job):
    mu = {}            # key -> dict
    ev = None
    phimin = {}        # (key, iseed, sector) -> phimin
    for line in open(path):
        f = line.split()
        if not f:
            continue
        t = f[0]
        if t == "E":
            ev = int(f[1]); continue
        if ev is None:
            continue
        tp = int(f[1]); key = (job, ev, tp)
        m = mu.setdefault(key, {"stubs": defaultdict(list), "C": set(), "F": defaultdict(list), "P": [], "A": set(),
                                "M": [], "B": [], "phimin": {}})
        if t == "T":
            m["pt"], m["eta"], m["phi"], m["vz"] = float(f[3]), float(f[4]), float(f[5]), float(f[6]); m["type"] = int(f[2])
        elif t == "S":
            ld = int(f[2])
            m["stubs"][ld + 1].append((float(f[4]), float(f[5]), float(f[6])))   # layer (1-6 barrel, 7-11 disks)
        elif t == "C":
            iseed, sec = int(f[2]), int(f[3]); m["C"].add((iseed, sec)); m["phimin"][(iseed, sec)] = float(f[4])
            m.setdefault("rinv", {})[(iseed, sec)] = float(f[5])
        elif t == "F":
            m["F"][(int(f[2]), int(f[3]))].append(f[4])
        elif t == "P":
            m["P"].append((int(f[2]), int(f[3]), int(f[4]), int(f[5]), float(f[6]), float(f[7]), float(f[8]),
                           float(f[9]), float(f[10]), float(f[11]), float(f[12])))
        elif t == "A":
            if int(f[5]):
                m["A"].add((int(f[2]), int(f[3]), int(f[4])))
        elif t == "M":
            m["M"].append((int(f[2]), float(f[3]), int(f[4]), int(f[5]), float(f[6]), float(f[7]), float(f[8]),
                           float(f[9]), int(f[10]), float(f[11]), float(f[12])))
        elif t == "B":
            m["B"].append((int(f[2]), int(f[3]), int(f[4]), int(f[5]), int(f[6])))
    return {k: v for k, v in mu.items() if "pt" in v}

def ntuple_found(root_path):
    import uproot
    a = uproot.open(root_path)["L1TrackHitNtupleMaker/eventTree"].arrays(["tp_pt", "tp_eta", "tp_phi", "tp_nmatch"], library="np")
    out = []
    for i in range(len(a["tp_pt"])):
        out.append(list(zip(a["tp_pt"][i], a["tp_eta"][i], a["tp_phi"][i], a["tp_nmatch"][i])))
    return out

def main(rundir, out):
    muons = {}
    found_by_event = {}
    for jd in sorted(glob.glob(os.path.join(rundir, "job*"))):
        job = os.path.basename(jd)
        tr = os.path.join(jd, "truthtrace.txt")
        if not os.path.exists(tr):
            continue
        muons.update(read(tr, job))
        nt = os.path.join(jd, "diag_ntuple.root")
        if os.path.exists(nt):
            for iev, tps in enumerate(ntuple_found(nt)):
                found_by_event[(job, iev + 1)] = tps
    # link to the ntuple's hybrid efficiency (tp_nmatch > 0) by matching pT/eta/phi
    for (job, ev, tp), m in muons.items():
        m["found"] = None
        for pt, eta, phi, nm in found_by_event.get((job, ev), []):
            if abs(pt - m["pt"]) < 1e-3 * m["pt"] and abs(eta - m["eta"]) < 1e-3 and abs(wrap(phi - m["phi"])) < 1e-3:
                m["found"] = nm > 0
    rows = []
    for key, m in muons.items():
        if abs(m["type"]) != 13:
            continue
        lay = set(m["stubs"])
        seedable = [s for s, (a, b) in BARREL_SEEDS.items() if a in lay and b in lay]
        tried = {s for s, _ in m["C"]}
        passed = {(s, sec) for (s, sec) in m["C"] if (s, sec) not in m["F"]}
        b_uniq = max([b[3] for b in m["B"]], default=-1)
        rows.append(dict(key=key, pt=m["pt"], eta=m["eta"], nlay_barrel=len([l for l in lay if l <= 6]),
                         seedable=bool(seedable), tried=bool(tried & set(seedable)) or bool(tried),
                         passed=bool(passed), reached_builder=b_uniq >= 0, builder_ok=b_uniq >= 2,
                         found=m["found"], m=m))
    return rows, muons

if __name__ == "__main__":
    rows, muons = main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else ".")
    import pickle
    pickle.dump(rows, open(os.path.join(sys.argv[2] if len(sys.argv) > 2 else ".", "trace_rows.pkl"), "wb"))
    print(len(rows), "muons")
