#!/usr/bin/env python3
"""Train the interaction-network edge classifier on graphs from build_graphs.py.

Events are split at random (fixed --seed) into train / validation / test; the
split is stored in the checkpoint so evaluate.py uses the same test events.
Loss is binary cross-entropy on both directions of every edge. The working
point delta is the edge-weight threshold where TPR = TNR on the validation set
(DeZoort et al., Sec. 4.3). The checkpoint with the lowest validation loss is kept.

Usage:
  python3 train.py --graphs ../output/graphs_qedee --out ../output/model_qedee
  python3 train.py --graphs ../output/graphs_qedee --out ../output/model_qedee_bend \
      --node-features r,phi,z,bend        # DefaultStub samples only
"""
import argparse
import json
import os
import time

import numpy as np
import torch

from gnn_common import EDGE_FEATURES, NODE_FEATURES, SIGNED_EDGE_FEATURES, load_shards, split_events
from model import InteractionNetwork, make_bidirectional, edge_scores


def roc_curve(scores, y):
    """Thresholds (descending) with the TPR and TNR of `score > threshold`."""
    order = np.argsort(-scores, kind="stable")
    s, t = scores[order], y[order].astype(bool)
    tp = np.cumsum(t)
    fp = np.cumsum(~t)
    last = np.r_[np.diff(s) != 0, True]  # one point per distinct score
    npos, nneg = max(t.sum(), 1), max((~t).sum(), 1)
    return s[last], tp[last] / npos, 1 - fp[last] / nneg


def auc(tpr, tnr):
    fpr = 1 - tnr
    return float(np.trapezoid(np.r_[0, tpr], np.r_[0, fpr]))


def balanced_threshold(scores, y):
    thr, tpr, tnr = roc_curve(scores, y)
    i = np.argmin(np.abs(tpr - tnr))
    return float(thr[i]), float(tpr[i]), float(tnr[i]), auc(tpr, tnr)


class GraphSet:
    """Events with selected node columns, normalised; yields concatenated batches."""

    def __init__(self, events, node_cols, edge_cols, norm):
        self.events = events
        self.node_cols = node_cols
        self.edge_cols = edge_cols
        self.norm = norm

    def batch(self, idx, device):
        xs, eis, eas, ys, off = [], [], [], [], 0
        for i in idx:
            ev = self.events[i]
            xs.append(ev["x"][:, self.node_cols])
            eis.append(ev["edge_index"].astype(np.int64) + off)
            eas.append(ev["edge_attr"][:, self.edge_cols])
            ys.append(ev["y"])
            off += len(ev["x"])
        x = (np.concatenate(xs) - self.norm["node_mean"]) / self.norm["node_std"]
        ea = np.concatenate(eas) / self.norm["edge_rms"]
        t = lambda a, dt=torch.float32: torch.as_tensor(np.ascontiguousarray(a), dtype=dt, device=device)
        return t(x), t(np.concatenate(eis, axis=1), torch.long), t(ea), t(np.concatenate(ys))

    def batches(self, idx, batch_size, device, shuffle=False, rng=None):
        idx = np.array(idx)
        if shuffle:
            idx = idx[rng.permutation(len(idx))]
        for k in range(0, len(idx), batch_size):
            yield self.batch(idx[k:k + batch_size], device)


def compute_norm(events, idx, node_cols, edge_cols):
    x = np.concatenate([events[i]["x"][:, node_cols] for i in idx]).astype(np.float64)
    ea = np.concatenate([events[i]["edge_attr"][:, edge_cols] for i in idx]).astype(np.float64)
    std = x.std(axis=0)
    rms = np.sqrt((ea ** 2).mean(axis=0))  # no centring: reversing an edge flips the sign
    return {"node_mean": x.mean(axis=0).astype(np.float32),
            "node_std": np.where(std > 0, std, 1).astype(np.float32),
            "edge_rms": np.where(rms > 0, rms, 1).astype(np.float32)}


@torch.no_grad()
def predict(model, data, idx, batch_size, device, signed_cols):
    """Edge scores and labels for the events idx, in event order."""
    model.eval()
    scores, ys = [], []
    for x, ei, ea, y in data.batches(idx, batch_size, device):
        if ei.shape[1] == 0:
            continue
        scores.append(edge_scores(model, x, ei, ea, signed_cols).cpu().numpy())
        ys.append(y.cpu().numpy())
    return np.concatenate(scores), np.concatenate(ys)


def run_epoch(model, data, idx, batch_size, device, signed_cols, pos_weight, opt=None, rng=None):
    train = opt is not None
    model.train(train)
    tot, n = 0.0, 0
    lossf = torch.nn.BCEWithLogitsLoss(pos_weight=torch.tensor(pos_weight, device=device), reduction="sum")
    with torch.set_grad_enabled(train):
        for x, ei, ea, y in data.batches(idx, batch_size, device, shuffle=train, rng=rng):
            if ei.shape[1] == 0:
                continue
            ei2, ea2 = make_bidirectional(ei, ea, signed_cols)
            loss = lossf(model(x, ei2, ea2), torch.cat([y, y]))
            if train:
                opt.zero_grad()
                (loss / ei2.shape[1]).backward()
                opt.step()
            tot += loss.item()
            n += ei2.shape[1]
    return tot / max(n, 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--graphs", required=True, help="directory written by build_graphs.py")
    ap.add_argument("--out", required=True, help="output directory for checkpoint and plots")
    ap.add_argument("--node-features", default="r,phi,z", help=f"comma list from {NODE_FEATURES}")
    ap.add_argument("--edge-features", default="dr,dphi,dz,deta,dR",
                    help=f"comma list from {EDGE_FEATURES} (z0, curv: pair's beamline intercept and q/pT)")
    ap.add_argument("--hidden", type=int, default=40)
    ap.add_argument("--iterations", type=int, default=1, help="message-passing steps (paper: 1)")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--lr", type=float, default=5e-3)
    ap.add_argument("--step", type=int, default=10, help="StepLR: decay the learning rate every N epochs")
    ap.add_argument("--gamma", type=float, default=0.95, help="StepLR decay factor (paper: 0.95)")
    ap.add_argument("--batch-size", type=int, default=64, help="graphs (events) per batch")
    ap.add_argument("--pos-weight", type=float, default=1.0, help="BCE weight of true edges")
    ap.add_argument("--val-frac", type=float, default=0.1)
    ap.add_argument("--test-frac", type=float, default=0.1)
    ap.add_argument("--max-events", type=int, default=0, help="use at most N events (0 = all)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--exclude-events", help=".npy of event indices (in --graphs order) kept out of training and "
                                             "validation, e.g. events that are evaluated elsewhere")
    ap.add_argument("--device", default="cpu", help="cpu, cuda or mps")
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.out, exist_ok=True)

    events = split_events(load_shards(args.graphs))
    if args.max_events:
        events = events[: args.max_events]
    perm = rng.permutation(len(events))
    n_test, n_val = int(args.test_frac * len(events)), int(args.val_frac * len(events))
    test_idx, val_idx, train_idx = perm[:n_test], perm[n_test:n_test + n_val], perm[n_test + n_val:]
    excluded = set(np.load(args.exclude_events).tolist()) if args.exclude_events else set()
    if excluded:
        print(f"[data] {len(excluded)} events excluded from training and validation")
    has_edges = lambda idx: [i for i in idx if events[i]["edge_index"].shape[1] > 0 and i not in excluded]
    train_idx, val_idx = has_edges(train_idx), has_edges(val_idx)
    print(f"[data] {len(events)} events; with edges: {len(train_idx)} train, {len(val_idx)} validation; "
          f"{n_test} test")

    names = args.node_features.split(",")
    node_cols = [NODE_FEATURES.index(n) for n in names]
    if "bend" in names and not any(np.any(events[i]["x"][:, NODE_FEATURES.index("bend")]) for i in train_idx):
        print("[data] WARNING: bend is zero everywhere (DummyStub?) -- it carries no information")
    edge_names = args.edge_features.split(",")
    edge_cols = [EDGE_FEATURES.index(n) for n in edge_names]
    signed_cols = [k for k, n in enumerate(edge_names) if n in SIGNED_EDGE_FEATURES]
    norm = compute_norm(events, train_idx, node_cols, edge_cols)
    data = GraphSet(events, node_cols, edge_cols, norm)
    ytr = np.concatenate([events[i]["y"] for i in train_idx])
    print(f"[data] train edges: {len(ytr)}, true fraction {ytr.mean():.4f}")

    device = torch.device(args.device)
    model = InteractionNetwork(len(node_cols), len(edge_cols), hidden=args.hidden,
                               iterations=args.iterations).to(device)
    print(f"[model] {sum(p.numel() for p in model.parameters())} parameters")
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=args.step, gamma=args.gamma)

    history, best = [], None
    for epoch in range(1, args.epochs + 1):
        t0 = time.time()
        tr = run_epoch(model, data, train_idx, args.batch_size, device, signed_cols, args.pos_weight, opt, rng)
        va = run_epoch(model, data, val_idx, args.batch_size, device, signed_cols, args.pos_weight)
        s, y = predict(model, data, val_idx, args.batch_size, device, signed_cols)
        delta, tpr, tnr, a = balanced_threshold(s, y)
        acc = float(((s > delta) == y.astype(bool)).mean())
        sched.step()
        history.append({"epoch": epoch, "train_loss": tr, "val_loss": va, "delta": delta,
                        "val_tpr": tpr, "val_tnr": tnr, "val_acc": acc, "val_auc": a})
        flag = ""
        if best is None or va < best["val_loss"]:
            best = history[-1]
            torch.save({"model": model.state_dict(), "args": vars(args), "node_features": names,
                        "node_cols": node_cols, "edge_features": edge_names, "edge_cols": edge_cols,
                        "signed_cols": signed_cols,
                        "norm": {k: v.tolist() for k, v in norm.items()}, "delta": delta,
                        "epoch": epoch, "test_idx": test_idx.tolist(), "val_idx": list(map(int, val_idx)),
                        "graphs": os.path.abspath(args.graphs)},
                       os.path.join(args.out, "model.pt"))
            flag = "  *saved"
        print(f"epoch {epoch:3d}  train {tr:.5f}  val {va:.5f}  delta {delta:.4f}  TPR {tpr:.4f}  "
              f"TNR {tnr:.4f}  acc {acc:.4f}  AUC {a:.5f}  ({time.time() - t0:.1f}s){flag}")

    with open(os.path.join(args.out, "history.json"), "w") as f:
        json.dump(history, f, indent=1)
    print(f"[done] best epoch {best['epoch']}: val loss {best['val_loss']:.5f}, delta {best['delta']:.4f}, "
          f"accuracy {best['val_acc']:.4f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ep = [h["epoch"] for h in history]
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.5))
    ax[0].plot(ep, [h["train_loss"] for h in history], label="train")
    ax[0].plot(ep, [h["val_loss"] for h in history], label="validation")
    ax[0].set_xlabel("epoch"); ax[0].set_ylabel("BCE loss"); ax[0].set_yscale("log"); ax[0].legend()
    ax[1].plot(ep, [h["val_acc"] for h in history])
    ax[1].set_xlabel("epoch"); ax[1].set_ylabel(r"validation accuracy at $\delta$")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "training.png"), dpi=120)


if __name__ == "__main__":
    main()
