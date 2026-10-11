#!/usr/bin/env python3
"""Overlay offline, several hybrid variants and the GNN on the same events.

  python3 overlay_variants.py --out fig.png --offline-ntuple DIR [--offline-pdg 11] \
      --run "hybrid (v2 cuts)=DIR1" --run "hybrid (ME 6, window x4)=DIR2"

Each DIR is a compare_hybrid.py output (rows.npz) for the same graphs, so the GNN
columns are identical and only the hybrid column differs. The GNN curve (GNN + one
crossing per layer, >= 4) is taken from the first run.
Panels: efficiency vs pT (all eta, and barrel |eta| < 1), vs eta and vs phi for
pT > --pt-min.
"""
import argparse, os
import numpy as np
from compare_hybrid import offline_rows, PT_BINS, ETA_BINS, PHI_BINS
from evaluate import binned_eff
from gnn_common import list_root_files


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="label=compare_hybrid output dir")
    ap.add_argument("--offline-ntuple")
    ap.add_argument("--offline-pdg", type=int, default=0)
    ap.add_argument("--pt-min", type=float, default=1.0)
    ap.add_argument("--title", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    runs = [(r.split("=", 1)[0], np.load(os.path.join(r.split("=", 1)[1], "rows.npz"))["rows"]) for r in a.run]
    # rows: pt, eta, gnn3, gnn4, ideal3, hybrid, n, noclf4, phi, gnn_split4, noclf_split4
    curves = [("GNN + 1 crossing/layer", runs[0][1], 9, "o", "C0")]
    for (lab, rows), (mk, c) in zip(runs, [("s", "C3"), ("D", "C1"), ("v", "C2")]):
        curves.append((lab, rows, 5, mk, c))
    off = None
    if a.offline_ntuple:
        off = offline_rows(list_root_files(a.offline_ntuple), a.offline_pdg)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    xpt = np.sqrt(PT_BINS[1:] * PT_BINS[:-1])
    panels = [(axs[0, 0], "pt", lambda r: np.ones(len(r), bool), "all |eta| < 2.4"),
              (axs[0, 1], "pt", lambda r: np.abs(r[:, 1]) < 1.0, "barrel |eta| < 1"),
              (axs[1, 0], "eta", lambda r: r[:, 0] > a.pt_min, f"pT > {a.pt_min:g} GeV"),
              (axs[1, 1], "phi", lambda r: r[:, 0] > a.pt_min, f"pT > {a.pt_min:g} GeV")]
    for ax, var, sel, title in panels:
        bins, col_x = {"pt": (PT_BINS, 0), "eta": (ETA_BINS, 1), "phi": (PHI_BINS, 8)}[var]
        x = xpt if var == "pt" else 0.5 * (bins[1:] + bins[:-1])
        if off is not None and var != "phi":
            m = sel(off) if var != "eta" else off[:, 0] > a.pt_min
            e, _ = binned_eff(off[m, 0 if var == "pt" else 1], off[m, 2], bins)
            ax.plot(x, e, "-", color="k", lw=1.5, label="offline generalTracks (highPurity)", zorder=0)
        for lab, rows, col, mk, c in curves:
            m = sel(rows)
            e, err = binned_eff(rows[m, col_x], rows[m, col], bins)
            ax.errorbar(x, e, yerr=err, fmt=mk, ms=4, color=c, label=lab)
        ax.set_ylim(0, 1.05); ax.set_ylabel("tracking efficiency"); ax.set_title(title, fontsize=9)
        ax.set_xlabel({"pt": r"particle $p_T$ [GeV]", "eta": r"particle $\eta$", "phi": r"particle $\phi$"}[var])
        if var == "pt":
            ax.set_xscale("log"); ax.set_xticks([0.3, 0.5, 1, 2, 5, 10]); ax.set_xticklabels(["0.3", "0.5", "1", "2", "5", "10"])
        ax.legend(fontsize=8, loc="lower right")
    if a.title:
        fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    fig.savefig(os.path.splitext(a.out)[0] + ".pdf")
    print("wrote", a.out)


if __name__ == "__main__":
    main()
