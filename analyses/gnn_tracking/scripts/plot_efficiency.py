#!/usr/bin/env python3
"""Efficiency vs pT, eta and phi: hybrid L1 tracks vs GNN track building with >= 4 and >= 3 layers.

  python3 plot_efficiency.py --sample "QED mu+mu-=DIR" --sample "QED e+e-=DIR" --out ../plots/eff_hybrid_vs_gnn

Each DIR is a compare_hybrid.py output (rows.npz, --denominator tp --count layers). The GNN
curves are the one-crossing-per-layer union-find at that run's --split-threshold, with
candidates needing stubs on >= 4 or >= 3 distinct layers. Error bars: Clopper-Pearson 68%.
--offline DIR:PDG (one per sample, same order) adds the offline efficiency as a line:
primary TPs (tp_ngenpart > 0) with |eta| < 2.4 and |pdgId| = PDG from the
OfflineTrackNtuple files in DIR, found if matched by >= 1 highPurity track
(offline_perf.C definition; different events from the L1 samples).
Writes <out>_<sample>.png/.pdf (1 x 3 panels per sample), <out>.png/.pdf (all samples),
and <out>.csv with the numbers.
"""
import argparse, csv, os, re
import numpy as np
from scipy.stats import beta

PT_BINS = np.array([0.3, 0.4, 0.5, 0.6, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0])
ETA_BINS = np.linspace(-2.4, 2.4, 17)
PHI_BINS = np.linspace(-np.pi, np.pi, 17)
# categorical slots 1-3 of the validated reference palette, in fixed order; markers as a second encoding
SERIES = [("Hybrid L1 tracks", 5, "#2a78d6", "s"),
          ("GNN, ≥ 4 layers", 9, "#eb6834", "o"),
          ("GNN, ≥ 3 layers", 11, "#1baf7a", "^")]
INK, MUTED, GRID = "#1f1f1e", "#5d5c55", "#e6e5df"


def eff(x, found, bins, min_n=10):
    k = np.histogram(x, bins, weights=found)[0]
    n = np.histogram(x, bins)[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        e = k / n
        lo = np.where(k > 0, beta.ppf(0.16, k, n - k + 1), 0.0)
        hi = np.where(k < n, beta.ppf(0.84, k + 1, n - k), 1.0)
    bad = n < min_n
    e[bad] = lo[bad] = hi[bad] = np.nan
    return e, e - lo, hi - e, n, k


OFFLINE_LABEL = "Offline (generalTracks, highPurity)"


def panels(axs, rows, pt_min, title=None, off=None):
    pt, eta, phi = rows[:, 0], rows[:, 1], rows[:, 8]
    hi = pt > pt_min
    specs = [(axs[0], pt, np.ones(len(pt), bool), PT_BINS, r"particle $p_T$ [GeV]", True),
             (axs[1], eta, hi, ETA_BINS, r"particle $\eta$", False),
             (axs[2], phi, hi, PHI_BINS, r"particle $\phi$", False)]
    if off is not None:
        ohi = off[:, 0] > pt_min
        ospecs = [(off[:, 0], np.ones(len(off), bool)), (off[:, 1], ohi), (off[:, 3], ohi)]
    for i, (ax, x, sel, bins, xlab, logx) in enumerate(specs):
        xc = np.sqrt(bins[1:] * bins[:-1]) if logx else 0.5 * (bins[1:] + bins[:-1])
        if off is not None:
            ox, osel = ospecs[i]
            e = eff(ox[osel], off[osel, 2], bins)[0]
            ax.plot(xc, e, "-", color=INK, lw=1.6, label=OFFLINE_LABEL, zorder=2)
        for j, (lab, col, color, mk) in enumerate(SERIES):
            e, el, eh, _, _ = eff(x[sel], rows[sel, col], bins)
            dx = (j - 1) * (0.012 if logx else 0.04) * (1 if not logx else 1)
            xs = xc * (1 + dx) if logx else xc + dx * (bins[1] - bins[0]) / 0.3 * 0.3
            ax.errorbar(xs, e, yerr=[el, eh], fmt=mk, ms=5.5, color=color, mec="white", mew=0.6,
                        elinewidth=1.2, capsize=0, label=lab, zorder=3 + j)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel(xlab, color=INK)
        ax.grid(True, color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(MUTED)
        ax.tick_params(colors=MUTED, labelcolor=INK)
        if logx:
            ax.set_xscale("log")
            ax.set_xticks([0.3, 0.5, 1, 2, 5, 10]); ax.set_xticklabels(["0.3", "0.5", "1", "2", "5", "10"])
            ax.set_xlim(0.28, 11)
        else:
            ax.set_xlim(bins[0], bins[-1])
    axs[0].set_ylabel("tracking efficiency", color=INK)
    axs[1].set_title(rf"$p_T > {pt_min:g}$ GeV", fontsize=10, color=MUTED, loc="left")
    axs[2].set_title(rf"$p_T > {pt_min:g}$ GeV", fontsize=10, color=MUTED, loc="left")
    if title:
        axs[0].set_title(title, fontsize=11.5, color=INK, loc="left", fontweight="bold")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", action="append", required=True, help="label=compare_hybrid output dir")
    ap.add_argument("--pt-min", type=float, default=1.0, help="pT cut for the eta and phi panels")
    ap.add_argument("--note", action="append", default=[], help="per-sample footnote, same order as --sample")
    ap.add_argument("--offline", action="append", default=[], help="DIR:PDG of OfflineTrackNtuple files, per sample")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10.5, "axes.labelsize": 11, "legend.fontsize": 9.5})
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    samples = [(s.split("=", 1)[0], np.load(os.path.join(s.split("=", 1)[1], "rows.npz"))["rows"]) for s in a.sample]
    offs = []
    if a.offline:
        from compare_hybrid import offline_rows
        from gnn_common import list_root_files
        for spec in a.offline:
            d, pdg = spec.rsplit(":", 1)
            offs.append(offline_rows(list_root_files(d), int(pdg)))
    offs += [None] * (len(samples) - len(offs))
    for r in [r for _, r in samples]:
        assert r.shape[1] >= 12, "rows.npz has no GNN >= 3 layers column; rerun compare_hybrid.py"

    def legend(fig):
        h, l = fig.axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper center", ncol=len(h), frameon=False, bbox_to_anchor=(0.5, 1.0))

    def slug(s):
        for a, b in (("μ", "mu"), ("⁺", "p"), ("⁻", "m"), ("+", "p"), ("-", "m")):
            s = s.replace(a, b)
        return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")
    for i, (lab, rows) in enumerate(samples):
        fig, axs = plt.subplots(1, 3, figsize=(13.5, 4.2))
        panels(axs, rows, a.pt_min, lab, offs[i])
        legend(fig)
        note = a.note[i] if i < len(a.note) else ""
        if note:
            fig.text(0.01, 0.005, note, fontsize=8.5, color=MUTED)
        fig.tight_layout(rect=(0, 0.04 if note else 0, 1, 0.93))
        for ext in ("png", "pdf"):
            fig.savefig(f"{a.out}_{slug(lab)}.{ext}", dpi=160)
        plt.close(fig)
    fig, axs = plt.subplots(len(samples), 3, figsize=(13.5, 4.0 * len(samples)), squeeze=False)
    for i, (lab, rows) in enumerate(samples):
        panels(axs[i], rows, a.pt_min, lab, offs[i])
    legend(fig)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext}", dpi=160)
    with open(f"{a.out}.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "variable", "bin_lo", "bin_hi", "n_particles"] + [s[0] for s in SERIES]
                   + ["offline_n_particles", OFFLINE_LABEL])
        for (lab, rows), off in zip(samples, offs):
            hi = rows[:, 0] > a.pt_min
            ohi = off[:, 0] > a.pt_min if off is not None else None
            for j, (var, x, sel, bins) in enumerate([("pT", rows[:, 0], np.ones(len(rows), bool), PT_BINS),
                                                     (f"eta (pT>{a.pt_min:g})", rows[:, 1], hi, ETA_BINS),
                                                     (f"phi (pT>{a.pt_min:g})", rows[:, 8], hi, PHI_BINS)]):
                cols = [eff(x[sel], rows[sel, c], bins)[0] for _, c, _, _ in SERIES]
                n = np.histogram(x[sel], bins)[0]
                if off is not None:
                    ox = off[:, [0, 1, 3][j]]
                    osel = np.ones(len(off), bool) if j == 0 else ohi
                    oe, on = eff(ox[osel], off[osel, 2], bins)[0], np.histogram(ox[osel], bins)[0]
                for b in range(len(bins) - 1):
                    extra = ([int(on[b]), "" if np.isnan(oe[b]) else round(float(oe[b]), 4)] if off is not None else ["", ""])
                    w.writerow([lab, var, round(bins[b], 3), round(bins[b + 1], 3), int(n[b])] +
                               ["" if np.isnan(c[b]) else round(float(c[b]), 4) for c in cols] + extra)
    print("wrote", f"{a.out}.png/.pdf/.csv", "and", ", ".join(f"{a.out}_{slug(l)}.png" for l, _ in samples))
    for lab, rows in samples:
        pt = rows[:, 0]
        line = f"{lab:14s}"
        for lo, hi_ in [(0.3, 0.5), (0.5, 0.75), (0.75, 1), (1, 2), (2, 100)]:
            m = (pt >= lo) & (pt < hi_)
            line += f" | {lo}-{hi_ if hi_ < 100 else '':<4} " + " ".join(f"{rows[m, c].mean():.2f}" for _, c, _, _ in SERIES)
        print(line, "  (hybrid / GNN>=4 / GNN>=3)")


if __name__ == "__main__":
    main()
