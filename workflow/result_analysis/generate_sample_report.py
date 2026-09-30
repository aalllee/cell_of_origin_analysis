#!/usr/bin/env python3
"""
Minimal per-sample report: how the cell types rank, and how cleanly the categories separate.

For each results/<sample> folder it scores every gene as the mean FFT-WPS intensity over the 193-199 bp band
(exactly collapse_fft_wps.py), correlates that with the mean expression of each Tabula Sapiens (cell type,
tissue) group (Pearson over all shared genes = the pipeline's correlate.py), and writes:

  ranking_<sample>.png     strongest 10 and weakest 10 groups, coloured immune / other / negative control
  spread_<sample>.png      the whole distribution of correlations, split by those three labels
  report_<sample>.csv      the 447 correlations with their labels
  summary.csv              one row per sample (printed as well)

Given several sample folders it adds a side-by-side comparison of the three distributions, which is the quick
way to see how much cleaner a deep sample separates than a shallow one.

Negative controls: erythrocytes and platelets (anucleate, all tissues), retinal cells (behind the blood-retina
barrier) and cardiac / skeletal muscle (post-mitotic). None of them should rank as a strong contributor.

Usage
  python generate_sample_report.py ../results/IH02
  python generate_sample_report.py ../results/IH02 ../results/healthy -o report_out
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
WORKFLOW_DIR = HERE.parent

CONTROL_TYPES = {"erythrocyte", "platelet",                                        # anucleate
                 "eye photoreceptor cell", "retinal bipolar neuron", "retina horizontal cell",
                 "retinal ganglion cell", "retinal pigment epithelial cell",       # behind the blood-retina barrier
                 "cardiac muscle cell", "slow muscle cell", "fast muscle cell",
                 "cell of skeletal muscle"}                                        # post-mitotic
LABELS = ["immune", "other", "negative control"]
COLORS = {"immune": "tab:red", "other": "tab:gray", "negative control": "black"}
INK = "#333333"


def fft_matrix(sample_dir):
    sample = sample_dir.name
    d = sample_dir / "proj" / "fft_summaries"
    for name in (f"fft_{sample}_WPS.tsv.gz", f"fft_{sample}_WPS.tsv"):
        if (d / name).exists():
            m = pd.read_csv(d / name, sep="\t", index_col=0)
            m.columns = m.columns.astype(float)
            return sample, m[sorted(m.columns)]
    raise SystemExit(f"no FFT matrix for '{sample}' in {d}")


def labels_for(groups, comp):
    def one(g):
        if g.split("_", 1)[0] in CONTROL_TYPES:
            return "negative control"
        return "immune" if comp[g] == "immune" else "other"
    return pd.Series({g: one(g) for g in groups})


def ranking_figure(r, lab, sample, out_png, n=10):
    sel = pd.concat([r.iloc[:n], r.iloc[-n:]])
    shown = list(range(1, n + 1)) + list(range(len(r) - n + 1, len(r) + 1))
    y = np.arange(len(sel))[::-1].astype(float)
    y[n:] -= 0.8

    fig, ax = plt.subplots(figsize=(8.5, 0.34 * len(sel) + 1.8))
    ax.barh(y, sel.values, color=[COLORS[lab[g]] for g in sel.index], alpha=0.9, height=0.75)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{k}. {g.split('_', 1)[0]} ({g.split('_', 1)[1].replace('_', ' ')})"
                        for k, g in zip(shown, sel.index)], fontsize=8)
    longest = sel.abs().max()
    for yi, v in zip(y, sel.values):
        if abs(v) >= 0.35 * longest:
            ax.annotate(f"{v:+.3f}", (v, yi), xytext=(6 if v < 0 else -6, 0), textcoords="offset points",
                        va="center", ha="left" if v < 0 else "right", fontsize=8, color="white", fontweight="bold")
        else:
            ax.annotate(f"{v:+.3f}", (max(v, 0.0), yi), xytext=(7, 0), textcoords="offset points",
                        va="center", ha="left", fontsize=8, color=INK, clip_on=False)
    ax.axvline(0, color="k", lw=0.7)
    ax.axhline(y[n - 1] - 0.9, color="0.8", lw=1, ls="--")
    pad = 0.10 * (sel.max() - sel.min())
    ax.set_xlim(sel.min() - pad, sel.max() + pad)
    ax.set_xlabel("correlation with expression (r)")
    ax.set_title(f"{sample}: strongest {n} and weakest {n} of {len(r)} groups", fontsize=11)
    ax.legend(handles=[Patch(color=c, label=k) for k, c in COLORS.items()], fontsize=8,
              loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout(rect=(0, 0, 0.88, 1))
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def spread_panel(ax, r, lab, title):
    rng = np.random.default_rng(0)
    for i, k in enumerate(LABELS):
        vals = r[lab[r.index] == k].values
        bp = ax.boxplot(vals, positions=[i], widths=0.55, patch_artist=True, showfliers=False)
        bp["boxes"][0].set(facecolor=COLORS[k], alpha=0.35)
        bp["medians"][0].set(color=COLORS[k], lw=2)
        ax.scatter(rng.normal(i, 0.07, len(vals)), vals, s=7, color=COLORS[k], alpha=0.45, zorder=3)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(range(len(LABELS)))
    ax.set_xticklabels([f"{k}\n(n={int((lab == k).sum())})" for k in LABELS], fontsize=9)
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_dir", type=Path, nargs="+", help="one or more results/<sample> folders")
    ap.add_argument("-o", "--output", type=Path, help="output directory (default: <sample_dir>/report)")
    ap.add_argument("--ref", type=Path, default=WORKFLOW_DIR / "resources/avg_expr_by_celltype_tissue.csv")
    ap.add_argument("--compartments", type=Path, default=WORKFLOW_DIR / "resources/celltype_compartment.csv")
    ap.add_argument("--band", type=float, nargs=2, default=(193, 199), metavar=("LOW", "HIGH"))
    args = ap.parse_args()

    ref = pd.read_csv(args.ref, index_col=0).drop(columns=["gene_name"], errors="ignore")
    comp = pd.read_csv(args.compartments, index_col=0)["compartment"].reindex(ref.columns)
    lab = labels_for(ref.columns, comp)

    results, rows = {}, []
    for sample_dir in args.sample_dir:
        sample, fft = fft_matrix(sample_dir)
        out_dir = args.output or (sample_dir / "report")
        out_dir.mkdir(parents=True, exist_ok=True)

        band_cols = [c for c in fft.columns if args.band[0] <= c <= args.band[1]]
        score = fft[band_cols].mean(axis=1)
        shared = ref.index.intersection(score.index)
        r = ref.loc[shared].corrwith(score.loc[shared]).sort_values()
        results[sample] = r

        pd.DataFrame({"correlation": r, "label": lab[r.index], "compartment": comp[r.index]}) \
            .to_csv(out_dir / f"report_{sample}.csv")
        ranking_figure(r, lab, sample, out_dir / f"ranking_{sample}.png")

        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        spread_panel(ax, r, lab, f"{sample}: correlations by label")
        ax.set_ylabel("correlation with expression (r)")
        plt.tight_layout()
        fig.savefig(out_dir / f"spread_{sample}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

        med = {k: r[lab[r.index] == k].median() for k in LABELS}
        top20 = lab[r.index[:20]]
        rows.append({"sample": sample, "genes": len(shared), "median gene score": round(score.median(), 1),
                     "strongest r": r.min(), "top group": r.index[0],
                     "median r immune": med["immune"], "median r other": med["other"],
                     "median r control": med["negative control"],
                     "immune - other": med["immune"] - med["other"],
                     "immune - control": med["immune"] - med["negative control"],
                     "immune in top 20": int((top20 == "immune").sum()),
                     "controls in top 20": int((top20 == "negative control").sum())})
        print(f"wrote {out_dir / f'ranking_{sample}.png'} and {out_dir / f'spread_{sample}.png'}")

    summary = pd.DataFrame(rows).set_index("sample")
    out_dir = args.output or args.sample_dir[0].parent
    summary.to_csv(out_dir / "summary.csv", float_format="%.4f")
    print("\n" + summary.round(3).to_string())

    if len(results) > 1:     # side by side: a deep sample should separate the labels far more cleanly
        fig, axes = plt.subplots(1, len(results), figsize=(4.6 * len(results), 4.6), sharey=True)
        for ax, (sample, r) in zip(np.atleast_1d(axes), results.items()):
            spread_panel(ax, r, lab, sample)
        np.atleast_1d(axes)[0].set_ylabel("correlation with expression (r)")
        fig.suptitle("Correlation spread by label, per sample", fontsize=12)
        plt.tight_layout()
        png = out_dir / ("spread_compare_" + "_".join(results) + ".png")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f"\nwrote {png}")


if __name__ == "__main__":
    main()
