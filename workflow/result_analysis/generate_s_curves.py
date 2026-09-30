#!/usr/bin/env python3
"""
S curves for one sample: correlation between per-gene FFT intensity and cell-type expression, at every period.

This is the Python version of the first figure in cfDNA/expression/plots.R ("allFreq_correlation_plot.pdf"):
for each period column of the FFT-WPS matrix (120-280 bp), correlate the per-gene intensity at that period with
each (cell type, tissue) group's mean expression, then plot correlation against period. Thin grey lines are the
individual groups; bold lines are category medians (blood/immune, the four negative controls, everything else).

Contributing cell types should dip at the nucleosome repeat length (193-199 bp); the controls should stay flat:
  erythrocytes and platelets are anucleate, and retinal / muscle cells sit behind a barrier or barely divide.

Inputs
  results/<sample>/           the pipeline's per-sample results folder; only
                              proj/fft_summaries/fft_<sample>_WPS.tsv.gz is read
  --ref                       Tabula Sapiens mean expression per group (workflow/resources/...csv)
  --compartments              compartment per group, for the blood/immune category

Usage
  python generate_s_curves.py ../results/IH02
  python generate_s_curves.py ../results/IH02 -o scurves_IH02 --highlight "classical monocyte_Blood"
  python generate_s_curves.py ../results/healthy ../results/cancer -o scurves_compare   # + overlay figure
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
WORKFLOW_DIR = HERE.parent

# four negative-control classes: colour, marker-free (lines only), and the cell types they contain
ERY, PLT, BARRIER, LOWTURN = ("erythrocyte (anucleate)", "platelet (anucleate)",
                              "behind blood barrier", "low turnover")
CONTROL_CATS = [ERY, PLT, BARRIER, LOWTURN]
CONTROL_TYPES = {
    ERY: {"erythrocyte"},
    PLT: {"platelet"},
    BARRIER: {"eye photoreceptor cell", "retinal bipolar neuron", "retina horizontal cell",
              "retinal ganglion cell", "retinal pigment epithelial cell"},
    LOWTURN: {"cardiac muscle cell", "slow muscle cell", "fast muscle cell", "cell of skeletal muscle"},
}
TYPE_TO_CAT = {t: k for k, types in CONTROL_TYPES.items() for t in types}
CAT_COLORS = {"blood/immune": "tab:red", "other": "dimgray",
              ERY: "#2a78d6", PLT: "#e87ba4", BARRIER: "#1baf7a", LOWTURN: "#eda100"}


def fft_matrix(sample_dir):
    """genes x periods FFT-WPS matrix written by convert_files.py, for the sample in this folder"""
    sample = sample_dir.name
    d = sample_dir / "proj" / "fft_summaries"
    for name in (f"fft_{sample}_WPS.tsv.gz", f"fft_{sample}_WPS.tsv"):
        if (d / name).exists():
            m = pd.read_csv(d / name, sep="\t", index_col=0)
            m.columns = m.columns.astype(float)
            return sample, m[sorted(m.columns)]
    raise SystemExit(f"no FFT matrix for '{sample}' in {d}")


def categorise(groups, comp):
    def one(g):
        cell_type = g.split("_", 1)[0]
        if cell_type in TYPE_TO_CAT:
            return TYPE_TO_CAT[cell_type]
        return "blood/immune" if comp[g] == "immune" else "other"
    return pd.Series({g: one(g) for g in groups})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample_dir", type=Path, nargs="+", help="one or more results/<sample> folders")
    ap.add_argument("-o", "--output", type=Path, help="output directory (default: <sample_dir>/s_curves)")
    ap.add_argument("--ref", type=Path, default=WORKFLOW_DIR / "resources/avg_expr_by_celltype_tissue.csv",
                    help="mean expression per (cell type, tissue) group")
    ap.add_argument("--compartments", type=Path, default=WORKFLOW_DIR / "resources/celltype_compartment.csv",
                    help="compartment per group (immune/epithelial/endothelial/stromal)")
    ap.add_argument("--band", type=float, nargs=2, default=(193, 199), metavar=("LOW", "HIGH"),
                    help="period band highlighted and summarised (default 193 199)")
    ap.add_argument("--highlight", default=None,
                    help="one group drawn in black with markers, like NB-4 in plots.R")
    args = ap.parse_args()

    ref = pd.read_csv(args.ref, index_col=0).drop(columns=["gene_name"], errors="ignore")
    comp = pd.read_csv(args.compartments, index_col=0)["compartment"].reindex(ref.columns)
    cat = categorise(ref.columns, comp)
    if args.highlight and args.highlight not in ref.columns:
        raise SystemExit(f"--highlight '{args.highlight}' is not a group in the reference")

    all_curves = {}
    for sample_dir in args.sample_dir:
        sample, fft = fft_matrix(sample_dir)
        out_dir = args.output or (sample_dir / "s_curves")
        out_dir.mkdir(parents=True, exist_ok=True)

        shared = ref.index.intersection(fft.index)
        print(f"\n{sample}: {fft.shape[0]} genes x {fft.shape[1]} periods | {len(shared)} shared with the reference")
        curves = pd.DataFrame({p: ref.loc[shared].corrwith(fft.loc[shared, p]) for p in fft.columns}).T
        curves.index.name = "period_bp"
        curves.to_csv(out_dir / f"s_curves_{sample}.csv.gz")
        all_curves[sample] = curves

        periods = curves.index.values
        in_band = (periods >= args.band[0]) & (periods <= args.band[1])
        rows = []
        for k in ["blood/immune", "other"] + CONTROL_CATS:
            med = curves.loc[:, cat[curves.columns] == k].median(axis=1)
            rows.append({"category": k, "groups": int((cat == k).sum()), "r at band": med[in_band].mean(),
                         "min r": med.min(), "period of min r": med.idxmin()})
        summary = pd.DataFrame(rows).set_index("category")
        summary.to_csv(out_dir / f"s_curves_summary_{sample}.csv")
        print(summary.round(3).to_string())

        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.axvspan(*args.band, color="gold", alpha=0.3, zorder=0)
        ax.axhline(0, color="k", lw=0.6)
        ax.plot(periods, curves.values, color="0.8", lw=0.4, alpha=0.5, zorder=1)
        for k in ["blood/immune"] + CONTROL_CATS:
            med = curves.loc[:, cat[curves.columns] == k].median(axis=1)
            ax.plot(periods, med, color=CAT_COLORS[k], lw=2.4, zorder=3, label=f"{k} (n={(cat == k).sum()})")
        ax.plot(periods, curves.loc[:, cat[curves.columns] == "other"].median(axis=1), color=CAT_COLORS["other"],
                lw=1.6, ls=":", zorder=2, label=f"other (n={(cat == 'other').sum()})")
        if args.highlight:
            ax.plot(periods, curves[args.highlight], color="k", lw=2, marker="o", ms=3.5, zorder=4,
                    label=args.highlight)
        ax.set_xlabel("period (bp)"); ax.set_ylabel("correlation with expression (r)")
        ax.set_title(f"{sample}: correlation of FFT-WPS intensity with expression, per period\n"
                     f"thin grey = each of the {len(cat)} (cell type, tissue) groups, bold = category median",
                     fontsize=11)
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
        plt.tight_layout(rect=(0, 0, 0.82, 1))
        fig.savefig(out_dir / f"s_curves_{sample}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out_dir / f's_curves_{sample}.png'}")

    # with several samples: one overlay of the category medians, line style per sample
    if len(all_curves) > 1:
        out_dir = args.output or args.sample_dir[0].parent
        out_dir.mkdir(parents=True, exist_ok=True)
        styles = ["-", "--", ":", "-."]
        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.axvspan(*args.band, color="gold", alpha=0.3, zorder=0)
        ax.axhline(0, color="k", lw=0.6)
        for (sample, curves), ls in zip(all_curves.items(), styles):
            for k in ["blood/immune", "other"] + CONTROL_CATS:
                med = curves.loc[:, cat[curves.columns] == k].median(axis=1)
                ax.plot(curves.index.values, med, color=CAT_COLORS[k], lw=2.2, ls=ls, alpha=0.9)
        handles = ([Line2D([], [], color=CAT_COLORS[k], lw=2.2, label=k)
                    for k in ["blood/immune", "other"] + CONTROL_CATS]
                   + [Line2D([], [], color="k", lw=2.2, ls=ls, label=s)
                      for s, ls in zip(all_curves, styles)])
        ax.legend(handles=handles, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
        ax.set_xlabel("period (bp)"); ax.set_ylabel("median correlation with expression (r)")
        ax.set_title("Category medians, " + " vs ".join(all_curves), fontsize=12)
        plt.tight_layout(rect=(0, 0, 0.82, 1))
        png = out_dir / ("s_curves_compare_" + "_".join(all_curves) + ".png")
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f"\nwrote {png}")


if __name__ == "__main__":
    main()
