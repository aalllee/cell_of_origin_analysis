#!/usr/bin/env python3
"""
Correlate a sample's per-gene FFT scores (output of collapse_fft_wps.py)
against a precomputed (cell_type, tissue) reference expression matrix.

Two reference layouts are accepted, detected from the file itself:

  avg_expr_by_celltype_tissue.csv   comma separated, an `ensemblid` index plus a
                                    `gene_name` column, values on a log1p scale
                                    (output of average_expression_by_celltype_tissue.py)
  geneExpression_Stanley.tsv        tab separated, a `feature_id` index, no gene
                                    name column, values are per-cell-type
                                    proportions (every column sums to 1)

Writes cell types ranked by correlation strength (ascending -- most
negative correlation first, matching the expected signal direction in the
193-199bp band), one row per (cell_type, tissue) group.

--transform controls what is correlated:

  none  (default)  correlate the values as they are. Correct for a reference that
                   is already on a log scale, such as avg_expr_by_celltype_tissue.csv.
                   Snyder et al. 2016 (`cfDNA/expression/plots.R`) apply log2, but to
                   Protein Atlas FPKM, a linear per-transcript abundance; logging an
                   already-logged quantity is not the same operation.
  log              log2 the values first. Correct for a linear abundance such as the
                   proportions in geneExpression_Stanley.tsv, where a handful of
                   high-abundance genes otherwise dominate every correlation -- this
                   is the case Snyder et al. were handling.

A centred log ratio is not offered because it would make no difference here: CLR
subtracts the per-column mean of the logs, a constant within each column, and
Pearson correlation is invariant to a shift of either variable. CLR and `log`
give identical correlations against a fixed score vector.

The script warns if a compositional reference is correlated without a transform.

usage: python3 correlate.py <reference> <gene_scores_tsv> -o <output_csv> [--transform log]
"""
import argparse

import numpy as np
import pandas as pd


def read_reference(path):
    """load the reference matrix, sniffing the delimiter and dropping label columns"""
    with open(path) as fh:
        header = fh.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    labels = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
    if labels:
        df = df.drop(columns=labels)
    return df, sep, labels


def is_compositional(df, tol=1e-3):
    """true when every column sums to 1, i.e. the values are proportions"""
    return bool(((df.sum(axis=0) - 1.0).abs() < tol).all())


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "avg_expr_csv", help="reference matrix: avg_expr_by_celltype_tissue.csv or geneExpression_Stanley.tsv"
    )
    parser.add_argument("gene_scores", help="path to this sample's gene_scores.tsv")
    parser.add_argument(
        "-o", "--output", required=True, help="output csv: cell_type, correlation per row"
    )
    parser.add_argument(
        "--transform", choices=["none", "log"], default="none",
        help="transform applied to expression before correlating (default none; use log "
             "for a linear reference such as geneExpression_Stanley.tsv)",
    )
    parser.add_argument(
        "--pseudocount", type=float, default=0.0,
        help="added before log2 when --transform log (default 0; the Stanley reference has "
             "no zeros so none is needed)",
    )
    args = parser.parse_args()

    avg_expr_df, sep, labels = read_reference(args.avg_expr_csv)
    print(f"reference matrix: {avg_expr_df.shape[0]} genes x {avg_expr_df.shape[1]} groups "
          f"({'tab' if sep == chr(9) else 'comma'} separated"
          f"{', dropped ' + ', '.join(labels) if labels else ''})")

    compositional = is_compositional(avg_expr_df)
    if compositional:
        print("reference is compositional (every column sums to 1)")
    if compositional and args.transform == "none":
        print("WARNING: correlating proportions without --transform log. A few high-abundance "
              "genes will dominate every correlation. Pass --transform log unless you mean this.")

    gene_scores = pd.read_csv(args.gene_scores, sep="\t")
    score_col = gene_scores.columns[1]  # per gene fft scores
    scores = gene_scores.set_index("gene_id")[score_col]

    common_genes = avg_expr_df.index.intersection(scores.index)
    print(f"{len(common_genes)} shared genes")
    if len(common_genes) == 0:
        raise SystemExit(
            "no shared genes: check that the reference and the annotation used for the WPS "
            "run describe the same gene IDs (unversioned ENSG, same genome build)"
        )

    avg_expr_common = avg_expr_df.loc[common_genes]
    scores_common = scores.loc[common_genes]

    if args.transform == "log":
        if args.pseudocount == 0.0 and (avg_expr_common.to_numpy() <= 0).any():
            raise SystemExit(
                "--transform log on a reference containing zeros; rerun with --pseudocount"
            )
        avg_expr_common = np.log2(avg_expr_common + args.pseudocount)
        print(f"applied log2 (pseudocount {args.pseudocount})")

    correlations = avg_expr_common.corrwith(scores_common)

    corr_df = correlations.reset_index()
    corr_df.columns = ["cell_type", "correlation"]
    corr_df = corr_df.sort_values("correlation")  # ascending most negative first

    corr_df.to_csv(args.output, index=False)
    print(f"wrote {len(corr_df)} rows to {args.output}")


if __name__ == "__main__":
    main()
