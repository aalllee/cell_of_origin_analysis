#!/usr/bin/env python3
"""
Correlate a sample's per-gene FFT scores (output of collapse_fft_wps.py)
against the precomputed Tabula Sapiens average-expression-by-(cell_type,
tissue) reference matrix (output of average_expression_by_celltype_tissue.py).

Writes cell types ranked by correlation strength (ascending -- most
negative correlation first, matching the expected signal direction in the
193-199bp band), one row per (cell_type, tissue) group.

usage: python3 correlate.py <avg_expr_csv> <gene_scores_tsv> -o <output_csv>
"""
import argparse

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "avg_expr_csv", help="path to precomputed avg_expr_by_celltype_tissue.csv"
    )
    parser.add_argument("gene_scores", help="path to this sample's gene_scores.tsv")
    parser.add_argument(
        "-o", "--output", required=True, help="output csv: cell_type, correlation per row"
    )
    args = parser.parse_args()

    avg_expr_df = pd.read_csv(args.avg_expr_csv, index_col=0)
    print(f"reference matrix: {avg_expr_df.shape[0]} genes x {avg_expr_df.shape[1] - 1} groups")

    gene_scores = pd.read_csv(args.gene_scores, sep="\t")
    score_col = gene_scores.columns[1]  # per gene fft scores
    scores = gene_scores.set_index("gene_id")[score_col]

    common_genes = avg_expr_df.index.intersection(scores.index)
    print(f"{len(common_genes)} shared genes")

    avg_expr_common = avg_expr_df.loc[common_genes].drop(columns=["gene_name"])
    scores_common = scores.loc[common_genes]

    correlations = avg_expr_common.corrwith(scores_common)

    corr_df = correlations.reset_index()
    corr_df.columns = ["cell_type", "correlation"]
    corr_df = corr_df.sort_values("correlation")  # ascending -- most negative first

    corr_df.to_csv(args.output, index=False)
    print(f"wrote {len(corr_df)} rows to {args.output}")


if __name__ == "__main__":
    main()
