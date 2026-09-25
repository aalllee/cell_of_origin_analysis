#!/usr/bin/env python3
"""
Build the Tabula Sapiens average expression reference matrix: average
normalized expression per (cell_type, tissue_in_publication) group,
restricted to 10x 3' v3 assay cells and excluding germ cells,
matching the methodology in Stanley et al. 2024 (Nat Commun).

usage: python3 average_expression_by_celltype_tissue.py <h5ad_path> <output_csv>
"""
import sys

import numpy as np
import pandas as pd
import anndata as ad


def main():
    if len(sys.argv) != 3:
        sys.exit(
            "usage: python3 average_expression_by_celltype_tissue.py "
            "<h5ad_path> <output_csv>"
        )
    h5ad_path, out_csv = sys.argv[1], sys.argv[2]

    
    adata = ad.read_h5ad(h5ad_path, backed="r")

    #filter to 10x 3' v3 assay cells and exlude germ cells
    mask = (
        (adata.obs["assay"] == "10x 3' v3")
        & (~adata.obs["cell_type"].str.contains("sperm", case=False, na=False))
    )
    adata_filtered = adata[mask, :]
    print(f"{mask.sum()} / {len(mask)} cells kept after assay/sperm filtering")

    cell_groups = adata_filtered.obs.groupby(
        ["cell_type", "tissue_in_publication"], observed=True
    ).indices
    print(f"{len(cell_groups)} (cell_type, tissue) groups found")

    result = {}
    for i, ((cell_type, tissue), idx) in enumerate(cell_groups.items(), 1):
        idx = np.sort(idx)
        chunk = adata_filtered[idx, :].X  # reads only this group's rows from disk
        mean_expr = (
            np.asarray(chunk.mean(axis=0)).ravel()
            if hasattr(chunk, "toarray")
            else np.asarray(chunk).mean(axis=0)
        )
        col_name = f"{cell_type}_{tissue}"
        result[col_name] = mean_expr
        if i % 50 == 0:
            print(f"[{i}/{len(cell_groups)}] groups done", flush=True)

    avg_expr_df = pd.DataFrame(result, index=adata_filtered.var_names)
    avg_expr_df.insert(0, "gene_name", adata_filtered.var["feature_name"].values)
    avg_expr_df.to_csv(out_csv)
    print(
        f"Wrote {avg_expr_df.shape[0]} genes x {avg_expr_df.shape[1] - 1} "
        f"groups to {out_csv}"
    )


if __name__ == "__main__":
    main()
