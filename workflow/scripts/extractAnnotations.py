"""
Extract protein-coding genes from a GENCODE GTF, one row per gene.

Keeps each gene's Ensembl canonical transcript when it also carries a CCDS
tag, and writes the transcript span as tab-separated rows (no header):
    gene_id    chrom    start    end    strand

Usage:
    python extractAnnotations.py gencode.v50.primary_assembly.annotation.gtf.gz -o transcriptAnno-GRCh38.tsv --min-length 10000
"""
import argparse
import gzip
import re
import pandas as pd

ATTR_RE = re.compile(r'(\w+) "([^"]*)"')


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def parse_attrs(field):
    return dict(ATTR_RE.findall(field))


def parse_tags(field):
    return set(m for k, m in ATTR_RE.findall(field) if k == "tag")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("annot_gtf", help="path to GENCODE annotation file (.gtf or .gtf.gz)")
    parser.add_argument("-o", "--output", required=True, help="output tsv: gene_id, chrom, start, end, strand")
    parser.add_argument("--min-length", default=0, type=int,
                        help="keep only genes at least this many bp long, start to end (default 0 = keep all)")
    args = parser.parse_args()

    annotation_rows = []

    with open_maybe_gzip(args.annot_gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")

            if fields[2] != "transcript":
                continue

            chrom, _, _, start, end, _, strand, _, attr_field = fields

            attrs = parse_attrs(attr_field)
            tags = parse_tags(attr_field)

            if attrs.get("gene_type") != "protein_coding":
                continue
            if "CCDS" not in tags or "Ensembl_canonical" not in tags:  #"Ensembl_canonical" not in tags
                continue

            start, end = int(start), int(end)
            if end - start + 1 < args.min_length:
                continue

            annotation_rows.append({
                "gene_id": attrs["gene_id"].split(".")[0],  # drop version: ENSG00000186092.7 -> ENSG00000186092
                "chrom": chrom.removeprefix("chr"),          # chr1 -> 1; drop this if your BAMs use chr names
                "start": start,
                "end": end,
                "strand": strand,
            })

    annotation_df = pd.DataFrame(annotation_rows).sort_values("gene_id").reset_index(drop=True)
    annotation_df.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"wrote {len(annotation_df)} genes to {args.output}")


if __name__ == "__main__":
    main()
