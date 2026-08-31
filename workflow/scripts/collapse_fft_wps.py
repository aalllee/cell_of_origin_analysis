#!/usr/bin/env python3
"""Collapse a genes x frequency FFT-WPS matrix (output of convert_files.py)
to one score per gene, by averaging the frequency columns in the 193-199bp
band, the "mean FFT intensity at the 193-199 frequency range" used in
Snyder et al. 2016 and the 2024 cell-of-origin paper to correlate cfDNA
nucleosome periodicity with gene expression.

Columns are matched by numeric value falling in [--low, --high], not by
exact string match since the frequency bins produced by fft_path.R depend on
the FFT resolution of the input data and won't necessarily include the
literal values 193/196/199 as separate columns.

Usage:
    python collapse_fft_wps.py fft_mysample_WPS.tsv.gz -o mysample_gene_scores.tsv
"""
import argparse
import csv
import gzip
import sys


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", help="fft_<sample>_WPS.tsv.gz (or .tsv) from convert_files.py")
    parser.add_argument("-o", "--output", required=True, help="output tsv: gene_id, score (.gz to compress)")
    parser.add_argument("--low", type=float, default=193, help="lower bound of frequency band, bp (default 193)")
    parser.add_argument("--high", type=float, default=199, help="upper bound of frequency band, bp (default 199)")
    args = parser.parse_args()

    with open_maybe_gzip(args.input) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        band_idx = []
        for i, col in enumerate(header[1:], start=1):
            try:
                freq = float(col)
            except ValueError:
                continue
            if args.low <= freq <= args.high:
                band_idx.append(i)

        if not band_idx:
            sys.exit(
                f"no columns fall within [{args.low}, {args.high}] bp -- "
                f"available frequencies: {header[1:]}"
            )
        sys.stderr.write(
            "averaging %d column(s) in [%g, %g]bp: %s\n"
            % (len(band_idx), args.low, args.high, [header[i] for i in band_idx])
        )

        with open_maybe_gzip(args.output, "wt") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["gene_id", "fft_wps_%g_%g" % (args.low, args.high)])
            n_genes = 0
            for row in reader:
                if not row:
                    continue
                gene_id = row[0]
                vals = []
                for i in band_idx:
                    if i >= len(row):
                        continue
                    try:
                        vals.append(float(row[i]))
                    except ValueError:
                        continue
                if not vals:
                    continue
                score = sum(vals) / len(vals)
                writer.writerow([gene_id, score])
                n_genes += 1

        sys.stderr.write("wrote %d genes to %s\n" % (n_genes, args.output))


if __name__ == "__main__":
    main()
