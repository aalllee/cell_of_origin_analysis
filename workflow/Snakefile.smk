import os
import sys

# bam_dir is required, passed via --config, e.g.:
# snakemake --config bam_dir="BAM DATA PATH"
# annotation defaults to the real gene annotation, override the same way if needed
if "bam_dir" not in config:
    sys.exit(
        "ERROR: bam_dir is required. Run with:\n"
        "  snakemake --use-conda -c4 -p --config bam_dir=<path>"
    )
BAM_DIR = config["bam_dir"]
ANNOTATION = config.get("annotation", "../annotation/transcriptAnno_GRCh38.body.tsv")

EXTRACT_SCRIPT = f"../cfDNA/expression/extractReadStartsFromBAM_Region_WPS.py"
FFT_SCRIPT = f"../cfDNA/expression/fft_path.R"
CONVERT_SCRIPT = f"../cfDNA/expression/convert_files.py"
COLLAPSE_SCRIPT = "scripts/collapse_fft_wps.py"
AVG_EXPR_SCRIPT = "scripts/average_expression_by_celltype_tissue.py"
CORRELATE_SCRIPT = "scripts/correlate.py"

# tabula_h5ad points at the downloaded Tabula Sapiens h5ad, override via --config
# snakemake --config bam_dir=<path> tabula_h5ad=<path>
TABULA_H5AD = config.get("tabula_h5ad", "../data/tabula_sapiens.h5ad")




#return a list of matching bam file names [sample1, sample2,...]
SAMPLES, = glob_wildcards(os.path.join(BAM_DIR, "{sample}.bam"))


#set custom directories, dispatch individual jobs per sample
rule all:
    input:
        expand("results/{sample}/correlation.csv", sample=SAMPLES),


# flag PCR duplicates (same outer fragment coordinates) before WPS, like Picard in Stanley et al. 2024;
# the WPS script skips reads *flagged* as duplicates, and input BAMs may carry no flags
rule mark_duplicates:
    input:
        os.path.join(BAM_DIR, "{sample}.bam"),
    output:
        bam="results/{sample}/{sample}.markdup.bam",
        bai="results/{sample}/{sample}.markdup.bam.bai",
        stats="results/{sample}/{sample}.markdup_stats.txt",
    threads: 4
    conda:
        "envs/cfdna.yaml"
    shell:
        """
        samtools collate -@ {threads} -O -u {input} {output.bam}.collate \
          | samtools fixmate -@ {threads} -m -u - - \
          | samtools sort -@ {threads} -u -T {output.bam}.sort - \
          | samtools markdup -@ {threads} -s -f {output.stats} -T {output.bam}.markdup - {output.bam}
        samtools index {output.bam}
        """


# the Python 2 scripts are called via $CONDA_PREFIX so an active venv/pyenv python can't shadow the env's python
rule extract_wps:
    input:
        bam="results/{sample}/{sample}.markdup.bam",
        bai="results/{sample}/{sample}.markdup.bam.bai",
        annotation=ANNOTATION,
    output:
        directory("results/{sample}/body"),
    params:
        script=EXTRACT_SCRIPT,
        min_insert=120,
        max_insert=180,
    conda:
        "envs/cfdna.yaml"
    shell:
        """
        mkdir -p {output}
        "$CONDA_PREFIX/bin/python" {params.script} \
          --minInsert={params.min_insert} --maxInsert={params.max_insert} \
          -i {input.annotation} \
          -o '{output}/block_%s.tsv.gz' \
          {input.bam}
        """


rule fft:
    input:
        "results/{sample}/body",
    output:
        directory("results/{sample}/proj/{sample}/fft"),
    params:
        script=FFT_SCRIPT,
    shell:
        """
        mkdir -p {output}
        ( cd {input}; ls block_*.tsv.gz ) | xargs -n 500 Rscript {params.script} $(pwd)/{input} $(pwd)/{output}
        """


rule convert_summary:
    input:
        fft="results/{sample}/proj/{sample}/fft",
        annotation=ANNOTATION,
    output:
        wps="results/{sample}/proj/fft_summaries/fft_{sample}_WPS.tsv.gz",
        cov="results/{sample}/proj/fft_summaries/fft_{sample}_cov.tsv.gz",
        starts="results/{sample}/proj/fft_summaries/fft_{sample}_starts.tsv.gz",
    params:
        script=CONVERT_SCRIPT,
    conda:
        "envs/cfdna.yaml"
    shell:
        """
        mkdir -p results/{wildcards.sample}/proj/fft_summaries
        "$CONDA_PREFIX/bin/python" {params.script} \
          -a {input.annotation} \
          -t results/{wildcards.sample} \
          -r results/{wildcards.sample} \
          -p proj \
          -i {wildcards.sample}
        """


rule collapse_scores:
    input:
        "results/{sample}/proj/fft_summaries/fft_{sample}_WPS.tsv.gz",
    output:
        "results/{sample}/gene_scores.tsv",
    params:
        script=COLLAPSE_SCRIPT,
        low=193,
        high=199,
    shell:
        "python3 {params.script} {input} -o {output} --low {params.low} --high {params.high}"


# built once. reused as input by every correlate job below
rule build_avg_expression:
    output:
        "resources/avg_expr_by_celltype_tissue.csv",
    params:
        script=AVG_EXPR_SCRIPT,
        h5ad=TABULA_H5AD,
    shell:
        "python3 {params.script} {params.h5ad} {output}"


rule correlate:
    input:
        avg_expr="resources/avg_expr_by_celltype_tissue.csv",
        scores="results/{sample}/gene_scores.tsv",
    output:
        "results/{sample}/correlation.csv",
    params:
        script=CORRELATE_SCRIPT,
    shell:
        "python3 {params.script} {input.avg_expr} {input.scores} -o {output}"