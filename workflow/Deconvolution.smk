import os
import sys

# bam_dir is required, passed via --config, e.g.:
# snakemake --config bam_dir="BAM DATA PATH"
# annotation and expression default to the Ensembl canonical GRCh38 pair below,
# override either the same way if needed, e.g. for a GRCh37 sample:
#   --config bam_dir=<path> annotation=resources/annotation/transcriptAnno-GRCh37.75.body.tsv
if "bam_dir" not in config:
    sys.exit(
        "ERROR: bam_dir is required. Run with:\n"
        "  snakemake --use-conda -c4 -p --config bam_dir=<path>"
    )
BAM_DIR = config["bam_dir"]

# the annotation and the expression reference are a matched pair: all 19,620 genes of
# geneExpression_Stanley.tsv are present in Ensemble_canonical_GRCh38.body.tsv. The older
# transcriptAnno_GRCh38.body.tsv misses 588 of them, so do not mix the two generations.
ANNOTATION = config.get("annotation", "resources/annotation/Ensemble_canonical_GRCh38.body.tsv")
EXPRESSION = config.get("expression", "resources/geneExpression_Stanley.tsv")

# geneExpression_Stanley.tsv holds per-cell-type proportions (every column sums to 1), a
# linear abundance, so it is logged before correlating. Set expr_transform=none for a
# reference that is already on a log scale, such as avg_expr_by_celltype_tissue.csv.
EXPR_TRANSFORM = config.get("expr_transform", "log")

EXTRACT_SCRIPT = f"../cfDNA/expression/extractReadStartsFromBAM_Region_WPS.py"
FFT_SCRIPT = f"../cfDNA/expression/fft_path.R"
CONVERT_SCRIPT = f"../cfDNA/expression/convert_files.py"
COLLAPSE_SCRIPT = "scripts/collapse_fft_wps.py"
CORRELATE_SCRIPT = "scripts/correlate.py"




#return a list of matching bam file names [sample1, sample2,...]
SAMPLES, = glob_wildcards(os.path.join(BAM_DIR, "{sample}.bam"))

localrules: all

#set custom directories, dispatch individual jobs per sample
rule all:
    input:
        expand("results/{sample}/correlation.csv", sample=SAMPLES),


# runs directly on the input BAMs (no duplicate marking); they must already be indexed
# the Python 2 scripts are called via $CONDA_PREFIX so an active venv/pyenv python can't shadow the env's python
rule extract_wps:
    input:
        bam=os.path.join(BAM_DIR, "{sample}.bam"),
        annotation=ANNOTATION,
    output:
        directory("results/{sample}/body"),
    threads: 1
    resources:
        time=180,
    params:
        sample = "{sample}",
        proj = "Project",
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
    threads: 1
    resources:
        time=180,
    params:
        sample = "{sample}",
        proj = "Project",
        script=FFT_SCRIPT,
    conda:
        "envs/r.yaml"
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
    threads: 1
    resources:
        time=60,
    params:
        sample = "{sample}",
        proj = "Project",
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
    threads: 1
    resources:
        time=60,
    params:
        sample = "{sample}",
        proj = "Project",
        script=COLLAPSE_SCRIPT,
        low=193,
        high=199,
    conda:
        "envs/python.yaml"
    shell:
        "python3 {params.script} {input} -o {output} --low {params.low} --high {params.high}"


# EXPRESSION is a shipped file, not a built one: the old build_avg_expression rule derived
# resources/avg_expr_by_celltype_tissue.csv from the Tabula Sapiens h5ad and is gone with it.
rule correlate:
    input:
        avg_expr=EXPRESSION,
        scores="results/{sample}/gene_scores.tsv",
    output:
        "results/{sample}/correlation.csv",
    threads: 1
    resources:
        time=60,
    params:
        sample = "{sample}",
        proj = "Project",
        script=CORRELATE_SCRIPT,
        transform=EXPR_TRANSFORM,
    conda:
        "envs/python.yaml"
    shell:
        "python3 {params.script} {input.avg_expr} {input.scores} -o {output}"
        " --transform {params.transform}"