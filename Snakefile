# =============================================================================
# germline-flow Snakefile
# Top-level entry point. Includes the per-stage rule files under workflow/rules.
# Run with:  snakemake --use-conda --cores N
#       or:  snakemake --profile profiles/slurm --use-conda
# =============================================================================

from pathlib import Path
import pandas as pd

# ---- config ------------------------------------------------------------------

configfile: "config/config.yaml"

SAMPLES_TSV = config["samples"]
samples_df = pd.read_csv(SAMPLES_TSV, sep="\t", dtype=str).set_index("sample", drop=False)
SAMPLES = samples_df["sample"].tolist()

# Sanity checks at parse time — fail loudly, fail early.
required_cols = {"sample", "fastq_1", "fastq_2"}
missing = required_cols - set(samples_df.columns)
if missing:
    raise ValueError(
        f"sample sheet {SAMPLES_TSV} is missing required columns: {sorted(missing)}"
    )
if samples_df.index.duplicated().any():
    dups = samples_df.index[samples_df.index.duplicated()].tolist()
    raise ValueError(f"duplicate sample IDs in {SAMPLES_TSV}: {dups}")

# Result directories — kept short for cluster log readability.
RESULTS  = Path("results")
ALIGN    = RESULTS / "align"
DEDUP    = RESULTS / "dedup"
BQSR     = RESULTS / "bqsr"
GVCF     = RESULTS / "gvcf"
JOINT    = RESULTS / "joint"
QC       = RESULTS / "qc"
PLINK    = RESULTS / "plink"
REPORT   = RESULTS / "report"

# ---- include rule modules ----------------------------------------------------
# common.smk must be included first because every other rule file references
# its helper functions.

include: "workflow/rules/common.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/markdup.smk"
include: "workflow/rules/bqsr.smk"
include: "workflow/rules/callvariants.smk"
include: "workflow/rules/jointgeno.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/export.smk"

# ---- default target ----------------------------------------------------------

rule all:
    """
    The 'all' target asks for the final artefacts: the joint VCF, the PLINK2
    export, and the rendered HTML cohort report. Everything else is pulled in
    as a transitive dependency.
    """
    input:
        JOINT  / "cohort.filtered.vcf.gz",
        PLINK  / "cohort.qc.pgen",
        PLINK  / "cohort.qc.pvar",
        PLINK  / "cohort.qc.psam",
        REPORT / "cohort_qc.html",
