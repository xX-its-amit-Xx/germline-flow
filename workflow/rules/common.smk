# =============================================================================
# common.smk — helper functions used by multiple rule files.
#
# Kept separate from the rule files because `snakemake --lint` flags any
# .smk file that mixes rule definitions with helper functions (the linter
# wants either lambda-on-the-spot or a dedicated module like this one).
# =============================================================================


def get_fastqs(wildcards):
    """Return the paired FASTQs for a given sample."""
    row = samples_df.loc[wildcards.sample]
    return {"r1": row["fastq_1"], "r2": row["fastq_2"]}


def get_sex(sample: str) -> str:
    """Return reported sex from the sample sheet, or 'U' if absent."""
    if "sex" in samples_df.columns:
        return str(samples_df.loc[sample, "sex"] or "U").upper()
    return "U"


def known_sites_args(_wildcards):
    """Build the `--known-sites ...` flag list for GATK BaseRecalibrator."""
    return " ".join(f"--known-sites {p}" for p in config["reference"]["known_sites"])


def hard_filter_args(kind: str) -> str:
    """Render the GATK VariantFiltration `--filter-name/--filter-expression` flags."""
    parts = []
    for f in config["joint_calling"]["hard_filters"][kind]:
        parts.append(f"--filter-name {f['name']} --filter-expression \"{f['expression']}\"")
    return " ".join(parts)


def sample_map_format_string() -> str:
    """
    Render the GenomicsDBImport sample-name-map as a printf format string with
    literal `\\t` and `\\n` escapes. printf interprets the escapes at runtime.
    """
    return "\\n".join(f"{s}\\t{GVCF}/{s}.g.vcf.gz" for s in SAMPLES)
