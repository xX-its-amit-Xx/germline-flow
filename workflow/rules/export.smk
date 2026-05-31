# =============================================================================
# export.smk — PLINK2 export, sample/variant QC, PCA, KING relatedness, report.
#
# Outputs the analysis-ready PLINK2 fileset (.pgen/.pvar/.psam) plus the
# auxiliary tables consumed by the HTML report:
#   - PCA eigenvectors / eigenvalues (ancestry, outlier detection)
#   - KING-robust kinship (cryptic relatedness, sample swaps)
# =============================================================================

rule plink_import:
    """Convert the joint-filtered VCF into PLINK2 pgen format."""
    input:
        vcf = JOINT / "cohort.filtered.vcf.gz",
    output:
        pgen = PLINK / "cohort.raw.pgen",
        pvar = PLINK / "cohort.raw.pvar",
        psam = PLINK / "cohort.raw.psam",
    log:
        "logs/plink/import.log",
    params:
        prefix = lambda wc, output: str(PLINK / "cohort.raw"),
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb  = 8000,
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        plink2 \
            --vcf {input.vcf} \
            --vcf-filter \
            --vcf-half-call missing \
            --double-id \
            --max-alleles 2 \
            --make-pgen \
            --out {params.prefix} \
            2> {log}
        """


rule plink_qc_filter:
    """
    Apply sample/variant QC: call rate, MAF, HWE. Thresholds in config.yaml.
    These cuts are reasonable defaults; tune per study (rare-variant cohorts
    will want maf_min: 0 and a separate gnomAD-based filter, for example).
    """
    input:
        pgen = PLINK / "cohort.raw.pgen",
        pvar = PLINK / "cohort.raw.pvar",
        psam = PLINK / "cohort.raw.psam",
    output:
        pgen = PLINK / "cohort.qc.pgen",
        pvar = PLINK / "cohort.qc.pvar",
        psam = PLINK / "cohort.qc.psam",
    log:
        "logs/plink/qc_filter.log",
    params:
        in_prefix  = lambda wc, input: str(input.pgen)[:-5],
        out_prefix = lambda wc, output: str(output.pgen)[:-5],
        call_rate  = config["plink_export"]["call_rate_min"],
        maf        = config["plink_export"]["maf_min"],
        hwe        = config["plink_export"]["hwe_p_min"],
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb  = 8000,
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        # Compute mind/geno as (1 - call_rate)
        mind_geno=$(awk -v c={params.call_rate} 'BEGIN{{printf "%.6f", 1 - c}}')
        plink2 \
            --pfile {params.in_prefix} \
            --mind ${{mind_geno}} \
            --geno ${{mind_geno}} \
            --maf {params.maf} \
            --hwe {params.hwe} \
            --make-pgen \
            --out {params.out_prefix} \
            2> {log}
        """


rule plink_pca:
    """PCA on QC'd genotypes. Used for ancestry / outlier visualisation."""
    input:
        pgen = PLINK / "cohort.qc.pgen",
        pvar = PLINK / "cohort.qc.pvar",
        psam = PLINK / "cohort.qc.psam",
    output:
        eigenvec = PLINK / "cohort.qc.eigenvec",
        eigenval = PLINK / "cohort.qc.eigenval",
    log:
        "logs/plink/pca.log",
    params:
        in_prefix  = lambda wc, input: str(input.pgen)[:-5],
        out_prefix = lambda wc, output: str(output.eigenvec).rsplit(".", 1)[0],
        ncomp      = config["plink_export"]["pca_components"],
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb  = 8000,
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        # --bad-freqs lets PLINK2 impute allele frequencies from < 50 samples.
        # It's a no-op on real cohorts (>=50 samples); only matters for tiny
        # test sets where imputed MAFs are necessarily inexact anyway.
        plink2 \
            --pfile {params.in_prefix} \
            --pca {params.ncomp} \
            --bad-freqs \
            --out {params.out_prefix} \
            2> {log}
        """


rule plink_king:
    """KING-robust kinship — flags cryptic relatedness & sample swaps."""
    input:
        pgen = PLINK / "cohort.qc.pgen",
        pvar = PLINK / "cohort.qc.pvar",
        psam = PLINK / "cohort.qc.psam",
    output:
        king_in  = PLINK / "cohort.qc.king.cutoff.in.id",
        king_out = PLINK / "cohort.qc.king.cutoff.out.id",
        kinship  = PLINK / "cohort.qc.kin0",
    log:
        "logs/plink/king.log",
    params:
        in_prefix  = lambda wc, input: str(input.pgen)[:-5],
        out_prefix = lambda wc, input: str(input.pgen)[:-5],
        cutoff     = config["plink_export"]["king_cutoff"],
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb  = 8000,
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        # --king-cutoff writes <prefix>.king.cutoff.in.id and .king.cutoff.out.id
        plink2 \
            --pfile {params.in_prefix} \
            --king-cutoff {params.cutoff} \
            --out {params.out_prefix} 2> {log}
        # Full pairwise table (small cohorts only — N^2). Emits <prefix>.kin0.
        plink2 \
            --pfile {params.in_prefix} \
            --make-king-table \
            --out {params.out_prefix} 2>> {log}
        # PLINK2 may emit empty .kin0 (no pairs above its internal threshold);
        # ensure the file exists so downstream rules don't fail on a stat.
        touch {output.kinship}
        """


rule relatedness_report_table:
    """Post-process the PLINK KING output into a tidy table for the report."""
    input:
        kinship = PLINK / "cohort.qc.kin0",
        psam    = PLINK / "cohort.qc.psam",
    output:
        tsv  = PLINK / "cohort.relatedness.tsv",
        json = PLINK / "cohort.relatedness.json",
    log:
        "logs/plink/relatedness_table.log",
    params:
        cutoff = config["plink_export"]["king_cutoff"],
    conda:
        "../envs/qc.yaml"
    resources:
        mem_mb = 2000,
        runtime = 15,
    script:
        "../scripts/relatedness.py"


rule report_html:
    """
    Build the cohort QC HTML report (plotly + Jinja2). All inputs are tables
    or JSON manifests produced by upstream rules; the script never touches BAMs.
    """
    input:
        cohort_tsv   = QC / "cohort_qc.tsv",
        cohort_json  = QC / "cohort_qc.json",
        eigenvec     = PLINK / "cohort.qc.eigenvec",
        eigenval     = PLINK / "cohort.qc.eigenval",
        related_tsv  = PLINK / "cohort.relatedness.tsv",
        related_json = PLINK / "cohort.relatedness.json",
        cohort_sum   = QC / "cohort" / "cohort.variant_calling_summary_metrics",
    output:
        html = REPORT / "cohort_qc.html",
    log:
        "logs/report/cohort_qc.log",
    params:
        title    = "germline-flow cohort QC report",
        samples  = SAMPLES,
        thresholds = config["qc"],
    conda:
        "../envs/report.yaml"
    resources:
        mem_mb = 2000,
        runtime = 15,
    script:
        "../scripts/make_report.py"
