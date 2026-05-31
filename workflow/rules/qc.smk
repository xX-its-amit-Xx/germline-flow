# =============================================================================
# qc.smk — Per-sample and cohort-level QC.
#
# Per sample:
#   - mosdepth: coverage distribution, % bases >= N x
#   - samtools stats: alignment summary stats
#   - bcftools stats: per-sample variant counts, Ti/Tv, het/hom
#   - sex_check.py: chrX / chrY coverage ratio -> inferred sex
# Cohort:
#   - GATK CollectVariantCallingMetrics (vs dbSNP truth set)
#   - cohort_qc.py: aggregate everything into one TSV + a JSON manifest
#   - VerifyBamID stub: contamination check (see note below)
#
# A note on contamination: a real run should use VerifyBamID2 with a
# population-allele-frequency resource (e.g. 1000G.phase3.100k.b38.vcf.gz.dat.*).
# That resource is large and licensed differently per region, so this pipeline
# leaves it as a clearly-labelled TODO rather than silently shipping a half-
# working check. See workflow/scripts/cohort_qc.py for where the value is read.
# =============================================================================

rule mosdepth:
    """Per-sample coverage at thresholds defined by qc.coverage_threshold."""
    input:
        bam = BQSR / "{sample}.bqsr.bam",
        bai = BQSR / "{sample}.bqsr.bam.bai",
    output:
        summary  = QC / "mosdepth" / "{sample}.mosdepth.summary.txt",
        regions  = QC / "mosdepth" / "{sample}.regions.bed.gz",
        thresh   = QC / "mosdepth" / "{sample}.thresholds.bed.gz",
    log:
        "logs/qc/mosdepth_{sample}.log",
    params:
        prefix    = lambda wc: str(QC / "mosdepth" / wc.sample),
        threshold = config["qc"]["coverage_threshold"],
    threads: 4
    resources:
        mem_mb  = 4000,
        runtime = 60,
    conda:
        "../envs/qc.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        mosdepth -t {threads} --no-per-base --by 1000 \
            --thresholds 1,10,{params.threshold},30 \
            {params.prefix} {input.bam} 2> {log}
        """


rule samtools_stats:
    input:
        bam = BQSR / "{sample}.bqsr.bam",
    output:
        stats = QC / "samtools" / "{sample}.stats.txt",
    log:
        "logs/qc/samtools_stats_{sample}.log",
    conda:
        "../envs/align.yaml"
    resources:
        mem_mb = 2000,
        runtime = 30,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.stats})
        samtools stats {input.bam} > {output.stats} 2> {log}
        """


rule bcftools_stats_per_sample:
    input:
        vcf = JOINT / "cohort.filtered.vcf.gz",
    output:
        stats = QC / "bcftools" / "{sample}.stats.txt",
    log:
        "logs/qc/bcftools_stats_{sample}.log",
    conda:
        "../envs/align.yaml"
    resources:
        mem_mb = 4000,
        runtime = 30,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.stats})
        bcftools stats -s {wildcards.sample} -f PASS {input.vcf} > {output.stats} 2> {log}
        """


rule sex_check:
    """Infer genetic sex from chrX vs chrY coverage and compare to the sample sheet."""
    input:
        summary = QC / "mosdepth" / "{sample}.mosdepth.summary.txt",
    output:
        json = QC / "sex" / "{sample}.sex.json",
    log:
        "logs/qc/sex_check_{sample}.log",
    params:
        reported_sex = lambda wc: get_sex(wc.sample),
    conda:
        "../envs/qc.yaml"
    resources:
        mem_mb = 1000,
        runtime = 10,
    script:
        "../scripts/sex_check.py"


rule collect_variant_calling_metrics:
    input:
        vcf   = JOINT / "cohort.filtered.vcf.gz",
        dbsnp = config["reference"]["known_sites"][-1],   # last known-sites is dbSNP
    output:
        summary = QC / "cohort" / "cohort.variant_calling_summary_metrics",
        detail  = QC / "cohort" / "cohort.variant_calling_detail_metrics",
    log:
        "logs/qc/collect_variant_calling_metrics.log",
    params:
        prefix = lambda wc, output: str(QC / "cohort" / "cohort"),
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb = 8000,
        runtime = 60,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {params.prefix})
        gatk --java-options "-Xmx{resources.mem_mb}m" CollectVariantCallingMetrics \
            -I {input.vcf} \
            --DBSNP {input.dbsnp} \
            -O {params.prefix} \
            2> {log}
        """


rule cohort_qc_aggregate:
    """Merge per-sample QC artefacts into a single tidy TSV + a JSON manifest."""
    input:
        mosdepth   = expand(str(QC / "mosdepth" / "{sample}.mosdepth.summary.txt"), sample=SAMPLES),
        samtools   = expand(str(QC / "samtools" / "{sample}.stats.txt"), sample=SAMPLES),
        bcftools   = expand(str(QC / "bcftools" / "{sample}.stats.txt"), sample=SAMPLES),
        sex        = expand(str(QC / "sex" / "{sample}.sex.json"), sample=SAMPLES),
        cohort_sum = QC / "cohort" / "cohort.variant_calling_summary_metrics",
    output:
        tsv  = QC / "cohort_qc.tsv",
        json = QC / "cohort_qc.json",
    log:
        "logs/qc/cohort_qc_aggregate.log",
    params:
        samples            = SAMPLES,
        coverage_threshold = config["qc"]["coverage_threshold"],
        titv_min           = config["qc"]["titv_min"],
        titv_max           = config["qc"]["titv_max"],
        het_hom_min        = config["qc"]["het_hom_min"],
        het_hom_max        = config["qc"]["het_hom_max"],
    conda:
        "../envs/qc.yaml"
    resources:
        mem_mb = 4000,
        runtime = 30,
    script:
        "../scripts/cohort_qc.py"
