# =============================================================================
# bqsr.smk — Base Quality Score Recalibration.
#
# Sequencers systematically over- or under-estimate base call qualities at
# certain cycles, contexts, and dinucleotides. BQSR builds a per-sample model
# from sites that should be invariant (known sites: dbSNP + Mills indels) and
# then rewrites the per-base qualities so downstream callers see calibrated
# evidence.
#
# Inputs : results/dedup/{sample}.dedup.bam, known-sites VCFs, reference FASTA
# Outputs: results/bqsr/{sample}.bqsr.bam (+ .bai), recal table
# Helpers: known_sites_args() — defined in common.smk
# =============================================================================


rule base_recalibrator:
    """First pass: learn the recalibration model from non-variant sites."""
    input:
        bam   = DEDUP / "{sample}.dedup.bam",
        fasta = config["reference"]["fasta"],
        ks    = config["reference"]["known_sites"],
        intervals = config["reference"]["calling_intervals"],
    output:
        table = BQSR / "{sample}.recal.table",
    log:
        "logs/bqsr/{sample}.recal.log",
    params:
        known_sites = known_sites_args,
    threads: config["resources"]["bqsr"]["threads"]
    resources:
        mem_mb  = config["resources"]["bqsr"]["mem_mb"],
        runtime = config["resources"]["bqsr"]["runtime"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        gatk --java-options "-Xmx{resources.mem_mb}m" BaseRecalibrator \
            -I {input.bam} \
            -R {input.fasta} \
            -L {input.intervals} \
            {params.known_sites} \
            -O {output.table} \
            2> {log}
        """


rule apply_bqsr:
    """Second pass: rewrite the BAM with the recalibrated quality scores."""
    input:
        bam   = DEDUP / "{sample}.dedup.bam",
        fasta = config["reference"]["fasta"],
        table = BQSR / "{sample}.recal.table",
        intervals = config["reference"]["calling_intervals"],
    output:
        bam = BQSR / "{sample}.bqsr.bam",
        bai = BQSR / "{sample}.bqsr.bam.bai",
    log:
        "logs/bqsr/{sample}.apply.log",
    threads: config["resources"]["bqsr"]["threads"]
    resources:
        mem_mb  = config["resources"]["bqsr"]["mem_mb"],
        runtime = config["resources"]["bqsr"]["runtime"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        gatk --java-options "-Xmx{resources.mem_mb}m" ApplyBQSR \
            -I {input.bam} \
            -R {input.fasta} \
            -L {input.intervals} \
            --bqsr-recal-file {input.table} \
            -O {output.bam} \
            2> {log}
        # GATK ApplyBQSR writes foo.bai next to foo.bam. Downstream rules want
        # the *.bam.bai form; create it as a copy if it's missing.
        bam_path="{output.bam}"
        gatk_bai="${{bam_path%.bam}}.bai"
        if [ -f "$gatk_bai" ] && [ ! -f "{output.bai}" ]; then
            cp "$gatk_bai" "{output.bai}"
        fi
        """
