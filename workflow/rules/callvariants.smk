# =============================================================================
# callvariants.smk — Per-sample HaplotypeCaller in GVCF mode.
#
# HaplotypeCaller in -ERC GVCF mode emits records at every site (variant and
# non-variant) so that joint genotyping later can distinguish "reference call"
# from "no data". The output gVCFs are the input to GenomicsDBImport.
#
# Inputs : results/bqsr/{sample}.bqsr.bam, reference FASTA, calling intervals
# Outputs: results/gvcf/{sample}.g.vcf.gz (+ .tbi)
# =============================================================================

rule haplotypecaller:
    input:
        bam       = BQSR / "{sample}.bqsr.bam",
        fasta     = config["reference"]["fasta"],
        intervals = config["reference"]["calling_intervals"],
    output:
        gvcf = GVCF / "{sample}.g.vcf.gz",
        tbi  = GVCF / "{sample}.g.vcf.gz.tbi",
    log:
        "logs/haplotypecaller/{sample}.log",
    threads: config["resources"]["haplotypecaller"]["threads"]
    resources:
        mem_mb  = config["resources"]["haplotypecaller"]["mem_mb"],
        runtime = config["resources"]["haplotypecaller"]["runtime"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        gatk --java-options "-Xmx{resources.mem_mb}m" HaplotypeCaller \
            -R {input.fasta} \
            -I {input.bam} \
            -L {input.intervals} \
            -O {output.gvcf} \
            -ERC GVCF \
            --native-pair-hmm-threads {threads} \
            2> {log}
        """
