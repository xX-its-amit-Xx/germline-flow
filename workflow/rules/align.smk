# =============================================================================
# align.smk — BWA-MEM2 alignment of paired FASTQs to a sorted, indexed BAM.
#
# Inputs : per-sample R1/R2 FASTQs (paths from samples.tsv)
# Outputs: results/align/{sample}.sorted.bam (+ .bai)
# Notes  : The @RG read group is set per sample because GATK requires it for
#          every downstream step. SM and ID both default to the sample name.
# =============================================================================

rule index_reference:
    """
    One-time BWA-MEM2 + samtools + GATK index build for the reference FASTA.
    Skipped after the first successful run because all outputs are present.
    """
    input:
        fasta = config["reference"]["fasta"],
    output:
        bwa_idx = config["reference"]["fasta"] + ".bwt.2bit.64",
        fai     = config["reference"]["fasta"] + ".fai",
        dict    = config["reference"]["fasta"].rsplit(".", 1)[0] + ".dict",
    log:
        "logs/index_reference.log",
    threads: 4
    resources:
        mem_mb  = 16000,
        runtime = 120,
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        (
          bwa-mem2 index {input.fasta}
          samtools faidx {input.fasta}
          gatk CreateSequenceDictionary -R {input.fasta} -O {output.dict}
        ) &> {log}
        """


rule bwa_mem2_align:
    """
    Align paired FASTQs with bwa-mem2 and stream into samtools sort.
    Produces a coordinate-sorted, indexed BAM tagged with the sample read group.
    """
    input:
        unpack(get_fastqs),
        fasta   = config["reference"]["fasta"],
        bwa_idx = config["reference"]["fasta"] + ".bwt.2bit.64",
    output:
        bam = ALIGN / "{sample}.sorted.bam",
        bai = ALIGN / "{sample}.sorted.bam.bai",
    log:
        "logs/align/{sample}.log",
    params:
        rg = lambda wc: (
            r"@RG\tID:{s}.L1\tSM:{s}\tLB:{s}.{lib}\tPL:{pl}".format(
                s   = wc.sample,
                lib = config["alignment"]["read_group_library_suffix"],
                pl  = config["alignment"]["read_group_platform"],
            )
        ),
    threads: config["resources"]["align"]["threads"]
    resources:
        mem_mb  = config["resources"]["align"]["mem_mb"],
        runtime = config["resources"]["align"]["runtime"],
    conda:
        "../envs/align.yaml"
    shell:
        r"""
        set -euo pipefail
        bwa-mem2 mem \
            -t {threads} \
            -R '{params.rg}' \
            -K 100000000 -Y \
            {input.fasta} {input.r1} {input.r2} 2> {log} \
        | samtools sort -@ {threads} -m 1G -O bam -o {output.bam} - 2>> {log}
        samtools index -@ {threads} {output.bam} 2>> {log}
        """
