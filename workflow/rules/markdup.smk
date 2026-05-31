# =============================================================================
# markdup.smk — flag PCR / optical duplicates.
#
# Two interchangeable implementations gated on config["dedup"]["tool"]:
#   - "samtools"   : samtools markdup, single-node, fast for moderate cohorts.
#   - "gatk_spark" : GATK MarkDuplicatesSpark, parallelisable on multi-core
#                    workers, slightly different metrics but equivalent flags.
#
# Inputs : results/align/{sample}.sorted.bam
# Outputs: results/dedup/{sample}.dedup.bam (+ .bai), per-sample metrics file
# =============================================================================

DEDUP_TOOL = config.get("dedup", {}).get("tool", "samtools").lower()
if DEDUP_TOOL not in {"samtools", "gatk_spark"}:
    raise ValueError(
        f"config.dedup.tool must be 'samtools' or 'gatk_spark', got {DEDUP_TOOL!r}"
    )


if DEDUP_TOOL == "samtools":

    rule markdup_samtools:
        """
        samtools markdup pipeline: sort by name -> fixmate -> sort by coord -> markdup.
        Flag (not remove) duplicates so downstream filtering can still see them.
        """
        input:
            bam = ALIGN / "{sample}.sorted.bam",
        output:
            bam     = DEDUP / "{sample}.dedup.bam",
            bai     = DEDUP / "{sample}.dedup.bam.bai",
            metrics = DEDUP / "{sample}.markdup_metrics.txt",
        log:
            "logs/markdup/{sample}.log",
        threads: config["resources"]["markdup"]["threads"]
        resources:
            mem_mb  = config["resources"]["markdup"]["mem_mb"],
            runtime = config["resources"]["markdup"]["runtime"],
        conda:
            "../envs/align.yaml"
        shell:
            r"""
            set -euo pipefail
            tmp=$(mktemp -d)
            trap "rm -rf $tmp" EXIT
            samtools sort -@ {threads} -n -T $tmp/qn -o $tmp/qn.bam {input.bam} 2> {log}
            samtools fixmate -@ {threads} -m $tmp/qn.bam $tmp/fix.bam 2>> {log}
            samtools sort -@ {threads} -T $tmp/co -o $tmp/co.bam $tmp/fix.bam 2>> {log}
            samtools markdup -@ {threads} -s -f {output.metrics} $tmp/co.bam {output.bam} 2>> {log}
            samtools index -@ {threads} {output.bam} 2>> {log}
            """

else:  # gatk_spark

    rule markdup_spark:
        """
        GATK MarkDuplicatesSpark: combines sort + markdup in one Spark job.
        Writes its own .bai. Spark runs locally with --local-cores by default.
        """
        input:
            bam   = ALIGN / "{sample}.sorted.bam",
            fasta = config["reference"]["fasta"],
        output:
            bam     = DEDUP / "{sample}.dedup.bam",
            bai     = DEDUP / "{sample}.dedup.bam.bai",
            metrics = DEDUP / "{sample}.markdup_metrics.txt",
        log:
            "logs/markdup/{sample}.log",
        threads: config["resources"]["markdup"]["threads"]
        resources:
            mem_mb  = config["resources"]["markdup"]["mem_mb"],
            runtime = config["resources"]["markdup"]["runtime"],
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            set -euo pipefail
            gatk --java-options "-Xmx{resources.mem_mb}m" MarkDuplicatesSpark \
                -I {input.bam} \
                -O {output.bam} \
                -M {output.metrics} \
                --spark-master local[{threads}] \
                --conf 'spark.local.dir=$TMPDIR' \
                2> {log}
            """
