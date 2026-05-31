# =============================================================================
# jointgeno.smk — GenomicsDBImport + GenotypeGVCFs + VQSR / hard-filter.
#
# Joint genotyping uses information across all samples to call low-confidence
# variants that a per-sample caller would have missed, and to assign the correct
# genotype (0/0 vs ./.) at sites where one sample has data and another doesn't.
#
# Filtering branches on config["joint_calling"]["filter_mode"]:
#   - "vqsr"        : statistical filtering trained on truth/training sets
#                     (requires ~30+ WGS or ~100+ WES samples)
#   - "hard_filter" : the GATK-recommended hard-filter expressions, suitable
#                     for small cohorts, single-sample, or callsets where VQSR
#                     fails to converge.
# =============================================================================

FILTER_MODE = config["joint_calling"]["filter_mode"].lower()
if FILTER_MODE not in {"vqsr", "hard_filter"}:
    raise ValueError(
        f"joint_calling.filter_mode must be 'vqsr' or 'hard_filter', got {FILTER_MODE!r}"
    )


rule write_sample_map:
    """
    Generate the tab-separated sample-name-map that GenomicsDBImport reads.
    Kept as its own rule (rather than inlined into genomicsdb_import) so the
    file ends up under JOINT/ where it's visible for debugging.
    """
    input:
        gvcfs = expand(str(GVCF / "{sample}.g.vcf.gz"), sample=SAMPLES),
    output:
        sample_map = JOINT / "sample_map.tsv",
    log:
        "logs/jointgeno/write_sample_map.log",
    params:
        # printf-format string with literal \t and \n escapes; printf interprets
        # them at runtime. See workflow/rules/common.smk for the builder.
        fmt = sample_map_format_string(),
    conda:
        # The rule only needs coreutils, but snakemake --lint requires every
        # rule to pin an environment. Reuse the lightweight align env.
        "../envs/align.yaml"
    resources:
        mem_mb  = 500,
        runtime = 5,
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.sample_map})
        printf '{params.fmt}\n' > {output.sample_map} 2> {log}
        """


rule genomicsdb_import:
    """
    Build a GenomicsDB workspace from all per-sample gVCFs.
    The workspace replaces the older CombineGVCFs step and scales to
    thousands of samples without n^2 memory blowup.
    """
    input:
        gvcfs     = expand(str(GVCF / "{sample}.g.vcf.gz"), sample=SAMPLES),
        tbis      = expand(str(GVCF / "{sample}.g.vcf.gz.tbi"), sample=SAMPLES),
        intervals = config["reference"]["calling_intervals"],
        sample_map = JOINT / "sample_map.tsv",
    output:
        workspace = directory(JOINT / "genomicsdb"),
    log:
        "logs/jointgeno/genomicsdb_import.log",
    threads: config["resources"]["jointgeno"]["threads"]
    resources:
        mem_mb  = config["resources"]["jointgeno"]["mem_mb"],
        runtime = config["resources"]["jointgeno"]["runtime"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        # GenomicsDBImport refuses to overwrite an existing workspace dir.
        rm -rf {output.workspace}
        gatk --java-options "-Xmx{resources.mem_mb}m" GenomicsDBImport \
            --sample-name-map {input.sample_map} \
            --genomicsdb-workspace-path {output.workspace} \
            -L {input.intervals} \
            --reader-threads {threads} \
            --batch-size 50 \
            2> {log}
        """


rule genotype_gvcfs:
    """Joint-genotype all samples from the GenomicsDB workspace."""
    input:
        workspace = JOINT / "genomicsdb",
        fasta     = config["reference"]["fasta"],
        intervals = config["reference"]["calling_intervals"],
    output:
        vcf = JOINT / "cohort.raw.vcf.gz",
        tbi = JOINT / "cohort.raw.vcf.gz.tbi",
    log:
        "logs/jointgeno/genotype_gvcfs.log",
    threads: config["resources"]["jointgeno"]["threads"]
    resources:
        mem_mb  = config["resources"]["jointgeno"]["mem_mb"],
        runtime = config["resources"]["jointgeno"]["runtime"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        set -euo pipefail
        gatk --java-options "-Xmx{resources.mem_mb}m" GenotypeGVCFs \
            -R {input.fasta} \
            -V gendb://{input.workspace} \
            -L {input.intervals} \
            -O {output.vcf} \
            2> {log}
        """


# ---- filtering branch --------------------------------------------------------

if FILTER_MODE == "hard_filter":

    rule select_and_hardfilter_snps:
        input:
            vcf   = JOINT / "cohort.raw.vcf.gz",
            fasta = config["reference"]["fasta"],
        output:
            vcf = JOINT / "cohort.snps.filtered.vcf.gz",
            tbi = JOINT / "cohort.snps.filtered.vcf.gz.tbi",
        log:
            "logs/jointgeno/hardfilter_snps.log",
        params:
            filters = hard_filter_args("snp"),
        conda:
            "../envs/gatk.yaml"
        resources:
            mem_mb  = 8000,
            runtime = 60,
        shell:
            r"""
            set -euo pipefail
            tmp=$(mktemp --suffix=.vcf.gz)
            gatk --java-options "-Xmx{resources.mem_mb}m" SelectVariants \
                -R {input.fasta} -V {input.vcf} --select-type-to-include SNP -O $tmp 2> {log}
            gatk --java-options "-Xmx{resources.mem_mb}m" VariantFiltration \
                -R {input.fasta} -V $tmp {params.filters} -O {output.vcf} 2>> {log}
            rm -f $tmp $tmp.tbi
            """

    rule select_and_hardfilter_indels:
        input:
            vcf   = JOINT / "cohort.raw.vcf.gz",
            fasta = config["reference"]["fasta"],
        output:
            vcf = JOINT / "cohort.indels.filtered.vcf.gz",
            tbi = JOINT / "cohort.indels.filtered.vcf.gz.tbi",
        log:
            "logs/jointgeno/hardfilter_indels.log",
        params:
            filters = hard_filter_args("indel"),
        conda:
            "../envs/gatk.yaml"
        resources:
            mem_mb  = 8000,
            runtime = 60,
        shell:
            r"""
            set -euo pipefail
            tmp=$(mktemp --suffix=.vcf.gz)
            gatk --java-options "-Xmx{resources.mem_mb}m" SelectVariants \
                -R {input.fasta} -V {input.vcf} \
                --select-type-to-include INDEL --select-type-to-include MIXED \
                -O $tmp 2> {log}
            gatk --java-options "-Xmx{resources.mem_mb}m" VariantFiltration \
                -R {input.fasta} -V $tmp {params.filters} -O {output.vcf} 2>> {log}
            rm -f $tmp $tmp.tbi
            """

    rule merge_filtered:
        """Concat the SNP + INDEL filtered callsets back into a single VCF."""
        input:
            snps   = JOINT / "cohort.snps.filtered.vcf.gz",
            indels = JOINT / "cohort.indels.filtered.vcf.gz",
        output:
            vcf = JOINT / "cohort.filtered.vcf.gz",
            tbi = JOINT / "cohort.filtered.vcf.gz.tbi",
        log:
            "logs/jointgeno/merge_filtered.log",
        conda:
            "../envs/align.yaml"   # bcftools lives in the align env
        resources:
            mem_mb  = 4000,
            runtime = 30,
        shell:
            r"""
            set -euo pipefail
            bcftools concat -a -Oz -o {output.vcf} {input.snps} {input.indels} 2> {log}
            bcftools index -t {output.vcf} 2>> {log}
            """

else:  # vqsr

    rule vqsr_snp:
        input:
            vcf       = JOINT / "cohort.raw.vcf.gz",
            fasta     = config["reference"]["fasta"],
            hapmap    = config["reference"]["vqsr_resources"]["snp"]["hapmap"],
            omni      = config["reference"]["vqsr_resources"]["snp"]["omni"],
            g1k       = config["reference"]["vqsr_resources"]["snp"]["g1k"],
            dbsnp     = config["reference"]["vqsr_resources"]["snp"]["dbsnp"],
        output:
            recal  = JOINT / "snp.recal",
            tranches = JOINT / "snp.tranches",
        log:
            "logs/jointgeno/vqsr_snp.log",
        threads: config["resources"]["vqsr"]["threads"]
        resources:
            mem_mb  = config["resources"]["vqsr"]["mem_mb"],
            runtime = config["resources"]["vqsr"]["runtime"],
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            set -euo pipefail
            gatk --java-options "-Xmx{resources.mem_mb}m" VariantRecalibrator \
                -R {input.fasta} -V {input.vcf} \
                --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {input.hapmap} \
                --resource:omni,known=false,training=true,truth=true,prior=12.0 {input.omni} \
                --resource:1000G,known=false,training=true,truth=false,prior=10.0 {input.g1k} \
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
                -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
                -mode SNP \
                -O {output.recal} --tranches-file {output.tranches} \
                2> {log}
            """

    rule vqsr_indel:
        input:
            vcf   = JOINT / "cohort.raw.vcf.gz",
            fasta = config["reference"]["fasta"],
            mills = config["reference"]["vqsr_resources"]["indel"]["mills"],
            dbsnp = config["reference"]["vqsr_resources"]["indel"]["dbsnp"],
        output:
            recal    = JOINT / "indel.recal",
            tranches = JOINT / "indel.tranches",
        log:
            "logs/jointgeno/vqsr_indel.log",
        threads: config["resources"]["vqsr"]["threads"]
        resources:
            mem_mb  = config["resources"]["vqsr"]["mem_mb"],
            runtime = config["resources"]["vqsr"]["runtime"],
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            set -euo pipefail
            gatk --java-options "-Xmx{resources.mem_mb}m" VariantRecalibrator \
                -R {input.fasta} -V {input.vcf} \
                --resource:mills,known=false,training=true,truth=true,prior=12.0 {input.mills} \
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {input.dbsnp} \
                -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
                -mode INDEL --max-gaussians 4 \
                -O {output.recal} --tranches-file {output.tranches} \
                2> {log}
            """

    rule apply_vqsr:
        input:
            vcf            = JOINT / "cohort.raw.vcf.gz",
            fasta          = config["reference"]["fasta"],
            snp_recal      = JOINT / "snp.recal",
            snp_tranches   = JOINT / "snp.tranches",
            indel_recal    = JOINT / "indel.recal",
            indel_tranches = JOINT / "indel.tranches",
        output:
            vcf = JOINT / "cohort.filtered.vcf.gz",
            tbi = JOINT / "cohort.filtered.vcf.gz.tbi",
        log:
            "logs/jointgeno/apply_vqsr.log",
        params:
            snp_ts   = config["joint_calling"]["vqsr"]["snp_truth_sensitivity"],
            indel_ts = config["joint_calling"]["vqsr"]["indel_truth_sensitivity"],
        threads: config["resources"]["vqsr"]["threads"]
        resources:
            mem_mb  = config["resources"]["vqsr"]["mem_mb"],
            runtime = config["resources"]["vqsr"]["runtime"],
        conda:
            "../envs/gatk.yaml"
        shell:
            r"""
            set -euo pipefail
            tmp=$(mktemp --suffix=.vcf.gz)
            gatk --java-options "-Xmx{resources.mem_mb}m" ApplyVQSR \
                -R {input.fasta} -V {input.vcf} \
                --recal-file {input.indel_recal} \
                --tranches-file {input.indel_tranches} \
                --truth-sensitivity-filter-level {params.indel_ts} \
                --create-output-variant-index true \
                -mode INDEL -O $tmp 2> {log}
            gatk --java-options "-Xmx{resources.mem_mb}m" ApplyVQSR \
                -R {input.fasta} -V $tmp \
                --recal-file {input.snp_recal} \
                --tranches-file {input.snp_tranches} \
                --truth-sensitivity-filter-level {params.snp_ts} \
                --create-output-variant-index true \
                -mode SNP -O {output.vcf} 2>> {log}
            rm -f $tmp $tmp.tbi
            """
