# Recipe 1 — Trio analysis (NA12878 / NA12891 / NA12892)

The CEPH/Utah pedigree 1463 is the most-sequenced human family on the planet
and its members are the de-facto truth set for human germline variant callers
(the Genome in a Bottle consortium, "GIAB"). This recipe runs germline-flow on
the trio's chr20 slice and benchmarks the resulting VCF against the GIAB v4.2.1
high-confidence truth set with `rtg vcfeval`.

## What you'll have at the end

* A jointly-called, hard-filtered VCF for the three-sample trio on GRCh38 chr20
  (~1.5 MB).
* A PLINK2 fileset with the trio's chr20 genotypes.
* A cohort QC HTML report — PCA should put the three on top of each other
  (they're EUR), and KING should report two 1st-degree pairs
  (mother-child, father-child) and one ~0 pair (parents are unrelated).
* An `rtg vcfeval` summary giving precision / recall vs the GIAB truth on
  the high-confidence regions of chr20.

## Data acquisition

Pin a working directory and pull everything into it. Total: ~6 GB.

```bash
mkdir -p data/trio/{fastq,reference,truth} && cd data/trio

# ---- GRCh38 reference + GATK known-sites bundle -----------------------------
# Use the Broad's no-alt analysis-set assembly. Mirrors:
#   gs://gcp-public-data--broad-references/hg38/v0/
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/

curl -sSL -o reference/GRCh38.primary_assembly.genome.fa.gz \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip reference/GRCh38.primary_assembly.genome.fa.gz
mv reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference/GRCh38.primary_assembly.genome.fa 2>/dev/null || true

# Broad bundle — Mills indels + dbSNP + Hapmap/Omni/1000G for VQSR.
for f in \
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
  Homo_sapiens_assembly38.dbsnp138.vcf \
  Homo_sapiens_assembly38.dbsnp138.vcf.idx \
  wgs_calling_regions.hg38.interval_list \
  hapmap_3.3.hg38.vcf.gz hapmap_3.3.hg38.vcf.gz.tbi \
  1000G_omni2.5.hg38.vcf.gz 1000G_omni2.5.hg38.vcf.gz.tbi \
  1000G_phase1.snps.high_confidence.hg38.vcf.gz 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
do
  curl -sSL -o "reference/${f}" "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/${f}"
done

# ---- NA12878 / NA12891 / NA12892 FASTQs from the GIAB 30x downsampled set --
# (Full 30x BAMs would be ~80 GB each; for a single-chromosome benchmark
# you'll want to subset to chr20 first — see "Pre-subset to chr20" below.)
# Direct FASTQ from the 1000 Genomes Project Illumina runs:
GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data"
curl -sSL -o fastq/NA12878_R1.fastq.gz "${GIAB_BASE}/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam"
# (Use samtools fastq to convert BAM->FASTQ; see Broad's "downsampled-30x" page.)

# ---- GIAB v4.2.1 high-confidence truth for NA12878 -------------------------
TRUTH_BASE="https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38"
curl -sSL -o truth/HG001_GRCh38.high_confidence.vcf.gz \
  "${TRUTH_BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
curl -sSL -o truth/HG001_GRCh38.high_confidence.vcf.gz.tbi \
  "${TRUTH_BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi"
curl -sSL -o truth/HG001_GRCh38.high_confidence.bed \
  "${TRUTH_BASE}/HG001_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
```

### Pre-subset to chr20 (optional but recommended for laptop runs)

If you'd rather run on just chr20 (~1/16 of the genome, ~30 min on 8 cores
instead of overnight):

```bash
# Subset the truth set and BED.
bcftools view -r chr20 -Oz -o truth/HG001_chr20.vcf.gz truth/HG001_GRCh38.high_confidence.vcf.gz
tabix -p vcf truth/HG001_chr20.vcf.gz
awk '$1 == "chr20"' truth/HG001_GRCh38.high_confidence.bed > truth/HG001_chr20.bed

# Subset the calling intervals.
grep -E '^@|^chr20\b' reference/wgs_calling_regions.hg38.interval_list \
  > reference/wgs_calling_regions.chr20.interval_list

# Subset the FASTQs by re-aligning a chr20-only BAM with samtools view + fastq.
# (Skipped here; see https://github.com/Illumina/Pisces for an end-to-end script.)
```

## Config diff

Edit `config/samples.tsv`:

```tsv
sample    fastq_1                                fastq_2                                sex
NA12878   data/trio/fastq/NA12878_R1.fastq.gz    data/trio/fastq/NA12878_R2.fastq.gz    F
NA12891   data/trio/fastq/NA12891_R1.fastq.gz    data/trio/fastq/NA12891_R2.fastq.gz    M
NA12892   data/trio/fastq/NA12892_R1.fastq.gz    data/trio/fastq/NA12892_R2.fastq.gz    F
```

Edit `config/config.yaml`:

```yaml
reference:
  fasta: data/trio/reference/GRCh38.primary_assembly.genome.fa
  known_sites:
    - data/trio/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    - data/trio/reference/Homo_sapiens_assembly38.dbsnp138.vcf
  calling_intervals: data/trio/reference/wgs_calling_regions.chr20.interval_list

joint_calling:
  filter_mode: hard_filter        # 3 samples is far below the VQSR cohort minimum

plink_export:
  call_rate_min: 0.90             # trio is small; loosen slightly
  maf_min: 0.0                    # keep singletons for trio analysis
```

## Run

```bash
snakemake --use-conda --cores 8
```

Expected wall-clock on chr20 only: ~30-45 min on 8 cores with mamba env caches warm.

## Benchmark against GIAB truth (NA12878 only)

```bash
# rtg vcfeval lives in the GATK env via the rtg-tools bioconda pkg.
mamba install -n base -c bioconda rtg-tools=3.12.1

# Build the SDF (sequence-dictionary format) for the reference once.
rtg format -o data/trio/reference/GRCh38.sdf data/trio/reference/GRCh38.primary_assembly.genome.fa

# Extract NA12878 from the joint VCF.
bcftools view -s NA12878 -Oz -o results/joint/NA12878.chr20.vcf.gz \
  results/joint/cohort.filtered.vcf.gz
tabix -p vcf results/joint/NA12878.chr20.vcf.gz

# Run the comparison, restricted to the high-confidence chr20 regions.
rtg vcfeval \
  --baseline=data/trio/truth/HG001_chr20.vcf.gz \
  --bed-regions=data/trio/truth/HG001_chr20.bed \
  --calls=results/joint/NA12878.chr20.vcf.gz \
  --template=data/trio/reference/GRCh38.sdf \
  --output=results/benchmark/na12878_chr20 \
  --sample=NA12878
cat results/benchmark/na12878_chr20/summary.txt
```

For chr20 only with the hard-filter branch you should see precision &gt; 0.99
and recall in the 0.97-0.99 range. Numbers below that usually mean coverage
is uneven (check the report's per-sample coverage figure) or the calling
intervals didn't get restricted correctly.

## Interpret the report

Open `results/report/cohort_qc.html`. You should see:

* **Coverage** &mdash; ~30x for all three samples; all bars green.
* **Ti/Tv** &mdash; all three near 2.05.
* **PCA** &mdash; three points clustered tightly; if a point is far from the
  others, your FASTQs are mis-labelled or contaminated.
* **Relatedness** &mdash; two pairs in `1st_degree` (NA12878-NA12891,
  NA12878-NA12892) and one pair `unrelated` (NA12891-NA12892). If you see a
  3rd 1st-degree pair, the samples are mis-labelled.

## Going further

* For Mendelian-error analysis, feed the joint VCF to
  [`rtg mendelian`](https://github.com/RealTimeGenomics/rtg-tools) with a
  PED file. Mendelian-error rate should be &lt; 1% on PASS calls.
* For a stricter benchmark, restrict to the GIAB "Tier1" stratifications
  (`benchmark/genome-stratifications`).
