# Recipe 2 — Small-cohort GWAS pilot (1000 Genomes chr22)

A toy end-to-end GWAS using 50 unrelated EUR samples from 1000 Genomes Phase 3
on chr22, simulating a quantitative phenotype with a single causal SNP. The
point isn't the science (the phenotype is made up); the point is to show how to
take germline-flow's PLINK2 output into `plink2 --glm` for an association sketch
and how to QC the result.

## What you'll have at the end

* Joint VCF for 50 samples on chr22 (~80 MB).
* PLINK2 pgen/pvar/psam (after sample/variant QC).
* Cohort QC report — PCA should separate EUR sub-populations (CEU, GBR, FIN…).
* A simulated phenotype TSV.
* PLINK2 GWAS summary statistics + a Manhattan plot.

## Data acquisition

```bash
mkdir -p data/g1k/{fastq,reference,pheno} && cd data/g1k

# Pick 50 unrelated EUR samples from the 1000G phase 3 pedigree.
# (See https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel)
cat > samples.tsv <<'EOF'
HG00096
HG00097
HG00099
HG00100
HG00101
HG00102
HG00103
HG00105
HG00106
HG00107
HG00108
HG00109
HG00110
HG00111
HG00112
HG00113
HG00114
HG00115
HG00116
HG00117
HG00118
HG00119
HG00120
HG00121
HG00122
HG00123
HG00125
HG00126
HG00127
HG00128
HG00129
HG00130
HG00131
HG00132
HG00133
HG00136
HG00137
HG00138
HG00139
HG00140
HG00141
HG00142
HG00143
HG00145
HG00146
HG00148
HG00149
HG00150
HG00151
HG00154
EOF

# Each sample's FASTQs live under
# https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/data/<POP>/<SAMPLE>/
# Subset to chr22 with samtools view + samtools fastq, OR start from the
# 1000G phase 3 high-coverage CRAMs and convert.
```

Reference and known sites: same Broad bundle as
[recipe 1](01-trio-analysis.md), restrict the calling intervals to chr22 only.

## Config diff

```yaml
# config/config.yaml
reference:
  fasta: data/g1k/reference/GRCh38.primary_assembly.genome.fa
  known_sites:
    - data/g1k/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    - data/g1k/reference/Homo_sapiens_assembly38.dbsnp138.vcf
  calling_intervals: data/g1k/reference/wgs_calling_regions.chr22.interval_list

joint_calling:
  filter_mode: hard_filter        # 50 samples is still well below the VQSR threshold

plink_export:
  call_rate_min: 0.97             # population study: stricter
  maf_min: 0.05                   # GWAS standard
  hwe_p_min: 1.0e-6
  king_cutoff: 0.0884
  pca_components: 10
```

## Run the pipeline

```bash
snakemake --use-conda --cores 16
```

Expected wall-clock: ~3-5 h on 16 cores; faster on a SLURM cluster.

## Simulate a phenotype

Pick a single chr22 SNP as the causal variant. Make height = 170 + 2.5 * dosage + N(0, 5).

```bash
mkdir -p data/g1k/pheno

# Pick a common chr22 variant with reasonable MAF in the cohort.
plink2 \
  --pfile results/plink/cohort.qc \
  --chr 22 --from-bp 30000000 --to-bp 31000000 \
  --maf 0.20 \
  --write-snplist \
  --out data/g1k/pheno/candidate_snps
CAUSAL=$(head -n 1 data/g1k/pheno/candidate_snps.snplist)
echo "causal SNP: ${CAUSAL}"

# Extract dosage for the causal SNP.
plink2 \
  --pfile results/plink/cohort.qc \
  --snp "${CAUSAL}" \
  --export A \
  --out data/g1k/pheno/causal_dosage

# Simulate phenotype: height = 170 + 2.5 * dosage + N(0, 5).
python3 - data/g1k/pheno/causal_dosage.raw data/g1k/pheno/height.pheno <<'PY'
import csv, random, sys
random.seed(42)
infile, outfile = sys.argv[1], sys.argv[2]
with open(infile) as fh:
    reader = csv.reader(fh, delimiter="\t")
    header = next(reader)
    dose_col = -1  # last column is the dosage of the chosen SNP
    rows = list(reader)
with open(outfile, "w") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["#FID", "IID", "height"])
    for r in rows:
        fid, iid = r[0], r[1]
        dose = float(r[dose_col]) if r[dose_col] != "NA" else 1.0
        h = 170 + 2.5 * dose + random.gauss(0, 5)
        w.writerow([fid, iid, f"{h:.3f}"])
PY
```

## Run the association test

```bash
mkdir -p results/gwas

# Use the PCs from the report as covariates (controls for population structure).
plink2 \
  --pfile results/plink/cohort.qc \
  --pheno data/g1k/pheno/height.pheno --pheno-name height \
  --covar results/plink/cohort.qc.eigenvec --covar-col-nums 3-7 \
  --glm hide-covar \
  --out results/gwas/height
```

`results/gwas/height.height.glm.linear` is the per-variant association summary
(p-value, effect size, SE). Your causal SNP should be the most significant hit.

## Manhattan + QQ plot

```bash
mamba install -n base -c bioconda r-qqman
Rscript - results/gwas/height.height.glm.linear results/gwas <<'R'
library(qqman)
args <- commandArgs(trailingOnly = TRUE)
gwas_file <- args[1]; outdir <- args[2]
gwas <- read.table(gwas_file, header = TRUE, comment.char = "")
gwas <- gwas[!is.na(gwas$P), ]
# Adapt column names: PLINK2 emits #CHROM POS ID REF ALT ... P
gwas$CHR <- as.integer(sub("chr", "", gwas[["X.CHROM"]]))
gwas$BP  <- gwas$POS
gwas$SNP <- gwas$ID
png(file.path(outdir, "manhattan.png"), width = 1200, height = 500)
manhattan(gwas, chr = "CHR", bp = "BP", p = "P", snp = "SNP",
          main = "Simulated height association (chr22)")
dev.off()
png(file.path(outdir, "qq.png"), width = 500, height = 500)
qq(gwas$P, main = "QQ plot")
dev.off()
R
```

## Interpret

* The cohort QC report's PCA should show the EUR sub-populations spread
  along PC1/PC2. If PC1 doesn't capture any of that — your cohort is too
  small or the variant filter dropped too many sites.
* The QQ plot should hug the diagonal everywhere except at the causal SNP
  region. A QQ that's lifted across the board means residual population
  structure — add more PCs as covariates.
* Don't draw conclusions from a 50-sample GWAS. This is a *workflow* sanity
  check; real association studies need thousands of samples.

## Where to go from here

* For polygenic-score work, hand the pgen off to `prsice` or `LDpred2`.
* For colocalisation / fine-mapping, use `coloc` or `SuSiE`.
* For mixed-model GWAS (handles cryptic relatedness instead of dropping
  relateds), use `regenie` or `SAIGE`.
