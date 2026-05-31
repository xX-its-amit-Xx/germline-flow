# Recipe 5 — Singleton "clinical-style" workflow

You have one proband, no parents, no cohort. Maybe a clinician handed you a
FASTQ and asked "any obvious findings?". This recipe runs germline-flow on a
single sample, restricts to a gene panel, and hands off to OpenCRAVAT for
ACMG-style annotation.

> **Read this first.** Single-sample variant calling without family or cohort
> context is *much* less accurate than joint-calling. This recipe is for
> research/education only. **It is not a clinical workflow.** Findings must be
> revalidated in a CAP/CLIA (or local equivalent) accredited laboratory before
> they're used to make any medical decision. See the project [disclaimer](../README.md#disclaimer).

## What you'll have at the end

* A single-sample joint-calling-style VCF restricted to your gene panel
  (~few hundred PASS variants).
* An OpenCRAVAT spreadsheet with ACMG classification, ClinVar significance,
  gnomAD frequencies, and CADD/REVEL scores.
* A coverage report showing which panel exons are under-covered (the
  most common cause of a "we didn't find the variant" false negative).

## Acquire a gene panel

Pick a published panel for the indication. For hereditary cardiomyopathy, the
[ClinGen HCM panel](https://search.clinicalgenome.org/kb/gene-validity?page=1&size=25&search=HCM)
genes are a reasonable starting point. Build a BED of their coding exons +
splice sites from Ensembl Biomart.

```bash
mkdir -p data/singleton/panel && cd data/singleton/panel

# Example: 8-gene HCM panel BED. Generate via Biomart or ucsc-table-browser.
cat > hcm_panel.genes <<'EOF'
MYH7
MYBPC3
TNNT2
TNNI3
TPM1
ACTC1
MYL2
MYL3
EOF

# Convert to a padded exon BED (±20 bp for splice sites). Use mysql or
# the UCSC table browser to pull exon coordinates; here's a one-liner with
# the UCSC API:
for gene in $(cat hcm_panel.genes); do
  curl -sSL "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeq;chrom=auto;name=${gene}"
done > hcm_panel.raw.json

# (In practice: use 'gencode' or 'mane' from biomaRt or
# https://www.ensembl.org/biomart/martview/ for the cleanest exon BED.
# Pad ±20 bp for splice acceptors / donors.)

# Convert to GATK .interval_list (assumes hcm_panel_padded20.bed already exists):
gatk BedToIntervalList \
    -I hcm_panel_padded20.bed \
    -O hcm_panel_padded20.interval_list \
    -SD ../../trio/reference/GRCh38.primary_assembly.genome.dict
```

## Config diff

```yaml
# config/config.yaml
samples: config/singleton.tsv

reference:
  fasta: data/trio/reference/GRCh38.primary_assembly.genome.fa
  known_sites:
    - data/trio/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    - data/trio/reference/Homo_sapiens_assembly38.dbsnp138.vcf
  calling_intervals: data/singleton/panel/hcm_panel_padded20.interval_list

joint_calling:
  # Hard-filter is mandatory for single-sample (VQSR needs cohort).
  filter_mode: hard_filter

qc:
  coverage_threshold: 30          # WES-style: target ~50-100x, alert at 30x
  titv_min: 2.5                   # WES-restricted panel; expect ~3.0
  titv_max: 3.5

plink_export:
  # PLINK QC is unhelpful on N=1. Disable variant filters but keep the export
  # for downstream PRS tooling.
  call_rate_min: 0.0
  maf_min: 0.0
  hwe_p_min: 0.0
  king_cutoff: 1.0                # effectively disable --king-cutoff
```

`config/singleton.tsv`:

```tsv
sample	fastq_1	fastq_2	sex
PROBAND01	data/singleton/fastq/PROBAND01_R1.fastq.gz	data/singleton/fastq/PROBAND01_R2.fastq.gz	U
```

## Run

```bash
snakemake --use-conda --cores 8
```

## Coverage QC on the panel

Open `results/report/cohort_qc.html`. The "coverage" figure tells you mean
coverage across the panel. To get per-exon coverage (essential for clinical
panels — you must show every coding base was &ge; 20x or report it as a gap):

```bash
mosdepth -t 4 --no-per-base --by data/singleton/panel/hcm_panel_padded20.bed \
    --thresholds 1,10,20,30 \
    results/qc/mosdepth_panel/PROBAND01 results/bqsr/PROBAND01.bqsr.bam

# Find any panel exon < 20x.
zcat results/qc/mosdepth_panel/PROBAND01.thresholds.bed.gz \
  | awk -F'\t' 'NR>1 && $5/($3-$2) < 0.95 { print $1, $2, $3, $4 }' \
  > results/qc/PROBAND01.panel_gaps.bed
```

A gap region usually means:

* the capture probe failed (re-design or re-sequence);
* the region is in a pseudogene family (rerun with a paralog-aware caller like
  `dragen-os` or `octopus`);
* coverage is low everywhere (resequence at higher depth).

## Hand off to OpenCRAVAT

```bash
mamba create -n cravat -c bioconda -c conda-forge open-cravat=2.4.2
mamba run -n cravat oc module install-base
mamba run -n cravat oc module install clinvar gnomad3 cadd revel \
    omim acmg cancer_genome_interpreter spliceai

mamba run -n cravat oc run results/joint/cohort.filtered.vcf.gz \
    -l hg38 \
    -a clinvar gnomad3 cadd revel omim acmg spliceai \
    -t excel csv vcf \
    -d results/cravat
```

Open `results/cravat/cohort.filtered.vcf.gz.cravat.html` and filter to
`ACMG = "Pathogenic" or "Likely pathogenic"` and `gnomAD AF < 0.001`. That
shortlist is the *starting point* for a clinical review — not the answer.

## Interpret

Every actionable variant should be:

1. Visually confirmed in IGV against the BAM (`results/bqsr/PROBAND01.bqsr.bam`).
2. Re-checked in a second tool (e.g., DeepVariant or DRAGEN) on the same FASTQ.
3. Confirmed by orthogonal assay (Sanger, ddPCR) in an accredited lab.
4. Interpreted by a clinical geneticist with the patient's phenotype and
   family history in front of them.

Anything short of those four steps is not a clinical finding.

## What this recipe does NOT do (intentionally)

* Trio-based de novo calling — use `gatk CalculateGenotypePosteriors` with a
  pedigree and a population-AF VCF, or DeepTrio.
* Structural variant calling — use `manta` or `delly`.
* Copy-number / aneuploidy detection — use `gatk gCNV` or `XHMM`.
* Mosaicism detection — needs &ge;500x targeted sequencing and a specialised
  caller like `Mutect2` in tumour-only mode.

For any of those, germline-flow is the wrong tool. Use the linked specialised
caller and join the results to germline-flow's output on (CHROM, POS).
