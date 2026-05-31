# Recipe 3 — Exome (WES) cohort study

Switching from WGS to WES is a four-line config change plus a couple of QC
threshold tweaks. The point of this recipe is to make the changes explicit so
you don't have to learn them by failing CI three times.

## What changes for WES

* **Calling intervals.** Use the capture-kit BED that came with the kit
  (Twist, Agilent SureSelect, IDT xGen, Illumina TruSeq, …). Pad ±100 bp on
  each side because reads at the edge of a capture region are still useful.
* **Coverage threshold.** WES protocols target 100x mean on-target; the
  germline-flow default of 20x is far too lenient. Use 50-100x.
* **Ti/Tv expected band.** WES is enriched for coding (CpG-rich) regions, so
  expected Ti/Tv is ~3.0-3.3, not the WGS ~2.0-2.1.
* **VQSR.** Needs &ge; 100 WES samples to converge; below that, hard-filter.
* **PLINK MAF.** WES rare-variant cohorts are often run with `maf_min: 0` and
  a separate gnomAD-based annotation step rather than a frequency filter.

## Acquire the capture BED

For a Twist Comprehensive Exome (the most common 2024 default):

```bash
mkdir -p data/wes/reference && cd data/wes/reference

# Twist publishes BEDs on github; pick the one matching your assembly.
# (For other kits, get the BED from the manufacturer's support page.)
curl -sSL -o twist_comprehensive_exome.GRCh38.bed.gz \
  https://www.twistbioscience.com/sites/default/files/resources/2022-12/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed
# (URL is illustrative — Twist's CDN URL changes; check their resources page.)
gunzip -f twist_comprehensive_exome.GRCh38.bed.gz

# Pad ±100 bp.
bedtools slop -i twist_comprehensive_exome.GRCh38.bed \
  -g ../../trio/reference/GRCh38.primary_assembly.genome.fa.fai \
  -b 100 \
  > twist_exome_padded100.bed

# GATK wants .interval_list, not .bed; convert.
gatk BedToIntervalList \
  -I twist_exome_padded100.bed \
  -O twist_exome_padded100.interval_list \
  -SD ../../trio/reference/GRCh38.primary_assembly.genome.dict
```

## Config diff

```yaml
# config/config.yaml
reference:
  fasta: data/trio/reference/GRCh38.primary_assembly.genome.fa
  known_sites:
    - data/trio/reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    - data/trio/reference/Homo_sapiens_assembly38.dbsnp138.vcf
  calling_intervals: data/wes/reference/twist_exome_padded100.interval_list

joint_calling:
  # 20 samples -> below the WES VQSR minimum. Use hard filters.
  # For >=100 WES samples, switch to "vqsr".
  filter_mode: hard_filter

qc:
  coverage_threshold: 50          # WES target: 100x mean, so 50x as a floor is reasonable
  titv_min: 2.8                   # WES expected ~3.0-3.3
  titv_max: 3.6
  het_hom_min: 1.4
  het_hom_max: 2.5

plink_export:
  call_rate_min: 0.95
  maf_min: 0.0                    # keep singletons for rare-variant analysis
  hwe_p_min: 1.0e-10              # very loose; WES rare variants violate HWE often
  king_cutoff: 0.0884
```

## Run

```bash
snakemake --use-conda --cores 16
```

Expected wall-clock for 20 WES samples on 16 cores: ~6-10 h, mostly in
HaplotypeCaller.

## Interpret the report

* **Coverage** &mdash; the bar plot shows the mean on-target coverage. The
  `frac_bases_above_threshold` column tells you what fraction of the capture
  was covered at &ge; 50x. Aim for &ge; 0.90 on a well-prepped Twist kit.
* **Ti/Tv** &mdash; should be 2.9-3.3 across all samples. A sample at 2.0
  has off-target reads dominating the call set (capture failed for that
  sample); a sample at 4.0 has aggressive filtering or low-coverage noise.
* **PCA** &mdash; WES PCs are noisier than WGS because of capture-region
  ascertainment; don't over-interpret PC3+.

## When to switch to VQSR

If your cohort grows past 100 WES samples, flip `filter_mode: vqsr` and run.
You'll see the report's filter column populated with "PASS" / "VQSRTrancheSNP\*"
tranches instead of the hard-filter labels. VQSR almost always recovers true
indels that hard-filtering drops — but only on big-enough cohorts.

## Going further

* For rare-variant burden tests, annotate with VEP (see recipe 4) then run
  `regenie` or `SAIGE-GENE+` aggregating by gene.
* For Mendelian-disease prioritisation, hand off to
  [OpenCRAVAT](https://opencravat.org/) — see recipe 4.
