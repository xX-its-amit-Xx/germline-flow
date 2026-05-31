# Recipe 4 — Integrating germline-flow output with downstream tools

germline-flow stops at "analysis-ready joint VCF + PLINK2 pgen". This recipe
shows how to plug that output into the most common things you'll want next.

```
                     +------------------- VEP / snpEff (functional annotation)
                     |
                     +------------------- OpenCRAVAT (clinical prioritisation,
                     |                                ACMG, ClinVar, gnomAD)
results/joint/       +------------------- Hail / glow (scalable analysis at
  cohort.filtered.vcf.gz                              cohort-of-thousands scale)
                     |
                     +------------------- rtg vcfeval (benchmark vs GIAB)
                     |
                     +------------------- bcftools / vcfanno (lightweight
                                                              annotation)

results/plink/cohort.qc.{pgen,pvar,psam}
                     +------------------- plink2 --glm / regenie / SAIGE (GWAS)
                     +------------------- LDpred2 / PRScs / PRSice (PRS)
                     +------------------- ADMIXTURE / fastSTRUCTURE (ancestry)
```

---

## 4.1 VEP — Ensembl Variant Effect Predictor

Adds gene/transcript/consequence/CADD/SIFT/PolyPhen/etc annotations to each
variant. The annotated VCF is the standard input for variant-level filtering.

```bash
# Install via the bioconda channel.
mamba install -n vep -c bioconda ensembl-vep=112
mamba run -n vep vep_install \
    --AUTO cfp --SPECIES homo_sapiens --ASSEMBLY GRCh38 \
    --CACHEDIR data/vep_cache --PLUGINS CADD,REVEL,SpliceAI

mamba run -n vep vep \
    -i results/joint/cohort.filtered.vcf.gz \
    -o results/annot/cohort.vep.vcf.gz \
    --vcf --compress_output bgzip \
    --cache --dir_cache data/vep_cache \
    --assembly GRCh38 --fasta data/trio/reference/GRCh38.primary_assembly.genome.fa \
    --everything \
    --plugin CADD,/path/to/CADD/whole_genome_SNVs.tsv.gz,/path/to/CADD/InDels.tsv.gz
```

Use `--pick` for one consequence per variant (recommended for clinical
filtering); omit it for full transcript-level output (better for research).

---

## 4.2 snpEff — alternative annotator

Lighter-weight than VEP, no cache to manage. Good for batch annotation in
container pipelines where VEP's plugin ecosystem is overkill.

```bash
mamba install -n snpeff -c bioconda snpeff=5.2a
mamba run -n snpeff snpEff -v -dataDir data/snpeff_data GRCh38.105 \
    results/joint/cohort.filtered.vcf.gz \
    > results/annot/cohort.snpeff.vcf
bgzip results/annot/cohort.snpeff.vcf && tabix -p vcf results/annot/cohort.snpeff.vcf.gz
```

snpEff also emits a per-cohort HTML report — `snpEff_summary.html` — that
complements germline-flow's cohort QC report nicely.

---

## 4.3 OpenCRAVAT — clinical / Mendelian prioritisation

Annotates with ClinVar, ACMG criteria, gene-disease relationships, PharmGKB,
and ~100 other databases through a plugin system. Used in clinical-research
workflows for Mendelian disease.

```bash
mamba create -n cravat -c bioconda -c conda-forge open-cravat=2.4.2
mamba run -n cravat oc module install-base
mamba run -n cravat oc module install clinvar gnomad3 cadd revel \
    omim acmg cancer_genome_interpreter

mamba run -n cravat oc run results/joint/cohort.filtered.vcf.gz \
    -l hg38 \
    -a clinvar gnomad3 cadd revel omim acmg \
    -t excel csv vcf \
    -d results/cravat
```

Open `results/cravat/cohort.filtered.vcf.gz.cravat.html` for the interactive
viewer — sorts by ACMG classification, lets you filter by inheritance
pattern, etc.

---

## 4.4 Hail — scale to cohorts of thousands

When the cohort gets large enough that PLINK2 starts to feel slow (~10k+
samples), move to [Hail](https://hail.is/) which runs on Spark.

```python
import hail as hl
hl.init(default_reference="GRCh38")
mt = hl.import_vcf(
    "results/joint/cohort.filtered.vcf.gz",
    reference_genome="GRCh38",
    force_bgz=True,
).write("results/hail/cohort.mt", overwrite=True)

mt = hl.read_matrix_table("results/hail/cohort.mt")
mt = hl.sample_qc(mt)
mt = hl.variant_qc(mt)
mt.cols().show(5)
```

Hail's `sample_qc` and `variant_qc` produce a superset of what germline-flow's
PLINK QC computes; the report metrics line up exactly so you can cross-check.

---

## 4.5 vcfeval — benchmarking against a truth set

If you have a sample with a published truth VCF (any GIAB sample), vcfeval is
the gold-standard for precision/recall.

```bash
mamba install -n rtg -c bioconda rtg-tools=3.12.1

# One-time: build the SDF for the reference.
mamba run -n rtg rtg format \
    -o data/trio/reference/GRCh38.sdf \
    data/trio/reference/GRCh38.primary_assembly.genome.fa

# Per-sample evaluation.
mamba run -n rtg rtg vcfeval \
    --baseline=data/giab/HG001_GRCh38.high_confidence.vcf.gz \
    --bed-regions=data/giab/HG001_GRCh38.high_confidence.bed \
    --calls=results/joint/cohort.filtered.vcf.gz \
    --template=data/trio/reference/GRCh38.sdf \
    --output=results/benchmark/HG001 \
    --sample=NA12878
```

Compare the `summary.txt` precision/recall against the GA4GH benchmarking
toolkit's reference numbers at
https://github.com/ga4gh/benchmarking-tools/blob/master/resources/expected-results/.

---

## 4.6 bcftools / vcfanno — lightweight annotation in the pipeline

If you don't want VEP's heavyweight cache, you can annotate inline with
`bcftools annotate` or `vcfanno` (single-pass merge of many BED/VCF tracks).
This is the fastest way to add gnomAD AF, ClinVar significance, and conservation
scores during the pipeline run itself.

```bash
mamba install -n vcfanno -c bioconda vcfanno=0.3.5

cat > vcfanno.toml <<'EOF'
[[annotation]]
file    = "data/annot/gnomad.v3.1.2.sites.AF.vcf.gz"
fields  = ["AF"]
names   = ["gnomad_AF"]
ops     = ["self"]

[[annotation]]
file    = "data/annot/clinvar.vcf.gz"
fields  = ["CLNSIG", "CLNDN"]
names   = ["clinvar_sig", "clinvar_dn"]
ops     = ["self", "self"]
EOF

vcfanno vcfanno.toml results/joint/cohort.filtered.vcf.gz \
    | bgzip > results/annot/cohort.annot.vcf.gz
tabix -p vcf results/annot/cohort.annot.vcf.gz
```

---

## 4.7 ADMIXTURE — global ancestry inference

```bash
mamba install -n admixture -c bioconda admixture=1.3.0

# ADMIXTURE wants PLINK1 bed/bim/fam.
plink2 --pfile results/plink/cohort.qc \
       --make-bed --out results/admixture/cohort.qc

mamba run -n admixture admixture --cv \
    results/admixture/cohort.qc.bed 5 \
    > results/admixture/cohort.K5.log
```

The `.Q` matrix has one row per sample and K columns = ancestry fractions.
Drop into the report by post-processing.

---

## 4.8 LDpred2 / PRScs — polygenic risk scores

Hand the PLINK2 pgen to your favourite PRS tool. See the
[LDpred2 vignette](https://privefl.github.io/bigsnpr/articles/LDpred2.html)
for the canonical example.

---

## Putting it all together

A typical research-lab "Monday morning" downstream stack:

1. germline-flow → joint VCF + pgen.
2. VEP + vcfanno → annotated VCF with gnomAD/ClinVar.
3. OpenCRAVAT → clinical/ACMG triage spreadsheet for the wet-lab collaborator.
4. PLINK2 → PCA + KING (already in the report).
5. ADMIXTURE → global ancestry fractions.
6. plink2 --glm / regenie → association statistics for the phenotype of interest.
7. rtg vcfeval → benchmark against any GIAB controls in the cohort.

Each of those is one command and one input. The point of stopping germline-flow
at the joint VCF / pgen step is to keep that interface clean.
