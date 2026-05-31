# Cookbook

End-to-end worked examples with real public data. Each recipe is fully
reproducible — every URL, every command, every config diff is in the file.

| # | Recipe | Data | Why you'd read it |
|---|--------|------|-------------------|
| 1 | [Trio analysis (NA12878 family)](01-trio-analysis.md) | GIAB NA12878 / NA12891 / NA12892 chr20 | Run a parent-child-parent trio, check Mendelian-inheritance consistency, benchmark against the GIAB high-confidence truth set with `rtg vcfeval`. |
| 2 | [Small-cohort GWAS pilot](02-small-cohort-gwas.md) | 1000 Genomes phase 3 chr22, 50 EUR samples | Joint-call a small cohort, run PCA, drop relateds, hand the PLINK2 pgen off to `plink2 --glm` for an association sketch. |
| 3 | [Exome cohort (WES) study](03-exome-cohort.md) | Synthetic 20-sample WES with Twist exome BED | Switch the intervals file to a capture-kit BED, tune QC thresholds for WES Ti/Tv, choose between VQSR and hard-filtering. |
| 4 | [Integrations: VEP, snpEff, Hail, OpenCRAVAT, vcfeval](04-integrations.md) | Output of any of the above | Wire the germline-flow output into the most common downstream tools — annotation, scalable analysis, clinical prioritisation, accuracy benchmarking. |
| 5 | [Singleton / clinical-style case](05-singleton-clinical.md) | One proband, no family | The "I got a single sample on my desk" workflow — disable VQSR, restrict to gene panels, hand off to OpenCRAVAT for ACMG annotations. Includes the obligatory "this is not clinical software" reminder. |

## How to use a recipe

Each recipe is a self-contained Markdown file with four sections:

1. **What you'll have at the end** — concrete outputs, with expected size.
2. **Data acquisition** — `curl`/`wget`/`s3` commands for the public data,
   including SHA256s where the upstream publishes them.
3. **Config diff** — exactly what to change in `config/config.yaml` and
   `config/samples.tsv`. Diffs only, not the whole file.
4. **Run + interpret** — the `snakemake` invocation, what the report should
   look like, and what to do if it doesn't.

## Want to contribute a recipe?

Open a PR adding `cookbook/NN-your-recipe.md` and a row in the table above.
The bar is: a graduate student should be able to follow it end-to-end on a
laptop or a single SLURM partition without asking anybody questions.
