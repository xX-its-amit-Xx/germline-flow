# test/

A self-contained end-to-end smoke test for germline-flow. CI runs it on every
push and PR via [.github/workflows/ci.yml](../.github/workflows/ci.yml). You can
also run it locally:

```bash
bash test/scripts/build_test_data.sh    # ~3 min, ~150 MB of intermediate files
make test CORES=4                       # ~8-12 min on a laptop
bash test/scripts/check_outputs.sh
```

## What's in here

| File | Purpose |
|------|---------|
| `config.yaml`              | Test-only config that points at a tiny chr20 slice + the synthetic FASTQ set. |
| `samples.tsv`              | Three pseudo-samples (`SIM01-03`) used for the test cohort. |
| `scripts/build_test_data.sh` | Downloads a 5 Mb chr20 slice of GRCh38, builds the BWA-MEM2 / samtools / GATK indexes, and uses `wgsim` to simulate paired-end reads at 5x coverage per sample. Idempotent. |
| `scripts/check_outputs.sh` | Asserts that the final VCF, PLINK2 fileset, and HTML report exist and are non-empty. |
| `data/`                    | Created by `build_test_data.sh`. **gitignored** because it's ~150 MB. |

## Why synthetic data?

Public real data (NA12878 FASTQ, GIAB high-confidence VCF, dbSNP) are large
(>>2 GB) and live on mirrors that occasionally rate-limit or go down. We use
`wgsim` to simulate three paired-end read sets from a small slice of the real
GRCh38 reference. That gives us:

* deterministic CI runtime,
* no external mirror dependency,
* enough variation across the three samples to exercise joint genotyping,
  PCA, and KING relatedness,
* output we can fully validate against the reference we simulated from.

## Want to test against real data?

See [`cookbook/01-trio-analysis.md`](../cookbook/01-trio-analysis.md). It walks
through running the pipeline on the GIAB NA12878 trio.
