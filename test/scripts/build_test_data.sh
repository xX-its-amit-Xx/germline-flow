#!/usr/bin/env bash
# =============================================================================
# build_test_data.sh
#
# Idempotently builds the test fixtures used by `make test` and CI:
#   1. A 5 Mb slice of GRCh38 chr20 (chr20:1-5,000,000) FASTA, plus all the
#      BWA-MEM2 / samtools / GATK index files derived from it.
#   2. Two tiny known-sites VCFs (dbSNP placeholder + indels placeholder),
#      generated from a handful of randomly chosen positions in the slice.
#   3. Three paired-end FASTQ sets (~5x coverage) simulated with `wgsim`.
#
# The script reuses any existing artefact, so re-running is cheap.
# Total runtime cold: ~2-3 min on a laptop. ~150 MB on disk.
#
# Required tools (installed via the workflow/envs/lint.yaml mamba env, which
# pulls in samtools / bcftools / wgsim transitively; otherwise install them
# from bioconda):
#   - samtools / bcftools / tabix
#   - bwa-mem2
#   - gatk4 (for CreateSequenceDictionary)
#   - wgsim
#   - curl
# =============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/../.." && pwd)"
DATA="${ROOT}/test/data"
REF_DIR="${DATA}/reference"
FQ_DIR="${DATA}/fastq"

mkdir -p "${REF_DIR}" "${FQ_DIR}"

REGION="chr20:1-5000000"
SLICE_FA="${REF_DIR}/chr20_slice.fa"
SLICE_FAI="${SLICE_FA}.fai"
SLICE_DICT="${REF_DIR}/chr20_slice.dict"
SLICE_BWA="${SLICE_FA}.bwt.2bit.64"
SLICE_INTERVALS="${REF_DIR}/chr20_slice.interval_list"

# Public GRCh38 primary assembly FASTA (~3 GB). We only need chr20 so we stream
# samtools faidx against a remote-indexed copy when available, but most reliably
# we fall back to a small chr20-only mirror.
GRCH38_CHR20_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz"

# Tools we need on PATH.
need() { command -v "$1" >/dev/null 2>&1 || { echo "missing tool: $1" >&2; exit 1; }; }
need samtools
need bcftools
need bwa-mem2
need gatk
need wgsim
need curl
need tabix

echo "==> [1/4] chr20 reference slice"

if [[ ! -s "${SLICE_FA}" ]]; then
  tmp_chr20="${REF_DIR}/chr20.fa.gz"
  if [[ ! -s "${tmp_chr20}" ]]; then
    echo "    downloading ${GRCH38_CHR20_URL}"
    curl -sSL --retry 3 "${GRCH38_CHR20_URL}" -o "${tmp_chr20}"
  fi
  echo "    extracting ${REGION}"
  # Decompress, index, slice.
  gunzip -k -f "${tmp_chr20}"
  samtools faidx "${REF_DIR}/chr20.fa"
  samtools faidx "${REF_DIR}/chr20.fa" "${REGION}" \
    | awk 'BEGIN{first=1} /^>/{if(first){print ">chr20"; first=0; next} else next} {print}' \
    > "${SLICE_FA}"
  rm -f "${REF_DIR}/chr20.fa" "${REF_DIR}/chr20.fa.fai" "${tmp_chr20}"
fi

if [[ ! -s "${SLICE_FAI}" ]]; then
  samtools faidx "${SLICE_FA}"
fi

if [[ ! -s "${SLICE_DICT}" ]]; then
  gatk CreateSequenceDictionary -R "${SLICE_FA}" -O "${SLICE_DICT}" > /dev/null
fi

if [[ ! -s "${SLICE_BWA}" ]]; then
  echo "    bwa-mem2 index"
  bwa-mem2 index "${SLICE_FA}" 2>/dev/null
fi

if [[ ! -s "${SLICE_INTERVALS}" ]]; then
  # GATK .interval_list = SAM header + 1-based intervals. We can derive both
  # from the dict (which is itself a SAM header).
  {
    cat "${SLICE_DICT}"
    awk 'BEGIN{OFS="\t"} {print $1, 1, $2, "+", "chr20_slice"}' "${SLICE_FAI}"
  } > "${SLICE_INTERVALS}"
fi

echo "==> [2/4] placeholder known-sites VCFs"
# Pick 50 SNP positions and 10 indel positions deterministically. These exist
# only so BQSR has a non-empty known-sites file; they don't need to be real
# polymorphisms.
DBSNP_VCF="${REF_DIR}/dbsnp.chr20_slice.vcf.gz"
INDEL_VCF="${REF_DIR}/known_indels.chr20_slice.vcf.gz"

if [[ ! -s "${DBSNP_VCF}" ]]; then
  python3 - "${SLICE_FA}" "${REF_DIR}/dbsnp.chr20_slice.vcf" <<'PY'
import random, sys, textwrap
from pathlib import Path

fa = Path(sys.argv[1]).read_text().splitlines()
seq = "".join(line for line in fa if not line.startswith(">"))
random.seed(42)
positions = sorted(random.sample(range(1000, len(seq) - 1000, 1000), 50))
header = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership placeholder">
    ##contig=<ID=chr20,length={L}>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
""".format(L=len(seq)))
records = []
swap = {"A": "G", "C": "T", "G": "A", "T": "C", "N": "A"}
for i, p in enumerate(positions, 1):
    ref = seq[p - 1].upper()
    if ref == "N":
        continue
    alt = swap[ref]
    records.append(f"chr20\t{p}\trs{i}\t{ref}\t{alt}\t.\tPASS\tDB")
Path(sys.argv[2]).write_text(header + "\n".join(records) + "\n")
PY
  bgzip -f "${REF_DIR}/dbsnp.chr20_slice.vcf"
  tabix -p vcf "${DBSNP_VCF}"
fi

if [[ ! -s "${INDEL_VCF}" ]]; then
  python3 - "${SLICE_FA}" "${REF_DIR}/known_indels.chr20_slice.vcf" <<'PY'
import random, sys, textwrap
from pathlib import Path

fa = Path(sys.argv[1]).read_text().splitlines()
seq = "".join(line for line in fa if not line.startswith(">"))
random.seed(7)
positions = sorted(random.sample(range(2000, len(seq) - 2000, 4000), 10))
header = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##contig=<ID=chr20,length={L}>
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
""".format(L=len(seq)))
records = []
for i, p in enumerate(positions, 1):
    ref = seq[p - 1:p + 1].upper()
    if "N" in ref:
        continue
    # 1-bp deletion: REF=AB, ALT=A
    alt = ref[0]
    records.append(f"chr20\t{p}\tindel{i}\t{ref}\t{alt}\t.\tPASS\t.")
Path(sys.argv[2]).write_text(header + "\n".join(records) + "\n")
PY
  bgzip -f "${REF_DIR}/known_indels.chr20_slice.vcf"
  tabix -p vcf "${INDEL_VCF}"
fi

echo "==> [3/4] simulating paired-end FASTQs (wgsim)"
# ~5x coverage on 5 Mb = ~167k reads of 150 bp per sample. wgsim is fast.
N_READS=170000
READ_LEN=150
INSERT=400

for i in 01 02 03; do
  R1="${FQ_DIR}/SIM${i}_R1.fastq"
  R2="${FQ_DIR}/SIM${i}_R2.fastq"
  if [[ -s "${R1}.gz" && -s "${R2}.gz" ]]; then
    continue
  fi
  seed=$((100 + 10#$i))
  wgsim \
      -1 "${READ_LEN}" -2 "${READ_LEN}" \
      -N "${N_READS}" \
      -d "${INSERT}" -s 50 \
      -r 0.001 -R 0.10 -X 0.30 \
      -S "${seed}" \
      "${SLICE_FA}" "${R1}" "${R2}" \
      > "${FQ_DIR}/SIM${i}.wgsim.log" 2>&1
  gzip -f "${R1}" "${R2}"
done

echo "==> [4/4] done"
ls -lh "${REF_DIR}" "${FQ_DIR}" | sed 's/^/    /'
