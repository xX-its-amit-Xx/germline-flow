#!/usr/bin/env bash
# Assert the pipeline produced its core artefacts and they're non-empty.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"

REQUIRED=(
  "results/joint/cohort.filtered.vcf.gz"
  "results/joint/cohort.filtered.vcf.gz.tbi"
  "results/plink/cohort.qc.pgen"
  "results/plink/cohort.qc.pvar"
  "results/plink/cohort.qc.psam"
  "results/plink/cohort.qc.eigenvec"
  "results/qc/cohort_qc.tsv"
  "results/qc/cohort_qc.json"
  "results/report/cohort_qc.html"
)

fail=0
for f in "${REQUIRED[@]}"; do
  if [[ ! -s "${f}" ]]; then
    echo "MISSING or empty: ${f}" >&2
    fail=1
  else
    printf 'OK  %s  (%s)\n' "${f}" "$(du -h "${f}" | cut -f1)"
  fi
done

# Light sanity on the report: must mention each sample.
if [[ -s "results/report/cohort_qc.html" ]]; then
  for s in SIM01 SIM02 SIM03; do
    if ! grep -q "${s}" "results/report/cohort_qc.html"; then
      echo "report missing sample ${s}" >&2
      fail=1
    fi
  done
fi

if [[ "${fail}" -ne 0 ]]; then
  echo "FAIL: required outputs missing or report incomplete" >&2
  exit 1
fi
echo "PASS: all expected outputs present"
