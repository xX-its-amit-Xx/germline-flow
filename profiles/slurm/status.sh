#!/usr/bin/env bash
# SLURM job status helper for Snakemake's --cluster-status option.
# Maps `sacct` job states onto Snakemake's {success, failed, running} vocabulary.
set -euo pipefail

JOBID="$1"

# Poll sacct; -P=pipe-delimited, -n=no header, -o State only.
state=$(sacct -j "${JOBID}" --format=State --noheader --parsable2 2>/dev/null \
        | head -n 1 | awk -F'+' '{print $1}' | xargs)

case "${state}" in
  ""|PENDING|RUNNING|REQUEUED|CONFIGURING|RESIZING|SUSPENDED)
    echo "running"
    ;;
  COMPLETED)
    echo "success"
    ;;
  BOOT_FAIL|CANCELLED|DEADLINE|FAILED|NODE_FAIL|OUT_OF_MEMORY|PREEMPTED|REVOKED|STOPPED|TIMEOUT)
    echo "failed"
    ;;
  *)
    # Unknown state: be conservative and report failed so Snakemake doesn't hang.
    echo "failed"
    ;;
esac
