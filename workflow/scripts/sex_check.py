"""
sex_check.py — infer genetic sex from chrX and chrY coverage.

Reads the mosdepth `*.mosdepth.summary.txt` file, computes mean coverage on
chrX and chrY (handling both "chr"-prefixed and unprefixed contig naming),
and compares the chrY/chrX ratio to a simple cutoff:

    ratio < 0.1   -> XX (female)
    ratio > 0.3   -> XY (male)
    otherwise     -> ambiguous (sex chromosome aneuploidy, low coverage, or
                    cross-sample contamination)

A "reported_sex" from the sample sheet, if provided, is recorded for
later mismatch flagging in the cohort report.
"""

import json
from pathlib import Path
from typing import Any

# Snakemake's `script:` directive injects a globals block at the top of the
# script, which means `from __future__ import annotations` would land after
# it and fail with SyntaxError. We target Python 3.11+ (qc.yaml pins it), so
# PEP 604 `X | Y` unions and PEP 585 generic builtins work natively without
# the future import.

# Snakemake injects the global `snakemake` object when running via `script:`.
snakemake: Any  # type: ignore[no-redef]

XX_MAX_YX_RATIO = 0.1
XY_MIN_YX_RATIO = 0.3


def _norm_contig(name: str) -> str:
    n = name.lower()
    if n.startswith("chr"):
        n = n[3:]
    return n


def parse_mosdepth_summary(path: Path) -> dict[str, float]:
    """
    Return {contig: mean_coverage} for whole-chromosome rows only.

    mosdepth's summary file has rows like:
        chrom  length      bases       mean  min  max
        chr1   248956422   123456789   8.40  0    250
        chr1_region 248956422 ...
    We keep only rows whose chrom name has no "_" suffix (chr1, not chr1_region).
    """
    out: dict[str, float] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            mean_idx = header.index("mean")
        except ValueError as exc:
            raise RuntimeError(f"unexpected mosdepth header: {header}") from exc
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if not parts or len(parts) < len(header):
                continue
            chrom = parts[0]
            # Skip per-region rows: mosdepth appends "_region" / "_total".
            if chrom.endswith(("_region", "_total")) or chrom == "total":
                continue
            try:
                mean = float(parts[mean_idx])
            except ValueError:
                continue
            out[_norm_contig(chrom)] = mean
    return out


def infer_sex(cov: dict[str, float]) -> tuple[str, float | None]:
    x = cov.get("x")
    y = cov.get("y")
    if x is None or x <= 0:
        return "U", None
    if y is None:
        return "F", 0.0
    ratio = y / x
    if ratio < XX_MAX_YX_RATIO:
        inferred = "F"
    elif ratio > XY_MIN_YX_RATIO:
        inferred = "M"
    else:
        inferred = "U"
    return inferred, ratio


def main() -> None:
    summary_path = Path(snakemake.input.summary)
    out_path = Path(snakemake.output.json)
    reported = str(snakemake.params.reported_sex or "U").upper()
    sample = snakemake.wildcards.sample

    cov = parse_mosdepth_summary(summary_path)
    inferred, ratio = infer_sex(cov)
    mismatch = (reported in {"F", "M"}) and (inferred in {"F", "M"}) and (inferred != reported)

    payload = {
        "sample": sample,
        "reported_sex": reported,
        "inferred_sex": inferred,
        "chrX_mean_coverage": cov.get("x"),
        "chrY_mean_coverage": cov.get("y"),
        "chrY_over_chrX_ratio": ratio,
        "mismatch": mismatch,
        "notes": (
            "ratio<0.1 => XX; ratio>0.3 => XY; otherwise ambiguous "
            "(low coverage, sex-chromosome aneuploidy, or contamination)"
        ),
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
