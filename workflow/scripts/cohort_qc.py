"""
cohort_qc.py — aggregate per-sample QC into a single tidy table + JSON manifest.

Inputs (provided by Snakemake via the rule's `input:` block):
    mosdepth   list[Path]  per-sample mosdepth summary files
    samtools   list[Path]  per-sample `samtools stats` files
    bcftools   list[Path]  per-sample `bcftools stats -s SAMPLE` files
    sex        list[Path]  per-sample JSON from sex_check.py
    cohort_sum Path        GATK CollectVariantCallingMetrics summary

Outputs:
    tsv   Path             one row per sample, columns documented in COLUMNS
    json  Path              structured manifest consumed by make_report.py

Verifies pass/fail against thresholds in qc:* config keys.

NOTE on contamination: the `contamination_freemix` column is intentionally
emitted as NaN and flagged with a TODO. A production deployment should add a
VerifyBamID2 rule and patch its FREEMIX output in here.
"""

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

# No `from __future__ import annotations` here: Snakemake's `script:` directive
# prepends a globals injection block, which would push the future import past
# line 1 and raise SyntaxError. We target Python 3.11+ (qc.yaml pins it) so
# PEP 604 `X | Y` and PEP 585 generic builtins work natively.

snakemake: Any  # type: ignore[no-redef]


# ---- mosdepth ---------------------------------------------------------------

def _parse_mosdepth_summary(path: Path, coverage_threshold: int) -> dict[str, float | None]:
    """
    Pulls the genome-wide mean coverage and the fraction of bases at or above
    `coverage_threshold` from mosdepth's summary text file.
    """
    with path.open() as fh:
        lines = fh.read().splitlines()
    if not lines:
        return {"mean_coverage": None, f"frac_bases_ge_{coverage_threshold}x": None}

    header = lines[0].split("\t")
    mean_idx = header.index("mean") if "mean" in header else None
    rows = [ln.split("\t") for ln in lines[1:] if ln]

    # mosdepth emits 'total' (whole-genome) and 'total_region' rows. Prefer 'total'.
    total_row = next((r for r in rows if r[0] == "total"), None)
    mean = float(total_row[mean_idx]) if (total_row and mean_idx is not None) else None

    # Fraction of bases above threshold: read from the *.thresholds.bed.gz sibling
    # if mosdepth was run with --thresholds. We approximate from the summary
    # row's 'mean' as a coarse proxy if the thresholds file isn't readable here.
    # (The make_report stage uses the bed.gz directly; this column is a quick read.)
    return {
        "mean_coverage": mean,
        f"frac_bases_ge_{coverage_threshold}x": None,  # filled below from thresholds bed
    }


def _parse_mosdepth_thresholds(path: Path, threshold: int) -> float | None:
    """
    Reads <prefix>.thresholds.bed.gz produced by mosdepth --thresholds and
    returns the fraction of windows whose coverage met `threshold`.

    The thresholds file has columns:
        chrom start end region 1X 10X 20X 30X ...
    where each NX column is the number of bases in that window covered >= NX.
    Fraction is sum(NX) / sum(end - start) across all windows.
    """
    import gzip
    if not path.exists():
        return None
    total_bases = 0
    above = 0
    col_name = f"{threshold}X"
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            col_idx = header.index(col_name)
        except ValueError:
            return None
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            try:
                start = int(parts[1])
                end   = int(parts[2])
                bases_above = int(parts[col_idx])
            except (ValueError, IndexError):
                continue
            total_bases += end - start
            above += bases_above
    return (above / total_bases) if total_bases else None


# ---- samtools stats ----------------------------------------------------------

_SN_LINE = re.compile(r"^SN\s+([^:]+):\s+(\S+)")


def _parse_samtools_stats(path: Path) -> dict[str, float]:
    """Extract the SN ('summary numbers') section into a flat dict."""
    out: dict[str, float] = {}
    with path.open() as fh:
        for line in fh:
            m = _SN_LINE.match(line)
            if not m:
                continue
            key, val = m.group(1).strip(), m.group(2)
            try:
                out[key] = float(val)
            except ValueError:
                # non-numeric SN line (rare); skip.
                pass
    return out


# ---- bcftools stats ----------------------------------------------------------

def _parse_bcftools_stats(path: Path, sample: str) -> dict[str, float | None]:
    """
    bcftools stats sections used here:
      'PSC' (per-sample counts): nRefHom, nNonRefHom, nHets, nIndels, nSingletons
      'SN'  (cohort): number of SNPs/indels (we keep ts/tv from 'TSTV')
      'TSTV': cohort Ti/Tv -> we recompute per-sample as nTs/nTv from PSC if avail
    """
    out: dict[str, float | None] = {
        "n_snps": None, "n_indels": None,
        "n_het": None, "n_hom_alt": None, "n_hom_ref": None, "n_singletons": None,
        "het_hom_ratio": None, "titv": None,
    }
    with path.open() as fh:
        for line in fh:
            if line.startswith("PSC\t"):
                parts = line.rstrip("\n").split("\t")
                # PSC  id  sample  nRefHom  nNonRefHom  nHets  nTransitions  nTransversions  nIndels  ...
                if len(parts) < 9 or parts[2] != sample:
                    continue
                try:
                    n_hom_ref = int(parts[3])
                    n_hom_alt = int(parts[4])
                    n_het     = int(parts[5])
                    n_ts      = int(parts[6])
                    n_tv      = int(parts[7])
                    n_indels  = int(parts[8])
                except (ValueError, IndexError):
                    continue
                out["n_hom_ref"] = n_hom_ref
                out["n_hom_alt"] = n_hom_alt
                out["n_het"]     = n_het
                out["n_indels"]  = n_indels
                out["n_snps"]    = n_ts + n_tv
                out["titv"]      = (n_ts / n_tv) if n_tv else None
                out["het_hom_ratio"] = (n_het / n_hom_alt) if n_hom_alt else None
                if len(parts) > 9:
                    try:
                        out["n_singletons"] = int(parts[9])
                    except ValueError:
                        pass
    return out


# ---- GATK cohort metrics ----------------------------------------------------

def _parse_gatk_cohort_summary(path: Path) -> dict[str, float]:
    """
    GATK CollectVariantCallingMetrics writes a small TSV with a header row of
    metric names and one row of values. Returns {metric: value}.
    """
    out: dict[str, float] = {}
    with path.open() as fh:
        header: list[str] | None = None
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("METRICS CLASS"):
                # skip Picard-style metric class banner and blanks
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            if len(parts) != len(header):
                continue
            for k, v in zip(header, parts, strict=True):
                try:
                    out[k] = float(v)
                except ValueError:
                    pass
            break  # only the first data row matters
    return out


# ---- main -------------------------------------------------------------------

COLUMNS = [
    "sample",
    "mean_coverage",
    "frac_bases_above_threshold",
    "total_reads", "mapped_reads", "mapped_pct",
    "duplicate_pct", "insert_size_mean",
    "n_snps", "n_indels", "n_het", "n_hom_alt", "n_singletons",
    "titv", "het_hom_ratio",
    "inferred_sex", "reported_sex", "sex_mismatch",
    "contamination_freemix",
    "pass_coverage", "pass_titv", "pass_het_hom", "pass_sex", "pass_overall",
]


def _build_row(
    sample: str,
    mosdepth_summary: Path,
    mosdepth_thresh: Path,
    samtools_stats: Path,
    bcftools_stats: Path,
    sex_json: Path,
    coverage_threshold: int,
    titv_min: float,
    titv_max: float,
    het_hom_min: float,
    het_hom_max: float,
) -> dict[str, Any]:
    md = _parse_mosdepth_summary(mosdepth_summary, coverage_threshold)
    md_frac = _parse_mosdepth_thresholds(mosdepth_thresh, coverage_threshold)
    st = _parse_samtools_stats(samtools_stats)
    bc = _parse_bcftools_stats(bcftools_stats, sample)
    sx = json.loads(sex_json.read_text())

    total_reads = st.get("raw total sequences")
    mapped      = st.get("reads mapped")
    mapped_pct  = (mapped / total_reads * 100) if (total_reads and mapped) else None
    dup_pct     = (st.get("reads duplicated", 0) / total_reads * 100) if total_reads else None

    pass_cov = md["mean_coverage"] is not None and md["mean_coverage"] >= coverage_threshold
    pass_titv = bc["titv"] is not None and titv_min <= bc["titv"] <= titv_max
    pass_hh   = bc["het_hom_ratio"] is not None and het_hom_min <= bc["het_hom_ratio"] <= het_hom_max
    pass_sex  = not sx.get("mismatch", False)

    row = {
        "sample": sample,
        "mean_coverage": md["mean_coverage"],
        "frac_bases_above_threshold": md_frac,
        "total_reads": total_reads,
        "mapped_reads": mapped,
        "mapped_pct": mapped_pct,
        "duplicate_pct": dup_pct,
        "insert_size_mean": st.get("insert size average"),
        "n_snps": bc["n_snps"],
        "n_indels": bc["n_indels"],
        "n_het": bc["n_het"],
        "n_hom_alt": bc["n_hom_alt"],
        "n_singletons": bc["n_singletons"],
        "titv": bc["titv"],
        "het_hom_ratio": bc["het_hom_ratio"],
        "inferred_sex": sx.get("inferred_sex"),
        "reported_sex": sx.get("reported_sex"),
        "sex_mismatch": sx.get("mismatch"),
        # TODO: replace with VerifyBamID2 FREEMIX value when that rule is added.
        "contamination_freemix": None,
        "pass_coverage": pass_cov,
        "pass_titv":     pass_titv,
        "pass_het_hom":  pass_hh,
        "pass_sex":      pass_sex,
        "pass_overall":  bool(pass_cov and pass_titv and pass_hh and pass_sex),
    }
    return row


def main() -> None:
    samples: list[str]        = list(snakemake.params.samples)
    mosdepth_paths            = [Path(p) for p in snakemake.input.mosdepth]
    samtools_paths            = [Path(p) for p in snakemake.input.samtools]
    bcftools_paths            = [Path(p) for p in snakemake.input.bcftools]
    sex_paths                 = [Path(p) for p in snakemake.input.sex]
    cohort_sum                = Path(snakemake.input.cohort_sum)

    coverage_threshold = int(snakemake.params.coverage_threshold)
    titv_min           = float(snakemake.params.titv_min)
    titv_max           = float(snakemake.params.titv_max)
    het_hom_min        = float(snakemake.params.het_hom_min)
    het_hom_max        = float(snakemake.params.het_hom_max)

    # Build {sample: path} maps from the parallel input lists. Snakemake
    # preserves the iteration order of `expand`, so positional zipping is safe.
    # strict=True makes a length mismatch fail loudly here rather than silently
    # dropping a sample later.
    md_summary = dict(zip(samples, mosdepth_paths, strict=True))
    md_thresh  = {
        s: p.with_name(p.name.replace(".mosdepth.summary.txt", ".thresholds.bed.gz"))
        for s, p in md_summary.items()
    }
    st_map = dict(zip(samples, samtools_paths, strict=True))
    bc_map = dict(zip(samples, bcftools_paths, strict=True))
    sx_map = dict(zip(samples, sex_paths,      strict=True))

    rows = [
        _build_row(
            s,
            md_summary[s], md_thresh[s], st_map[s], bc_map[s], sx_map[s],
            coverage_threshold, titv_min, titv_max, het_hom_min, het_hom_max,
        )
        for s in samples
    ]
    df = pd.DataFrame(rows, columns=COLUMNS)

    tsv_out = Path(snakemake.output.tsv)
    json_out = Path(snakemake.output.json)
    tsv_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(tsv_out, sep="\t", index=False)

    cohort_summary = _parse_gatk_cohort_summary(cohort_sum)
    manifest = {
        "samples": samples,
        "thresholds": {
            "coverage_threshold": coverage_threshold,
            "titv_min": titv_min, "titv_max": titv_max,
            "het_hom_min": het_hom_min, "het_hom_max": het_hom_max,
        },
        "cohort_variant_calling_metrics": cohort_summary,
        "n_pass": int(df["pass_overall"].sum()),
        "n_fail": int((~df["pass_overall"].astype(bool)).sum()),
    }
    json_out.write_text(json.dumps(manifest, indent=2, default=float) + "\n")


if __name__ == "__main__":
    main()
