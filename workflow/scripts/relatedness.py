"""
relatedness.py — turn PLINK2 KING output into a report-friendly table.

PLINK2 emits:
  * <prefix>.kin0          tab-separated pairwise KING-robust kinship
  * <prefix>.king.cutoff.{in,out}.id   IDs kept / removed by --king-cutoff

KING kinship coefficient interpretation (Manichaikul et al. 2010):
    > 0.354     duplicate / monozygotic twin
    0.177-0.354 1st-degree (parent-child, full siblings)
    0.0884-0.177 2nd-degree (half-sib, grandparent, avuncular)
    0.0442-0.0884 3rd-degree (first cousins)
    < 0.0442    unrelated
"""

import json
from pathlib import Path
from typing import Any

import pandas as pd

# `from __future__ import annotations` cannot be used here because Snakemake's
# `script:` directive prepends its globals injection, pushing any __future__
# imports past line 1. We target Python 3.11+ so the modern type syntax works
# natively without it.

snakemake: Any  # type: ignore[no-redef]


def classify(kinship: float) -> str:
    if kinship > 0.354:
        return "duplicate_or_MZ_twin"
    if kinship > 0.177:
        return "1st_degree"
    if kinship > 0.0884:
        return "2nd_degree"
    if kinship > 0.0442:
        return "3rd_degree"
    return "unrelated"


def main() -> None:
    kin_path = Path(snakemake.input.kinship)
    tsv_out  = Path(snakemake.output.tsv)
    json_out = Path(snakemake.output.json)
    cutoff   = float(snakemake.params.cutoff)
    # snakemake.input.psam is declared so the rule depends on it being current,
    # but the kinship table already carries the IDs we need.

    tsv_out.parent.mkdir(parents=True, exist_ok=True)

    if not kin_path.exists() or kin_path.stat().st_size == 0:
        # PLINK2 omits the file entirely when no pair exceeds its internal
        # threshold; treat that as "all unrelated".
        empty = pd.DataFrame(columns=["ID1", "ID2", "KINSHIP", "relationship", "above_cutoff"])
        empty.to_csv(tsv_out, sep="\t", index=False)
        json_out.write_text(json.dumps({
            "cutoff": cutoff,
            "n_pairs_above_cutoff": 0,
            "pairs": [],
        }, indent=2) + "\n")
        return

    df = pd.read_csv(kin_path, sep="\t")
    # PLINK2 column casing varies across versions; normalise.
    df.columns = [c.replace("#", "").strip() for c in df.columns]
    # Standard PLINK2 --make-king-table columns: FID1 IID1 FID2 IID2 NSNP HETHET IBS0 KINSHIP
    id1_col = "IID1" if "IID1" in df.columns else df.columns[1]
    id2_col = "IID2" if "IID2" in df.columns else df.columns[3]
    kin_col = "KINSHIP" if "KINSHIP" in df.columns else df.columns[-1]

    out = pd.DataFrame({
        "ID1":      df[id1_col],
        "ID2":      df[id2_col],
        "KINSHIP":  df[kin_col].astype(float),
    })
    out["relationship"]  = out["KINSHIP"].map(classify)
    out["above_cutoff"]  = out["KINSHIP"] >= cutoff
    out = out.sort_values("KINSHIP", ascending=False)
    out.to_csv(tsv_out, sep="\t", index=False)

    flagged = out[out["above_cutoff"]]
    json_out.write_text(json.dumps({
        "cutoff": cutoff,
        "n_pairs_above_cutoff": int(len(flagged)),
        "pairs": flagged.head(50).to_dict(orient="records"),
    }, indent=2, default=float) + "\n")


if __name__ == "__main__":
    main()
