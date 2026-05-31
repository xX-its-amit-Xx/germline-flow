"""
make_report.py — build the cohort QC HTML report.

Reads the tidy TSV / JSON produced by cohort_qc.py and the PLINK2 PCA + KING
outputs, renders three plotly figures, and embeds them in a Jinja2 template.
Each figure ships a one-paragraph "what + why" caption so a student reading
the report knows what they're looking at without leaving the page.

Figures:
  1. Per-sample coverage distribution (bar) with the QC threshold drawn in.
  2. Per-sample Ti/Tv ratio (scatter) with the expected WGS / WES bands.
  3. PCA scatter (PC1 vs PC2) coloured by inferred sex — first-pass ancestry /
     outlier detection.
  4. Relatedness summary panel (table) from PLINK2 KING.

The template is inlined to avoid a fourth file under workflow/scripts/.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from jinja2 import Template
from plotly.io import to_html

snakemake: Any  # type: ignore[no-redef]


# ---- captions (the teaching tone) -------------------------------------------

CAPTIONS = {
    "coverage": (
        "<b>What this shows.</b> Mean sequencing depth for each sample across the "
        "calling intervals, computed by mosdepth on the post-BQSR BAM. The dashed "
        "line is the coverage threshold from config (<code>qc.coverage_threshold</code>); "
        "any bar below it is flagged in the per-sample QC table."
        "<br><b>Why it matters.</b> Germline variant calls become unreliable below "
        "~15-20x because GATK can't distinguish a heterozygous site from a stretch "
        "of reference reads with one sequencing error. Low-coverage samples produce "
        "false negatives at heterozygous sites and inflated singleton counts."
    ),
    "titv": (
        "<b>What this shows.</b> Transition/transversion (Ti/Tv) ratio per sample, "
        "from the per-sample slice of <code>bcftools stats</code>. The shaded band "
        "is the expected range for the assay type (set in <code>qc.titv_min/max</code>)."
        "<br><b>Why it matters.</b> Real germline SNVs are enriched for transitions "
        "(A&lt;-&gt;G, C&lt;-&gt;T) at ~2.0-2.1 for WGS or ~3.0-3.3 for WES (CpG-rich "
        "coding regions). A Ti/Tv too close to 0.5 means the callset is dominated by "
        "random sequencing errors. A Ti/Tv far above the expected band suggests "
        "aggressive filtering has stripped the harder-to-call transversions."
    ),
    "pca": (
        "<b>What this shows.</b> PCA on the QC'd PLINK2 genotypes "
        "(<code>plink2 --pca</code>). Each point is one sample; PC1 vs PC2 captures "
        "the largest two axes of genetic variance in the cohort. Points are coloured "
        "by inferred genetic sex."
        "<br><b>Why it matters.</b> PCs 1-3 typically separate continental ancestries "
        "in human cohorts (1000 Genomes-style). Outliers along these axes are usually "
        "either samples from a different population than the rest (which you'll want "
        "to control for in GWAS) or technical artefacts (contamination, library "
        "prep batch effects)."
    ),
    "relatedness": (
        "<b>What this shows.</b> KING-robust kinship (Manichaikul et al., 2010) for "
        "every sample pair, computed by <code>plink2 --make-king-table</code>. The "
        "table lists pairs whose kinship coefficient exceeds "
        "<code>plink_export.king_cutoff</code> (default 0.0884 = 2nd-degree relatives)."
        "<br><b>Why it matters.</b> Unintended close relatives in a population study "
        "inflate type-I error in association testing and bias allele-frequency "
        "estimates. Duplicates or MZ twins (kinship > 0.354) are usually sample "
        "labelling errors that need to be resolved before downstream analysis."
    ),
}


# ---- figure builders --------------------------------------------------------

def fig_coverage(df: pd.DataFrame, threshold: float) -> go.Figure:
    fig = px.bar(
        df,
        x="sample",
        y="mean_coverage",
        color="pass_coverage",
        color_discrete_map={True: "#2a9d8f", False: "#e76f51"},
        labels={"mean_coverage": "Mean coverage (x)", "sample": "Sample"},
        title="Per-sample mean coverage",
    )
    fig.add_hline(
        y=threshold,
        line_dash="dash",
        line_color="black",
        annotation_text=f"threshold = {threshold}x",
        annotation_position="top left",
    )
    fig.update_layout(showlegend=False, height=400, margin=dict(l=40, r=20, t=40, b=40))
    return fig


def fig_titv(df: pd.DataFrame, titv_min: float, titv_max: float) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["sample"], y=df["titv"], mode="markers",
        marker=dict(size=12, color=df["pass_titv"].map({True: "#2a9d8f", False: "#e76f51"})),
        text=df["sample"], name="Ti/Tv",
    ))
    fig.add_hrect(
        y0=titv_min, y1=titv_max, fillcolor="#a8dadc", opacity=0.25, line_width=0,
        annotation_text=f"expected ({titv_min}-{titv_max})",
        annotation_position="top left",
    )
    fig.update_layout(
        title="Per-sample Ti/Tv ratio",
        xaxis_title="Sample",
        yaxis_title="Ti / Tv",
        height=400, margin=dict(l=40, r=20, t=40, b=40),
    )
    return fig


def fig_pca(eigenvec: pd.DataFrame, sex_map: dict[str, str]) -> go.Figure:
    eigenvec = eigenvec.copy()
    eigenvec["sex"] = eigenvec["IID"].map(sex_map).fillna("U")
    fig = px.scatter(
        eigenvec, x="PC1", y="PC2", color="sex", hover_name="IID",
        color_discrete_map={"F": "#e76f51", "M": "#2a9d8f", "U": "#999999"},
        title="PCA of cohort genotypes (PC1 vs PC2)",
    )
    fig.update_traces(marker=dict(size=10))
    fig.update_layout(height=450, margin=dict(l=40, r=20, t=40, b=40))
    return fig


def relatedness_table_html(related_json: dict) -> str:
    pairs = related_json.get("pairs", [])
    cutoff = related_json.get("cutoff")
    if not pairs:
        return (
            f"<p><i>No sample pairs exceed the KING cutoff "
            f"(&ge;&nbsp;{cutoff}). All pairs are inferred unrelated "
            "at the configured threshold.</i></p>"
        )
    rows = "".join(
        f"<tr><td>{p['ID1']}</td><td>{p['ID2']}</td>"
        f"<td>{p['KINSHIP']:.4f}</td><td>{p['relationship']}</td></tr>"
        for p in pairs
    )
    return (
        "<table class='related'>"
        "<thead><tr><th>ID1</th><th>ID2</th><th>KINSHIP</th><th>Inferred relationship</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
    )


# ---- Jinja2 template --------------------------------------------------------

TEMPLATE = Template(r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{{ title }}</title>
<style>
  :root { --fg:#222; --bg:#fafafa; --muted:#666; --accent:#264653; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
         color: var(--fg); background: var(--bg); margin: 0; padding: 2rem 3rem; max-width: 1200px; }
  h1 { color: var(--accent); border-bottom: 3px solid var(--accent); padding-bottom: .3rem; }
  h2 { color: var(--accent); margin-top: 2.5rem; }
  .meta { color: var(--muted); font-size: .9em; margin-bottom: 1.5rem; }
  .caption { background: #eef3f5; border-left: 4px solid var(--accent);
             padding: .75rem 1rem; margin: .5rem 0 1.5rem 0; font-size: .92em; line-height: 1.45; }
  table { border-collapse: collapse; margin: 1rem 0; font-size: .9em; }
  th, td { border: 1px solid #ddd; padding: .35rem .6rem; text-align: left; }
  th { background: #e9ecef; }
  table.summary td.pass { background: #d8f3dc; }
  table.summary td.fail { background: #ffd6d6; }
  table.related th { background: #f1faee; }
  .disclaimer { background: #fff3cd; border: 1px solid #ffe69c; padding: 1rem 1.25rem;
                margin: 2rem 0; font-size: .9em; border-radius: 4px; }
  code { background: #f1f1f1; padding: .1rem .3rem; border-radius: 3px; }
  footer { margin-top: 3rem; color: var(--muted); font-size: .8em; }
</style>
</head>
<body>

<h1>{{ title }}</h1>
<div class="meta">Generated by germline-flow &middot; {{ n_samples }} sample(s)
&middot; {{ n_pass }} pass / {{ n_fail }} fail overall QC</div>

<div class="disclaimer">
  <b>Research / education tooling.</b> This report is not validated for clinical use.
  Validate any actionable finding in an accredited laboratory before making clinical
  decisions. See the <a href="https://github.com/your-org/germline-flow#disclaimer">disclaimer</a> in the project README.
</div>

<h2>1. Per-sample QC summary</h2>
{{ summary_table | safe }}

<h2>2. Coverage</h2>
<div class="caption">{{ caption_coverage | safe }}</div>
{{ fig_coverage | safe }}

<h2>3. Ti/Tv ratio</h2>
<div class="caption">{{ caption_titv | safe }}</div>
{{ fig_titv | safe }}

<h2>4. Principal-components analysis</h2>
<div class="caption">{{ caption_pca | safe }}</div>
{{ fig_pca | safe }}

<h2>5. Relatedness (KING-robust kinship)</h2>
<div class="caption">{{ caption_relatedness | safe }}</div>
{{ relatedness_table | safe }}

<h2>6. Cohort variant-calling metrics (GATK)</h2>
<div class="caption">Per-cohort metrics from <code>gatk CollectVariantCallingMetrics</code>
against the dbSNP truth set. Use <code>DBSNP_TI_TV_RATIO</code> and
<code>NOVEL_TI_TV_RATIO</code> together: known-site Ti/Tv should be close to the
expected band; novel Ti/Tv well below it suggests the novel calls are enriched for
error.</div>
{{ cohort_metrics_table | safe }}

<footer>
  Built with Snakemake &middot; figures by plotly &middot; report styled with the
  germline-flow default Jinja2 template.<br>
  Cite the GATK Best Practices (Van der Auwera &amp; O'Connor, 2020) when using
  the output in published work.
</footer>

</body>
</html>
""")


def _summary_table_html(df: pd.DataFrame) -> str:
    cols = [
        "sample", "mean_coverage", "frac_bases_above_threshold",
        "n_snps", "n_indels", "titv", "het_hom_ratio",
        "inferred_sex", "reported_sex", "pass_overall",
    ]
    def _fmt(v: Any) -> str:
        if v is None or (isinstance(v, float) and pd.isna(v)):
            return "NA"
        if isinstance(v, float):
            return f"{v:.3f}"
        return str(v)
    head = "<tr>" + "".join(f"<th>{c}</th>" for c in cols) + "</tr>"
    body = []
    for _, row in df.iterrows():
        pass_cls = "pass" if row["pass_overall"] else "fail"
        cells = "".join(
            f'<td class="{pass_cls if c == "pass_overall" else ""}">{_fmt(row[c])}</td>'
            for c in cols
        )
        body.append(f"<tr>{cells}</tr>")
    return f'<table class="summary"><thead>{head}</thead><tbody>{"".join(body)}</tbody></table>'


def _cohort_metrics_table_html(metrics: dict[str, float]) -> str:
    if not metrics:
        return "<p><i>No cohort metrics produced.</i></p>"
    rows = "".join(f"<tr><td>{k}</td><td>{v:g}</td></tr>" for k, v in metrics.items())
    return (
        "<table><thead><tr><th>Metric</th><th>Value</th></tr></thead>"
        f"<tbody>{rows}</tbody></table>"
    )


def _read_eigenvec(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = [c.replace("#", "").strip() for c in df.columns]
    # PLINK2 emits IID as the second column; older versions used FID then IID.
    if "IID" not in df.columns:
        df = df.rename(columns={df.columns[1]: "IID"})
    return df


def main() -> None:
    cohort_tsv   = Path(snakemake.input.cohort_tsv)
    cohort_json  = Path(snakemake.input.cohort_json)
    eigenvec     = Path(snakemake.input.eigenvec)
    related_json = Path(snakemake.input.related_json)
    out_html     = Path(snakemake.output.html)

    qc_df       = pd.read_csv(cohort_tsv, sep="\t")
    cohort_meta = json.loads(cohort_json.read_text())
    related     = json.loads(related_json.read_text())
    pca_df      = _read_eigenvec(eigenvec)

    thresholds = cohort_meta["thresholds"]

    sex_map = dict(zip(qc_df["sample"], qc_df["inferred_sex"].fillna("U"), strict=True))

    f_cov  = fig_coverage(qc_df, thresholds["coverage_threshold"])
    f_titv = fig_titv(qc_df, thresholds["titv_min"], thresholds["titv_max"])
    f_pca  = fig_pca(pca_df, sex_map)

    rendered = TEMPLATE.render(
        title=str(snakemake.params.title),
        n_samples=len(qc_df),
        n_pass=cohort_meta["n_pass"],
        n_fail=cohort_meta["n_fail"],
        summary_table=_summary_table_html(qc_df),
        caption_coverage=CAPTIONS["coverage"],
        caption_titv=CAPTIONS["titv"],
        caption_pca=CAPTIONS["pca"],
        caption_relatedness=CAPTIONS["relatedness"],
        fig_coverage=to_html(f_cov, include_plotlyjs="cdn", full_html=False),
        fig_titv=to_html(f_titv, include_plotlyjs=False, full_html=False),
        fig_pca=to_html(f_pca, include_plotlyjs=False, full_html=False),
        relatedness_table=relatedness_table_html(related),
        cohort_metrics_table=_cohort_metrics_table_html(
            cohort_meta.get("cohort_variant_calling_metrics", {})
        ),
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(rendered)


if __name__ == "__main__":
    main()
