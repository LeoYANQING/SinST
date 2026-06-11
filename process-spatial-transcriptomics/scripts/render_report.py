#!/usr/bin/env python3
"""Render a compact HTML and PNG QC report for a processed h5ad file."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import anndata as ad

from pipeline_utils import json_safe


def _category_colors(values: pd.Series) -> tuple[np.ndarray, list[tuple[str, Any]]]:
    categorical = values.astype(str).astype("category")
    codes = categorical.cat.codes.to_numpy()
    palette = plt.get_cmap("tab20")
    colors = palette(codes % 20)
    legend = [
        (str(category), palette(index % 20))
        for index, category in enumerate(categorical.cat.categories)
    ]
    return colors, legend


def _save_qc(adata: ad.AnnData, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5))
    axes[0].hist(adata.obs["total_counts"], bins=50, color="#4C78A8")
    axes[0].set_title("Total counts")
    axes[1].hist(adata.obs["n_genes_by_counts"], bins=50, color="#59A14F")
    axes[1].set_title("Detected genes")
    mt = adata.obs.get("pct_counts_mt")
    if mt is not None and mt.notna().any():
        axes[2].hist(mt.dropna(), bins=50, color="#E15759")
        axes[2].set_title("Mitochondrial %")
    else:
        axes[2].text(0.5, 0.5, "No mitochondrial annotation", ha="center", va="center")
        axes[2].set_axis_off()
    fig.tight_layout()
    fig.savefig(output, dpi=140)
    plt.close(fig)


def _save_umap(adata: ad.AnnData, output: Path) -> None:
    embedding = np.asarray(adata.obsm["X_umap"])
    colors, legend = _category_colors(adata.obs["sample_id"])
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    ax.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=colors,
        s=3,
        linewidths=0,
        rasterized=True,
    )
    ax.set_title("UMAP by sample")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    if len(legend) <= 20:
        handles = [
            plt.Line2D(
                [0], [0], marker="o", linestyle="", color=color, label=label
            )
            for label, color in legend
        ]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5), fontsize=7)
    fig.tight_layout()
    fig.savefig(output, dpi=140)
    plt.close(fig)


def _save_spatial(adata: ad.AnnData, output: Path, max_panels: int = 12) -> None:
    samples = list(pd.unique(adata.obs["sample_id"].astype(str)))[:max_panels]
    ncols = min(4, len(samples))
    nrows = int(np.ceil(len(samples) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows), squeeze=False)
    coords = np.asarray(adata.obsm["spatial"])
    sample_values = adata.obs["sample_id"].astype(str).to_numpy()
    for ax, sample in zip(axes.ravel(), samples):
        mask = sample_values == sample
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=2,
            linewidths=0,
            c="#4C78A8",
            rasterized=True,
        )
        ax.set_title(f"{sample} (n={int(mask.sum()):,})", fontsize=9)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_axis_off()
    for ax in axes.ravel()[len(samples):]:
        ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output, dpi=140)
    plt.close(fig)


def render_report(
    adata: ad.AnnData,
    output_dir: str | Path,
    *,
    audit: dict[str, Any] | None = None,
    validation: dict[str, Any] | None = None,
    report_name: str = "report.html",
) -> Path:
    output_dir = Path(output_dir)
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    _save_qc(adata, figures / "qc_distributions.png")
    if "X_umap" in adata.obsm:
        _save_umap(adata, figures / "umap_by_sample.png")
    if "spatial" in adata.obsm:
        _save_spatial(adata, figures / "spatial_by_sample.png")

    sample_table = (
        adata.obs.groupby("sample_id", observed=True)
        .agg(
            n_obs=("total_counts", "size"),
            median_counts=("total_counts", "median"),
            median_genes=("n_genes_by_counts", "median"),
        )
        .reset_index()
    )
    provenance = adata.uns.get("pipeline_provenance", {})
    warning = (
        (audit or {}).get("canonical_confounding", {}).get("warning")
        or (audit or {}).get("confounding", {}).get("warning")
        or provenance.get("confounding", {}).get("warning")
    )
    warning_html = (
        f'<div class="warning"><strong>Confounding warning:</strong> {html.escape(warning)}</div>'
        if warning
        else ""
    )
    validation_status = (validation or {}).get("status", "not run")
    document = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Spatial transcriptomics preprocessing report</title>
<style>
body {{ font-family: sans-serif; max-width: 1100px; margin: 2rem auto; color: #222; }}
table {{ border-collapse: collapse; width: 100%; }}
th, td {{ border: 1px solid #ddd; padding: 0.4rem; text-align: left; }}
th {{ background: #f3f5f7; }}
img {{ max-width: 100%; margin: 0.5rem 0 1.5rem; }}
pre {{ background: #f6f8fa; padding: 1rem; overflow-x: auto; }}
.warning {{ background: #fff3cd; border: 1px solid #e5c96c; padding: 1rem; }}
</style></head><body>
<h1>Spatial transcriptomics preprocessing report</h1>
<p><strong>Shape:</strong> {adata.n_obs:,} observations x {adata.n_vars:,} genes</p>
<p><strong>Validation:</strong> {html.escape(validation_status)}</p>
{warning_html}
<h2>Samples</h2>
{sample_table.to_html(index=False, border=0)}
<h2>QC distributions</h2>
<img src="figures/qc_distributions.png" alt="QC distributions">
<h2>UMAP</h2>
<img src="figures/umap_by_sample.png" alt="UMAP by sample">
<h2>Spatial overview</h2>
<img src="figures/spatial_by_sample.png" alt="Spatial overview">
<h2>Provenance</h2>
<pre>{html.escape(json.dumps(json_safe(provenance), indent=2, sort_keys=True))}</pre>
</body></html>"""
    output = output_dir / report_name
    output.write_text(document)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--audit", type=Path, default=None)
    parser.add_argument("--validation", type=Path, default=None)
    parser.add_argument("--report-name", default="report.html")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    audit = json.loads(args.audit.read_text()) if args.audit else None
    validation = json.loads(args.validation.read_text()) if args.validation else None
    output = render_report(
        adata,
        args.output_dir,
        audit=audit,
        validation=validation,
        report_name=args.report_name,
    )
    print(output)


if __name__ == "__main__":
    main()
