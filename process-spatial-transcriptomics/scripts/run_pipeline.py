#!/usr/bin/env python3
"""Run the standardized spatial transcriptomics preprocessing pipeline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

from pipeline_utils import (
    audit_adata,
    build_spatial_graph,
    calculate_qc,
    canonicalize_adata,
    confounding_audit,
    filter_qc,
    load_config,
    package_versions,
    sample_summary,
    utc_now,
    write_json,
    write_yaml,
)

import scanpy as sc

from render_report import render_report
from validate_output import validate_adata


def _copy_counts(adata: ad.AnnData, source: str) -> None:
    if source == "X":
        counts = adata.X.copy()
    else:
        counts = adata.layers[source].copy()
    if sparse.issparse(counts):
        counts = counts.tocsr()
    adata.layers["counts"] = counts


def _run_expression_pipeline(
    adata: ad.AnnData, config: dict[str, Any]
) -> dict[str, Any]:
    seed = int(config["runtime"]["random_seed"])
    expression = config["expression"]
    embedding = config["embedding"]

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=float(expression["target_sum"]))
    sc.pp.log1p(adata)

    hvg = expression["hvg"]
    n_top = min(int(hvg["n_top_genes"]), adata.n_vars)
    batch_key = hvg.get("batch_key")
    if batch_key not in adata.obs or adata.obs[batch_key].nunique() < 2:
        batch_key = None
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top,
        flavor=hvg.get("flavor", "seurat"),
        batch_key=batch_key,
        subset=False,
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    if n_hvg < 2:
        raise ValueError(f"Only {n_hvg} highly variable genes were selected")

    pca_input = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(
        pca_input,
        zero_center=False,
        max_value=float(embedding["scale_max_value"]),
    )
    max_pcs = min(
        int(embedding["n_pcs"]),
        pca_input.n_obs - 1,
        pca_input.n_vars - 1,
    )
    if max_pcs < 2:
        raise ValueError("Not enough observations or HVGs to calculate PCA")
    sc.tl.pca(
        pca_input,
        n_comps=max_pcs,
        zero_center=True,
        svd_solver="arpack",
        random_state=seed,
    )
    adata.obsm["X_pca"] = np.asarray(pca_input.obsm["X_pca"], dtype=np.float32)
    adata.uns["pca"] = pca_input.uns["pca"]
    adata.varm["PCs"] = np.zeros((adata.n_vars, max_pcs), dtype=np.float32)
    adata.varm["PCs"][adata.var["highly_variable"].to_numpy()] = np.asarray(
        pca_input.varm["PCs"], dtype=np.float32
    )

    neighbor_cfg = embedding["neighbors"]
    n_neighbors = min(int(neighbor_cfg["n_neighbors"]), adata.n_obs - 1)
    n_pcs = min(int(neighbor_cfg["n_pcs"]), max_pcs)
    sc.pp.neighbors(
        adata,
        n_neighbors=n_neighbors,
        n_pcs=n_pcs,
        use_rep="X_pca",
        random_state=seed,
    )
    sc.tl.umap(
        adata,
        min_dist=float(embedding["umap"]["min_dist"]),
        random_state=seed,
    )

    leiden_cfg = embedding["leiden"]
    leiden_status: dict[str, Any] = {"requested": bool(leiden_cfg["enabled"])}
    if leiden_cfg["enabled"]:
        try:
            sc.tl.leiden(
                adata,
                resolution=float(leiden_cfg["resolution"]),
                random_state=seed,
                key_added="leiden",
                flavor=leiden_cfg.get("flavor", "igraph"),
                n_iterations=int(leiden_cfg.get("n_iterations", 2)),
                directed=False,
            )
            leiden_status["completed"] = True
        except (ImportError, ModuleNotFoundError) as exc:
            leiden_status.update({"completed": False, "warning": str(exc)})
    return {
        "n_highly_variable_genes": n_hvg,
        "n_pcs": max_pcs,
        "n_neighbors": n_neighbors,
        "leiden": leiden_status,
    }


def run_pipeline(
    input_path: Path,
    output_dir: Path,
    config: dict[str, Any],
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    adata = ad.read_h5ad(input_path)
    audit = audit_adata(adata, config, input_path)
    write_json(audit, output_dir / "input_audit.json")
    if not audit["counts_source"]:
        raise ValueError(
            "No integer-like raw count matrix was found. Set expression.counts_layer "
            "to a valid raw-count layer before running preprocessing."
        )

    _copy_counts(adata, audit["counts_source"])
    canonical = canonicalize_adata(adata, config, audit)
    qc_before = calculate_qc(adata, config)
    adata, filtering = filter_qc(adata, config)
    if adata.n_obs < 3 or adata.n_vars < 3:
        raise ValueError(
            f"QC left too little data for embedding: shape={adata.shape}"
        )
    expression = _run_expression_pipeline(adata, config)
    spatial = build_spatial_graph(adata, config)
    confounding = confounding_audit(
        adata.obs, canonical["batch_key"], canonical["time_key"]
    )
    audit["canonical_confounding"] = confounding
    write_json(audit, output_dir / "input_audit.json")

    sample_summary(adata).to_csv(output_dir / "sample_qc_summary.csv", index=False)
    provenance = {
        "pipeline": "process-spatial-transcriptomics",
        "pipeline_schema_version": "1.0",
        "started_at_utc": started,
        "completed_at_utc": utc_now(),
        "source": audit["source"],
        "profile": config["profile"],
        "canonical_keys": canonical,
        "qc_before_filtering": qc_before,
        "filtering": filtering,
        "expression": expression,
        "spatial_graph": spatial,
        "confounding": confounding,
        "integration": config["integration"],
        "random_seed": int(config["runtime"]["random_seed"]),
        "package_versions": package_versions(
            (
                "anndata",
                "scanpy",
                "numpy",
                "pandas",
                "scipy",
                "scikit-learn",
                "umap-learn",
                "igraph",
                "leidenalg",
            )
        ),
    }
    adata.uns["pipeline_provenance"] = provenance
    adata.uns["pipeline_config"] = config

    output_h5ad = output_dir / config["output"]["h5ad_name"]
    write_yaml(config, output_dir / "resolved_config.yaml")
    adata.write_h5ad(output_h5ad, compression="gzip")
    validation = validate_adata(adata)
    write_json(validation, output_dir / "validation.json")
    report = render_report(
        adata,
        output_dir,
        audit=audit,
        validation=validation,
        report_name=config["output"]["report_name"],
    )
    result = {
        "output_h5ad": str(output_h5ad),
        "report": str(report),
        "validation": validation,
        "shape": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
    }
    write_json(result, output_dir / "run_summary.json")
    if validation["status"] != "pass":
        raise RuntimeError(
            "Pipeline completed but output validation failed: "
            + "; ".join(validation["errors"])
        )
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Input .h5ad file")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument(
        "--profile",
        choices=("generic", "axolotl-stereoseq", "visium"),
        default=None,
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config, args.profile)
    result = run_pipeline(args.input, args.output_dir, config)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
