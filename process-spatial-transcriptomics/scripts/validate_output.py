#!/usr/bin/env python3
"""Validate the canonical output contract for a processed h5ad file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

from pipeline_utils import json_safe, matrix_stats, write_json


def validate_adata(adata: ad.AnnData) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []

    for key in ("sample_id", "batch", "total_counts", "n_genes_by_counts"):
        if key not in adata.obs:
            errors.append(f"Missing obs[{key!r}]")
    if "counts" not in adata.layers:
        errors.append("Missing layers['counts']")
    if "spatial" not in adata.obsm:
        errors.append("Missing obsm['spatial']")
    for key in ("X_pca", "X_umap"):
        if key not in adata.obsm:
            errors.append(f"Missing obsm[{key!r}]")
    for key in ("connectivities", "distances"):
        if key not in adata.obsp:
            errors.append(f"Missing expression graph obsp[{key!r}]")
    for key in ("spatial_connectivities", "spatial_distances"):
        if key not in adata.obsp:
            errors.append(f"Missing spatial graph obsp[{key!r}]")
    if "highly_variable" not in adata.var:
        errors.append("Missing var['highly_variable']")
    if "pipeline_provenance" not in adata.uns:
        errors.append("Missing uns['pipeline_provenance']")

    matrix_checks: dict[str, Any] = {}
    if adata.n_obs == 0 or adata.n_vars == 0:
        errors.append("Output has zero observations or variables")
    else:
        matrix_checks["X"] = matrix_stats(adata.X, adata.shape)
        if not matrix_checks["X"]["finite"]:
            errors.append("X contains non-finite sampled values")
        if matrix_checks["X"]["sample_min"] is not None and matrix_checks["X"]["sample_min"] < 0:
            warnings.append("X contains negative values; expected log-normalized expression")

        if "counts" in adata.layers:
            matrix_checks["counts"] = matrix_stats(
                adata.layers["counts"], adata.shape
            )
            if not matrix_checks["counts"]["integer_like"]:
                errors.append("layers['counts'] is not integer-like")
            if not matrix_checks["counts"]["finite"]:
                errors.append("layers['counts'] contains non-finite sampled values")

    graph_checks: dict[str, Any] = {}
    for key in ("connectivities", "spatial_connectivities"):
        if key not in adata.obsp:
            continue
        graph = adata.obsp[key]
        graph_checks[key] = {
            "shape": list(graph.shape),
            "nnz": int(graph.nnz) if sparse.issparse(graph) else int(np.count_nonzero(graph)),
        }
        if graph.shape != (adata.n_obs, adata.n_obs):
            errors.append(f"obsp[{key!r}] has invalid shape {graph.shape}")

    if "spatial_connectivities" in adata.obsp and "sample_id" in adata.obs:
        graph = sparse.coo_matrix(adata.obsp["spatial_connectivities"])
        samples = adata.obs["sample_id"].astype(str).to_numpy()
        cross = samples[graph.row] != samples[graph.col]
        n_cross = int(cross.sum())
        graph_checks["spatial_cross_sample_edges"] = n_cross
        if n_cross:
            errors.append(f"Spatial graph contains {n_cross} cross-sample edges")

    if "spatial" in adata.obsm:
        spatial = np.asarray(adata.obsm["spatial"])
        if spatial.shape[0] != adata.n_obs or spatial.ndim != 2:
            errors.append(f"obsm['spatial'] has invalid shape {spatial.shape}")
        if not np.isfinite(spatial).all():
            errors.append("obsm['spatial'] contains non-finite values")

    if "X_pca" in adata.obsm and np.asarray(adata.obsm["X_pca"]).shape[0] != adata.n_obs:
        errors.append("obsm['X_pca'] observation dimension does not match")
    if "X_umap" in adata.obsm and np.asarray(adata.obsm["X_umap"]).shape != (adata.n_obs, 2):
        errors.append("obsm['X_umap'] must have shape (n_obs, 2)")

    return {
        "status": "pass" if not errors else "fail",
        "shape": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
        "errors": errors,
        "warnings": warnings,
        "matrix_checks": matrix_checks,
        "graph_checks": graph_checks,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    result = validate_adata(adata)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        write_json(result, args.output)
    print(json.dumps(json_safe(result), indent=2, sort_keys=True))
    if result["status"] != "pass":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
