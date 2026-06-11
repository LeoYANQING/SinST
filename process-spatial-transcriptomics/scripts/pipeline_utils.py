#!/usr/bin/env python3
"""Shared utilities for the spatial transcriptomics preprocessing skill."""

from __future__ import annotations

import copy
import json
import os
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("NUMBA_CACHE_DIR", f"/tmp/st-pipeline-numba-{os.getuid()}")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/st-pipeline-matplotlib-{os.getuid()}")
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.spatial import Delaunay, cKDTree


SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
DEFAULT_CONFIG = SKILL_DIR / "assets" / "default-config.yaml"


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def load_config(path: str | Path | None, profile: str | None) -> dict[str, Any]:
    config_path = Path(path) if path else DEFAULT_CONFIG
    with config_path.open() as handle:
        raw = yaml.safe_load(handle) or {}

    profiles = raw.pop("profiles", {})
    selected = profile or raw.get("profile", "generic")
    if selected not in profiles:
        raise ValueError(
            f"Unknown profile {selected!r}; choose one of {sorted(profiles)}"
        )
    resolved = deep_merge(raw, profiles[selected] or {})
    resolved["profile"] = selected
    resolved["_config_source"] = str(config_path.resolve())
    return resolved


def write_yaml(data: dict[str, Any], path: str | Path) -> None:
    with Path(path).open("w") as handle:
        yaml.safe_dump(json_safe(data), handle, sort_keys=False)


def write_json(data: Any, path: str | Path) -> None:
    with Path(path).open("w") as handle:
        json.dump(json_safe(data), handle, indent=2, sort_keys=True)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [json_safe(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, pd.Timestamp):
        return value.isoformat()
    if pd.isna(value) if np.isscalar(value) and not isinstance(value, str) else False:
        return None
    return value


def package_versions(names: Iterable[str]) -> dict[str, str | None]:
    result: dict[str, str | None] = {}
    for name in names:
        try:
            result[name] = version(name)
        except PackageNotFoundError:
            result[name] = None
    return result


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def choose_existing(
    columns: Iterable[str], requested: str | None, candidates: Iterable[str]
) -> str | None:
    columns = set(columns)
    if requested and requested != "auto":
        return requested if requested in columns else None
    return next((candidate for candidate in candidates if candidate in columns), None)


def _dense_block(matrix: Any, rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    try:
        block = matrix[rows, :][:, cols]
    except (TypeError, IndexError):
        block = matrix[np.ix_(rows, cols)]
    if sparse.issparse(block):
        block = block.toarray()
    return np.asarray(block)


def sample_matrix(
    matrix: Any,
    shape: tuple[int, int],
    *,
    max_rows: int = 256,
    max_cols: int = 1024,
    seed: int = 0,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n_rows, n_cols = shape
    rows = np.sort(rng.choice(n_rows, min(max_rows, n_rows), replace=False))
    cols = np.sort(rng.choice(n_cols, min(max_cols, n_cols), replace=False))
    return _dense_block(matrix, rows, cols)


def matrix_stats(matrix: Any, shape: tuple[int, int], seed: int = 0) -> dict[str, Any]:
    values = sample_matrix(matrix, shape, seed=seed)
    finite = values[np.isfinite(values)]
    nonzero = finite[finite != 0]
    integer_like = bool(
        nonzero.size == 0
        or np.allclose(nonzero, np.round(nonzero), rtol=0, atol=1e-6)
    )
    return {
        "dtype": str(values.dtype),
        "sample_size": int(values.size),
        "sample_nonzero": int(nonzero.size),
        "sample_min": float(finite.min()) if finite.size else None,
        "sample_max": float(finite.max()) if finite.size else None,
        "sample_min_nonzero": float(nonzero.min()) if nonzero.size else None,
        "integer_like": integer_like,
        "finite": bool(finite.size == values.size),
    }


def compare_log1p(
    expression: Any,
    counts: Any,
    shape: tuple[int, int],
    seed: int = 0,
) -> bool:
    rng = np.random.default_rng(seed)
    n_rows, n_cols = shape
    rows = np.sort(rng.choice(n_rows, min(128, n_rows), replace=False))
    cols = np.sort(rng.choice(n_cols, min(512, n_cols), replace=False))
    x = _dense_block(expression, rows, cols)
    c = _dense_block(counts, rows, cols)
    return bool(np.allclose(np.expm1(x), c, rtol=1e-5, atol=1e-5))


def resolve_keys(adata: ad.AnnData, config: dict[str, Any]) -> dict[str, str | None]:
    canonical = config["canonical"]
    return {
        "sample_source": choose_existing(
            adata.obs.columns,
            canonical.get("sample_key"),
            ("sample_id", "library_id", "Batch", "batch"),
        ),
        "batch_source": choose_existing(
            adata.obs.columns,
            canonical.get("batch_key"),
            ("batch", "Batch", "library_id", "sample_id"),
        ),
        "time_source": choose_existing(
            adata.obs.columns,
            canonical.get("time_key"),
            ("time", "Time", "stage", "Stage", "dpi"),
        ),
        "gene_symbol_source": choose_existing(
            adata.var.columns,
            canonical.get("gene_symbol_column"),
            ("gene_symbol", "gene_symbols", "symbol", "Gene"),
        ),
        "spatial_source": choose_existing(
            adata.obsm.keys(),
            canonical.get("spatial_key"),
            ("spatial", "X_spatial"),
        ),
    }


def resolve_counts_layer(
    adata: ad.AnnData, config: dict[str, Any]
) -> tuple[str | None, Any | None]:
    requested = config["expression"].get("counts_layer", "counts")
    candidates = [requested, "counts", "raw_counts"]
    seen: set[str] = set()
    for key in candidates:
        if not key or key in seen:
            continue
        seen.add(key)
        if key in adata.layers:
            stats = matrix_stats(adata.layers[key], adata.shape)
            if stats["integer_like"] and stats["finite"]:
                return key, adata.layers[key]
    x_stats = matrix_stats(adata.X, adata.shape)
    if x_stats["integer_like"] and x_stats["finite"]:
        return "X", adata.X
    return None, None


def confounding_audit(
    obs: pd.DataFrame, batch_key: str | None, time_key: str | None
) -> dict[str, Any]:
    result = {
        "available": False,
        "batch_nested_in_time": None,
        "time_nested_in_batch": None,
        "warning": None,
    }
    if not batch_key or not time_key or batch_key not in obs or time_key not in obs:
        return result

    frame = obs[[batch_key, time_key]].dropna().astype(str).drop_duplicates()
    if frame.empty:
        return result
    batch_to_time = frame.groupby(batch_key, observed=True)[time_key].nunique()
    time_to_batch = frame.groupby(time_key, observed=True)[batch_key].nunique()
    result.update(
        {
            "available": True,
            "n_batches": int(frame[batch_key].nunique()),
            "n_times": int(frame[time_key].nunique()),
            "batch_nested_in_time": bool((batch_to_time == 1).all()),
            "time_nested_in_batch": bool((time_to_batch == 1).all()),
            "max_times_per_batch": int(batch_to_time.max()),
            "max_batches_per_time": int(time_to_batch.max()),
        }
    )
    if (
        result["batch_nested_in_time"]
        and result["n_batches"] > 1
        and result["n_times"] > 1
    ):
        result["warning"] = (
            "Each batch occurs at only one time point. Batch and temporal effects "
            "are not identifiable without additional experimental replication."
        )
    return result


def derive_time_values(
    obs: pd.DataFrame,
    config: dict[str, Any],
    *,
    sample_values: pd.Series | None = None,
) -> tuple[pd.Series | None, dict[str, Any] | None]:
    canonical = config["canonical"]
    mapping = canonical.get("time_mapping") or {}
    if mapping:
        source = sample_values
        if source is None:
            if "sample_id" not in obs:
                return None, None
            source = obs["sample_id"].astype(str)
        mapped = source.astype(str).map(
            {str(key): value for key, value in mapping.items()}
        )
        return mapped, {
            "method": "mapping",
            "mapped_values": int(mapped.notna().sum()),
        }

    pattern = canonical.get("time_regex")
    if not pattern:
        return None, None
    source_key = canonical.get("time_regex_source", "sample_id")
    if source_key == "sample_id" and sample_values is not None:
        source = sample_values.astype(str)
    elif source_key in obs:
        source = obs[source_key].astype(str)
    else:
        return None, {
            "method": "regex",
            "pattern": pattern,
            "source": source_key,
            "error": "source column unavailable",
        }
    extracted = source.str.extract(pattern, expand=False)
    if isinstance(extracted, pd.DataFrame):
        extracted = extracted.iloc[:, 0]
    numeric = pd.to_numeric(extracted, errors="coerce")
    return numeric, {
        "method": "regex",
        "pattern": pattern,
        "source": source_key,
        "mapped_values": int(numeric.notna().sum()),
    }


def audit_adata(
    adata: ad.AnnData, config: dict[str, Any], source: str | Path
) -> dict[str, Any]:
    keys = resolve_keys(adata, config)
    counts_key, counts = resolve_counts_layer(adata, config)
    x_stats = matrix_stats(adata.X, adata.shape)
    counts_stats = matrix_stats(counts, adata.shape) if counts is not None else None

    if counts is None:
        expression_state = "unknown_no_integer_counts"
    elif counts_key == "X":
        expression_state = "raw_counts_in_X"
    elif compare_log1p(adata.X, counts, adata.shape):
        expression_state = "X_equals_log1p_counts"
    else:
        expression_state = "counts_layer_available_X_other"

    source_path = Path(source)
    stat = source_path.stat()
    time_derivation = None
    if keys["time_source"]:
        confounding = confounding_audit(
            adata.obs, keys["batch_source"], keys["time_source"]
        )
    else:
        sample_values = (
            adata.obs[keys["sample_source"]].astype(str)
            if keys["sample_source"]
            else pd.Series("sample_1", index=adata.obs_names)
        )
        derived_time, time_derivation = derive_time_values(
            adata.obs, config, sample_values=sample_values
        )
        if derived_time is not None and keys["batch_source"]:
            derived_frame = pd.DataFrame(
                {
                    "_batch": adata.obs[keys["batch_source"]].astype(str).to_numpy(),
                    "_time": derived_time.to_numpy(),
                },
                index=adata.obs_names,
            )
            confounding = confounding_audit(
                derived_frame, "_batch", "_time"
            )
        else:
            confounding = confounding_audit(adata.obs, None, None)

    audit = {
        "source": str(source_path.resolve()),
        "source_size_bytes": stat.st_size,
        "source_mtime_utc": datetime.fromtimestamp(
            stat.st_mtime, tz=timezone.utc
        ).isoformat(),
        "profile": config["profile"],
        "shape": {"n_obs": int(adata.n_obs), "n_vars": int(adata.n_vars)},
        "duplicate_obs_names": int(adata.obs_names.duplicated().sum()),
        "duplicate_var_names": int(adata.var_names.duplicated().sum()),
        "keys": keys,
        "obs_columns": list(map(str, adata.obs.columns)),
        "var_columns": list(map(str, adata.var.columns)),
        "layers": list(map(str, adata.layers.keys())),
        "obsm": list(map(str, adata.obsm.keys())),
        "obsp": list(map(str, adata.obsp.keys())),
        "uns": list(map(str, adata.uns.keys())),
        "counts_source": counts_key,
        "expression_state": expression_state,
        "X_stats": x_stats,
        "counts_stats": counts_stats,
        "time_derivation": time_derivation,
        "confounding": confounding,
        "audited_at_utc": utc_now(),
    }
    if keys["spatial_source"]:
        spatial = np.asarray(adata.obsm[keys["spatial_source"]])
        audit["spatial"] = {
            "key": keys["spatial_source"],
            "shape": list(spatial.shape),
            "finite": bool(np.isfinite(spatial).all()),
        }
    else:
        audit["spatial"] = None
    return audit


def canonicalize_adata(
    adata: ad.AnnData, config: dict[str, Any], audit: dict[str, Any]
) -> dict[str, str | None]:
    keys = audit["keys"]
    sample_source = keys["sample_source"]
    if sample_source:
        adata.obs["sample_id"] = adata.obs[sample_source].astype(str).to_numpy()
    else:
        adata.obs["sample_id"] = "sample_1"

    batch_source = keys["batch_source"]
    if batch_source:
        adata.obs["batch"] = adata.obs[batch_source].astype(str).to_numpy()
    else:
        adata.obs["batch"] = adata.obs["sample_id"].astype(str).to_numpy()

    time_source = keys["time_source"]
    if time_source:
        adata.obs["time"] = adata.obs[time_source].to_numpy()
    else:
        derived_time, _ = derive_time_values(
            adata.obs,
            config,
            sample_values=adata.obs["sample_id"].astype(str),
        )
        if derived_time is not None and derived_time.notna().any():
            adata.obs["time"] = derived_time.to_numpy()

    gene_source = keys["gene_symbol_source"]
    if gene_source:
        adata.var["gene_symbol"] = adata.var[gene_source].astype(str).to_numpy()

    spatial_source = keys["spatial_source"]
    if spatial_source and spatial_source != "spatial":
        adata.obsm["spatial"] = np.asarray(adata.obsm[spatial_source]).copy()

    if adata.obs_names.has_duplicates:
        original = adata.obs_names.astype(str)
        adata.obs["source_obs_name"] = original
        adata.obs_names = pd.Index(
            adata.obs["sample_id"].astype(str) + "::" + original
        )
        adata.obs_names_make_unique()
    if adata.var_names.has_duplicates:
        adata.var["source_var_name"] = adata.var_names.astype(str)
        adata.var_names_make_unique()

    return {
        "sample_key": "sample_id",
        "batch_key": "batch",
        "time_key": "time" if "time" in adata.obs else None,
        "spatial_key": "spatial" if "spatial" in adata.obsm else None,
        "gene_symbol_column": "gene_symbol" if "gene_symbol" in adata.var else None,
    }


def mito_mask(adata: ad.AnnData, config: dict[str, Any]) -> np.ndarray:
    qc = config["qc"]
    sources = [pd.Series(adata.var_names.astype(str), index=adata.var_names)]
    for column in ("gene_symbol", "Gene", "symbol", "gene_symbols"):
        if column in adata.var:
            sources.append(adata.var[column].astype(str))
    patterns = tuple(str(x).upper() for x in qc.get("mitochondrial_prefixes", []))
    exact = {str(x).upper() for x in qc.get("mitochondrial_genes", [])}
    mask = np.zeros(adata.n_vars, dtype=bool)
    for source in sources:
        values = source.astype(str).str.upper()
        if patterns:
            mask |= values.str.startswith(patterns).to_numpy()
        if exact:
            mask |= values.isin(exact).to_numpy()
    return mask


def calculate_qc(adata: ad.AnnData, config: dict[str, Any]) -> dict[str, Any]:
    counts = adata.layers["counts"]
    n_counts = np.asarray(counts.sum(axis=1)).ravel()
    n_genes = np.asarray((counts > 0).sum(axis=1)).ravel()
    adata.obs["total_counts"] = n_counts.astype(np.float32)
    adata.obs["n_genes_by_counts"] = n_genes.astype(np.int32)

    mt = mito_mask(adata, config)
    adata.var["mt"] = mt
    if mt.any():
        mt_counts = np.asarray(counts[:, mt].sum(axis=1)).ravel()
        adata.obs["pct_counts_mt"] = (
            100.0 * mt_counts / np.maximum(n_counts, 1.0)
        ).astype(np.float32)
    else:
        adata.obs["pct_counts_mt"] = np.nan
    return {
        "mitochondrial_genes_detected": int(mt.sum()),
        "median_total_counts": float(np.median(n_counts)),
        "median_genes_by_counts": float(np.median(n_genes)),
    }


def filter_qc(
    adata: ad.AnnData, config: dict[str, Any]
) -> tuple[ad.AnnData, dict[str, Any]]:
    qc = config["qc"]
    keep = adata.obs["n_genes_by_counts"].to_numpy() >= int(qc["min_genes"])
    min_counts = qc.get("min_counts")
    if min_counts is not None:
        keep &= adata.obs["total_counts"].to_numpy() >= float(min_counts)
    max_genes = qc.get("max_genes")
    if max_genes is not None:
        keep &= adata.obs["n_genes_by_counts"].to_numpy() <= float(max_genes)
    max_counts = qc.get("max_counts")
    if max_counts is not None:
        keep &= adata.obs["total_counts"].to_numpy() <= float(max_counts)
    max_pct_mt = qc.get("max_pct_mt")
    if max_pct_mt is not None and adata.var["mt"].any():
        keep &= adata.obs["pct_counts_mt"].to_numpy() <= float(max_pct_mt)

    before_cells = adata.n_obs
    filtered = adata[keep].copy()
    gene_ncells = np.asarray((filtered.layers["counts"] > 0).sum(axis=0)).ravel()
    gene_keep = gene_ncells >= int(qc["min_cells"])
    before_genes = filtered.n_vars
    filtered = filtered[:, gene_keep].copy()
    return filtered, {
        "cells_before": int(before_cells),
        "cells_after": int(filtered.n_obs),
        "cells_removed": int(before_cells - filtered.n_obs),
        "genes_before": int(before_genes),
        "genes_after": int(filtered.n_vars),
        "genes_removed": int(before_genes - filtered.n_vars),
    }


def _knn_edges(coords: np.ndarray, n_neighbors: int) -> tuple[np.ndarray, np.ndarray]:
    k = min(max(1, n_neighbors), max(1, len(coords) - 1))
    distances, indices = cKDTree(coords).query(coords, k=k + 1)
    rows = np.repeat(np.arange(len(coords)), k)
    cols = indices[:, 1:].reshape(-1)
    dists = distances[:, 1:].reshape(-1)
    return np.column_stack([rows, cols]), dists


def _delaunay_edges(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    triangulation = Delaunay(coords)
    simplices = triangulation.simplices
    edge_set: set[tuple[int, int]] = set()
    for simplex in simplices:
        for i in range(len(simplex)):
            for j in range(i + 1, len(simplex)):
                a, b = sorted((int(simplex[i]), int(simplex[j])))
                edge_set.add((a, b))
    edges = np.asarray(sorted(edge_set), dtype=int)
    distances = np.linalg.norm(coords[edges[:, 0]] - coords[edges[:, 1]], axis=1)
    return edges, distances


def build_spatial_graph(
    adata: ad.AnnData, config: dict[str, Any]
) -> dict[str, Any]:
    spatial_cfg = config["spatial"]
    if "spatial" not in adata.obsm:
        if spatial_cfg.get("required", True):
            raise KeyError("No spatial coordinates found in adata.obsm")
        return {"built": False, "reason": "no_spatial_coordinates"}

    coords_all = np.asarray(adata.obsm["spatial"], dtype=float)
    if coords_all.ndim != 2 or coords_all.shape[1] not in (2, 3):
        raise ValueError("obsm['spatial'] must have shape (n_obs, 2) or (n_obs, 3)")
    if not np.isfinite(coords_all).all():
        raise ValueError("obsm['spatial'] contains non-finite coordinates")

    samples = adata.obs["sample_id"].astype(str).to_numpy()
    all_rows: list[np.ndarray] = []
    all_cols: list[np.ndarray] = []
    all_distances: list[np.ndarray] = []
    per_sample: dict[str, Any] = {}

    for sample in pd.unique(samples):
        global_index = np.flatnonzero(samples == sample)
        coords = coords_all[global_index]
        if len(coords) < 2:
            per_sample[sample] = {"n_obs": int(len(coords)), "n_edges": 0}
            continue

        method = spatial_cfg.get("method", "delaunay")
        used = method
        try:
            if method == "delaunay" and len(coords) >= coords.shape[1] + 2:
                edges, distances = _delaunay_edges(coords)
            else:
                used = "knn"
                edges, distances = _knn_edges(
                    coords, int(spatial_cfg.get("n_neighbors", 6))
                )
        except Exception:
            used = "knn"
            edges, distances = _knn_edges(
                coords, int(spatial_cfg.get("n_neighbors", 6))
            )

        prune_quantile = spatial_cfg.get("distance_prune_quantile")
        if prune_quantile is not None and len(distances):
            threshold = float(np.quantile(distances, float(prune_quantile)))
            keep = distances <= threshold
            edges, distances = edges[keep], distances[keep]
        else:
            threshold = None

        rows = global_index[edges[:, 0]]
        cols = global_index[edges[:, 1]]
        all_rows.extend([rows, cols])
        all_cols.extend([cols, rows])
        all_distances.extend([distances, distances])
        per_sample[sample] = {
            "n_obs": int(len(coords)),
            "n_undirected_edges": int(len(edges)),
            "method": used,
            "distance_prune_threshold": threshold,
        }

    if all_rows:
        rows = np.concatenate(all_rows)
        cols = np.concatenate(all_cols)
        distances = np.concatenate(all_distances)
    else:
        rows = cols = np.asarray([], dtype=int)
        distances = np.asarray([], dtype=float)

    shape = (adata.n_obs, adata.n_obs)
    connectivity = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)), shape=shape
    )
    distance = sparse.csr_matrix(
        (distances.astype(np.float32), (rows, cols)), shape=shape
    )
    connectivity.eliminate_zeros()
    distance.eliminate_zeros()
    adata.obsp["spatial_connectivities"] = connectivity
    adata.obsp["spatial_distances"] = distance
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "distances_key": "spatial_distances",
        "params": {
            "method": spatial_cfg.get("method", "delaunay"),
            "n_neighbors": int(spatial_cfg.get("n_neighbors", 6)),
            "per_sample": True,
        },
    }
    return {
        "built": True,
        "directed_nonzero_entries": int(connectivity.nnz),
        "per_sample": per_sample,
    }


def sample_summary(adata: ad.AnnData) -> pd.DataFrame:
    columns = ["sample_id", "batch"]
    if "time" in adata.obs:
        columns.append("time")
    metric_columns = ["total_counts", "n_genes_by_counts"]
    if "pct_counts_mt" in adata.obs:
        metric_columns.append("pct_counts_mt")
    frame = adata.obs[columns + metric_columns].copy()
    grouped = frame.groupby(columns, observed=True, dropna=False)
    summary = grouped.agg(
        n_obs=("total_counts", "size"),
        median_total_counts=("total_counts", "median"),
        median_genes_by_counts=("n_genes_by_counts", "median"),
    ).reset_index()
    if "pct_counts_mt" in adata.obs and adata.obs["pct_counts_mt"].notna().any():
        mt = grouped["pct_counts_mt"].median().rename("median_pct_counts_mt")
        summary = summary.merge(mt.reset_index(), on=columns, how="left")
    return summary
