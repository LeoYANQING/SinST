# Output Contract

## AnnData Fields

Required output:

| Location | Key | Meaning |
|---|---|---|
| `X` | | library-size normalized, `log1p` full-gene expression |
| `layers` | `counts` | original integer-like counts after cell/gene filtering |
| `obs` | `sample_id` | tissue section or biological sample |
| `obs` | `batch` | technical batch/library label |
| `obs` | `time` | temporal value when supplied |
| `obs` | `total_counts` | count total before normalization |
| `obs` | `n_genes_by_counts` | detected genes before normalization |
| `var` | `highly_variable` | HVG mask used for PCA |
| `obsm` | `spatial` | source physical coordinates |
| `obsm` | `X_pca` | PCA coordinates |
| `obsm` | `X_umap` | UMAP coordinates |
| `obsp` | `connectivities`, `distances` | expression-space neighbor graph |
| `obsp` | `spatial_connectivities`, `spatial_distances` | physical graph |
| `uns` | `pipeline_provenance` | source, parameters, versions, warnings |
| `uns` | `pipeline_config` | fully resolved configuration |

`obs["leiden"]` is required only when Leiden dependencies are installed and the
configuration enables clustering.

## Run Artifacts

Each output directory contains:

- `processed.h5ad` or the configured filename;
- `input_audit.json`;
- `resolved_config.yaml`;
- `sample_qc_summary.csv`;
- `validation.json`;
- `run_summary.json`;
- `report.html`;
- `figures/qc_distributions.png`;
- `figures/umap_by_sample.png`;
- `figures/spatial_by_sample.png`.

## Validation Invariants

- counts are integer-like and finite;
- `X` is finite and nonnegative;
- embedding rows equal `n_obs`;
- graph matrices are square with shape `(n_obs, n_obs)`;
- the physical graph contains zero cross-sample edges;
- coordinates are finite and two- or three-dimensional;
- provenance is present.
