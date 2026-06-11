# SinST

SinST is an audit-first spatial transcriptomics preprocessing pipeline and
Codex skill for AnnData datasets.

It inspects the input state before preprocessing, preserves raw counts, produces
standard PCA/UMAP outputs, constructs physical graphs separately for each tissue
section, detects batch-time confounding, and validates the resulting `.h5ad`.

## Features

- Detect raw counts, log-transformed matrices, and ambiguous input states.
- Canonicalize sample, batch, time, gene symbol, and spatial-coordinate fields.
- Calculate QC metrics and apply configurable cell/gene filters.
- Run normalization, `log1p`, HVG selection, PCA, neighbors, UMAP, and Leiden.
- Keep all filtered genes in the final AnnData while using HVGs for PCA.
- Construct Delaunay or kNN spatial graphs within each sample.
- Detect non-identifiable batch-time experimental designs.
- Generate a processed AnnData, validation results, figures, and an HTML report.
- Provide profiles for generic AnnData, axolotl Stereo-seq, and Visium data.

## Requirements

Python 3.10 or later is recommended.

Create an isolated environment and install the dependencies:

```bash
conda create -n sinst python=3.10 -y
conda activate sinst

pip install \
  anndata \
  scanpy \
  numpy \
  pandas \
  scipy \
  matplotlib \
  pyyaml \
  scikit-learn \
  umap-learn \
  python-igraph \
  leidenalg \
  h5py
```

Clone the repository:

```bash
git clone https://github.com/LeoYANQING/SinST.git
cd SinST
```

All commands below are run from the repository root.

## Quick Start

### 1. Inspect the input

Always audit a dataset before preprocessing:

```bash
python process-spatial-transcriptomics/scripts/inspect_dataset.py \
  /path/to/input.h5ad \
  --profile generic \
  --output results/input_audit.json
```

Check the audit:

```bash
cat results/input_audit.json
```

Important fields include:

- `counts_source`: location of the integer-like raw counts.
- `expression_state`: detected state of `X`.
- `duplicate_obs_names`: number of duplicate observation names.
- `spatial`: coordinate key, shape, and finite-value check.
- `confounding.warning`: batch-time identifiability warning.

Do not continue if `counts_source` is `null`. Configure the correct raw-count
layer instead of treating normalized expression as counts.

### 2. Run the complete pipeline

```bash
python process-spatial-transcriptomics/scripts/run_pipeline.py \
  /path/to/input.h5ad \
  --profile generic \
  --output-dir results
```

The pipeline executes:

```text
input audit
  -> metadata canonicalization
  -> QC and filtering
  -> normalize_total
  -> log1p
  -> HVG selection
  -> PCA
  -> expression neighbors
  -> UMAP and Leiden
  -> per-sample spatial graph
  -> validation and report
```

The source `.h5ad` is never overwritten.

### 3. Check the output

```bash
cat results/validation.json
```

The output is ready for downstream analysis when it contains:

```json
{
  "status": "pass"
}
```

Open `results/report.html` to inspect the QC, UMAP, spatial overview, provenance,
and confounding warnings.

## Profiles

Available profiles:

| Profile | Intended input |
|---|---|
| `generic` | AnnData with raw counts and spatial coordinates |
| `axolotl-stereoseq` | ARTISTA-style axolotl Stereo-seq data |
| `visium` | Visium AnnData with `obsm["spatial"]` |

Select a profile with:

```bash
--profile generic
--profile axolotl-stereoseq
--profile visium
```

## ARTISTA Example

Audit an ARTISTA regeneration dataset:

```bash
python process-spatial-transcriptomics/scripts/inspect_dataset.py \
  /home/data/ARTISTA/Regeneration.h5ad \
  --profile axolotl-stereoseq \
  --output results/artista/input_audit.json
```

Run the complete pipeline:

```bash
python process-spatial-transcriptomics/scripts/run_pipeline.py \
  /home/data/ARTISTA/Regeneration.h5ad \
  --profile axolotl-stereoseq \
  --output-dir results/artista
```

The axolotl profile:

- reads sample identity from `obs["Batch"]`;
- copies display annotation from `var["Gene"]`;
- derives numeric DPI values from ARTISTA sample names;
- disables mitochondrial filtering unless a trusted gene mapping is supplied;
- constructs a Delaunay graph independently for each section;
- reports batch-time confounding instead of silently applying integration.

Large ARTISTA files require substantial memory and runtime. Run the audit first
and consider testing a single section before processing the combined dataset.

## Output Files

Each run creates:

```text
results/
├── processed.h5ad
├── input_audit.json
├── resolved_config.yaml
├── sample_qc_summary.csv
├── validation.json
├── run_summary.json
├── report.html
└── figures/
    ├── qc_distributions.png
    ├── umap_by_sample.png
    └── spatial_by_sample.png
```

The processed AnnData contains:

```text
X                                  normalized log1p full-gene expression
layers["counts"]                   filtered integer-like raw counts
obs["sample_id"], obs["batch"]
obs["time"]                        when time metadata is available
obs["total_counts"]
obs["n_genes_by_counts"]
obs["leiden"]                      when Leiden is enabled
var["highly_variable"]
obsm["spatial"]
obsm["X_pca"]
obsm["X_umap"]
obsp["connectivities"]             expression-space graph
obsp["spatial_connectivities"]     physical per-sample graph
uns["pipeline_provenance"]
uns["pipeline_config"]
```

## Custom Configuration

Copy the default configuration:

```bash
cp \
  process-spatial-transcriptomics/assets/default-config.yaml \
  sinst-config.yaml
```

Edit `sinst-config.yaml`, then run:

```bash
python process-spatial-transcriptomics/scripts/run_pipeline.py \
  /path/to/input.h5ad \
  --config sinst-config.yaml \
  --profile generic \
  --output-dir results/custom
```

Common settings:

```yaml
expression:
  target_sum: 10000
  hvg:
    n_top_genes: 2000

qc:
  min_genes: 100
  min_cells: 3
  max_pct_mt: 20

embedding:
  n_pcs: 50
  neighbors:
    n_neighbors: 15
    n_pcs: 30

spatial:
  method: delaunay
  distance_prune_quantile: 0.99
```

Use an explicit `canonical.time_mapping` when time cannot be reliably parsed
from an existing metadata column or sample name.

## Independent Validation

Validate a processed dataset again after downstream modifications:

```bash
python process-spatial-transcriptomics/scripts/validate_output.py \
  results/processed.h5ad \
  --output results/validation.json
```

Regenerate the report:

```bash
python process-spatial-transcriptomics/scripts/render_report.py \
  results/processed.h5ad \
  --output-dir results
```

## Use as a Codex Skill

The skill is located at:

```text
process-spatial-transcriptomics/
```

It can be copied or linked into the Codex skills directory:

```bash
mkdir -p "${CODEX_HOME:-$HOME/.codex}/skills"

ln -s "$(pwd)/process-spatial-transcriptomics" \
  "${CODEX_HOME:-$HOME/.codex}/skills/process-spatial-transcriptomics"
```

Example request:

```text
Use $process-spatial-transcriptomics to audit and preprocess
/path/to/input.h5ad with the axolotl-stereoseq profile.
```

## Interpretation Notes

- Expression and physical neighbors are different graphs.
- Spatial graphs must not connect separate tissue sections before registration.
- Batch correction is intentionally not automatic.
- If each batch occurs at only one time point, batch and time effects are not
  identifiable from the dataset alone.
- Preserve source coordinates; store display rotations or registered coordinates
  under additional `obsm` keys.

For detailed policies, see:

- [`SKILL.md`](process-spatial-transcriptomics/SKILL.md)
- [`input-contract.md`](process-spatial-transcriptomics/references/input-contract.md)
- [`output-contract.md`](process-spatial-transcriptomics/references/output-contract.md)
- [`qc-policy.md`](process-spatial-transcriptomics/references/qc-policy.md)
- [`multi-sample-policy.md`](process-spatial-transcriptomics/references/multi-sample-policy.md)
- [`species-axolotl.md`](process-spatial-transcriptomics/references/species-axolotl.md)

## Current Scope

- Input format: AnnData `.h5ad`.
- No automatic doublet detection.
- No automatic batch integration.
- No spatial autocorrelation or Moran's I stage yet.

## License

SinST is released under the [MIT License](LICENSE).
