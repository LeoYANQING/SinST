---
name: process-spatial-transcriptomics
description: Audit, standardize, preprocess, validate, and report spatial transcriptomics datasets stored as AnnData .h5ad files. Use when Codex receives a new spatial transcriptomics dataset and needs a reproducible pipeline covering raw-count detection, metadata canonicalization, QC, normalization, HVG selection, PCA, expression neighbors, UMAP, Leiden clustering, per-sample spatial graph construction, batch-time confounding checks, canonical AnnData export, or an axolotl Stereo-seq/Visium preprocessing profile.
---

# Process Spatial Transcriptomics

Convert a spatial transcriptomics `.h5ad` into a validated, reproducible AnnData
artifact. Treat inspection and provenance as required pipeline stages, not optional
notebook commentary.

## Workflow

1. Inspect the input before changing it.
2. Select a profile and review dataset-specific metadata.
3. Resolve raw counts and canonical keys.
4. Run preprocessing into a new output directory.
5. Inspect warnings, figures, and sample metrics.
6. Validate the output contract before using it downstream.

Never overwrite the source dataset. Never infer that a dataset is fully processed
only because PCA, UMAP, or annotation fields already exist.

## Inspect

Run:

```bash
python scripts/inspect_dataset.py INPUT.h5ad \
  --profile generic \
  --output OUTPUT/input_audit.json
```

Read [input-contract.md](references/input-contract.md) when counts, coordinate,
or metadata fields are ambiguous. Stop before preprocessing if no integer-like raw
counts can be identified.

Use `axolotl-stereoseq` for ARTISTA-style axolotl data and read
[species-axolotl.md](references/species-axolotl.md). Use `visium` for a Visium
AnnData with `obsm["spatial"]`; otherwise start with `generic`.

## Configure

Copy `assets/default-config.yaml` into the run output or project configuration
area only when defaults need adjustment. Pass the edited file with `--config`.

Review [qc-policy.md](references/qc-policy.md) before changing filtering thresholds.
Keep upper QC cutoffs explicit; do not silently remove a fixed top percentile.

For multi-sample or temporal data, read
[multi-sample-policy.md](references/multi-sample-policy.md). Leave
`integration.method: none` unless the experimental design makes batch correction
identifiable and the user explicitly selects a method.

## Run

```bash
python scripts/run_pipeline.py INPUT.h5ad \
  --profile axolotl-stereoseq \
  --output-dir OUTPUT
```

The runner:

- reconstructs `X` from confirmed raw counts;
- preserves integer counts in `layers["counts"]`;
- normalizes and log-transforms full-gene `X`;
- selects HVGs without permanently subsetting the output;
- scales only a temporary HVG matrix for PCA;
- builds expression neighbors, UMAP, and optional Leiden clusters;
- builds a separate spatial graph within each `sample_id`;
- records parameters, versions, source metadata, and confounding warnings;
- writes a processed `.h5ad`, QC tables, figures, HTML report, and validation JSON.

If Leiden dependencies are unavailable, preserve the successful embedding output
and report the skipped clustering step. Do not hide the warning.

## Validate

The main runner validates automatically. Re-run validation independently after
any downstream edit:

```bash
python scripts/validate_output.py OUTPUT/processed.h5ad \
  --output OUTPUT/validation.json
```

Require `status: pass` before downstream modeling. Read
[output-contract.md](references/output-contract.md) for required fields and
their meaning.

Regenerate a report after modifying a processed AnnData:

```bash
python scripts/render_report.py OUTPUT/processed.h5ad \
  --output-dir OUTPUT
```

## Interpretation Rules

- Distinguish expression neighbors in `obsp["connectivities"]` from physical
  neighbors in `obsp["spatial_connectivities"]`.
- Build physical graphs per sample or tissue section. Cross-section spatial edges
  are invalid unless an explicit registration method creates them.
- Treat `batch` nested within `time` as non-identifiable confounding. Report it;
  do not claim that automated integration recovered pure temporal biology.
- Preserve raw coordinates. Put display rotations or registered coordinates in
  additional `obsm` keys instead of silently replacing `obsm["spatial"]`.
- Keep random seeds and resolved configuration with every output.
