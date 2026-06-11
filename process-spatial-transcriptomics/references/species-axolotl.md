# Axolotl Stereo-seq Profile

## ARTISTA Characteristics

ARTISTA-style files use:

- `obs["Batch"]` for section/library identity;
- `var["Gene"]` for mixed human-readable annotation;
- stable IDs containing `AMEX60DD...`;
- `obsm["spatial"]` and sometimes `obsm["X_spatial"]`;
- integer counts in `layers["counts"]`;
- source `X` that may equal `log1p(counts)` without library-size normalization.

Therefore, rebuild normalized `X` from `layers["counts"]` rather than trusting
the presence of precomputed PCA or a float-valued `X`.

## Gene Annotation

Do not assume human or mouse naming conventions. Entries can combine an axolotl
stable ID with inferred ortholog labels. Preserve `var_names` as stable feature
identifiers and copy `var["Gene"]` to `var["gene_symbol"]` for display.

`MT-` prefix matching can detect nothing in axolotl files. Skip mitochondrial
filtering unless a trusted axolotl mitochondrial gene mapping is supplied through
the configuration.

## Coordinates

Preserve raw chip coordinates in `obsm["spatial"]`. ARTISTA `uns["angle_dict"]`
contains display rotations; applying those rotations is visualization metadata,
not biological registration.

Store rotated or registered coordinates under separate names such as:

- `obsm["spatial_display"]`;
- `obsm["spatial_registered"]`.

## Samples and Names

ARTISTA combined files can contain duplicated observation names. Canonicalize to
`sample_id::source_obs_name` and preserve the original name.

Construct Delaunay graphs separately for every `Batch`. Pooled files do not
necessarily retain source spatial graph matrices even when individual section
files do.

## Confounding

In ARTISTA regeneration data, sequencing chip and time can be strongly or fully
confounded. Treat the pipeline's nesting warning as a hard interpretation
constraint. Do not enable generic batch correction and describe the result as
pure regeneration dynamics.

Use explicit time mappings for values such as 2, 5, 10, 15, 20, 30, and 60 DPI.
The profile records and applies a DPI regular expression to `sample_id`; replace
it with `canonical.time_mapping` when sample names are irregular. Do not rely on
lexical ordering of strings such as `10DPI` and `5DPI`.
