# Input Contract

## Supported Input

The deterministic runner accepts one AnnData `.h5ad` file. Convert other platform
formats to AnnData before invoking the runner.

Required biological data:

- observations in rows and genes in columns;
- integer-like raw counts in `layers["counts"]`, another configured layer, or `X`;
- physical coordinates in `obsm["spatial"]` or a configured alternative;
- a sample/section identifier for concatenated datasets.

Recommended metadata:

- batch or library identifier;
- time point for temporal data;
- gene symbols in `var`;
- cell type or region annotation when available.

## State Detection

The inspector samples matrices and classifies common states:

- `raw_counts_in_X`: `X` is integer-like and can seed `layers["counts"]`;
- `X_equals_log1p_counts`: `expm1(X)` matches an integer count layer;
- `counts_layer_available_X_other`: raw counts exist but `X` has another state;
- `unknown_no_integer_counts`: processing must stop until raw counts are supplied.

Do not reverse log transformation to manufacture counts when no raw count layer
exists. `expm1(X)` is acceptable only as a diagnostic equality check against an
existing count matrix.

## Canonical Keys

The pipeline writes:

- `obs["sample_id"]`;
- `obs["batch"]`;
- `obs["time"]` when available;
- `var["gene_symbol"]` when available;
- `obsm["spatial"]`.

When observation names are duplicated, it preserves the original value in
`obs["source_obs_name"]` and prefixes names with `sample_id`.

Time can be supplied through an existing `obs` column, an explicit
`canonical.time_mapping`, or a configured `canonical.time_regex`. The resolved
configuration and provenance record which rule was used.

## Refusal Conditions

Stop rather than guess when:

- no integer-like raw counts are available;
- spatial coordinates contain NaN or have unsupported dimensions;
- no sample identifier exists for a concatenated multi-section object;
- QC removes nearly all observations or genes;
- metadata mappings assign the same cell to contradictory samples or times.
