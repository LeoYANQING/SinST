# Multi-Sample and Temporal Policy

## Graph Construction

Treat each tissue section as a separate physical coordinate system. Build
`spatial_connectivities` block by block using `sample_id`. Do not connect cells
across samples before explicit registration.

Expression neighbors may span samples after normalization because they describe
transcriptional similarity, not physical adjacency.

## Gene Coordination

Resolve gene identifiers before concatenation. Preserve stable IDs and add gene
symbols as annotation. Do not merge genes solely because display symbols match.

Use a documented join:

- inner join for a shared modeling feature space;
- outer join only when missing genes are represented deliberately.

## HVG Selection

For multiple batches, select HVGs with a batch key and inspect
`highly_variable_nbatches`. A gene variable in only one batch may represent
biology, platform bias, or both.

Do not automatically require intersection HVGs across every batch; report the
stability distribution and choose a threshold appropriate to the experiment.

## Batch Correction

Default to no integration. Batch correction is a scientific model choice, not a
generic preprocessing cleanup.

Before Harmony, scVI, scANVI, or another method, inspect the design matrix:

- whether every time has multiple technical batches;
- whether a batch spans multiple biological conditions;
- whether donor, section, chip, and time are separable;
- which biological covariates must be preserved.

If each batch appears at only one time, batch is nested in time. Technical and
temporal effects are not identifiable from this dataset alone. Produce sensitivity
analyses, retain uncorrected results, and avoid claims of complete deconfounding.

## Temporal Metadata

Store both:

- an ordered human-readable stage label when useful;
- a numeric time value for trajectory or dynamics models.

Do not parse numeric time silently from arbitrary batch strings. Use an explicit
mapping when filenames encode complex stages.
