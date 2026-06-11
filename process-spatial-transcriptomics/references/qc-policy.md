# QC Policy

## Principles

Calculate QC from raw counts, not normalized `X`. Preserve metrics before
filtering so removals remain auditable.

Default filters are deliberately conservative:

- minimum detected genes per observation;
- minimum observations per gene;
- mitochondrial percentage only when mitochondrial genes are actually identified.

Upper count/gene thresholds default to `null`. Set them from observed
distributions or validated platform expectations. Do not remove the top fixed
percentile merely because it is convenient.

## Mitochondrial Genes

Use species-appropriate gene annotation. If no trustworthy mitochondrial mapping
exists, store `pct_counts_mt` as missing and skip mitochondrial filtering. Never
write zero percent as though it were measured.

## Doublets

High counts or detected genes are screening signals, not proof of a doublet.
Add a platform-appropriate doublet method as a separate, reported stage when
needed. Do not mix inferred doublet calls into generic QC without recording the
method and score.

## Spatial QC

Inspect each section independently for:

- isolated coordinates;
- duplicated coordinates;
- disconnected tissue fragments;
- unusually long graph edges;
- section orientation and coordinate units.

Distance pruning removes unusually long Delaunay edges but can also disconnect
real narrow structures. Review the spatial figure and per-sample graph metrics.

## Threshold Changes

For each changed cutoff, record:

- biological/platform justification;
- number and percentage removed per sample;
- whether removal is concentrated in one batch, time, or cell type;
- before/after QC figures.
