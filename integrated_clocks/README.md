# Integrated clocks

These workflows use one expression vector per biological sample after intersecting genes across bulk RNA-seq and whole-sample scRNA pseudobulk profiles.

- `shared/scripts/` contains metadata preparation, leakage-safe feature selection, model training, benchmark diagnostics, frozen external projection, validation, and publication-table builders used by both integrated models.
- `elasticnet/scripts/` contains direct ElasticNet-specific gene interpretation and sensitivity analyses.
- `pca_glmnet/scripts/` contains PCA/glmnet component interpretation, GSEA, and contribution analyses.
- `cross_species/` contains the mouse-to-human direction-conservation, zero-shot projection, GTEx adequacy, and final reporting workflow used for Supplementary Figure 7.

The canonical selector lock is model-family-specific: ElasticNet uses the 2,000-gene HVG baseline, while PCA/glmnet uses the 2,000-gene training-only Version B selector. Validation selected the family-specific candidates; the held-out test set did not participate in selector choice.

Prepared data should be placed outside Git and exposed through `OVARY_CLOCK_DATA_ROOT`. Shared results default to `integrated_clocks/shared/results/`; set `CLOCK_RESULTS_DIR` to redirect them.
