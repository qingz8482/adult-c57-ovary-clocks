# Metacell clocks

This directory contains the active July 2026 cell-type-specific metacell ElasticNet workflow.

The primary evaluation unit is a held-out `metacell_id`. Biological samples can contribute distinct metacells to multiple folds, and this overlap is explicitly audited. The resulting performance is therefore metacell-level internal out-of-fold CV, not animal-level independent validation.

## Inputs

Set `OVARY_CLOCK_DATA_ROOT` to a prepared data directory and `OVARY_METACELL_PREPROCESS_ROOT` to a directory containing:

- `step_02_preprocessing_filtering/filtered_metacell_counts.rds`
- `step_02_preprocessing_filtering/filtered_metacell_metadata.tsv`

Set `OVARY_CLOCK_OUTPUT_ROOT` for generated outputs. The canonical run is written to `$OVARY_CLOCK_OUTPUT_ROOT/metacell_clocks`.

## Order

1. `run_scRNA_metacell_elasticnet_metacell_cv.R`
2. `build_prjna863443_stress_test_contrasts.R`
3. `project_GSE235021_to_metacell_clocks.R` when its prepared Seurat inputs are available
4. `run_scRNA_metacell_elasticnet_coefficient_gsea.R`
5. `organize_scRNA_metacell_cv_figure_package.R`

PRJNA863443 and GSE235021 are projection-only stress tests. Their samples do not enter clock fitting or internal CV.
