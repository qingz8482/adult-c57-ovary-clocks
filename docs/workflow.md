# Workflow

## Integrated clocks

1. Prepare/harmonize metadata with `integrated_clocks/shared/scripts/01_prepare_bulk_pseudobulk_metadata.R`.
2. Prepare and approve sample/group alignment with scripts `02b`–`02e`.
3. Fit both model families with `03_train_bulk_pseudobulk_clocks.R`, which sources the Version B selector and performs leakage-safe family-specific selection.
4. Run scripts `04`–`06` for locked diagnostics, primary prediction plots, frozen deployment projection, alignment audits, matched-control verification, and simplified projection classes.
5. Run ElasticNet interpretation scripts under `integrated_clocks/elasticnet/scripts/`.
6. Run PCA/glmnet interpretation and GSEA scripts under `integrated_clocks/pca_glmnet/scripts/`.
7. Build publication methods and supplementary tables with shared scripts `12` and `13`.

Feature scoring, scaling, PCA, glmnet tuning, and calibration must be refitted inside each training fold or split. Validation and test data must never influence those operations.

## Metacell clocks

Follow the order in `metacell_clocks/README.md`. Eligibility requires confirmed/probable labels, label confidence of at least 0.70, exclusion of Unknown and Erythroid, and at least 25 cells per metacell. All reported ages and errors use weeks.

## External projections

Deployment models are frozen before projection. Projection does not refit selected genes, scaling, PCA, glmnet coefficients, or calibration. The 300-sample external collection consists of 288 analysis-domain samples and 12 developmental/OOD placements. GSE106339 is classified by biological role, while its prepubertal contrast remains descriptive rather than part of the 34 primary adult perturbation contrasts.

## Cross-species analysis

Run `integrated_clocks/cross_species/run_all.sh`, then the GTEx adequacy and final-reporting runners in their numbered subdirectories. This workflow consumes frozen integrated-clock artifacts and produces the Supplementary Figure 7 analyses without refitting the mouse clocks.
