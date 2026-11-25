# Model development results

Execution status: **complete**. Completed controller phases: A, B,C, final, sensitivities, A_repeated.

All results here use the repaired-input derivative, training-only expression filtering, full nested study-grouped validation and no explicit batch correction. The six procedures were prespecified; this report does not select a winning procedure.

## Main Analysis A: 136 profiles, four whole-study folds

| Recipe | Procedure | Pooled MAE (weeks) | Mean study MAE (weeks) | Predictive R² | Mean bias (weeks) |
|---|---|---:|---:|---:|---:|
| HVG–ElasticNet | Fixed 2,000 | 16.96 | 17.93 | -0.02 | 1.36 |
| HVG–ElasticNet | Nested count tuning | 16.10 | 16.17 | 0.07 | 1.81 |
| Version B–ElasticNet | Fixed 2,000 | 15.75 | 17.77 | 0.08 | -0.28 |
| Version B–ElasticNet | Nested count tuning | 15.00 | 16.79 | 0.16 | -0.31 |
| Version B–PCA/glmnet | Fixed 2,000 | 16.23 | 16.01 | 0.03 | 0.60 |
| Version B–PCA/glmnet | Nested count tuning | 15.73 | 15.44 | 0.10 | -0.67 |

These are study-held-out estimates for the specified modeling procedures. They should not be compared as if interchangeable with the earlier sample-level repeated-CV MAEs. Negative predictive R² means squared prediction error exceeds the observed cohort-mean reference used in the R² denominator.

## Independent study holdouts

| Recipe | Procedure | NIHMS B MAE (weeks) | Winkler C MAE (weeks) |
|---|---|---:|---:|
| HVG–ElasticNet | fixed_2000 | 24.44 | 19.53 |
| HVG–ElasticNet | count_tuning | 23.18 | 17.65 |
| Version B–ElasticNet | fixed_2000 | 23.40 | 17.12 |
| Version B–ElasticNet | count_tuning | 23.39 | 17.12 |
| Version B–PCA/glmnet | fixed_2000 | 20.92 | 17.73 |
| Version B–PCA/glmnet | count_tuning | 21.30 | 17.73 |

Each holdout has its own training-only tuning and model. Neither uses the final all-study model.

## Repeated study partitions

| Recipe | Procedure | Completed repeats | Median pooled MAE | IQR | Median study MAE | IQR |
|---|---|---:|---:|---|---:|---|
| HVG–ElasticNet | fixed_2000 | 30 | 16.48 | 15.56–16.68 | 16.40 | 15.87–16.81 |
| HVG–ElasticNet | count_tuning | 30 | 15.96 | 15.26–17.07 | 15.95 | 15.09–16.95 |
| Version B–ElasticNet | fixed_2000 | 30 | 15.75 | 15.34–16.13 | 16.41 | 16.13–16.96 |
| Version B–ElasticNet | count_tuning | 30 | 15.74 | 14.70–16.29 | 16.19 | 15.80–16.71 |
| Version B–PCA/glmnet | fixed_2000 | 30 | 16.15 | 15.54–16.72 | 15.36 | 14.80–15.87 |
| Version B–PCA/glmnet | count_tuning | 30 | 15.70 | 15.19–16.93 | 14.86 | 14.13–15.67 |

These IQRs describe dependence on study grouping. The repeated predictions are not independent biological replicates, and the IQRs are not confidence intervals.

## Prespecified sensitivities: main-A pooled MAE in weeks

| Recipe | Procedure | Inverse-stratum weights | Omit single-age training studies | Log-age target | Verified inputs (78 profiles) |
|---|---|---:|---:|---:|---:|
| HVG–ElasticNet | fixed_2000 | 15.99 | 11.52 | 16.40 | 15.39 |
| HVG–ElasticNet | count_tuning | 14.78 | 9.90 | 15.90 | 14.97 |
| Version B–ElasticNet | fixed_2000 | 15.29 | 14.12 | 16.16 | 16.25 |
| Version B–ElasticNet | count_tuning | 14.62 | 14.13 | 14.88 | 16.21 |
| Version B–PCA/glmnet | fixed_2000 | 15.40 | 16.11 | 16.05 | 16.92 |
| Version B–PCA/glmnet | count_tuning | 15.60 | 15.23 | 15.76 | 16.66 |

Historical inverse-stratum weights are missing for 34 profiles; this sensitivity recomputes inverse training-stratum frequencies from the existing stratum labels within each training split. The verified-input sensitivity has different study, age and modality coverage. Its MAE is not a paired comparison on the full 136-profile population. Each other sensitivity changes one aspect of the primary recipe.

## Final model settings

| Recipe | Procedure | N genes | PCs | Alpha | Lambda | Nonzero effective genes |
|---|---|---:|---:|---:|---:|---:|
| HVG–ElasticNet | fixed_2000 | 2000 | 0 | 0.1 | 38.54 | 46 |
| HVG–ElasticNet | count_tuning | 3000 | 0 | 0.5 | 0.9112 | 61 |
| Version B–ElasticNet | fixed_2000 | 2000 | 0 | 0.75 | 7.565 | 8 |
| Version B–ElasticNet | count_tuning | 5000 | 0 | 0.1 | 53.37 | 64 |
| Version B–PCA/glmnet | fixed_2000 | 2000 | 10 | 0.75 | 73.91 | 2000 |
| Version B–PCA/glmnet | count_tuning | 3000 | 10 | 1 | 53.37 | 3000 |

Final all-study models are saved for deployment; their training errors are deliberately not presented as performance. A lambda.1se-style procedure can select an intercept-only model when study-transfer uncertainty is large; zero effective genes are reported explicitly.

## Validation and provenance

- Each complete Analysis A has exactly one prediction for every intended profile. All inner and outer splits keep each study intact.
- The independent checks compare Version B with gene-wise regressions, perturb held-out data, verify the common 100-lambda grid, and reproduce predictions from saved models and original-gene coefficients.
- Input recovery agrees with existing filtered exports. The CPM denominator now uses each profile’s full available pre-filter total; see [normalization comparison](shared_inputs/normalization_comparison.tsv). Existing input files and historical models were preserved.
- The benchmark retains 78 verified-count profiles, 26 processed-expression profiles and 32 unresolved profiles. No statistical transformation changes those provenance labels.
- Initial validation fits affected by floating-point lambda grouping are isolated under `validation_superseded_lambda_key_bug/` and `validation_first_production_fold/`. They must not be used as results. Current production restores the exact prespecified lambda vector and checks complete cross-fold grid aggregation.
- LOWESS age gaps are descriptive out-of-fold residualizations using observed ages, not a deployment calibration or independent biological-age validation.

See [comparison tables](comparison/), [pipeline documentation](README.md), [unit-check results](tests/test_results.csv) and [saved-fit checks](tests/completed_run_checks.csv).

Completed validation: 29 implementation checks, 858 saved-model replay checks, 28,494 held-out predictions and 4,431 study-disjoint inner partitions.
