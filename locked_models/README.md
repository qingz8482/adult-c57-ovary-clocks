# Locked models

[`hvg_exclusion_vbpca_20261012/`](hvg_exclusion_vbpca_20261012/README.md) holds the two integrated clocks reported in the manuscript. They were locked on 2026-10-12 and have not been refitted since.

| Folder | Clock | Training set | Settings | Non-zero predictors |
|---|---|---|---|---|
| `HVG-ElasticNet-exclude-single-age/` | EN | 110 profiles, 9 studies, 8–83 wk | N = 3,000, α = 0.5, λ = 0.1789 | 79 genes |
| `VersionB-PCA-glmnet/` | PCA/EN | 136 profiles, 14 studies, 7–83 wk | N = 3,000, 10 PCs, α = 1, λ = 53.367 | PC1, PC2, PC5, PC8 |

Each model folder contains:

- `model.rds`: the deployable object used by `scripts/predict_frozen.R`;
- `genes.txt`, `scaler.csv`, `coefficients.csv` and `intercept.txt`: the model as plain text;
- for PCA/EN only, `pca_loadings.csv` and `pc_coefficients.csv`;
- `provenance/`: the inner-CV candidates.

For both clocks, the prediction in weeks is `intercept.txt` + Σ `coefficient` × x, where x is ln(1 + CPM) and unmeasured genes take their `center` value. `coefficient` is in weeks per ln(1 + CPM) unit, and `standardized_coefficient` is the same weight per training SD. For PCA/EN, the gene coefficients equal `pca_loadings.csv` × `pc_coefficients.csv` / scale. Supplementary Data 3c gives the same back-projected weights on the standardized scale (loadings × PC coefficients). `tests/test_locked_models.R` checks that `model.rds` and the plain-text export give the same predictions.
