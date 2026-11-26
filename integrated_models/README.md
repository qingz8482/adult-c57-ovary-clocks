# Integrated clocks

The two integrated clocks reported in the manuscript. They were locked on 2026-10-12 (internal release `hvg_exclusion_vbpca_20261012`, as named in `manifest.json`) and have not been refitted since. Every model file is byte-identical to its completed final fit.

| Folder | Clock | Training rule | Training set | Settings | Non-zero predictors |
|---|---|---|---|---|---|
| `HVG-ElasticNet-exclude-single-age/` | EN | Studies with one distinct age excluded from training | 110 profiles, 9 studies, 8–83 wk | N = 3,000, α = 0.5, λ = 0.1789 | 79 genes |
| `VersionB-PCA-glmnet/` | PCA/EN | All studies | 136 profiles, 14 studies, 7–83 wk | N = 3,000, 10 PCs, α = 1, λ = 53.367 | PC1, PC2, PC5, PC8 |

Both use count tuning: N is chosen from 500/1k/2k/3k/5k by inner leave-one-study-out MAE (study-weighted), λ by the one-SE rule within each path, then the lowest-scoring path wins.

The EN clock excludes these studies from training: Morris (16 profiles, 7 wk), CRA017492 (3, 8.7 wk), GSE252791 (3, 12 wk), Zhao (1, 12.9 wk) and GSE194196 (3, 16 wk). The same rule was applied inside every outer training set of its evaluation.

## Files

Each model folder contains:

- `model.rds`: the deployable object used by `scripts/predict_frozen.R`;
- `genes.txt`, `scaler.csv`, `coefficients.csv` and `intercept.txt`: the model as plain text;
- for PCA/EN only, `pca_loadings.csv` and `pc_coefficients.csv`;
- `provenance/`: the inner-CV candidates.

For both clocks, the prediction in weeks is `intercept.txt` + Σ `coefficient` × x, where x is ln(1 + CPM) and unmeasured genes take their `center` value. `coefficient` is in weeks per ln(1 + CPM) unit, and `standardized_coefficient` is the same weight per training SD. For PCA/EN, the gene coefficients equal `pca_loadings.csv` × `pc_coefficients.csv` / scale. Supplementary Data 3c gives the same back-projected weights on the standardized scale (loadings × PC coefficients). `tests/test_integrated_models.R` checks that `model.rds` and the plain-text export give the same predictions.

## Relationship to the previous release

The earlier workspace release (`rebutal3 .../locked_models/count_tuning_20261011`) is unchanged. These clocks supersede it for deployment:

- **HVG–ElasticNet** is replaced by the exclusion model.
- **Version B–PCA/glmnet** is carried over unchanged. The files are byte-identical; see `manifest.json` sources.
- **Version B–ElasticNet** is not carried into this release.

## Performance (held-out, from model development, MAE in weeks)

| Model | Stage A pooled (4 folds) | A study-mean | B NIHMS | C Winkler | 30 repeated splits, median [IQR] |
|---|---|---|---|---|---|
| HVG–EN exclusion | 9.90 | 10.07 | 9.33 | 16.91 | 10.98 [10.6–11.3] |
| Version B–PCA | 15.73 | 15.44 | 21.30 | 17.73 | 15.70 [15.2–16.9] |

These outer estimates assess each full procedure with fold-specific models. The deployment fits have no independent test set in this cohort.

The exclusion rule was adopted after seeing the rebutal3 sensitivity results and the rebutal4 outer results. Its estimates therefore describe model development, not independent validation; external data is required for that. The two models were trained on different training sets, so their numbers are not a like-for-like comparison of recipes. See `provenance/rebutal4_comparison_vs_rebutal3.md`.

## Applying the frozen models

From this directory:

```bash
python3 verify_checksums.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript scripts/predict_frozen.R HVG-ElasticNet-exclude-single-age/model.rds new_logcpm.csv predictions.csv
```

Substitute `VersionB-PCA-glmnet/model.rds` for the PCA clock. Input is a genes-by-samples CSV/TSV: gene identifiers in the first column, then expression as ln(1 + 10⁶ × gene quantity / full available pre-filter sample total), with the documented gene mapping. Compute the denominator before subsetting to clock genes. The model supplies frozen centring, scaling, PCA projection where applicable, and training-mean imputation for missing genes. No feature selection, batch correction or calibration is applied at prediction. Output ages are in weeks.

The HVG–EN model's youngest training sample is 8 weeks. Treat predictions for samples outside 8–83 weeks as extrapolation.

## Verification and provenance

- `validation/replay_checks.csv`: each locked RDS is identical to its source object and reproduces its source predictions exactly on its own training profiles (110 + 136 replay predictions). The exported gene-scale coefficients reproduce the predictions within 1e-13 weeks. These replays are reproducibility checks, not performance estimates.
- `manifest.json`: source paths, SHA-256 of every copied file, settings, training rules and selection scope.
- `checksums.tsv`: SHA-256 of every release file that is byte-identical to the original release, which includes all model artifacts. Check with `verify_checksums.py`. This README was rewritten for the repository, and five text files (`manifest.json`, two `sessionInfo.txt`, `provenance/sample_metadata.csv`, `provenance/source_checksums.tsv`) had cluster file paths removed; their original checksums are listed in `checksums_sanitised_files.tsv`.
- `scripts/core.R` and `scripts/predict_frozen.R` are byte-identical across rebutal3, its locked release and rebutal4. `scripts/build_release.py` assembled this release, and `scripts/validate_release.R` produced the replay checks.
- `provenance/`: per-model final inner-CV candidates and splits, run signatures, configs and R sessions; shared sample metadata and input checksums; the rebutal4 design README and comparison; the rebutal3 results report; and the previous release's manifest.
