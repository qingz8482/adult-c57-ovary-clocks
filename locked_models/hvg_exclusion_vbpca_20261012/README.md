# Locked clocks — release hvg_exclusion_vbpca_20261012

Locked 2026-10-12 (Asia/Singapore) by user selection. This release contains two deployment models, copied byte-for-byte from their completed final fits. Nothing was refitted. Do not overwrite this directory; create a new versioned release for any future change.

| Model | Training rule | Training set | Input genes | Nonzero predictors | Saved model |
|---|---|---|---:|---|---|
| HVG–ElasticNet (single-age exclusion) | studies with one distinct age excluded from training | 110 profiles, 9 studies, 8–83 wk | 3,000 | 79 genes | [model.rds](HVG-ElasticNet-exclude-single-age/model.rds) |
| Version B–PCA/glmnet | all studies (unchanged from release count_tuning_20261011) | 136 profiles, 14 studies, 7–83 wk | 3,000 | 4 of 10 PC inputs: PC1, PC2, PC5, PC8 | [model.rds](VersionB-PCA-glmnet/model.rds) |

Settings: HVG–EN N = 3,000, α = 0.5, λ = 0.1789. Version B–PCA N = 3,000, 10 PCs, α = 1, λ = 53.367. Both use count tuning: N is chosen from 500/1k/2k/3k/5k by inner leave-one-study-out MAE (study-weighted), λ by the one-SE rule within each path, then the lowest-scoring path wins.

The HVG–EN model excludes these studies from training: Morris (16 profiles, 7 wk), CRA017492 (3, 8.7 wk), GSE252791 (3, 12 wk), Zhao (1, 12.9 wk) and GSE194196 (3, 16 wk). The same rule was applied inside every outer training set of its evaluation.

## Relationship to the previous release

`rebutal3 .../locked_models/count_tuning_20261011` is unchanged. This release supersedes it for deployment:

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
- `checksums.tsv`: SHA-256 of every release file that is byte-identical to the original release, which includes all model artifacts. Check with `verify_checksums.py`. Six text files (`README.md`, `manifest.json`, two `sessionInfo.txt`, `provenance/sample_metadata.csv`, `provenance/source_checksums.tsv`) had cluster file paths removed for publication; their original checksums are listed in `checksums_sanitised_files.tsv`.
- `scripts/core.R` and `scripts/predict_frozen.R` are byte-identical across rebutal3, its locked release and rebutal4. `scripts/build_release.py` assembled this release, and `scripts/validate_release.R` produced the replay checks.
- `provenance/`: per-model final inner-CV candidates and splits, run signatures, configs and R sessions; shared sample metadata and input checksums; the rebutal4 design README and comparison; the rebutal3 results report; and the previous release's manifest.
