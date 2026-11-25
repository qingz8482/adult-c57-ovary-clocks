# Adult C57 ovary clocks

Code, locked models and metadata for three mouse ovarian transcriptomic ageing clocks:

1. **EN clock**: an integrated bulk RNA-seq plus scRNA-seq pseudobulk clock (HVG selector → ElasticNet), trained without single-age studies;
2. **PCA/EN clock**: an integrated clock (Version B age/dataset-penalized selector → PCA → ElasticNet); and
3. **cell-type-specific metacell clocks**: ElasticNet clocks for eight harmonized ovarian cell classes, built on v2 cell labels.

This release matches the revised manuscript. Version `v0.1-initial-submission` holds the superseded code from the original submission, which used a fixed 83/14/20 train/validation/test split.

## Study design in brief

- **Model development.** 136 control-like profiles from 14 studies (69 bulk RNA-seq profiles from seven studies, 67 scRNA-seq pseudobulk profiles from seven studies; 7–83 weeks). Development used study-grouped nested cross-validation:
  - **Stage A:** four age- and modality-balanced folds, plus 30 repeated near-balanced partitions.
  - **Stages B and C:** NIHMS1554107 or the Winkler atlas held out.
  - **Inner loop:** leave-one-study-out tuning with the one-SE λ rule.
- **Training rule.** The EN clock excludes the five single-age studies (26 profiles) from every training set and from the final fit (110 profiles, 9 studies, 8–83 weeks). The PCA/EN clock uses all 136 profiles.
- **Locked clocks.** Release `hvg_exclusion_vbpca_20261012`:
  - EN: N = 3,000 genes, α = 0.5, λ = 0.179, 79 non-zero genes.
  - PCA/EN: N = 3,000 genes, 10 PC inputs, α = 1, λ = 53.4; active PCs are PC1, PC2, PC5 and PC8.
  - No calibration or batch correction is applied.
- **Applications.**
  - Frozen projection of intervention cohorts, with matched-control Welch tests and phenotype benchmarking.
  - Comparison with the multi-tissue tAge clocks.
  - Two-stage cross-species evaluation in human ovary (GTEx and GSE202601).
- **Metacell clocks.**
  - Cells from eight scRNA-seq datasets were re-annotated with one procedure that combines marker panel v2.2, reference transfer and consensus rules.
  - Metacells were formed from 50 cells (minimum 25 retained).
  - Each class had a five-fold metacell-level nested CV; 7,049 metacells were modelled across eight classes.

## Repository map

| Folder | Contents |
|---|---|
| [`locked_models/`](locked_models/README.md) | The two deployed integrated clocks: `model.rds`, exported coefficients, scalers, PCA loadings, provenance, replay checks and checksums. |
| [`integrated_clocks/`](integrated_clocks/README.md) | Input preparation, nested-CV model development, PCA/EN interpretation, external projection, tAge comparison and cross-species analyses. |
| [`metacell_clocks/`](metacell_clocks/README.md) | v1 label audit, v2 labelling and metacell clock fitting, and the v2 figure and table release, including stress-test projections. |
| [`figures/`](figures/README.md) | Assembly scripts for Fig. 1, Fig. 4e and Supplementary Figs 1–3. |
| [`metadata/`](metadata/README.md) | Supplementary Tables 1–5 as TSV, plus the sample, fold and contrast tables the pipelines read. |
| [`docs/`](docs/workflow.md) | Workflow, figure-to-code map, data access, environment and reproducibility boundary. |
| [`tests/`](tests/README.md) | Repository validation, metadata checks and a locked-model smoke test. |

## Applying the locked clocks

```bash
cd locked_models/hvg_exclusion_vbpca_20261012
python3 verify_checksums.py
Rscript scripts/predict_frozen.R HVG-ElasticNet-exclude-single-age/model.rds new_logcpm.csv predictions.csv
Rscript scripts/predict_frozen.R VersionB-PCA-glmnet/model.rds new_logcpm.csv predictions.csv
```

The input is a genes-by-samples CSV/TSV with gene symbols in the first column. Values are ln(1 + CPM), with CPM computed from the full pre-filter library total. Missing model genes are imputed with the training mean. Output ages are in weeks. Predictions outside the training age range (EN: 8–83 weeks; PCA/EN: 7–83 weeks) are extrapolations. See [`locked_models/hvg_exclusion_vbpca_20261012/README.md`](locked_models/hvg_exclusion_vbpca_20261012/README.md).

## Installation

```bash
conda env create -f environment.yml
conda activate adult-c57-ovary-clocks
```

Alternatively, restore the R dependencies with `renv::restore()`. Versions used for the locked fits are recorded in each model's `provenance/sessionInfo.txt` (R 4.4.3, glmnet 5.0, data.table 1.17.8). See [`docs/environment.md`](docs/environment.md).

## Running the analyses

Expression matrices, Seurat objects and generated results are not distributed. The scripts were run inside a fixed analysis workspace. Set:

```bash
export OVARY_CLOCK_ROOT=/path/to/workspace              # folder holding adult_c57_ovary_clocks/
export OVARY_CLOCK_UPSTREAM_ROOT=/path/to/upstream_counts # legacy model-ready count tables
export RSCRIPT=Rscript
```

[`docs/workflow.md`](docs/workflow.md) gives the run order and the workspace folder each component expects. [`docs/data_access.md`](docs/data_access.md) lists the public accessions.

Validate a clone:

```bash
bash tests/validate_repository.sh
bash tests/run_smoke_tests.sh
```

## Citation and license

Citation metadata are in [`CITATION.cff`](CITATION.cff). The code is released under the [MIT License](LICENSE). Public datasets, MSigDB/KEGG gene sets and the third-party tAge models keep their original licenses and terms.
