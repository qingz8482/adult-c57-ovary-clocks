# Workflow

## Environment variables

| Variable | Meaning |
|---|---|
| `OVARY_CLOCK_ROOT` | Workspace folder that contains `adult_c57_ovary_clocks/`, which holds the prepared inputs (`data_input/`, standardized Seurat objects, metacell preprocessing) |
| `OVARY_CLOCK_UPSTREAM_ROOT` | Folder holding the upstream model-ready count tables (`outputs/model_ready_counts/<study>.model_ready_counts.csv`) and FASTQ-derived count metadata |
| `RSCRIPT` | R interpreter used by the shell runners (default `Rscript`) |
| `TAGE_PYTHON` | Python environment with the official tAge package (tAge comparison only) |

Some scripts read sibling outputs by relative path, for example `../lockedmodel` or `../../rebutal3_batcheffects removal_CV_LOO/...`. Each repository folder must therefore sit at its workspace location, given below relative to `$OVARY_CLOCK_ROOT/adult_c57_ovary_clocks/rebutal/`. The workspace folder names are kept as used, so the provenance records inside the outputs stay valid.

| Repository folder | Workspace location |
|---|---|
| `integrated_clocks/01_input_preparation/scripts` | `rebutal3_batcheffects removal_CV_LOO/model_development/scripts` |
| `integrated_clocks/02_model_development` (EN) | `rebutal4/model_development` |
| `integrated_clocks/02_model_development` with `production_rebutal3_reference.json` as `config/production.json` (PCA/EN) | `rebutal3_batcheffects removal_CV_LOO/model_development` |
| `integrated_clocks/03_pca_interpretation` (except `backprojection/`) | `rebutal3_batcheffects removal_CV_LOO/model_development/VersionB-PCA/diagnostic` |
| `integrated_clocks/03_pca_interpretation/backprojection` | `rebutal4/backprojection_PCA` |
| `integrated_models` | `rebutal4/lockedmodel` |
| `integrated_clocks/04_external_projection` | `rebutal4/lockedmodel_application` |
| `integrated_clocks/05_tage_comparison` | `rebutal4/compare_tAge` |
| `integrated_clocks/06_cross_species` | `rebutal4/generalization_tohuman` |
| `metacell_clocks/00_label_audit` | `rebutal4/metacell_review` |
| `metacell_clocks/01_labelling_v2` | `rebutal4/metacell_review/v2` |
| `metacell_clocks/02_release_v2` | `rebutal4/metacell_review/reproduce` |
| `integrated_clocks/07_en_ora/scripts` | `rebutal4/figure1_model_development/scripts` |

## Run order

1. **Inputs.** Run `integrated_clocks/01_input_preparation/scripts/prepare_inputs.R`.
2. **Model development.** Run `integrated_clocks/02_model_development/scripts/run_all.py --phases='A;B,C;final;A_repeated'`, once with each config, then `summarize_results.R`. The final fits form `integrated_models`. Its `scripts/build_release.py` and `scripts/validate_release.R` assembled and replay-checked the release.
3. **PCA/EN interpretation.** Run `integrated_clocks/03_pca_interpretation/scripts/run_all.sh`, which calls `supple2_PCA/scripts/run_all.sh`.
4. **External projection.** Run `integrated_clocks/04_external_projection/scripts/run_all.sh`, then the figure scripts `08_clock_concordance_no_developmental_ood.R` and `14_…`–`17_…`, then `detectable_shift/detectable_shift.R`.
5. **tAge comparison.** Run `integrated_clocks/05_tage_comparison/scripts/run_all.sh`; it needs the projection outputs from step 4.
6. **Cross-species.** Run `integrated_clocks/06_cross_species/scripts/run_all.sh`.
7. **Metacell clocks.** Follow [`../metacell_clocks/README.md`](../metacell_clocks/README.md): run `01_labelling_v2`, then `02_release_v2`.
8. **EN non-zero-gene ORA.** Run `integrated_clocks/07_en_ora/scripts/enrich_elasticnet_nonzero_genes.R` in an R environment with clusterProfiler and org.Mm.eg.db.

[`manuscript_code_map.tsv`](manuscript_code_map.tsv) lists the analysis code behind each figure and supplementary table. Figure-assembly scripts are not included.

## Invariants

- All data-dependent steps run inside the training studies of each outer split: gene filtering, feature scoring and selection, centring and scaling, PCA and glmnet tuning. Held-out studies never influence fitting or selection.
- No age calibration, clipping or batch correction is applied, in development or at projection.
- Projection uses the frozen genes, centres, scales, PCA rotation and coefficients. Missing genes take the training mean, which contributes zero after scaling.
- Ages and errors are in weeks. Delta age is predicted minus chronological age.
- Metacell CV holds out metacells, not biological samples. The expected sample overlap between folds is audited.
