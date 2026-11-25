# Metadata

## `manuscript_tables/`: the published Supplementary Data

TSV exports of the final Supplementary Data workbook (Supplementary Data 1–5). Title and header rows are kept as in the workbook. `README.tsv` describes each table.

| File | Content |
|---|---|
| `SD1_datasets.tsv` | Dataset summary (30 dataset-modality rows plus the human GSE202601 cohort), dataset roles and age-bin definitions |
| `SD2a_model_development_profiles.tsv` | The 136 development profiles with Stage A fold, Stage B/C holdout and final-fit membership for each clock |
| `SD2b_external_projection_profiles.tsv` | The 281 external projection profiles with projection class |
| `SD3a_EN_nonzero_coefficients.tsv` | The 79 non-zero EN coefficients |
| `SD3b_PCA_EN_pc_coefficients.tsv` | PCA/EN PC coefficients and variance explained |
| `SD3c_PCA_EN_backprojected_gene_weights.tsv` | Back-projected PCA/EN gene weights for the 3,000 input genes |
| `SD4_benchmarked_contrasts.tsv` | The 36 benchmarked contrasts: source-study phenotype, expected direction and both clocks' delta-age effect, Welch P and BH q |
| `SD5_marker_panel_v2.tsv` | Marker panel v2.2 used for cell-class annotation |

## `pipeline_inputs/`: the tables the code reads

| File | Read by | Content |
|---|---|---|
| `integrated_development_samples.csv` | `integrated_clocks/01–02` | 136 development profiles (14 studies). Local file-path columns were removed. |
| `integrated_fold_assignments.csv` | `integrated_clocks/02` | Stage A fold (`cv_fold`) and the 30 repeated partitions (`repeat_01`–`repeat_30`) |
| `integrated_fold_balance_audit.tsv` | — | Train/test size, age range and modality balance for every fold and repeat |
| `external_projection_samples.tsv` | `integrated_clocks/04` | All 300 candidate projection rows before final exclusions, including reference controls and developmental OOD placements |
| `external_projection_contrasts.tsv`, `prespecified_contrasts.tsv` | `integrated_clocks/04–05` | All 41 contrast definitions, including placement-only and developmental contrasts. Supplementary Data 4 holds the 36 that were benchmarked. |
| `metacell_fold_assignments.tsv` | `metacell_clocks/01–02` | 7,049 eligible metacells with v2 class, label confidence, luteal fraction and outer fold |
| `metacell_cv_eligibility.tsv` | — | Per-class eligibility (eight eligible classes; oocyte excluded) |

The pipeline tables are a superset of the manuscript tables. Where counts differ, the manuscript tables reflect the final reporting filters. `sample_key` joins the integrated tables, and `metacell_id` and `sample_key` join the metacell tables. Empty values mean unavailable or not applicable.

## `dataset_accessions.tsv`

Dataset-to-repository and accession lookup.
