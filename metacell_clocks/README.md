# Cell-type-specific metacell clocks (v2 labels)

| Folder | What it does |
|---|---|
| `00_label_audit/` | Audits the labels behind the first-submission (v1) metacell clocks. Only Morris and Isola carried author labels; about 85% of clock cells had been labelled by one of three different in-house marker procedures. Luteal-like cells were also assigned inconsistently between granulosa and theca/stroma across datasets. This audit motivated the v2 relabelling. `tables/03_cluster_review_flags.tsv` holds the per-cluster review flags it reads. |
| `01_labelling_v2/` | One labelling procedure applied to all eight datasets, followed by metacell construction and the metacell clocks. |
| `02_release_v2/` | Figures and tables rebuilt on the v2 labels and clocks (Fig. 4, Fig. 5, Supplementary Figs 5–6). |

## 01_labelling_v2: run order

1. `01_marker_v2_scoring.R`: per-cluster marker-signature scores (marker panel v2.2, `reference_tables/01_marker_panel_v2.tsv`) z-scored across clusters, with an absolute detection gate.
2. `02a_build_reference.R`, `02b_reference_transfer.R`: Seurat anchor-based label transfer from the Morris and Isola atlases. Confidence is calibrated by cross-atlas prediction.
3. `03_benchmark.R`: benchmark against the author labels (`reference_tables/03_author_label_map.tsv`).
4. `04_consensus_v2.R`: clusters are confirmed when marker scoring and reference mapping agree; discordant clusters go to manual review. Luteal state within granulosa is assigned only when both evidence sources agree. The ontology is in `reference_tables/00_ontology_v2.tsv`.
5. `05_cross_dataset_consistency.R`, `06_luteal_theca_gene_check.R`: consistency checks.
6. `run_scRNA_metacell_elasticnet_build_v2.R`: builds metacells. Cells with confidence ≥ 0.70 that are not erythroid or unknown are split at random into groups of 50 within each dataset × sample × class, and raw counts are summed; groups with fewer than 25 cells are dropped.
7. `run_scRNA_metacell_elasticnet_metacell_cv_v2.R`: class eligibility (≥ 30 metacells, ≥ 3 age points, ≥ 2 age bins), then five-fold metacell-level nested CV balanced across age bin and dataset. Feature caps are 20/40/80/160/320 genes, and α/λ are tuned inside each outer training fold. The script also fits the final deployment model for each class.
8. `07_compare_v1_v2.R`, `08_plot_granulosa_luteal_sensitivity.R`: v1-versus-v2 comparison and the granulosa luteal-fraction sensitivity analysis. The `_sens` runners repeat steps 6–7 for that sensitivity analysis.

Folds are assigned by metacell, so metacells from one biological sample can fall in different folds. The fit audits this expected sample overlap. Fold assignments are in [`../metadata/pipeline_inputs/metacell_fold_assignments.tsv`](../metadata/pipeline_inputs/metacell_fold_assignments.tsv).

## 02_release_v2: run order

| # | Script | Output |
|---|---|---|
| 1 | `build_prjna863443_stress_test_contrasts.R` | Matched VCD/Foxl2 contrasts |
| 2 | `run_scRNA_metacell_elasticnet_coefficient_gsea.R` | Coefficient-ranked GSEA per class (GO BP and Reactome) |
| 3 | `organize_scRNA_metacell_cv_figure_package.R` | Overview, eligibility, internal CV, coefficient and GSEA figures and tables |
| 4 | `project_GSE235021_to_metacell_clocks.R` (with `v2_external_labeling.R`) | GSE235021 projection with v2 labels for the external libraries |
| 5 | `12_add_stress_projection_pvalues.R`, `13_make_horizontal_stress_projection_metacell_delta_composite.R` | Stress-test projection tests and composite |
| 6 | `refresh_filtered_elasticnet_metacell_internal_grouped_cv_plot.R` | Internal CV plot |
| 7 | `step_02_preprocessing_filtering/*.R` | Annotation evidence, marker panel, preprocessing and metacell summaries, and composition versus age (raw and dataset-adjusted) |
| 8 | `22_add_stress_heatmap_significance.R` | VCD/Foxl2 and GSE235021 delta heatmaps with Welch tests, plus the formal 6→9-month positive-control rule |
| 9 | `23_plot_gse235021_six_clock_faceted_welch.R` | GSE235021 per-clock diet contrasts (BH over clocks × contrasts) |
| 10 | `25_gse235021_alias_corrected_reprojection.R` | Independent check of the GSE235021 alias fix |
| 11 | `24_sfig6_deployment_fit_scatter.R`, `25_sfig6_vcd_foxl2_merged.R` | Deployment-fit scatter source tables (used by 26/27) |
| 12 | `26_…`, `27_…` | Per-clock fit-and-delta composites for PRJNA863443 and GSE235021 (Fig. 5, Supplementary Fig. 6C–P) |
| 13 | `28_…`, `29_…` | Stress-test cell-class composition |
| 14 | `plot_sf5_gsea_nes_artboard.R` | Supplementary Fig. 5 GSEA NES panels |

After any change to the GSE235021 projection, rerun scripts 4 → 8 → 9 → 10 → 11 → 12.

### GSE235021 identifier handling

The training counts are not symbol-harmonized across datasets. A few model features are non-zero in only one training dataset: old mm10 symbols and R-sanitized mitochondrial names from Morris. `project_GSE235021_to_metacell_clocks.R` zero-fills alias-resolved features that are non-zero in at most one training dataset (171 features). It reports coverage over deployable features, meaning those non-zero in at least two training datasets. This removes a near-constant offset in absolute predicted age. Within-dataset contrasts are essentially unchanged.
