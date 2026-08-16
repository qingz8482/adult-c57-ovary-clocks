#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
results_dir <- Sys.getenv("CLOCK_RESULTS_DIR", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = FALSE)
diag_dir <- file.path(results_dir, "integrated_diagnostics")
stopifnot(dir.exists(diag_dir))

layout <- data.table(
  section = c(
    "00_summary_index",
    "01_primary_predictions",
    "02_locked_clock_selection",
    "03_prediction_scatter_residuals",
    "04_dataset_bias",
    "05_stratified_performance",
    "06_naive_baseline",
    "07_paired_clock_comparison",
    "08_dynamic_range_calibration",
    "09_selected_gene_diagnostics",
    "10_pca_diagnostics",
    "11_projection_readiness",
    "12_PCA deployment model GSEA"
  ),
  description = c(
    "Narrative summary, warnings, plot manifest, and organization manifest.",
    "Example-style primary prediction plots for the two frozen integrated clocks, with training, validation, and held-out test facets.",
    "Family-specific locked selectors, candidate selector validation/test MAE, and prediction-row audit.",
    "Observed-vs-predicted and residual-vs-age diagnostics for locked clocks.",
    "Dataset-level residual bias and absolute-error diagnostics.",
    "Performance stratified by modality, biological age bin, and split age stratum.",
    "Naive train-mean and train-median baseline comparison.",
    "Paired sample-level comparison between retained ElasticNet and PCA/glmnet clocks.",
    "Dynamic range and calibration summaries.",
    "Integrated selected-gene and Version B feature-score diagnostics.",
    "PCA variance, coefficient, and score diagnostics.",
    "External projection readiness audit.",
    "PCA deployment model gene-set enrichment analysis and leading-edge summaries."
  )
)

for (section in layout$section) {
  dir.create(file.path(diag_dir, section), recursive = TRUE, showWarnings = FALSE)
}

rules <- list(
  "00_summary_index" = c(
    "integrated_clock_diagnostic_summary.md",
    "diagnostic_warnings.tsv",
    "plot_manifest.tsv"
  ),
  "01_primary_predictions" = c(
    "integrated_primary_predictions.pdf",
    "integrated_primary_predictions.png",
    "integrated_primary_predictions.svg",
    "integrated_primary_predictions_fit_stats.tsv",
    "integrated_primary_predictions_source.tsv"
  ),
  "02_locked_clock_selection" = c(
    "01_family_specific_locked_clocks.pdf",
    "01_family_specific_locked_clocks.png",
    "02_selector_validation_test_mae.pdf",
    "02_selector_validation_test_mae.png",
    "family_specific_locked_clock_summary.tsv",
    "prediction_row_audit.tsv"
  ),
  "03_prediction_scatter_residuals" = c(
    "03_locked_predicted_vs_observed.pdf",
    "03_locked_predicted_vs_observed.png",
    "04_locked_residual_vs_age.pdf",
    "04_locked_residual_vs_age.png",
    "04_locked_residual_vs_age_fit_stats.tsv",
    "04_locked_residual_vs_age_source.tsv",
    "worst_residual_samples.tsv"
  ),
  "04_dataset_bias" = c(
    "05_locked_abs_error_by_dataset.pdf",
    "05_locked_abs_error_by_dataset.png",
    "05_locked_residual_by_dataset.pdf",
    "05_locked_residual_by_dataset.png",
    "residual_bias_by_dataset.tsv"
  ),
  "05_stratified_performance" = c(
    "06_locked_error_by_modality.pdf",
    "06_locked_error_by_modality.png",
    "07_locked_error_by_biological_age_bin.pdf",
    "07_locked_error_by_biological_age_bin.png",
    "08_locked_error_by_split_age_stratum.pdf",
    "08_locked_error_by_split_age_stratum.png",
    "performance_by_biological_age_bin.tsv",
    "performance_by_modality.tsv",
    "performance_by_split_age_stratum.tsv"
  ),
  "06_naive_baseline" = c(
    "09_naive_baseline_comparison.pdf",
    "09_naive_baseline_comparison.png",
    "naive_baseline_performance.tsv"
  ),
  "07_paired_clock_comparison" = c(
    "10_paired_retained_clock_comparison_by_age_bin.pdf",
    "10_paired_retained_clock_comparison_by_age_bin.png",
    "10_paired_retained_clock_comparison_by_dataset.pdf",
    "10_paired_retained_clock_comparison_by_dataset.png",
    "10_paired_retained_clock_comparison_by_modality.pdf",
    "10_paired_retained_clock_comparison_by_modality.png",
    "paired_retained_clock_sample_comparison.tsv",
    "paired_retained_clock_summary.tsv",
    "paired_retained_clock_summary_by_age_bin.tsv",
    "paired_retained_clock_summary_by_dataset.tsv",
    "paired_retained_clock_summary_by_modality.tsv"
  ),
  "08_dynamic_range_calibration" = c(
    "11_dynamic_range_calibration.pdf",
    "11_dynamic_range_calibration.png",
    "calibration_dynamic_range_summary.tsv"
  ),
  "09_selected_gene_diagnostics" = c(
    "12_top30_versionB_genes_integrated_elasticnet.pdf",
    "12_top30_versionB_genes_integrated_elasticnet.png",
    "12_top30_versionB_genes_integrated_pca_glmnet.pdf",
    "12_top30_versionB_genes_integrated_pca_glmnet.png",
    "12_versionB_partial_r2_integrated_elasticnet.pdf",
    "12_versionB_partial_r2_integrated_elasticnet.png",
    "12_versionB_partial_r2_integrated_pca_glmnet.pdf",
    "12_versionB_partial_r2_integrated_pca_glmnet.png",
    "selected_gene_summary.tsv",
    "top30_versionB_genes.tsv"
  ),
  "10_pca_diagnostics" = c(
    "13_pca_diagnostics_pc1_pc2_scores.pdf",
    "13_pca_diagnostics_pc1_pc2_scores.png",
    "13_pca_diagnostics_pc_coefficients.pdf",
    "13_pca_diagnostics_pc_coefficients.png",
    "13_pca_diagnostics_variance_explained.pdf",
    "13_pca_diagnostics_variance_explained.png",
    "pca_glmnet_pc_coefficients.tsv",
    "pca_score_summary_from_rds.tsv",
    "pca_variance_summary.tsv"
  ),
  "11_projection_readiness" = c(
    "14_projection_readiness_audit.pdf",
    "14_projection_readiness_audit.png",
    "projection_readiness_audit.tsv"
  )
)

move_log <- data.table(
  file = character(),
  section = character(),
  old_path = character(),
  new_path = character(),
  status = character()
)

for (section in names(rules)) {
  for (file in rules[[section]]) {
    old_path <- file.path(diag_dir, file)
    new_path <- file.path(diag_dir, section, file)
    source_paths <- old_path
    if (section != "00_summary_index") {
      source_paths <- c(source_paths, file.path(diag_dir, "00_summary_index", file))
    }
    source_paths <- source_paths[file.exists(source_paths)]
    source_paths <- source_paths[normalizePath(source_paths, mustWork = FALSE) != normalizePath(new_path, mustWork = FALSE)]
    source_path <- if (length(source_paths) > 0) source_paths[[1]] else NA_character_
    if (!is.na(source_path) && file.exists(new_path)) {
      ok <- file.copy(source_path, new_path, overwrite = TRUE)
      if (ok) unlink(source_path)
      status <- if (ok) "replaced_existing" else "replace_failed"
    } else if (!is.na(source_path)) {
      ok <- file.rename(source_path, new_path)
      status <- if (ok) "moved" else "move_failed"
    } else if (file.exists(new_path)) {
      status <- "already_in_section"
    } else {
      status <- "not_found"
    }
    move_log <- rbind(
      move_log,
      data.table(file = file, section = section, old_path = old_path, new_path = new_path, status = status)
    )
  }
}

top_level_files <- list.files(diag_dir, full.names = FALSE, recursive = FALSE)
top_level_files <- top_level_files[file.info(file.path(diag_dir, top_level_files))$isdir == FALSE]
extra_files <- setdiff(top_level_files, c("README.md", "organization_manifest.tsv"))
if (length(extra_files) > 0) {
  section <- "00_summary_index"
  for (file in extra_files) {
    old_path <- file.path(diag_dir, file)
    new_path <- file.path(diag_dir, section, file)
    ok <- file.rename(old_path, new_path)
    move_log <- rbind(
      move_log,
      data.table(
        file = file,
        section = section,
        old_path = old_path,
        new_path = new_path,
        status = if (ok) "moved_unassigned_to_summary" else "unassigned_move_failed"
      )
    )
  }
}

all_section_files <- list.files(diag_dir, recursive = TRUE, full.names = TRUE)
all_section_files <- all_section_files[file.info(all_section_files)$isdir == FALSE]
manifest <- data.table(
  section = basename(dirname(all_section_files)),
  file = basename(all_section_files),
  relative_path = sub(paste0("^", diag_dir, "/?"), "", all_section_files),
  absolute_path = normalizePath(all_section_files, mustWork = FALSE)
)
manifest <- manifest[!file %in% c("README.md", "organization_manifest.tsv")]
setorder(manifest, section, file)
fwrite(manifest, file.path(diag_dir, "organization_manifest.tsv"), sep = "\t")

readme <- c(
  "# Integrated diagnostics folder",
  "",
  "This directory is organized by diagnostic question. Each section folder keeps the relevant plots and source/summary tables together.",
  "",
  "## Sections",
  "",
  paste0("- `", layout$section, "`: ", layout$description),
  "",
  "## Manifest",
  "",
  "- `organization_manifest.tsv` lists every organized file and its section.",
  "- `00_summary_index/plot_manifest.tsv` is the original plot manifest from the diagnostics script.",
  "- `00_summary_index/integrated_clock_diagnostic_summary.md` is the narrative diagnostic summary.",
  "",
  "## Note",
  "",
  "If `04_integrated_clock_diagnostics.R` or `05_plot_integrated_primary_predictions.R` is rerun, new files may be written at the top level again. Rerun `scripts/06_organize_integrated_diagnostics.R` to restore this layout."
)
writeLines(readme, file.path(diag_dir, "README.md"))

cat("Organized integrated diagnostics into section folders under:\n")
cat(diag_dir, "\n")
cat("Files in manifest:", nrow(manifest), "\n")
cat("Move statuses:\n")
print(move_log[, .N, by = status][order(status)])
