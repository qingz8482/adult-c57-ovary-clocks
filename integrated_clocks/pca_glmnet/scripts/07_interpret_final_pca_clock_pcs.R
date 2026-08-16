#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
results_dir <- Sys.getenv("CLOCK_RESULTS_DIR", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = FALSE)
out_dir <- file.path(results_dir, "integrated_diagnostics", "10_pca_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

locked_path <- file.path(results_dir, "locked_model_selection.tsv")
if (!file.exists(locked_path)) stop("Missing locked model selection table: ", locked_path)
locked <- fread(locked_path, sep = "\t")
pca_selector <- locked[model_family == "pca_glmnet", selector][1]
if (is.na(pca_selector) || !nzchar(pca_selector)) stop("Could not determine locked PCA/glmnet selector.")

benchmark_path <- file.path(
  results_dir,
  "pca_glmnet_clock",
  "integrated_bulk_pseudobulk",
  pca_selector,
  "locked_train_validation_to_test_pca_glmnet_model.rds"
)
deployment_path <- file.path(
  results_dir,
  "model_artifacts",
  "deployment",
  paste0("integrated_bulk_pseudobulk__pca_glmnet__", pca_selector, "__deployment_all_control_like.rds")
)
metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
feature_score_path <- file.path(
  results_dir,
  "model_artifacts",
  "deployment",
  paste0("pca_glmnet__", pca_selector),
  "selected_gene_table.tsv"
)
if (!file.exists(feature_score_path)) {
  feature_score_path <- file.path(
    results_dir,
    "pca_glmnet_clock",
    "integrated_bulk_pseudobulk",
    pca_selector,
    "selected_gene_table.tsv"
  )
}
stopifnot(file.exists(benchmark_path), file.exists(deployment_path), file.exists(metadata_path))

meta <- fread(metadata_path)
meta[, age_week := as.numeric(age_week)]
deployment_artifact_for_features <- readRDS(deployment_path)
if (file.exists(feature_score_path)) {
  feature_scores <- fread(feature_score_path, sep = "\t")
} else {
  feature_scores <- data.table(gene_id = deployment_artifact_for_features$selected_genes)
  feature_scores[, `:=`(
    feature_score = NA_real_,
    age_beta = NA_real_,
    feature_score_source = paste0("fallback_selected_genes_for_", pca_selector)
  )]
}
if (!"gene_id" %in% names(feature_scores)) stop("Feature-score table lacks gene_id.")
if (!"feature_score" %in% names(feature_scores)) feature_scores[, feature_score := NA_real_]
if (!"age_beta" %in% names(feature_scores)) feature_scores[, age_beta := NA_real_]
feature_scores[, version_b_age_direction := fifelse(
  is.na(age_beta), NA_character_,
  fifelse(age_beta >= 0, "age-increasing", "age-decreasing")
)]

extract_pca_model <- function(path, artifact_scope) {
  obj <- readRDS(path)
  if (artifact_scope == "benchmark_locked_train_validation") {
    pca <- obj$pca
    genes <- obj$genes
    selected_pc_names <- paste0("PC", seq_len(obj$pc_count))
    cc <- as.matrix(coef(obj$model))
    pc_coef <- data.table(pc = rownames(cc), raw_pc_coefficient = as.numeric(cc[, 1]))
    pc_coef <- pc_coef[pc %in% selected_pc_names]
    calib_slope <- obj$calibration$slope %||% obj$calibration[["slope"]] %||% 1
    pc_coef[, calibrated_pc_coefficient := raw_pc_coefficient * as.numeric(calib_slope)]
    pc_var <- as.data.table(obj$pc_variance_explained)
    if (!"pc" %in% names(pc_var)) pc_var[, pc := seq_len(.N)]
    if (!"variance_explained" %in% names(pc_var)) setnames(pc_var, names(pc_var)[2], "variance_explained")
    reference_meta <- meta[sample_key %in% rownames(pca$x)]
  } else {
    pca <- obj$model$pca_fit
    genes <- obj$selected_genes
    selected_pc_names <- obj$model$selected_pc_names
    pc_coef <- as.data.table(obj$model$pc_coefficient_table)
    pc_var <- as.data.table(obj$model$pc_variance_explained)
    reference_meta <- as.data.table(obj$reference_metadata)
  }

  if (!"sample_key" %in% names(reference_meta)) {
    reference_meta[, sample_key := rownames(reference_meta)]
  }
  reference_meta[, age_week := as.numeric(age_week)]
  score_dt <- as.data.table(pca$x, keep.rownames = "sample_key")
  score_dt <- merge(score_dt, reference_meta[, .(sample_key, age_week, modality, dataset_folder)], by = "sample_key", all.x = TRUE)

  pc_summary <- rbindlist(lapply(selected_pc_names, function(pc_name) {
    pc_num <- as.integer(sub("^PC", "", pc_name))
    score <- score_dt[[pc_name]]
    ok <- is.finite(score) & is.finite(score_dt$age_week)
    data.table(
      artifact_scope = artifact_scope,
      pc = pc_name,
      pc_number = pc_num,
      n_training_samples_for_pca = sum(ok),
      pc_age_pearson = if (sum(ok) >= 3 && uniqueN(score[ok]) > 1 && uniqueN(score_dt$age_week[ok]) > 1) {
        suppressWarnings(cor(score[ok], score_dt$age_week[ok], method = "pearson"))
      } else NA_real_,
      pc_age_spearman = if (sum(ok) >= 3 && uniqueN(score[ok]) > 1 && uniqueN(score_dt$age_week[ok]) > 1) {
        suppressWarnings(cor(score[ok], score_dt$age_week[ok], method = "spearman"))
      } else NA_real_
    )
  }), fill = TRUE)

  pc_var[, pc_number := as.integer(pc)]
  pc_var[, pc := paste0("PC", pc_number)]
  pc_summary <- merge(pc_summary, pc_var[, .(pc, variance_explained)], by = "pc", all.x = TRUE)
  pc_summary <- merge(pc_summary, pc_coef[, .(pc, raw_pc_coefficient, calibrated_pc_coefficient)], by = "pc", all.x = TRUE)
  pc_summary[, cumulative_variance_through_pc := cumsum(variance_explained[order(pc_number)])[rank(pc_number, ties.method = "first")]]
  pc_summary[, pc_score_age_axis := fifelse(pc_age_pearson >= 0, "higher_PC_score_older", "higher_PC_score_younger")]
  pc_summary[, clock_coefficient_direction := fifelse(calibrated_pc_coefficient >= 0, "higher_PC_score_increases_predicted_age", "higher_PC_score_decreases_predicted_age")]
  setorder(pc_summary, pc_number)

  rotation <- as.data.table(pca$rotation[, selected_pc_names, drop = FALSE], keep.rownames = "gene_id")
  gene_pc <- melt(rotation, id.vars = "gene_id", variable.name = "pc", value.name = "pc_loading")
  gene_pc <- merge(gene_pc, pc_summary[, .(artifact_scope, pc, pc_number, pc_age_pearson, calibrated_pc_coefficient)], by = "pc")
  gene_pc <- merge(gene_pc, feature_scores, by = "gene_id", all.x = TRUE)
  gene_pc[, abs_pc_loading := abs(pc_loading)]
  gene_pc[, derived_gene_pc_coefficient := pc_loading * calibrated_pc_coefficient]
  gene_pc[, abs_derived_gene_pc_coefficient := abs(derived_gene_pc_coefficient)]
  gene_pc[, pc_loading_age_direction := fifelse(pc_loading * pc_age_pearson >= 0, "aged-aligned-loading", "young-aligned-loading")]
  gene_pc[, clock_contribution_direction := fifelse(derived_gene_pc_coefficient >= 0, "positive-derived-coefficient", "negative-derived-coefficient")]
  gene_pc[, pc_loading_rank := frank(-abs_pc_loading, ties.method = "first"), by = .(artifact_scope, pc)]
  setorder(gene_pc, artifact_scope, pc_number, pc_loading_rank)

  gene_counts <- gene_pc[, .(
    n_genes_in_pc = .N,
    n_genes_nonzero_loading = sum(pc_loading != 0, na.rm = TRUE),
    n_aged_predictor_by_clock_direction = sum(clock_contribution_direction == "aged-predictor", na.rm = TRUE),
    n_young_predictor_by_clock_direction = sum(clock_contribution_direction == "young-predictor", na.rm = TRUE),
    n_age_increasing_by_version_b_beta = sum(version_b_age_direction == "age-increasing", na.rm = TRUE),
    n_age_decreasing_by_version_b_beta = sum(version_b_age_direction == "age-decreasing", na.rm = TRUE)
  ), by = .(artifact_scope, pc, pc_number)]
  gene_counts[, total_unique_genes_across_selected_pcs := uniqueN(gene_pc$gene_id), by = artifact_scope]

  list(pc_summary = pc_summary, gene_pc = gene_pc, gene_counts = gene_counts)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

benchmark <- extract_pca_model(benchmark_path, "benchmark_locked_train_validation")
deployment <- extract_pca_model(deployment_path, "deployment_all_control_like_external_projection")

pc_summary <- rbindlist(list(benchmark$pc_summary, deployment$pc_summary), fill = TRUE)
gene_pc <- rbindlist(list(benchmark$gene_pc, deployment$gene_pc), fill = TRUE)
gene_counts <- rbindlist(list(benchmark$gene_counts, deployment$gene_counts), fill = TRUE)
top_gene_pc <- gene_pc[pc_loading_rank <= 50]

fwrite(pc_summary, file.path(out_dir, "final_pca_clock_selected_pc_summary.tsv"), sep = "\t")
fwrite(gene_counts, file.path(out_dir, "final_pca_clock_gene_counts_by_pc.tsv"), sep = "\t")
fwrite(gene_pc, file.path(out_dir, "final_pca_clock_gene_pc_loadings_with_feature_scores.tsv"), sep = "\t")
fwrite(top_gene_pc, file.path(out_dir, "final_pca_clock_top50_gene_pc_loadings.tsv"), sep = "\t")

deployment_pc_summary <- pc_summary[artifact_scope == "deployment_all_control_like_external_projection"]

theme_pca <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

p <- ggplot(deployment_pc_summary, aes(x = factor(pc, levels = paste0("PC", seq_len(max(pc_number)))), y = variance_explained)) +
  geom_col(fill = "#4C78A8", color = "grey25", width = 0.68) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * variance_explained)),
    vjust = -0.35, size = 3
  ) +
  labs(
    title = "Final deployment PCA/glmnet selected PCs: variance explained",
    x = "Selected PC",
    y = "Variance explained"
  ) +
  theme_pca()
ggsave(file.path(out_dir, "final_pca_clock_selected_pc_variance.pdf"), p, width = 10, height = 5.8)
ggsave(file.path(out_dir, "final_pca_clock_selected_pc_variance.png"), p, width = 10, height = 5.8, dpi = 300, bg = "white")

p <- ggplot(deployment_pc_summary, aes(x = pc_age_pearson, y = calibrated_pc_coefficient, label = pc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(size = 3, color = "#4C78A8") +
  geom_text(nudge_y = 0.04, size = 3.4, show.legend = FALSE) +
  labs(
    title = "Final deployment PCA/glmnet selected PCs: age correlation and clock coefficient",
    x = "PC score Pearson correlation with age",
    y = "Calibrated PC coefficient"
  ) +
  theme_pca()
ggsave(file.path(out_dir, "final_pca_clock_pc_age_correlation_vs_coefficient.pdf"), p, width = 8, height = 6)
ggsave(file.path(out_dir, "final_pca_clock_pc_age_correlation_vs_coefficient.png"), p, width = 8, height = 6, dpi = 300, bg = "white")

plot_top <- top_gene_pc[artifact_scope == "deployment_all_control_like_external_projection" & pc_loading_rank <= 20]
p <- ggplot(plot_top, aes(x = reorder(gene_id, abs_derived_gene_pc_coefficient), y = abs_derived_gene_pc_coefficient, fill = clock_contribution_direction)) +
  geom_col(color = "grey25") +
  coord_flip() +
  facet_wrap(~pc, scales = "free_y", ncol = 2) +
  scale_fill_manual(
    values = c("positive-derived-coefficient" = "#2F6F9F", "negative-derived-coefficient" = "#B5524A"),
    labels = c(
      "positive-derived-coefficient" = "Positive derived coefficient",
      "negative-derived-coefficient" = "Negative derived coefficient"
    )
  ) +
  labs(
    title = "Deployment PCA/glmnet final clock: top per-PC gene coefficients",
    subtitle = "Derived per-PC gene coefficient = PC loading x calibrated PC coefficient.",
    x = "Gene",
    y = "Absolute derived per-PC coefficient",
    fill = NULL
  ) +
  theme_pca()
ggsave(file.path(out_dir, "final_pca_clock_deployment_top20_genes_per_pc.pdf"), p, width = 12, height = 13)
ggsave(file.path(out_dir, "final_pca_clock_deployment_top20_genes_per_pc.png"), p, width = 12, height = 13, dpi = 300, bg = "white")

cat("Wrote final PCA clock PC interpretation outputs under:\n")
cat(out_dir, "\n")
cat("Deployment selected PCs:", deployment$pc_summary[, .N], "\n")
cat("Benchmark selected PCs:", benchmark$pc_summary[, .N], "\n")
cat("Deployment selected genes:", deployment$gene_counts[, unique(total_unique_genes_across_selected_pcs)], "\n")
