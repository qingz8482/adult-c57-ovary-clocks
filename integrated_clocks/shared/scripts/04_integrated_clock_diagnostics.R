#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(glmnet)
})

options(stringsAsFactors = FALSE)

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
out_dir <- file.path(results_dir, "integrated_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

expected_scope <- "integrated_bulk_pseudobulk"
retained <- data.table(model_family = character(), locked_selector = character())

warnings_log <- character()
plot_manifest <- data.table(prefix = character(), png_path = character(), pdf_path = character())

add_warning <- function(msg) {
  warnings_log <<- c(warnings_log, msg)
  message("WARNING: ", msg)
}

detect_col <- function(dt, candidates, context, required = TRUE) {
  nms <- names(dt)
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[[1]])
  hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit) > 0) return(hit[[1]])
  msg <- paste0("Could not detect ", context, " column. Tried: ", paste(candidates, collapse = ", "))
  if (required) stop(msg, call. = FALSE)
  add_warning(msg)
  NA_character_
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

selector_label <- function(x) {
  fifelse(
    x == "hvg_baseline", "HVG",
    fifelse(x == "version_b_age_dataset_penalized", "Version B", as.character(x))
  )
}

family_label <- function(x) {
  fifelse(x == "elasticnet", "ElasticNet", fifelse(x == "pca_glmnet", "PCA/glmnet", as.character(x)))
}

split_label <- function(x) {
  fifelse(x == "validation", "Validation", fifelse(x == "test", "Test", fifelse(x == "train", "Train", as.character(x))))
}

bio_bin_from_age <- function(age_week) {
  fifelse(
    is.na(age_week), NA_character_,
    fifelse(
      age_week >= 7 & age_week < 16, "young_adult",
      fifelse(
        age_week >= 16 & age_week < 40, "mature_adult",
        fifelse(age_week >= 40 & age_week < 65, "advanced_maternal_age",
          fifelse(age_week >= 65, "post_reproductive_ageing", NA_character_)
        )
      )
    )
  )
}

recode_bio_age_bin <- function(x) {
  fifelse(
    x == "mid_adult", "mature_adult",
    fifelse(
      x == "aged", "advanced_maternal_age",
      fifelse(x == "old_aged", "post_reproductive_ageing", x)
    )
  )
}

split_stratum_from_age <- function(age_week) {
  fifelse(
    is.na(age_week), NA_character_,
    fifelse(
      age_week >= 7 & age_week < 13, "early_adult",
      fifelse(age_week >= 13 & age_week < 39, "reproductive_adult", fifelse(age_week >= 39, "aged_adult", NA_character_))
    )
  )
}

corr_or_na <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || uniqueN(x[ok]) < 2 || uniqueN(y[ok]) < 2) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

metric_dt <- function(dt, by_cols = character()) {
  out <- dt[, .(
    n = .N,
    mae_week = mean(abs_error_week, na.rm = TRUE),
    rmse_week = sqrt(mean(residual_week^2, na.rm = TRUE)),
    median_abs_error_week = stats::median(abs_error_week, na.rm = TRUE),
    mean_residual_week = mean(residual_week, na.rm = TRUE),
    median_residual_week = stats::median(residual_week, na.rm = TRUE),
    residual_sd_week = stats::sd(residual_week, na.rm = TRUE),
    pearson = corr_or_na(age_week, predicted_age_week, "pearson"),
    spearman = corr_or_na(age_week, predicted_age_week, "spearman")
  ), by = by_cols]
  out[]
}

save_plot <- function(p, prefix, width = 10, height = 7) {
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, p, width = width, height = height, bg = "white", limitsize = FALSE)
  plot_manifest <<- rbind(
    plot_manifest,
    data.table(prefix = prefix, png_path = png_path, pdf_path = pdf_path)
  )
}

theme_diag <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

pred <- read_tsv(file.path(results_dir, "all_clock_predictions.tsv"))
perf <- read_tsv(file.path(results_dir, "all_clock_performance.tsv"))
locked <- read_tsv(file.path(results_dir, "locked_model_selection.tsv"))
retained <- unique(locked[, .(model_family, locked_selector = selector)])
retained <- retained[model_family %in% c("elasticnet", "pca_glmnet")]
setorder(retained, model_family)
retained[, clock_label := paste0(family_label(model_family), "\n", selector_label(locked_selector), " selector")]
locked_map <- retained[, .(model_family, selector = locked_selector, clock_label)]
locked_keys <- paste(locked_map$model_family, locked_map$selector, sep = "\r")

sample_col <- detect_col(pred, c("sample_key", "sample_id", "sample"), "sample key")
scope_col <- detect_col(pred, c("training_scope", "model_scope", "scope"), "training scope")
family_col <- detect_col(pred, c("model_family", "model_type", "clock_model"), "model family")
selector_col <- detect_col(pred, c("selector", "model_version", "feature_selector"), "selector")
split_col <- detect_col(pred, c("split_role", "split", "role"), "split role")
modality_col <- detect_col(pred, c("modality", "sample_modality"), "modality")
dataset_col <- detect_col(pred, c("dataset_folder", "dataset"), "dataset folder")
age_col <- detect_col(pred, c("age_week", "observed_age_week", "observed_age"), "age")
pred_col <- detect_col(pred, c("predicted_age_week", "pred_age_week", "predicted_age"), "prediction")
fit_stage_col <- detect_col(pred, c("fit_stage", "stage"), "fit stage", required = FALSE)

setnames(pred, c(sample_col, scope_col, family_col, selector_col, split_col, modality_col, dataset_col, age_col, pred_col),
  c("sample_key", "training_scope", "model_family", "selector", "split_role", "modality", "dataset_folder", "age_week", "predicted_age_week"),
  skip_absent = TRUE
)
if (!is.na(fit_stage_col) && fit_stage_col != "fit_stage") setnames(pred, fit_stage_col, "fit_stage")
pred[, age_week := safe_num(age_week)]
pred[, predicted_age_week := safe_num(predicted_age_week)]
if (!"residual_week" %in% names(pred)) pred[, residual_week := predicted_age_week - age_week]
if (!"abs_error_week" %in% names(pred)) pred[, abs_error_week := abs(residual_week)]
if (!"biological_age_bin" %in% names(pred)) pred[, biological_age_bin := bio_bin_from_age(age_week)]
pred[, biological_age_bin := recode_bio_age_bin(biological_age_bin)]
if ("bio_age_bin" %in% names(pred)) pred[, bio_age_bin := recode_bio_age_bin(bio_age_bin)]
if (!"split_age_stratum" %in% names(pred)) pred[, split_age_stratum := split_stratum_from_age(age_week)]
pred[, residual_week := safe_num(residual_week)]
pred[, abs_error_week := safe_num(abs_error_week)]
pred[, model_family_label := family_label(model_family)]
pred[, selector_label := selector_label(selector)]
pred[, split_label := split_label(split_role)]

essential_missing <- !complete.cases(pred[, .(sample_key, training_scope, model_family, selector, split_role, modality, dataset_folder, age_week, predicted_age_week)])
if (any(essential_missing)) {
  add_warning(paste0("Dropped ", sum(essential_missing), " prediction rows with missing essential values."))
  pred <- pred[!essential_missing]
}

if (uniqueN(pred$training_scope) != 1L || unique(pred$training_scope) != expected_scope) {
  stop("Expected exactly one training_scope named ", expected_scope, call. = FALSE)
}

selected_gene_files <- list.files(results_dir, pattern = "selected.*genes.*\\.tsv$", recursive = TRUE, full.names = TRUE)
bad_selected <- selected_gene_files[grepl("bulk_rnaseq|scrna_pseudobulk|pseudobulk_only|bulk_only", basename(selected_gene_files))]
if (length(bad_selected) > 0) {
  stop("Selected-gene files appear modality-specific: ", paste(bad_selected, collapse = "; "), call. = FALSE)
}

all_samples_n <- uniqueN(pred$sample_key)
prediction_audit <- pred[, .(
  n_predictions = .N,
  n_unique_samples = uniqueN(sample_key),
  n_bulk = uniqueN(sample_key[modality == "bulk_rnaseq"]),
  n_pseudobulk = uniqueN(sample_key[modality == "scrna_pseudobulk"]),
  expected_n_if_all_117_samples = 117L
), by = .(model_family, selector, split_role)]
prediction_audit[, missing_or_extra_flag := fifelse(n_unique_samples == 117L, "all_117_present", "not_all_117_expected_for_split_or_locked_stage")]
fwrite(prediction_audit[order(model_family, selector, split_role)], file.path(out_dir, "prediction_row_audit.tsv"), sep = "\t")

locked_expected <- merge(retained, locked, by.x = c("model_family", "locked_selector"), by.y = c("model_family", "selector"), all.x = TRUE)
if (any(is.na(locked_expected$mae_week))) {
  add_warning("Locked selection table did not fully match the expected family-specific locked selectors.")
}

locked_eval <- merge(pred, locked_map, by = c("model_family", "selector"))
locked_eval <- locked_eval[fit_stage %in% c("train_to_validation", "locked_train_validation_to_test")]
locked_eval <- locked_eval[split_role %in% c("validation", "test")]
locked_eval[, split_role := factor(split_role, levels = c("validation", "test"))]
locked_eval[, model_family := factor(model_family, levels = c("elasticnet", "pca_glmnet"))]

locked_metrics <- metric_dt(locked_eval, c("training_scope", "model_family", "selector", "split_role"))
locked_wide <- dcast(
  locked_metrics,
  training_scope + model_family + selector ~ split_role,
  value.var = "mae_week"
)
setnames(locked_wide, c("validation", "test"), c("validation_mae_week", "test_mae_week"), skip_absent = TRUE)
locked_wide[, validation_to_test_mae_gap := test_mae_week - validation_mae_week]
locked_wide[, retained_clock_statement := "Both ElasticNet and PCA/glmnet are retained as family-specific integrated clocks."]
setnames(locked_wide, "selector", "locked_selector")
fwrite(locked_wide[order(model_family)], file.path(out_dir, "family_specific_locked_clock_summary.tsv"), sep = "\t")

# A. Family-specific locked clock summary.
summary_long <- melt(
  locked_wide,
  id.vars = c("training_scope", "model_family", "locked_selector"),
  measure.vars = c("validation_mae_week", "test_mae_week"),
  variable.name = "metric",
  value.name = "mae_week"
)
summary_long[, split := fifelse(metric == "validation_mae_week", "Validation", "Test")]
summary_long[, model_family_label := family_label(model_family)]
summary_long[, selector_lab := selector_label(locked_selector)]
p <- ggplot(summary_long, aes(x = model_family_label, y = mae_week, fill = split)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64, color = "grey25") +
  geom_text(aes(label = sprintf("%.2f", mae_week)), position = position_dodge(width = 0.72), vjust = -0.35, size = 4) +
  geom_text(data = unique(summary_long[, .(model_family_label, selector_lab)]),
    aes(x = model_family_label, y = 0, label = paste("Locked:", selector_lab)),
    inherit.aes = FALSE, vjust = 1.25, size = 3.7
  ) +
  labs(
    title = "Integrated clock benchmark: both model families retained",
    x = "Model family",
    y = "MAE (weeks)",
    fill = "Split"
  ) +
  expand_limits(y = max(summary_long$mae_week, na.rm = TRUE) * 1.18) +
  theme_diag()
save_plot(p, "01_family_specific_locked_clocks", 9, 6)

# B. Candidate selector validation/test MAE.
candidate_eval <- pred[split_role %in% c("validation", "test") & fit_stage %in% c("train_to_validation", "locked_train_validation_to_test")]
candidate_metrics <- metric_dt(candidate_eval, c("model_family", "selector", "split_role"))
candidate_metrics[, locked := paste(model_family, selector, sep = "\r") %in% locked_keys]
candidate_metrics[, candidate_label := paste(family_label(model_family), selector_label(selector), sep = " | ")]
p <- ggplot(candidate_metrics, aes(x = candidate_label, y = mae_week, fill = split_label(split_role), alpha = locked)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.67, color = "grey25") +
  geom_text(aes(label = sprintf("%.2f", mae_week)), position = position_dodge(width = 0.75), vjust = -0.35, size = 3.4) +
  facet_wrap(~family_label(model_family), scales = "free_x") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.55), guide = "none") +
  labs(
    title = "Candidate selector validation/test MAE",
    subtitle = "Locked selectors are shown at full opacity; test MAE is available for locked configurations.",
    x = NULL,
    y = "MAE (weeks)",
    fill = "Split"
  ) +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_plot(p, "02_selector_validation_test_mae", 12, 6.5)

annot_metrics <- metric_dt(locked_eval, c("model_family", "split_role"))
lims <- range(c(locked_eval$age_week, locked_eval$predicted_age_week), na.rm = TRUE)
pad <- diff(lims) * 0.06
if (!is.finite(pad) || pad == 0) pad <- 1
lims <- c(lims[1] - pad, lims[2] + pad)
annot_metrics[, label := sprintf("n=%d\nMAE=%.2f\nRMSE=%.2f\nPearson=%s\nSpearman=%s",
  n, mae_week, rmse_week,
  fifelse(is.na(pearson), "NA", sprintf("%.2f", pearson)),
  fifelse(is.na(spearman), "NA", sprintf("%.2f", spearman))
)]
annot_metrics[, x_pos := lims[1] + diff(lims) * 0.04]
annot_metrics[, y_pos := lims[2] - diff(lims) * 0.04]

# C. Predicted vs observed.
p <- ggplot(locked_eval, aes(x = age_week, y = predicted_age_week, color = modality, shape = dataset_folder)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35") +
  geom_point(size = 2.6, alpha = 0.9) +
  geom_label(data = annot_metrics, aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.2, fill = "white", linewidth = 0.2
  ) +
  facet_grid(split_label(split_role) ~ family_label(model_family)) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(title = "Locked integrated clocks: predicted vs observed age", x = "Observed age (weeks)", y = "Predicted age (weeks)", color = "Modality", shape = "Dataset") +
  theme_diag() +
  guides(shape = guide_legend(ncol = 1))
save_plot(p, "03_locked_predicted_vs_observed", 13, 8)

primary_eval <- merge(pred, locked_map, by = c("model_family", "selector"))
primary_eval <- primary_eval[
  (split_role == "train" & fit_stage == "train_fit") |
    (split_role == "validation" & fit_stage == "train_to_validation") |
    (split_role == "test" & fit_stage == "locked_train_validation_to_test")
]
if (nrow(primary_eval) > 0) {
  primary_eval[, split_label := fifelse(
    split_role == "train", "training model",
    fifelse(split_role == "validation", "validation", "held-out test")
  )]
  primary_eval[, split_label := factor(split_label, levels = c("training model", "validation", "held-out test"))]
  primary_eval[, clock_label := factor(clock_label, levels = locked_map$clock_label)]
  primary_eval[, dataset_label := gsub("_", " ", dataset_folder)]
  primary_eval[, dataset_label := gsub("2024 NatureAging Isola aging mouse ovary atlas", "Isola ageing atlas", dataset_label, fixed = TRUE)]
  primary_eval[, dataset_label := gsub("2025 GeroScience Zhang multiomic aging mouse ovary", "Zhang multiomic ageing", dataset_label, fixed = TRUE)]
  primary_eval[, dataset_label := gsub("2026 CommBio Galligos ovarian inflammaging", "Galligos inflammaging", dataset_label, fixed = TRUE)]
  primary_eval[, dataset_label := gsub("CRA003645", "CRA003645 bulk timecourse", dataset_label, fixed = TRUE)]
  primary_eval[, dataset_label := gsub("GSE252791", "GSE252791 bulk validation", dataset_label, fixed = TRUE)]
  primary_eval[, biological_age_bin := factor(
    biological_age_bin,
    levels = c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
  )]
  primary_eval[, split_age_stratum := factor(split_stratum_from_age(age_week), levels = c("early_adult", "reproductive_adult", "aged_adult"))]
}

# D. Residual vs age.
residual_plot_df <- copy(primary_eval)
if (nrow(residual_plot_df) == 0) {
  add_warning("Skipped residual vs age training/validation/test plot: no frozen integrated predictions were found.")
} else {
  dataset_levels <- sort(unique(residual_plot_df$dataset_label))
  dataset_palette <- grDevices::hcl.colors(length(dataset_levels), palette = "Dark 3")
  names(dataset_palette) <- dataset_levels

  residual_stats <- residual_plot_df[, .(
    n = .N,
    mae_week = mean(abs_error_week, na.rm = TRUE),
    mean_residual_week = mean(residual_week, na.rm = TRUE),
    residual_age_pearson = corr_or_na(age_week, residual_week, "pearson")
  ), by = .(clock_label, split_label)]
  residual_stats[, label := sprintf("n=%d\nMAE=%.1f\nmean resid.=%.1f\nr(age,resid.)=%s",
    n, mae_week, mean_residual_week,
    fifelse(is.na(residual_age_pearson), "NA", sprintf("%.2f", residual_age_pearson))
  )]

  age_lims <- range(residual_plot_df$age_week, na.rm = TRUE)
  age_pad <- diff(age_lims) * 0.06
  if (!is.finite(age_pad) || age_pad == 0) age_pad <- 1
  age_lims <- c(floor(age_lims[1] - age_pad), ceiling(age_lims[2] + age_pad))
  residual_lims <- range(residual_plot_df$residual_week, na.rm = TRUE)
  residual_pad <- diff(residual_lims) * 0.12
  if (!is.finite(residual_pad) || residual_pad == 0) residual_pad <- 1
  residual_lims <- c(floor(residual_lims[1] - residual_pad), ceiling(residual_lims[2] + residual_pad))

  p <- ggplot(residual_plot_df, aes(x = age_week, y = residual_week, color = dataset_label, shape = biological_age_bin)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(size = 2.8, alpha = 0.9) +
    geom_text(
      data = residual_stats,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.05,
      vjust = 1.15,
      size = 3
    ) +
    facet_grid(clock_label ~ split_label) +
    coord_cartesian(xlim = age_lims, ylim = residual_lims) +
    scale_color_manual(values = dataset_palette, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
    labs(
      title = "Integrated bulk+pseudobulk clocks: residuals for training fit, validation, and primary held-out test",
      x = "Observed age (weeks)",
      y = "Residual (predicted - observed, weeks)",
      color = "Dataset",
      shape = "Biological age bin"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  smooth_data <- residual_plot_df[, if (.N >= 6) .SD, by = .(clock_label, split_label)]
  if (nrow(smooth_data) > 0) {
    p <- p + geom_smooth(
      data = smooth_data,
      aes(x = age_week, y = residual_week, group = interaction(clock_label, split_label)),
      inherit.aes = FALSE,
      method = "loess",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.55
    )
  } else {
    add_warning("Skipped residual smooth: no residual model/split facet had at least six samples.")
  }
  save_plot(p, "04_locked_residual_vs_age", 15, 7)
  fwrite(residual_plot_df, file.path(out_dir, "04_locked_residual_vs_age_source.tsv"), sep = "\t")
  fwrite(residual_stats, file.path(out_dir, "04_locked_residual_vs_age_fit_stats.tsv"), sep = "\t")
}

# E-H. Stratified errors.
residual_bias_by_dataset <- metric_dt(locked_eval, c("model_family", "split_role", "dataset_folder"))
fwrite(residual_bias_by_dataset[order(model_family, split_role, dataset_folder)], file.path(out_dir, "residual_bias_by_dataset.tsv"), sep = "\t")
p <- ggplot(locked_eval, aes(x = reorder(dataset_folder, residual_week, FUN = median), y = residual_week, color = modality)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 2.2) +
  coord_flip() +
  facet_grid(family_label(model_family) ~ split_label(split_role), scales = "free_y", space = "free_y") +
  labs(title = "Locked residual bias by dataset", x = "Dataset", y = "Residual (weeks)", color = "Modality") +
  theme_diag()
save_plot(p, "05_locked_residual_by_dataset", 13, 9)
p <- ggplot(locked_eval, aes(x = reorder(dataset_folder, abs_error_week, FUN = median), y = abs_error_week, color = modality)) +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.75, size = 2.2) +
  coord_flip() +
  facet_grid(family_label(model_family) ~ split_label(split_role), scales = "free_y", space = "free_y") +
  labs(title = "Locked absolute error by dataset", x = "Dataset", y = "Absolute error (weeks)", color = "Modality") +
  theme_diag()
save_plot(p, "05_locked_abs_error_by_dataset", 13, 9)

performance_by_modality <- metric_dt(primary_eval, c("model_family", "split_role", "modality"))
fwrite(performance_by_modality[order(model_family, split_role, modality)], file.path(out_dir, "performance_by_modality.tsv"), sep = "\t")
p <- ggplot(primary_eval, aes(x = modality, y = abs_error_week, color = modality)) +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 2.3) +
  facet_grid(clock_label ~ split_label) +
  labs(title = "Locked absolute error by modality", x = "Modality", y = "Absolute error (weeks)", color = "Modality") +
  theme_diag() +
  theme(legend.position = "right")
save_plot(p, "06_locked_error_by_modality", 15, 7)

bio_levels <- c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
locked_eval[, biological_age_bin := factor(bio_bin_from_age(age_week), levels = bio_levels)]
performance_by_bio <- metric_dt(primary_eval, c("model_family", "split_role", "biological_age_bin"))
fwrite(performance_by_bio[order(model_family, split_role, biological_age_bin)], file.path(out_dir, "performance_by_biological_age_bin.tsv"), sep = "\t")
p <- ggplot(primary_eval, aes(x = biological_age_bin, y = abs_error_week, color = biological_age_bin)) +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 2.3) +
  facet_grid(clock_label ~ split_label) +
  labs(title = "Locked absolute error by biological age bin", x = "Biological age bin", y = "Absolute error (weeks)", color = "Age bin") +
  theme_diag() +
  theme(legend.position = "right")
save_plot(p, "07_locked_error_by_biological_age_bin", 15, 7)

stratum_levels <- c("early_adult", "reproductive_adult", "aged_adult")
locked_eval[, split_age_stratum := factor(split_stratum_from_age(age_week), levels = stratum_levels)]
performance_by_stratum <- metric_dt(primary_eval, c("model_family", "split_role", "split_age_stratum"))
fwrite(performance_by_stratum[order(model_family, split_role, split_age_stratum)], file.path(out_dir, "performance_by_split_age_stratum.tsv"), sep = "\t")
p <- ggplot(primary_eval, aes(x = split_age_stratum, y = abs_error_week, color = split_age_stratum)) +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.8, size = 2.3) +
  facet_grid(clock_label ~ split_label) +
  labs(title = "Locked absolute error by balanced split age stratum", x = "Split age stratum", y = "Absolute error (weeks)", color = "Split age stratum") +
  theme_diag() +
  theme(legend.position = "right")
save_plot(p, "08_locked_error_by_split_age_stratum", 15, 7)

# I. Naive baseline comparison.
train_rows <- pred[split_role == "train" & fit_stage == "train_fit"]
if (nrow(train_rows) == 0) {
  add_warning("Training rows unavailable; skipped naive baseline comparison.")
  fwrite(data.table(warning = "Training rows unavailable."), file.path(out_dir, "naive_baseline_performance.tsv"), sep = "\t")
} else {
  train_age <- unique(train_rows[, .(sample_key, modality, age_week)])
  mean_age <- mean(train_age$age_week, na.rm = TRUE)
  median_age <- stats::median(train_age$age_week, na.rm = TRUE)
  baseline_rows <- unique(locked_eval[, .(sample_key, modality, dataset_folder, split_role, age_week)])
  baseline_pred <- rbindlist(list(
    copy(baseline_rows)[, `:=`(method = "Train mean baseline", predicted_age_week = mean_age)],
    copy(baseline_rows)[, `:=`(method = "Train median baseline", predicted_age_week = median_age)]
  ))
  baseline_pred[, `:=`(
    residual_week = predicted_age_week - age_week,
    abs_error_week = abs(predicted_age_week - age_week)
  )]
  model_perf <- metric_dt(locked_eval, c("model_family", "split_role"))
  model_perf[, method := paste(family_label(model_family), "locked clock")]
  model_perf[, group := "overall"]
  model_perf_mod <- metric_dt(locked_eval, c("model_family", "split_role", "modality"))
  model_perf_mod[, method := paste(family_label(model_family), "locked clock")]
  model_perf_mod[, group := modality]
  baseline_perf <- metric_dt(baseline_pred, c("method", "split_role"))
  baseline_perf[, `:=`(model_family = NA_character_, group = "overall")]
  baseline_perf_mod <- metric_dt(baseline_pred, c("method", "split_role", "modality"))
  baseline_perf_mod[, `:=`(model_family = NA_character_, group = modality)]
  naive_perf <- rbindlist(list(model_perf, model_perf_mod, baseline_perf, baseline_perf_mod), fill = TRUE)
  fwrite(naive_perf, file.path(out_dir, "naive_baseline_performance.tsv"), sep = "\t")
  p <- ggplot(naive_perf[group == "overall"], aes(x = method, y = mae_week, fill = split_label(split_role))) +
    geom_col(position = position_dodge(width = 0.75), color = "grey25") +
    coord_flip() +
    labs(title = "Locked clocks vs naive age baselines", x = NULL, y = "MAE (weeks)", fill = "Split") +
    theme_diag()
  save_plot(p, "09_naive_baseline_comparison", 11, 6.5)
}

# J. Paired retained clock comparison.
paired <- dcast(
  primary_eval,
  sample_key + modality + dataset_folder + dataset_label + split_role + split_label + age_week + biological_age_bin + split_age_stratum ~ model_family,
  value.var = "abs_error_week"
)
setnames(paired, c("elasticnet", "pca_glmnet"), c("abs_error_elasticnet", "abs_error_pca_glmnet"), skip_absent = TRUE)
paired <- paired[!is.na(abs_error_elasticnet) & !is.na(abs_error_pca_glmnet)]
paired[, delta_abs_error := abs_error_pca_glmnet - abs_error_elasticnet]
fwrite(paired, file.path(out_dir, "paired_retained_clock_sample_comparison.tsv"), sep = "\t")
paired_summary <- paired[, {
  p_value <- NA_real_
  if (.N >= 5) p_value <- suppressWarnings(stats::wilcox.test(delta_abs_error, mu = 0)$p.value)
  .(
    n = .N,
    median_delta_abs_error = stats::median(delta_abs_error, na.rm = TRUE),
    mean_delta_abs_error = mean(delta_abs_error, na.rm = TRUE),
    wilcoxon_p_value_descriptive = p_value,
    retained_clock_statement = "Both clocks are retained; paired deltas are diagnostic only."
  )
}, by = .(split_role)]
fwrite(paired_summary, file.path(out_dir, "paired_retained_clock_summary.tsv"), sep = "\t")

paired_delta_summary <- function(dt, by_cols) {
  dt[, {
    p_value <- NA_real_
    if (.N >= 5) p_value <- suppressWarnings(stats::wilcox.test(delta_abs_error, mu = 0)$p.value)
    .(
      n = .N,
      median_delta_abs_error = stats::median(delta_abs_error, na.rm = TRUE),
      mean_delta_abs_error = mean(delta_abs_error, na.rm = TRUE),
      wilcoxon_p_value_descriptive = p_value
    )
  }, by = by_cols]
}
paired_summary_by_dataset <- paired_delta_summary(paired, c("split_role", "dataset_label", "modality"))
paired_summary_by_modality <- paired_delta_summary(paired, c("split_role", "modality"))
paired_summary_by_age_bin <- paired_delta_summary(paired, c("split_role", "biological_age_bin"))
fwrite(paired_summary_by_dataset[order(split_role, dataset_label, modality)], file.path(out_dir, "paired_retained_clock_summary_by_dataset.tsv"), sep = "\t")
fwrite(paired_summary_by_modality[order(split_role, modality)], file.path(out_dir, "paired_retained_clock_summary_by_modality.tsv"), sep = "\t")
fwrite(paired_summary_by_age_bin[order(split_role, biological_age_bin)], file.path(out_dir, "paired_retained_clock_summary_by_age_bin.tsv"), sep = "\t")

paired[, split_label := factor(split_label, levels = c("training model", "validation", "held-out test"))]
delta_subtitle <- "Delta = PCA/glmnet absolute error - ElasticNet absolute error; positive values favor ElasticNet."
p <- ggplot(paired, aes(x = dataset_label, y = delta_abs_error, color = modality)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.85, size = 2.4) +
  coord_flip() +
  facet_grid(. ~ split_label, scales = "free_y") +
  labs(
    title = "Paired retained-clock delta absolute error by dataset",
    subtitle = delta_subtitle,
    x = "Dataset",
    y = "Delta absolute error (weeks)",
    color = "Modality"
  ) +
  theme_diag() +
  theme(legend.position = "right")
save_plot(p, "10_paired_retained_clock_comparison_by_dataset", 15, 7)
p <- ggplot(paired, aes(x = modality, y = delta_abs_error, color = modality)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.85, size = 2.4) +
  facet_grid(. ~ split_label) +
  labs(
    title = "Paired retained-clock delta absolute error by modality",
    subtitle = delta_subtitle,
    x = "Modality",
    y = "Delta absolute error (weeks)",
    color = "Modality"
  ) +
  theme_diag() +
  theme(legend.position = "right")
save_plot(p, "10_paired_retained_clock_comparison_by_modality", 15, 7)
p <- ggplot(paired, aes(x = biological_age_bin, y = delta_abs_error, color = biological_age_bin)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, color = "grey30") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.85, size = 2.4) +
  facet_grid(. ~ split_label) +
  labs(
    title = "Paired retained-clock delta absolute error by age bin",
    subtitle = delta_subtitle,
    x = "Biological age bin",
    y = "Delta absolute error (weeks)",
    color = "Age bin"
  ) +
  theme_diag() +
  theme(legend.position = "right")
save_plot(p, "10_paired_retained_clock_comparison_by_age_bin", 15, 7)

# K. Dynamic range and calibration.
calib <- primary_eval[, {
  fit <- tryCatch(lm(predicted_age_week ~ age_week), error = function(e) NULL)
  .(
    n = .N,
    observed_age_min = min(age_week, na.rm = TRUE),
    observed_age_max = max(age_week, na.rm = TRUE),
    observed_age_range = diff(range(age_week, na.rm = TRUE)),
    predicted_age_min = min(predicted_age_week, na.rm = TRUE),
    predicted_age_max = max(predicted_age_week, na.rm = TRUE),
    predicted_age_range = diff(range(predicted_age_week, na.rm = TRUE)),
    calibration_intercept = if (is.null(fit)) NA_real_ else unname(coef(fit)[1]),
    calibration_slope = if (is.null(fit)) NA_real_ else unname(coef(fit)[2]),
    residual_correlation_with_observed_age = corr_or_na(age_week, residual_week, "pearson")
  )
}, by = .(model_family, selector, split_role)]
fwrite(calib, file.path(out_dir, "calibration_dynamic_range_summary.tsv"), sep = "\t")
dist_long <- rbindlist(list(
  primary_eval[, .(clock_label, split_label, age_type = "Observed", age_week = age_week)],
  primary_eval[, .(clock_label, split_label, age_type = "Predicted", age_week = predicted_age_week)]
))
p <- ggplot(dist_long, aes(x = age_type, y = age_week, fill = age_type)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "grey25") +
  geom_jitter(width = 0.12, alpha = 0.65, size = 1.8) +
  facet_grid(clock_label ~ split_label) +
  labs(title = "Dynamic range and calibration: observed vs predicted distributions", x = NULL, y = "Age (weeks)", fill = "Age type") +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 0), legend.position = "right")
save_plot(p, "11_dynamic_range_calibration", 15, 7)

# L. Integrated selected-gene diagnostics.
selected_summary <- data.table()
top30_all <- data.table()
feature_paths <- list.files(results_dir, pattern = "feature_score_table\\.tsv$", recursive = TRUE, full.names = TRUE)
selected_paths <- selected_gene_files
for (path in selected_paths) {
  dt <- tryCatch(fread(path, sep = "\t"), error = function(e) NULL)
  if (is.null(dt)) next
  rel <- sub(paste0("^", results_dir, "/?"), "", path)
  family <- fifelse(grepl("elasticnet", rel), "elasticnet", fifelse(grepl("pca_glmnet", rel), "pca_glmnet", NA_character_))
  selector <- fifelse(grepl("versionB|version_b", rel), "version_b_age_dataset_penalized", fifelse(grepl("hvg", rel), "hvg_baseline", NA_character_))
  selected_summary <- rbind(selected_summary, data.table(
    training_scope = expected_scope,
    model_family = family,
    selector = selector,
    file_path = path,
    n_genes = nrow(dt),
    has_feature_score = "feature_score" %in% names(dt)
  ), fill = TRUE)
}
for (path in feature_paths) {
  dt <- tryCatch(fread(path, sep = "\t"), error = function(e) NULL)
  if (is.null(dt) || !"gene_id" %in% names(dt) || !"feature_score" %in% names(dt)) next
  rel <- sub(paste0("^", results_dir, "/?"), "", path)
  family <- fifelse(grepl("elasticnet", rel), "elasticnet", fifelse(grepl("pca_glmnet", rel), "pca_glmnet", NA_character_))
  top <- head(dt[order(-feature_score)], 30)
  top[, `:=`(training_scope = expected_scope, model_family = family, selector = "version_b_age_dataset_penalized", source_file = path)]
  top30_all <- rbind(top30_all, top, fill = TRUE)
  p <- ggplot(top[order(feature_score)], aes(x = reorder(gene_id, feature_score), y = feature_score)) +
    geom_col(fill = "#4C78A8", color = "grey25") +
    coord_flip() +
    labs(title = paste("Top 30 integrated Version B genes:", family_label(family)), x = "Gene", y = "Feature score") +
    theme_diag()
  save_plot(p, paste0("12_top30_versionB_genes_integrated_", family), 9, 8)
  batch_col <- if ("partial_r2_source_batch" %in% names(dt)) "partial_r2_source_batch" else if ("partial_r2_dataset" %in% names(dt)) "partial_r2_dataset" else NA_character_
  if (!is.na(batch_col) && "partial_r2_age" %in% names(dt)) {
    label_dt <- head(dt[order(-feature_score)], 12)
    p <- ggplot(dt, aes(x = partial_r2_age, y = .data[[batch_col]])) +
      geom_point(alpha = 0.35, size = 1.2, color = "#4C78A8") +
      geom_text(data = label_dt, aes(label = gene_id), size = 3, check_overlap = TRUE) +
      labs(title = paste("Version B age vs source-batch partial R2:", family_label(family)), x = "Partial R2 age", y = batch_col) +
      theme_diag()
    save_plot(p, paste0("12_versionB_partial_r2_integrated_", family), 8, 6)
  }
}
fwrite(selected_summary[order(model_family, selector, file_path)], file.path(out_dir, "selected_gene_summary.tsv"), sep = "\t")
fwrite(top30_all, file.path(out_dir, "top30_versionB_genes.tsv"), sep = "\t")

# M. PCA diagnostics.
pca_var_paths <- list.files(results_dir, pattern = "pca_variance\\.tsv$", recursive = TRUE, full.names = TRUE)
pca_model_paths <- list.files(results_dir, pattern = "pca_glmnet_model\\.rds$", recursive = TRUE, full.names = TRUE)
if (length(pca_var_paths) == 0 && length(pca_model_paths) == 0) {
  add_warning("No PCA diagnostic files found; skipped PCA diagnostics.")
  fwrite(data.table(warning = "No PCA diagnostic files found."), file.path(out_dir, "pca_diagnostics_warnings.tsv"), sep = "\t")
} else {
  pca_var <- rbindlist(lapply(pca_var_paths, function(path) {
    dt <- fread(path, sep = "\t")
    rel <- sub(paste0("^", results_dir, "/?"), "", path)
    dt[, `:=`(
      source_file = path,
      selector = fifelse(grepl("version_b_age_dataset_penalized", rel), "version_b_age_dataset_penalized", "hvg_baseline"),
      fit_stage = fifelse(grepl("locked_train_validation_to_test", basename(path)), "locked_train_validation_to_test",
        fifelse(grepl("train_to_validation", basename(path)), "train_to_validation", "train_fit")
      )
    )]
    dt
  }), fill = TRUE)
  fwrite(pca_var, file.path(out_dir, "pca_variance_summary.tsv"), sep = "\t")
  p <- ggplot(pca_var[pc <= 20], aes(x = pc, y = variance_explained, color = selector)) +
    geom_line() +
    geom_point() +
    facet_wrap(~fit_stage) +
    labs(title = "PCA diagnostics: variance explained", x = "PC", y = "Variance explained", color = "Selector") +
    theme_diag()
  save_plot(p, "13_pca_diagnostics_variance_explained", 10, 6)

  coef_dt <- data.table()
  score_dt <- data.table()
  for (path in pca_model_paths) {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
    if (is.null(obj)) next
    rel <- sub(paste0("^", results_dir, "/?"), "", path)
    selector <- fifelse(grepl("version_b_age_dataset_penalized", rel), "version_b_age_dataset_penalized", "hvg_baseline")
    fit_stage <- fifelse(grepl("locked_train_validation_to_test", basename(path)), "locked_train_validation_to_test",
      fifelse(grepl("train_to_validation", basename(path)), "train_to_validation", "train_fit")
    )
    if (!is.null(obj$model)) {
      cc <- tryCatch(as.matrix(coef(obj$model)), error = function(e) NULL)
      if (!is.null(cc)) {
        tmp <- data.table(term = rownames(cc), coefficient = as.numeric(cc[, 1]))
        tmp <- tmp[term != "(Intercept)"]
        tmp[, `:=`(selector = selector, fit_stage = fit_stage, source_file = path)]
        coef_dt <- rbind(coef_dt, tmp, fill = TRUE)
      }
    }
    if (!is.null(obj$pca) && !is.null(obj$pca$x) && ncol(obj$pca$x) >= 2) {
      tmp <- data.table(sample_key = rownames(obj$pca$x), PC1 = obj$pca$x[, 1], PC2 = obj$pca$x[, 2])
      tmp <- merge(tmp, unique(pred[, .(sample_key, modality, dataset_folder, age_week)]), by = "sample_key", all.x = TRUE)
      tmp[, `:=`(selector = selector, fit_stage = fit_stage, source_file = path)]
      score_dt <- rbind(score_dt, tmp, fill = TRUE)
    }
  }
  if (nrow(coef_dt) > 0) {
    fwrite(coef_dt, file.path(out_dir, "pca_glmnet_pc_coefficients.tsv"), sep = "\t")
    p <- ggplot(coef_dt, aes(x = term, y = coefficient, fill = selector)) +
      geom_col(color = "grey25") +
      facet_wrap(~fit_stage, scales = "free_x") +
      labs(title = "PCA/glmnet PC coefficients", x = "PC", y = "Coefficient", fill = "Selector") +
      theme_diag()
    save_plot(p, "13_pca_diagnostics_pc_coefficients", 11, 6)
  } else {
    add_warning("No PCA/glmnet PC coefficients could be extracted.")
  }
  if (nrow(score_dt) > 0) {
    fwrite(score_dt, file.path(out_dir, "pca_score_summary_from_rds.tsv"), sep = "\t")
    score_dt[, pca_score_panel := paste(selector_label(selector), fit_stage)]
    score_dt[, pca_score_panel := factor(
      pca_score_panel,
      levels = c(
        "HVG train_fit",
        "HVG train_to_validation",
        "HVG locked_train_validation_to_test",
        "Version B train_fit",
        "Version B train_to_validation",
        "Version B locked_train_validation_to_test"
      )
    )]
    p <- ggplot(score_dt, aes(x = PC1, y = PC2, color = age_week, shape = modality)) +
      geom_point(size = 2.4, alpha = 0.9) +
      facet_wrap(~pca_score_panel, scales = "free", ncol = 3, drop = FALSE) +
      scale_color_viridis_c(option = "C", na.value = "grey70") +
      labs(title = "PCA score diagnostics extracted from model RDS", color = "Age (weeks)", shape = "Modality") +
      theme_diag()
    save_plot(p, "13_pca_diagnostics_pc1_pc2_scores", 12, 8)
  } else {
    add_warning("No PCA score matrix could be extracted from model RDS files.")
  }
}

# N. Projection readiness audit.
all_files <- list.files(results_dir, recursive = TRUE, full.names = TRUE)
all_files <- all_files[!startsWith(normalizePath(all_files, mustWork = FALSE), normalizePath(out_dir, mustWork = FALSE))]
artifact_present <- function(patterns) {
  any(vapply(patterns, function(pat) any(grepl(pat, all_files)), logical(1)))
}
projection_audit <- rbindlist(list(
  data.table(model_family = "elasticnet", requirement = "selected genes", present = artifact_present(c("elasticnet.*selected_genes\\.tsv", "selected_genes/integrated_bulk_pseudobulk_elasticnet"))),
  data.table(model_family = "elasticnet", requirement = "gene centers/scales or preprocessing object", present = artifact_present(c("elasticnet_model\\.rds"))),
  data.table(model_family = "elasticnet", requirement = "glmnet coefficients/model object", present = artifact_present(c("elasticnet_coefficients\\.tsv", "elasticnet_model\\.rds"))),
  data.table(model_family = "elasticnet", requirement = "calibration model", present = artifact_present(c("calibration.*rds", "calibration.*tsv"))),
  data.table(model_family = "pca_glmnet", requirement = "selected genes", present = artifact_present(c("pca_glmnet.*selected_genes\\.tsv", "selected_genes/integrated_bulk_pseudobulk_pca_glmnet"))),
  data.table(model_family = "pca_glmnet", requirement = "gene centers/scales", present = artifact_present(c("pca_glmnet_model\\.rds"))),
  data.table(model_family = "pca_glmnet", requirement = "PCA rotation/loadings", present = artifact_present(c("pca_glmnet_model\\.rds"))),
  data.table(model_family = "pca_glmnet", requirement = "selected PC list", present = artifact_present(c("pca_glmnet_model\\.rds", "pca_glmnet_cv_grid\\.tsv"))),
  data.table(model_family = "pca_glmnet", requirement = "glmnet model/coefficients", present = artifact_present(c("pca_glmnet_model\\.rds", "pca_glmnet.*coefficients"))),
  data.table(model_family = "pca_glmnet", requirement = "calibration model", present = artifact_present(c("calibration.*rds", "calibration.*tsv"))),
  data.table(model_family = "pca_glmnet", requirement = "PC contribution table", present = artifact_present(c("pc_contribution", "contribution_table")))
), fill = TRUE)
projection_audit[, missing_action := fifelse(
  present,
  "available",
  fifelse(
    requirement == "calibration model",
    "Save calibration intercept/slope or calibration lm object from the training script.",
    fifelse(
      requirement == "PC contribution table",
      "Save per-PC and per-gene contribution tables during PCA/glmnet training if interpretability is needed.",
      paste("Save", requirement, "from the training script.")
    )
  )
)]
projection_audit[, projection_ready := all(present), by = model_family]
fwrite(projection_audit, file.path(out_dir, "projection_readiness_audit.tsv"), sep = "\t")
p <- ggplot(projection_audit, aes(x = requirement, y = model_family, fill = present)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = fifelse(present, "present", "missing")), size = 3.6) +
  scale_fill_manual(values = c("TRUE" = "#4C78A8", "FALSE" = "#F58518")) +
  labs(title = "Projection readiness audit", x = "Required artifact", y = "Model family", fill = "Present") +
  theme_diag() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_plot(p, "14_projection_readiness_audit", 12, 5.5)

# Summary markdown.
worst <- locked_eval[order(-abs_error_week)][1:min(.N, 12), .(
  sample_key, model_family, selector, split_role, modality, dataset_folder, age_week, predicted_age_week, residual_week, abs_error_week
)]
fwrite(worst, file.path(out_dir, "worst_residual_samples.tsv"), sep = "\t")

dataset_signal <- residual_bias_by_dataset[order(-mae_week)][1:min(.N, 6)]
modality_signal <- performance_by_modality[order(-mae_week)][1:min(.N, 6)]
bio_signal <- performance_by_bio[order(-mae_week)][1:min(.N, 6)]
naive_overall <- if (file.exists(file.path(out_dir, "naive_baseline_performance.tsv"))) {
  fread(file.path(out_dir, "naive_baseline_performance.tsv"), sep = "\t")
} else data.table()
projection_status <- projection_audit[, .(projection_ready = all(present)), by = model_family]
locked_summary_line <- function(family) {
  row <- locked_wide[as.character(model_family) == family]
  if (nrow(row) == 0) return(paste0("- ", family_label(family), ": locked selector unavailable."))
  paste0(
    "- ", family_label(family), ": locked selector `", row$locked_selector[[1]],
    "`; validation MAE ", sprintf("%.2f", row$validation_mae_week[[1]]),
    " weeks; test MAE ", sprintf("%.2f", row$test_mae_week[[1]]), " weeks."
  )
}
pca_locked_selector <- locked_wide[as.character(model_family) == "pca_glmnet", locked_selector][1]
elastic_locked_selector <- locked_wide[as.character(model_family) == "elasticnet", locked_selector][1]
version_b_statement <- if (identical(pca_locked_selector, "version_b_age_dataset_penalized")) {
  "- Version B is useful for PCA/glmnet in this benchmark because it reduced validation MAE relative to the PCA/glmnet HVG baseline and was the locked PCA/glmnet selector."
} else {
  paste0("- PCA/glmnet locked selector was `", pca_locked_selector, "` in this run; Version B plots remain descriptive feature-selection diagnostics, not the retained PCA/glmnet deployment selector.")
}

md <- c(
  "# Integrated Clock Diagnostic Summary",
  "",
  "## Run Summary",
  paste0("- Training scope: `", expected_scope, "`."),
  paste0("- Prediction rows analyzed: ", nrow(pred), "."),
  paste0("- Unique samples in prediction table: ", all_samples_n, "."),
  "- Both ElasticNet and PCA/glmnet are retained as family-specific integrated clocks.",
  "",
  "## Locked Selectors And MAE",
  locked_summary_line("elasticnet"),
  locked_summary_line("pca_glmnet"),
  "",
  "## Naive Baseline",
  if (nrow(naive_overall) > 0 && "method" %in% names(naive_overall)) {
    paste0("- Naive baseline table written to `naive_baseline_performance.tsv`; locked clocks are compared against train mean and train median predictors overall and by modality.")
  } else {
    "- Naive baseline comparison was skipped because training rows were unavailable."
  },
  "",
  "## Worst Residual Samples",
  paste0("- Worst residual sample table written to `worst_residual_samples.tsv`; maximum locked absolute error was ",
    sprintf("%.2f", max(worst$abs_error_week, na.rm = TRUE)), " weeks."),
  "",
  "## Error Structure",
  paste0("- Highest dataset-level MAE observed in: ", paste(dataset_signal[, paste0(model_family, "/", split_role, "/", dataset_folder, "=", sprintf("%.2f", mae_week))], collapse = "; "), "."),
  paste0("- Highest modality-level MAE observed in: ", paste(modality_signal[, paste0(model_family, "/", split_role, "/", modality, "=", sprintf("%.2f", mae_week))], collapse = "; "), "."),
  paste0("- Highest biological-age-bin MAE observed in: ", paste(bio_signal[, paste0(model_family, "/", split_role, "/", biological_age_bin, "=", sprintf("%.2f", mae_week))], collapse = "; "), "."),
  "",
  "## Feature-Selector Utility",
  version_b_statement,
  paste0("- ElasticNet locked selector was `", elastic_locked_selector, "`."),
  "",
  "## Projection Readiness",
  paste0("- ElasticNet projection ready: ", projection_status[model_family == "elasticnet", projection_ready], "."),
  paste0("- PCA/glmnet projection ready: ", projection_status[model_family == "pca_glmnet", projection_ready], "."),
  "- Missing artifacts and exact save requirements are listed in `projection_readiness_audit.tsv`.",
  "",
  "## Recommended Next Step",
  "- Save calibration intercept/slope objects in the training workflow before external projection, then run excluded-sample projection with both retained locked clocks."
)
writeLines(md, file.path(out_dir, "integrated_clock_diagnostic_summary.md"))

if (length(warnings_log) > 0) {
  fwrite(data.table(warning = warnings_log), file.path(out_dir, "diagnostic_warnings.tsv"), sep = "\t")
} else {
  fwrite(data.table(warning = character()), file.path(out_dir, "diagnostic_warnings.tsv"), sep = "\t")
}
fwrite(plot_manifest, file.path(out_dir, "plot_manifest.tsv"), sep = "\t")

cat("Integrated diagnostics output directory:", out_dir, "\n")
cat("Plots created:", nrow(plot_manifest), "\n")
cat("Warnings:", length(warnings_log), "\n")
