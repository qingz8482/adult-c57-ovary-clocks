#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(scales)
})

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)

options(stringsAsFactors = FALSE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
plot_dir <- file.path(results_dir, "readable_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

locked_path <- file.path(results_dir, "locked_model_selection.tsv")
pred_path <- file.path(results_dir, "all_clock_predictions.tsv")
perf_path <- file.path(results_dir, "all_clock_performance.tsv")
stopifnot(file.exists(locked_path), file.exists(pred_path), file.exists(perf_path))

warnings_log <- character()
plot_manifest <- tibble()

add_warning <- function(msg) {
  warnings_log <<- c(warnings_log, msg)
  message("WARNING: ", msg)
}

detect_col <- function(df, candidates, context, required = TRUE) {
  nms <- names(df)
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[[1]])
  hit_lower <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(hit_lower) > 0) return(hit_lower[[1]])
  msg <- paste0("Could not detect ", context, " column. Tried: ", paste(candidates, collapse = ", "))
  if (required) stop(msg, call. = FALSE)
  add_warning(msg)
  NA_character_
}

optional_col <- function(df, candidates, context) {
  detect_col(df, candidates, context, required = FALSE)
}

numeric_or_na <- function(x) suppressWarnings(as.numeric(x))

selector_label <- function(x) {
  case_when(
    str_detect(tolower(x), "hvg") ~ "HVG",
    str_detect(tolower(x), "version_b|dataset_penalized|age_dataset") ~ "Version B",
    TRUE ~ as.character(x)
  )
}

family_label <- function(x) {
  case_when(
    x == "pca_glmnet" ~ "PCA/glmnet",
    x == "elasticnet" ~ "ElasticNet",
    TRUE ~ as.character(x)
  )
}

modality_label <- function(x) {
  case_when(
    x == "bulk_rnaseq" ~ "Bulk RNA-seq",
    x == "scrna_pseudobulk" ~ "scRNA pseudobulk",
    TRUE ~ as.character(x)
  )
}

bio_bin_from_age <- function(age_week) {
  case_when(
    is.na(age_week) ~ NA_character_,
    age_week >= 7 & age_week < 16 ~ "young_adult",
    age_week >= 16 & age_week < 40 ~ "mature_adult",
    age_week >= 40 & age_week < 65 ~ "advanced_maternal_age",
    age_week >= 65 ~ "post_reproductive_ageing",
    TRUE ~ NA_character_
  )
}

split_stratum_from_age <- function(age_week) {
  case_when(
    is.na(age_week) ~ NA_character_,
    age_week >= 7 & age_week < 13 ~ "early_adult",
    age_week >= 13 & age_week < 39 ~ "reproductive_adult",
    age_week >= 39 ~ "aged_adult",
    TRUE ~ NA_character_
  )
}

corr_or_na <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) return(NA_real_)
  suppressWarnings(stats::cor(x[ok], y[ok], method = method))
}

metric_table <- function(df, by_cols) {
  df %>%
    group_by(across(all_of(by_cols))) %>%
    summarise(
      n = n(),
      mae_week = mean(abs_error_week, na.rm = TRUE),
      rmse_week = sqrt(mean(residual_week^2, na.rm = TRUE)),
      median_abs_error_week = median(abs_error_week, na.rm = TRUE),
      mean_residual_week = mean(residual_week, na.rm = TRUE),
      median_residual_week = median(residual_week, na.rm = TRUE),
      residual_sd_week = sd(residual_week, na.rm = TRUE),
      pearson = corr_or_na(age_week, predicted_age_week, "pearson"),
      spearman = corr_or_na(age_week, predicted_age_week, "spearman"),
      .groups = "drop"
    )
}

save_plot <- function(plot, prefix, title, width = 10, height = 7) {
  png_path <- file.path(plot_dir, paste0(prefix, ".png"))
  pdf_path <- file.path(plot_dir, paste0(prefix, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, bg = "white", limitsize = FALSE)
  plot_manifest <<- bind_rows(
    plot_manifest,
    tibble(
      prefix = prefix,
      title = title,
      png_path = png_path,
      pdf_path = pdf_path,
      width = width,
      height = height
    )
  )
}

theme_clock <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey30"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

read_main_tsv <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

standardize_predictions <- function(df) {
  sample_col <- detect_col(df, c("sample_key", "sample_id", "sample"), "sample key")
  sample_id_col <- optional_col(df, c("sample_id", "sample", "sample_key"), "sample display id")
  modality_col <- detect_col(df, c("modality", "data_type"), "modality")
  family_col <- detect_col(df, c("model_family", "model_type", "clock_model", "method"), "model family")
  selector_col <- detect_col(df, c("selector", "feature_selector", "gene_selector", "feature_selection", "model_version"), "selector")
  split_col <- detect_col(df, c("split", "split_role", "role", "data_split"), "split")
  age_col <- detect_col(df, c("age_week", "observed_age_week", "observed_age"), "age")
  pred_col <- detect_col(df, c("pred_age_week", "predicted_age_week", "predicted_age", "prediction"), "prediction")
  dataset_col <- optional_col(df, c("dataset_folder", "dataset_label", "dataset"), "dataset")
  fit_stage_col <- optional_col(df, c("fit_stage", "stage", "training_stage"), "fit stage")
  raw_pred_col <- optional_col(df, c("raw_predicted_age_week", "raw_prediction", "raw_pred"), "raw prediction")

  out <- tibble(
    sample_key = as.character(df[[sample_col]]),
    sample_id = if (!is.na(sample_id_col)) as.character(df[[sample_id_col]]) else as.character(df[[sample_col]]),
    modality = as.character(df[[modality_col]]),
    model_family = as.character(df[[family_col]]),
    selector = as.character(df[[selector_col]]),
    split = as.character(df[[split_col]]),
    age_week = numeric_or_na(df[[age_col]]),
    predicted_age_week = numeric_or_na(df[[pred_col]]),
    dataset_folder = if (!is.na(dataset_col)) as.character(df[[dataset_col]]) else NA_character_,
    fit_stage = if (!is.na(fit_stage_col)) as.character(df[[fit_stage_col]]) else NA_character_,
    raw_predicted_age_week = if (!is.na(raw_pred_col)) numeric_or_na(df[[raw_pred_col]]) else NA_real_
  )

  essential_missing <- is.na(out$sample_key) | is.na(out$modality) | is.na(out$model_family) |
    is.na(out$selector) | is.na(out$split) | is.na(out$age_week) | is.na(out$predicted_age_week)
  if (any(essential_missing)) {
    add_warning(paste0("Dropped ", sum(essential_missing), " prediction rows because essential columns were missing."))
  }

  out %>%
    filter(!essential_missing) %>%
    mutate(
      residual_week = predicted_age_week - age_week,
      abs_error_week = abs(residual_week),
      bio_age_bin = bio_bin_from_age(age_week),
      split_age_stratum = split_stratum_from_age(age_week),
      selector_label = selector_label(selector),
      model_family_label = family_label(model_family),
      modality_label = modality_label(modality),
      model_id = paste(modality_label, model_family_label, selector_label, sep = " | "),
      split = factor(split, levels = unique(c("train", "validation", "test", sort(unique(split)))))
    )
}

standardize_locked <- function(df) {
  modality_col <- detect_col(df, c("modality", "data_type"), "locked modality")
  family_col <- detect_col(df, c("model_family", "model_type", "clock_model", "method"), "locked model family")
  selector_col <- detect_col(df, c("selector", "feature_selector", "gene_selector", "feature_selection", "model_version"), "locked selector")
  n_col <- optional_col(df, c("n_samples", "n", "sample_count"), "locked n")
  mae_col <- optional_col(df, c("mae_week", "mae", "mean_absolute_error"), "locked MAE")
  rmse_col <- optional_col(df, c("rmse_week", "rmse"), "locked RMSE")
  pearson_col <- optional_col(df, c("pearson", "pearson_r"), "locked Pearson")
  spearman_col <- optional_col(df, c("spearman", "spearman_rho"), "locked Spearman")

  tibble(
    modality = as.character(df[[modality_col]]),
    model_family = as.character(df[[family_col]]),
    selector = as.character(df[[selector_col]]),
    validation_n = if (!is.na(n_col)) numeric_or_na(df[[n_col]]) else NA_real_,
    validation_mae_week = if (!is.na(mae_col)) numeric_or_na(df[[mae_col]]) else NA_real_,
    validation_rmse_week = if (!is.na(rmse_col)) numeric_or_na(df[[rmse_col]]) else NA_real_,
    validation_pearson = if (!is.na(pearson_col)) numeric_or_na(df[[pearson_col]]) else NA_real_,
    validation_spearman = if (!is.na(spearman_col)) numeric_or_na(df[[spearman_col]]) else NA_real_
  ) %>%
    distinct(modality, model_family, selector, .keep_all = TRUE) %>%
    mutate(
      selector_label = selector_label(selector),
      model_family_label = family_label(model_family),
      modality_label = modality_label(modality),
      model_id = paste(modality_label, model_family_label, selector_label, sep = " | ")
    )
}

parse_result_path <- function(paths) {
  rel <- str_remove(paths, paste0("^", stringr::fixed(results_dir), "/"))
  parts <- str_split(rel, "/", simplify = TRUE)
  tibble(
    path = paths,
    result_family_dir = parts[, 1],
    modality = parts[, 2],
    selector = parts[, 3],
    model_family = case_when(
      parts[, 1] == "pca_glmnet_clock" ~ "pca_glmnet",
      parts[, 1] == "elasticnet_clock" ~ "elasticnet",
      TRUE ~ parts[, 1]
    ),
    modality_label = modality_label(parts[, 2]),
    model_family_label = family_label(case_when(
      parts[, 1] == "pca_glmnet_clock" ~ "pca_glmnet",
      parts[, 1] == "elasticnet_clock" ~ "elasticnet",
      TRUE ~ parts[, 1]
    )),
    selector_label = selector_label(parts[, 3])
  )
}

locked_raw <- read_main_tsv(locked_path)
pred_raw <- read_main_tsv(pred_path)
perf_raw <- read_main_tsv(perf_path)

pred <- standardize_predictions(pred_raw)
locked <- standardize_locked(locked_raw)
locked_keys <- locked %>% select(modality, model_family, selector)
locked_pred <- pred %>% semi_join(locked_keys, by = c("modality", "model_family", "selector"))
locked_test <- locked_pred %>% filter(split == "test")

if (nrow(locked_test) == 0) {
  add_warning("No locked test predictions were found. Main test plots will be skipped.")
}

cohort <- pred %>%
  distinct(sample_key, modality, modality_label, split, dataset_folder, age_week, bio_age_bin, split_age_stratum)

write_tsv(
  locked %>%
    arrange(modality, model_family) %>%
    select(modality, modality_label, model_family, model_family_label, selector, selector_label, starts_with("validation")),
  file.path(plot_dir, "locked_model_selection_clean.tsv")
)

# A. Cohort and split overview.
cohort_dataset <- cohort %>%
  count(modality_label, split, dataset_folder, name = "n") %>%
  mutate(dataset_folder_wrapped = str_wrap(dataset_folder, width = 28))

p <- ggplot(cohort_dataset, aes(x = split, y = n, fill = dataset_folder_wrapped)) +
  geom_col(color = "grey20", linewidth = 0.15) +
  facet_wrap(~modality_label, scales = "free_y") +
  labs(
    title = "Cohort Overview by Dataset",
    subtitle = "Unique samples counted once per split and modality.",
    x = "Split",
    y = "Number of samples",
    fill = "Dataset"
  ) +
  theme_clock(14) +
  guides(fill = guide_legend(ncol = 1))
save_plot(p, "00_cohort_overview_by_dataset", "Cohort overview by dataset", width = 12, height = 8)

bio_levels <- c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
cohort_bio <- cohort %>%
  mutate(bio_age_bin = factor(bio_age_bin, levels = bio_levels)) %>%
  count(modality_label, split, bio_age_bin, name = "n")

p <- ggplot(cohort_bio, aes(x = split, y = n, fill = bio_age_bin)) +
  geom_col(color = "grey20", linewidth = 0.15) +
  facet_wrap(~modality_label, scales = "free_y") +
  labs(
    title = "Cohort Overview by Biological Age Bin",
    subtitle = "Biological bins are derived from continuous age in weeks.",
    x = "Split",
    y = "Number of samples",
    fill = "Biological age bin"
  ) +
  theme_clock(14)
save_plot(p, "00_cohort_overview_by_biological_age_bin", "Cohort overview by biological age bin", width = 11, height = 7)

stratum_levels <- c("early_adult", "reproductive_adult", "aged_adult")
cohort_stratum <- cohort %>%
  mutate(split_age_stratum = factor(split_age_stratum, levels = stratum_levels)) %>%
  count(modality_label, split, split_age_stratum, name = "n")

p <- ggplot(cohort_stratum, aes(x = split, y = n, fill = split_age_stratum)) +
  geom_col(color = "grey20", linewidth = 0.15) +
  facet_wrap(~modality_label, scales = "free_y") +
  labs(
    title = "Cohort Overview by Balanced Split Age Stratum",
    subtitle = "Split strata are derived from continuous age in weeks.",
    x = "Split",
    y = "Number of samples",
    fill = "Split age stratum"
  ) +
  theme_clock(14)
save_plot(p, "00_cohort_overview_by_split_age_stratum", "Cohort overview by split age stratum", width = 11, height = 7)

# B. Locked model selection summary.
locked_table <- locked %>%
  mutate(
    tile_label = sprintf(
      "%s\nvalidation MAE: %s wk\nn = %s",
      selector_label,
      if_else(is.na(validation_mae_week), "NA", sprintf("%.2f", validation_mae_week)),
      if_else(is.na(validation_n), "NA", as.character(as.integer(validation_n)))
    )
  )

p <- ggplot(locked_table, aes(x = model_family_label, y = modality_label, fill = selector_label)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = tile_label), size = 4.2, lineheight = 0.95) +
  scale_fill_manual(values = c("HVG" = "#4C78A8", "Version B" = "#F58518"), na.value = "grey75") +
  labs(
    title = "Locked Model Selection",
    subtitle = "Locked choices are read directly from locked_model_selection.tsv; test results are not used for selection.",
    x = "Model family",
    y = "Modality",
    fill = "Selector"
  ) +
  theme_clock(15) +
  theme(axis.text.x = element_text(angle = 0))
save_plot(p, "01_locked_model_selection", "Locked model selection", width = 9, height = 5.5)

# C. Locked test MAE.
locked_test_metrics <- metric_table(
  locked_test,
  c("modality", "modality_label", "model_family", "model_family_label", "selector", "selector_label")
) %>%
  mutate(model_plot_label = paste(modality_label, model_family_label, selector_label, sep = "\n")) %>%
  arrange(modality, model_family)

write_tsv(locked_test_metrics, file.path(plot_dir, "02_locked_test_mae_metrics.tsv"))

p <- ggplot(locked_test_metrics, aes(x = fct_reorder(model_plot_label, mae_week), y = mae_week, fill = modality_label)) +
  geom_col(width = 0.72, color = "grey20") +
  geom_text(aes(label = sprintf("%.2f wk", mae_week)), vjust = -0.35, size = 4.2) +
  labs(
    title = "Locked Test MAE",
    subtitle = "Frozen test split only; lower is better.",
    x = "Locked model",
    y = "MAE (weeks)",
    fill = "Modality"
  ) +
  expand_limits(y = max(locked_test_metrics$mae_week, na.rm = TRUE) * 1.18) +
  theme_clock(14)
save_plot(p, "02_locked_test_mae", "Locked test MAE", width = 10.5, height = 6.5)

# D. Predicted vs observed age.
if (nrow(locked_test) > 0) {
  lims <- range(c(locked_test$age_week, locked_test$predicted_age_week), na.rm = TRUE)
  pad <- diff(lims) * 0.06
  if (!is.finite(pad) || pad == 0) pad <- 1
  lims <- c(lims[1] - pad, lims[2] + pad)
  facet_metrics <- locked_test_metrics %>%
    mutate(
      label = sprintf("n=%d\nMAE=%.2f\nRMSE=%.2f\nPearson=%s\nSpearman=%s",
        n, mae_week, rmse_week,
        if_else(is.na(pearson), "NA", sprintf("%.2f", pearson)),
        if_else(is.na(spearman), "NA", sprintf("%.2f", spearman))
      ),
      x_pos = lims[1] + diff(lims) * 0.03,
      y_pos = lims[2] - diff(lims) * 0.04
    )
  worst <- locked_test %>%
    group_by(modality_label, model_family_label) %>%
    slice_max(abs_error_week, n = 3, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(sample_label = str_trunc(sample_id, width = 24))

  p <- ggplot(locked_test, aes(x = age_week, y = predicted_age_week)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey35") +
    geom_point(aes(color = dataset_folder), size = 2.8, alpha = 0.9) +
    geom_label(
      data = facet_metrics,
      aes(x = x_pos, y = y_pos, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.4,
      linewidth = 0.2,
      fill = "white"
    ) +
    facet_grid(modality_label ~ model_family_label) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = "Locked Test Predictions vs Observed Age",
      subtitle = "Labels mark the three largest absolute residuals per facet.",
      x = "Observed age (weeks)",
      y = "Predicted age (weeks)",
      color = "Dataset"
    ) +
    theme_clock(14) +
    guides(color = guide_legend(ncol = 1))
  if (has_ggrepel && nrow(worst) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = worst,
      aes(label = sample_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  } else {
    add_warning("ggrepel is not installed; skipped worst-residual labels for predicted-vs-observed plot.")
  }
  save_plot(p, "03_locked_predicted_vs_observed_test", "Locked predicted vs observed test age", width = 12, height = 8)
}

# E. Residual vs observed age.
if (nrow(locked_test) > 0) {
  smooth_data <- locked_test %>%
    group_by(modality_label, model_family_label) %>%
    filter(n() >= 6) %>%
    ungroup()
  p <- ggplot(locked_test, aes(x = age_week, y = residual_week)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
    geom_point(aes(color = dataset_folder), size = 2.8, alpha = 0.9) +
    facet_grid(modality_label ~ model_family_label) +
    labs(
      title = "Locked Test Residuals vs Observed Age",
      subtitle = "Residual = predicted age - observed age; negative values indicate underprediction.",
      x = "Observed age (weeks)",
      y = "Residual (weeks)",
      color = "Dataset"
    ) +
    theme_clock(14) +
    guides(color = guide_legend(ncol = 1))
  if (nrow(smooth_data) > 0) {
    p <- p + geom_smooth(data = smooth_data, method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.9, color = "black")
  } else {
    add_warning("Skipped residual smooth: no locked test facet had at least six samples.")
  }
  save_plot(p, "04_locked_residual_vs_age_test", "Locked residual vs observed age", width = 12, height = 8)
}

# F. Residual by dataset.
eval_locked <- locked_pred %>% filter(split %in% c("validation", "test"))
if (nrow(eval_locked) == 0) {
  eval_locked <- locked_pred %>% filter(split != "train")
  add_warning("No validation/test rows found for residual-by-dataset plot; used all non-training predictions instead.")
}

dataset_bias <- metric_table(
  eval_locked,
  c("modality", "modality_label", "model_family", "model_family_label", "selector", "selector_label", "split", "dataset_folder")
) %>%
  arrange(modality, model_family, split, dataset_folder)
write_tsv(dataset_bias, file.path(plot_dir, "residual_bias_by_dataset.tsv"))

p <- eval_locked %>%
  mutate(dataset_plot = fct_reorder(str_wrap(dataset_folder, 28), residual_week, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = dataset_plot, y = residual_week, color = split)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey35") +
  geom_boxplot(outlier.shape = NA, color = "grey25") +
  geom_jitter(width = 0.18, height = 0, size = 2, alpha = 0.75) +
  coord_flip() +
  facet_grid(modality_label ~ model_family_label, scales = "free_y", space = "free_y") +
  labs(
    title = "Locked Residual Bias by Dataset",
    subtitle = "Validation and test predictions for locked model configurations.",
    x = "Dataset",
    y = "Residual (weeks)",
    color = "Split"
  ) +
  theme_clock(13)
save_plot(p, "05_locked_residual_by_dataset", "Locked residual by dataset", width = 13, height = 9)

p <- eval_locked %>%
  mutate(dataset_plot = fct_reorder(str_wrap(dataset_folder, 28), abs_error_week, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = dataset_plot, y = abs_error_week, color = split)) +
  geom_boxplot(outlier.shape = NA, color = "grey25") +
  geom_jitter(width = 0.18, height = 0, size = 2, alpha = 0.75) +
  coord_flip() +
  facet_grid(modality_label ~ model_family_label, scales = "free_y", space = "free_y") +
  labs(
    title = "Locked Absolute Error by Dataset",
    subtitle = "Validation and test predictions for locked model configurations.",
    x = "Dataset",
    y = "Absolute error (weeks)",
    color = "Split"
  ) +
  theme_clock(13)
save_plot(p, "05_locked_abs_error_by_dataset", "Locked absolute error by dataset", width = 13, height = 9)

# G. Error by biological age bin.
bio_summary <- metric_table(
  locked_test %>% mutate(bio_age_bin = factor(bio_age_bin, levels = bio_levels)),
  c("modality", "modality_label", "model_family", "model_family_label", "selector", "selector_label", "split", "bio_age_bin")
) %>%
  arrange(modality, model_family, bio_age_bin)
write_tsv(bio_summary, file.path(plot_dir, "error_by_biological_age_bin.tsv"))

p <- locked_test %>%
  mutate(bio_age_bin = factor(bio_age_bin, levels = bio_levels)) %>%
  ggplot(aes(x = bio_age_bin, y = abs_error_week, color = bio_age_bin)) +
  geom_boxplot(outlier.shape = NA, color = "grey25") +
  geom_jitter(width = 0.18, height = 0, size = 2.4, alpha = 0.8) +
  facet_grid(modality_label ~ model_family_label) +
  labs(
    title = "Locked Test Error by Biological Age Bin",
    subtitle = "Absolute error in weeks; bins are derived from continuous age.",
    x = "Biological age bin",
    y = "Absolute error (weeks)",
    color = "Age bin"
  ) +
  theme_clock(13)
save_plot(p, "06_locked_error_by_biological_age_bin", "Locked error by biological age bin", width = 12, height = 8)

# H. Error by balanced split stratum.
stratum_summary <- metric_table(
  locked_test %>% mutate(split_age_stratum = factor(split_age_stratum, levels = stratum_levels)),
  c("modality", "modality_label", "model_family", "model_family_label", "selector", "selector_label", "split", "split_age_stratum")
) %>%
  arrange(modality, model_family, split_age_stratum)
write_tsv(stratum_summary, file.path(plot_dir, "error_by_split_age_stratum.tsv"))

p <- locked_test %>%
  mutate(split_age_stratum = factor(split_age_stratum, levels = stratum_levels)) %>%
  ggplot(aes(x = split_age_stratum, y = abs_error_week, color = split_age_stratum)) +
  geom_boxplot(outlier.shape = NA, color = "grey25") +
  geom_jitter(width = 0.18, height = 0, size = 2.4, alpha = 0.8) +
  facet_grid(modality_label ~ model_family_label) +
  labs(
    title = "Locked Test Error by Balanced Split Age Stratum",
    subtitle = "Absolute error in weeks; split strata are derived from continuous age.",
    x = "Balanced split age stratum",
    y = "Absolute error (weeks)",
    color = "Split age stratum"
  ) +
  theme_clock(13)
save_plot(p, "07_locked_error_by_split_age_stratum", "Locked error by split age stratum", width = 12, height = 8)

# I. Validation vs test MAE for candidate models.
candidate_metrics <- metric_table(
  pred %>% filter(split %in% c("validation", "test")),
  c("modality", "modality_label", "model_family", "model_family_label", "selector", "selector_label", "split")
)

candidate_vt <- candidate_metrics %>%
  select(modality, modality_label, model_family, model_family_label, selector, selector_label, split, mae_week) %>%
  mutate(split = as.character(split)) %>%
  pivot_wider(names_from = split, values_from = mae_week, names_prefix = "mae_") %>%
  left_join(locked_keys %>% mutate(is_locked = TRUE), by = c("modality", "model_family", "selector")) %>%
  mutate(is_locked = if_else(is.na(is_locked), FALSE, is_locked))
write_tsv(candidate_vt, file.path(plot_dir, "08_candidate_validation_vs_test_mae.tsv"))

missing_test_n <- sum(!is.na(candidate_vt$mae_validation) & is.na(candidate_vt$mae_test))
if (missing_test_n > 0) {
  add_warning(paste0("Candidate validation-vs-test plot has ", missing_test_n, " candidate configurations without test MAE because only locked models were evaluated on test."))
}

candidate_vt_plot <- candidate_vt %>% filter(!is.na(mae_validation), !is.na(mae_test))
if (nrow(candidate_vt_plot) > 0) {
  lims <- range(c(candidate_vt_plot$mae_validation, candidate_vt_plot$mae_test), na.rm = TRUE)
  lims <- c(0, max(lims, na.rm = TRUE) * 1.15)
  p <- ggplot(candidate_vt_plot, aes(x = mae_validation, y = mae_test, color = modality_label, shape = model_family_label)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(size = 3.4) +
    facet_wrap(~selector_label) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(
      title = "Candidate Validation vs Test MAE",
      subtitle = "Only locked configurations have test MAE in the first-pass output.",
      x = "Validation MAE (weeks)",
      y = "Test MAE (weeks)",
      color = "Modality",
      shape = "Model family"
    ) +
    theme_clock(14)
  locked_labels <- candidate_vt_plot %>%
    filter(is_locked) %>%
    mutate(label = paste(modality_label, model_family_label, sep = "\n"))
  if (has_ggrepel && nrow(locked_labels) > 0) {
    p <- p + ggrepel::geom_text_repel(data = locked_labels, aes(label = label), size = 3.4, show.legend = FALSE)
  } else if (nrow(locked_labels) > 0) {
    p <- p + geom_text(data = locked_labels, aes(label = label), size = 3.2, vjust = -0.8, show.legend = FALSE)
  }
  save_plot(p, "08_candidate_validation_vs_test_mae", "Candidate validation vs test MAE", width = 10.5, height = 6.5)
} else {
  add_warning("Skipped candidate validation-vs-test plot because no candidate had both validation and test MAE.")
}

# J. HVG vs Version B comparison.
delta_source <- candidate_metrics %>%
  mutate(selector_class = selector_label) %>%
  filter(selector_class %in% c("HVG", "Version B")) %>%
  select(modality, modality_label, model_family, model_family_label, split, selector_class, mae_week)

delta_mae <- delta_source %>%
  mutate(split = as.character(split)) %>%
  pivot_wider(names_from = selector_class, values_from = mae_week) %>%
  mutate(delta_MAE = `Version B` - HVG)

missing_delta <- delta_mae %>% filter(is.na(delta_MAE))
if (nrow(missing_delta) > 0) {
  add_warning(paste0("HVG-vs-Version-B delta MAE has ", nrow(missing_delta), " incomplete modality/model/split combinations. This is expected for test if non-locked models were not test-evaluated."))
}
write_tsv(delta_mae, file.path(plot_dir, "hvg_vs_versionB_delta_mae.tsv"))

delta_plot <- delta_mae %>% filter(!is.na(delta_MAE))
if (nrow(delta_plot) > 0) {
  p <- delta_plot %>%
    mutate(
      split = factor(split, levels = c("validation", "test")),
      model_plot_label = paste(modality_label, model_family_label, sep = "\n")
    ) %>%
    ggplot(aes(x = model_plot_label, y = delta_MAE, fill = modality_label)) +
    geom_hline(yintercept = 0, color = "grey30") +
    geom_col(color = "grey20", width = 0.72) +
    facet_wrap(~split, drop = FALSE) +
    labs(
      title = "Version B minus HVG Delta MAE",
      subtitle = "Negative values mean Version B has lower MAE than HVG.",
      x = "Modality and model family",
      y = "Delta MAE (weeks)",
      fill = "Modality"
    ) +
    theme_clock(14)
  save_plot(p, "09_hvg_vs_versionB_delta_mae", "HVG vs Version B delta MAE", width = 11, height = 6.5)
} else {
  add_warning("Skipped HVG-vs-Version-B delta plot because no complete HVG/Version-B pairs were available.")
}

# K. Prediction dynamic range.
dynamic_summary <- locked_test %>%
  group_by(modality, modality_label, model_family, model_family_label, selector, selector_label) %>%
  summarise(
    n = n(),
    observed_min = min(age_week, na.rm = TRUE),
    observed_max = max(age_week, na.rm = TRUE),
    observed_range = observed_max - observed_min,
    predicted_min = min(predicted_age_week, na.rm = TRUE),
    predicted_max = max(predicted_age_week, na.rm = TRUE),
    predicted_range = predicted_max - predicted_min,
    .groups = "drop"
  )
write_tsv(dynamic_summary, file.path(plot_dir, "10_prediction_dynamic_range_summary.tsv"))

range_long <- dynamic_summary %>%
  mutate(model_plot_label = paste(modality_label, model_family_label, selector_label, sep = "\n")) %>%
  select(model_plot_label, modality_label, observed_range, predicted_range) %>%
  pivot_longer(c(observed_range, predicted_range), names_to = "range_type", values_to = "age_range_week") %>%
  mutate(range_type = recode(range_type, observed_range = "Observed age range", predicted_range = "Predicted age range"))

p_range <- ggplot(range_long, aes(x = model_plot_label, y = age_range_week, fill = range_type)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68, color = "grey20") +
  labs(
    title = "Observed vs Predicted Age Dynamic Range",
    subtitle = "Range compression suggests prediction collapse toward a mean age.",
    x = "Locked model",
    y = "Age range (weeks)",
    fill = "Range"
  ) +
  theme_clock(13)

dist_long <- locked_test %>%
  select(modality_label, model_family_label, selector_label, age_week, predicted_age_week) %>%
  mutate(model_plot_label = paste(modality_label, model_family_label, selector_label, sep = "\n")) %>%
  pivot_longer(c(age_week, predicted_age_week), names_to = "age_type", values_to = "age_week") %>%
  mutate(age_type = recode(age_type, age_week = "Observed", predicted_age_week = "Predicted"))

p_dist <- ggplot(dist_long, aes(x = age_type, y = age_week, fill = age_type)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "grey20") +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.65) +
  facet_wrap(~model_plot_label) +
  labs(
    title = "Observed and Predicted Test Age Distributions",
    x = NULL,
    y = "Age (weeks)",
    fill = "Age type"
  ) +
  theme_clock(13) +
  theme(axis.text.x = element_text(angle = 0))

if (has_patchwork) {
  p <- p_range / p_dist + patchwork::plot_layout(heights = c(0.95, 1.15))
  save_plot(p, "10_prediction_dynamic_range", "Prediction dynamic range", width = 12, height = 11)
} else {
  add_warning("patchwork is not installed; saved prediction dynamic range as separate range and distribution plots.")
  save_plot(p_range, "10_prediction_dynamic_range_ranges", "Prediction dynamic range ranges", width = 12, height = 6.5)
  save_plot(p_dist, "10_prediction_dynamic_range_distribution", "Prediction dynamic range distribution", width = 12, height = 7.5)
}

# L. Naive baseline comparison.
train_unique <- pred %>% filter(split == "train") %>% distinct(modality, modality_label, sample_key, age_week)
if (nrow(train_unique) == 0) {
  add_warning("Training rows are not available in all_clock_predictions.tsv; skipped naive baseline comparison.")
} else {
  baseline_values <- train_unique %>%
    group_by(modality, modality_label) %>%
    summarise(train_mean_age = mean(age_week, na.rm = TRUE), train_median_age = median(age_week, na.rm = TRUE), .groups = "drop")
  test_samples <- locked_test %>%
    distinct(modality, modality_label, sample_key, age_week)
  baseline_perf <- test_samples %>%
    left_join(baseline_values, by = c("modality", "modality_label")) %>%
    group_by(modality, modality_label) %>%
    summarise(
      split = "test",
      train_mean_baseline_mae = mean(abs(train_mean_age - age_week), na.rm = TRUE),
      train_median_baseline_mae = mean(abs(train_median_age - age_week), na.rm = TRUE),
      .groups = "drop"
    )
  locked_model_perf <- locked_test_metrics %>%
    transmute(
      modality,
      modality_label,
      model_family,
      model_family_label,
      selector,
      selector_label,
      split = "test",
      method = "Locked model",
      mae_week
    )
  baseline_rows <- locked_test_metrics %>%
    select(modality, modality_label, model_family, model_family_label, selector, selector_label) %>%
    left_join(baseline_perf, by = c("modality", "modality_label")) %>%
    pivot_longer(c(train_mean_baseline_mae, train_median_baseline_mae), names_to = "method", values_to = "mae_week") %>%
    mutate(method = recode(method, train_mean_baseline_mae = "Train mean baseline", train_median_baseline_mae = "Train median baseline"))
  baseline_all <- bind_rows(locked_model_perf, baseline_rows) %>%
    mutate(model_plot_label = paste(modality_label, model_family_label, selector_label, sep = "\n"))
  write_tsv(baseline_all, file.path(plot_dir, "naive_baseline_performance.tsv"))

  p <- ggplot(baseline_all, aes(x = model_plot_label, y = mae_week, fill = method)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.7, color = "grey20") +
    geom_text(aes(label = sprintf("%.1f", mae_week)), position = position_dodge(width = 0.78), vjust = -0.3, size = 3.6) +
    labs(
      title = "Locked Models vs Naive Age Baselines",
      subtitle = "Baselines predict the training-set mean or median age within each modality.",
      x = "Locked model",
      y = "Test MAE (weeks)",
      fill = "Prediction method"
    ) +
    expand_limits(y = max(baseline_all$mae_week, na.rm = TRUE) * 1.18) +
    theme_clock(13)
  save_plot(p, "11_naive_baseline_comparison", "Naive baseline comparison", width = 12, height = 7)
}

# M. Optional PCA diagnostics.
pca_root <- file.path(results_dir, "pca_glmnet_clock")
pca_var_files <- list.files(pca_root, pattern = "pca_variance.tsv$", recursive = TRUE, full.names = TRUE)
if (length(pca_var_files) > 0) {
  pca_var <- bind_rows(lapply(pca_var_files, function(path) {
    meta <- parse_result_path(path)
    read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
      mutate(path = path, fit_stage = str_remove(basename(path), "_pca_variance.tsv$")) %>%
      bind_cols(meta %>% select(modality, modality_label, model_family, model_family_label, selector, selector_label))
  }))
  p <- pca_var %>%
    mutate(facet_label = paste(modality_label, selector_label, fit_stage, sep = "\n")) %>%
    ggplot(aes(x = pc, y = variance_explained)) +
    geom_line(color = "#4C78A8", linewidth = 0.8) +
    geom_point(color = "#4C78A8", size = 1.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    facet_wrap(~facet_label, scales = "free_y") +
    labs(
      title = "PCA/glmnet PC Variance Explained",
      subtitle = "Found pca_variance.tsv files under results/pca_glmnet_clock.",
      x = "Principal component",
      y = "Variance explained"
    ) +
    theme_clock(12)
  save_plot(p, "12_pca_variance_explained", "PCA variance explained", width = 13, height = 9)
} else {
  add_warning("No PCA variance TSV files found under results/pca_glmnet_clock.")
}

pca_cv_files <- list.files(pca_root, pattern = "pca_glmnet_cv_grid.tsv$", recursive = TRUE, full.names = TRUE)
if (length(pca_cv_files) > 0) {
  pca_cv <- bind_rows(lapply(pca_cv_files, function(path) {
    meta <- parse_result_path(path)
    read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
      mutate(path = path, fit_stage = str_remove(basename(path), "_pca_glmnet_cv_grid.tsv$")) %>%
      bind_cols(meta %>% select(modality, modality_label, model_family, model_family_label, selector, selector_label))
  }))
  if ("pc_count" %in% names(pca_cv)) {
    p <- pca_cv %>%
      mutate(facet_label = paste(modality_label, selector_label, fit_stage, sep = "\n")) %>%
      ggplot(aes(x = pc_count, y = cv_mae_week, color = factor(alpha))) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.9) +
      facet_wrap(~facet_label, scales = "free_y") +
      labs(
        title = "PCA/glmnet CV MAE by PC Count",
        subtitle = "Lower CV MAE is better; alpha is glmnet mixing parameter.",
        x = "PC count",
        y = "CV MAE (weeks)",
        color = "alpha"
      ) +
      theme_clock(12)
    save_plot(p, "12_pca_cv_mae_by_pc_count", "PCA CV MAE by PC count", width = 13, height = 9)
  }
} else {
  add_warning("No PCA/glmnet CV grid TSV files found for optional PC-count diagnostics.")
}
add_warning("No PC score, PC age-association, or PC dataset-R2 diagnostic files were found; skipped those optional PCA plots.")

# N. Optional selected-gene diagnostics.
selected_files <- list.files(results_dir, pattern = "selected_genes.tsv$", recursive = TRUE, full.names = TRUE)
feature_score_files <- list.files(results_dir, pattern = "feature_score_table.tsv$", recursive = TRUE, full.names = TRUE)
selected_table_files <- list.files(results_dir, pattern = "selected_gene_table.tsv$", recursive = TRUE, full.names = TRUE)

if (length(selected_files) == 0) {
  add_warning("No selected_genes.tsv files found; skipped selected-gene diagnostics.")
} else {
  selected_counts <- bind_rows(lapply(selected_files, function(path) {
    meta <- parse_result_path(path)
    dt <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
    tibble(path = path, n_selected_genes = nrow(dt)) %>%
      bind_cols(meta %>% select(modality, modality_label, model_family, model_family_label, selector, selector_label))
  }))
  write_tsv(selected_counts, file.path(plot_dir, "selected_gene_counts.tsv"))

  p <- selected_counts %>%
    mutate(model_plot_label = paste(modality_label, model_family_label, selector_label, sep = "\n")) %>%
    ggplot(aes(x = model_plot_label, y = n_selected_genes, fill = selector_label)) +
    geom_col(color = "grey20", width = 0.72) +
    geom_text(aes(label = n_selected_genes), vjust = -0.3, size = 3.8) +
    labs(
      title = "Number of Selected Genes",
      subtitle = "Counts are read from selected_genes.tsv files.",
      x = "Model",
      y = "Selected genes",
      fill = "Selector"
    ) +
    expand_limits(y = max(selected_counts$n_selected_genes, na.rm = TRUE) * 1.12) +
    theme_clock(13)
  save_plot(p, "13_selected_gene_counts", "Selected gene counts", width = 12, height = 6.5)

  selected_genes <- bind_rows(lapply(selected_files, function(path) {
    meta <- parse_result_path(path)
    dt <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
    gene_col <- detect_col(dt, c("gene_id", "gene", "symbol"), paste0("gene column in ", path))
    tibble(gene_id = as.character(dt[[gene_col]]), path = path) %>%
      bind_cols(meta %>% select(modality, modality_label, model_family, model_family_label, selector, selector_label))
  }))

  overlap_rows <- selected_genes %>%
    filter(selector_label %in% c("HVG", "Version B")) %>%
    group_by(modality, modality_label, model_family, model_family_label) %>%
    summarise(
      hvg_genes = list(unique(gene_id[selector_label == "HVG"])),
      version_b_genes = list(unique(gene_id[selector_label == "Version B"])),
      .groups = "drop"
    ) %>%
    mutate(
      hvg_n = lengths(hvg_genes),
      version_b_n = lengths(version_b_genes),
      overlap_n = mapply(function(a, b) length(intersect(a, b)), hvg_genes, version_b_genes),
      overlap_fraction_of_hvg = if_else(hvg_n > 0, overlap_n / hvg_n, NA_real_),
      overlap_fraction_of_version_b = if_else(version_b_n > 0, overlap_n / version_b_n, NA_real_)
    ) %>%
    select(-hvg_genes, -version_b_genes)
  write_tsv(overlap_rows, file.path(plot_dir, "selected_gene_overlap_hvg_vs_versionB.tsv"))

  if (nrow(overlap_rows) > 0) {
    p <- overlap_rows %>%
      mutate(model_plot_label = paste(modality_label, model_family_label, sep = "\n")) %>%
      ggplot(aes(x = model_plot_label, y = overlap_n, fill = modality_label)) +
      geom_col(color = "grey20", width = 0.72) +
      geom_text(aes(label = overlap_n), vjust = -0.3, size = 4) +
      labs(
        title = "Selected-Gene Overlap: HVG vs Version B",
        subtitle = "Overlap is computed within each modality and model family.",
        x = "Model",
        y = "Overlap size",
        fill = "Modality"
      ) +
      expand_limits(y = max(overlap_rows$overlap_n, na.rm = TRUE) * 1.15) +
      theme_clock(13)
    save_plot(p, "13_selected_gene_overlap_hvg_vs_versionB", "Selected gene overlap HVG vs Version B", width = 10.5, height = 6.5)
  }
}

if (length(feature_score_files) > 0) {
  feature_scores <- bind_rows(lapply(feature_score_files, function(path) {
    meta <- parse_result_path(path)
    read_tsv(path, show_col_types = FALSE, progress = FALSE) %>%
      mutate(path = path) %>%
      bind_cols(meta %>% select(modality, modality_label, model_family, model_family_label, selector, selector_label))
  }))
  if (all(c("partial_r2_age", "partial_r2_dataset") %in% names(feature_scores))) {
    p <- feature_scores %>%
      ggplot(aes(x = partial_r2_dataset, y = partial_r2_age, color = modality_label)) +
      geom_point(alpha = 0.22, size = 0.9) +
      facet_grid(model_family_label ~ modality_label) +
      labs(
        title = "Version B Feature Scores: Age vs Dataset Partial R2",
        subtitle = "Genes toward high age R2 and low dataset R2 are more age-informative and less dataset-dominated.",
        x = "Dataset partial R2",
        y = "Age partial R2",
        color = "Modality"
      ) +
      theme_clock(13)
    save_plot(p, "13_feature_score_age_vs_dataset_partial_r2", "Feature score age vs dataset partial R2", width = 12, height = 8)
  } else {
    add_warning("Feature score files do not contain both partial_r2_age and partial_r2_dataset; skipped feature-score scatter.")
  }
} else {
  add_warning("No feature_score_table.tsv files found; skipped feature-score diagnostics.")
}

if (length(selected_table_files) > 0) {
  top_tables <- bind_rows(lapply(selected_table_files, function(path) {
    meta <- parse_result_path(path)
    dt <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
    gene_col <- detect_col(dt, c("gene_id", "gene", "symbol"), paste0("gene column in ", path))
    score_col <- optional_col(dt, c("feature_score", "score", "age_score"), paste0("score column in ", path))
    if (is.na(score_col)) return(tibble())
    tibble(gene_id = as.character(dt[[gene_col]]), feature_score = numeric_or_na(dt[[score_col]]), path = path) %>%
      bind_cols(meta %>% select(modality, modality_label, model_family, model_family_label, selector, selector_label))
  }))
  if (nrow(top_tables) > 0) {
    write_tsv(top_tables, file.path(plot_dir, "versionB_selected_gene_tables_combined.tsv"))
    for (key in unique(paste(top_tables$modality, top_tables$model_family, sep = "||"))) {
      sub <- top_tables %>%
        filter(paste(modality, model_family, sep = "||") == key) %>%
        arrange(desc(feature_score)) %>%
        slice_head(n = 30)
      if (nrow(sub) == 0) next
      prefix <- paste0(
        "13_top30_versionB_genes_",
        str_replace_all(unique(sub$modality), "[^A-Za-z0-9]+", "_"),
        "_",
        str_replace_all(unique(sub$model_family), "[^A-Za-z0-9]+", "_")
      )
      p <- sub %>%
        mutate(gene_id = fct_reorder(gene_id, feature_score)) %>%
        ggplot(aes(x = gene_id, y = feature_score, fill = modality_label)) +
        geom_col(color = "grey20", width = 0.72) +
        coord_flip() +
        labs(
          title = paste("Top 30 Version B Genes:", unique(sub$modality_label), unique(sub$model_family_label)),
          subtitle = "Ranked by feature_score from selected_gene_table.tsv.",
          x = "Gene",
          y = "Feature score",
          fill = "Modality"
        ) +
        theme_clock(12)
      save_plot(p, prefix, paste("Top 30 Version B genes", unique(sub$modality), unique(sub$model_family)), width = 9, height = 10)
    }
  }
} else {
  add_warning("No selected_gene_table.tsv files found; skipped top selected-gene score plots.")
}

# Quality checks against expected first-pass values.
expected_mae <- tibble(
  modality = c("bulk_rnaseq", "bulk_rnaseq", "scrna_pseudobulk", "scrna_pseudobulk"),
  model_family = c("elasticnet", "pca_glmnet", "pca_glmnet", "elasticnet"),
  selector_label = c("HVG", "HVG", "HVG", "Version B"),
  expected_mae_week = c(18.53, 17.83, 22.57, 19.22)
)

mae_check <- locked_test_metrics %>%
  left_join(expected_mae, by = c("modality", "model_family", "selector_label")) %>%
  mutate(abs_diff = abs(mae_week - expected_mae_week))

cat("\nLocked test MAEs computed by plotting script:\n")
print(mae_check %>% select(modality, model_family, selector_label, mae_week, expected_mae_week, abs_diff))

bad_mae <- mae_check %>% filter(!is.na(expected_mae_week), abs_diff > 0.1)
if (nrow(bad_mae) > 0) {
  add_warning(paste0("Locked test MAE mismatch >0.1 weeks for: ", paste(bad_mae$model_id, collapse = "; ")))
}
if (any(is.na(mae_check$expected_mae_week))) {
  add_warning("Some locked test models had no expected MAE entry for first-pass tolerance check.")
}

write_tsv(plot_manifest, file.path(plot_dir, "plot_manifest.tsv"))

if (length(warnings_log) == 0) warnings_log <- "No warnings."
writeLines(warnings_log, file.path(plot_dir, "plot_warnings.txt"))

summary_lines <- c(
  "# Bulk + Pseudobulk Clock Readable Plot Summary",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Inputs",
  paste0("- `", locked_path, "`"),
  paste0("- `", pred_path, "`"),
  paste0("- `", perf_path, "`"),
  "",
  "## Locked Model Selection",
  paste0(
    "- ", locked$modality_label, " / ", locked$model_family_label, ": ",
    locked$selector_label, " (validation MAE ",
    ifelse(is.na(locked$validation_mae_week), "NA", sprintf("%.2f", locked$validation_mae_week)),
    " weeks)"
  ),
  "",
  "## Locked Test Performance",
  paste0(
    "- ", locked_test_metrics$modality_label, " / ", locked_test_metrics$model_family_label,
    " / ", locked_test_metrics$selector_label, ": MAE ",
    sprintf("%.2f", locked_test_metrics$mae_week),
    " weeks, RMSE ", sprintf("%.2f", locked_test_metrics$rmse_week),
    " weeks, n=", locked_test_metrics$n
  ),
  "",
  "## Notes",
  "- Locked model plots use only the combinations listed in `locked_model_selection.tsv`.",
  "- Test performance was not used to reselect models.",
  "- Candidate test comparison is limited because first-pass outputs contain test predictions only for locked models.",
  "",
  "## Warnings",
  paste0("- ", warnings_log)
)
writeLines(summary_lines, file.path(plot_dir, "results_plot_summary.md"))

png_n <- length(list.files(plot_dir, pattern = "\\.png$", full.names = TRUE))
pdf_n <- length(list.files(plot_dir, pattern = "\\.pdf$", full.names = TRUE))
cat("\nReadable plots saved to:", plot_dir, "\n")
cat("PNG files:", png_n, "\n")
cat("PDF files:", pdf_n, "\n")
cat("Plot manifest:", file.path(plot_dir, "plot_manifest.tsv"), "\n")
cat("Summary:", file.path(plot_dir, "results_plot_summary.md"), "\n")
cat("Warnings:", file.path(plot_dir, "plot_warnings.txt"), "\n")
