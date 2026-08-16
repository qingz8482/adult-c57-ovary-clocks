#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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

pred_path <- file.path(results_dir, "all_clock_predictions.tsv")
locked_path <- file.path(results_dir, "locked_model_selection.tsv")
stopifnot(file.exists(pred_path))
stopifnot(file.exists(locked_path))

pred <- fread(pred_path, sep = "\t")
locked <- fread(locked_path, sep = "\t")
required_cols <- c(
  "training_scope", "model_family", "selector", "fit_stage", "sample_key",
  "modality", "dataset_folder", "split_role", "age_week", "predicted_age_week"
)
missing_cols <- setdiff(required_cols, names(pred))
if (length(missing_cols) > 0) {
  stop("Missing required prediction columns: ", paste(missing_cols, collapse = ", "))
}

expected_scope <- "integrated_bulk_pseudobulk"
if (uniqueN(pred$training_scope) != 1L || unique(pred$training_scope) != expected_scope) {
  stop("Expected exactly one training_scope named ", expected_scope)
}

pred[, age_week := as.numeric(age_week)]
pred[, predicted_age_week := as.numeric(predicted_age_week)]
if (!"residual_week" %in% names(pred)) pred[, residual_week := predicted_age_week - age_week]
if (!"abs_error_week" %in% names(pred)) pred[, abs_error_week := abs(residual_week)]
recode_bio_age_bin <- function(x) {
  fifelse(
    x == "mid_adult", "mature_adult",
    fifelse(
      x == "aged", "advanced_maternal_age",
      fifelse(x == "old_aged", "post_reproductive_ageing", x)
    )
  )
}
if (!"biological_age_bin" %in% names(pred)) {
  if ("bio_age_bin" %in% names(pred)) {
    pred[, biological_age_bin := bio_age_bin]
  } else {
    pred[, biological_age_bin := fifelse(
      age_week >= 7 & age_week < 16, "young_adult",
      fifelse(
        age_week >= 16 & age_week < 40, "mature_adult",
        fifelse(age_week >= 40 & age_week < 65, "advanced_maternal_age", "post_reproductive_ageing")
      )
    )]
  }
}
pred[, biological_age_bin := recode_bio_age_bin(biological_age_bin)]
if ("bio_age_bin" %in% names(pred)) pred[, bio_age_bin := recode_bio_age_bin(bio_age_bin)]

selector_label <- function(x) {
  fifelse(
    x == "hvg_baseline", "HVG",
    fifelse(x == "version_b_age_dataset_penalized", "Version B", as.character(x))
  )
}
family_label <- function(x) {
  fifelse(x == "elasticnet", "ElasticNet", fifelse(x == "pca_glmnet", "PCA/glmnet", as.character(x)))
}
locked_map <- unique(locked[model_family %in% c("elasticnet", "pca_glmnet"), .(model_family, selector)])
setorder(locked_map, model_family, selector)
locked_map[, clock_label := paste0("Integrated ", family_label(model_family), "\n", selector_label(selector), " selector")]

plot_df <- merge(pred, locked_map, by = c("model_family", "selector"))
plot_df <- plot_df[
  (split_role == "train" & fit_stage == "train_fit") |
    (split_role == "validation" & fit_stage == "train_to_validation") |
    (split_role == "test" & fit_stage == "locked_train_validation_to_test")
]
if (nrow(plot_df) == 0) {
  stop("No frozen integrated training/validation/test predictions were found.")
}

exclude_cra003645_53w <- identical(
  toupper(Sys.getenv("EXCLUDE_CRA003645_53W", unset = "FALSE")),
  "TRUE"
)
exclusion_audit <- plot_df[0]
if (exclude_cra003645_53w) {
  exclusion_audit <- plot_df[
    split_role == "test" &
      dataset_folder == "CRA003645" &
      age_week == 53
  ]
  expected_excluded_samples <- c("Ovary1Y_R1", "Ovary1Y_R2", "Ovary1Y_R3")
  observed_excluded_samples <- sort(unique(exclusion_audit$sample_id))
  if (!identical(observed_excluded_samples, expected_excluded_samples)) {
    stop(
      "CRA003645 53-week exclusion did not resolve the three expected samples. Observed: ",
      paste(observed_excluded_samples, collapse = ", ")
    )
  }
  plot_df <- plot_df[
    !(split_role == "test" & dataset_folder == "CRA003645" & age_week == 53)
  ]
}

plot_df[, split_label := fifelse(
  split_role == "train", "training model",
  fifelse(
    split_role == "validation",
    "validation",
    fifelse(
      exclude_cra003645_53w,
      "held-out test\n(53-wk CRA excluded)",
      "held-out test"
    )
  )
)]
test_split_label <- if (exclude_cra003645_53w) {
  "held-out test\n(53-wk CRA excluded)"
} else {
  "held-out test"
}
plot_df[, split_label := factor(
  split_label,
  levels = c("training model", "validation", test_split_label)
)]
plot_df[, clock_label := factor(clock_label, levels = locked_map$clock_label)]
plot_df[, dataset_label := gsub("_", " ", dataset_folder)]
plot_df[, dataset_label := gsub("2024 NatureAging Isola aging mouse ovary atlas", "Isola ageing atlas", dataset_label, fixed = TRUE)]
plot_df[, dataset_label := gsub("2025 GeroScience Zhang multiomic aging mouse ovary", "Zhang multiomic ageing", dataset_label, fixed = TRUE)]
plot_df[, dataset_label := gsub("2026 CommBio Galligos ovarian inflammaging", "Galligos inflammaging", dataset_label, fixed = TRUE)]
plot_df[, dataset_label := gsub("CRA003645", "CRA003645 bulk timecourse", dataset_label, fixed = TRUE)]
plot_df[, dataset_label := gsub("GSE252791", "GSE252791 bulk validation", dataset_label, fixed = TRUE)]
plot_df[, biological_age_bin := factor(
  biological_age_bin,
  levels = c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
)]
dataset_levels <- sort(unique(plot_df$dataset_label))
dataset_palette <- grDevices::hcl.colors(length(dataset_levels), palette = "Dark 3")
names(dataset_palette) <- dataset_levels

fit_metric <- function(x, y, metric) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 2 || uniqueN(x[ok]) < 2 || uniqueN(y[ok]) < 2) return(NA_real_)
  fit <- lm(y[ok] ~ x[ok])
  if (metric == "slope") return(unname(coef(fit)[2]))
  if (metric == "r2") return(unname(summary(fit)$r.squared))
  if (metric == "pearson") return(unname(cor(x[ok], y[ok], method = "pearson")))
  NA_real_
}

fit_stats <- plot_df[, .(
  n = .N,
  mae_week = mean(abs(predicted_age_week - age_week), na.rm = TRUE),
  slope = fit_metric(age_week, predicted_age_week, "slope"),
  r2 = fit_metric(age_week, predicted_age_week, "r2")
), by = .(clock_label, split_label)]
fit_stats[, label := sprintf("n=%d\nMAE=%.1f\nslope=%.2f\nR2=%.2f", n, mae_week, slope, r2)]

lims <- range(c(plot_df$age_week, plot_df$predicted_age_week), na.rm = TRUE)
pad <- diff(lims) * 0.06
if (!is.finite(pad) || pad == 0) pad <- 1
lims <- c(floor(lims[1] - pad), ceiling(lims[2] + pad))

p_primary <- ggplot(plot_df, aes(x = age_week, y = predicted_age_week, color = dataset_label, shape = biological_age_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
  geom_smooth(
    data = plot_df[split_role == "train"],
    aes(x = age_week, y = predicted_age_week, group = 1),
    inherit.aes = FALSE,
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_text(
    data = fit_stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.15,
    size = 3
  ) +
  facet_grid(clock_label ~ split_label) +
  coord_cartesian(xlim = lims, ylim = lims) +
  scale_color_manual(values = dataset_palette, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "right",
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Observed age (weeks)",
    y = "Predicted age (weeks)",
    color = "Dataset",
    shape = "Biological age bin",
    title = if (exclude_cra003645_53w) {
      "Integrated bulk+pseudobulk clocks: post-hoc test sensitivity analysis"
    } else {
      "Integrated bulk+pseudobulk clocks: training fit, validation, and primary held-out test"
    },
    subtitle = if (exclude_cra003645_53w) {
      "Three 53-week CRA003645 samples excluded from the displayed test subset; models were not refitted"
    } else {
      NULL
    }
  )

stem <- file.path(out_dir, "integrated_primary_predictions")
ggsave(paste0(stem, ".svg"), p_primary, width = 15, height = 7, device = grDevices::svg, limitsize = FALSE)
ggsave(paste0(stem, ".pdf"), p_primary, width = 15, height = 7, device = grDevices::pdf, limitsize = FALSE)
ggsave(paste0(stem, ".png"), p_primary, width = 15, height = 7, dpi = 300, bg = "white", limitsize = FALSE)

fwrite(plot_df, file.path(out_dir, "integrated_primary_predictions_source.tsv"), sep = "\t")
fwrite(fit_stats, file.path(out_dir, "integrated_primary_predictions_fit_stats.tsv"), sep = "\t")
if (exclude_cra003645_53w) {
  exclusion_audit[, exclusion_reason := "post_hoc_age_restricted_sensitivity_view"]
  fwrite(
    exclusion_audit,
    file.path(out_dir, "integrated_primary_predictions_exclusion_audit.tsv"),
    sep = "\t"
  )
}

cat("Wrote integrated primary prediction plots:\n")
cat(paste0(stem, ".pdf\n"))
cat(paste0(stem, ".svg\n"))
cat(paste0(stem, ".png\n"))
