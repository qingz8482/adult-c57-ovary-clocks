#!/usr/bin/env Rscript

# Reproduce copy of Supplementary Figure 6 in its original faceted layout
# (SFig6_metacell_clock_deployment_fit_stress_test_samples): one panel per cell-type metacell clock,
# faceted VCD model | Foxl2 model | Rice diet. Grey points = deployment fit of the final all-metacell
# ElasticNet model (in-sample prediction of the frozen model used for projection) with its lm fit;
# coloured symbols = external stress-test samples (mean +/- SD of their projected metacells) at their
# age at collection. Adapted from rebutal4/supple_fig6/scripts/make_sfig6.R; inputs are the
# reproduce (v2 labels/clocks) step_02-06 outputs and the reproduce GSE235021 projection.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
  library(Matrix)
})

run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")
out_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection/scatter plot")
panel_dir <- file.path(out_dir, "panels")
src_dir <- file.path(out_dir, "source_data")
for (d in c(panel_dir, src_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
gse_table_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection/GSE235021/tables")

celltype_label <- c(
  B_Plasma = "B / plasma", Endothelial = "Endothelial", Epithelial = "Epithelial",
  Granulosa = "Granulosa", Immune_Myeloid = "Immune myeloid",
  Pericyte_Smooth_Muscle = "Pericyte / smooth muscle", T_NK = "T / NK",
  Theca_Stroma = "Theca / stroma"
)

# ---------------------------------------------------------------- deployment fit
models <- readRDS(file.path(run_root, "step_04_model_building/celltype_metacell_elasticnet/celltype_metacell_elasticnet_final_models.rds"))
folds <- fread(file.path(run_root, "step_03_metacell_cv_split/metacell_fold_assignments.tsv"))
cv_eval <- fread(file.path(run_root, "step_05_internal_cv_validation/celltype_clock_cv_evaluation_summary.tsv"))
counts <- readRDS(file.path(run_root, "step_02_preprocessing_filtering/filtered_metacell_counts.rds"))
lib_size <- Matrix::colSums(counts)

# train_prediction order = cell-type subset of the fold table; verified by re-predicting
order_check <- list()
train_fit <- rbindlist(lapply(names(models), function(ct) {
  m <- models[[ct]]
  meta <- folds[harmonized_cell_class == ct]
  stopifnot(nrow(meta) == length(m$train_prediction))
  sub <- counts[m$features, meta$metacell_id, drop = FALSE]
  x <- t(as.matrix(log1p(t(t(sub) / lib_size[meta$metacell_id]) * 1e6)))
  repred <- as.numeric(predict(m$fit, newx = x))
  order_check[[ct]] <<- data.table(harmonized_cell_class = ct, n = length(repred),
                                   max_abs_diff = max(abs(repred - m$train_prediction)),
                                   cor = cor(repred, m$train_prediction))
  meta[, .(harmonized_cell_class = ct, metacell_id, sample_key, sample_id, dataset_id,
           age_weeks, predicted_age_weeks = m$train_prediction)]
}))
order_check <- rbindlist(order_check)
print(order_check)
fwrite(order_check, file.path(src_dir, "deployment_fit_order_verification.csv"))
stopifnot(all(order_check$max_abs_diff < 1e-6))

fit_stats <- train_fit[, {
  lmf <- lm(predicted_age_weeks ~ age_weeks)
  .(n_metacells = .N, n_samples = uniqueN(sample_key),
    fit_MAE_weeks = mean(abs(predicted_age_weeks - age_weeks)),
    fit_R2 = cor(predicted_age_weeks, age_weeks)^2,
    fit_slope = unname(coef(lmf)[2]), fit_intercept = unname(coef(lmf)[1]))
}, by = harmonized_cell_class]
fit_stats <- merge(fit_stats, cv_eval[, .(harmonized_cell_class, oof_cv_MAE_weeks = mae_weeks, oof_cv_R2 = r2)],
                   by = "harmonized_cell_class", all.x = TRUE)

# ---------------------------------------------------------------- external samples
prj <- fread(file.path(run_root, "step_06_projection_stress_test/prjna863443_projection_metacell_predictions.tsv"))
prj[, `:=`(
  injection_age = fcase(grepl("_3m_", sample_id), "3m", grepl("_10m_", sample_id), "10m", default = NA_character_),
  post_injection_days = fcase(grepl("_30d_", sample_id), 30L, grepl("_90d_", sample_id), 90L, default = NA_integer_)
)]
# Age at collection: VCD/vehicle = injection age + post-injection days; Foxl2 young = 4 mo, old = 9 mo
prj[, age_at_collection_weeks := fcase(
  intervention_class %in% c("vcd_control", "vcd_model") & injection_age == "3m" & post_injection_days == 30L, 17.3,
  intervention_class %in% c("vcd_control", "vcd_model") & injection_age == "3m" & post_injection_days == 90L, 25.9,
  intervention_class %in% c("vcd_control", "vcd_model") & injection_age == "10m" & post_injection_days == 30L, 47.7,
  intervention_class %in% c("vcd_control", "vcd_model") & injection_age == "10m" & post_injection_days == 90L, 56.3,
  grepl("_young_", sample_id), 17.4,
  grepl("_old_", sample_id), 39.1
)]
stopifnot(!anyNA(prj$age_at_collection_weeks))
prj[, `:=`(
  group = fcase(intervention_class == "vcd_control", "Vehicle (VCD ctrl)", intervention_class == "vcd_model", "VCD",
                intervention_class == "foxl2_control", "Foxl2 WT", intervention_class == "foxl2_model", "Foxl2+/-"),
  study_arm = fifelse(grepl("^vcd", intervention_class), "VCD model", "Foxl2 model"),
  stratum = fifelse(grepl("^vcd", intervention_class), paste0("inj ", injection_age, " + ", post_injection_days, " d"),
                    fifelse(grepl("_young_", sample_id), "4 mo", "9 mo")),
  sample_label = sample_id
)]

gse <- fread(file.path(gse_table_dir, "GSE235021_metacell_clock_predictions.tsv"))
gse[, `:=`(
  age_at_collection_weeks = age_weeks,
  group = fcase(diet == "polished_rice", "Polished rice (ctrl)", diet == "black_rice", "Black rice",
                diet == "white_brown_rice", "White-brown rice"),
  study_arm = "Rice diet",
  stratum = paste0(age_month, " mo"),
  sample_label = sample_title
)]

keep <- c("study_arm", "group", "stratum", "sample_key", "sample_label", "harmonized_cell_class",
          "metacell_id", "n_cells", "age_at_collection_weeks", "predicted_age_weeks")
ext_mc <- rbind(prj[, ..keep], gse[, ..keep])[harmonized_cell_class %in% names(models)]
stopifnot(!anyNA(ext_mc$group))
ext_sample <- ext_mc[, .(n_metacells = .N, mean_predicted_age_weeks = mean(predicted_age_weeks),
                         sd_predicted_age_weeks = sd(predicted_age_weeks)),
                     by = .(study_arm, group, stratum, sample_key, sample_label, harmonized_cell_class,
                            age_at_collection_weeks)]

arm_levels <- c("VCD model", "Foxl2 model", "Rice diet")
group_levels <- c("Vehicle (VCD ctrl)", "VCD", "Foxl2 WT", "Foxl2+/-",
                  "Polished rice (ctrl)", "Black rice", "White-brown rice")
group_cols <- c("Vehicle (VCD ctrl)" = "#6F7F8F", "VCD" = "#D1495B", "Foxl2 WT" = "#7FB3D5", "Foxl2+/-" = "#1F5FA8",
                "Polished rice (ctrl)" = "#B9A37E", "Black rice" = "#3B2A4A", "White-brown rice" = "#E39B2E")
group_shapes <- c("Vehicle (VCD ctrl)" = 21, "VCD" = 24, "Foxl2 WT" = 21, "Foxl2+/-" = 23,
                  "Polished rice (ctrl)" = 21, "Black rice" = 22, "White-brown rice" = 25)
ext_sample[, `:=`(group = factor(group, levels = group_levels), study_arm = factor(study_arm, levels = arm_levels))]

# Horizontal dodge so groups/replicates sharing an age do not overlap
ext_sample[, x_plot := {
  g <- as.integer(droplevels(group)); off_g <- (g - (max(g) + 1) / 2) * 1.6
  rep_i <- frank(sample_label, ties.method = "dense")
  off_g + (rep_i - (max(rep_i) + 1) / 2) * 0.45
} + age_at_collection_weeks, by = .(harmonized_cell_class, study_arm, round(age_at_collection_weeks))]

# Training metacells are drawn in every facet as the shared reference
train_facet <- rbindlist(lapply(arm_levels, function(a) copy(train_fit)[, study_arm := factor(a, levels = arm_levels)]))
stratum_labels <- unique(ext_sample[, .(study_arm, stratum, age_at_collection_weeks)])

fwrite(train_fit, file.path(src_dir, "SFig6_deployment_fit_training_metacells.csv"))
fwrite(fit_stats, file.path(src_dir, "SFig6_deployment_fit_statistics.csv"))
fwrite(ext_mc, file.path(src_dir, "SFig6_external_metacell_predictions.csv"))
fwrite(ext_sample[order(harmonized_cell_class, study_arm, group, stratum, sample_label)],
       file.path(src_dir, "SFig6_external_sample_predictions.csv"))

# ---------------------------------------------------------------- plotting
legend_title <- "Sample (mean ± SD of its metacell predictions)"

plot_celltype <- function(ct, base_size = 8, show_legend = TRUE) {
  tf <- train_facet[harmonized_cell_class == ct]
  es <- ext_sample[harmonized_cell_class == ct]
  st <- fit_stats[harmonized_cell_class == ct]
  rng <- range(c(tf$age_weeks, tf$predicted_age_weeks, es$mean_predicted_age_weeks + fcoalesce(es$sd_predicted_age_weeks, 0),
                 es$mean_predicted_age_weeks - fcoalesce(es$sd_predicted_age_weeks, 0)), na.rm = TRUE)
  lims <- c(floor(min(0, rng[1]) / 25) * 25, ceiling(rng[2] / 25) * 25)
  sub <- sprintf("Deployment fit: n = %d metacells / %d samples\nMAE = %.1f wk, R² = %.2f, slope = %.2f  |  OOF CV: MAE = %.1f wk, R² = %.2f",
                 st$n_metacells, st$n_samples, st$fit_MAE_weeks, st$fit_R2, st$fit_slope, st$oof_cv_MAE_weeks, st$oof_cv_R2)
  ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_point(data = tf, aes(age_weeks, predicted_age_weeks), colour = "grey72", size = 0.5, alpha = 0.55, stroke = 0) +
    geom_smooth(data = tf, aes(age_weeks, predicted_age_weeks), method = "lm", formula = y ~ x, se = TRUE,
                colour = "black", fill = "grey80", linewidth = 0.5) +
    geom_vline(data = stratum_labels, aes(xintercept = age_at_collection_weeks), colour = "grey88", linewidth = 0.3) +
    geom_text(data = stratum_labels, aes(x = age_at_collection_weeks, y = lims[1] + 1, label = stratum),
              size = base_size / 3.9, angle = 90, hjust = 0, vjust = -0.4, colour = "grey35") +
    geom_errorbar(data = es, aes(x = x_plot, ymin = mean_predicted_age_weeks - sd_predicted_age_weeks,
                                 ymax = mean_predicted_age_weeks + sd_predicted_age_weeks, colour = group),
                  width = 0, linewidth = 0.35, na.rm = TRUE) +
    geom_point(data = es, aes(x_plot, mean_predicted_age_weeks, fill = group, shape = group),
               colour = "black", size = 1.6, stroke = 0.3) +
    facet_wrap(~ study_arm, nrow = 1, drop = FALSE) +
    scale_fill_manual(values = group_cols, limits = group_levels, drop = FALSE, name = legend_title) +
    scale_colour_manual(values = group_cols, limits = group_levels, drop = FALSE, name = legend_title) +
    scale_shape_manual(values = group_shapes, limits = group_levels, drop = FALSE, name = legend_title) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE) +
    labs(title = paste0(celltype_label[[ct]], " clock"), subtitle = sub,
         x = "Chronological age (weeks)", y = "Predicted age (weeks)") +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", color = "grey70"),
          plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = base_size - 1.8, colour = "grey25", lineheight = 0.95),
          aspect.ratio = 1, legend.key = element_blank(),
          legend.position = if (show_legend) "bottom" else "none", legend.title.position = "top") +
    (if (show_legend) guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.4)))
     else guides(fill = "none", colour = "none", shape = "none"))
}

stem <- "SFig6_metacell_clock_deployment_fit_stress_test_samples"
for (ct in names(models)) {
  f <- file.path(panel_dir, paste0("SFig6_", ct, "_deployment_fit_external_samples"))
  p <- plot_celltype(ct, base_size = 9)
  ggsave(paste0(f, ".pdf"), p, width = 9, height = 4.4, device = cairo_pdf)
  ggsave(paste0(f, ".png"), p, width = 9, height = 4.4, dpi = 300)
}
combined <- wrap_plots(lapply(names(models), function(ct) plot_celltype(ct, show_legend = ct == names(models)[1])),
                       ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
ggsave(file.path(out_dir, paste0(stem, ".pdf")), combined, width = 16, height = 17, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(stem, ".png")), combined, width = 16, height = 17, dpi = 200)
writeLines(capture.output(sessionInfo()), file.path(src_dir, "sessionInfo.txt"))
message("Wrote ", file.path(out_dir, paste0(stem, ".pdf")))
