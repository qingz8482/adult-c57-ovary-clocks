#!/usr/bin/env Rscript

# Variant of Supplementary Figure 6 (24_sfig6_deployment_fit_scatter.R) with the VCD and Foxl2
# PRJNA863443 arms merged into one panel per clock: facets "VCD + Foxl2 (PRJNA863443)" | "Rice diet".
# Reads the source tables written by script 24 (no re-prediction or refitting); styling, colours,
# symbols and statistics are unchanged. Run after 24_sfig6_deployment_fit_scatter.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")
out_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection/scatter plot")
src_dir <- file.path(out_dir, "source_data")

celltype_label <- c(
  B_Plasma = "B / plasma", Endothelial = "Endothelial", Epithelial = "Epithelial",
  Granulosa = "Granulosa", Immune_Myeloid = "Immune myeloid",
  Pericyte_Smooth_Muscle = "Pericyte / smooth muscle", T_NK = "T / NK",
  Theca_Stroma = "Theca / stroma"
)

train_fit <- fread(file.path(src_dir, "SFig6_deployment_fit_training_metacells.csv"))
fit_stats <- fread(file.path(src_dir, "SFig6_deployment_fit_statistics.csv"))
ext_sample <- fread(file.path(src_dir, "SFig6_external_sample_predictions.csv"))

merged_arm <- "VCD + Foxl2 (PRJNA863443)"
arm_levels <- c(merged_arm, "Rice diet")
ext_sample[study_arm %in% c("VCD model", "Foxl2 model"), study_arm := merged_arm]
stopifnot(all(ext_sample$study_arm %in% arm_levels))

group_levels <- c("Vehicle (VCD ctrl)", "VCD", "Foxl2 WT", "Foxl2+/-",
                  "Polished rice (ctrl)", "Black rice", "White-brown rice")
group_cols <- c("Vehicle (VCD ctrl)" = "#6F7F8F", "VCD" = "#D1495B", "Foxl2 WT" = "#7FB3D5", "Foxl2+/-" = "#1F5FA8",
                "Polished rice (ctrl)" = "#B9A37E", "Black rice" = "#3B2A4A", "White-brown rice" = "#E39B2E")
group_shapes <- c("Vehicle (VCD ctrl)" = 21, "VCD" = 24, "Foxl2 WT" = 21, "Foxl2+/-" = 23,
                  "Polished rice (ctrl)" = 21, "Black rice" = 22, "White-brown rice" = 25)
ext_sample[, `:=`(group = factor(group, levels = group_levels), study_arm = factor(study_arm, levels = arm_levels))]

# Same horizontal dodge as script 24, now across all groups sharing an age within the merged panel
# (VCD inj 3m + 30 d at 17.3 wk and Foxl2 4 mo at 17.4 wk dodge as one cluster)
ext_sample[, x_plot := {
  g <- as.integer(droplevels(group)); off_g <- (g - (max(g) + 1) / 2) * 1.6
  rep_i <- frank(sample_label, ties.method = "dense")
  off_g + (rep_i - (max(rep_i) + 1) / 2) * 0.45
} + age_at_collection_weeks, by = .(harmonized_cell_class, study_arm, round(age_at_collection_weeks))]

train_facet <- rbindlist(lapply(arm_levels, function(a) copy(train_fit)[, study_arm := factor(a, levels = arm_levels)]))
# One reference line per age; strata at (nearly) the same age share one label
stratum_labels <- unique(ext_sample[, .(study_arm, stratum, age_at_collection_weeks)])
stratum_labels <- stratum_labels[, .(age_at_collection_weeks = mean(age_at_collection_weeks),
                                     stratum = paste(unique(stratum), collapse = " / ")),
                                 by = .(study_arm, age_round = round(age_at_collection_weeks))]

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

celltypes <- fit_stats$harmonized_cell_class
stem <- "SFig6_metacell_clock_deployment_fit_stress_test_samples_VCD_Foxl2_merged"
combined <- wrap_plots(lapply(celltypes, function(ct) plot_celltype(ct, show_legend = ct == celltypes[1])),
                       ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(legend.position = "bottom")
# Two facets per clock instead of three: two-thirds of the original 16 in width keeps panel size
ggsave(file.path(out_dir, paste0(stem, ".pdf")), combined, width = 11.5, height = 17, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(stem, ".png")), combined, width = 11.5, height = 17, dpi = 200)
message("Wrote ", file.path(out_dir, paste0(stem, ".pdf")))
