#!/usr/bin/env Rscript

# GSE235021 version of the stress-test composite (script 26 layout): for each cell-type metacell clock,
#   (1) deployment-fit scatter with the GSE235021 rice-diet samples (SFig6 "Rice diet" arm, script 24),
#   (2) 6-month and (3) 9-month diet delta-age boxplots, as in
#       GSE235021/GSE235021_sample_delta_six_clocks_{6,9}m_faceted_welch_tests (script 23).
# Box panels are drawn from script 23's plot-data tables (metacell deltas, n annotations and the
# Welch / BH labels), so the statistics match it. 6-month boxes keep script 23's opacity; 9-month
# boxes are half-transparent. Clocks: script 23's seven (B_Plasma excluded there, too few metacells).
# Run after scripts 23 and 24.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")
stress_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection")
sfig6_src <- file.path(stress_dir, "scatter plot", "source_data")
gse_table_dir <- file.path(stress_dir, "GSE235021", "tables")
out_dir <- file.path(stress_dir, "composite")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Order and names as in the VCD/Foxl2 composite (script 26); B_Plasma has no GSE235021 box panels
celltypes <- c("Theca_Stroma", "Granulosa", "T_NK", "Pericyte_Smooth_Muscle",
               "Immune_Myeloid", "Endothelial", "Epithelial")
scatter_title <- c(Theca_Stroma = "Theca / stroma clock", Granulosa = "Granulosa clock", T_NK = "T / NK clock",
                   Pericyte_Smooth_Muscle = "Smooth muscle clock", Immune_Myeloid = "Immune myeloid clock",
                   Endothelial = "Endothelial clock", Epithelial = "Epithelial clock")
strip_label <- c(Theca_Stroma = "Theca / Stroma", Granulosa = "Granulosa", T_NK = "T / NK",
                 Pericyte_Smooth_Muscle = "Smooth Muscle", Immune_Myeloid = "Immune / Myeloid",
                 Endothelial = "Endothelial", Epithelial = "Epithelial")
base_size <- 7.5

theme_clock <- function(base_size = 8.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.25),
      strip.text = element_text(face = "bold"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "plain")
    )
}

# ---------------------------------------------------------------- scatter (deployment fit, rice diet)
train_fit <- fread(file.path(sfig6_src, "SFig6_deployment_fit_training_metacells.csv"))
ext_sample <- fread(file.path(sfig6_src, "SFig6_external_sample_predictions.csv"))[study_arm == "Rice diet"]
group_levels <- c("Polished rice (ctrl)", "Black rice", "White-brown rice")
group_cols <- c("Polished rice (ctrl)" = "#B9A37E", "Black rice" = "#3B2A4A", "White-brown rice" = "#E39B2E")
group_shapes <- c("Polished rice (ctrl)" = 21, "Black rice" = 22, "White-brown rice" = 25)
ext_sample[, group := factor(group, levels = group_levels)]
stopifnot(!anyNA(ext_sample$group))
# Dodge as in scripts 24/25
ext_sample[, x_plot := {
  g <- as.integer(droplevels(group)); off_g <- (g - (max(g) + 1) / 2) * 1.6
  rep_i <- frank(sample_label, ties.method = "dense")
  off_g + (rep_i - (max(rep_i) + 1) / 2) * 0.45
} + age_at_collection_weeks, by = .(harmonized_cell_class, round(age_at_collection_weeks))]
stratum_labels <- unique(ext_sample[, .(stratum, age_at_collection_weeks)])[
  , .(age_at_collection_weeks = mean(age_at_collection_weeks), stratum = paste(unique(stratum), collapse = " / ")),
  by = .(age_round = round(age_at_collection_weeks))]
sample_legend <- "Sample (mean ± SD of its metacell predictions)"

plot_scatter <- function(ct) {
  tf <- train_fit[harmonized_cell_class == ct]
  es <- ext_sample[harmonized_cell_class == ct]
  rng <- range(c(tf$age_weeks, tf$predicted_age_weeks, es$mean_predicted_age_weeks + fcoalesce(es$sd_predicted_age_weeks, 0),
                 es$mean_predicted_age_weeks - fcoalesce(es$sd_predicted_age_weeks, 0)), na.rm = TRUE)
  lims <- c(floor(min(0, rng[1]) / 25) * 25, ceiling(rng[2] / 25) * 25)
  ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
    geom_point(data = tf, aes(age_weeks, predicted_age_weeks), colour = "grey72", size = 0.4, alpha = 0.55, stroke = 0) +
    geom_smooth(data = tf, aes(age_weeks, predicted_age_weeks), method = "lm", formula = y ~ x, se = TRUE,
                colour = "black", fill = "grey80", linewidth = 0.5) +
    geom_vline(data = stratum_labels, aes(xintercept = age_at_collection_weeks), colour = "grey88", linewidth = 0.3) +
    geom_text(data = stratum_labels, aes(x = age_at_collection_weeks, y = lims[1] + 1, label = stratum),
              size = base_size / 3.9, angle = 90, hjust = 0, vjust = -0.4, colour = "grey35") +
    geom_errorbar(data = es, aes(x = x_plot, ymin = mean_predicted_age_weeks - sd_predicted_age_weeks,
                                 ymax = mean_predicted_age_weeks + sd_predicted_age_weeks, colour = group),
                  width = 0, linewidth = 0.35, na.rm = TRUE) +
    geom_point(data = es, aes(x_plot, mean_predicted_age_weeks, fill = group, shape = group),
               colour = "black", size = 1.5, stroke = 0.3) +
    scale_fill_manual(values = group_cols, limits = group_levels, name = sample_legend) +
    scale_colour_manual(values = group_cols, limits = group_levels, name = sample_legend) +
    scale_shape_manual(values = group_shapes, limits = group_levels, name = sample_legend) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE) +
    labs(title = scatter_title[[ct]], x = "Chronological age (weeks)", y = "Predicted age (weeks)") +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          aspect.ratio = 1, legend.key = element_blank(), legend.title.position = "top") +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.2)))
}

# ---------------------------------------------------------------- boxplots (script 23 plot data)
diet_levels <- c("Polished rice", "Black rice", "White/brown rice")
diet_palette <- c("Polished rice" = "#9E9E9E", "Black rice" = "#3B4CC0", "White/brown rice" = "#B85C38")
diet_outline <- c("Polished rice" = "#616161", "Black rice" = "#253494", "White/brown rice" = "#7F3B22")
box_src <- lapply(c(`6` = 6, `9` = 9), function(age) {
  f <- file.path(gse_table_dir, sprintf("GSE235021_sample_delta_six_clocks_%dm_faceted_welch_tests_plot_data.tsv", age))
  fread(file = f, sep = "\t")
})

plot_box <- function(age, ct) {
  src <- box_src[[as.character(age)]][harmonized_cell_class == ct]
  pts <- src[record_type == "projected_metacell"]
  stopifnot(nrow(pts) > 0)
  pts[, `:=`(diet = factor(diet, levels = diet_levels), strip = strip_label[[ct]])]
  y_max <- max(pts$delta_predicted_age_weeks); span <- max(y_max - min(pts$delta_predicted_age_weeks), 1)
  # n labels and brackets placed as in script 23
  n_lab <- src[record_type == "metacell_level_n_annotation"][, `:=`(diet = factor(diet, levels = diet_levels),
                                                                   label_y = y_max + span * 0.06, strip = strip_label[[ct]])]
  brk <- src[record_type == "metacell_level_significance_bracket"][, `:=`(case_diet = factor(diet, levels = diet_levels),
                                                                         strip = strip_label[[ct]])]
  brk[, `:=`(x_start = 1, x_end = as.numeric(case_diet),
             bracket_y = y_max + span * fifelse(case_diet == "Black rice", 0.20, 0.36))]
  brk[, `:=`(tick_y = bracket_y - span * 0.035, label_y = bracket_y + span * 0.015, label_x = (x_start + x_end) / 2)]
  stopifnot(nrow(n_lab) == 3L, nrow(brk) == 2L)
  box_alpha <- if (age == 9) 0.5 else 0.88
  legend_name <- sprintf("%d-month diet delta age", age)
  ggplot(pts, aes(diet, delta_predicted_age_weeks, fill = diet)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = "black", linewidth = 0.26, alpha = box_alpha) +
    geom_point(aes(color = diet), position = position_jitter(width = 0.12, height = 0, seed = 20260728 + age),
               shape = 21, size = 0.4, alpha = 0.50, stroke = 0.14) +
    geom_text(data = n_lab, aes(x = diet, y = label_y, label = paste0("n=", n_metacells), color = diet),
              inherit.aes = FALSE, size = 1.6, vjust = 0, show.legend = FALSE) +
    geom_segment(data = brk, aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y, color = case_diet),
                 inherit.aes = FALSE, linewidth = 0.3, show.legend = FALSE) +
    geom_segment(data = brk, aes(x = x_start, xend = x_start, y = tick_y, yend = bracket_y, color = case_diet),
                 inherit.aes = FALSE, linewidth = 0.3, show.legend = FALSE) +
    geom_segment(data = brk, aes(x = x_end, xend = x_end, y = tick_y, yend = bracket_y, color = case_diet),
                 inherit.aes = FALSE, linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = brk, aes(x = label_x, y = label_y, label = significance_label, color = case_diet),
              inherit.aes = FALSE, size = 1.8, fontface = "bold", vjust = 0, show.legend = FALSE) +
    facet_wrap(~ strip) +
    scale_fill_manual(values = diet_palette, drop = FALSE, name = legend_name) +
    scale_color_manual(values = diet_outline, drop = FALSE, name = legend_name) +
    scale_x_discrete(labels = c("Polished rice" = "Polished", "Black rice" = "Black", "White/brown rice" = "White/brown")) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    labs(x = sprintf("%d months", age), y = "Delta age (weeks)") +
    theme_clock(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = base_size - 1),
          axis.title.x = element_text(size = base_size - 0.5))
}

# ---------------------------------------------------------------- assemble
plots <- list(); tags <- character()
tag_letters <- letters[3:(2 + 2 * length(celltypes))]   # c.. as in the VCD/Foxl2 composite
for (i in seq_along(celltypes)) {
  ct <- celltypes[[i]]
  plots <- c(plots, list(plot_scatter(ct), plot_box(6, ct), plot_box(9, ct)))
  tags <- c(tags, tag_letters[2 * i - 1], tag_letters[2 * i], "")
}
# Odd number of clocks: fill the last row's right half with empty space
if (length(celltypes) %% 2 == 1) { plots <- c(plots, rep(list(plot_spacer()), 3)) }
combined <- wrap_plots(plots, ncol = 6, widths = c(1, 0.8, 0.8, 1, 0.8, 0.8), guides = "collect") +
  plot_annotation(tag_levels = list(tags)) &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 11))

stem <- "GSE235021_clock_fit_and_delta_composite"
n_rows <- ceiling(length(celltypes) / 2)
ggsave(file.path(out_dir, paste0(stem, ".pdf")), combined, width = 16, height = 3.5 * n_rows, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(stem, ".png")), combined, width = 16, height = 3.5 * n_rows, dpi = 200, bg = "white")
fwrite(ext_sample[harmonized_cell_class %in% celltypes], file.path(out_dir, paste0(stem, "_scatter_sample_source.tsv")), sep = "\t")
message("Wrote ", file.path(out_dir, paste0(stem, ".pdf")))
