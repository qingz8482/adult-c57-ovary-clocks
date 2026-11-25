#!/usr/bin/env Rscript

# Composite figure (layout of the user's mock-up, panels c-r): for each cell-type metacell clock,
#   (1) deployment-fit scatter with the PRJNA863443 VCD + Foxl2 samples in one panel
#       (as SFig6_..._VCD_Foxl2_merged, script 25; rice diet omitted),
#   (2) VCD delta-age boxplots by matched contrast and (3) Foxl2 delta-age boxplots by age stratum
#       (as stress_projection_pvalue/{vcd,foxl2}_delta_distribution_by_celltype_combined, script 12).
# Boxplots reuse script 12's own helpers (Welch tests, n labels, brackets), so statistics match it.
# Inputs: SFig6 source tables (script 24) and the reproduce step_08 plot data. Run after 12, 24, 25.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")
script12 <- file.path(run_root, "scripts", "12_add_stress_projection_pvalues.R")
stress_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection")
sfig6_src <- file.path(stress_dir, "scatter plot", "source_data")

# Script 12 definitions only (constants, helpers, plot_matched_distribution); stop before its run loop
s12 <- readLines(script12)
s12_end <- grep("^all_tests <- list\\(\\)", s12) - 1L
stopifnot(length(s12_end) == 1L)
eval(parse(text = s12[seq_len(s12_end)]))
# Defined after script 12's block, which sets its own out_dir
out_dir <- file.path(stress_dir, "composite")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Order and names follow the mock-up
celltypes <- c("Theca_Stroma", "Granulosa", "T_NK", "Pericyte_Smooth_Muscle",
               "Immune_Myeloid", "Endothelial", "B_Plasma", "Epithelial")
scatter_title <- c(Theca_Stroma = "Theca / stroma clock", Granulosa = "Granulosa clock", T_NK = "T / NK clock",
                   Pericyte_Smooth_Muscle = "Smooth muscle clock", Immune_Myeloid = "Immune myeloid clock",
                   Endothelial = "Endothelial clock", B_Plasma = "B cells clock", Epithelial = "Epithelial clock")
strip_label <- c(Theca_Stroma = "Theca / Stroma", Granulosa = "Granulosa", T_NK = "T / NK",
                 Pericyte_Smooth_Muscle = "Smooth Muscle", Immune_Myeloid = "Immune / Myeloid",
                 Endothelial = "Endothelial", B_Plasma = "B Plasma", Epithelial = "Epithelial")
base_size <- 7.5

# ---------------------------------------------------------------- scatter (deployment fit)
train_fit <- fread(file.path(sfig6_src, "SFig6_deployment_fit_training_metacells.csv"))
ext_sample <- fread(file.path(sfig6_src, "SFig6_external_sample_predictions.csv"))[study_arm %in% c("VCD model", "Foxl2 model")]
group_levels <- c("Vehicle (VCD ctrl)", "VCD", "Foxl2 WT", "Foxl2+/-")
group_cols <- c("Vehicle (VCD ctrl)" = "#6F7F8F", "VCD" = "#D1495B", "Foxl2 WT" = "#7FB3D5", "Foxl2+/-" = "#1F5FA8")
group_shapes <- c("Vehicle (VCD ctrl)" = 21, "VCD" = 24, "Foxl2 WT" = 21, "Foxl2+/-" = 23)
ext_sample[, group := factor(group, levels = group_levels)]
stopifnot(!anyNA(ext_sample$group))
# Dodge as in scripts 24/25, across VCD and Foxl2 groups sharing an age
ext_sample[, x_plot := {
  g <- as.integer(droplevels(group)); off_g <- (g - (max(g) + 1) / 2) * 1.6
  rep_i <- frank(sample_label, ties.method = "dense")
  off_g + (rep_i - (max(rep_i) + 1) / 2) * 0.45
} + age_at_collection_weeks, by = .(harmonized_cell_class, round(age_at_collection_weeks))]
stratum_labels <- unique(ext_sample[, .(stratum, age_at_collection_weeks)])[
  , .(age_at_collection_weeks = mean(age_at_collection_weeks), stratum = paste(unique(stratum), collapse = " / ")),
  by = .(age_round = round(age_at_collection_weeks))]

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
    scale_fill_manual(values = group_cols, limits = group_levels, name = "Sample (mean ± SD of its metacell predictions)") +
    scale_colour_manual(values = group_cols, limits = group_levels, name = "Sample (mean ± SD of its metacell predictions)") +
    scale_shape_manual(values = group_shapes, limits = group_levels, name = "Sample (mean ± SD of its metacell predictions)") +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE) +
    labs(title = scatter_title[[ct]], x = "Chronological age (weeks)", y = "Predicted age (weeks)") +
    theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          aspect.ratio = 1, legend.key = element_blank(), legend.title.position = "top") +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.2)))
}

# ---------------------------------------------------------------- boxplots (script 12)
box_data <- function(family) {
  dt <- read_plot_data(paste0(tolower(family), "_delta_distribution_by_celltype_combined_plot_data.tsv"))
  if (family == "VCD") {
    dt[, contrast_label := sub("^Injection (\\w+); post (\\w+)$", "\\1 \\2", contrast_label)]
  } else {
    dt[, contrast_label := sub("^Age ", "", contrast_label)]
  }
  dt
}
vcd_dt <- box_data("VCD")
foxl2_dt <- box_data("Foxl2")
stopifnot(all(celltypes %in% vcd_dt$harmonized_cell_class), all(celltypes %in% foxl2_dt$harmonized_cell_class))

all_box_tests <- list()
plot_box <- function(dt_all, family, ct) {
  dt <- dt_all[harmonized_cell_class == ct]
  # Strip text = display name (the facet variable of script 12's combined layout)
  dt[, harmonized_cell_class := strip_label[[ct]]]
  res <- plot_matched_distribution(dt, family, figure_id = paste0("composite_", family, "_", ct), title = NULL, combined = TRUE)
  all_box_tests[[length(all_box_tests) + 1L]] <<- res$tests[, `:=`(source_cell_class = ct)]
  p <- res$plot
  # VCD boxes keep script 12's original opacity; Foxl2 boxes are half-transparent to set them apart
  p$layers[[2]]$aes_params$alpha <- if (family == "Foxl2") 0.5 else 0.88
  case_lab <- c(vcd_control = "Vehicle control", vcd = "VCD",
                foxl2_control = "Foxl2 WT", foxl2 = "Foxl2+/-")
  legend_name <- paste(family, "delta age")
  p +
    scale_fill_manual(values = projection_group_palette, labels = case_lab, name = legend_name) +
    scale_color_manual(values = projection_group_outline_palette, labels = case_lab, name = legend_name) +
    labs(title = NULL, subtitle = NULL, x = NULL, y = "Delta age (weeks)") +
    theme(aspect.ratio = NULL, axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = base_size - 1),
          plot.margin = margin(5.5, 5.5, 5.5, 5.5))
}

# ---------------------------------------------------------------- assemble
plots <- list(); tags <- character()
letters_used <- letters[3:18]   # c..r as in the mock-up
for (i in seq_along(celltypes)) {
  ct <- celltypes[[i]]
  plots <- c(plots, list(plot_scatter(ct), plot_box(vcd_dt, "VCD", ct), plot_box(foxl2_dt, "Foxl2", ct)))
  tags <- c(tags, letters_used[2 * i - 1], letters_used[2 * i], "")
}
combined <- wrap_plots(plots, ncol = 6, widths = c(1, 1.05, 0.7, 1, 1.05, 0.7), guides = "collect") +
  plot_annotation(tag_levels = list(tags)) &
  theme(legend.position = "bottom", plot.tag = element_text(face = "bold", size = 11))

stem <- "stress_test_clock_fit_and_delta_composite"
ggsave(file.path(out_dir, paste0(stem, ".pdf")), combined, width = 16, height = 14, device = cairo_pdf)
ggsave(file.path(out_dir, paste0(stem, ".png")), combined, width = 16, height = 14, dpi = 200, bg = "white")

# Source data: scatter samples and box-panel Welch tests
fwrite(ext_sample[harmonized_cell_class %in% celltypes], file.path(out_dir, paste0(stem, "_scatter_sample_source.tsv")), sep = "\t")
fwrite(rbindlist(all_box_tests, fill = TRUE), file.path(out_dir, paste0(stem, "_box_welch_tests.tsv")), sep = "\t")
message("Wrote ", file.path(out_dir, paste0(stem, ".pdf")))
