#!/usr/bin/env Rscript

# Redraw the GSE235021 six-clock diet-effect figures (6 and 9 months) in the faceted layout of
# foxl2_delta_distribution_by_celltype_combined: one free-y panel per cell-type metacell clock,
# boxes + jittered metacells per diet, coloured n labels, dashed zero line.
# Inputs are the plot-data tables written by project_GSE235021_to_metacell_clocks.R (projected
# metacells, n annotations and Welch/BH significance rows); no values are recomputed.
# Reproduce copy (v2 labels/clocks): only root_dir/run_root changed. Run from the repo root.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

root_dir <- Sys.getenv("OVARY_CLOCK_ROOT")  # reproduce: fixed root
run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")  # reproduce: v2 labels/clocks
gse_dir <- file.path(run_root, "step_08_figures_tables_release", "figures",
                     "05_stress_test_projection", "GSE235021")
table_dir <- file.path(gse_dir, "tables")

diet_levels <- c("Polished rice", "Black rice", "White/brown rice")
diet_palette <- c("Polished rice" = "#9E9E9E", "Black rice" = "#3B4CC0", "White/brown rice" = "#B85C38")
diet_outline <- c("Polished rice" = "#616161", "Black rice" = "#253494", "White/brown rice" = "#7F3B22")
# Panel order follows the Foxl2/VCD combined figures (clocks by training metacell count)
panel_levels <- c("Theca_Stroma", "Granulosa", "Immune_Myeloid", "Endothelial",
                  "Pericyte_Smooth_Muscle", "Epithelial", "T_NK")
# T_NK was excluded from the upstream six-clock tables; its metacell deltas are taken from the
# all-clock delta table and all seven clocks are re-tested here with the upstream rule (Welch t-test of
# case vs matched-control metacell deltas per contrast; BH across all 7 clocks x 4 contrasts = 28 tests).
# B_Plasma stays excluded (1-2 metacells per library; 9-month contrasts not estimable).
extra_classes <- "T_NK"
diet_code_labels <- c(polished_rice = "Polished rice", black_rice = "Black rice", white_brown_rice = "White/brown rice")

format_test_value <- function(value) {
  fcase(is.na(value), "NE", value < 0.001, "<0.001", default = sprintf("%.3f", value))
}

contrast_metacells <- fread(file = file.path(table_dir, "GSE235021_contrast_metacell_plot_source.tsv"), sep = "\t")
six_clock_tests <- fread(file = file.path(table_dir, "GSE235021_metacell_level_welch_t_tests_six_clocks.tsv"), sep = "\t")
all_metacell_delta <- fread(file = file.path(table_dir, "GSE235021_metacell_delta_to_matched_controls.tsv"), sep = "\t")

seven_clock_tests <- contrast_metacells[harmonized_cell_class %in% panel_levels, {
  case_v <- delta_predicted_age_weeks[contrast_group == "case"]
  ctrl_v <- delta_predicted_age_weeks[contrast_group == "matched_control"]
  p <- if (length(case_v) >= 2 && length(ctrl_v) >= 2) {
    tryCatch(t.test(case_v, ctrl_v, var.equal = FALSE)$p.value, error = function(e) NA_real_)
  } else NA_real_
  .(n_case_metacells = length(case_v), n_control_metacells = length(ctrl_v), welch_t_p_value = p)
}, by = .(contrast_id, harmonized_cell_class, age_month, diet_code = sub("_[0-9]+m_vs_polished$", "",
                                                                         sub("^GSE235021_", "", contrast_id)))]
stopifnot(nrow(seven_clock_tests) == 28L, all(is.finite(seven_clock_tests$welch_t_p_value)))
seven_clock_tests[, welch_t_p_value_bh := p.adjust(welch_t_p_value, method = "BH")]
seven_clock_tests[, significance_label := fcase(
  welch_t_p_value_bh < 0.001, "***", welch_t_p_value_bh < 0.01, "**", welch_t_p_value_bh < 0.05, "*",
  default = paste0("p=", format_test_value(welch_t_p_value)))]
# The six upstream clocks must reproduce the release p-values exactly (BH changes with the larger family)
chk <- merge(seven_clock_tests, six_clock_tests[, .(contrast_id, harmonized_cell_class, p_release = welch_t_p_value)],
             by = c("contrast_id", "harmonized_cell_class"))
stopifnot(nrow(chk) == 24L, max(abs(chk$welch_t_p_value - chk$p_release) / pmax(chk$p_release, 1e-300)) < 1e-8)
fwrite(seven_clock_tests, file.path(table_dir, "GSE235021_metacell_level_welch_t_tests_seven_clocks.tsv"), sep = "\t")

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

plot_age <- function(age_value) {
  stem_in <- sprintf("GSE235021_sample_delta_six_clocks_%dm_horizontal_welch_tests", age_value)
  src <- fread(file = file.path(table_dir, paste0(stem_in, "_plot_data.tsv")), sep = "\t")

  pts <- src[record_type == "projected_metacell"]
  extra <- all_metacell_delta[harmonized_cell_class %in% extra_classes & age_month == age_value,
                              .(record_type = "projected_metacell", metacell_id, sample_id,
                                diet = unname(diet_code_labels[diet]), harmonized_cell_class, delta_predicted_age_weeks)]
  stopifnot(!anyNA(extra$diet), nrow(extra) > 0)
  pts <- rbind(pts, extra, fill = TRUE)
  pts[, diet := factor(diet, levels = diet_levels)]
  pts[, panel := factor(harmonized_cell_class, levels = panel_levels)]
  stopifnot(!anyNA(pts$diet), !anyNA(pts$panel))

  ranges <- pts[, .(y_min = min(delta_predicted_age_weeks), y_max = max(delta_predicted_age_weeks)), by = panel]
  ranges[, span := pmax(y_max - y_min, 1)]

  n_lab <- pts[, .N, by = .(panel, diet)][ranges, on = "panel"]
  n_lab[, `:=`(n_label = paste0("n=", N), label_y = y_max + span * 0.06)]

  brk <- seven_clock_tests[age_month == age_value,
                           .(harmonized_cell_class, diet_code, welch_t_p_value, welch_t_p_value_bh, significance_label)]
  brk[, panel := factor(harmonized_cell_class, levels = panel_levels)]
  brk[, case_diet := factor(fifelse(diet_code == "black_rice", "Black rice", "White/brown rice"), levels = diet_levels)]
  brk <- brk[ranges, on = "panel"]
  # Polished (x=1) -> Black (x=2) lower bracket; Polished -> White/brown (x=3) upper bracket
  brk[, `:=`(
    x_start = 1,
    x_end = as.numeric(case_diet),
    bracket_y = y_max + span * fifelse(case_diet == "Black rice", 0.20, 0.36)
  )]
  brk[, `:=`(tick_y = bracket_y - span * 0.035, label_y = bracket_y + span * 0.015,
             label_x = (x_start + x_end) / 2)]

  title <- sprintf("GSE235021 %d-month diet delta-age distributions by cell-type clock", age_value)
  p <- ggplot(pts, aes(diet, delta_predicted_age_weeks, fill = diet)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_boxplot(width = 0.55, outlier.shape = NA, color = "black", linewidth = 0.26, alpha = 0.88) +
    geom_point(
      aes(color = diet),
      position = position_jitter(width = 0.12, height = 0, seed = 20260728 + age_value),
      shape = 21, size = 0.45, alpha = 0.50, stroke = 0.14
    ) +
    geom_text(data = n_lab, aes(x = diet, y = label_y, label = n_label, color = diet),
              inherit.aes = FALSE, size = 1.9, vjust = 0, show.legend = FALSE) +
    geom_segment(data = brk, aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y, color = case_diet),
                 inherit.aes = FALSE, linewidth = 0.3, show.legend = FALSE) +
    geom_segment(data = brk, aes(x = x_start, xend = x_start, y = tick_y, yend = bracket_y, color = case_diet),
                 inherit.aes = FALSE, linewidth = 0.3, show.legend = FALSE) +
    geom_segment(data = brk, aes(x = x_end, xend = x_end, y = tick_y, yend = bracket_y, color = case_diet),
                 inherit.aes = FALSE, linewidth = 0.3, show.legend = FALSE) +
    geom_text(data = brk, aes(x = label_x, y = label_y, label = significance_label, color = case_diet),
              inherit.aes = FALSE, size = 2.1, fontface = "bold", vjust = 0, show.legend = FALSE) +
    facet_wrap(~ panel, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = diet_palette, drop = FALSE) +
    scale_color_manual(values = diet_outline, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.10))) +
    theme_clock(base_size = 7.4) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 6.2),
      panel.spacing.x = grid::unit(0.34, "lines"),
      panel.spacing.y = grid::unit(0.38, "lines"),
      legend.position = "right"
    ) +
    labs(
      x = "Diet",
      y = "Metacell predicted-age delta from matched-control mean (weeks)",
      fill = "Diet",
      color = "Diet",
      title = title,
      subtitle = paste0(
        "Each point is one projected metacell; brackets compare each diet with age-matched polished-rice controls.\n",
        "Stars use BH q across 7 clocks x 4 contrasts (* <0.05, ** <0.01, *** <0.001); non-significant labels show Welch p. T_NK added; B_Plasma excluded (too few metacells)."
      )
    )

  stem_out <- sprintf("GSE235021_sample_delta_six_clocks_%dm_faceted_welch_tests", age_value)
  pdf_path <- file.path(gse_dir, paste0(stem_out, ".pdf"))
  png_path <- file.path(gse_dir, paste0(stem_out, ".png"))
  ggsave(pdf_path, p, width = 10.4, height = 6.4)
  # Rasterize the PDF so PNG and PDF share font metrics
  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    system2(pdftoppm, c("-r", "240", "-png", "-singlefile", shQuote(pdf_path), shQuote(sub("\\.png$", "", png_path))))
  } else {
    ggsave(png_path, p, width = 10.4, height = 6.4, dpi = 240, bg = "white")
  }
  out <- rbind(
    pts[, .(record_type, metacell_id, sample_id, diet = as.character(diet), harmonized_cell_class,
            delta_predicted_age_weeks)],
    n_lab[, .(record_type = "metacell_level_n_annotation", diet = as.character(diet),
              harmonized_cell_class = as.character(panel), n_metacells = N)],
    brk[, .(record_type = "metacell_level_significance_bracket", diet = as.character(case_diet),
            harmonized_cell_class, welch_t_p_value, welch_t_p_value_bh, significance_label)],
    fill = TRUE
  )
  fwrite(out, file.path(table_dir, paste0(stem_out, "_plot_data.tsv")), sep = "\t")
  message("Wrote ", pdf_path)
}

invisible(lapply(c(6, 9), plot_age))
