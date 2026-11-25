#!/usr/bin/env Rscript

# Redraw the VCD (by matched contrast) and Foxl2 (by age stratum) projection median-delta heatmaps,
# and draw the matching GSE235021 diet-contrast heatmap, with per-tile significance (Welch t-test, case vs control metacell deltas within each
# contrast/stratum and cell-type clock; from 12_add_stress_projection_pvalues.R) and half-height
# tiles. GSE235021 tests are recomputed here for all eight clocks with the same metacell-level Welch
# rule (the GSE235021 release table covers six clocks only; its p-values are reproduced exactly).
# All heatmaps use the same absolute tile (row) height: half of the original VCD row height.
# Reproduce copy (v2 labels/clocks): only root_dir/run_root changed. Run from the repo root after organize_scRNA_metacell_cv_figure_package.R and
# 12_add_stress_projection_pvalues.R.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

root_dir <- Sys.getenv("OVARY_CLOCK_ROOT")  # reproduce: fixed root
run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")  # reproduce: v2 labels/clocks
release_dir <- file.path(run_root, "step_08_figures_tables_release")
section <- "05_stress_test_projection"
plot_data_dir <- file.path(release_dir, "plot_data", section)
figure_dir <- file.path(release_dir, "figures", section)
pvalue_path <- file.path(figure_dir, "stress_projection_pvalue", "stress_projection_delta_distribution_welch_t_tests.tsv")
gse_dir <- file.path(figure_dir, "GSE235021")
gse_table_dir <- file.path(gse_dir, "tables")

# Original figure geometry from organize_scRNA_metacell_cv_figure_package.R (projection_heatmap_by_contrast)
orig_width <- 8.6
orig_height <- 5.4
tile_height_scale <- 0.5

heatmaps <- list(
  vcd = list(
    figure_id = "vcd_projection_delta_heatmap_by_contrast",
    title = "VCD projection median delta heatmap by matched contrast",
    pvalue_figure_id = "vcd_delta_distribution_by_celltype_combined",
    case_name = "VCD",
    y_label = "Matched contrast",
    loader = function(spec) load_heatmap_data(spec)
  ),
  foxl2 = list(
    figure_id = "foxl2_projection_delta_heatmap_by_age_stratum",
    title = "Foxl2 projection median delta heatmap by age stratum",
    pvalue_figure_id = "foxl2_delta_distribution_by_celltype_combined",
    case_name = "Foxl2",
    y_label = "Matched contrast",
    loader = function(spec) load_heatmap_data(spec)
  ),
  gse235021 = list(
    figure_id = "GSE235021_projection_delta_heatmap_by_contrast",
    title = "GSE235021 projection median delta heatmap by matched diet contrast",
    case_name = "diet (black or white/brown rice)",
    control_name = "age-matched polished-rice",
    y_label = "Matched contrast",
    out_figure_dir = gse_dir,
    out_plot_data_dir = gse_table_dir,
    loader = function(spec) load_gse235021_data(spec),
    positive_control = TRUE
  )
)

base_subtitle <- "Each tile summarizes case-sample deltas computed only against controls from the same matched contrast and cell type."

interpretation_shapes <- c(
  confident_projection_context = 16,
  cautious_projection_context = 17,
  diagnostic_only_projection_context = 4
)

interpretation_labels <- c(
  confident_projection_context = "Confident",
  cautious_projection_context = "Cautious",
  diagnostic_only_projection_context = "Diagnostic only"
)

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

pvals_all <- fread(file = pvalue_path, sep = "\t")

load_heatmap_data <- function(spec) {
  heatmap_path <- file.path(plot_data_dir, paste0(spec$figure_id, "_plot_data.tsv"))
  dat <- fread(file = heatmap_path, sep = "\t")
  added_cols <- c("n_case_metacells_test", "n_control_metacells_test", "p_welch_t",
                  "p_welch_t_fdr_within_figure", "p_star", "sig_label")
  dat[, intersect(added_cols, names(dat)) := NULL]  # idempotent re-runs
  contrast_order <- unique(dat$contrast_label)
  pvals <- pvals_all[
    figure_id == spec$pvalue_figure_id,
    .(contrast_label, harmonized_cell_class, n_case_metacells_test = n_case_metacells,
      n_control_metacells_test = n_control_metacells, p_welch_t, p_welch_t_fdr_within_figure, p_star)
  ]
  dat <- merge(dat, pvals, by = c("contrast_label", "harmonized_cell_class"), all.x = TRUE, sort = FALSE)
  missing <- dat[is.na(p_star)]
  if (nrow(missing) > 0) {
    warning(spec$figure_id, ": no p-value row for ",
            paste(missing$contrast_label, missing$harmonized_cell_class, sep = " / ", collapse = ", "))
  }
  # Stars for p<0.05; "NA" where the test could not run (<2 metacells in a group); blank = ns
  dat[, sig_label := fcase(
    is.na(p_star) | p_star == "NA", "NA",
    p_star == "ns", "",
    default = p_star
  )]
  # Preserve original axis ordering: cell types by median delta, contrasts in file order
  celltype_order <- dat[, .(m = median(median_delta_predicted_age_weeks, na.rm = TRUE)),
                        by = harmonized_cell_class][order(m), harmonized_cell_class]
  dat[, celltype_label := factor(harmonized_cell_class, levels = celltype_order)]
  dat[, contrast_label := factor(contrast_label, levels = contrast_order)]
  dat
}

p_to_star <- function(p) {
  fcase(is.na(p), "NA", p < 1e-4, "****", p < 1e-3, "***", p < 1e-2, "**", p < 0.05, "*", default = "ns")
}

# GSE235021: tile value = median of case-library deltas vs pooled age-matched polished-rice controls
# (same quantity as the VCD/Foxl2 median_delta_predicted_age_weeks); test = Welch t-test of case vs
# control metacell deltas, BH across all tested tiles in the figure.
load_gse235021_data <- function(spec) {
  contrast_def <- fread(file = file.path(gse_table_dir, "GSE235021_contrast_definition.tsv"), sep = "\t")
  summary_tbl <- fread(file = file.path(gse_table_dir, "GSE235021_contrast_summary_by_celltype.tsv"), sep = "\t")
  metacells <- fread(file = file.path(gse_table_dir, "GSE235021_contrast_metacell_plot_source.tsv"), sep = "\t")
  six_clock_tests <- fread(file = file.path(gse_table_dir, "GSE235021_metacell_level_welch_t_tests_six_clocks.tsv"), sep = "\t")
  # Clock quality is a property of the frozen clock, shared with the PRJNA863443 projections
  clock_quality <- unique(fread(
    file = file.path(plot_data_dir, "vcd_projection_delta_heatmap_by_contrast_plot_data.tsv"), sep = "\t"
  )[, .(harmonized_cell_class, model_status, projection_interpretation_flag, cv_MAE_weeks, cv_R2, cv_slope)])
  clock_classes <- sort(clock_quality$harmonized_cell_class)

  grid <- CJ(contrast_id = contrast_def$contrast_id, harmonized_cell_class = clock_classes)
  grid <- merge(grid, contrast_def[, .(contrast_id, contrast_label, age_month, case_diet, control_diet)], by = "contrast_id")
  case_summary <- summary_tbl[contrast_group == "case", .(
    contrast_id, harmonized_cell_class, n_case_samples = n_samples, n_case_metacells = n_metacells,
    median_delta_predicted_age_weeks = median_sample_delta_weeks, mean_delta_predicted_age_weeks = mean_sample_delta_weeks
  )]
  control_summary <- summary_tbl[contrast_group == "matched_control", .(
    contrast_id, harmonized_cell_class, n_control_samples = n_samples, n_control_metacells = n_metacells
  )]
  dat <- merge(grid, case_summary, by = c("contrast_id", "harmonized_cell_class"), all.x = TRUE)
  dat <- merge(dat, control_summary, by = c("contrast_id", "harmonized_cell_class"), all.x = TRUE)
  dat[, effect_estimable := !is.na(median_delta_predicted_age_weeks) & !is.na(n_control_metacells)]
  dat[!(effect_estimable), median_delta_predicted_age_weeks := NA_real_]

  tests <- metacells[is.finite(delta_predicted_age_weeks), {
    case_v <- delta_predicted_age_weeks[contrast_group == "case"]
    ctrl_v <- delta_predicted_age_weeks[contrast_group == "matched_control"]
    p <- if (length(case_v) >= 2 && length(ctrl_v) >= 2) {
      tryCatch(t.test(case_v, ctrl_v, var.equal = FALSE)$p.value, error = function(e) NA_real_)
    } else NA_real_
    .(n_case_metacells_test = length(case_v), n_control_metacells_test = length(ctrl_v), p_welch_t = p)
  }, by = .(contrast_id, harmonized_cell_class)]
  check <- merge(tests, six_clock_tests[, .(contrast_id, harmonized_cell_class, welch_t_p_value)],
                 by = c("contrast_id", "harmonized_cell_class"))
  if (nrow(check) != nrow(six_clock_tests) ||
      any(abs(check$p_welch_t - check$welch_t_p_value) > 1e-10 * pmax(1, check$welch_t_p_value), na.rm = TRUE)) {
    stop("Recomputed GSE235021 metacell Welch p-values do not reproduce the six-clock release table")
  }
  dat <- merge(dat, tests, by = c("contrast_id", "harmonized_cell_class"), all.x = TRUE)
  dat[!(effect_estimable), p_welch_t := NA_real_]
  dat[, p_welch_t_fdr_within_figure := p.adjust(p_welch_t, method = "BH")]
  dat[, p_star := p_to_star(p_welch_t)]
  dat <- merge(dat, clock_quality, by = "harmonized_cell_class", all.x = TRUE)
  dat[, `:=`(stress_test_family = "GSE235021_diet", inference_unit = "projected_metacell")]
  dat[, sig_label := fcase(
    !effect_estimable, "NE",
    p_star == "NA", "NA",
    p_star == "ns", "",
    default = p_star
  )]

  celltype_order <- dat[, .(m = median(median_delta_predicted_age_weeks, na.rm = TRUE)),
                        by = harmonized_cell_class][order(m), harmonized_cell_class]
  dat[, celltype_label := factor(harmonized_cell_class, levels = celltype_order)]
  dat[, contrast_label := factor(contrast_label, levels = contrast_def$contrast_label)]
  setorder(dat, contrast_label, celltype_label)
  dat
}

# GSE235021 within-dataset age positive control: each clock should place 9-month libraries ~13 weeks
# older than 6-month libraries of the same diet. Unit = library median predicted age per cell class.
# PASS = 9M > 6M in every diet AND diet-adjusted age effect one-sided p < 0.05 (lm: age + diet);
# NT = a diet lacks libraries at either age. Note: age coincides with sequencing flowcell in this
# dataset (all 9M libraries on one flowcell), so PASS shows age resolution but cannot exclude batch.
posctrl_expected_weeks <- 13
diet_row_labels <- c(polished_rice = "9 vs 6 mo: polished", white_brown_rice = "9 vs 6 mo: white/brown",
                     black_rice = "9 vs 6 mo: black")

load_gse235021_posctrl <- function(celltype_levels) {
  s <- fread(file = file.path(gse_table_dir, "GSE235021_sample_celltype_prediction_summary.tsv"), sep = "\t")
  s <- s[is.finite(median_predicted_age_weeks)]
  per_diet <- s[, {
    a9 <- median_predicted_age_weeks[age_month == 9]
    a6 <- median_predicted_age_weeks[age_month == 6]
    .(n_6m_libraries = length(a6), n_9m_libraries = length(a9),
      diff_9m_minus_6m_weeks = if (length(a6) && length(a9)) mean(a9) - mean(a6) else NA_real_,
      p_welch_t = if (length(a6) >= 2 && length(a9) >= 2) t.test(a9, a6)$p.value else NA_real_)
  }, by = .(harmonized_cell_class, diet)]
  pooled <- s[, {
    tab <- table(factor(diet, levels = names(diet_row_labels)), factor(age_month, levels = c(6, 9)))
    if (all(tab > 0)) {
      f <- lm(median_predicted_age_weeks ~ factor(age_month) + diet)
      co <- summary(f)$coefficients["factor(age_month)9", ]
      .(diff_9m_minus_6m_weeks = co[["Estimate"]],
        p_one_sided = pt(co[["t value"]], df = f$df.residual, lower.tail = FALSE))
    } else .(diff_9m_minus_6m_weeks = NA_real_, p_one_sided = NA_real_)
  }, by = harmonized_cell_class]
  verdict <- merge(
    per_diet[, .(all_diets_positive = all(diff_9m_minus_6m_weeks > 0), any_missing = anyNA(diff_9m_minus_6m_weeks)),
             by = harmonized_cell_class],
    pooled, by = "harmonized_cell_class")
  verdict[, positive_control := fcase(
    any_missing | is.na(p_one_sided), "NT",
    all_diets_positive & p_one_sided < 0.05, "PASS",
    default = "FAIL")]

  rows <- rbind(
    per_diet[, .(harmonized_cell_class, row_label = diet_row_labels[diet], value = diff_9m_minus_6m_weeks,
                 label = ifelse(is.na(diff_9m_minus_6m_weeks), "NT", sprintf("%+.1f", diff_9m_minus_6m_weeks)))],
    pooled[, .(harmonized_cell_class, row_label = "Pooled, diet-adjusted", value = diff_9m_minus_6m_weeks,
               label = ifelse(is.na(diff_9m_minus_6m_weeks), "NT",
                              sprintf("%+.1f (p=%s)", diff_9m_minus_6m_weeks, formatC(p_one_sided, digits = 2, format = "g"))))],
    verdict[, .(harmonized_cell_class, row_label = "Positive control", value = NA_real_, label = positive_control)]
  )
  rows[, celltype_label := factor(harmonized_cell_class, levels = celltype_levels)]
  rows[, row_label := factor(row_label, levels = rev(c("Positive control", "Pooled, diet-adjusted", diet_row_labels)))]
  table_out <- merge(
    dcast(per_diet, harmonized_cell_class ~ diet, value.var = c("diff_9m_minus_6m_weeks", "p_welch_t", "n_6m_libraries", "n_9m_libraries")),
    verdict[, .(harmonized_cell_class, pooled_diff_9m_minus_6m_weeks = diff_9m_minus_6m_weeks,
                pooled_p_one_sided = p_one_sided, positive_control)],
    by = "harmonized_cell_class")
  table_out[, expected_diff_weeks := posctrl_expected_weeks]
  list(rows = rows, table = table_out)
}

verdict_colors <- c(PASS = "#1B7837", FAIL = "#B2182B", NT = "grey40")

build_posctrl_plot <- function(pc) {
  r <- pc$rows
  ggplot(r, aes(celltype_label, row_label)) +
    geom_tile(aes(fill = value), color = "white", linewidth = 0.3) +
    geom_text(data = r[row_label != "Positive control"], aes(label = label), size = 2.2, color = "black") +
    geom_text(data = r[row_label == "Positive control"], aes(label = label, color = label),
              size = 2.6, fontface = "bold", show.legend = FALSE) +
    scale_fill_gradient2(low = "#762A83", mid = "white", high = "#1B7837", midpoint = 0,
                         limits = c(-25, 25), oob = scales::squish, na.value = "grey96",
                         breaks = c(-20, -10, 0, 13, 20)) +
    scale_color_manual(values = verdict_colors) +
    scale_x_discrete(drop = FALSE) +
    guides(fill = guide_colorbar(order = 3, title.vjust = 0.8, barwidth = unit(1.8, "in"), barheight = unit(0.12, "in"))) +
    theme_clock() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
          plot.subtitle = element_text(size = 8)) +
    labs(y = NULL, fill = "9M - 6M predicted\nage (weeks)",
         subtitle = sprintf("a  Age positive control: 9- vs 6-month libraries within each diet (expected about +%d weeks)",
                            posctrl_expected_weeks))
}

build_plot <- function(d, spec, annotate = TRUE) {
  if (!"effect_estimable" %in% names(d)) d <- copy(d)[, effect_estimable := TRUE]
  has_ne <- any(!d$effect_estimable)
  p <- ggplot(d, aes(celltype_label, contrast_label, fill = median_delta_predicted_age_weeks)) +
    geom_tile(color = "white", linewidth = 0.3)
  if (annotate) {
    p <- p +
      geom_point(data = d[(effect_estimable)], aes(shape = projection_interpretation_flag), color = "black", size = 1.5,
                 position = position_nudge(x = -0.28)) +
      geom_text(data = d[!sig_label %in% c("NA", "NE")], aes(label = sig_label), size = 2.6, fontface = "bold",
                color = "black", position = position_nudge(x = 0.12), vjust = 0.75) +
      geom_text(data = d[sig_label == "NA"], aes(label = sig_label), size = 2.2,
                color = "black", position = position_nudge(x = 0.12)) +
      geom_text(data = d[sig_label == "NE"], aes(label = sig_label), size = 2.2, color = "grey35")
  } else {
    p <- p + geom_point(aes(shape = projection_interpretation_flag), color = "black", size = 1.7)
  }
  p +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey88",
                         limits = spec$fill_limits) +
    scale_shape_manual(values = interpretation_shapes, labels = interpretation_labels) +
    guides(
      fill = guide_colorbar(order = 1, title.vjust = 0.8,
                            barwidth = unit(1.8, "in"), barheight = unit(0.12, "in")),
      shape = guide_legend(order = 2, title.vjust = 0.5, nrow = 1)
    ) +
    theme_clock() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7),
      legend.margin = margin(0, 12, 0, 12),
      legend.box.spacing = unit(10, "pt"),
      plot.caption = element_text(hjust = 0, size = 7, color = "grey25"),
      plot.caption.position = "plot"
    ) +
    labs(
      x = "Cell type",
      y = spec$y_label,
      fill = "Median case-control delta (weeks)",
      shape = "Projection context",
      title = spec$title,
      subtitle = base_subtitle,
      caption = if (annotate) paste0(
        "Stars: Welch t-test of ", spec$case_name, " vs ", if (is.null(spec$control_name)) "control" else spec$control_name,
        " metacell deltas within each matched contrast and cell-type clock ",
        "\n(* p<0.05, ** p<0.01, *** p<0.001, **** p<1e-4; blank = ns; NA = <2 metacells in a group",
        if (has_ne) "; NE = not estimable, no case or matched-control metacells" else "",
        ").\nRaw p-values and within-figure BH-FDR are in the plot-data table."
      ) else NULL
    )
}

panel_rows_of <- function(g) unique(g$layout$t[grepl("^panel", g$layout$name)])

fixed_height_in <- function(g) {
  rows <- panel_rows_of(g)
  sum(vapply(seq_along(g$heights), function(i) {
    if (i %in% rows) 0 else convertHeight(g$heights[i], "in", valueOnly = TRUE)
  }, numeric(1)))
}

pdf(NULL)

# Reference tile height: original VCD heatmap (8.6 x 5.4 in) row height, halved
vcd_ref <- load_heatmap_data(heatmaps$vcd)
vcd_ref_g <- ggplotGrob(
  build_plot(vcd_ref, heatmaps$vcd, annotate = FALSE) +
    guides(fill = guide_colorbar(), shape = guide_legend()) +
    theme(legend.position = "right", legend.title = element_text(), legend.text = element_text()) +
    labs(fill = "Median case-control\ndelta (weeks)", caption = NULL)
)
orig_panel_in <- orig_height - fixed_height_in(vcd_ref_g)
row_height_in <- orig_panel_in / nlevels(vcd_ref$contrast_label) * tile_height_scale

# Foxl2 uses the VCD heatmap's colour bar (same limits = data range of the VCD tiles), so colours compare
vcd_fill_limits <- range(heatmaps$vcd$loader(heatmaps$vcd)$median_delta_predicted_age_weeks, na.rm = TRUE)
heatmaps$vcd$fill_limits <- vcd_fill_limits
heatmaps$foxl2$fill_limits <- vcd_fill_limits
foxl2_range <- range(heatmaps$foxl2$loader(heatmaps$foxl2)$median_delta_predicted_age_weeks, na.rm = TRUE)
stopifnot(foxl2_range[1] >= vcd_fill_limits[1], foxl2_range[2] <= vcd_fill_limits[2])
message(sprintf("VCD/Foxl2 colour bar: %.2f to %.2f weeks", vcd_fill_limits[1], vcd_fill_limits[2]))

rendered <- lapply(heatmaps, function(spec) {
  dat <- spec$loader(spec)
  p <- build_plot(dat, spec)
  g <- ggplotGrob(p)
  panel_in <- row_height_in * nlevels(dat$contrast_label)
  g$heights[panel_rows_of(g)] <- unit(panel_in, "in")
  height_in <- fixed_height_in(g) + panel_in
  list(spec = spec, dat = dat, g = g, panel_in = panel_in, height_in = height_in)
})

# Equal tile widths across heatmaps: give every figure the widest left axis / y-title columns,
# so panels (same number of cell-type columns) get identical widths at the same page width
left_cols <- c("axis-l", "ylab-l")
col_of <- function(g, name) unique(g$layout$l[g$layout$name == name])
left_width_in <- sapply(left_cols, function(nm) max(vapply(rendered, function(r) {
  convertWidth(sum(r$g$widths[col_of(r$g, nm)]), "in", valueOnly = TRUE)
}, numeric(1))))
rendered <- lapply(rendered, function(r) {
  for (nm in left_cols) r$g$widths[col_of(r$g, nm)] <- unit(left_width_in[[nm]], "in")
  r
})

invisible(dev.off())

for (r in rendered) {
  out_figure_dir <- if (is.null(r$spec$out_figure_dir)) figure_dir else r$spec$out_figure_dir
  out_plot_data_dir <- if (is.null(r$spec$out_plot_data_dir)) plot_data_dir else r$spec$out_plot_data_dir
  pdf_path <- file.path(out_figure_dir, paste0(r$spec$figure_id, ".pdf"))
  png_path <- file.path(out_figure_dir, paste0(r$spec$figure_id, ".png"))
  if (isTRUE(r$spec$positive_control)) {
    # Stack the age positive-control panel (a) above the diet-contrast heatmap (b); same columns.
    pc <- load_gse235021_posctrl(levels(r$dat$celltype_label))
    fwrite(pc$table, file.path(out_plot_data_dir, paste0(r$spec$figure_id, "_positive_control.tsv")), sep = "\t")
    heat <- build_plot(r$dat, r$spec) +
      labs(title = NULL, subtitle = "b  Diet contrasts: median case-library delta vs age-matched polished-rice controls",
           caption = NULL) +
      theme(plot.subtitle = element_text(size = 8))
    n_pc_rows <- nlevels(pc$rows$row_label)
    combined <- (build_posctrl_plot(pc) / heat) +
      plot_layout(heights = unit(row_height_in * c(n_pc_rows, nlevels(r$dat$contrast_label)), "in"), guides = "collect") +
      plot_annotation(
        title = r$spec$title,
        caption = paste0(
          "a: library-median predicted age per cell-class clock; diet rows = mean(9M) - mean(6M); pooled = age coefficient from ",
          "lm(age + diet), one-sided p.\nPASS = older at 9M in every diet and pooled p < 0.05; FAIL otherwise; NT = a diet lacks ",
          "libraries at one age. All 9M libraries share one sequencing flowcell,\nso PASS shows age resolution but cannot exclude batch. ",
          "Diet contrasts (b) are interpretable only for PASS clocks.\n",
          "b stars: Welch t-test of ", r$spec$case_name, " vs ", r$spec$control_name,
          " metacell deltas (* p<0.05, ** p<0.01, *** p<0.001, **** p<1e-4; blank = ns;\nNA = <2 metacells in a group; ",
          "NE = not estimable). Raw p-values, BH-FDR and positive-control statistics are in the plot-data tables."),
        theme = theme(plot.title = element_text(face = "bold", size = 10.5),
                      plot.caption = element_text(hjust = 0, size = 6.6, color = "grey25"))) &
      theme(legend.position = "bottom", legend.box = "vertical", legend.box.just = "left")
    r$height_in <- r$height_in + row_height_in * n_pc_rows + 1.5
    r$g <- combined
  }
  ggsave(pdf_path, r$g, width = orig_width, height = r$height_in)
  # Rasterize the PDF so the PNG uses the same font metrics the layout was measured with
  # (the default PNG device substitutes a wider system font and can clip long axis labels)
  pdftoppm <- Sys.which("pdftoppm")
  if (nzchar(pdftoppm)) {
    png_stem <- sub("\\.png$", "", png_path)
    system2(pdftoppm, c("-r", "240", "-png", "-singlefile", shQuote(pdf_path), shQuote(png_stem)))
  } else {
    ggsave(png_path, r$g, width = orig_width, height = r$height_in, dpi = 240, bg = "white")
  }

  out_data <- copy(r$dat)
  out_data[, `:=`(celltype_label = as.character(celltype_label), contrast_label = as.character(contrast_label))]
  fwrite(out_data, file.path(out_plot_data_dir, paste0(r$spec$figure_id, "_plot_data.tsv")), sep = "\t")

  message(sprintf("%s: tile row height %.2f in, panel %.2f in, figure %.1f x %.2f in",
                  r$spec$figure_id, row_height_in, r$panel_in, orig_width, r$height_in))
  message("Wrote ", pdf_path)
}
