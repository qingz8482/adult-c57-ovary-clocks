#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

root_candidates <- c(
  normalizePath(getwd(), mustWork = FALSE),
  Sys.getenv("OVARY_CLOCK_ROOT"),
  Sys.getenv("OVARY_CLOCK_ROOT")
)
root_dir <- Sys.getenv("OVARY_CLOCK_ROOT")  # reproduce: fixed root
run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")  # reproduce: v2 labels/clocks
release_dir <- file.path(run_root, "step_08_figures_tables_release")
section <- "05_stress_test_projection"
plot_data_dir <- file.path(release_dir, "plot_data", section)
figure_dir <- file.path(release_dir, "figures", section)
out_dir <- file.path(figure_dir, "stress_projection_pvalue")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

projection_group_palette <- c(
  vcd = "#F8766D",
  vcd_control = "#00BFC4",
  foxl2 = "#F8766D",
  foxl2_control = "#00BFC4"
)
projection_group_outline_palette <- c(
  vcd = "#B63A32",
  vcd_control = "#00787B",
  foxl2 = "#B63A32",
  foxl2_control = "#00787B"
)
stress_box_width <- 0.42
combined_panel_y_to_x_ratio <- 221 / 209

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

safe_id <- function(x) {
  x <- tolower(gsub("[^A-Za-z0-9]+", "_", x))
  gsub("^_|_$", "", x)
}

format_p <- function(p) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.4f", p))
  )
}

p_to_star <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(p < 1e-4, "****",
      ifelse(p < 1e-3, "***",
        ifelse(p < 1e-2, "**",
          ifelse(p < 0.05, "*", "ns")
        )
      )
    )
  )
}

p_label <- function(p) {
  star <- p_to_star(p)
  ifelse(star == "ns", paste0("p=", format_p(p)), star)
}

save_plot <- function(p, stem, width, height) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  invisible(c(pdf_path, png_path))
}

read_plot_data <- function(filename) {
  path <- file.path(plot_data_dir, filename)
  if (!file.exists(path)) stop("Missing plot-data table: ", path, call. = FALSE)
  dt <- fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
  dt[, source_plot_data := path]
  dt[, delta_predicted_age_weeks := as.numeric(delta_predicted_age_weeks)]
  dt
}

group_order_for_family <- function(family) {
  if (identical(family, "VCD")) c("vcd_control", "vcd") else c("foxl2_control", "foxl2")
}

case_control_groups <- function(family) {
  if (identical(family, "VCD")) {
    list(control = "vcd_control", case = "vcd")
  } else {
    list(control = "foxl2_control", case = "foxl2")
  }
}

t_test_one <- function(dt, family, by_cols, figure_id) {
  groups <- case_control_groups(family)
  if (length(by_cols) == 0) {
    keys <- data.table(.all = "all")
    split_dt <- list(dt)
  } else {
    keys <- unique(dt[, ..by_cols])
    split_dt <- lapply(seq_len(nrow(keys)), function(ii) {
      keep <- rep(TRUE, nrow(dt))
      for (cc in by_cols) {
        keep <- keep & as.character(dt[[cc]]) == as.character(keys[[cc]][ii])
      }
      dt[keep]
    })
  }
  out <- lapply(seq_along(split_dt), function(ii) {
    x <- split_dt[[ii]]
    vals_case <- x[projection_group == groups$case & is.finite(delta_predicted_age_weeks), delta_predicted_age_weeks]
    vals_control <- x[projection_group == groups$control & is.finite(delta_predicted_age_weeks), delta_predicted_age_weeks]
    p <- if (length(vals_case) >= 2 && length(vals_control) >= 2) {
      tryCatch(t.test(vals_case, vals_control, var.equal = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    base <- if (length(by_cols) == 0) data.table() else keys[ii, ..by_cols]
    cbind(
      data.table(
        figure_id = figure_id,
        stress_test_family = family,
        case_group = groups$case,
        control_group = groups$control,
        n_case_metacells = length(vals_case),
        n_control_metacells = length(vals_control),
        mean_case_delta_week = mean(vals_case, na.rm = TRUE),
        mean_control_delta_week = mean(vals_control, na.rm = TRUE),
        mean_difference_week = mean(vals_case, na.rm = TRUE) - mean(vals_control, na.rm = TRUE),
        median_case_delta_week = median(vals_case, na.rm = TRUE),
        median_control_delta_week = median(vals_control, na.rm = TRUE),
        p_welch_t = p,
        p_star = p_to_star(p),
        p_plot_label = p_label(p)
      ),
      base
    )
  })
  rbindlist(out, fill = TRUE)
}

make_count_labels <- function(dt, x_col, facet_cols = character()) {
  by_cols <- unique(c(facet_cols, x_col, "projection_group"))
  labels <- dt[, .(
    n_metacells = .N,
    y_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    y_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = by_cols]
  label_by <- unique(c(facet_cols, x_col))
  ranges <- dt[, .(
    panel_y_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_y_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = label_by]
  ranges[, y_span := pmax(panel_y_max - panel_y_min, 1)]
  labels <- merge(labels, ranges, by = label_by, all.x = TRUE)
  labels[, `:=`(
    label_y = panel_y_max + y_span * 0.08,
    n_label = paste0("n=", n_metacells)
  )]
  labels
}

make_annotation <- function(dt, tests, x_col, facet_cols = character()) {
  label_by <- unique(c(facet_cols, x_col))
  ranges <- dt[, .(
    panel_y_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_y_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = label_by]
  ranges[, y_span := pmax(panel_y_max - panel_y_min, 1)]
  ann <- merge(tests, ranges, by = label_by, all.x = TRUE)
  ann[, `:=`(
    label_y = panel_y_max + y_span * 0.20,
    bracket_y = panel_y_max + y_span * 0.16,
    tip_y = panel_y_max + y_span * 0.12,
    x_pos = as.numeric(factor(get(x_col), levels = unique(dt[[x_col]]))),
    x_left = as.numeric(factor(get(x_col), levels = unique(dt[[x_col]]))) - 0.18,
    x_right = as.numeric(factor(get(x_col), levels = unique(dt[[x_col]]))) + 0.18
  )]
  ann
}

plot_matched_distribution <- function(dt, family, figure_id, title, combined = FALSE) {
  dt <- copy(dt)
  dt[, projection_group := factor(as.character(projection_group), levels = group_order_for_family(family))]
  dt[, contrast_label := factor(as.character(contrast_label), levels = unique(as.character(contrast_label)))]
  if ("harmonized_cell_class" %in% names(dt)) {
    dt[, harmonized_cell_class := factor(as.character(harmonized_cell_class), levels = unique(as.character(harmonized_cell_class)))]
  }
  facet_cols <- if (combined) "harmonized_cell_class" else character()
  tests <- t_test_one(dt, family, unique(c(facet_cols, "contrast_label")), figure_id)
  count_labels <- make_count_labels(dt, "contrast_label", facet_cols)
  ann <- make_annotation(dt, tests, "contrast_label", facet_cols)

  # Half-transparent box fills for both groups in the Foxl2 combined figure
  box_alpha <- if (figure_id == "foxl2_delta_distribution_by_celltype_combined") 0.5 else 0.88
  p <- ggplot(dt, aes(contrast_label, delta_predicted_age_weeks, fill = projection_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_boxplot(
      width = stress_box_width,
      outlier.shape = NA,
      color = "black",
      linewidth = if (combined) 0.26 else 0.30,
      alpha = box_alpha,
      position = position_dodge2(width = 0.68, preserve = "single")
    ) +
    geom_point(
      aes(fill = projection_group, color = projection_group),
      position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.68, seed = if (combined) 23 else 17),
      shape = 21,
      size = if (combined) 0.34 else 0.42,
      alpha = if (combined) 0.50 else 0.55,
      stroke = if (combined) 0.14 else 0.16
    ) +
    geom_text(
      data = count_labels,
      aes(x = contrast_label, y = label_y, label = n_label, color = projection_group, group = projection_group),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.68),
      size = if (combined) 1.45 else 2.0,
      vjust = 0,
      show.legend = FALSE
    ) +
    geom_segment(
      data = ann,
      aes(x = x_left, xend = x_right, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.28
    ) +
    geom_segment(
      data = ann,
      aes(x = x_left, xend = x_left, y = tip_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.28
    ) +
    geom_segment(
      data = ann,
      aes(x = x_right, xend = x_right, y = tip_y, yend = bracket_y),
      inherit.aes = FALSE,
      linewidth = 0.28
    ) +
    geom_text(
      data = ann,
      aes(x = x_pos, y = label_y, label = p_plot_label),
      inherit.aes = FALSE,
      size = if (combined) 1.65 else 2.4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = projection_group_palette, drop = TRUE) +
    scale_color_manual(values = projection_group_outline_palette, drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.34))) +
    theme_clock(base_size = if (combined) 7.4 else 8.2) +
    theme(
      axis.text.x = element_text(angle = if (combined) 35 else 30, hjust = 1, vjust = 1, size = if (combined) 5.8 else rel(1)),
      legend.position = "right"
    ) +
    labs(
      x = if (family == "VCD") "VCD matched contrast" else "Foxl2 matched age stratum",
      y = "Delta predicted age vs matched controls (weeks)",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = "Welch t-test labels compare case vs matched-control metacell delta distributions within each displayed contrast/cell type."
    )

  if (combined) {
    p <- p +
      facet_wrap(~ harmonized_cell_class, scales = "free_y", ncol = 4) +
      theme(
        aspect.ratio = combined_panel_y_to_x_ratio,
        panel.spacing.x = grid::unit(0.34, "lines"),
        panel.spacing.y = grid::unit(0.38, "lines"),
        plot.margin = margin(5.5, 5.5, 5.5, 12, unit = "pt")
      )
  }

  list(plot = p, tests = tests)
}

make_legacy_count_labels <- function(dt) {
  labels <- dt[, .(
    n_metacells = .N,
    x_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    x_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = .(celltype_label, projection_group)]
  ranges <- dt[, .(
    panel_x_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_x_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = celltype_label]
  ranges[, x_span := pmax(panel_x_max - panel_x_min, 1)]
  labels <- merge(labels, ranges, by = "celltype_label", all.x = TRUE)
  labels[, `:=`(
    label_x = panel_x_max + x_span * 0.05,
    n_label = paste0("n=", n_metacells)
  )]
  labels
}

plot_legacy_distribution <- function(dt, family, figure_id, title) {
  dt <- copy(dt)
  dt[, projection_group := factor(as.character(projection_group), levels = group_order_for_family(family))]
  dt[, celltype_label := factor(celltype_label, levels = rev(unique(celltype_label[order(harmonized_cell_class)])))]
  tests <- t_test_one(dt, family, "celltype_label", figure_id)
  count_labels <- make_legacy_count_labels(dt)
  ranges <- dt[, .(
    panel_x_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_x_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = celltype_label]
  ranges[, x_span := pmax(panel_x_max - panel_x_min, 1)]
  ann <- merge(tests, ranges, by = "celltype_label", all.x = TRUE)
  ann[, `:=`(
    label_x = panel_x_max + x_span * 0.18,
    bracket_x = panel_x_max + x_span * 0.14,
    tip_x = panel_x_max + x_span * 0.10
  )]
  x_limits <- range(c(dt$delta_predicted_age_weeks, count_labels$label_x, ann$label_x), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04

  p <- ggplot(dt, aes(delta_predicted_age_weeks, celltype_label, fill = projection_group)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_boxplot(
      width = stress_box_width,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.30,
      alpha = 0.88,
      position = position_dodge2(width = 0.68, preserve = "single")
    ) +
    geom_point(
      aes(fill = projection_group, color = projection_group),
      position = position_jitterdodge(jitter.width = 0, jitter.height = 0.12, dodge.width = 0.68, seed = 31),
      shape = 21,
      size = 0.48,
      alpha = 0.50,
      stroke = 0.18
    ) +
    geom_text(
      data = count_labels,
      aes(x = label_x, y = celltype_label, label = n_label, color = projection_group, group = projection_group),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.68),
      size = 2.0,
      hjust = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = ann,
      aes(x = label_x, y = celltype_label, label = p_plot_label),
      inherit.aes = FALSE,
      size = 2.25,
      fontface = "bold",
      hjust = 0
    ) +
    scale_fill_manual(values = projection_group_palette, drop = TRUE) +
    scale_color_manual(values = projection_group_outline_palette, drop = TRUE) +
    scale_x_continuous(limits = c(x_limits[1] - x_pad, x_limits[2] + x_pad), expand = expansion(mult = c(0.02, 0.02))) +
    theme_clock(base_size = 8.2) +
    theme(legend.position = "right") +
    labs(
      x = "Predicted - observed age (weeks)",
      y = "Cell type and internal CV quality",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = "Welch t-test labels compare pooled case vs control metacell delta distributions within each cell type."
    )
  list(plot = p, tests = tests)
}

all_tests <- list()
append_tests <- function(x) {
  all_tests[[length(all_tests) + 1L]] <<- x
}

matched_files <- list.files(
  plot_data_dir,
  pattern = "^(vcd|foxl2)_delta_distribution_.*_plot_data\\.tsv$",
  full.names = FALSE
)
matched_files <- matched_files[!grepl("metacell_delta_distribution", matched_files)]

for (filename in sort(matched_files)) {
  dt <- read_plot_data(filename)
  family <- if (startsWith(filename, "vcd_")) "VCD" else "Foxl2"
  figure_id <- sub("_plot_data\\.tsv$", "", filename)
  stem <- figure_id
  combined <- grepl("_by_celltype_combined$", figure_id)
  celltype <- unique(as.character(dt$harmonized_cell_class))
  title <- if (combined) {
    if (family == "VCD") "VCD delta-age distribution by cell type with Welch t-test labels" else "Foxl2 delta-age distribution by cell type with Welch t-test labels"
  } else {
    paste0(family, " delta-age distribution with Welch t-test labels: ", celltype[[1]])
  }
  res <- plot_matched_distribution(dt, family, figure_id, title, combined = combined)
  save_plot(res$plot, stem, width = if (combined) 11.2 else 7.2, height = if (combined) 7.6 else 4.8)
  append_tests(res$tests)
}

legacy_files <- c(
  "vcd_metacell_delta_distribution_by_celltype_plot_data.tsv",
  "foxl2_metacell_delta_distribution_by_celltype_plot_data.tsv"
)
for (filename in legacy_files) {
  dt <- read_plot_data(filename)
  family <- if (startsWith(filename, "vcd_")) "VCD" else "Foxl2"
  figure_id <- sub("_plot_data\\.tsv$", "", filename)
  title <- if (family == "VCD") {
    "VCD pooled metacell delta-age distribution with Welch t-test labels"
  } else {
    "Foxl2 pooled metacell delta-age distribution with Welch t-test labels"
  }
  res <- plot_legacy_distribution(dt, family, figure_id, title)
  save_plot(res$plot, figure_id, width = 9.2, height = 6.6)
  append_tests(res$tests)
}

test_table <- rbindlist(all_tests, fill = TRUE)
test_table[, p_welch_t_fdr_within_figure := p.adjust(p_welch_t, method = "BH"), by = figure_id]
setcolorder(
  test_table,
  c(
    "figure_id", "stress_test_family", "harmonized_cell_class", "celltype_label",
    "contrast_label", "case_group", "control_group", "n_case_metacells",
    "n_control_metacells", "mean_case_delta_week", "mean_control_delta_week",
    "mean_difference_week", "median_case_delta_week", "median_control_delta_week",
    "p_welch_t", "p_welch_t_fdr_within_figure", "p_star", "p_plot_label"
  )[c(
    "figure_id", "stress_test_family", "harmonized_cell_class", "celltype_label",
    "contrast_label", "case_group", "control_group", "n_case_metacells",
    "n_control_metacells", "mean_case_delta_week", "mean_control_delta_week",
    "mean_difference_week", "median_case_delta_week", "median_control_delta_week",
    "p_welch_t", "p_welch_t_fdr_within_figure", "p_star", "p_plot_label"
  ) %in% names(test_table)]
)
fwrite(
  test_table,
  file.path(out_dir, "stress_projection_delta_distribution_welch_t_tests.tsv"),
  sep = "\t",
  na = "NA"
)

message("Wrote stress projection p-value plots and t-test table to: ", out_dir)
