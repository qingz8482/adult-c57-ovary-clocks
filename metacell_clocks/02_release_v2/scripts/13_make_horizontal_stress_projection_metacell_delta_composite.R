#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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
plot_data_dir <- file.path(release_dir, "plot_data", "05_stress_test_projection")
out_dir <- file.path(
  release_dir,
  "figures",
  "05_stress_test_projection",
  "stress_projection_pvalue"
)
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

t_test_by_celltype <- function(dt, family) {
  groups <- case_control_groups(family)
  dt[, {
    vals_case <- delta_predicted_age_weeks[
      projection_group == groups$case & is.finite(delta_predicted_age_weeks)
    ]
    vals_control <- delta_predicted_age_weeks[
      projection_group == groups$control & is.finite(delta_predicted_age_weeks)
    ]
    p <- if (length(vals_case) >= 2L && length(vals_control) >= 2L) {
      tryCatch(t.test(vals_case, vals_control, var.equal = FALSE)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    .(p_plot_label = p_label(p))
  }, by = celltype_label]
}

make_count_labels <- function(dt) {
  labels <- dt[, .(
    n_metacells = .N,
    panel_y_max_group = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = .(celltype_label, projection_group)]
  ranges <- dt[, .(
    panel_y_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_y_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = celltype_label]
  ranges[, y_span := pmax(panel_y_max - panel_y_min, 1)]
  labels <- merge(labels, ranges, by = "celltype_label", all.x = TRUE)
  labels[, `:=`(
    label_y = panel_y_max + y_span * 0.05,
    n_label = paste0("n=", n_metacells)
  )]
  labels
}

make_legacy_plot <- function(dt, family, title) {
  dt <- copy(dt)
  dt[, delta_predicted_age_weeks := as.numeric(delta_predicted_age_weeks)]
  dt[, projection_group := factor(
    as.character(projection_group),
    levels = group_order_for_family(family)
  )]
  dt[, celltype_label := factor(
    celltype_label,
    levels = unique(celltype_label[order(harmonized_cell_class)])
  )]

  count_labels <- make_count_labels(dt)
  ranges <- dt[, .(
    panel_y_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_y_max = max(delta_predicted_age_weeks, na.rm = TRUE)
  ), by = celltype_label]
  ranges[, y_span := pmax(panel_y_max - panel_y_min, 1)]
  ann <- merge(t_test_by_celltype(dt, family), ranges, by = "celltype_label", all.x = TRUE)
  ann[, label_y := panel_y_max + y_span * 0.18]

  y_limits <- range(
    c(dt$delta_predicted_age_weeks, count_labels$label_y, ann$label_y),
    na.rm = TRUE
  )
  y_pad <- diff(y_limits) * 0.04

  ggplot(dt, aes(celltype_label, delta_predicted_age_weeks, fill = projection_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_boxplot(
      width = 0.42,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.30,
      alpha = 0.88,
      position = position_dodge2(width = 0.68, preserve = "single")
    ) +
    geom_point(
      aes(color = projection_group),
      position = position_jitterdodge(
        jitter.width = 0.12,
        jitter.height = 0,
        dodge.width = 0.68,
        seed = 31
      ),
      shape = 21,
      size = 0.48,
      alpha = 0.50,
      stroke = 0.18
    ) +
    geom_text(
      data = count_labels,
      aes(
        x = celltype_label,
        y = label_y,
        label = n_label,
        color = projection_group,
        group = projection_group
      ),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.68),
      size = 2.0,
      vjust = 0,
      show.legend = FALSE
    ) +
    geom_text(
      data = ann,
      aes(x = celltype_label, y = label_y, label = p_plot_label),
      inherit.aes = FALSE,
      size = 2.25,
      fontface = "bold",
      vjust = 0
    ) +
    scale_fill_manual(values = projection_group_palette, drop = TRUE) +
    scale_color_manual(values = projection_group_outline_palette, drop = TRUE) +
    scale_y_continuous(
      limits = c(y_limits[1] - y_pad, y_limits[2] + y_pad),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    theme_clock(base_size = 8.2) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 6.3),
      legend.position = "right"
    ) +
    labs(
      x = "Cell type and internal CV quality",
      y = "Predicted - observed age (weeks)",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = paste(
        "Welch t-test labels compare pooled case vs control metacell delta",
        "distributions within each cell type."
      )
    )
}

read_plot_data <- function(stem) {
  path <- file.path(plot_data_dir, paste0(stem, "_plot_data.tsv"))
  if (!file.exists(path)) stop("Missing plot-data table: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

vcd_plot <- make_legacy_plot(
  read_plot_data("vcd_metacell_delta_distribution_by_celltype"),
  "VCD",
  "VCD pooled metacell delta-age distribution with Welch t-test labels"
)
foxl2_plot <- make_legacy_plot(
  read_plot_data("foxl2_metacell_delta_distribution_by_celltype"),
  "Foxl2",
  "Foxl2 pooled metacell delta-age distribution with Welch t-test labels"
)

horizontal_composite <- vcd_plot | foxl2_plot
output_stem <- file.path(
  out_dir,
  "vcd_foxl2_metacell_delta_distribution_by_celltype_horizontal"
)

ggsave(
  paste0(output_stem, ".pdf"),
  horizontal_composite,
  width = 18.4,
  height = 7.2,
  units = "in",
  device = grDevices::cairo_pdf,
  limitsize = FALSE
)
ggsave(
  paste0(output_stem, ".png"),
  horizontal_composite,
  width = 18.4,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

message("Wrote horizontal stress-projection composite to: ", output_stem, ".{pdf,png}")
