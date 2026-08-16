#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
projection_dir <- file.path(project_dir, "results", "external_projection")
raw_dir <- file.path(projection_dir, "diagnostics", "raw delta_age by dataset_with t test")
relative_dir <- file.path(projection_dir, "diagnostics", "newly added with pvalue")
out_dir <- file.path(projection_dir, "diagnostics", "combined three-panel delta age")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

middle_dataset <- Sys.getenv("DELTA_MIDDLE_DATASET", unset = "GSE279446")
allowed_middle_datasets <- c("GSE279446", "GSE273724")
if (!middle_dataset %chin% allowed_middle_datasets) {
  stop(
    "DELTA_MIDDLE_DATASET must be one of: ",
    paste(allowed_middle_datasets, collapse = ", "),
    call. = FALSE
  )
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

format_p <- function(p) {
  fifelse(
    is.na(p),
    "p = NA",
    fifelse(p < 0.001, paste0("p = ", formatC(p, format = "e", digits = 2)),
            paste0("p = ", sprintf("%.4f", p)))
  )
}

model_labels <- c(elasticnet = "ElasticNet", pca_glmnet = "PCA/glmnet")
class_palette <- c(
  intervention = "#4E79A7",
  genetic_model = "#E15759",
  reference_control = "#2A9D8F"
)
class_labels <- c(
  intervention = "Intervention",
  genetic_model = "Genetic model",
  reference_control = "Reference control"
)

raw_predictions <- read_tsv(file.path(
  projection_dir, "simplified_classes", "external_projection_predictions_simplified.tsv"
))
if ("plot_include_primary" %in% names(raw_predictions)) {
  raw_predictions <- raw_predictions[truthy(plot_include_primary)]
}

raw_data <- raw_predictions[
  dataset_folder %chin% c("CRA017492", middle_dataset) & is.finite(delta_age_week),
  .(
    dataset_folder,
    model_family,
    sample_key,
    plot_group = projection_group,
    projection_class,
    delta_age_week = as.numeric(delta_age_week)
  )
]
raw_data <- unique(
  raw_data,
  by = c("dataset_folder", "model_family", "sample_key", "plot_group", "delta_age_week")
)

relative_data <- read_tsv(file.path(
  relative_dir, "newly_added7_relative_delta_plot_source_data.tsv"
))[
  dataset_folder == "GSE98497",
  .(
    dataset_folder,
    model_family,
    sample_key,
    plot_group,
    projection_class,
    delta_age_week = as.numeric(relative_delta_age_for_plot)
  )
]

raw_stats <- read_tsv(file.path(
  raw_dir, "raw_delta_age_condition_vs_control_tests_by_dataset_model.tsv"
))[
  dataset_folder %chin% c("CRA017492", middle_dataset) &
    n_condition > 0 & n_control > 0,
  .(
    dataset_folder,
    model_family,
    contrast_id = contrast_label,
    control_group,
    condition_group,
    n_control,
    n_condition,
    p_value = as.numeric(p_welch_t)
  )
]

relative_stats <- read_tsv(file.path(
  relative_dir, "newly_added7_relative_delta_welch_ttests.tsv"
))[
  dataset_folder == "GSE98497",
  .(
    dataset_folder,
    model_family,
    contrast_id,
    control_group = "WT",
    condition_group = "Mtor cKO",
    n_control,
    n_condition = n_experimental,
    p_value = as.numeric(p_value)
  )
]

middle_panel_specs <- list(
  GSE279446 = list(
    levels = c(
      "young WT", "young Zeb1_200H", "young Zeb1_200M",
      "adult WT", "adult Zeb1_200H", "adult Zeb1_200M"
    ),
    labels = c(
      "young WT" = "Young WT",
      "young Zeb1_200H" = "Young Zeb1-200H",
      "young Zeb1_200M" = "Young Zeb1-200M",
      "adult WT" = "Adult WT",
      "adult Zeb1_200H" = "Adult Zeb1-200H",
      "adult Zeb1_200M" = "Adult Zeb1-200M"
    ),
    title = "GSE279446",
    subtitle = "Raw delta age"
  ),
  GSE273724 = list(
    levels = c(
      "adult 27-28 week vehicle", "adult 27-28 week plus HEP14",
      "aged 73-74 week vehicle", "aged 73-74 week plus HEP14",
      "reference control"
    ),
    labels = c(
      "adult 27-28 week vehicle" = "Adult vehicle",
      "adult 27-28 week plus HEP14" = "Adult + HEP14",
      "aged 73-74 week vehicle" = "Aged vehicle",
      "aged 73-74 week plus HEP14" = "Aged + HEP14",
      "reference control" = "Reference control"
    ),
    title = "GSE273724",
    subtitle = "Raw delta age"
  )
)

panel_specs <- list(
  CRA017492 = list(
    levels = c("aged water control", "aged NMN", "young control"),
    labels = c(
      "aged water control" = "Aged water",
      "aged NMN" = "Aged NMN",
      "young control" = "Young control"
    ),
    title = "CRA017492",
    subtitle = "Raw delta age"
  ),
  GSE98497 = list(
    levels = c("WT", "Mtor cKO"),
    labels = c("WT" = "WT", "Mtor cKO" = "Mtor cKO"),
    title = "GSE98497",
    subtitle = "Matched-control relative delta age"
  )
)
panel_specs[[middle_dataset]] <- middle_panel_specs[[middle_dataset]]
panel_specs <- panel_specs[c("CRA017492", middle_dataset, "GSE98497")]

all_data <- rbindlist(list(raw_data, relative_data), use.names = TRUE)
all_stats <- rbindlist(list(raw_stats, relative_stats), use.names = TRUE)
all_data[, model_family := factor(model_family, levels = names(model_labels))]
all_stats[, model_family := factor(model_family, levels = names(model_labels))]

prepare_panel <- function(dataset_id) {
  spec <- panel_specs[[dataset_id]]
  dt <- copy(all_data[dataset_folder == dataset_id])
  st <- copy(all_stats[dataset_folder == dataset_id])
  if (!nrow(dt)) stop("No plotting rows for ", dataset_id, call. = FALSE)
  if (!nrow(st)) stop("No statistical rows for ", dataset_id, call. = FALSE)

  dt[, plot_group := factor(plot_group, levels = spec$levels)]
  positions <- setNames(seq_along(spec$levels), spec$levels)
  st[, `:=`(
    x_start = as.numeric(positions[control_group]),
    x_end = as.numeric(positions[condition_group]),
    p_label = format_p(p_value)
  )]
  if (any(!is.finite(st$x_start)) || any(!is.finite(st$x_end))) {
    stop("Unresolved comparison positions for ", dataset_id, call. = FALSE)
  }

  maxima <- dt[, .(panel_max = max(delta_age_week)), by = model_family]
  st <- merge(st, maxima, by = "model_family", all.x = TRUE)
  st[, comparison_level := seq_len(.N), by = .(model_family, control_group)]
  st[, `:=`(
    bracket_y = panel_max + 5.0 + 7.0 * (comparison_level - 1L),
    tick_y = panel_max + 4.0 + 7.0 * (comparison_level - 1L),
    label_y = panel_max + 7.0 + 7.0 * (comparison_level - 1L)
  )]
  list(data = dt, stats = st, spec = spec)
}

prepared <- lapply(names(panel_specs), prepare_panel)
names(prepared) <- names(panel_specs)
annotation_max <- max(vapply(prepared, function(x) max(x$stats$label_y), numeric(1)))
data_min <- min(all_data$delta_age_week)
y_limits <- c(
  floor((data_min - 4) / 10) * 10,
  ceiling((annotation_max + 4) / 10) * 10
)
y_breaks <- seq(y_limits[[1]], y_limits[[2]], by = 20)

theme_delta <- function() {
  theme_bw(base_size = 8.5, base_family = "sans") +
    theme(
      panel.border = element_rect(colour = "grey30", fill = NA, linewidth = 0.45),
      panel.grid.major = element_line(colour = "grey88", linewidth = 0.35),
      panel.grid.minor = element_line(colour = "grey94", linewidth = 0.25),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1, size = 7.2),
      axis.text.y = element_text(size = 7.5),
      axis.title.y = element_text(size = 8.5),
      strip.background = element_rect(fill = "grey94", colour = "grey40", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 7.8),
      plot.title = element_text(face = "bold", size = 10, margin = margin(b = 1.5)),
      plot.subtitle = element_text(size = 7.7, colour = "grey25", margin = margin(b = 4)),
      plot.margin = margin(6, 7, 5, 7),
      legend.position = "bottom",
      panel.spacing.x = unit(3.5, "mm")
    )
}

make_plot <- function(dataset_id) {
  obj <- prepared[[dataset_id]]
  dt <- obj$data
  st <- obj$stats
  spec <- obj$spec
  panel_y_limits <- if (dataset_id == "GSE98497") {
    c(-40, y_limits[[2]])
  } else {
    y_limits
  }
  panel_y_breaks <- seq(panel_y_limits[[1]], panel_y_limits[[2]], by = 20)
  legend_data <- data.table(
    plot_group = factor(rep(spec$levels[[1]], length(class_palette)), levels = spec$levels),
    delta_age_week = rep(panel_y_limits[[1]], length(class_palette)),
    projection_class = factor(names(class_palette), levels = names(class_palette))
  )

  ggplot(dt, aes(plot_group, delta_age_week, colour = projection_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.35) +
    geom_boxplot(
      outlier.shape = NA, width = 0.62, linewidth = 0.45, alpha = 0.10,
      show.legend = FALSE
    ) +
    geom_jitter(
      width = 0.13, height = 0, size = 1.55, alpha = 0.82,
      show.legend = FALSE
    ) +
    geom_point(
      data = legend_data,
      aes(plot_group, delta_age_week, colour = projection_class),
      inherit.aes = FALSE,
      alpha = 0,
      show.legend = TRUE
    ) +
    geom_segment(
      data = st,
      aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.38
    ) +
    geom_segment(
      data = st,
      aes(x = x_start, xend = x_start, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.38
    ) +
    geom_segment(
      data = st,
      aes(x = x_end, xend = x_end, y = tick_y, yend = bracket_y),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.38
    ) +
    geom_text(
      data = st,
      aes(x = (x_start + x_end) / 2, y = label_y, label = p_label),
      inherit.aes = FALSE,
      colour = "black",
      size = 2.35,
      fontface = "bold"
    ) +
    facet_wrap(
      ~model_family,
      nrow = 1,
      labeller = as_labeller(model_labels)
    ) +
    scale_x_discrete(labels = spec$labels, drop = FALSE) +
    scale_y_continuous(
      limits = panel_y_limits,
      breaks = panel_y_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    scale_colour_manual(
      values = class_palette,
      limits = names(class_palette),
      breaks = names(class_palette),
      labels = class_labels,
      drop = FALSE
    ) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.3))) +
    labs(
      title = spec$title,
      subtitle = spec$subtitle,
      x = NULL,
      y = "Delta age (weeks)",
      colour = NULL
    ) +
    theme_delta()
}

set.seed(20260811)
p_cra <- make_plot("CRA017492")
p_middle <- make_plot(middle_dataset)
p_gse98497 <- make_plot("GSE98497")

caption <- paste(
  "Brackets show unadjusted two-sided Welch t-tests evaluated separately within each model family.",
  paste0("CRA017492 and ", middle_dataset, " show raw delta age; GSE98497 shows matched-control relative delta age.")
)

combined <- (p_cra | p_middle | p_gse98497) +
  plot_layout(widths = c(1.05, 1.75, 0.95), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    caption = caption,
    theme = theme(
      plot.tag = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 7.2, colour = "grey30", hjust = 0),
      plot.margin = margin(5, 7, 5, 7)
    )
  ) &
  theme(legend.position = "bottom")

output_prefix <- paste0("CRA017492_", middle_dataset, "_GSE98497")
stem <- file.path(out_dir, paste0(output_prefix, "_delta_age_shared_y"))
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 18, height = 5.8, family = "sans")
print(combined)
dev.off()

ragg::agg_png(
  paste0(stem, ".png"),
  width = 18,
  height = 5.8,
  units = "in",
  res = 300,
  background = "white"
)
print(combined)
dev.off()

source_export <- rbindlist(lapply(names(prepared), function(dataset_id) {
  dt <- copy(prepared[[dataset_id]]$data)
  dt[, delta_definition := if (dataset_id == "GSE98497") {
    "matched-control relative delta age"
  } else {
    "raw delta age"
  }]
  dt
}), use.names = TRUE)
fwrite(
  source_export,
  file.path(out_dir, paste0(output_prefix, "_delta_age_source_data.tsv")),
  sep = "\t",
  na = "NA"
)
fwrite(
  all_stats,
  file.path(out_dir, paste0(output_prefix, "_delta_age_statistics.tsv")),
  sep = "\t",
  na = "NA"
)
fwrite(
  data.table(
    dataset_folder = c("CRA017492", middle_dataset, "GSE98497"),
    y_min_week = c(y_limits[[1]], y_limits[[1]], -40),
    y_max_week = rep(y_limits[[2]], 3),
    y_break_step_week = rep(20, 3)
  ),
  file.path(out_dir, paste0(output_prefix, "_axis_contract.tsv")),
  sep = "\t"
)

message("Wrote combined three-panel figure and source tables to: ", out_dir)
