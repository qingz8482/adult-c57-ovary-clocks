#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  supplied_path <- sub("^--file=", "", script_arg)
  supplied_path <- gsub("~\\+~", " ", supplied_path)
  normalizePath(supplied_path, mustWork = TRUE)
} else {
  normalizePath("build_dataset_adjusted_composition_plots.R", mustWork = TRUE)
}

out_dir <- dirname(script_path)
source_dir <- dirname(out_dir)
individual_dir <- file.path(out_dir, "individual_celltypes")
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

composition_path <- file.path(source_dir, "sample_celltype_composition_by_age.tsv")
reference_statistics_path <- file.path(source_dir, "celltype_age_association_statistics.tsv")
if (!file.exists(composition_path) || !file.exists(reference_statistics_path)) {
  stop("Missing source composition or statistics table", call. = FALSE)
}

composition <- fread(composition_path) |> as_tibble()
reference_statistics <- fread(reference_statistics_path) |> as_tibble()

required_columns <- c(
  "dataset_id", "sample_key", "sample_id", "age_weeks",
  "harmonized_cell_class", "cell_fraction", "celltype_label", "dataset_label"
)
missing_columns <- setdiff(required_columns, names(composition))
if (length(missing_columns) > 0L) {
  stop("Composition table is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

celltype_levels <- c(
  "B_Plasma", "Endothelial", "Epithelial", "Granulosa", "Immune_Myeloid",
  "Oocyte_Germline", "Pericyte_Smooth_Muscle", "T_NK", "Theca_Stroma"
)
celltype_labels <- c(
  B_Plasma = "B / plasma cells",
  Endothelial = "Endothelial",
  Epithelial = "Epithelial",
  Granulosa = "Granulosa",
  Immune_Myeloid = "Myeloid immune",
  Oocyte_Germline = "Oocyte / germline",
  Pericyte_Smooth_Muscle = "Pericyte / smooth muscle",
  T_NK = "T / NK cells",
  Theca_Stroma = "Theca / stroma"
)

dataset_levels <- c(
  "2022_eLife_Morris_cycling_murine_ovary_atlas",
  "2024_Cell_Winkler_FRT_scRNA_spatial",
  "2024_LifeMedicine_Jiang_whole_life_mouse_ovary",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary",
  "2025_CommBio_Zhao_developing_ovary_scRNA_spatial",
  "2026_CommBio_Galligos_ovarian_inflammaging",
  "2025_Kim_PRJNA863443_mouse_menopause_models"
)
dataset_labels <- c(
  "2022_eLife_Morris_cycling_murine_ovary_atlas" = "Morris 2022",
  "2024_Cell_Winkler_FRT_scRNA_spatial" = "Winkler 2024",
  "2024_LifeMedicine_Jiang_whole_life_mouse_ovary" = "Jiang 2024",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "Isola 2024",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary" = "Zhang 2025",
  "2025_CommBio_Zhao_developing_ovary_scRNA_spatial" = "Zhao 2025",
  "2026_CommBio_Galligos_ovarian_inflammaging" = "Galligos 2026",
  "2025_Kim_PRJNA863443_mouse_menopause_models" = "Kim 2025"
)
dataset_palette <- c(
  "Morris 2022" = "#4477AA",
  "Winkler 2024" = "#EE6677",
  "Jiang 2024" = "#228833",
  "Isola 2024" = "#AA3377",
  "Zhang 2025" = "#66CCEE",
  "Zhao 2025" = "#CCBB44",
  "Galligos 2026" = "#EE8866",
  "Kim 2025" = "#000000"
)

composition <- composition |>
  mutate(
    age_weeks = as.numeric(age_weeks),
    cell_fraction = as.numeric(cell_fraction),
    harmonized_cell_class = factor(harmonized_cell_class, levels = celltype_levels),
    celltype_label = factor(
      unname(celltype_labels[as.character(harmonized_cell_class)]),
      levels = unname(celltype_labels[celltype_levels])
    ),
    dataset_id = factor(dataset_id, levels = dataset_levels),
    dataset_label = factor(
      unname(dataset_labels[as.character(dataset_id)]),
      levels = unname(dataset_labels[dataset_levels])
    )
  )

if (anyNA(composition$harmonized_cell_class) || anyNA(composition$dataset_id)) {
  stop("Composition table contains an unexpected cell class or dataset", call. = FALSE)
}
if (anyNA(composition$age_weeks) || anyNA(composition$cell_fraction)) {
  stop("Composition table contains missing age or fraction values", call. = FALSE)
}

# For each cell class, fit the same inferential model used in the source analysis.
# The plotted line is an equal-dataset-weighted marginal prediction: at every age,
# model-matrix rows are created for all eight dataset fixed effects and averaged.
fit_adjusted_model <- function(df) {
  model_df <- df |>
    mutate(
      age_10_weeks = age_weeks / 10,
      dataset_id = factor(dataset_id, levels = dataset_levels)
    )
  fit <- lm(cell_fraction ~ age_10_weeks + dataset_id, data = model_df)
  age_grid <- seq(min(model_df$age_weeks), max(model_df$age_weeks), length.out = 240)
  terms_no_response <- delete.response(terms(fit))
  beta <- coef(fit)
  covariance <- vcov(fit)
  critical_value <- qt(0.975, df = df.residual(fit))

  prediction_rows <- lapply(age_grid, function(age_value) {
    newdata <- tibble(
      age_10_weeks = age_value / 10,
      dataset_id = factor(dataset_levels, levels = dataset_levels)
    )
    design <- model.matrix(
      terms_no_response,
      data = newdata,
      contrasts.arg = fit$contrasts
    )
    design <- design[, names(beta), drop = FALSE]
    marginal_design <- colMeans(design)
    estimate <- as.numeric(marginal_design %*% beta)
    standard_error <- sqrt(as.numeric(
      t(marginal_design) %*% covariance %*% marginal_design
    ))
    tibble(
      age_weeks = age_value,
      adjusted_fraction = estimate,
      adjusted_ci_low = estimate - critical_value * standard_error,
      adjusted_ci_high = estimate + critical_value * standard_error
    )
  }) |> bind_rows()

  age_row <- summary(fit)$coefficients["age_10_weeks", ]
  age_ci <- confint(fit, "age_10_weeks", level = 0.95)
  model_summary <- tibble(
    n_samples = nrow(model_df),
    n_datasets = n_distinct(model_df$dataset_id),
    age_min_weeks = min(model_df$age_weeks),
    age_max_weeks = max(model_df$age_weeks),
    adjusted_fraction_change_per_10_weeks = unname(age_row["Estimate"]),
    adjusted_fraction_change_ci_low = unname(age_ci[1]),
    adjusted_fraction_change_ci_high = unname(age_ci[2]),
    adjusted_age_p = unname(age_row["Pr(>|t|)"]),
    adjusted_model_r_squared = summary(fit)$r.squared
  )

  list(predictions = prediction_rows, statistics = model_summary)
}

model_results <- lapply(celltype_levels, function(celltype_id) {
  df <- composition |>
    filter(as.character(harmonized_cell_class) == celltype_id)
  result <- fit_adjusted_model(df)
  result$predictions <- result$predictions |>
    mutate(harmonized_cell_class = celltype_id)
  result$statistics <- result$statistics |>
    mutate(harmonized_cell_class = celltype_id)
  result
})

predictions <- bind_rows(lapply(model_results, `[[`, "predictions")) |>
  mutate(
    harmonized_cell_class = factor(harmonized_cell_class, levels = celltype_levels),
    celltype_label = factor(
      unname(celltype_labels[as.character(harmonized_cell_class)]),
      levels = unname(celltype_labels[celltype_levels])
    )
  )

statistics <- bind_rows(lapply(model_results, `[[`, "statistics")) |>
  mutate(
    adjusted_age_q_bh = p.adjust(adjusted_age_p, method = "BH"),
    adjusted_percentage_point_change_per_10_weeks = 100 * adjusted_fraction_change_per_10_weeks,
    adjusted_percentage_point_ci_low = 100 * adjusted_fraction_change_ci_low,
    adjusted_percentage_point_ci_high = 100 * adjusted_fraction_change_ci_high,
    celltype_label = unname(celltype_labels[harmonized_cell_class])
  ) |>
  arrange(match(harmonized_cell_class, celltype_levels))

# Verify that the refit exactly reproduces the source adjusted coefficients.
coefficient_check <- statistics |>
  select(harmonized_cell_class, refit = adjusted_fraction_change_per_10_weeks) |>
  left_join(
    reference_statistics |>
      select(harmonized_cell_class, reference = adjusted_fraction_change_per_10_weeks),
    by = "harmonized_cell_class"
  ) |>
  mutate(absolute_difference = abs(refit - reference))

if (anyNA(coefficient_check$reference) || max(coefficient_check$absolute_difference) > 1e-10) {
  stop("Refitted coefficients do not reproduce the source statistics", call. = FALSE)
}
theca_slope <- statistics |>
  filter(harmonized_cell_class == "Theca_Stroma") |>
  pull(adjusted_percentage_point_change_per_10_weeks)
if (length(theca_slope) != 1L || theca_slope >= 0) {
  stop("Theca / stroma adjusted slope is not negative", call. = FALSE)
}
if (any(!is.finite(unlist(predictions[c("adjusted_fraction", "adjusted_ci_low", "adjusted_ci_high")])))) {
  stop("Adjusted prediction table contains non-finite values", call. = FALSE)
}

fwrite(
  predictions |>
    mutate(
      harmonized_cell_class = as.character(harmonized_cell_class),
      celltype_label = as.character(celltype_label)
    ),
  file.path(out_dir, "dataset_adjusted_marginal_predictions.tsv"),
  sep = "\t"
)
fwrite(statistics, file.path(out_dir, "dataset_adjusted_model_statistics.tsv"), sep = "\t")
fwrite(coefficient_check, file.path(out_dir, "dataset_adjusted_coefficient_hard_check.tsv"), sep = "\t")

theme_composition <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.3, colour = "black"),
      plot.title = element_text(size = base_size + 0.8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#333333"),
      plot.caption = element_text(size = base_size - 1, colour = "#555555", hjust = 0),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", hjust = 0),
      legend.position = "bottom",
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.8),
      legend.key.height = grid::unit(2.2, "mm"),
      legend.key.width = grid::unit(3.2, "mm"),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

save_plot_bundle <- function(plot, stem, width_mm, height_mm, make_tiff = FALSE) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::svg(
    paste0(stem, ".svg"), width = width_in, height = height_in,
    family = "sans", pointsize = 7
  )
  print(plot)
  grDevices::dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = "sans", pointsize = 7
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()

  if (make_tiff) {
    ragg::agg_tiff(
      paste0(stem, ".tiff"), width = width_in, height = height_in,
      units = "in", res = 600, compression = "lzw", background = "white"
    )
    print(plot)
    grDevices::dev.off()
  }
}

dataset_guide <- guide_legend(
  title = "Dataset", nrow = 2, byrow = TRUE,
  override.aes = list(size = 2.4, alpha = 1)
)

wrap_text <- function(x, width) {
  paste(strwrap(x, width = width), collapse = "\n")
}

format_q <- function(x) {
  ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x))
}

make_single_plot <- function(celltype_id) {
  df <- composition |>
    filter(as.character(harmonized_cell_class) == celltype_id)
  pred <- predictions |>
    filter(as.character(harmonized_cell_class) == celltype_id)
  stat <- statistics |>
    filter(harmonized_cell_class == celltype_id)
  subtitle <- wrap_text(sprintf(
    "Dataset-adjusted change = %+.2f percentage points per 10 weeks (95%% CI, %+.2f to %+.2f; BH q = %s)",
    stat$adjusted_percentage_point_change_per_10_weeks,
    stat$adjusted_percentage_point_ci_low,
    stat$adjusted_percentage_point_ci_high,
    format_q(stat$adjusted_age_q_bh)
  ), width = 88)

  ggplot() +
    geom_ribbon(
      data = pred,
      aes(age_weeks, ymin = adjusted_ci_low, ymax = adjusted_ci_high),
      fill = "#9ECAE1", alpha = 0.32
    ) +
    geom_line(
      data = pred,
      aes(age_weeks, adjusted_fraction),
      colour = "#2F6B9A", linewidth = 0.68
    ) +
    geom_point(
      data = df,
      aes(age_weeks, cell_fraction, fill = dataset_label),
      shape = 21, size = 2.25, stroke = 0.3, colour = "white", alpha = 0.9,
      position = position_jitter(width = 0.55, height = 0, seed = 20260714)
    ) +
    scale_fill_manual(values = dataset_palette, drop = FALSE, guide = dataset_guide) +
    scale_x_continuous(
      breaks = seq(0, 90, 15), limits = c(4, 90),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 0.1),
      expand = expansion(mult = c(0.03, 0.05))
    ) +
    labs(
      x = "Chronological age (weeks)",
      y = "Fraction of confident cells",
      title = unname(celltype_labels[celltype_id]),
      subtitle = subtitle,
      caption = wrap_text(paste0(
        "Each point is one biological sample (n = ", nrow(df), "). ",
        "Line and 95% CI: equal-dataset-weighted marginal prediction from ",
        "cell fraction ~ age per 10 weeks + dataset."
      ), width = 82)
    ) +
    theme_composition()
}

single_plots <- setNames(lapply(celltype_levels, make_single_plot), celltype_levels)
for (celltype_id in celltype_levels) {
  stem <- file.path(
    individual_dir,
    paste0("dataset_adjusted_cell_fraction_vs_age_", tolower(gsub("_", "-", celltype_id)))
  )
  save_plot_bundle(single_plots[[celltype_id]], stem, width_mm = 120, height_mm = 92)
}

overview_plot <- ggplot() +
  geom_ribbon(
    data = predictions,
    aes(age_weeks, ymin = adjusted_ci_low, ymax = adjusted_ci_high),
    fill = "#9ECAE1", alpha = 0.32
  ) +
  geom_line(
    data = predictions,
    aes(age_weeks, adjusted_fraction),
    colour = "#2F6B9A", linewidth = 0.58
  ) +
  geom_point(
    data = composition,
    aes(age_weeks, cell_fraction, fill = dataset_label),
    shape = 21, size = 1.55, stroke = 0.22, colour = "white", alpha = 0.88,
    position = position_jitter(width = 0.55, height = 0, seed = 20260714)
  ) +
  facet_wrap(
    vars(celltype_label), ncol = 3, scales = "free_y",
    axes = "all_x", axis.labels = "all_x"
  ) +
  scale_fill_manual(values = dataset_palette, drop = FALSE, guide = dataset_guide) +
  scale_x_continuous(
    breaks = seq(0, 90, 30), limits = c(4, 90),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 0.1),
    expand = expansion(mult = c(0.03, 0.08))
  ) +
  labs(
    x = "Chronological age (weeks)",
    y = "Fraction of confident cells",
    title = "Dataset-adjusted ovarian cell-class composition across adult age",
    subtitle = "Natural-ageing samples; panels use independent y-axis ranges",
    caption = wrap_text(paste0(
      "Each point is one biological sample (n = ", n_distinct(composition$sample_key), "). ",
      "Lines and 95% CIs are equal-dataset-weighted marginal predictions from ",
      "sample-level linear models of cell fraction on age per 10 weeks and dataset."
    ), width = 126)
  ) +
  theme_composition(base_size = 7) +
  theme(
    panel.spacing = grid::unit(3.0, "mm"),
    strip.text = element_text(size = 7.2, face = "bold", hjust = 0),
    legend.position = "bottom"
  )

overview_stem <- file.path(
  out_dir,
  "dataset_adjusted_celltype_composition_fraction_vs_age_overview"
)
save_plot_bundle(
  overview_plot, overview_stem,
  width_mm = 183, height_mm = 190, make_tiff = TRUE
)

cat("Dataset-adjusted composition figures written to: ", out_dir, "\n", sep = "")
cat("Theca / stroma adjusted slope: ", sprintf("%.3f", theca_slope),
    " percentage points per 10 weeks\n", sep = "")
