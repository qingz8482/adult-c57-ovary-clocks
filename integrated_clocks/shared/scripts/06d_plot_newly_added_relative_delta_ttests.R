#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- sub("^--file=", "", script_arg[[1]])
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
projection_dir <- file.path(project_dir, "results", "external_projection")
input_path <- file.path(
  projection_dir,
  "simplified_classes",
  "external_projection_predictions_all_simplified.tsv"
)
official_relative_path <- file.path(
  projection_dir,
  "simplified_classes",
  "external_projection_relative_to_controls_all_simplified.tsv"
)
out_dir <- file.path(projection_dir, "diagnostics", "newly added with pvalue")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

pred <- read_tsv(input_path)
official_relative <- read_tsv(official_relative_path)

required_columns <- c(
  "sample_key", "dataset_folder", "model_family", "delta_age_week",
  "projection_class", "projection_group", "comparison_id"
)
missing_columns <- setdiff(required_columns, names(pred))
if (length(missing_columns)) {
  stop("Prediction input is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

contrast_spec <- data.table(
  dataset_folder = c(
    "GSE106339",
    "GSE128240",
    rep("GSE154890", 3),
    rep("GSE235021", 4),
    rep("GSE249151", 2),
    "GSE98497",
    "PRJNA787028"
  ),
  comparison_id = c(
    "GSE106339_prepubertal_GOAT_KO_vs_WT",
    "GSE128240_CTX_vs_control",
    rep("GSE154890_natural_age_external_placement", 3),
    rep("GSE235021_rice_diets_vs_polished_6m", 2),
    rep("GSE235021_rice_diets_vs_polished_9m", 2),
    rep("GSE249151_VCD_vs_vehicle", 2),
    "GSE98497_Mtor_cKO_vs_WT",
    "PRJNA787028_high_vs_low_housing_density"
  ),
  contrast_id = c(
    "GSE106339_GOAT_KO_vs_WT_3wk",
    "GSE128240_CTX_vs_control",
    "GSE154890_6m_vs_3m",
    "GSE154890_9m_vs_3m",
    "GSE154890_12m_vs_3m",
    "GSE235021_black_vs_polished_6m",
    "GSE235021_brown_vs_polished_6m",
    "GSE235021_black_vs_polished_9m",
    "GSE235021_brown_vs_polished_9m",
    "GSE249151_VCD_day15_vs_vehicle",
    "GSE249151_VCD_day30_vs_vehicle",
    "GSE98497_Mtor_cKO_vs_WT",
    "PRJNA787028_8pc_vs_2pc"
  ),
  control_group = c(
    "3-week WT",
    "same-study control",
    rep("13-week natural-age reference", 3),
    "6-month polished rice",
    "6-month polished rice",
    "9-month polished rice",
    "9-month polished rice",
    "sesame-oil vehicle control",
    "sesame-oil vehicle control",
    "WT control",
    "2 mice per cage"
  ),
  experimental_group = c(
    "3-week GOAT KO",
    "cyclophosphamide",
    "26.1-week natural-age reference",
    "39.1-week natural-age reference",
    "52.2-week natural-age reference",
    "black rice",
    "white brown rice",
    "black rice",
    "white brown rice",
    "VCD day 15",
    "VCD day 30",
    "Mtor oocyte-conditional KO",
    "8 mice per cage"
  ),
  control_label = c(
    "WT",
    "Control",
    rep("3 months", 3),
    "6m polished",
    "6m polished",
    "9m polished",
    "9m polished",
    "Vehicle",
    "Vehicle",
    "WT",
    "2/cage"
  ),
  experimental_label = c(
    "GOAT KO",
    "Cyclophosphamide",
    "6 months",
    "9 months",
    "12 months",
    "6m black",
    "6m brown",
    "9m black",
    "9m brown",
    "VCD day 15",
    "VCD day 30",
    "Mtor cKO",
    "8/cage"
  ),
  interpretation_flag = c(
    "developmental_OOD_descriptive",
    "matched_same-study_control",
    rep("natural_age_reference_placement", 3),
    rep("matched_same-study_age-stratified_control", 4),
    "matched_same-study_control",
    "control_age_confounded",
    "matched_same-study_control",
    "matched_same-study_control"
  )
)
contrast_spec[, contrast_order := seq_len(.N), by = dataset_folder]

new_datasets <- unique(contrast_spec$dataset_folder)
observed_datasets <- unique(pred[dataset_folder %in% new_datasets]$dataset_folder)
if (!setequal(new_datasets, observed_datasets)) {
  stop(
    "The seven-dataset input is incomplete: ",
    paste(setdiff(new_datasets, observed_datasets), collapse = ", "),
    call. = FALSE
  )
}

relative_lookup <- official_relative[, .(
  sample_key,
  model_family,
  official_relative_delta_age = as.numeric(relative_delta_age)
)]
relative_lookup <- unique(relative_lookup, by = c("sample_key", "model_family"))
pred <- merge(
  pred,
  relative_lookup,
  by = c("sample_key", "model_family"),
  all.x = TRUE,
  sort = FALSE
)

make_contrast_data <- function(spec_row) {
  dt <- pred[
    dataset_folder == spec_row$dataset_folder &
      comparison_id == spec_row$comparison_id &
      projection_group %chin% c(spec_row$control_group, spec_row$experimental_group)
  ]
  if (!nrow(dt)) {
    stop("No rows found for contrast ", spec_row$contrast_id, call. = FALSE)
  }
  dt[, `:=`(
    contrast_id = spec_row$contrast_id,
    contrast_order = spec_row$contrast_order,
    comparison_role = fifelse(
      projection_group == spec_row$control_group,
      "Control",
      "Experimental"
    ),
    plot_group = fifelse(
      projection_group == spec_row$control_group,
      spec_row$control_label,
      spec_row$experimental_label
    ),
    control_group = spec_row$control_group,
    experimental_group = spec_row$experimental_group,
    interpretation_flag = spec_row$interpretation_flag
  )]
  dt[, control_mean_delta_age_week := mean(
    delta_age_week[comparison_role == "Control"],
    na.rm = TRUE
  ), by = model_family]
  dt[, derived_relative_delta_age_week :=
       delta_age_week - control_mean_delta_age_week]
  dt[, relative_delta_age_for_plot := fifelse(
    is.finite(official_relative_delta_age),
    official_relative_delta_age,
    derived_relative_delta_age_week
  )]
  dt[, relative_delta_source := fifelse(
    is.finite(official_relative_delta_age),
    "official_matched_control_relative_delta",
    fifelse(
      spec_row$dataset_folder == "GSE106339",
      "descriptive_WT_centered_delta_prepubertal_OOD",
      "descriptive_3month_reference_centered_delta"
    )
  )]
  dt
}

contrast_data <- rbindlist(
  lapply(seq_len(nrow(contrast_spec)), function(i) make_contrast_data(contrast_spec[i])),
  use.names = TRUE,
  fill = TRUE
)
contrast_data[, display_projection_class := projection_class]
contrast_data[dataset_folder == "GSE106339", display_projection_class := fifelse(
  comparison_role == "Control",
  "reference_control",
  "genetic_model"
)]

duplicated_keys <- contrast_data[, .N, by = .(
  contrast_id, sample_key, model_family
)][N != 1L]
if (nrow(duplicated_keys)) {
  stop("Duplicated sample/model rows were created in contrast assembly.", call. = FALSE)
}
if (any(!is.finite(contrast_data$relative_delta_age_for_plot))) {
  stop("Non-finite plot-relative delta ages remain.", call. = FALSE)
}

welch_test <- function(dt) {
  control <- dt[comparison_role == "Control"]$relative_delta_age_for_plot
  experimental <- dt[comparison_role == "Experimental"]$relative_delta_age_for_plot
  if (length(control) < 2L || length(experimental) < 2L) {
    stop("Welch test requires at least two samples per group.", call. = FALSE)
  }
  result <- tryCatch(
    t.test(experimental, control, alternative = "two.sided", var.equal = FALSE),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    stop("Welch test failed: ", conditionMessage(result), call. = FALSE)
  }
  data.table(
    n_control = length(control),
    n_experimental = length(experimental),
    mean_control_relative_delta_age_week = mean(control),
    mean_experimental_relative_delta_age_week = mean(experimental),
    mean_difference_week = mean(experimental) - mean(control),
    t_statistic = unname(result$statistic),
    degrees_of_freedom = unname(result$parameter),
    p_value = result$p.value,
    confidence_interval_low_week = result$conf.int[[1]],
    confidence_interval_high_week = result$conf.int[[2]]
  )
}

stats <- contrast_data[, welch_test(.SD), by = .(
  dataset_folder,
  contrast_id,
  contrast_order,
  model_family,
  control_group,
  experimental_group,
  interpretation_flag,
  relative_delta_source
)]
stats[, p_value_bh_all_tests := p.adjust(p_value, method = "BH")]
stats[, significance_label := fifelse(
  p_value < 0.0001, "****",
  fifelse(
    p_value < 0.001, "***",
    fifelse(
      p_value < 0.01, "**",
      fifelse(
        p_value < 0.05, "*",
        paste0("p = ", formatC(p_value, format = "f", digits = 3))
      )
    )
  )
)]
stats[, test_method := "two-sided Welch t-test"]
stats[, significance_basis := "unadjusted two-sided p-value"]
setorder(stats, dataset_folder, contrast_order, model_family)

stats_path <- file.path(out_dir, "newly_added7_relative_delta_welch_ttests.tsv")
source_path <- file.path(out_dir, "newly_added7_relative_delta_plot_source_data.tsv")
fwrite(stats, stats_path, sep = "\t")
fwrite(
  contrast_data[, .(
    dataset_folder, contrast_id, contrast_order, sample_key, sample_id,
    model_family, projection_class, display_projection_class, comparison_role, plot_group,
    delta_age_week, official_relative_delta_age,
    derived_relative_delta_age_week, relative_delta_age_for_plot,
    relative_delta_source, interpretation_flag
  )],
  source_path,
  sep = "\t"
)

class_palette <- c(
  developmental_OOD = "#B07AA1",
  diet = "#F28E2B",
  genetic_model = "#E15759",
  intervention = "#4E79A7",
  reference_control = "#7F7F7F",
  stress = "#9C755F"
)

theme_projection <- function() {
  theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      axis.text.x = element_text(angle = 32, hjust = 1),
      legend.position = "bottom",
      plot.caption = element_text(size = 8, color = "grey35", hjust = 0)
    )
}

save_plot_pair <- function(plot, stem, width = 10, height = 6.5) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = 300)
  c(pdf_path, png_path)
}

plot_paths <- character()
for (dataset_id in new_datasets) {
  dt <- copy(contrast_data[dataset_folder == dataset_id])
  st <- copy(stats[dataset_folder == dataset_id])
  spec <- contrast_spec[dataset_folder == dataset_id]
  plot_palette <- class_palette
  if (dataset_id == "GSE106339") {
    plot_palette[c("reference_control", "genetic_model")] <- c("#00BFC4", "#F8766D")
  }

  group_levels <- unique(as.vector(rbind(
    spec$control_label,
    spec$experimental_label
  )))
  dt[, plot_group := factor(plot_group, levels = group_levels)]
  group_position <- setNames(seq_along(group_levels), group_levels)
  st <- merge(
    st,
    unique(spec[, .(
      contrast_id, control_label, experimental_label
    )]),
    by = "contrast_id",
    all.x = TRUE
  )
  st[, `:=`(
    x_start = as.numeric(group_position[control_label]),
    x_end = as.numeric(group_position[experimental_label])
  )]

  panel_range <- dt[, .(
    panel_min = min(relative_delta_age_for_plot),
    panel_max = max(relative_delta_age_for_plot)
  ), by = model_family]
  panel_range[, panel_span := pmax(panel_max - panel_min, 1)]
  st <- merge(st, panel_range, by = "model_family", all.x = TRUE)
  st[, y_position := panel_max + panel_span * (0.12 + 0.15 * (contrast_order - 1L))]
  st[, tick_height := panel_span * 0.035]

  dataset_caption <- paste(
    "Brackets show two-sided Welch t-tests. Stars denote unadjusted p < 0.05",
    "(* <0.05, ** <0.01, *** <0.001, **** <0.0001); otherwise the exact p-value is shown.",
    "BH-adjusted values are provided in the accompanying TSV."
  )
  if (dataset_id == "GSE106339") {
    dataset_caption <- paste(
      dataset_caption,
      "GSE106339 is a descriptive comparison centered on its same-study WT reference-control group."
    )
  }
  if (dataset_id == "GSE154890") {
    dataset_caption <- paste(
      dataset_caption,
      "GSE154890 is natural-age reference placement centered on the 3-month group, not a treatment effect."
    )
  }
  if (dataset_id == "GSE249151") {
    dataset_caption <- paste(
      dataset_caption,
      "The VCD day-30 comparison is age-confounded relative to the vehicle collection time."
    )
  }
  dataset_caption <- paste(strwrap(dataset_caption, width = 175), collapse = "\n")

  p <- ggplot(
    dt,
    aes(plot_group, relative_delta_age_for_plot, color = display_projection_class)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12) +
    geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
    geom_segment(
      data = st,
      aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_segment(
      data = st,
      aes(
        x = x_start, xend = x_start,
        y = y_position, yend = y_position - tick_height
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_segment(
      data = st,
      aes(
        x = x_end, xend = x_end,
        y = y_position, yend = y_position - tick_height
      ),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_text(
      data = st,
      aes(
        x = (x_start + x_end) / 2,
        y = y_position + tick_height,
        label = significance_label
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 3.2,
      vjust = 0
    ) +
    facet_wrap(~model_family, scales = "free_y") +
    scale_color_manual(values = plot_palette, drop = FALSE) +
    expand_limits(y = max(st$y_position + st$tick_height * 4)) +
    labs(
      title = paste0("Matched-control relative delta age: ", dataset_id),
      subtitle = "Significant comparisons use stars; non-significant comparisons show exact p-values",
      x = NULL,
      y = "Relative delta age (weeks)",
      color = "Projection class",
      caption = dataset_caption
    ) +
    theme_projection()

  plot_paths <- c(
    plot_paths,
    save_plot_pair(
      p,
      paste0("05_relative_delta_age_matched_controls_", dataset_id, "_pvalue")
    )
  )
}

expected_plot_count <- length(new_datasets) * 2L
if (length(plot_paths) != expected_plot_count || any(!file.exists(plot_paths))) {
  stop("Not all seven PDF/PNG plot pairs were written.", call. = FALSE)
}

readme <- c(
  "# Relative delta-age Welch tests for seven newly added datasets",
  "",
  "- Each model family was tested separately with a two-sided Welch t-test.",
  "- Stars use unadjusted two-sided p-values; non-significant comparisons show exact p-values.",
  "- BH-adjusted p-values across all tests are retained in the statistics table.",
  "- Existing official matched-control relative delta ages were used when available.",
  "- In GSE106339, GOAT-KO samples are genetic models and WT samples are reference controls; the WT-centered contrast remains descriptive because all samples are prepubertal/OOD.",
  "- GSE154890 is natural-age placement centered on the 3-month reference, not a treatment-effect comparison.",
  "- GSE249151 VCD day 30 versus vehicle is flagged as control-age-confounded.",
  "",
  paste0("- Statistics: `", basename(stats_path), "`"),
  paste0("- Plot source data: `", basename(source_path), "`")
)
writeLines(readme, file.path(out_dir, "README.md"))

message(
  "Wrote ", length(new_datasets), " PDF/PNG plot pairs and ",
  nrow(stats), " model-specific Welch tests to: ", out_dir
)
