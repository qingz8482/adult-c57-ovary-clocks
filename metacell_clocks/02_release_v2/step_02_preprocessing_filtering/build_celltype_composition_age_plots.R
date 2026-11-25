#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  supplied_path <- sub("^--file=", "", script_arg)
  supplied_path <- gsub("~\\+~", " ", supplied_path)
  normalizePath(supplied_path, mustWork = TRUE)
} else {
  normalizePath("build_celltype_composition_age_plots.R", mustWork = TRUE)
}
step_dir <- dirname(script_path)
out_dir <- file.path(step_dir, "celltype_composition_age")
individual_dir <- file.path(out_dir, "individual_celltypes")
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)

audit_path <- file.path(step_dir, "confident_celltype_filter_audit.tsv")
metadata_path <- file.path(step_dir, "filtered_metacell_metadata.tsv")
if (!file.exists(audit_path) || !file.exists(metadata_path)) {
  stop("Missing preprocessing inputs in ", step_dir, call. = FALSE)
}

audit <- fread(audit_path) |> as_tibble()
metacell_meta <- fread(metadata_path) |> as_tibble()

# The carried-forward cell audit uses the legacy germline label; harmonize it
# to the canonical class name used by the active metacell metadata.
audit <- audit |>
  mutate(harmonized_cell_class = recode(harmonized_cell_class, Oocyte_Germ = "Oocyte_Germline"))

required_audit <- c(
  "dataset_id", "sample_key", "harmonized_cell_class", "n_cells",
  "confident_celltype_filter_pass"
)
required_meta <- c(
  "dataset_id", "sample_key", "sample_id", "age_weeks", "reproductive_age_bin",
  "intervention_class", "cv_role"
)
if (length(setdiff(required_audit, names(audit))) > 0L) {
  stop("Audit table is missing: ", paste(setdiff(required_audit, names(audit)), collapse = ", "), call. = FALSE)
}
if (length(setdiff(required_meta, names(metacell_meta))) > 0L) {
  stop("Metacell metadata is missing: ", paste(setdiff(required_meta, names(metacell_meta)), collapse = ", "), call. = FALSE)
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

sample_meta <- metacell_meta |>
  filter(cv_role == "metacell_cv_candidate", intervention_class == "natural_aging") |>
  transmute(
    dataset_id,
    sample_key,
    sample_id,
    age_weeks = as.numeric(age_weeks),
    reproductive_age_bin = as.character(reproductive_age_bin)
  ) |>
  distinct()

if (anyDuplicated(sample_meta$sample_key)) {
  stop("Natural-aging sample metadata is not unique by sample_key", call. = FALSE)
}
if (anyNA(sample_meta$age_weeks)) {
  stop("Natural-aging sample metadata contains missing age_weeks", call. = FALSE)
}

unexpected_celltypes <- audit |>
  filter(confident_celltype_filter_pass %in% TRUE) |>
  distinct(harmonized_cell_class) |>
  pull(harmonized_cell_class) |>
  setdiff(celltype_levels)
if (length(unexpected_celltypes) > 0L) {
  stop("Unexpected confident cell classes: ", paste(unexpected_celltypes, collapse = ", "), call. = FALSE)
}

observed_counts <- audit |>
  filter(confident_celltype_filter_pass %in% TRUE) |>
  inner_join(sample_meta |> select(dataset_id, sample_key), by = c("dataset_id", "sample_key")) |>
  group_by(dataset_id, sample_key, harmonized_cell_class) |>
  summarise(n_cells = sum(as.numeric(n_cells)), .groups = "drop")

composition <- sample_meta |>
  crossing(harmonized_cell_class = celltype_levels) |>
  left_join(observed_counts, by = c("dataset_id", "sample_key", "harmonized_cell_class")) |>
  mutate(n_cells = replace_na(n_cells, 0)) |>
  group_by(sample_key) |>
  mutate(
    total_confident_cells = sum(n_cells),
    cell_fraction = n_cells / total_confident_cells
  ) |>
  ungroup() |>
  mutate(
    celltype_label = factor(
      unname(celltype_labels[harmonized_cell_class]),
      levels = unname(celltype_labels[celltype_levels])
    ),
    dataset_label = factor(
      unname(dataset_labels[dataset_id]),
      levels = unname(dataset_labels[names(dataset_labels)])
    )
  ) |>
  arrange(match(harmonized_cell_class, celltype_levels), age_weeks, dataset_id, sample_id)

if (any(composition$total_confident_cells <= 0)) {
  stop("At least one sample has no confident cells", call. = FALSE)
}
fraction_sums <- composition |>
  summarise(fraction_sum = sum(cell_fraction), .by = sample_key)
if (any(abs(fraction_sums$fraction_sum - 1) > 1e-10)) {
  stop("Cell fractions do not sum to one within sample", call. = FALSE)
}
if (anyNA(composition$dataset_label)) {
  stop("At least one dataset lacks a display label", call. = FALSE)
}

fit_one_celltype <- function(df) {
  spearman <- suppressWarnings(cor.test(df$age_weeks, df$cell_fraction, method = "spearman", exact = FALSE))
  model_df <- df |>
    mutate(age_10_weeks = age_weeks / 10, dataset_id = factor(dataset_id))
  adjusted_model <- lm(cell_fraction ~ age_10_weeks + dataset_id, data = model_df)
  model_coef <- summary(adjusted_model)$coefficients
  age_row <- model_coef["age_10_weeks", ]
  age_ci <- confint(adjusted_model, "age_10_weeks", level = 0.95)

  tibble(
    n_samples = nrow(df),
    n_datasets = n_distinct(df$dataset_id),
    age_min_weeks = min(df$age_weeks),
    age_max_weeks = max(df$age_weeks),
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    adjusted_fraction_change_per_10_weeks = unname(age_row["Estimate"]),
    adjusted_fraction_change_ci_low = unname(age_ci[1]),
    adjusted_fraction_change_ci_high = unname(age_ci[2]),
    adjusted_age_p = unname(age_row["Pr(>|t|)"]),
    adjusted_model_r_squared = summary(adjusted_model)$r.squared
  )
}

statistics <- composition |>
  group_by(harmonized_cell_class, celltype_label) |>
  group_modify(~fit_one_celltype(.x)) |>
  ungroup() |>
  mutate(
    spearman_q_bh = p.adjust(spearman_p, method = "BH"),
    adjusted_age_q_bh = p.adjust(adjusted_age_p, method = "BH"),
    adjusted_percentage_point_change_per_10_weeks = 100 * adjusted_fraction_change_per_10_weeks,
    adjusted_percentage_point_ci_low = 100 * adjusted_fraction_change_ci_low,
    adjusted_percentage_point_ci_high = 100 * adjusted_fraction_change_ci_high
  ) |>
  arrange(adjusted_age_q_bh, match(harmonized_cell_class, celltype_levels))

fwrite(
  composition |>
    mutate(celltype_label = as.character(celltype_label), dataset_label = as.character(dataset_label)),
  file.path(out_dir, "sample_celltype_composition_by_age.tsv"),
  sep = "\t"
)
fwrite(
  statistics |> mutate(celltype_label = as.character(celltype_label)),
  file.path(out_dir, "celltype_age_association_statistics.tsv"),
  sep = "\t"
)

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

save_plot_bundle <- function(plot, stem, width_mm, height_mm, tiff = FALSE) {
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

  if (tiff) {
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

make_single_plot <- function(celltype_id) {
  df <- composition |> filter(harmonized_cell_class == celltype_id)
  stat <- statistics |> filter(harmonized_cell_class == celltype_id)
  upper <- max(df$cell_fraction) * 1.12
  upper <- max(upper, 0.015)
  subtitle <- wrap_text(sprintf(
    "Pooled Spearman ρ = %.2f (BH q = %s); dataset-adjusted Δ = %+.2f pp per 10 weeks (BH q = %s)",
    stat$spearman_rho,
    pvalue(stat$spearman_q_bh, accuracy = 0.001, add_p = FALSE),
    stat$adjusted_percentage_point_change_per_10_weeks,
    pvalue(stat$adjusted_age_q_bh, accuracy = 0.001, add_p = FALSE)
  ), width = 82)

  ggplot(df, aes(age_weeks, cell_fraction)) +
    geom_smooth(
      method = "glm", formula = y ~ x,
      method.args = list(family = quasibinomial(link = "logit")), se = TRUE,
      colour = "#2F6B9A", fill = "#9ECAE1", linewidth = 0.65, alpha = 0.28
    ) +
    geom_point(
      aes(fill = dataset_label), shape = 21, size = 2.25,
      stroke = 0.3, colour = "white", alpha = 0.9,
      position = position_jitter(width = 0.55, height = 0, seed = 20260714)
    ) +
    scale_fill_manual(values = dataset_palette, drop = FALSE, guide = dataset_guide) +
    scale_x_continuous(
      breaks = seq(0, 90, 15), limits = c(4, 90),
      expand = expansion(mult = c(0, 0.01))
    ) +
    scale_y_continuous(
      labels = label_percent(accuracy = 0.1),
      expand = expansion(mult = c(0, 0.03))
    ) +
    coord_cartesian(ylim = c(0, upper)) +
    labs(
      x = "Chronological age (weeks)",
      y = "Fraction of confident cells",
      title = unname(celltype_labels[celltype_id]),
      subtitle = subtitle,
      caption = wrap_text(paste0(
        "Each point is one biological sample (n = ", nrow(df), "). ",
        "Line and 95% CI: pooled sample-level quasi-binomial logit trend; pp: percentage points."
      ), width = 105)
    ) +
    theme_composition()
}

single_plots <- setNames(lapply(celltype_levels, make_single_plot), celltype_levels)
for (celltype_id in celltype_levels) {
  stem <- file.path(
    individual_dir,
    paste0("cell_fraction_vs_age_", tolower(gsub("_", "-", celltype_id)))
  )
  save_plot_bundle(single_plots[[celltype_id]], stem, width_mm = 120, height_mm = 92)
}

combined_plot <- ggplot(composition, aes(age_weeks, cell_fraction)) +
  geom_smooth(
    method = "glm", formula = y ~ x,
    method.args = list(family = quasibinomial(link = "logit")), se = TRUE,
    colour = "#2F6B9A", fill = "#9ECAE1", linewidth = 0.55, alpha = 0.28
  ) +
  geom_point(
    aes(fill = dataset_label), shape = 21, size = 1.55,
    stroke = 0.22, colour = "white", alpha = 0.88,
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
  scale_y_continuous(labels = label_percent(accuracy = 0.1), expand = expansion(mult = c(0.02, 0.08))) +
  labs(
    x = "Chronological age (weeks)",
    y = "Fraction of confident cells",
    title = "Ovarian cell-class composition across adult age",
    subtitle = "Natural-aging samples; panels use independent y-axis ranges",
    caption = paste0(
      "Each point is one biological sample (n = ", n_distinct(composition$sample_key), "). ",
      "Lines and 95% CIs are pooled descriptive quasi-binomial logit trends; dataset-adjusted statistics are provided separately."
    )
  ) +
  theme_composition(base_size = 7) +
  theme(
    panel.spacing = grid::unit(3.0, "mm"),
    strip.text = element_text(size = 7.2, face = "bold", hjust = 0),
    legend.position = "bottom"
  )

save_plot_bundle(
  combined_plot,
  file.path(out_dir, "celltype_composition_fraction_vs_age_overview"),
  width_mm = 183,
  height_mm = 190,
  tiff = TRUE
)

hard_checks <- tibble(
  check = c(
    "one_row_per_sample_celltype", "all_expected_celltypes", "no_missing_age",
    "positive_sample_cell_totals", "fractions_sum_to_one", "all_individual_pngs_exist",
    "overview_export_bundle_exists"
  ),
  pass = c(
    nrow(composition) == nrow(sample_meta) * length(celltype_levels),
    setequal(unique(composition$harmonized_cell_class), celltype_levels),
    !anyNA(composition$age_weeks),
    all(composition$total_confident_cells > 0),
    all(abs(fraction_sums$fraction_sum - 1) <= 1e-10),
    all(file.exists(file.path(
      individual_dir,
      paste0("cell_fraction_vs_age_", tolower(gsub("_", "-", celltype_levels)), ".png")
    ))),
    all(file.exists(file.path(
      out_dir,
      paste0("celltype_composition_fraction_vs_age_overview", c(".svg", ".pdf", ".png", ".tiff"))
    )))
  ),
  value = c(
    paste0(nrow(composition), " rows for ", nrow(sample_meta), " samples x ", length(celltype_levels), " cell types"),
    paste(sort(unique(composition$harmonized_cell_class)), collapse = ";"),
    as.character(sum(is.na(composition$age_weeks))),
    as.character(min(composition$total_confident_cells)),
    format(max(abs(fraction_sums$fraction_sum - 1)), scientific = TRUE),
    as.character(length(celltype_levels)),
    "SVG;PDF;PNG;TIFF"
  )
)
fwrite(hard_checks, file.path(out_dir, "celltype_composition_age_hard_checks.tsv"), sep = "\t")
if (!all(hard_checks$pass)) {
  stop("One or more cell-composition hard checks failed", call. = FALSE)
}

writeLines(
  c(
    "# Cell-type composition fraction across age",
    "",
    "This analysis uses confident cells from natural-aging, metacell-CV-candidate samples. Fractions are computed within each biological sample after summing confirmed and probable cells; missing sample-cell-type combinations are represented as zero.",
    "",
    "Each cell type has a standalone sample-level plot in `individual_celltypes/`, and the overview presents all nine classes in a three-by-three grid. Points are colored by source dataset. Pooled quasi-binomial logit trends are descriptive; dataset-adjusted sample-level linear models and pooled Spearman correlations are reported in `celltype_age_association_statistics.tsv`, with Benjamini-Hochberg correction across the nine cell types.",
    "",
    "Because cell recovery depends on dissociation and capture, these figures show composition associations rather than direct evidence of in-vivo abundance changes. Dataset-adjusted results should be emphasized over pooled trends.",
    "",
    "Run from any directory with:",
    "",
    "```bash",
    "Rscript build_celltype_composition_age_plots.R",
    "```"
  ),
  file.path(out_dir, "README.md")
)

message("Wrote cell-type composition-versus-age outputs to: ", out_dir)
