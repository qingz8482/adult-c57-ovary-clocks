#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
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
figure_dir <- file.path(release_dir, "figures")
plot_data_dir <- file.path(release_dir, "plot_data")

paths <- list(
  oof = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_oof_predictions.tsv"),
  cv_eval = file.path(run_root, "step_05_internal_cv_validation", "celltype_clock_cv_evaluation_summary.tsv")
)
missing_paths <- paths[!file.exists(unlist(paths))]
if (length(missing_paths) > 0) stop("Missing input: ", paste(unlist(missing_paths), collapse = "; "), call. = FALSE)

age_bin_levels <- c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
age_bin_labels <- c(
  young_adult = "Young adult",
  mature_adult = "Mature adult",
  advanced_maternal_age = "Advanced maternal age",
  post_reproductive_ageing = "Post-reproductive"
)
age_bin_shapes <- c(
  young_adult = 16,
  mature_adult = 17,
  advanced_maternal_age = 15,
  post_reproductive_ageing = 18
)
dataset_palette <- c(
  "2022_eLife_Morris_cycling_murine_ovary_atlas" = "#0072B2",
  "2024_Cell_Winkler_FRT_scRNA_spatial" = "#D55E00",
  "2024_LifeMedicine_Jiang_whole_life_mouse_ovary" = "#009E73",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "#CC79A7",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary" = "#7F3C8D",
  "2025_CommBio_Zhao_developing_ovary_scRNA_spatial" = "#11A579",
  "2026_CommBio_Galligos_ovarian_inflammaging" = "#E69F00",
  "2025_Kim_PRJNA863443_mouse_menopause_models" = "#56B4E9"
)

theme_clock_plot <- function(base_size = 8.5) {
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

normalize_dataset_id <- function(dataset_id, dataset_folder) {
  dataset_id <- as.character(dataset_id)
  dataset_folder <- as.character(dataset_folder)
  if_else(is.na(dataset_id) | dataset_id == "", dataset_folder, dataset_id)
}

oof_tbl <- fread(paths$oof) %>% as_tibble()
cv_eval <- fread(paths$cv_eval) %>% as_tibble()

facet_metrics <- cv_eval %>%
  mutate(
    facet_label = harmonized_cell_class,
    metric_label = sprintf("n=%s\nMAE=%.1f wk\nslope=%.2f\nR2=%.2f", n_metacells, mae_weeks, slope, r2)
  ) %>%
  select(harmonized_cell_class, facet_label, metric_label)

plot_data <- oof_tbl %>%
  left_join(facet_metrics, by = "harmonized_cell_class") %>%
  mutate(
    dataset_id = normalize_dataset_id(dataset_id, dataset_folder),
    reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels),
    dataset_id = factor(dataset_id, levels = names(dataset_palette))
  )

fwrite(plot_data, file.path(plot_data_dir, "filtered_elasticnet_metacell_internal_grouped_cv_plot_data.tsv"), sep = "\t")

metric_label_data <- facet_metrics %>%
  mutate(age_weeks = -Inf, predicted_age_weeks = Inf)
facet_labels <- setNames(facet_metrics$facet_label, facet_metrics$harmonized_cell_class)

p_cv <- ggplot(plot_data, aes(age_weeks, predicted_age_weeks, color = dataset_id, shape = reproductive_age_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_point(alpha = 0.55, size = 1.05) +
  geom_text(
    data = metric_label_data,
    aes(age_weeks, predicted_age_weeks, label = metric_label),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.04,
    size = 2.35,
    lineheight = 0.92,
    color = "grey20"
  ) +
  facet_wrap(~ harmonized_cell_class, scales = "free", labeller = labeller(harmonized_cell_class = facet_labels)) +
  scale_color_manual(values = dataset_palette, drop = FALSE) +
  scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
  theme_clock_plot(base_size = 8.5) +
  labs(
    title = "Filtered ElasticNet metacell clocks: internal grouped CV",
    x = "Observed age (weeks)",
    y = "Predicted age (weeks)",
    color = "Dataset",
    shape = "Age bin"
  )

ggsave(file.path(figure_dir, "filtered_elasticnet_metacell_internal_grouped_cv.pdf"), p_cv, width = 11, height = 8)
ggsave(file.path(figure_dir, "filtered_elasticnet_metacell_internal_grouped_cv.png"), p_cv, width = 11, height = 8, dpi = 240)

message("Refreshed filtered ElasticNet metacell internal grouped CV figure in ", figure_dir)
