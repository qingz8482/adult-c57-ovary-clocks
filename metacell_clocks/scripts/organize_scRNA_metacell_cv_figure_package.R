#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(stringr)
})

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)][1])
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
data_root <- normalizePath(Sys.getenv("OVARY_CLOCK_DATA_ROOT", unset = file.path(repo_root, "data")), mustWork = FALSE)
output_root <- normalizePath(Sys.getenv("OVARY_CLOCK_OUTPUT_ROOT", unset = file.path(repo_root, "outputs")), mustWork = FALSE)
run_root <- file.path(output_root, "metacell_clocks")
release_dir <- file.path(run_root, "step_08_figures_tables_release")
figure_root <- file.path(release_dir, "figures")
plot_data_root <- file.path(release_dir, "plot_data")
table_root <- file.path(release_dir, "tables")

sections <- c(
  "01_overview_eligibility",
  "02_internal_cv_performance",
  "03_feature_selection",
  "04_gsea",
  "05_stress_test_projection",
  "99_combined_panels"
)
for (section in sections) {
  dir.create(file.path(figure_root, section), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(plot_data_root, section), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(table_root, section), recursive = TRUE, showWarnings = FALSE)
}

paths <- list(
  eligibility = file.path(run_root, "step_03_metacell_cv_split", "celltype_metacell_cv_eligibility.tsv"),
  oof = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_oof_predictions.tsv"),
  cv_eval = file.path(run_root, "step_05_internal_cv_validation", "celltype_clock_cv_evaluation_summary.tsv"),
  final_models = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_final_models.rds"),
  coefficients = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_coefficients_by_fold.tsv"),
  selected_features = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_selected_features_by_fold.tsv"),
  residual_age = file.path(run_root, "step_05_internal_cv_validation", "residual_bias_by_age_bin.tsv"),
  residual_dataset = file.path(run_root, "step_05_internal_cv_validation", "residual_bias_by_dataset.tsv"),
  projection_metacell = file.path(run_root, "step_06_projection_stress_test", "prjna863443_projection_metacell_predictions.tsv"),
  projection_sample = file.path(run_root, "step_06_projection_stress_test", "prjna863443_projection_sample_celltype_summary.tsv"),
  projection_delta_old = file.path(run_root, "step_06_projection_stress_test", "prjna863443_projection_delta_by_celltype.tsv"),
  vcd_summary_old = file.path(run_root, "step_06_projection_stress_test", "vcd_projection_summary_by_celltype.tsv"),
  foxl2_summary_old = file.path(run_root, "step_06_projection_stress_test", "foxl2_projection_summary_by_celltype.tsv"),
  projection_delta = file.path(run_root, "step_06_projection_stress_test", "prjna863443_projection_delta_by_contrast_celltype.tsv"),
  vcd_summary = file.path(run_root, "step_06_projection_stress_test", "vcd_projection_summary_by_contrast_celltype.tsv"),
  foxl2_summary = file.path(run_root, "step_06_projection_stress_test", "foxl2_projection_summary_by_contrast_celltype.tsv"),
  prjna863443_metadata = file.path(data_root, "scRNA_dataset", "PRJNA863443", "metadata", "PRJNA863443_sample_metadata.tsv"),
  gsea_rank = file.path(run_root, "step_09_gsea", "elasticnet_coefficient_gsea", "elasticnet_coefficient_rank_table.tsv"),
  gsea_results = file.path(run_root, "step_09_gsea", "elasticnet_coefficient_gsea", "elasticnet_coefficient_gsea_results.tsv"),
  gsea_display = file.path(run_root, "step_09_gsea", "elasticnet_coefficient_gsea", "elasticnet_coefficient_gsea_display_terms.tsv"),
  gsea_status = file.path(run_root, "step_09_gsea", "elasticnet_coefficient_gsea", "elasticnet_coefficient_gsea_status.tsv"),
  ora_results = file.path(run_root, "step_09_gsea", "selected_gene_ora", "selected_gene_ora_results.tsv"),
  ora_display = file.path(run_root, "step_09_gsea", "selected_gene_ora", "selected_gene_ora_display_terms.tsv"),
  gsea_method = file.path(run_root, "step_09_gsea", "GSEA_METHOD.md"),
  fold_summary = file.path(run_root, "step_03_metacell_cv_split", "fold_summary_by_celltype.tsv"),
  sample_overlap = file.path(run_root, "step_03_metacell_cv_split", "sample_overlap_expected_audit.tsv")
)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required table: ", path, call. = FALSE)
  fread(path) %>% as_tibble()
}

eligibility <- read_tsv(paths$eligibility)
oof <- read_tsv(paths$oof)
cv_eval <- read_tsv(paths$cv_eval)
coef_tbl <- read_tsv(paths$coefficients)
selected_features <- read_tsv(paths$selected_features)
residual_age <- read_tsv(paths$residual_age)
residual_dataset <- read_tsv(paths$residual_dataset)
projection_metacell <- read_tsv(paths$projection_metacell)
projection_sample <- read_tsv(paths$projection_sample)
projection_delta <- read_tsv(paths$projection_delta)
projection_delta_old <- read_tsv(paths$projection_delta_old)
vcd_summary <- read_tsv(paths$vcd_summary)
foxl2_summary <- read_tsv(paths$foxl2_summary)
vcd_summary_old <- read_tsv(paths$vcd_summary_old)
foxl2_summary_old <- read_tsv(paths$foxl2_summary_old)
prjna863443_metadata <- read_tsv(paths$prjna863443_metadata)

final_models <- readRDS(paths$final_models)
final_coef_tbl <- bind_rows(lapply(names(final_models), function(cell_class) {
  final_models[[cell_class]]$coefficient_table %>%
    as_tibble() %>%
    mutate(
      harmonized_cell_class = cell_class,
      model_source = "final_all_metacells_elasticnet"
    )
})) %>%
  transmute(
    harmonized_cell_class,
    model_source,
    gene,
    elasticnet_coefficient = as.numeric(coefficient),
    abs_elasticnet_coefficient = abs(as.numeric(coefficient)),
    nonzero_coefficient,
    coefficient_direction = case_when(
      elasticnet_coefficient > 0 ~ "positive",
      elasticnet_coefficient < 0 ~ "negative",
      TRUE ~ "zero"
    ),
    abs_coefficient_rank,
    mean_expression,
    expression_fraction,
    training_sd,
    rho,
    abs_rho,
    best_alpha,
    best_lambda,
    selected_feature_cap,
    feature_source
  )
final_coef_table_path <- file.path(table_root, "03_feature_selection", "celltype_metacell_elasticnet_final_coefficients.tsv")

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
direction_palette <- c(positive = "#FF7566", negative = "#23C7C7")
eligibility_palette <- c(`FALSE` = "#F8766D", `TRUE` = "#00BFC4")
status_palette <- c(
  strong_internal_clock = "#1B9E77",
  usable_internal_clock = "#4E79A7",
  weak_internal_clock = "#E69F00",
  diagnostic_only = "#8F8F8F"
)
interpretation_shapes <- c(
  confident_projection_context = 16,
  cautious_projection_context = 17,
  diagnostic_only_projection_context = 4
)
interpretation_alpha <- c(
  confident_projection_context = 0.95,
  cautious_projection_context = 0.72,
  diagnostic_only_projection_context = 0.48
)
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

clean_token <- function(x) {
  x <- as.character(x)
  x <- if_else(is.na(x) | x == "", NA_character_, x)
  x <- str_replace_all(x, "[^A-Za-z0-9+/-]+", "_")
  str_replace_all(x, "^_|_$", "")
}

rel <- function(path) {
  path <- as.character(path)
  prefix <- paste0(normalizePath(release_dir, mustWork = FALSE), "/")
  path_norm <- normalizePath(path, mustWork = FALSE)
  ifelse(startsWith(path_norm, prefix), substring(path_norm, nchar(prefix) + 1L), path)
}

normalize_dataset_id <- function(dataset_id, dataset_folder) {
  dataset_id <- as.character(dataset_id)
  dataset_folder <- as.character(dataset_folder)
  if_else(is.na(dataset_id) | dataset_id == "", dataset_folder, dataset_id)
}

manifest <- tibble(
  figure_id = character(),
  figure_title = character(),
  figure_path_png = character(),
  figure_path_pdf = character(),
  plot_data_path = character(),
  source_table_path = character(),
  section = character(),
  status = character()
)

add_manifest <- function(figure_id, title, section, png_path, pdf_path, plot_data_path, source_table_path, status = "complete") {
  manifest <<- bind_rows(
    manifest,
    tibble(
      figure_id = figure_id,
      figure_title = title,
      figure_path_png = rel(png_path),
      figure_path_pdf = rel(pdf_path),
      plot_data_path = rel(plot_data_path),
      source_table_path = paste(rel(source_table_path), collapse = ";"),
      section = section,
      status = status
    )
  )
}

save_plot_with_data <- function(plot, data, section, figure_id, title, source_table_path,
                                width = 8, height = 5, status = "complete") {
  png_path <- file.path(figure_root, section, paste0(figure_id, ".png"))
  pdf_path <- file.path(figure_root, section, paste0(figure_id, ".pdf"))
  data_path <- file.path(plot_data_root, section, paste0(figure_id, "_plot_data.tsv"))
  fwrite(data, data_path, sep = "\t")
  ggsave(pdf_path, plot, width = width, height = height)
  ggsave(png_path, plot, width = width, height = height, dpi = 240)
  add_manifest(figure_id, title, section, png_path, pdf_path, data_path, source_table_path, status)
}

copy_section_table <- function(path, section) {
  if (file.exists(path)) file.copy(path, file.path(table_root, section, basename(path)), overwrite = TRUE)
}

for (path in c(paths$eligibility, paths$fold_summary, paths$sample_overlap)) copy_section_table(path, "01_overview_eligibility")
for (path in c(paths$oof, paths$cv_eval, paths$residual_age, paths$residual_dataset)) copy_section_table(path, "02_internal_cv_performance")
for (path in c(paths$coefficients, paths$selected_features)) copy_section_table(path, "03_feature_selection")
fwrite(final_coef_tbl, final_coef_table_path, sep = "\t")
for (path in c(paths$gsea_rank, paths$gsea_results, paths$gsea_display, paths$gsea_status, paths$ora_results, paths$ora_display, paths$gsea_method)) copy_section_table(path, "04_gsea")
for (path in c(
  paths$projection_metacell,
  paths$projection_sample,
  paths$projection_delta,
  paths$vcd_summary,
  paths$foxl2_summary,
  paths$prjna863443_metadata,
  paths$projection_delta_old,
  paths$vcd_summary_old,
  paths$foxl2_summary_old
)) copy_section_table(path, "05_stress_test_projection")

# 01 overview and eligibility
elig_plot_data <- eligibility %>%
  mutate(
    clock_eligible = as.character(eligible_for_metacell_cv),
    harmonized_cell_class = reorder(harmonized_cell_class, n_metacells)
  ) %>%
  arrange(desc(n_metacells))
p_elig <- ggplot(elig_plot_data, aes(n_metacells, harmonized_cell_class, fill = clock_eligible)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = eligibility_palette, breaks = names(eligibility_palette), drop = FALSE) +
  theme_clock() +
  labs(x = "Metacells", y = "Cell type", fill = "Clock eligible", title = "Metacell availability and model eligibility")
save_plot_with_data(
  p_elig, elig_plot_data, "01_overview_eligibility", "metacell_availability_and_model_eligibility",
  "Metacell availability and model eligibility", paths$eligibility, width = 7.2, height = 4.8
)

nonzero_by_celltype <- coef_tbl %>%
  filter(nonzero_coefficient) %>%
  group_by(harmonized_cell_class) %>%
  summarise(n_nonzero_genes = n_distinct(gene), .groups = "drop") %>%
  right_join(eligibility %>% select(harmonized_cell_class), by = "harmonized_cell_class") %>%
  mutate(n_nonzero_genes = if_else(is.na(n_nonzero_genes), 0L, as.integer(n_nonzero_genes))) %>%
  arrange(desc(n_nonzero_genes), harmonized_cell_class) %>%
  mutate(harmonized_cell_class = reorder(harmonized_cell_class, n_nonzero_genes))
p_nonzero <- ggplot(nonzero_by_celltype, aes(n_nonzero_genes, harmonized_cell_class)) +
  geom_col(fill = eligibility_palette[["TRUE"]], width = 0.72) +
  theme_clock() +
  labs(x = "Non-zero ElasticNet genes", y = "Cell type", title = "Non-zero selected genes by cell type")
save_plot_with_data(
  p_nonzero, nonzero_by_celltype, "01_overview_eligibility", "nonzero_selected_genes_by_celltype",
  "Non-zero selected genes by cell type", paths$coefficients, width = 7.2, height = 4.8
)

# 02 internal OOF performance
facet_metrics <- cv_eval %>%
  mutate(
    facet_label = sprintf(
      "%s\nn=%s, MAE=%.1f wk\nslope=%.2f, R2=%.2f",
      harmonized_cell_class, n_metacells, mae_weeks, slope, r2
    )
  ) %>%
  select(harmonized_cell_class, facet_label)
facet_labels <- setNames(facet_metrics$facet_label, facet_metrics$harmonized_cell_class)
oof_plot_data <- oof %>%
  left_join(facet_metrics %>% mutate(facet_label_tsv = gsub("\n", "; ", facet_label)) %>% select(harmonized_cell_class, facet_label_tsv),
            by = "harmonized_cell_class") %>%
  mutate(
    dataset_id = normalize_dataset_id(dataset_id, dataset_folder),
    reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels),
    dataset_id = factor(dataset_id, levels = names(dataset_palette))
  )
p_oof <- ggplot(oof_plot_data, aes(age_weeks, predicted_age_weeks, color = dataset_id, shape = reproductive_age_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_point(alpha = 0.55, size = 1.05) +
  facet_wrap(~ harmonized_cell_class, scales = "free", labeller = labeller(harmonized_cell_class = facet_labels)) +
  scale_color_manual(values = dataset_palette, drop = FALSE) +
  scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
  theme_clock() +
  theme(aspect.ratio = 1) +
  labs(
    x = "Observed age (weeks)",
    y = "Predicted age (weeks)",
    color = "Dataset",
    shape = "Age bin",
    title = "Filtered ElasticNet metacell clocks: internal out-of-fold predictions"
  )
save_plot_with_data(
  p_oof, oof_plot_data, "02_internal_cv_performance", "internal_oof_metacell_predictions_facet",
  "Filtered ElasticNet metacell clocks: internal out-of-fold predictions",
  paths$oof, width = 11, height = 8
)

for (ct in sort(unique(oof$harmonized_cell_class))) {
  ct_id <- safe_id(ct)
  ct_data <- oof_plot_data %>% filter(harmonized_cell_class == ct)
  ct_metric <- cv_eval %>% filter(harmonized_cell_class == ct)
  title <- sprintf(
    "%s: internal OOF metacell predictions (n=%s, MAE=%.1f, R2=%.2f)",
    ct, ct_metric$n_metacells[[1]], ct_metric$mae_weeks[[1]], ct_metric$r2[[1]]
  )
  title_wrapped <- str_wrap(title, width = 34)
  p_ct <- ggplot(ct_data, aes(age_weeks, predicted_age_weeks, color = dataset_id, shape = reproductive_age_bin)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_point(alpha = 0.65, size = 1.25) +
    scale_color_manual(values = dataset_palette, drop = FALSE) +
    scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
    theme_clock(base_size = 9) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 10, face = "bold", lineheight = 0.95)
    ) +
    labs(x = "Observed age (weeks)", y = "Predicted age (weeks)", color = "Dataset", shape = "Age bin", title = title_wrapped)
  save_plot_with_data(
    p_ct, ct_data, "02_internal_cv_performance", paste0("internal_oof_", ct_id, "_pred_vs_obs"),
    title, paths$oof, width = 6.4, height = 5.0
  )
}

perf_data <- cv_eval %>% arrange(desc(mae_weeks)) %>% mutate(harmonized_cell_class = reorder(harmonized_cell_class, mae_weeks))
p_mae <- ggplot(perf_data, aes(mae_weeks, harmonized_cell_class)) +
  geom_col(fill = "#2f6f73", width = 0.72) +
  theme_clock() +
  labs(x = "OOF metacell MAE (weeks)", y = "Cell type", title = "Cell-type clock MAE summary")
save_plot_with_data(p_mae, perf_data, "02_internal_cv_performance", "celltype_clock_performance_summary",
                    "Cell-type clock MAE summary", paths$cv_eval, width = 7.2, height = 4.8)

r2_data <- cv_eval %>% arrange(r2) %>% mutate(harmonized_cell_class = reorder(harmonized_cell_class, r2))
p_r2 <- ggplot(r2_data, aes(r2, harmonized_cell_class)) +
  geom_col(fill = "#4E79A7", width = 0.72) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_clock() +
  labs(x = "R2", y = "Cell type", title = "Cell-type clock R2 summary")
save_plot_with_data(p_r2, r2_data, "02_internal_cv_performance", "celltype_clock_r2_summary",
                    "Cell-type clock R2 summary", paths$cv_eval, width = 7.2, height = 4.8)

slope_data <- cv_eval %>% arrange(slope) %>% mutate(harmonized_cell_class = reorder(harmonized_cell_class, slope))
p_slope <- ggplot(slope_data, aes(slope, harmonized_cell_class)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_col(fill = "#59A14F", width = 0.72) +
  theme_clock() +
  labs(x = "Predicted ~ observed slope", y = "Cell type", title = "Cell-type clock slope summary")
save_plot_with_data(p_slope, slope_data, "02_internal_cv_performance", "celltype_clock_slope_summary",
                    "Cell-type clock slope summary", paths$cv_eval, width = 7.2, height = 4.8)

age_bias_data <- oof %>%
  group_by(harmonized_cell_class, reproductive_age_bin) %>%
  summarise(
    n_metacells = n(),
    mean_delta_age_weeks = mean(delta_age_weeks, na.rm = TRUE),
    median_delta_age_weeks = median(delta_age_weeks, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels))
p_age_bias <- ggplot(age_bias_data, aes(reproductive_age_bin, harmonized_cell_class, fill = mean_delta_age_weeks)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  scale_x_discrete(labels = age_bin_labels, drop = FALSE) +
  theme_clock() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Age bin", y = "Cell type", fill = "Mean delta\n(weeks)", title = "Residual bias by age bin and cell type")
save_plot_with_data(p_age_bias, age_bias_data, "02_internal_cv_performance", "residual_bias_by_age_bin_and_celltype",
                    "Residual bias by age bin and cell type", paths$oof, width = 8.2, height = 5.2)

dataset_bias_data <- oof %>%
  mutate(dataset_id = normalize_dataset_id(dataset_id, dataset_folder)) %>%
  group_by(harmonized_cell_class, dataset_id) %>%
  summarise(
    n_metacells = n(),
    mean_delta_age_weeks = mean(delta_age_weeks, na.rm = TRUE),
    median_delta_age_weeks = median(delta_age_weeks, na.rm = TRUE),
    .groups = "drop"
  )
p_dataset_bias <- ggplot(dataset_bias_data, aes(dataset_id, harmonized_cell_class, fill = mean_delta_age_weeks)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  theme_clock() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(x = "Dataset", y = "Cell type", fill = "Mean delta\n(weeks)", title = "Residual bias by dataset and cell type")
save_plot_with_data(p_dataset_bias, dataset_bias_data, "02_internal_cv_performance", "residual_bias_by_dataset_and_celltype",
                    "Residual bias by dataset and cell type", paths$oof, width = 9.5, height = 5.2)

# 03 feature selection
direction_data <- coef_tbl %>%
  filter(nonzero_coefficient) %>%
  mutate(direction = if_else(coefficient > 0, "positive", "negative")) %>%
  group_by(harmonized_cell_class, direction) %>%
  summarise(n_nonzero_genes = n_distinct(gene), .groups = "drop") %>%
  tidyr::complete(harmonized_cell_class, direction = c("positive", "negative"), fill = list(n_nonzero_genes = 0L)) %>%
  mutate(plot_n_nonzero_genes = abs(n_nonzero_genes))
p_direction <- ggplot(direction_data, aes(plot_n_nonzero_genes, direction, fill = direction)) +
  geom_col(width = 0.72) +
  facet_wrap(~ harmonized_cell_class, scales = "free_x") +
  scale_fill_manual(values = direction_palette, drop = FALSE) +
  theme_clock() +
  labs(x = "Non-zero genes", y = NULL, fill = "Coefficient", title = "Selected gene direction by cell type")
save_plot_with_data(p_direction, direction_data, "03_feature_selection", "selected_gene_direction_by_celltype",
                    "Selected gene direction by cell type", paths$coefficients, width = 10, height = 6.8)

top_n_features <- 6L
top_feature_data <- final_coef_tbl %>%
  filter(nonzero_coefficient, coefficient_direction %in% c("positive", "negative")) %>%
  group_by(harmonized_cell_class) %>%
  arrange(desc(abs_elasticnet_coefficient), gene, .by_group = TRUE) %>%
  slice_head(n = top_n_features) %>%
  mutate(
    feature_rank = row_number(),
    signed_coefficient = elasticnet_coefficient,
    plot_abs_coefficient = abs_elasticnet_coefficient,
    gene_plot_id = paste(harmonized_cell_class, feature_rank, gene, sep = "___")
  ) %>%
  ungroup() %>%
  mutate(gene_plot_id = factor(gene_plot_id, levels = rev(unique(gene_plot_id))))
gene_plot_labels <- function(x) sub("^.*___[0-9]+___", "", x)
p_top_features <- ggplot(top_feature_data, aes(plot_abs_coefficient, gene_plot_id, fill = coefficient_direction)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
  geom_col(width = 0.72) +
  facet_wrap(~ harmonized_cell_class, scales = "free_y", ncol = 3) +
  scale_y_discrete(labels = gene_plot_labels) +
  scale_fill_manual(values = direction_palette, drop = FALSE) +
  theme_clock(base_size = 8.2) +
  theme(
    axis.text.y = element_text(size = 6.1, color = "black"),
    strip.text = element_text(size = 7.2, face = "bold"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank()
  ) +
  labs(
    x = "Absolute final ElasticNet coefficient",
    y = NULL,
    fill = "Coefficient direction",
    title = "Top ElasticNet gene features by cell-type clock"
  )
save_plot_with_data(
  p_top_features, top_feature_data, "03_feature_selection", "top_elasticnet_gene_features_by_celltype",
  "Top ElasticNet gene features by cell-type clock", final_coef_table_path, width = 9.5, height = 6.8
)

# 04 GSEA
gsea_available <- file.exists(paths$gsea_results) &&
  file.exists(paths$gsea_display) &&
  nrow(fread(file = paths$gsea_results, nrows = 1)) > 0

gsea_yield <- cv_eval %>%
  select(harmonized_cell_class, n_metacells) %>%
  mutate(
    n_gsea_terms_analyzed = if (file.exists(paths$gsea_results)) {
      result_counts <- fread(file = paths$gsea_results) %>%
        as_tibble() %>%
        count(harmonized_cell_class, name = "n_gsea_terms_analyzed")
      result_counts$n_gsea_terms_analyzed[match(as.character(harmonized_cell_class), result_counts$harmonized_cell_class)] %>%
        tidyr::replace_na(0L)
    } else {
      0L
    },
    n_nominally_significant_terms = if (file.exists(paths$gsea_results)) {
      nominal_counts <- fread(file = paths$gsea_results) %>%
        as_tibble() %>%
        filter(!is.na(p_value), p_value <= 0.05) %>%
        count(harmonized_cell_class, name = "n_nominally_significant_terms")
      nominal_counts$n_nominally_significant_terms[match(as.character(harmonized_cell_class), nominal_counts$harmonized_cell_class)] %>%
        tidyr::replace_na(0L)
    } else {
      0L
    },
    gsea_status = if_else(gsea_available, "complete_elasticnet_coefficient_gsea", "not_available_no_metacell_cv_gsea_source_table")
  ) %>%
  arrange(n_gsea_terms_analyzed, harmonized_cell_class) %>%
  mutate(
    harmonized_cell_class = factor(harmonized_cell_class, levels = harmonized_cell_class),
    term_count_label = as.character(n_gsea_terms_analyzed),
    term_count_label_x = n_gsea_terms_analyzed + max(n_gsea_terms_analyzed, 1L) * 0.012
  )
p_gsea_yield <- ggplot(gsea_yield, aes(n_gsea_terms_analyzed, harmonized_cell_class)) +
  geom_col(fill = "#6F94C6", width = 0.72) +
  geom_text(aes(x = term_count_label_x, label = term_count_label), hjust = 0, size = 3.0, color = "grey20") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_clock(base_size = 8) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(color = "black")
  ) +
  labs(
    x = "GSEA terms analyzed",
    y = "Cell type",
    title = "GSEA analyzed term yield"
  )
save_plot_with_data(
  p_gsea_yield, gsea_yield, "04_gsea", "gsea_eligibility_and_display_term_yield",
  "GSEA analyzed term yield",
  if (file.exists(paths$gsea_results)) paths$gsea_results else paths$coefficients,
  width = 5.4, height = 3.6,
  status = if (file.exists(paths$gsea_rank)) "complete" else "complete_no_gsea_source_table"
)

pretty_gsea_term <- function(x) {
  x <- sub("^GOBP_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  replacements <- c(
    "\\bDna\\b" = "DNA",
    "\\bRna\\b" = "RNA",
    "\\bMrna\\b" = "mRNA",
    "\\bAtp\\b" = "ATP",
    "\\bNad\\b" = "NAD",
    "\\bIl\\b" = "IL",
    "\\bTnf\\b" = "TNF",
    "\\bNfkb\\b" = "NFKB",
    "\\bP53\\b" = "P53",
    "\\bTgf\\b" = "TGF",
    "\\bMhc\\b" = "MHC"
  )
  for (pattern in names(replacements)) x <- gsub(pattern, replacements[[pattern]], x)
  x
}

if (gsea_available) {
  gsea_display <- fread(paths$gsea_display) %>%
    as_tibble() %>%
    mutate(
      harmonized_cell_class = factor(harmonized_cell_class, levels = cv_eval$harmonized_cell_class),
      is_status_row = display_status != "display_term" | is.na(NES),
      NES_direction_display = case_when(
        NES > 0 ~ "Positive NES = older-predictive clock-gene programs",
        NES < 0 ~ "Negative NES = younger-predictive clock-gene programs",
        TRUE ~ "GSEA status"
      ),
      term_display = if_else(is_status_row, term_name, pretty_gsea_term(term_name))
    ) %>%
    group_by(harmonized_cell_class) %>%
    arrange(is_status_row, NES, p_adjust, p_value, .by_group = TRUE) %>%
    mutate(
      display_order = row_number(),
      term_plot_id = paste(harmonized_cell_class, display_order, term_display, sep = "___"),
      label_x = NES + if_else(NES >= 0, -0.06, 0.06),
      label_hjust = if_else(NES >= 0, 1, 0)
    ) %>%
    ungroup()
  gsea_slot_count <- max(20L, max(gsea_display$display_order, na.rm = TRUE))
  gsea_padding <- gsea_display %>%
    count(harmonized_cell_class, name = "n_existing_slots") %>%
    filter(n_existing_slots < gsea_slot_count) %>%
    group_by(harmonized_cell_class, n_existing_slots) %>%
    summarise(display_order = list(seq.int(n_existing_slots + 1L, gsea_slot_count)), .groups = "drop") %>%
    tidyr::unnest(display_order) %>%
    mutate(
      gene_set_source = NA_character_,
      term_id = NA_character_,
      term_name = "",
      pathway_size = NA_integer_,
      NES = NA_real_,
      ES = NA_real_,
      p_value = NA_real_,
      p_adjust = NA_real_,
      leading_edge_genes = NA_character_,
      n_ranked_genes = NA_integer_,
      rank_source = NA_character_,
      NES_direction = "padding",
      evidence_label = "padding",
      p_label = "",
      display_rank = NA_integer_,
      display_status = "padding_empty_slot",
      is_status_row = FALSE,
      NES_direction_display = NA_character_,
      term_display = "",
      term_plot_id = paste(harmonized_cell_class, display_order, "empty_slot", sep = "___"),
      label_x = NA_real_,
      label_hjust = 0.5,
      is_padding_row = TRUE
    )
  gsea_display <- gsea_display %>%
    mutate(is_padding_row = FALSE) %>%
    bind_rows(gsea_padding) %>%
    arrange(harmonized_cell_class, display_order)
  gsea_display$term_plot_id <- factor(gsea_display$term_plot_id, levels = unique(gsea_display$term_plot_id))
  gsea_term_labels <- function(x) {
    labels <- sub("^.*___[0-9]+___", "", x)
    if_else(labels == "empty_slot", "", str_wrap(labels, width = 38))
  }

  p_gsea_nes <- ggplot(gsea_display, aes(NES, term_plot_id)) +
    geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
    geom_col(
      data = gsea_display %>% filter(!is_status_row, !is_padding_row),
      aes(fill = NES_direction_display),
      width = 0.68
    ) +
    geom_text(
      data = gsea_display %>% filter(!is_status_row, !is_padding_row),
      aes(x = label_x, label = p_label, hjust = label_hjust),
      size = 2.1,
      color = "grey20"
    ) +
    geom_text(
      data = gsea_display %>% filter(is_status_row, !is_padding_row),
      aes(x = 0, label = str_wrap(term_display, width = 38)),
      size = 2.5,
      color = "grey35"
    ) +
    facet_wrap(~ harmonized_cell_class, scales = "free_y", ncol = 2, drop = FALSE) +
    scale_y_discrete(labels = gsea_term_labels) +
    scale_fill_manual(values = c(
      "Positive NES = older-predictive clock-gene programs" = "#B2182B",
      "Negative NES = younger-predictive clock-gene programs" = "#2166AC"
    ), drop = FALSE) +
    theme_clock(base_size = 8) +
    theme(
      axis.text.y = element_text(size = 6.2),
      legend.position = "bottom",
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25)
    ) +
    labs(
      x = "NES",
      y = NULL,
      fill = "Coefficient-ranked enrichment direction",
      title = "ElasticNet coefficient-ranked GSEA by cell type"
    )
  save_plot_with_data(
    p_gsea_nes, gsea_display, "04_gsea", "elasticnet_coefficient_gsea_nes_panels_by_celltype",
    "ElasticNet coefficient-ranked GSEA by cell type",
    c(paths$gsea_results, paths$gsea_display),
    width = 9.6, height = 17.2, status = "complete"
  )
  file.copy(
    file.path(figure_root, "04_gsea", "elasticnet_coefficient_gsea_nes_panels_by_celltype.png"),
    file.path(figure_root, "04_gsea", "gsea_nes_panels_by_celltype.png"),
    overwrite = TRUE
  )
  file.copy(
    file.path(figure_root, "04_gsea", "elasticnet_coefficient_gsea_nes_panels_by_celltype.pdf"),
    file.path(figure_root, "04_gsea", "gsea_nes_panels_by_celltype.pdf"),
    overwrite = TRUE
  )
  file.copy(
    file.path(plot_data_root, "04_gsea", "elasticnet_coefficient_gsea_nes_panels_by_celltype_plot_data.tsv"),
    file.path(plot_data_root, "04_gsea", "gsea_nes_panels_by_celltype_plot_data.tsv"),
    overwrite = TRUE
  )
} else {
  gsea_status_data <- gsea_yield %>%
    transmute(
      harmonized_cell_class = as.character(harmonized_cell_class),
      term = "No metacell-CV GSEA source table",
      NES = 0,
      p_value = NA_real_,
      direction = "not_available",
      status = gsea_status
    )
  p_gsea_nes <- ggplot(gsea_status_data, aes(NES, harmonized_cell_class)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(color = "grey35", size = 2) +
    theme_clock() +
    labs(
      x = "NES",
      y = "Cell type",
      title = "GSEA NES panels by cell type",
      subtitle = "No canonical metacell-CV GSEA result table was found"
    )
  save_plot_with_data(
    p_gsea_nes, gsea_status_data, "04_gsea", "gsea_nes_panels_by_celltype",
    "GSEA NES panels by cell type", paths$coefficients,
    width = 7.2, height = 4.8, status = "not_available_no_source_table"
  )
}

# 05 projection stress-test plots
clock_quality <- cv_eval %>%
  transmute(
    harmonized_cell_class,
    cv_MAE_weeks = mae_weeks,
    cv_R2 = r2,
    cv_slope = slope,
    n_metacells_internal_cv = n_metacells,
    model_status = case_when(
      is.na(cv_R2) | is.na(cv_MAE_weeks) | n_metacells_internal_cv < 50 ~ "diagnostic_only",
      cv_R2 >= 0.80 & cv_MAE_weeks <= 8.5 & cv_slope >= 0.70 ~ "strong_internal_clock",
      cv_R2 >= 0.60 & cv_MAE_weeks <= 12.5 & cv_slope >= 0.50 ~ "usable_internal_clock",
      cv_R2 >= 0.45 & cv_MAE_weeks <= 14.5 ~ "weak_internal_clock",
      TRUE ~ "diagnostic_only"
    ),
    projection_interpretation_flag = case_when(
      model_status %in% c("strong_internal_clock", "usable_internal_clock") ~ "confident_projection_context",
      model_status == "weak_internal_clock" ~ "cautious_projection_context",
      TRUE ~ "diagnostic_only_projection_context"
    )
  )

ensure_model_quality <- function(tbl) {
  required_quality <- c(
    "cv_MAE_weeks", "cv_R2", "cv_slope", "n_metacells_internal_cv",
    "model_status", "projection_interpretation_flag"
  )
  if (all(required_quality %in% names(tbl))) {
    tbl
  } else {
    tbl %>% left_join(clock_quality, by = "harmonized_cell_class")
  }
}

add_projection_labels <- function(tbl) {
  ensure_model_quality(tbl) %>%
    mutate(
      status_marker = if_else(model_status %in% c("weak_internal_clock", "diagnostic_only"), "*", ""),
      mae_label = sprintf("MAE=%.1fw", cv_MAE_weeks),
      model_quality_label = sprintf(
        "%s%s | %s; R2=%.2f; slope=%.2f",
        harmonized_cell_class, status_marker, mae_label, cv_R2, cv_slope
      ),
      contrast_label = case_when(
        stress_test_family == "VCD" ~ sprintf("Injection %s; post %sd", injection_age, post_injection_days),
        stress_test_family == "Foxl2" ~ paste("Age", gsub("_", " ", age_group)),
        TRUE ~ matched_contrast_id
      )
    )
}

write_projection_quality_data <- function(summary_tbl, filename) {
  out_path <- file.path(plot_data_root, "05_stress_test_projection", filename)
  out_data <- add_projection_labels(summary_tbl) %>%
    arrange(stress_test_family, matched_contrast_id, desc(abs(delta_predicted_age_weeks)), harmonized_cell_class)
  fwrite(out_data, out_path, sep = "\t")
  out_data
}

stress_plot_summary <- function(summary_tbl) {
  add_projection_labels(summary_tbl) %>%
    group_by(stress_test_family, matched_contrast_id, contrast_label, harmonized_cell_class) %>%
    summarise(
      control_sample_ids = first(control_sample_ids),
      case_sample_ids = paste(sort(unique(case_sample_id)), collapse = ";"),
      age_group = first(age_group),
      injection_age = first(injection_age),
      post_injection_days = first(post_injection_days),
      mean_delta_predicted_age_weeks = mean(delta_predicted_age_weeks, na.rm = TRUE),
      median_delta_predicted_age_weeks = median(delta_predicted_age_weeks, na.rm = TRUE),
      n_case_samples = n_distinct(case_sample_id),
      n_control_samples = first(n_control_samples),
      n_case_metacells = sum(n_case_metacells, na.rm = TRUE),
      n_control_metacells = first(n_control_metacells),
      cv_MAE_weeks = first(cv_MAE_weeks),
      cv_R2 = first(cv_R2),
      cv_slope = first(cv_slope),
      n_metacells_internal_cv = first(n_metacells_internal_cv),
      model_status = first(model_status),
      projection_interpretation_flag = first(projection_interpretation_flag),
      mae_label = first(mae_label),
      model_quality_label = first(model_quality_label),
      .groups = "drop"
    ) %>%
    mutate(
      plot_scope = "all_cell_types_exploratory",
      signed_delta_label_x = median_delta_predicted_age_weeks + 0.6,
      label_hjust = 0
    )
}

vcd_delta_with_quality <- write_projection_quality_data(
  vcd_summary,
  "vcd_projection_delta_by_contrast_with_model_quality.tsv"
)
foxl2_delta_with_quality <- write_projection_quality_data(
  foxl2_summary,
  "foxl2_projection_delta_by_age_stratum_with_model_quality.tsv"
)

stress_sample_metadata <- prjna863443_metadata %>%
  mutate(
    age_group = clean_token(age_label),
    stress_test_family = case_when(
      model == "VCD" ~ "VCD",
      str_detect(model, "^Foxl2") ~ "Foxl2",
      TRUE ~ NA_character_
    ),
    case_control_status = case_when(
      model == "VCD" & condition == "VCD" ~ "case",
      model == "VCD" & condition == "CTL" ~ "control",
      str_detect(model, "^Foxl2") & genotype == "Foxl2+/-" ~ "case",
      str_detect(model, "^Foxl2") & genotype == "Foxl2+/+" ~ "control",
      TRUE ~ NA_character_
    ),
    treatment_or_genotype = case_when(
      model == "VCD" ~ condition,
      str_detect(model, "^Foxl2") ~ genotype,
      TRUE ~ treatment
    ),
    matched_contrast_id = case_when(
      model == "VCD" ~ paste0("VCD_injection_", clean_token(injection_age), "_post_", clean_token(post_injection_days)),
      str_detect(model, "^Foxl2") ~ paste0("Foxl2_age_", age_group),
      TRUE ~ NA_character_
    ),
    contrast_label = case_when(
      model == "VCD" ~ sprintf("VCD: injection %s; post %sd", injection_age, post_injection_days),
      str_detect(model, "^Foxl2") ~ paste("Foxl2: age", gsub("_", " ", age_group)),
      TRUE ~ NA_character_
    ),
    injection_age_num = suppressWarnings(as.numeric(str_extract(injection_age, "\\d+"))),
    post_injection_days_num = suppressWarnings(as.numeric(post_injection_days)),
    sample_display = paste0(sample_id, " (", case_control_status, ")")
  ) %>%
  filter(stress_test_family %in% c("VCD", "Foxl2"), case_control_status %in% c("case", "control")) %>%
  select(
    sample_id, sample_display, stress_test_family, case_control_status,
    treatment_or_genotype, age_weeks, age_group, injection_age,
    post_injection_days, matched_contrast_id, contrast_label,
    injection_age_num, post_injection_days_num
  )

stress_celltype_levels <- eligibility %>%
  arrange(desc(n_metacells)) %>%
  pull(harmonized_cell_class) %>%
  as.character()
stress_sample_levels <- stress_sample_metadata %>%
  mutate(
    family_order = match(stress_test_family, c("VCD", "Foxl2")),
    status_order = match(case_control_status, c("control", "case"))
  ) %>%
  arrange(family_order, injection_age_num, post_injection_days_num, age_weeks, status_order, sample_id) %>%
  pull(sample_display)
stress_contrast_levels <- stress_sample_metadata %>%
  distinct(stress_test_family, matched_contrast_id, contrast_label, injection_age_num, post_injection_days_num, age_weeks) %>%
  mutate(family_order = match(stress_test_family, c("VCD", "Foxl2"))) %>%
  arrange(family_order, injection_age_num, post_injection_days_num, age_weeks, matched_contrast_id) %>%
  pull(contrast_label)

stress_metacell_annotated <- projection_metacell %>%
  inner_join(
    stress_sample_metadata %>%
      transmute(
        sample_id,
        sample_display,
        stress_test_family,
        case_control_status,
        treatment_or_genotype,
        sample_age_weeks = age_weeks,
        age_group,
        injection_age,
        post_injection_days,
        matched_contrast_id,
        injection_age_num,
        post_injection_days_num
      ),
    by = "sample_id"
  ) %>%
  ensure_model_quality()

stress_metacell_control_means <- stress_metacell_annotated %>%
  filter(case_control_status == "control") %>%
  group_by(stress_test_family, matched_contrast_id, harmonized_cell_class) %>%
  summarise(
    matched_control_mean_predicted_age_weeks = mean(predicted_age_weeks, na.rm = TRUE),
    n_control_metacells = n(),
    n_control_samples = n_distinct(sample_id),
    control_sample_ids = paste(sort(unique(sample_id)), collapse = ";"),
    .groups = "drop"
  )

stress_metacell_delta_data <- stress_metacell_annotated %>%
  left_join(
    stress_metacell_control_means,
    by = c("stress_test_family", "matched_contrast_id", "harmonized_cell_class")
  ) %>%
  group_by(stress_test_family, matched_contrast_id, harmonized_cell_class) %>%
  mutate(
    n_case_metacells = sum(case_control_status == "case"),
    n_control_metacells = first(n_control_metacells),
    n_case_samples = n_distinct(sample_id[case_control_status == "case"]),
    n_control_samples = first(n_control_samples)
  ) %>%
  ungroup() %>%
  mutate(
    case_sample_id = sample_id,
    case_metacell_id = metacell_id,
    projection_metacell_id = metacell_id,
    predicted_age_metacell_weeks = predicted_age_weeks,
    predicted_age_case_metacell_weeks = if_else(case_control_status == "case", predicted_age_weeks, NA_real_),
    projection_group = case_when(
      stress_test_family == "VCD" & case_control_status == "case" ~ "vcd",
      stress_test_family == "VCD" & case_control_status == "control" ~ "vcd_control",
      stress_test_family == "Foxl2" & case_control_status == "case" ~ "foxl2",
      stress_test_family == "Foxl2" & case_control_status == "control" ~ "foxl2_control",
      TRUE ~ case_control_status
    ),
    delta_predicted_age_weeks = predicted_age_metacell_weeks - matched_control_mean_predicted_age_weeks,
    delta_definition = "predicted_age_metacell - mean(predicted_age_matched_control_metacells_within_same_matched_contrast_and_celltype)",
    plot_scope = "matched_contrast_specific_metacell_delta_with_controls"
  ) %>%
  add_projection_labels() %>%
  mutate(
    contrast_label = factor(contrast_label, levels = unique(contrast_label)),
    harmonized_cell_class = factor(harmonized_cell_class, levels = stress_celltype_levels),
    case_sample_id = factor(case_sample_id, levels = unique(case_sample_id)),
    projection_group = factor(projection_group, levels = c("vcd", "vcd_control", "foxl2", "foxl2_control"))
  ) %>%
  group_by(stress_test_family, matched_contrast_id, harmonized_cell_class, projection_group) %>%
  mutate(n_metacells_projection_group_celltype_contrast = n()) %>%
  group_by(stress_test_family, harmonized_cell_class, projection_group) %>%
  mutate(n_metacells_projection_group_celltype_total = n()) %>%
  ungroup() %>%
  arrange(stress_test_family, matched_contrast_id, harmonized_cell_class, case_sample_id, case_metacell_id)

availability_grid <- expand.grid(
  sample_id = stress_sample_metadata$sample_id,
  harmonized_cell_class = stress_celltype_levels,
  stringsAsFactors = FALSE
) %>%
  as_tibble()

stress_availability_data <- availability_grid %>%
  left_join(stress_sample_metadata, by = "sample_id") %>%
  left_join(
    projection_sample %>%
      select(sample_id, harmonized_cell_class, n_projected_metacells, n_projected_cells),
    by = c("sample_id", "harmonized_cell_class")
  ) %>%
  mutate(
    n_projected_metacells = if_else(is.na(n_projected_metacells), 0L, as.integer(n_projected_metacells)),
    n_projected_cells = if_else(is.na(n_projected_cells), 0L, as.integer(n_projected_cells)),
    metacell_available = n_projected_metacells > 0,
    sample_display = factor(sample_display, levels = rev(stress_sample_levels)),
    harmonized_cell_class = factor(harmonized_cell_class, levels = stress_celltype_levels),
    contrast_label = factor(contrast_label, levels = stress_contrast_levels)
  ) %>%
  arrange(stress_test_family, matched_contrast_id, sample_id, harmonized_cell_class)

stress_availability_summary <- stress_availability_data %>%
  group_by(stress_test_family, harmonized_cell_class) %>%
  summarise(
    n_projected_metacells = sum(n_projected_metacells, na.rm = TRUE),
    n_projected_cells = sum(n_projected_cells, na.rm = TRUE),
    n_stress_test_samples = n_distinct(sample_id),
    n_samples_with_metacells = n_distinct(sample_id[metacell_available]),
    .groups = "drop"
  ) %>%
  mutate(
    harmonized_cell_class = factor(harmonized_cell_class, levels = rev(stress_celltype_levels)),
    stress_test_family = factor(stress_test_family, levels = c("VCD", "Foxl2"))
  )

stress_availability_summary_labels <- stress_availability_summary %>%
  group_by(harmonized_cell_class) %>%
  summarise(
    total_projected_metacells = sum(n_projected_metacells, na.rm = TRUE),
    n_samples_with_metacells = sum(n_samples_with_metacells, na.rm = TRUE),
    n_stress_test_samples = sum(n_stress_test_samples, na.rm = TRUE),
    label = sprintf("%s/%s samples", n_samples_with_metacells, n_stress_test_samples),
    .groups = "drop"
  ) %>%
  mutate(label_x = total_projected_metacells + max(total_projected_metacells, na.rm = TRUE) * 0.02)
stress_availability_summary_xlim <- max(stress_availability_summary_labels$total_projected_metacells, na.rm = TRUE) * 1.36

p_stress_availability_summary <- ggplot(
  stress_availability_summary,
  aes(n_projected_metacells, harmonized_cell_class, fill = stress_test_family)
) +
  geom_col(width = 0.72) +
  geom_text(
    data = stress_availability_summary_labels,
    aes(x = label_x, y = harmonized_cell_class, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.65,
    color = "grey25"
  ) +
  scale_fill_manual(values = c(VCD = "#4E79A7", Foxl2 = "#59A14F"), drop = FALSE) +
  scale_x_continuous(limits = c(0, stress_availability_summary_xlim), expand = expansion(mult = c(0, 0))) +
  theme_clock() +
  labs(
    x = "Projected metacells across stress-test samples",
    y = "Cell type",
    fill = "Stress-test family",
    title = "Stress-test metacell availability summary",
    subtitle = "Summed PRJNA863443 VCD and Foxl2 projected metacells; labels show sample-cell-type coverage."
  )
save_plot_with_data(
  p_stress_availability_summary,
  stress_availability_summary,
  "05_stress_test_projection",
  "stress_test_metacell_availability_summary",
  "Stress-test metacell availability summary",
  c(paths$projection_sample, paths$prjna863443_metadata),
  width = 7.4,
  height = 4.9
)

p_stress_availability <- ggplot(stress_availability_data, aes(harmonized_cell_class, sample_display, fill = n_projected_metacells)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = if_else(n_projected_metacells > 0L, as.character(n_projected_metacells), "")),
            size = 2.05, color = "grey15") +
  facet_wrap(~ contrast_label, ncol = 3, scales = "free_y") +
  scale_fill_gradient(low = "#F2F2F2", high = "#2166AC", trans = "sqrt", breaks = c(0, 5, 25, 100)) +
  theme_clock(base_size = 7.8) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.spacing.x = grid::unit(0.45, "lines"),
    panel.spacing.y = grid::unit(0.55, "lines"),
    legend.position = "right"
  ) +
  labs(
    x = "Cell type",
    y = "Stress-test sample",
    fill = "Projected\nmetacells",
    title = "Stress-test metacell availability by sample and cell type",
    subtitle = "Matched PRJNA863443 VCD and Foxl2 samples; blank tiles indicate no projected metacells for that sample-cell-type pair."
  )
save_plot_with_data(
  p_stress_availability,
  stress_availability_data,
  "05_stress_test_projection",
  "stress_test_metacell_availability_by_sample",
  "Stress-test metacell availability by sample and cell type",
  c(paths$projection_sample, paths$prjna863443_metadata),
  width = 11.8,
  height = 8.8
)

metacell_n_labels <- function(dat, x_var, facet_vars = character()) {
  group_vars <- c(facet_vars, x_var, "projection_group")
  labels <- dat %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_metacells = n(),
      y_min = min(delta_predicted_age_weeks, na.rm = TRUE),
      y_max = max(delta_predicted_age_weeks, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(n_label = paste0("n=", n_metacells))
  if (length(facet_vars) > 0L) {
    labels %>%
      group_by(across(all_of(facet_vars))) %>%
      mutate(
        facet_y_min = min(y_min, na.rm = TRUE),
        facet_y_max = max(y_max, na.rm = TRUE),
        label_y = facet_y_max + pmax((facet_y_max - facet_y_min) * 0.08, 0.06)
      ) %>%
      ungroup()
  } else {
    facet_y_min <- min(labels$y_min, na.rm = TRUE)
    facet_y_max <- max(labels$y_max, na.rm = TRUE)
    labels %>%
      mutate(label_y = facet_y_max + max((facet_y_max - facet_y_min) * 0.08, 0.06, na.rm = TRUE))
  }
}

delta_age_by_pair_and_celltype <- function(dat, figure_id, title, source_path) {
  dat <- dat %>%
    mutate(
      delta_summary = sprintf(
        "control metacell mean=%0.1f; metacell=%0.1f",
        matched_control_mean_predicted_age_weeks,
        predicted_age_metacell_weeks
      )
    ) %>%
    arrange(matched_contrast_id, harmonized_cell_class, case_sample_id, case_metacell_id)
  count_labels <- metacell_n_labels(dat, "case_sample_id", c("harmonized_cell_class", "contrast_label"))
  p <- ggplot(dat, aes(case_sample_id, delta_predicted_age_weeks)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_boxplot(
      aes(fill = projection_group),
      width = stress_box_width,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.28,
      alpha = 0.88
    ) +
    geom_point(
      aes(fill = projection_group, color = projection_group),
      position = position_jitter(width = 0.12, height = 0, seed = 7),
      shape = 21,
      size = 0.42,
      alpha = 0.55,
      stroke = 0.16
    ) +
    geom_text(
      data = count_labels,
      aes(x = case_sample_id, y = label_y, label = n_label, color = projection_group),
      inherit.aes = FALSE,
      size = 1.45,
      vjust = 0,
      show.legend = FALSE
    ) +
    facet_grid(harmonized_cell_class ~ contrast_label, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = projection_group_palette, drop = TRUE) +
    scale_color_manual(values = projection_group_outline_palette, drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
    theme_clock(base_size = 6.8) +
    theme(
      axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1, size = 5.4),
      strip.text.y = element_text(angle = 0, size = 6.2),
      panel.spacing.x = grid::unit(0.32, "lines"),
      panel.spacing.y = grid::unit(0.24, "lines"),
      legend.position = "right"
    ) +
    labs(
      x = "Sample within matched contrast",
      y = "Delta predicted age vs matched controls (weeks)",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = "Each point is one projected metacell; control and case deltas use the same matched-control metacell mean within contrast and cell type."
    )
  save_plot_with_data(
    p,
    dat,
    "05_stress_test_projection",
    figure_id,
    title,
    source_path,
    width = 12.8,
    height = 9.2
  )
}

vcd_metacell_delta_data <- stress_metacell_delta_data %>%
  filter(stress_test_family == "VCD")
foxl2_metacell_delta_data <- stress_metacell_delta_data %>%
  filter(stress_test_family == "Foxl2")
metacell_delta_source_paths <- c(paths$projection_metacell, paths$prjna863443_metadata)

delta_age_by_pair_and_celltype(
  vcd_metacell_delta_data,
  "vcd_delta_age_by_pair_and_celltype",
  "VCD metacell delta predicted age by matched pair and cell type",
  metacell_delta_source_paths
)
delta_age_by_pair_and_celltype(
  foxl2_metacell_delta_data,
  "foxl2_delta_age_by_pair_and_celltype",
  "Foxl2 metacell delta predicted age by matched age stratum and cell type",
  metacell_delta_source_paths
)

projection_heatmap_by_contrast <- function(summary_tbl, figure_id, title, source_path) {
  dat <- stress_plot_summary(summary_tbl) %>%
    mutate(
      celltype_label = reorder(harmonized_cell_class, median_delta_predicted_age_weeks, median, na.rm = TRUE),
      contrast_label = factor(contrast_label, levels = unique(contrast_label))
    )
  p <- ggplot(dat, aes(celltype_label, contrast_label, fill = median_delta_predicted_age_weeks)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_point(aes(shape = projection_interpretation_flag), color = "black", size = 1.7) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    scale_shape_manual(values = interpretation_shapes, drop = FALSE) +
    theme_clock() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1)) +
    labs(
      x = "Cell type",
      y = "Matched contrast",
      fill = "Median case-control\ndelta (weeks)",
      shape = "Projection context",
      title = title,
      subtitle = "Each tile summarizes case-sample deltas computed only against controls from the same matched contrast and cell type."
    )
  save_plot_with_data(p, dat, "05_stress_test_projection", figure_id,
                      title, source_path, width = 8.6, height = 5.4)
}

projection_ranked_by_contrast <- function(summary_tbl, figure_id, title, source_path, confidence_filtered = FALSE) {
  dat <- stress_plot_summary(summary_tbl) %>%
    { if (confidence_filtered) filter(., model_status %in% c("strong_internal_clock", "usable_internal_clock")) else . } %>%
    mutate(
      plot_scope = if_else(confidence_filtered, "strong_or_usable_internal_clock_only", "all_cell_types_exploratory"),
      contrast_label = factor(contrast_label, levels = unique(contrast_label)),
      celltype_label = reorder(harmonized_cell_class, median_delta_predicted_age_weeks)
    )
  subtitle <- if (confidence_filtered) {
    "Confidence-filtered: only strong_internal_clock and usable_internal_clock models are shown."
  } else {
    "All cell types shown; exploratory. Point color, size, and labels display internal metacell-CV quality."
  }
  p <- ggplot(dat, aes(median_delta_predicted_age_weeks, celltype_label)) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.35) +
    geom_segment(aes(x = 0, xend = median_delta_predicted_age_weeks, yend = celltype_label),
                 color = "grey72", linewidth = 0.45) +
    geom_point(
      aes(size = cv_R2, color = model_status, shape = projection_interpretation_flag, alpha = projection_interpretation_flag),
      stroke = 0.75
    ) +
    geom_text(aes(x = signed_delta_label_x, label = mae_label, hjust = label_hjust), size = 2.25, color = "grey20") +
    facet_wrap(~ contrast_label, scales = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0.34, 0.42))) +
    scale_size_continuous(limits = c(0, 1), range = c(2.0, 5.0), breaks = c(0.5, 0.7, 0.9)) +
    scale_color_manual(values = status_palette, drop = FALSE) +
    scale_shape_manual(values = interpretation_shapes, drop = FALSE) +
    scale_alpha_manual(values = interpretation_alpha, drop = FALSE, guide = "none") +
    theme_clock(base_size = 8.2) +
    theme(legend.position = "right") +
    labs(
      x = "Median case-control delta predicted age (weeks)",
      y = "Cell type",
      size = "Internal CV R2",
      color = "Internal clock status",
      shape = "Projection context",
      title = title,
      subtitle = subtitle
    )
  save_plot_with_data(p, dat, "05_stress_test_projection", figure_id,
                      title, source_path, width = 9.4, height = 6.4)
}

projection_heatmap_by_contrast(
  vcd_summary,
  "vcd_projection_delta_heatmap_by_contrast",
  "VCD projection median delta heatmap by matched contrast",
  paths$vcd_summary
)
projection_heatmap_by_contrast(
  vcd_summary,
  "vcd_delta_heatmap_by_pair",
  "VCD median delta heatmap by matched pair",
  paths$vcd_summary
)

projection_heatmap_by_contrast(
  foxl2_summary,
  "foxl2_projection_delta_heatmap_by_age_stratum",
  "Foxl2 projection median delta heatmap by age stratum",
  paths$foxl2_summary
)
projection_heatmap_by_contrast(
  foxl2_summary,
  "foxl2_delta_heatmap_by_age_stratum",
  "Foxl2 median delta heatmap by age stratum",
  paths$foxl2_summary
)

delta_distribution_data <- function(metacell_delta_tbl) {
  metacell_delta_tbl %>%
    mutate(
      contrast_label = factor(contrast_label, levels = unique(contrast_label)),
      harmonized_cell_class = factor(harmonized_cell_class, levels = stress_celltype_levels),
      case_sample_id = factor(case_sample_id, levels = unique(case_sample_id))
    ) %>%
    arrange(matched_contrast_id, harmonized_cell_class, case_sample_id, case_metacell_id)
}

plot_delta_distribution <- function(dat, figure_id, title, source_path, width = 7.2, height = 4.8) {
  count_labels <- metacell_n_labels(dat, "contrast_label")
  p <- ggplot(dat, aes(contrast_label, delta_predicted_age_weeks, fill = projection_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
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
      position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.68, seed = 17),
      shape = 21,
      size = 0.42,
      alpha = 0.55,
      stroke = 0.16
    ) +
    geom_text(
      data = count_labels,
      aes(x = contrast_label, y = label_y, label = n_label, color = projection_group, group = projection_group),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.68),
      size = 2.2,
      vjust = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = projection_group_palette, drop = TRUE) +
    scale_color_manual(values = projection_group_outline_palette, drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
    theme_clock(base_size = 8.2) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      legend.position = "right"
    ) +
    labs(
      x = "Matched contrast / age stratum",
      y = "Delta predicted age vs matched controls (weeks)",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = "Each point is one projected metacell; boxes compare case and matched-control projection groups within each contrast."
    )
  save_plot_with_data(p, dat, "05_stress_test_projection", figure_id,
                      title, source_path, width = width, height = height)
}

plot_combined_delta_distribution <- function(metacell_delta_tbl, family, figure_id, title, source_path) {
  dat <- delta_distribution_data(metacell_delta_tbl) %>%
    mutate(harmonized_cell_class = factor(harmonized_cell_class, levels = stress_celltype_levels))
  count_labels <- metacell_n_labels(dat, "contrast_label", "harmonized_cell_class")
  p <- ggplot(dat, aes(contrast_label, delta_predicted_age_weeks, fill = projection_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_boxplot(
      width = stress_box_width,
      outlier.shape = NA,
      color = "black",
      linewidth = 0.26,
      alpha = 0.88,
      position = position_dodge2(width = 0.68, preserve = "single")
    ) +
    geom_point(
      aes(fill = projection_group, color = projection_group),
      position = position_jitterdodge(jitter.width = 0.10, jitter.height = 0, dodge.width = 0.68, seed = 23),
      shape = 21,
      size = 0.34,
      alpha = 0.50,
      stroke = 0.14
    ) +
    geom_text(
      data = count_labels,
      aes(x = contrast_label, y = label_y, label = n_label, color = projection_group, group = projection_group),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.68),
      size = 1.65,
      vjust = 0,
      show.legend = FALSE
    ) +
    facet_wrap(~ harmonized_cell_class, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = projection_group_palette, drop = TRUE) +
    scale_color_manual(values = projection_group_outline_palette, drop = TRUE) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
    theme_clock(base_size = 7.4) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 5.8),
      panel.spacing.x = grid::unit(0.34, "lines"),
      panel.spacing.y = grid::unit(0.38, "lines"),
      legend.position = "right"
    ) +
    labs(
      x = if_else(family == "VCD", "VCD matched contrast", "Foxl2 matched age stratum"),
      y = "Delta predicted age vs matched controls (weeks)",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = if_else(
        family == "VCD",
        "VCD and VCD-control metacell deltas across 4 injection-age/post-injection matched contrasts.",
        "Foxl2 and Foxl2-control metacell deltas across 2 matched age strata."
      )
    )
  save_plot_with_data(p, dat, "05_stress_test_projection", figure_id,
                      title, source_path, width = 11.2, height = 7.6)
}

write_per_celltype_delta_distributions <- function(metacell_delta_tbl, family, source_path) {
  all_dat <- delta_distribution_data(metacell_delta_tbl)
  celltypes <- sort(unique(as.character(all_dat$harmonized_cell_class)))
  for (celltype in celltypes) {
    dat <- all_dat %>%
      filter(as.character(harmonized_cell_class) == celltype) %>%
      mutate(harmonized_cell_class = as.character(harmonized_cell_class))
    if (family == "VCD") {
      figure_id <- paste0("vcd_delta_distribution_", safe_id(celltype), "_by_matched_contrast")
      title <- paste0("VCD delta-age distribution by matched contrast: ", celltype)
    } else {
      figure_id <- paste0("foxl2_delta_distribution_", safe_id(celltype), "_by_age_stratum")
      title <- paste0("Foxl2 delta-age distribution by age stratum: ", celltype)
    }
    plot_delta_distribution(dat, figure_id, title, source_path, width = 7.2, height = 4.8)
  }
}

write_per_celltype_delta_distributions(vcd_metacell_delta_data, "VCD", metacell_delta_source_paths)
write_per_celltype_delta_distributions(foxl2_metacell_delta_data, "Foxl2", metacell_delta_source_paths)
plot_combined_delta_distribution(
  vcd_metacell_delta_data,
  "VCD",
  "vcd_delta_distribution_by_celltype_combined",
  "VCD delta-age distributions by cell type and matched contrast",
  metacell_delta_source_paths
)
plot_combined_delta_distribution(
  foxl2_metacell_delta_data,
  "Foxl2",
  "foxl2_delta_distribution_by_celltype_combined",
  "Foxl2 delta-age distributions by cell type and matched age stratum",
  metacell_delta_source_paths
)

save_legacy_plot_with_data <- function(plot, data, section, figure_id, source_table_path,
                                       width = 8, height = 5) {
  png_path <- file.path(figure_root, section, paste0(figure_id, ".png"))
  pdf_path <- file.path(figure_root, section, paste0(figure_id, ".pdf"))
  data_path <- file.path(plot_data_root, section, paste0(figure_id, "_plot_data.tsv"))
  fwrite(data, data_path, sep = "\t")
  ggsave(pdf_path, plot, width = width, height = height)
  ggsave(png_path, plot, width = width, height = height, dpi = 240)
  invisible(list(png_path = png_path, pdf_path = pdf_path, data_path = data_path, source_table_path = source_table_path))
}

legacy_no_strata_metacell_data <- function(family) {
  raw_groups <- if (family == "VCD") c("vcd_model", "vcd_control") else c("foxl2_model", "foxl2_control")
  plot_groups <- if (family == "VCD") c("vcd", "vcd_control") else c("foxl2", "foxl2_control")
  projection_metacell %>%
    filter(intervention_class %in% raw_groups) %>%
    ensure_model_quality() %>%
    mutate(
      projection_group = case_when(
        intervention_class == "vcd_model" ~ "vcd",
        intervention_class == "foxl2_model" ~ "foxl2",
        TRUE ~ intervention_class
      ),
      projection_group = factor(projection_group, levels = plot_groups),
      projection_metacell_id = metacell_id,
      predicted_age_metacell_weeks = predicted_age_weeks,
      delta_predicted_age_weeks = delta_age_weeks,
      delta_definition = "predicted_age_metacell - observed/sample_age_weeks; pooled legacy no-strata display",
      plot_scope = "legacy_pooled_all_strata_no_matched_control_inference",
      celltype_label = paste0(
        harmonized_cell_class,
        if_else(projection_interpretation_flag == "diagnostic_only_projection_context", "*", ""),
        "\nMAE=", sprintf("%0.1f", cv_MAE_weeks),
        "w; R2=", sprintf("%0.2f", cv_R2),
        "; slope=", sprintf("%0.2f", cv_slope)
      ),
      harmonized_cell_class = factor(harmonized_cell_class, levels = stress_celltype_levels),
      celltype_label = factor(celltype_label, levels = rev(unique(celltype_label[order(harmonized_cell_class)])))
    ) %>%
    group_by(harmonized_cell_class, projection_group) %>%
    mutate(n_metacells_projection_group_celltype_total = n()) %>%
    ungroup() %>%
    arrange(harmonized_cell_class, projection_group, projection_metacell_id)
}

legacy_no_strata_count_labels <- function(dat) {
  dat %>%
    group_by(celltype_label, projection_group) %>%
    summarise(
      n_metacells = n(),
      x_min = min(delta_predicted_age_weeks, na.rm = TRUE),
      x_max = max(delta_predicted_age_weeks, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(celltype_label) %>%
    mutate(
      label_x = max(x_max, na.rm = TRUE) + pmax((max(x_max, na.rm = TRUE) - min(x_min, na.rm = TRUE)) * 0.05, 0.8),
      n_label = paste0("n=", n_metacells)
    ) %>%
    ungroup()
}

plot_legacy_no_strata_metacell_distribution <- function(dat, figure_id, title, subtitle, source_path) {
  count_labels <- legacy_no_strata_count_labels(dat)
  x_limits <- range(c(dat$delta_predicted_age_weeks, count_labels$label_x), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04
  p <- ggplot(dat, aes(delta_predicted_age_weeks, celltype_label, fill = projection_group)) +
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
      size = 2.05,
      hjust = 0,
      show.legend = FALSE
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
      subtitle = subtitle
    )
  save_legacy_plot_with_data(
    p, dat, "05_stress_test_projection", figure_id, source_path,
    width = 9.0, height = 6.6
  )
}

plot_legacy_no_strata_heatmap <- function(summary_tbl, figure_id, title, source_path) {
  dat <- summary_tbl %>%
    ensure_model_quality() %>%
    mutate(
      projection_group = case_when(
        intervention_class == "vcd_model" ~ "vcd",
        intervention_class == "foxl2_model" ~ "foxl2",
        TRUE ~ intervention_class
      ),
      projection_group = factor(projection_group, levels = c("vcd", "vcd_control", "foxl2", "foxl2_control")),
      harmonized_cell_class = factor(harmonized_cell_class, levels = stress_celltype_levels),
      label = sprintf("%0.1f\nn=%s", median_delta_age_weeks, n_projected_metacells),
      plot_scope = "legacy_pooled_all_strata_no_matched_control_inference"
    )
  p <- ggplot(dat, aes(harmonized_cell_class, projection_group, fill = median_delta_age_weeks)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label), size = 2.35, color = "grey15", lineheight = 0.9) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    theme_clock(base_size = 8.2) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.ticks.y = element_blank()
    ) +
    labs(
      x = "Cell type",
      y = "Projection group",
      fill = "Median predicted -\nobserved age (weeks)",
      title = title,
      subtitle = "Legacy pooled no-strata display; not used for matched stress-test inference."
    )
  save_legacy_plot_with_data(
    p, dat, "05_stress_test_projection", figure_id, source_path,
    width = 8.4, height = 3.6
  )
}

plot_legacy_vs_controls <- function(dat, figure_id, title, source_path) {
  plot_dat <- dat %>%
    mutate(
      group = if_else(grepl("_control$", as.character(projection_group)), "control", "case"),
      group = factor(group, levels = c("case", "control")),
      plot_scope = "legacy_pooled_all_strata_case_control_observed_age_delta"
    )
  count_labels <- plot_dat %>%
    group_by(celltype_label, group) %>%
    summarise(
      n_metacells = n(),
      x_min = min(delta_predicted_age_weeks, na.rm = TRUE),
      x_max = max(delta_predicted_age_weeks, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(celltype_label) %>%
    mutate(
      label_x = max(x_max, na.rm = TRUE) + pmax((max(x_max, na.rm = TRUE) - min(x_min, na.rm = TRUE)) * 0.05, 0.8),
      n_label = paste0("n=", n_metacells)
    ) %>%
    ungroup()
  x_limits <- range(c(plot_dat$delta_predicted_age_weeks, count_labels$label_x), na.rm = TRUE)
  x_pad <- diff(x_limits) * 0.04
  p <- ggplot(plot_dat, aes(delta_predicted_age_weeks, celltype_label, fill = group)) +
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
      aes(fill = group, color = group),
      position = position_jitterdodge(jitter.width = 0, jitter.height = 0.12, dodge.width = 0.68, seed = 37),
      shape = 21,
      size = 0.48,
      alpha = 0.50,
      stroke = 0.18
    ) +
    geom_text(
      data = count_labels,
      aes(x = label_x, y = celltype_label, label = n_label, color = group, group = group),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.68),
      size = 2.05,
      hjust = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c(case = "#F8766D", control = "#00BFC4"), drop = FALSE) +
    scale_color_manual(values = c(case = "#B63A32", control = "#00787B"), drop = FALSE) +
    scale_x_continuous(limits = c(x_limits[1] - x_pad, x_limits[2] + x_pad), expand = expansion(mult = c(0.02, 0.02))) +
    theme_clock(base_size = 8.2) +
    theme(legend.position = "right") +
    labs(
      x = "Predicted - observed age (weeks)",
      y = "Cell type and internal CV quality",
      fill = "Projection group",
      color = "Projection group",
      title = title,
      subtitle = "Legacy pooled no-strata case/control display; not used for matched stress-test inference."
    )
  save_legacy_plot_with_data(
    p, plot_dat, "05_stress_test_projection", figure_id, source_path,
    width = 9.0, height = 6.6
  )
}

legacy_vcd_no_strata_data <- legacy_no_strata_metacell_data("VCD")
legacy_foxl2_no_strata_data <- legacy_no_strata_metacell_data("Foxl2")
plot_legacy_no_strata_metacell_distribution(
  legacy_vcd_no_strata_data,
  "vcd_metacell_delta_distribution_by_celltype",
  "VCD metacell delta distribution by cell type",
  "Legacy no-strata display: each point is one projected metacell; delta is predicted age minus observed/sample age.",
  paths$projection_metacell
)
plot_legacy_no_strata_metacell_distribution(
  legacy_foxl2_no_strata_data,
  "foxl2_metacell_delta_distribution_by_celltype",
  "Foxl2 metacell delta distribution by cell type",
  "Legacy no-strata display: each point is one projected metacell; delta is predicted age minus observed/sample age.",
  paths$projection_metacell
)
plot_legacy_no_strata_heatmap(
  vcd_summary_old,
  "vcd_projection_delta_heatmap",
  "VCD projection delta heatmap",
  paths$vcd_summary_old
)
plot_legacy_no_strata_heatmap(
  foxl2_summary_old,
  "foxl2_projection_delta_heatmap",
  "Foxl2 projection delta heatmap",
  paths$foxl2_summary_old
)
plot_legacy_vs_controls(
  legacy_vcd_no_strata_data,
  "vcd_vs_natural_controls_by_celltype",
  "VCD versus controls by cell type",
  paths$projection_metacell
)
plot_legacy_vs_controls(
  legacy_foxl2_no_strata_data,
  "foxl2_vs_natural_controls_by_celltype",
  "Foxl2 versus controls by cell type",
  paths$projection_metacell
)

old_pooled_stress_figures <- tribble(
  ~figure_id, ~figure_title, ~source_table_path,
  "vcd_projection_delta_heatmap", "VCD projection delta heatmap - superseded pooled all-age plot", paths$vcd_summary_old,
  "vcd_ranked_delta_by_celltype", "VCD ranked projection delta by cell type - superseded pooled all-age plot", paths$vcd_summary_old,
  "vcd_ranked_delta_by_celltype_with_model_quality", "VCD ranked projection delta with internal clock quality - superseded pooled all-age plot", paths$vcd_summary_old,
  "vcd_ranked_delta_by_celltype_confidence_filtered", "VCD ranked projection delta, confidence-filtered - superseded pooled all-age plot", paths$vcd_summary_old,
  "vcd_metacell_delta_distribution_by_celltype", "VCD metacell delta distribution - superseded pooled all-age plot", paths$projection_metacell,
  "vcd_vs_natural_controls_by_celltype", "VCD versus natural controls - superseded invalid-control plot", paths$projection_metacell,
  "foxl2_projection_delta_heatmap", "Foxl2 projection delta heatmap - superseded pooled all-age plot", paths$foxl2_summary_old,
  "foxl2_ranked_delta_by_celltype", "Foxl2 ranked projection delta by cell type - superseded pooled all-age plot", paths$foxl2_summary_old,
  "foxl2_ranked_delta_by_celltype_with_model_quality", "Foxl2 ranked projection delta with internal clock quality - superseded pooled all-age plot", paths$foxl2_summary_old,
  "foxl2_ranked_delta_by_celltype_confidence_filtered", "Foxl2 ranked projection delta, confidence-filtered - superseded pooled all-age plot", paths$foxl2_summary_old,
  "foxl2_metacell_delta_distribution_by_celltype", "Foxl2 metacell delta distribution - superseded pooled all-age plot", paths$projection_metacell,
  "foxl2_vs_natural_controls_by_celltype", "Foxl2 versus natural controls - superseded invalid-control plot", paths$projection_metacell,
  "vcd_ranked_delta_by_celltype_by_contrast", "VCD ranked projection delta by cell type and matched contrast - superseded ranked-dot plot", paths$vcd_summary,
  "vcd_ranked_delta_by_celltype_by_contrast_confidence_filtered", "VCD ranked projection delta by matched contrast, confidence-filtered - superseded ranked-dot plot", paths$vcd_summary,
  "foxl2_ranked_delta_by_celltype_by_age_stratum", "Foxl2 ranked projection delta by cell type and age stratum - superseded ranked-dot plot", paths$foxl2_summary,
  "foxl2_ranked_delta_by_celltype_by_age_stratum_confidence_filtered", "Foxl2 ranked projection delta by age stratum, confidence-filtered - superseded ranked-dot plot", paths$foxl2_summary
)
for (i in seq_len(nrow(old_pooled_stress_figures))) {
  old_row <- old_pooled_stress_figures[i, ]
  old_status <- if (grepl("ranked_delta|lollipop", old_row$figure_id)) {
    "superseded_ranked_or_lollipop_plot"
  } else {
    "superseded_pooled_all_age_invalid"
  }
  add_manifest(
    old_row$figure_id,
    old_row$figure_title,
    "05_stress_test_projection",
    file.path(figure_root, "05_stress_test_projection", paste0(old_row$figure_id, ".png")),
    file.path(figure_root, "05_stress_test_projection", paste0(old_row$figure_id, ".pdf")),
    file.path(plot_data_root, "05_stress_test_projection", paste0(old_row$figure_id, "_plot_data.tsv")),
    old_row$source_table_path,
    status = old_status
  )
}

# 99 combined panels: a compact index figure for navigation.
combined_data <- cv_eval %>%
  select(harmonized_cell_class, n_metacells, mae_weeks, r2, slope) %>%
  arrange(mae_weeks)
p_combined <- ggplot(combined_data, aes(mae_weeks, reorder(harmonized_cell_class, mae_weeks))) +
  geom_point(aes(size = n_metacells, color = r2), alpha = 0.9) +
  scale_color_gradient(low = "#BDBDBD", high = "#2166AC", limits = c(0, 1)) +
  theme_clock() +
  labs(x = "OOF metacell MAE (weeks)", y = "Cell type", color = "R2", size = "Metacells",
       title = "Combined metacell clock performance overview")
save_plot_with_data(
  p_combined, combined_data, "99_combined_panels", "combined_metacell_clock_performance_overview",
  "Combined metacell clock performance overview", paths$cv_eval, width = 7.2, height = 4.8
)

manifest <- manifest %>% arrange(section, figure_id)
fwrite(manifest, file.path(release_dir, "figure_manifest.tsv"), sep = "\t")

active_pooled_stress <- manifest %>%
  filter(figure_id %in% old_pooled_stress_figures$figure_id, status == "complete")
active_ranked_or_lollipop_stress <- manifest %>%
  filter(
    section == "05_stress_test_projection",
    status == "complete",
    grepl("ranked_delta|lollipop", figure_id)
  )
superseded_stress <- manifest %>%
  filter(
    figure_id %in% old_pooled_stress_figures$figure_id,
    grepl("^superseded", status)
  )

contrast_samples <- bind_rows(vcd_delta_with_quality, foxl2_delta_with_quality) %>%
  group_by(stress_test_family, matched_contrast_id) %>%
  summarise(
    case_samples = paste(sort(unique(case_sample_id)), collapse = "; "),
    control_samples = paste(sort(unique(unlist(strsplit(paste(control_sample_ids, collapse = ";"), ";")))), collapse = "; "),
    n_case_samples = n_distinct(case_sample_id),
    n_control_samples = first(n_control_samples),
    .groups = "drop"
  )

contrast_celltypes <- bind_rows(vcd_delta_with_quality, foxl2_delta_with_quality) %>%
  group_by(stress_test_family, matched_contrast_id) %>%
  summarise(
    n_cell_types = n_distinct(harmonized_cell_class),
    cell_types = paste(sort(unique(harmonized_cell_class)), collapse = "; "),
    .groups = "drop"
  )

new_distribution_ids <- manifest %>%
  filter(
    section == "05_stress_test_projection",
    status == "complete",
    grepl("^(vcd|foxl2)_delta_distribution_", figure_id)
  )

hard_checks <- tibble(
  check = c(
    "VCD has exactly 4 matched contrasts",
    "Foxl2 has exactly 2 matched age strata",
    "No VCD delta is computed across unmatched injection_age/post_injection_days",
    "No Foxl2 delta is computed across unmatched age_group",
    "Natural ageing reference samples are not used as universal controls for VCD/Foxl2",
    "Every stress-test delta has a matched_contrast_id",
    "All new stress-test delta-distribution figures are active and use metacell-level projection-group plot data",
    "Ranked/lollipop stress-test plots are marked superseded or removed from the active manifest",
    "Pooled all-age stress-test plots are marked superseded or removed from the active manifest"
  ),
  passed = c(
    n_distinct(vcd_delta_with_quality$matched_contrast_id) == 4L,
    n_distinct(foxl2_delta_with_quality$matched_contrast_id) == 2L,
    all(!is.na(vcd_delta_with_quality$injection_age) & !is.na(vcd_delta_with_quality$post_injection_days)) &&
      all(grepl("^VCD_injection_.+_post_.+", vcd_delta_with_quality$matched_contrast_id)),
    all(!is.na(foxl2_delta_with_quality$age_group)) &&
      all(grepl("^Foxl2_age_.+", foxl2_delta_with_quality$matched_contrast_id)),
    !any(grepl("(^|;)(AC_|YF_|OF_)", bind_rows(vcd_delta_with_quality, foxl2_delta_with_quality)$control_sample_ids)),
    all(!is.na(bind_rows(vcd_delta_with_quality, foxl2_delta_with_quality)$matched_contrast_id)),
    nrow(new_distribution_ids) == 18L &&
      all(grepl("case_metacell_id", vapply(
        file.path(release_dir, new_distribution_ids$plot_data_path),
        function(path) paste(names(fread(path, nrows = 1)), collapse = ";"),
        character(1)
      ))) &&
      all(grepl("projection_group", vapply(
        file.path(release_dir, new_distribution_ids$plot_data_path),
        function(path) paste(names(fread(path, nrows = 1)), collapse = ";"),
        character(1)
      ))),
    nrow(active_ranked_or_lollipop_stress) == 0L,
    nrow(active_pooled_stress) == 0
  )
)

md_table <- function(tbl) {
  tbl <- as.data.frame(tbl)
  header <- paste(names(tbl), collapse = " | ")
  divider <- paste(rep("---", ncol(tbl)), collapse = " | ")
  rows <- apply(tbl, 1, function(row) paste(row, collapse = " | "))
  c(paste0("| ", header, " |"), paste0("| ", divider, " |"), paste0("| ", rows, " |"))
}

audit_lines <- c(
  "# PRJNA863443 Stress-Test Contrast Audit",
  "",
  "PRJNA863443 stress-test projections are interpreted within matched design strata. VCD contrasts are matched by injection_age and post_injection_days. Foxl2 contrasts are matched by age group. Active stress-test delta-distribution figures are metacell-level: each dot is one projected metacell from the case or matched-control projection group, compared with the mean predicted age of matched control metacells from the same matched contrast and cell type. Pooled all-age VCD/Foxl2 deltas are invalid and have been superseded.",
  "",
  "## Matched Contrast Counts",
  "",
  paste0("- VCD matched contrasts: ", n_distinct(vcd_delta_with_quality$matched_contrast_id)),
  paste0("- Foxl2 matched age strata: ", n_distinct(foxl2_delta_with_quality$matched_contrast_id)),
  "",
  "## Case/Control Samples Per Contrast",
  "",
  md_table(contrast_samples),
  "",
  "## Cell Types Available Per Contrast",
  "",
  md_table(contrast_celltypes),
  "",
  "## Pooled Plot Status",
  "",
  paste0("- Any pooled all-age stress-test plot remains active: ", if_else(nrow(active_pooled_stress) > 0, "YES", "NO")),
  paste0("- Any ranked/lollipop stress-test plot remains active: ", if_else(nrow(active_ranked_or_lollipop_stress) > 0, "YES", "NO")),
  "- Superseded stress-test plots:",
  paste0("  - `", superseded_stress$figure_id, "` (", superseded_stress$status, ")"),
  "",
  "## Hard Checks",
  "",
  md_table(hard_checks %>% mutate(status = if_else(passed, "PASS", "FAIL")) %>% select(status, check))
)
writeLines(audit_lines, file.path(release_dir, "STRESS_TEST_CONTRAST_AUDIT.md"))

overview_lines <- c(
  "# Figure Package Overview",
  "",
  "Canonical run: scRNA_metacell_elasticnet_metacell_cv_20260707",
  "",
  "Primary framing: this package summarizes metacell-level internal out-of-fold ElasticNet clock predictions. Held-out test means held-out CV fold metacells, not a biological-sample holdout.",
  "",
  "## Overview/Eligibility",
  "",
  "- `figures/01_overview_eligibility/metacell_availability_and_model_eligibility.*`",
  "- `figures/01_overview_eligibility/nonzero_selected_genes_by_celltype.*`",
  "",
  "## Internal Held-Out Performance",
  "",
  "- `figures/02_internal_cv_performance/internal_oof_metacell_predictions_facet.*`",
  "- `figures/02_internal_cv_performance/internal_oof_<celltype>_pred_vs_obs.*`",
  "- `figures/02_internal_cv_performance/celltype_clock_performance_summary.*`",
  "- `figures/02_internal_cv_performance/celltype_clock_r2_summary.*`",
  "- `figures/02_internal_cv_performance/celltype_clock_slope_summary.*`",
  "- `figures/02_internal_cv_performance/residual_bias_by_age_bin_and_celltype.*`",
  "- `figures/02_internal_cv_performance/residual_bias_by_dataset_and_celltype.*`",
  "",
  "## Feature Selection",
  "",
  "- `figures/03_feature_selection/selected_gene_direction_by_celltype.*`",
  "- `figures/03_feature_selection/top_elasticnet_gene_features_by_celltype.*`",
  "",
  "## GSEA",
  "",
  "- `figures/04_gsea/gsea_eligibility_and_display_term_yield.*`",
  if (gsea_available) c(
    "- `figures/04_gsea/elasticnet_coefficient_gsea_nes_panels_by_celltype.*`",
    "- Primary GSEA is ElasticNet coefficient-ranked and interprets clock genes; positive NES indicates older-predictive clock-gene programs and negative NES indicates younger-predictive clock-gene programs.",
    "- ORA tables are supplementary selected-gene enrichment and should not be mislabeled as NES/GSEA.",
    "- `tables/04_gsea/GSEA_METHOD.md` documents the GSEA method and interpretation limits."
  ) else c(
    "- `figures/04_gsea/gsea_nes_panels_by_celltype.*`",
    "- No canonical metacell-CV GSEA result table was found, so the NES panel is explicitly marked as not available in `figure_manifest.tsv`."
  ),
  "",
  "## Stress-Test Projection",
  "",
  "Stress-test inference is shown as metacell-level delta-age distributions within matched contrasts. Each dot is one projected metacell from the case or matched-control projection group, where delta_predicted_age_weeks = predicted_age_metacell - mean(predicted_age_matched_control_metacells_within_same_matched_contrast_and_celltype).",
  "Box labels report the projected metacell count (`n=`) for each projection group within the plotted cell-type clock and matched contrast/age stratum.",
  "PRJNA863443 stress-test projections are interpreted within matched design strata. VCD has 4 matched contrasts defined by injection_age and post_injection_days. Foxl2 has 2 matched age strata defined by age group.",
  "Pooled all-age stress-test plots and ranked/lollipop stress-test plots are no longer active; any remaining legacy files are marked superseded in `figure_manifest.tsv`.",
  "",
  "- `figures/05_stress_test_projection/stress_test_metacell_availability_summary.*`",
  "- `figures/05_stress_test_projection/stress_test_metacell_availability_by_sample.*`",
  "- `figures/05_stress_test_projection/vcd_delta_age_by_pair_and_celltype.*`",
  "- `figures/05_stress_test_projection/foxl2_delta_age_by_pair_and_celltype.*`",
  "- `figures/05_stress_test_projection/vcd_delta_heatmap_by_pair.*`",
  "- `figures/05_stress_test_projection/foxl2_delta_heatmap_by_age_stratum.*`",
  "- `figures/05_stress_test_projection/vcd_projection_delta_heatmap_by_contrast.*`",
  "- `figures/05_stress_test_projection/foxl2_projection_delta_heatmap_by_age_stratum.*`",
  "- `figures/05_stress_test_projection/vcd_delta_distribution_<celltype>_by_matched_contrast.*`",
  "- `figures/05_stress_test_projection/foxl2_delta_distribution_<celltype>_by_age_stratum.*`",
  "- `figures/05_stress_test_projection/vcd_delta_distribution_by_celltype_combined.*`",
  "- `figures/05_stress_test_projection/foxl2_delta_distribution_by_celltype_combined.*`",
  "- `plot_data/05_stress_test_projection/vcd_projection_delta_by_contrast_with_model_quality.tsv`",
  "- `plot_data/05_stress_test_projection/foxl2_projection_delta_by_age_stratum_with_model_quality.tsv`",
  "- `STRESS_TEST_CONTRAST_AUDIT.md`",
  "",
  "All figure plot-data TSVs are under matching `plot_data/` subfolders; copied source tables are under matching `tables/` subfolders."
)
writeLines(overview_lines, file.path(release_dir, "FIGURE_OVERVIEW.md"))

message("Figure package organized: ", release_dir)
message("Figures in manifest: ", nrow(manifest))
