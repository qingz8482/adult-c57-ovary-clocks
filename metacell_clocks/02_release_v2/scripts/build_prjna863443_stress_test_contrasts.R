#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
})

options(stringsAsFactors = FALSE)

root_candidates <- c(
  normalizePath(getwd(), mustWork = FALSE),
  Sys.getenv("OVARY_CLOCK_ROOT"),
  Sys.getenv("OVARY_CLOCK_ROOT")
)
root_dir <- Sys.getenv("OVARY_CLOCK_ROOT")  # reproduce: fixed root
run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")  # reproduce: v2 labels/clocks
projection_dir <- file.path(run_root, "step_06_projection_stress_test")
metadata_path <- file.path(
  root_dir,
  "adult_c57_ovary_clocks",
  "data_input",
  "scRNA_dataset",
  "PRJNA863443",
  "metadata",
  "PRJNA863443_sample_metadata.tsv"
)

paths <- list(
  projection_sample = file.path(projection_dir, "prjna863443_projection_sample_celltype_summary.tsv"),
  cv_eval = file.path(run_root, "step_05_internal_cv_validation", "celltype_clock_cv_evaluation_summary.tsv"),
  delta_contrast = file.path(projection_dir, "prjna863443_projection_delta_by_contrast_celltype.tsv"),
  vcd_contrast = file.path(projection_dir, "vcd_projection_summary_by_contrast_celltype.tsv"),
  foxl2_contrast = file.path(projection_dir, "foxl2_projection_summary_by_contrast_celltype.tsv")
)

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required table: ", path, call. = FALSE)
  fread(path, na.strings = c("", "NA")) %>% as_tibble()
}

clean_token <- function(x) {
  x <- as.character(x)
  x <- if_else(is.na(x) | x == "", NA_character_, x)
  x <- str_replace_all(x, "[^A-Za-z0-9+/-]+", "_")
  str_replace_all(x, "^_|_$", "")
}

metadata <- read_tsv(metadata_path) %>%
  mutate(
    age_group = clean_token(age_label),
    projection_branch = case_when(
      str_detect(model_detail, "^natural_aging") ~ "natural_aging_reference",
      model == "VCD" | str_detect(model, "^Foxl2") ~ "stress_test_projection",
      TRUE ~ "other_projection"
    ),
    stress_test_family = case_when(
      model == "VCD" ~ "VCD",
      str_detect(model, "^Foxl2") ~ "Foxl2",
      str_detect(model_detail, "^natural_aging") ~ "natural_aging_reference",
      TRUE ~ NA_character_
    ),
    case_control_status = case_when(
      model == "VCD" & condition == "VCD" ~ "case",
      model == "VCD" & condition == "CTL" ~ "control",
      str_detect(model, "^Foxl2") & genotype == "Foxl2+/-" ~ "case",
      str_detect(model, "^Foxl2") & genotype == "Foxl2+/+" ~ "control",
      str_detect(model_detail, "^natural_aging") ~ "reference",
      TRUE ~ NA_character_
    ),
    treatment_or_genotype = case_when(
      model == "VCD" ~ condition,
      str_detect(model, "^Foxl2") ~ genotype,
      str_detect(model_detail, "^natural_aging") ~ "WT_natural_aging_reference",
      TRUE ~ treatment
    ),
    matched_contrast_id = case_when(
      model == "VCD" ~ paste0("VCD_injection_", clean_token(injection_age), "_post_", clean_token(post_injection_days)),
      str_detect(model, "^Foxl2") ~ paste0("Foxl2_age_", age_group),
      str_detect(model_detail, "^natural_aging") ~ paste0("natural_aging_reference_age_", age_group),
      TRUE ~ NA_character_
    ),
    matched_control_key = case_when(
      model == "VCD" ~ paste0("CTL_injection_", clean_token(injection_age), "_post_", clean_token(post_injection_days)),
      str_detect(model, "^Foxl2") ~ paste0("Foxl2_WT_age_", age_group),
      str_detect(model_detail, "^natural_aging") ~ paste0("natural_aging_reference_age_", age_group),
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    sample_id,
    projection_branch,
    stress_test_family,
    case_control_status,
    treatment_or_genotype,
    metadata_age_weeks = age_weeks,
    age_group,
    injection_age,
    post_injection_days,
    matched_contrast_id,
    matched_control_key
  )

projection_sample <- read_tsv(paths$projection_sample)
cv_eval <- read_tsv(paths$cv_eval)

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

projection_enriched <- projection_sample %>%
  left_join(metadata, by = "sample_id") %>%
  mutate(
    age_weeks = coalesce(metadata_age_weeks, age_weeks),
    age_weeks = as.numeric(age_weeks),
    post_injection_days = as.character(post_injection_days)
  )

missing_metadata <- projection_enriched %>%
  filter(is.na(stress_test_family)) %>%
  distinct(sample_id) %>%
  pull(sample_id)
if (length(missing_metadata) > 0) {
  stop("Projection samples missing PRJNA863443 metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}

matched_projection <- projection_enriched %>%
  filter(stress_test_family %in% c("VCD", "Foxl2"), case_control_status %in% c("case", "control"))

control_means <- matched_projection %>%
  filter(case_control_status == "control") %>%
  group_by(stress_test_family, matched_contrast_id, matched_control_key, harmonized_cell_class) %>%
  summarise(
    control_sample_ids = paste(sort(unique(sample_id)), collapse = ";"),
    matched_control_mean_predicted_age_weeks = mean(predicted_age_weeks, na.rm = TRUE),
    n_control_metacells = sum(n_projected_metacells, na.rm = TRUE),
    n_control_samples = n_distinct(sample_id),
    .groups = "drop"
  )

case_rows <- matched_projection %>%
  filter(case_control_status == "case") %>%
  transmute(
    stress_test_family,
    matched_contrast_id,
    harmonized_cell_class,
    case_sample_id = sample_id,
    case_control_status,
    treatment_or_genotype,
    age_weeks,
    age_group,
    injection_age,
    post_injection_days,
    predicted_age_weeks,
    n_case_metacells = n_projected_metacells,
    n_case_samples = 1L
  )

delta_by_contrast <- case_rows %>%
  left_join(control_means, by = c("stress_test_family", "matched_contrast_id", "harmonized_cell_class")) %>%
  left_join(clock_quality, by = "harmonized_cell_class") %>%
  mutate(
    delta_predicted_age_weeks = predicted_age_weeks - matched_control_mean_predicted_age_weeks
  ) %>%
  select(
    stress_test_family,
    matched_contrast_id,
    harmonized_cell_class,
    case_sample_id,
    control_sample_ids,
    case_control_status,
    treatment_or_genotype,
    age_weeks,
    age_group,
    injection_age,
    post_injection_days,
    predicted_age_weeks,
    matched_control_mean_predicted_age_weeks,
    delta_predicted_age_weeks,
    n_case_metacells,
    n_control_metacells,
    n_case_samples,
    n_control_samples,
    cv_MAE_weeks,
    cv_R2,
    cv_slope,
    n_metacells_internal_cv,
    model_status,
    projection_interpretation_flag
  ) %>%
  arrange(stress_test_family, matched_contrast_id, harmonized_cell_class, case_sample_id)

unmatched_cases <- delta_by_contrast %>%
  filter(is.na(control_sample_ids) | is.na(matched_control_mean_predicted_age_weeks))
if (nrow(unmatched_cases) > 0) {
  stop("Matched controls missing for one or more projected case rows.", call. = FALSE)
}

if (any(delta_by_contrast$stress_test_family == "VCD" & (is.na(delta_by_contrast$injection_age) | is.na(delta_by_contrast$post_injection_days)))) {
  stop("VCD matched contrasts must include injection_age and post_injection_days.", call. = FALSE)
}
if (any(delta_by_contrast$stress_test_family == "Foxl2" & is.na(delta_by_contrast$age_group))) {
  stop("Foxl2 matched contrasts must include age_group.", call. = FALSE)
}
if (any(is.na(delta_by_contrast$matched_contrast_id) | delta_by_contrast$matched_contrast_id == "")) {
  stop("Every stress-test delta must have matched_contrast_id.", call. = FALSE)
}
if (any(str_detect(delta_by_contrast$control_sample_ids, "(^|;)(AC_|YF_|OF_)"))) {
  stop("Natural ageing reference samples cannot be used as VCD/Foxl2 controls.", call. = FALSE)
}

fwrite(delta_by_contrast, paths$delta_contrast, sep = "\t", na = "")
fwrite(filter(delta_by_contrast, stress_test_family == "VCD"), paths$vcd_contrast, sep = "\t", na = "")
fwrite(filter(delta_by_contrast, stress_test_family == "Foxl2"), paths$foxl2_contrast, sep = "\t", na = "")

message("Wrote matched stress-test projection contrast tables:")
message(" - ", paths$delta_contrast)
message(" - ", paths$vcd_contrast)
message(" - ", paths$foxl2_contrast)
