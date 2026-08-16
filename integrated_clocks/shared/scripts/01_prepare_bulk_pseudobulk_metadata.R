#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)][1])
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
data_input_dir <- normalizePath(
  Sys.getenv("OVARY_CLOCK_DATA_ROOT", unset = file.path(root_dir, "data_input")),
  mustWork = FALSE
)
out_dir <- file.path(project_dir, "metadata")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bulk_metadata_path <- file.path(data_input_dir, "bulkseq_dataset", "metadata_C57BL6_bulk_seq.csv")
scrna_metadata_path <- file.path(
  data_input_dir,
  "scRNA_pseudobulk",
  "clocks",
  "sample_level_pseudobulk",
  "step_00_metadata",
  "sample_pseudobulk_metadata.tsv"
)
existing_combined_metadata_path <- file.path(out_dir, "bulk_pseudobulk_sample_metadata.csv")
gse235021_scrna_metadata_path <- file.path(
  data_input_dir,
  "aggregated_pseudobulk",
  "GSE235021",
  "GSE235021_sample_metadata.tsv"
)

stopifnot(file.exists(bulk_metadata_path))
if (!file.exists(scrna_metadata_path) && !file.exists(existing_combined_metadata_path)) {
  stop(
    "Neither the historical scRNA pseudobulk metadata nor an existing combined metadata file exists.",
    call. = FALSE
  )
}

adult_min_age_week <- 7
expression_filter_rule <- "gene retained when count >= 10 in at least 20% of samples within each dataset/modality"
benchmark_validation_datasets <- c(
  "GSE252791",
  "GSE279573",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary"
)
benchmark_test_datasets <- c(
  "CRA003645",
  "2026_CommBio_Galligos_ovarian_inflammaging"
)
old_age_stress_test_dataset <- "NIHMS1554107"
scrna_atlas_stress_test_dataset <- "2024_Cell_Winkler_FRT_scRNA_spatial"

# Matched controls retained only to anchor post-lock external projections.
# These samples are intentionally excluded from the 117-sample model-development
# cohort even though their biological annotations are control-like.
projection_reference_control_sample_keys <- c(
  paste0("bulk_rnaseq__GSE273724__GSE273724__young_7_8w__SRR301255", 58:61),
  paste0("bulk_rnaseq__GSE279446__GSE279446__wt_ovary_subset__WT_ovary", 1:4),
  paste0(
    "bulk_rnaseq__GSE279446__GSE279446__wt_ovary_subset__Zeb1_WT_ovary_",
    1:4,
    "_youngCohort"
  ),
  paste0(
    "bulk_rnaseq__GSE279573__GSE279573__wt_chow_fertile__SRR310133",
    c(39, 44)
  ),
  paste0("bulk_rnaseq__GSE329024__Ova-K", 1:4),
  paste0("bulk_rnaseq__GSE329024__Ova-S", 1:4)
)

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

age_bin_from_week <- function(age_week) {
  fifelse(
    is.na(age_week), NA_character_,
    fifelse(
      age_week < adult_min_age_week, "developmental_or_prepubertal",
      fifelse(
        age_week < 16, "young_adult",
        fifelse(
          age_week < 40, "mature_adult",
          fifelse(age_week < 65, "advanced_maternal_age", "post_reproductive_ageing")
        )
      )
    )
  )
}

split_stratum_from_week <- function(age_week) {
  fifelse(
    is.na(age_week), NA_character_,
    fifelse(
      age_week < adult_min_age_week, "developmental_or_prepubertal",
      fifelse(
        age_week < 13, "early_adult",
        fifelse(age_week < 39, "reproductive_adult", "aged_adult")
      )
    )
  )
}

sc_compatible_age_bin <- function(age_week) {
  fifelse(
    is.na(age_week), NA_character_,
    fifelse(
      age_week < adult_min_age_week, "developmental_or_prepubertal",
      fifelse(
        age_week < 16, "young_adult",
        fifelse(age_week < 40, "mid_adult", fifelse(age_week < 65, "aged", "old_aged"))
      )
    )
  )
}

clean_strain <- function(strain) {
  x <- toupper(gsub("[^A-Za-z0-9]", "", na_chr(strain)))
  fifelse(
    is.na(x), NA_character_,
    fifelse(grepl("^C57", x) | x %chin% c("C57BL6", "C57BL6J", "C57BL"), "C57_compatible", "non_C57_or_unclear")
  )
}

classify_bulk_context <- function(genotype, treatment, sample_title, type) {
  txt <- tolower(paste(na_chr(genotype), na_chr(treatment), na_chr(sample_title), na_chr(type)))
  fifelse(
    grepl("external_projection_only", txt), "external_projection_only",
    fifelse(
      grepl("cums|stress", txt), "stress_model",
      fifelse(
        grepl("ko|knock|mutant|tr/tr|cko|heterozygous|zeb1|abcg5|bnc1|srsf1|rnf111", txt), "genetic_model",
        fifelse(
          grepl("nmn|hep14|high-fat|sucrose|hfhs|fmt|microbiota|diet|vehicle|water control|rescue", txt),
          "intervention_or_vehicle",
          fifelse(grepl("\\bwt\\b|control|natural aging|untreated", txt), "control_or_natural_aging", "unclassified")
        )
      )
    )
  )
}

dataset_display_name <- function(dataset_folder) {
  lookup <- c(
    "2022_eLife_Morris_cycling_murine_ovary_atlas" = "cycling_murine_ovary_atlas",
    "2024_Cell_Winkler_FRT_scRNA_spatial" = "female_reproductive_tract_aging_atlas",
    "2024_LifeMedicine_Jiang_whole_life_mouse_ovary" = "whole_life_mouse_ovary",
    "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "aging_mouse_ovary_atlas",
    "2025_CommBio_Zhao_developing_ovary_scRNA_spatial" = "developing_ovary_scRNA_spatial",
    "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary" = "multiomic_aging_mouse_ovary",
    "2026_CommBio_Galligos_ovarian_inflammaging" = "ovarian_inflammaging",
    "CRA003645" = "bulk_CRA003645_C57_ovary_timecourse",
    "CRA017492" = "bulk_CRA017492_NMN_aging_ovary",
    "GSE194196" = "bulk_GSE194196_Bnc1_model",
    "GSE228975" = "bulk_GSE228975_neonatal_Srsf1_model",
    "GSE240547" = "bulk_GSE240547_diet_ovary",
    "GSE252791" = "bulk_GSE252791_Rnf111_estrus",
    "GSE273724" = "bulk_GSE273724_HEP14_aging_ovary",
    "GSE279446" = "bulk_GSE279446_Zeb1_model",
    "GSE279573" = "bulk_GSE279573_Abcg5g8_model",
    "GSE329024" = "bulk_GSE329024_CUMS_ovary",
    "GSE235021" = "bulk_GSE235021_rice_diet_ovary",
    "GSE154890" = "bulk_GSE154890_natural_aging_ovary",
    "GSE249151" = "bulk_GSE249151_VCD_POI_ovary",
    "GSE106339" = "bulk_GSE106339_GOAT_model_ovary",
    "GSE128240" = "bulk_GSE128240_cyclophosphamide_ovary",
    "GSE98497" = "bulk_GSE98497_Mtor_CKO_ovary",
    "PRJNA787028" = "bulk_PRJNA787028_housing_density_ovary",
    "PRJNA924255" = "bulk_PRJNA924255_AAV_sh_Sema4C_ovary",
    "PRJNA1002222" = "bulk_PRJNA1002222_Cd38_KO_ovary",
    "NIHMS1554107" = "bulk_NIHMS1554107_natural_aging_timecourse",
    "PRJNA1076509" = "bulk_PRJNA1076509_FMT_ovary"
  )
  out <- unname(lookup[as.character(dataset_folder)])
  out[is.na(out)] <- as.character(dataset_folder)[is.na(out)]
  out
}

append_reason <- function(reason, flag, label) {
  flag[is.na(flag)] <- FALSE
  reason[flag & is.na(reason)] <- label
  add <- flag & !is.na(reason) & !grepl(label, reason, fixed = TRUE)
  reason[add] <- paste(reason[add], label, sep = ";")
  reason
}

bulk <- fread(bulk_metadata_path)
bulk[, modality := "bulk_rnaseq"]
bulk[, sample_id := na_chr(`Sample name`)]
bulk[, dataset_name := dataset_display_name(dataset_folder)]
bulk[, age_week := fifelse(`Chronological age unit` == "Week", as.numeric(`Chronological age`),
  fifelse(`Chronological age unit` == "Day", as.numeric(`Chronological age`) / 7, NA_real_)
)]
bulk[, age_month_equivalent := age_week / 4.34524]
bulk[, strain_canonical := clean_strain(fcoalesce(na_chr(Strain), na_chr(Breed)))]
bulk[, tissue_canonical := tolower(na_chr(Tissue))]
bulk[, sex_canonical := tolower(na_chr(Sex))]
bulk[, reproductive_age_bin := age_bin_from_week(age_week)]
bulk[, bio_age_bin := reproductive_age_bin]
bulk[, split_age_stratum := split_stratum_from_week(age_week)]
bulk[, scRNA_compatible_age_bin := sc_compatible_age_bin(age_week)]
bulk[, sample_context := classify_bulk_context(Genotype, Treatment, `Sample title`, Type)]
bulk[, adult_age_range_include := age_week >= adult_min_age_week]
bulk[, exclusion_reason := NA_character_]
bulk[, exclusion_reason := append_reason(exclusion_reason, is.na(age_week), "missing_age_week")]
bulk[, exclusion_reason := append_reason(exclusion_reason, !is.na(age_week) & age_week < adult_min_age_week, "younger_than_7_weeks_or_developmental")]
bulk[, exclusion_reason := append_reason(exclusion_reason, strain_canonical != "C57_compatible", "non_C57_or_unclear_strain")]
bulk[, exclusion_reason := append_reason(exclusion_reason, tissue_canonical != "ovary", "non_ovary_tissue")]
bulk[, exclusion_reason := append_reason(exclusion_reason, sex_canonical != "female", "non_female_or_unclear_sex")]
bulk[, age_range_clock_candidate := is.na(exclusion_reason)]
bulk[, control_like_clock_candidate := age_range_clock_candidate & sample_context == "control_or_natural_aging"]
bulk[, expression_counts_path := file.path(
  root_dir, "data_input", "bulkseq_dataset", dataset_folder, "input", paste0(dataset_folder, "_gene_counts_filtered.csv")
)]
bulk[, expression_logcpm_path := file.path(
  root_dir, "data_input", "bulkseq_dataset", dataset_folder, "input", paste0(dataset_folder, "_gene_logcpm_filtered.csv")
)]
bulk[, expression_annotation_path := file.path(
  root_dir, "data_input", "bulkseq_dataset", dataset_folder, "input", paste0(dataset_folder, "_gene_annotation.csv")
)]
bulk[, metadata_source_path := bulk_metadata_path]
bulk[, expression_filter_rule := expression_filter_rule]
bulk[, n_cells_used := NA_integer_]
bulk[, n_classes_used := NA_integer_]
bulk[, included_harmonized_cell_classes := NA_character_]
bulk[, curation_status_summary := NA_character_]

bulk_out <- bulk[, .(
  modality,
  sample_key = paste(modality, dataset_folder, sample_id, sep = "__"),
  sample_id,
  dataset_folder,
  dataset_name,
  dataset_id,
  accession = Accession,
  project_accession = `Project accession`,
  sample_title = `Sample title`,
  age_week,
  age_month_equivalent,
  reproductive_age_bin,
  bio_age_bin,
  split_age_stratum,
  scRNA_compatible_age_bin,
  strain = Strain,
  strain_canonical,
  sex = Sex,
  tissue = Tissue,
  genotype = Genotype,
  treatment = Treatment,
  sample_context,
  dev_stage = `Dev stage`,
  dataformat,
  age_raw = `Age raw`,
  age_curation_confidence = `Age curation confidence`,
  n_cells_used,
  n_classes_used,
  included_harmonized_cell_classes,
  curation_status_summary,
  age_range_clock_candidate,
  control_like_clock_candidate,
  exclusion_reason,
  expression_counts_path,
  expression_logcpm_path,
  expression_annotation_path,
  metadata_source_path,
  expression_filter_rule
)]

if (file.exists(scrna_metadata_path)) {
  scrna <- fread(scrna_metadata_path)
  scrna[, modality := "scrna_pseudobulk"]
  scrna[, dataset_name := dataset_display_name(dataset_folder)]
  scrna[, age_month_equivalent := age_week / 4.34524]
  scrna[, reproductive_age_bin := age_bin_from_week(age_week)]
  scrna[, bio_age_bin := reproductive_age_bin]
  scrna[, split_age_stratum := split_stratum_from_week(age_week)]
  scrna[, scRNA_compatible_age_bin := sc_compatible_age_bin(age_week)]
  scrna[, sample_context := "control_or_natural_aging"]
  scrna[, strain := "C57BL/6-compatible"]
  scrna[, strain_canonical := "C57_compatible"]
  scrna[, sex := "female"]
  scrna[, tissue := "ovary"]
  scrna[, genotype := "WT_or_unspecified_C57"]
  scrna[, treatment := reproductive_state]
  scrna[, dev_stage := reproductive_age_bin]
  scrna[, dataformat := "sample_gene_counts_and_logcpm"]
  scrna[, age_raw := as.character(age_week)]
  scrna[, age_curation_confidence := "inherited_from_scRNA_clock_metadata"]
  scrna[, age_range_clock_candidate := age_week >= adult_min_age_week]
  scrna[, control_like_clock_candidate := age_range_clock_candidate]
  scrna[, exclusion_reason := fifelse(age_range_clock_candidate, NA_character_, "younger_than_7_weeks_or_developmental")]
  scrna[, expression_counts_path := file.path(
    root_dir, "data_input", "aggregated_pseudobulk", dataset_folder, "input_genelevel", "sample_gene_counts_after_expression_filter.tsv"
  )]
  scrna[, expression_logcpm_path := file.path(
    root_dir, "data_input", "aggregated_pseudobulk", dataset_folder, "input_genelevel", "sample_gene_logcpm_after_expression_filter.tsv"
  )]
  scrna[, expression_annotation_path := file.path(
    root_dir, "data_input", "aggregated_pseudobulk", dataset_folder, "input_genelevel", "gene_aggregation_map.tsv"
  )]
  scrna[, metadata_source_path := scrna_metadata_path]
  scrna[, expression_filter_rule := expression_filter_rule]

  scrna_out <- scrna[, .(
    modality,
    sample_key,
    sample_id,
    dataset_folder,
    dataset_name,
    dataset_id = dataset_folder,
    accession = NA_character_,
    project_accession = NA_character_,
    sample_title = sample_id,
    age_week,
    age_month_equivalent,
    reproductive_age_bin,
    bio_age_bin,
    split_age_stratum,
    scRNA_compatible_age_bin,
    strain,
    strain_canonical,
    sex,
    tissue,
    genotype,
    treatment,
    sample_context,
    dev_stage,
    dataformat,
    age_raw,
    age_curation_confidence,
    n_cells_used,
    n_classes_used,
    included_harmonized_cell_classes,
    curation_status_summary,
    age_range_clock_candidate,
    control_like_clock_candidate,
    exclusion_reason,
    expression_counts_path,
    expression_logcpm_path,
    expression_annotation_path,
    metadata_source_path,
    expression_filter_rule
  )]
} else {
  message(
    "Historical scRNA pseudobulk metadata was removed; preserving scRNA rows from ",
    existing_combined_metadata_path
  )
  existing_combined <- fread(existing_combined_metadata_path)
  scrna_out <- existing_combined[modality == "scrna_pseudobulk"]
  missing_columns <- setdiff(names(bulk_out), names(scrna_out))
  if (length(missing_columns) > 0) {
    stop(
      "Existing combined metadata is missing scRNA columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  scrna_out <- scrna_out[, names(bulk_out), with = FALSE]
  scrna_out[, metadata_source_path := existing_combined_metadata_path]
}

if (file.exists(gse235021_scrna_metadata_path)) {
  gse235021_scrna <- fread(gse235021_scrna_metadata_path)
  gse235021_scrna[, `:=`(
    modality = "scrna_pseudobulk",
    sample_key = paste("scrna_pseudobulk", dataset_folder, sample_id, sep = "__"),
    dataset_name = "scrna_pseudobulk_GSE235021_rice_diet_ovary",
    dataset_id = dataset_folder,
    accession = sample_id,
    project_accession = bioproject,
    age_month_equivalent = age_month,
    reproductive_age_bin = age_bin_from_week(age_week),
    bio_age_bin = age_bin_from_week(age_week),
    split_age_stratum = split_stratum_from_week(age_week),
    scRNA_compatible_age_bin = sc_compatible_age_bin(age_week),
    strain_canonical = "C57_compatible",
    genotype = "WT",
    treatment = gsub("_", " ", diet, fixed = TRUE),
    sample_context = "external_projection_only",
    dev_stage = age_bin_from_week(age_week),
    dataformat = "scRNA-derived sample pseudobulk gene counts and logCPM",
    age_raw = paste0(age_month, " months"),
    age_curation_confidence = "repository_sample_metadata",
    n_cells_used = n_cells,
    n_classes_used = NA_integer_,
    included_harmonized_cell_classes = "all_Cell_Ranger_filtered_cells",
    curation_status_summary = "sample_level_raw_pseudobulk",
    age_range_clock_candidate = age_week >= adult_min_age_week,
    control_like_clock_candidate = FALSE,
    exclusion_reason = NA_character_,
    expression_counts_path = file.path(
      root_dir, "data_input", "aggregated_pseudobulk", dataset_folder,
      "input_genelevel", "sample_gene_counts_after_expression_filter.tsv"
    ),
    expression_logcpm_path = file.path(
      root_dir, "data_input", "aggregated_pseudobulk", dataset_folder,
      "input_genelevel", "sample_gene_logcpm_after_expression_filter.tsv"
    ),
    expression_annotation_path = file.path(
      root_dir, "data_input", "aggregated_pseudobulk", dataset_folder,
      "input_genelevel", "gene_aggregation_map.tsv"
    ),
    metadata_source_path = gse235021_scrna_metadata_path,
    expression_filter_rule = expression_filter_rule
  )]
  gse235021_scrna_out <- gse235021_scrna[, .(
    modality, sample_key, sample_id, dataset_folder, dataset_name, dataset_id,
    accession, project_accession, sample_title, age_week, age_month_equivalent,
    reproductive_age_bin, bio_age_bin, split_age_stratum,
    scRNA_compatible_age_bin, strain, strain_canonical, sex, tissue, genotype,
    treatment, sample_context, dev_stage, dataformat, age_raw,
    age_curation_confidence, n_cells_used, n_classes_used,
    included_harmonized_cell_classes, curation_status_summary,
    age_range_clock_candidate, control_like_clock_candidate, exclusion_reason,
    expression_counts_path, expression_logcpm_path, expression_annotation_path,
    metadata_source_path, expression_filter_rule
  )]
  scrna_out <- scrna_out[dataset_folder != "GSE235021"]
  scrna_out <- rbindlist(
    list(scrna_out, gse235021_scrna_out),
    use.names = TRUE,
    fill = FALSE
  )
}

combined <- rbindlist(list(scrna_out, bulk_out), use.names = TRUE)
setorder(combined, modality, dataset_folder, age_week, sample_id)

combined[, projection_reference_control :=
  sample_key %chin% projection_reference_control_sample_keys]
missing_projection_reference_keys <- setdiff(
  projection_reference_control_sample_keys,
  combined$sample_key
)
if (length(missing_projection_reference_keys) > 0L) {
  stop(
    "Projection-reference controls are missing from harmonized metadata: ",
    paste(missing_projection_reference_keys, collapse = ", "),
    call. = FALSE
  )
}
projection_reference_counts <- combined[
  projection_reference_control == TRUE,
  .N,
  by = dataset_folder
]
expected_projection_reference_counts <- data.table(
  dataset_folder = c("GSE273724", "GSE279446", "GSE279573", "GSE329024"),
  N = c(4L, 8L, 2L, 8L)
)
setorder(projection_reference_counts, dataset_folder)
setorder(expected_projection_reference_counts, dataset_folder)
if (!identical(projection_reference_counts, expected_projection_reference_counts)) {
  stop("Projection-reference control registry does not match expected 4/8/2/8 counts.", call. = FALSE)
}
combined[projection_reference_control == TRUE, `:=`(
  sample_context = "external_projection_reference_control",
  control_like_clock_candidate = FALSE
)]
combined[, exclusion_reason := append_reason(
  exclusion_reason,
  projection_reference_control,
  "external_projection_reference_control_not_model_building"
)]
combined[, model_building_eligible := control_like_clock_candidate]

combined[, split_role := fifelse(
  control_like_clock_candidate & dataset_folder %chin% benchmark_validation_datasets,
  "validation",
  fifelse(
    control_like_clock_candidate & dataset_folder %chin% benchmark_test_datasets,
    "test",
    fifelse(control_like_clock_candidate, "train", "not_primary_control_like")
  )
)]
combined[projection_reference_control == TRUE, split_role := "external_projection_only"]
combined[, primary_training_weight := fifelse(split_role == "train", 1, NA_real_)]
combined[, inverse_stratum_weight_sensitivity := NA_real_]
training_n <- combined[split_role == "train", .N]
training_n_strata <- combined[split_role == "train", uniqueN(split_age_stratum)]
combined[split_role == "train", inverse_stratum_weight_sensitivity := {
  training_n / (training_n_strata * .N)
}, by = split_age_stratum]
combined[, old_age_stress_test_role := fifelse(
  control_like_clock_candidate & dataset_folder == old_age_stress_test_dataset,
  "stress_test",
  fifelse(control_like_clock_candidate, "stress_train", "not_primary_control_like")
)]
combined[, scrna_atlas_stress_test_role := fifelse(
  control_like_clock_candidate & dataset_folder == scrna_atlas_stress_test_dataset,
  "stress_test",
  fifelse(control_like_clock_candidate, "stress_train", "not_primary_control_like")
)]
combined[projection_reference_control == TRUE, `:=`(
  old_age_stress_test_role = "external_projection_only",
  scrna_atlas_stress_test_role = "external_projection_only"
)]

combined[, counts_file_exists := file.exists(expression_counts_path)]
combined[, logcpm_file_exists := file.exists(expression_logcpm_path)]
combined[, annotation_file_exists := file.exists(expression_annotation_path)]

summary_by_dataset <- combined[, .(
  n_samples_total = .N,
  n_age_range_clock_candidates = sum(age_range_clock_candidate, na.rm = TRUE),
  n_control_like_clock_candidates = sum(control_like_clock_candidate, na.rm = TRUE),
  n_model_building_eligible = sum(model_building_eligible, na.rm = TRUE),
  n_projection_reference_controls = sum(projection_reference_control, na.rm = TRUE),
  min_age_week = suppressWarnings(min(age_week, na.rm = TRUE)),
  max_age_week = suppressWarnings(max(age_week, na.rm = TRUE)),
  n_counts_files_present = sum(counts_file_exists),
  n_logcpm_files_present = sum(logcpm_file_exists),
  sample_contexts = paste(sort(unique(na.omit(sample_context))), collapse = ";")
), by = .(modality, dataset_folder, dataset_name)]
setorder(summary_by_dataset, modality, dataset_folder)

summary_by_age_bin <- combined[age_range_clock_candidate == TRUE, .N, by = .(modality, reproductive_age_bin)]
setorder(summary_by_age_bin, modality, reproductive_age_bin)

age_bin_ranges <- data.table(
  reproductive_age_bin = c(
    "young_adult",
    "mature_adult",
    "advanced_maternal_age",
    "post_reproductive_ageing"
  ),
  age_range_week = c(
    "7 to <16",
    "16 to <40",
    "40 to <65",
    ">=65"
  ),
  age_range_month = c(
    "1.6 to <3.7",
    "3.7 to <9.2",
    "9.2 to <15.0",
    ">=15.0"
  ),
  age_range_note = c(
    "adult young-cycling/fertile window",
    "mature adult before advanced maternal-age bin",
    "advanced maternal-age to late reproductive-age window",
    "post-reproductive/estropausal-age window"
  )
)
summary_by_age_bin <- merge(summary_by_age_bin, age_bin_ranges, by = "reproductive_age_bin", all.x = TRUE)
setcolorder(summary_by_age_bin, c(
  "modality",
  "reproductive_age_bin",
  "age_range_week",
  "age_range_month",
  "age_range_note",
  "N"
))
setorder(summary_by_age_bin, modality, reproductive_age_bin)

split_balance_summary <- combined[
  control_like_clock_candidate == TRUE,
  .N,
  by = .(split_role, modality, split_age_stratum, bio_age_bin)
]
setorder(split_balance_summary, split_role, modality, split_age_stratum, bio_age_bin)

split_assignment <- combined[, .(
  sample_key,
  sample_id,
  modality,
  dataset_folder,
  dataset_name,
  age_week,
  bio_age_bin,
  split_age_stratum,
  control_like_clock_candidate,
  model_building_eligible,
  projection_reference_control,
  split_role,
  primary_training_weight,
  inverse_stratum_weight_sensitivity,
  old_age_stress_test_role,
  scrna_atlas_stress_test_role
)]
setorder(split_assignment, split_role, modality, dataset_folder, age_week, sample_id)

old_age_stress_balance_summary <- combined[
  control_like_clock_candidate == TRUE,
  .N,
  by = .(
    stress_role = old_age_stress_test_role,
    modality,
    split_age_stratum,
    bio_age_bin
  )
]
old_age_stress_balance_summary[, stress_test := "old_age_bulk_generalization"]

scrna_atlas_stress_balance_summary <- combined[
  control_like_clock_candidate == TRUE,
  .N,
  by = .(
    stress_role = scrna_atlas_stress_test_role,
    modality,
    split_age_stratum,
    bio_age_bin
  )
]
scrna_atlas_stress_balance_summary[, stress_test := "scrna_atlas_generalization"]

stress_test_balance_summary <- rbindlist(list(
  old_age_stress_balance_summary,
  scrna_atlas_stress_balance_summary
), use.names = TRUE)
setcolorder(stress_test_balance_summary, c(
  "stress_test",
  "stress_role",
  "modality",
  "split_age_stratum",
  "bio_age_bin",
  "N"
))
setorder(stress_test_balance_summary, stress_test, stress_role, modality, split_age_stratum, bio_age_bin)

fwrite(combined, file.path(out_dir, "bulk_pseudobulk_sample_metadata.csv"))
fwrite(summary_by_dataset, file.path(out_dir, "bulk_pseudobulk_dataset_summary.csv"))
fwrite(summary_by_age_bin, file.path(out_dir, "bulk_pseudobulk_age_bin_summary.csv"))
fwrite(split_assignment, file.path(out_dir, "bulk_pseudobulk_split_assignment.csv"))
fwrite(split_balance_summary, file.path(out_dir, "bulk_pseudobulk_split_balance_summary.csv"))
fwrite(stress_test_balance_summary, file.path(out_dir, "bulk_pseudobulk_stress_test_balance_summary.csv"))

cat("Wrote combined metadata:", file.path(out_dir, "bulk_pseudobulk_sample_metadata.csv"), "\n")
cat("Rows:", nrow(combined), "\n")
cat("Age-range clock candidates:", combined[age_range_clock_candidate == TRUE, .N], "\n")
cat("Control-like clock candidates:", combined[control_like_clock_candidate == TRUE, .N], "\n")
cat("Model-building eligible samples:", combined[model_building_eligible == TRUE, .N], "\n")
cat("Projection-only matched reference controls:", combined[projection_reference_control == TRUE, .N], "\n")
cat("Wrote split assignment:", file.path(out_dir, "bulk_pseudobulk_split_assignment.csv"), "\n")
cat("Wrote split balance summary:", file.path(out_dir, "bulk_pseudobulk_split_balance_summary.csv"), "\n")
cat("Wrote stress-test balance summary:", file.path(out_dir, "bulk_pseudobulk_stress_test_balance_summary.csv"), "\n")
