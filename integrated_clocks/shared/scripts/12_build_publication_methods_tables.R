#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
data_input_dir <- file.path(root_dir, "data_input")
results_dir <- file.path(project_dir, "results")
metadata_dir <- file.path(project_dir, "metadata")
out_dir <- file.path(project_dir, "methods_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

read_optional <- function(path, sep = NULL) {
  if (!file.exists(path)) return(data.table())
  if (is.null(sep)) sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  fread(path, sep = sep, na.strings = c("", "NA", "NaN"))
}

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  trimws(x)
}

collapse_unique <- function(x) {
  x <- sort(unique(na_chr(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

add_missing_cols <- function(dt, cols) {
  for (col in setdiff(cols, names(dt))) {
    dt[, (col) := NA]
  }
  dt
}

write_sheet <- function(dt, sheet_name, description) {
  path <- file.path(out_dir, paste0(sheet_name, ".tsv"))
  fwrite(dt, path, sep = "\t", na = "")
  data.table(
    sheet = sheet_name,
    file = basename(path),
    n_rows = nrow(dt),
    n_columns = ncol(dt),
    description = description
  )
}

file_info_dt <- function(paths, source_class) {
  paths <- unique(na_chr(paths))
  paths <- paths[!is.na(paths)]
  if (length(paths) == 0) {
    return(data.table())
  }
  info <- file.info(paths)
  data.table(
    source_class = source_class,
    path = paths,
    file_exists = file.exists(paths),
    file_size_bytes = as.numeric(info$size),
    modified_time = as.character(info$mtime),
    file_name = basename(paths)
  )
}

matrix_shape <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(data.table(n_rows = NA_integer_, n_sample_columns = NA_integer_))
  }
  header <- tryCatch(fread(path, nrows = 0), error = function(e) NULL)
  one_col <- tryCatch(fread(path, select = 1, showProgress = FALSE), error = function(e) NULL)
  data.table(
    n_rows = if (is.null(one_col)) NA_integer_ else nrow(one_col),
    n_sample_columns = if (is.null(header)) NA_integer_ else max(0L, ncol(header) - 1L)
  )
}

metadata_path <- file.path(metadata_dir, "bulk_pseudobulk_sample_metadata.csv")
dataset_summary_path <- file.path(metadata_dir, "bulk_pseudobulk_dataset_summary.csv")
split_balance_path <- file.path(metadata_dir, "bulk_pseudobulk_split_balance_summary.csv")
stress_balance_path <- file.path(metadata_dir, "bulk_pseudobulk_stress_test_balance_summary.csv")
age_bin_summary_path <- file.path(metadata_dir, "bulk_pseudobulk_age_bin_summary.csv")
pred_path <- file.path(results_dir, "integrated_bulk_pseudobulk_all_clock_predictions.tsv")
alignment_path <- file.path(results_dir, "integrated_bulk_pseudobulk_sample_expression_alignment.tsv")
alignment_issue_path <- file.path(results_dir, "integrated_bulk_pseudobulk_alignment_issues.tsv")
common_gene_path <- file.path(results_dir, "integrated_bulk_pseudobulk_common_gene_summary.tsv")
external_pred_path <- file.path(
  results_dir, "external_projection", "simplified_classes",
  "external_projection_predictions_all_simplified.tsv"
)
tage_contrast_manifest_path <- file.path(
  results_dir, "external_projection", "tAge", "manifests", "tage_contrast_manifest.tsv"
)
tage_sample_manifest_path <- file.path(
  results_dir, "external_projection", "tAge", "manifests", "tage_sample_manifest.tsv"
)
bulk_source_path <- file.path(data_input_dir, "bulkseq_dataset", "metadata_C57BL6_bulk_seq.csv")
scrna_source_path <- file.path(data_input_dir, "metadata_scRNA.tsv")
source_intake_path <- file.path(data_input_dir, "standardized_seurat", "00_metadata", "00_source_intake.tsv")
standardized_input_path <- file.path(data_input_dir, "standardized_seurat", "00_metadata", "standardized_input_manifest.tsv")
pseudobulk_manifest_path <- file.path(data_input_dir, "aggregated_pseudobulk", "aggregated_sample_pseudobulk_manifest.tsv")
clean_pseudobulk_manifest_path <- file.path(data_input_dir, "aggregated_pseudobulk", "clean_sample_pseudobulk_clock_input_manifest.tsv")
prjna_meta_path <- file.path(data_input_dir, "aggregated_pseudobulk", "PRJNA863443", "PRJNA863443_sample_metadata.tsv")

meta <- fread(metadata_path)
meta[, age_week := as.numeric(age_week)]
meta[, age_month_equivalent := as.numeric(age_month_equivalent)]
meta[, age_range_clock_candidate := truthy(age_range_clock_candidate)]
meta[, control_like_clock_candidate := truthy(control_like_clock_candidate)]
meta[, model_building_eligible := truthy(model_building_eligible)]
meta[, projection_reference_control := truthy(projection_reference_control)]
dataset_summary <- read_optional(dataset_summary_path)
split_balance <- read_optional(split_balance_path)
stress_balance <- read_optional(stress_balance_path)
age_bin_summary <- read_optional(age_bin_summary_path)
pred <- read_optional(pred_path)
alignment <- read_optional(alignment_path)
alignment_issues <- read_optional(alignment_issue_path)
common_gene_summary <- read_optional(common_gene_path)
external_pred <- read_optional(external_pred_path)
tage_contrast_manifest <- read_optional(tage_contrast_manifest_path)
tage_sample_manifest <- read_optional(tage_sample_manifest_path)
bulk_source <- read_optional(bulk_source_path, sep = ",")
scrna_source <- read_optional(scrna_source_path)
source_intake <- read_optional(source_intake_path)
standardized_input <- read_optional(standardized_input_path)
pseudobulk_manifest <- read_optional(pseudobulk_manifest_path)
clean_pseudobulk_manifest <- read_optional(clean_pseudobulk_manifest_path)
prjna_meta <- read_optional(prjna_meta_path)

added7_datasets <- c(
  "GSE235021", "GSE154890", "GSE249151", "GSE106339",
  "GSE128240", "GSE98497", "PRJNA787028"
)
prjna1002222_dataset <- "PRJNA1002222"

benchmark_samples <- if (nrow(pred) > 0) {
  unique(pred[, .(
    sample_key, split_role, modality, dataset_folder, bio_age_bin,
    split_age_stratum, age_week
  )])
} else {
  data.table()
}
if (nrow(benchmark_samples) > 0) {
  benchmark_pred_counts <- pred[, .(
    n_prediction_rows = .N,
    model_families = collapse_unique(model_family),
    selectors = collapse_unique(selector),
    fit_stages = collapse_unique(fit_stage)
  ), by = sample_key]
  benchmark_samples <- merge(
    benchmark_samples,
    benchmark_pred_counts,
    by = "sample_key",
    all.x = TRUE
  )
}

benchmark_keys <- if (nrow(benchmark_samples) > 0) benchmark_samples$sample_key else character()
external_unique <- if (nrow(external_pred) > 0) {
  unique(external_pred[, .(
    sample_key, dataset_folder, modality, age_week, bio_age_bin,
    projection_role, projection_class, projection_reason,
    matched_reference_control_flag, matched_reference_control_source,
    used_in_deployment_training_flag, primary_external_projection_include,
    primary_effect_interpretation_include, matched_control_allowed,
    reference_control_group_id, reference_control_group_label,
    treatment, genotype, sample_context, sample_title, sample_id, accession,
    projection_group, comparison_id, genotype_class, cohort,
    control_group_label, natural_aging_plot_group
  )])
} else {
  data.table()
}
external_samples <- copy(external_unique)
if (nrow(external_samples) > 0) {
  external_samples[, age_week := as.numeric(age_week)]
  external_samples[, split_age_stratum := NA_character_]
  external_samples[, split_role := projection_role]
  external_model_counts <- if (nrow(external_pred) > 0) {
    external_pred[, .(
      n_external_model_rows = .N,
      projected_model_families = collapse_unique(model_family),
      projected_selectors = collapse_unique(selector),
      min_gene_coverage_fraction = safe_min(gene_coverage_fraction),
      max_gene_coverage_fraction = safe_max(gene_coverage_fraction)
    ), by = sample_key]
  } else {
    data.table()
  }
  external_samples <- merge(external_samples, external_model_counts, by = "sample_key", all.x = TRUE)
  external_samples[, external_projection_added7 :=
    dataset_folder %in% added7_datasets &
      !(dataset_folder == "GSE235021" & modality == "scrna_pseudobulk")]
  external_samples[, external_projection_gse235021_scrna_addition :=
    dataset_folder == "GSE235021" & modality == "scrna_pseudobulk"]
  external_samples[, external_projection_prjna1002222_addition :=
    dataset_folder == prjna1002222_dataset]
  external_samples[, frozen_v1_model_ready := sample_key %in% benchmark_keys]
  external_samples[, projection_benchmark_boundary := "post_hoc_projection_does_not_modify_frozen_v1_benchmark"]
  external_samples[, model_building_control_scope := fcase(
    dataset_folder %in% c(
      "GSE273724", "GSE279446", "GSE279573", "PRJNA863443",
      added7_datasets
    ),
      "external_projection_or_sensitivity_only_not_primary_model_building",
    frozen_v1_model_ready == TRUE,
      "sample_also_has_frozen_v1_benchmark_role_projection_does_not_change_it",
    default = "external_projection_or_sensitivity_only"
  )]
  external_samples[, relative_control_limitation := fcase(
    dataset_folder == "GSE273724",
      "Only four young_7_8w samples are available as same-study relative controls; age bins are unbalanced, so these rows are sensitivity/projection context rather than frozen benchmark controls.",
    dataset_folder == "GSE279446",
      "Heterozygous and homozygous Zeb1 genetic models are stress-test projections compared with same-cohort WT controls.",
    dataset_folder == "GSE279573",
      "WT chow fertile samples provide same-study projection context, but the dataset is not treated as an eligible primary model-building control set.",
    dataset_folder == "PRJNA863443",
      "Added as a post hoc pseudobulk external projection cohort; not part of the frozen v1 training/test plan.",
    dataset_folder == "GSE235021",
      "Rice-diet effects are evaluated within age using dataset-level polished-rice controls pooled across bulk and scRNA-derived pseudobulk; 6- and 9-month contrasts are not pooled, and mixed-modality effects are modality/batch-confounded by design.",
    dataset_folder == "GSE154890",
      "Natural-age placement is evaluated against the 3-month reference and is not interpreted as a perturbation effect.",
    dataset_folder == "GSE249151",
      "VCD day 15 is compared with sesame-oil vehicle; the day-30 comparison uses the available earlier vehicle group and is control-age-confounded.",
    dataset_folder == "GSE106339",
      "GOAT-knockout samples are genetic models and WT samples are reference controls within the analysis-domain sample set; all are prepubertal (3 weeks), and the contrast is interpreted descriptively rather than as a primary adult contrast.",
    dataset_folder == "GSE128240",
      "Cyclophosphamide-treated ovaries are compared with same-study controls.",
    dataset_folder == "GSE98497",
      "Mtor oocyte-conditional knockout is compared with WT; PMSG exposure is shared by both groups.",
    dataset_folder == "PRJNA787028",
      "Eight-mice-per-cage housing is compared with two-mice-per-cage housing as the same-study reference.",
    dataset_folder == prjna1002222_dataset,
      "Cd38-knockout ovaries are compared with same-age WT controls separately at 8 and 12 months; the 2-month WT group is a placement-only young reference and is not part of either knockout contrast.",
    default = NA_character_
  )]
}

# Backward-compatible fallback for older tAge manifests. In the current build,
# PRJNA1002222 and both matched contrasts are already present in tAge.
prjna1002222_contrast_rows <- if (nrow(external_samples) > 0 &&
    !prjna1002222_dataset %chin% tage_contrast_manifest$dataset) {
  unique(external_samples[
    dataset_folder == prjna1002222_dataset &
      projection_group %chin% c("WT_8M", "Cd38KO_8M", "WT_12M", "Cd38KO_12M"),
    .(
      sample_key,
      dataset = prjna1002222_dataset,
      tage_contrast_id = comparison_id,
      tage_contrast_name = fcase(
        comparison_id == "PRJNA1002222_Cd38KO_vs_WT_8M",
          "PRJNA1002222 Cd38 knockout versus WT at 8 months",
        comparison_id == "PRJNA1002222_Cd38KO_vs_WT_12M",
          "PRJNA1002222 Cd38 knockout versus WT at 12 months"
      ),
      tage_contrast_type = "genetic_model",
      tage_contrast_role = fifelse(grepl("^WT_", projection_group), "control", "experimental"),
      tage_reference_group_id = paste0(comparison_id, "__controls"),
      tage_placement_only = FALSE,
      tage_confounding_flag = NA_character_,
      control_matching_age = fifelse(grepl("_8M$", comparison_id), "8 months", "12 months"),
      control_matching_genotype = "same-age WT reference",
      control_matching_treatment = "none; diestrus collection",
      control_matching_timing = "same study and collection age",
      control_matching_duration = NA_character_,
      control_matching_batch = "same study",
      control_matching_status = "prespecified_and_present",
      inclusion_reason = fifelse(grepl("^WT_", projection_group), "declared_reference", "declared_target"),
      exclusion_reason = NA_character_,
      contrast_warning = NA_character_
    )
  ])
} else data.table()
if (nrow(prjna1002222_contrast_rows) > 0) {
  tage_contrast_manifest <- tage_contrast_manifest[dataset != prjna1002222_dataset]
  tage_contrast_manifest <- rbindlist(
    list(tage_contrast_manifest, prjna1002222_contrast_rows),
    use.names = TRUE,
    fill = TRUE
  )
}

tage_contrast_definitions <- if (nrow(tage_contrast_manifest) > 0) {
  tage_contrast_manifest[, .(
    tage_contrast_name = collapse_unique(tage_contrast_name),
    tage_contrast_type = collapse_unique(tage_contrast_type),
    tage_reference_group_id = collapse_unique(tage_reference_group_id),
    tage_placement_only = all(tage_placement_only == TRUE),
    tage_confounding_flag = collapse_unique(tage_confounding_flag),
    control_matching_age = collapse_unique(control_matching_age),
    control_matching_genotype = collapse_unique(control_matching_genotype),
    control_matching_treatment = collapse_unique(control_matching_treatment),
    control_matching_timing = collapse_unique(control_matching_timing),
    control_matching_duration = collapse_unique(control_matching_duration),
    control_matching_batch = collapse_unique(control_matching_batch),
    control_matching_status = collapse_unique(control_matching_status),
    inclusion_reason = collapse_unique(inclusion_reason),
    exclusion_reason = collapse_unique(exclusion_reason),
    contrast_warning = collapse_unique(contrast_warning),
    n_reference_or_control_samples = uniqueN(
      sample_key[tage_contrast_role %in% c("control", "placement_reference")]
    ),
    n_experimental_or_target_samples = uniqueN(
      sample_key[tage_contrast_role %in% c("experimental", "placement_target")]
    ),
    n_samples_total = uniqueN(sample_key),
    external_projection_added7 = dataset %in% added7_datasets,
    external_projection_prjna1002222_addition = dataset == prjna1002222_dataset
  ), by = .(dataset, tage_contrast_id)]
} else {
  data.table()
}
if (nrow(tage_contrast_definitions) > 0) {
  setorder(tage_contrast_definitions, dataset, tage_contrast_id)
}
if (nrow(tage_contrast_manifest) > 0) {
  setorder(tage_contrast_manifest, dataset, tage_contrast_id, tage_contrast_role, sample_key)
}

dataset_summary_publication <- copy(dataset_summary)
if (
  nrow(dataset_summary_publication) > 0 &&
  nrow(external_samples[dataset_folder == "PRJNA863443"]) > 0 &&
  !("PRJNA863443" %in% dataset_summary_publication$dataset_folder)
) {
  prjna_after <- clean_pseudobulk_manifest[
    dataset_folder == "PRJNA863443" & stage == "after"
  ]
  prjna_external <- external_samples[dataset_folder == "PRJNA863443"]
  prjna_n_samples <- uniqueN(prjna_external$sample_key)
  prjna_counts_present <- nrow(prjna_after) == 1L && file.exists(prjna_after$counts_tsv)
  prjna_logcpm_present <- nrow(prjna_after) == 1L && file.exists(prjna_after$logcpm_tsv)
  dataset_summary_publication <- rbindlist(list(
    dataset_summary_publication,
    data.table(
      modality = "scrna_pseudobulk",
      dataset_folder = "PRJNA863443",
      dataset_name = "PRJNA863443_VCD_Foxl2_natural_ageing_ovary_pseudobulk",
      n_samples_total = prjna_n_samples,
      n_age_range_clock_candidates = 0L,
      n_control_like_clock_candidates = 0L,
      min_age_week = safe_min(prjna_external$age_week),
      max_age_week = safe_max(prjna_external$age_week),
      n_counts_files_present = ifelse(prjna_counts_present, prjna_n_samples, 0L),
      n_logcpm_files_present = ifelse(prjna_logcpm_present, prjna_n_samples, 0L),
      sample_contexts = paste(
        c(
          "external_projection_only", "control_or_natural_aging",
          "genetic_model", "intervention_or_vehicle"
        ),
        collapse = ";"
      )
    )
  ), use.names = TRUE, fill = TRUE)
}
if (nrow(dataset_summary_publication) > 0) {
  setorder(dataset_summary_publication, modality, dataset_folder)
}

source_intake_pub <- if (nrow(source_intake) > 0) {
  copy(source_intake)
} else {
  data.table()
}
if (nrow(source_intake_pub) > 0) {
  setnames(source_intake_pub, "dataset_id", "dataset_folder", skip_absent = TRUE)
}

pseudobulk_gene_after <- if (nrow(clean_pseudobulk_manifest) > 0) {
  clean_pseudobulk_manifest[stage == "after", .(
    dataset_folder,
    pseudobulk_gene_rows_after_filter = n_gene_rows,
    pseudobulk_input_samples = n_samples,
    pseudobulk_counts_tsv = counts_tsv,
    pseudobulk_logcpm_tsv = logcpm_tsv
  )]
} else {
  data.table()
}

model_ready_by_dataset <- if (nrow(benchmark_samples) > 0) {
  benchmark_samples[, .(
    n_model_ready_samples = uniqueN(sample_key),
    model_ready_split_roles = collapse_unique(split_role),
    model_ready_min_age_week = safe_min(age_week),
    model_ready_max_age_week = safe_max(age_week)
  ), by = .(dataset_folder, modality)]
} else {
  data.table()
}

external_by_dataset <- if (nrow(external_samples) > 0) {
  external_samples[, .(
    n_external_projection_samples = uniqueN(sample_key),
    external_projection_classes = collapse_unique(projection_class),
    external_projection_groups = collapse_unique(projection_group),
    external_min_age_week = safe_min(age_week),
    external_max_age_week = safe_max(age_week)
  ), by = .(dataset_folder, modality)]
} else {
  data.table()
}

tage_by_dataset <- if (nrow(tage_contrast_definitions) > 0) {
  tage_contrast_definitions[, .(
    n_tage_contrasts = uniqueN(tage_contrast_id),
    tage_contrast_ids = collapse_unique(tage_contrast_id),
    tage_contrast_types = collapse_unique(tage_contrast_type),
    tage_confounding_flags = collapse_unique(tage_confounding_flag)
  ), by = .(dataset_folder = dataset)]
} else {
  data.table()
}

dataset_overview <- meta[, .(
  dataset_name = collapse_unique(dataset_name),
  dataset_id = collapse_unique(dataset_id),
  project_accession = collapse_unique(project_accession),
  sample_accessions = collapse_unique(accession),
  n_harmonized_metadata_rows = .N,
  n_age_range_clock_candidates = sum(age_range_clock_candidate, na.rm = TRUE),
  n_control_like_clock_candidates = sum(control_like_clock_candidate, na.rm = TRUE),
  n_model_building_eligible = sum(model_building_eligible, na.rm = TRUE),
  n_projection_reference_controls = sum(projection_reference_control, na.rm = TRUE),
  min_age_week = safe_min(age_week),
  max_age_week = safe_max(age_week),
  biological_age_bins = collapse_unique(bio_age_bin),
  split_age_strata = collapse_unique(split_age_stratum),
  split_roles = collapse_unique(split_role),
  sample_contexts = collapse_unique(sample_context),
  strain_labels = collapse_unique(strain),
  strain_canonical = collapse_unique(strain_canonical),
  sex_labels = collapse_unique(sex),
  tissue_labels = collapse_unique(tissue),
  expression_filter_rule = collapse_unique(expression_filter_rule),
  counts_file_exists_all = all(counts_file_exists == TRUE, na.rm = TRUE),
  logcpm_file_exists_all = all(logcpm_file_exists == TRUE, na.rm = TRUE),
  annotation_file_exists_all = all(annotation_file_exists == TRUE, na.rm = TRUE)
), by = .(modality, dataset_folder)]

dataset_overview <- merge(dataset_overview, model_ready_by_dataset, by = c("dataset_folder", "modality"), all.x = TRUE)
dataset_overview <- merge(dataset_overview, external_by_dataset, by = c("dataset_folder", "modality"), all = TRUE)
if (nrow(tage_by_dataset) > 0) {
  dataset_overview <- merge(dataset_overview, tage_by_dataset, by = "dataset_folder", all.x = TRUE)
}
dataset_overview[, external_projection_added7 := dataset_folder %in% added7_datasets]
if (nrow(source_intake_pub) > 0) {
  dataset_overview <- merge(
    dataset_overview,
    source_intake_pub[, .(
      dataset_folder,
      paper_short = paper_short %||% NA_character_,
      paper = paper %||% NA_character_,
      source_accession = accession %||% NA_character_,
      input_class = input_class %||% NA_character_,
      source_format = source_format %||% NA_character_,
      raw_count_available = raw_count_available %||% NA_character_,
      source_analysis_use = analysis_use %||% NA_character_,
      source_limitations = limitations %||% NA_character_
    )],
    by = "dataset_folder",
    all.x = TRUE
  )
}
if (nrow(pseudobulk_gene_after) > 0) {
  dataset_overview <- merge(dataset_overview, pseudobulk_gene_after, by = "dataset_folder", all.x = TRUE)
}
setorder(dataset_overview, modality, dataset_folder)

sample_columns <- c(
  "sample_key", "modality", "dataset_folder", "dataset_name", "dataset_id",
  "sample_id", "accession", "project_accession", "sample_title",
  "age_week", "age_month_equivalent", "age_raw", "age_curation_confidence",
  "bio_age_bin", "reproductive_age_bin", "split_age_stratum",
  "strain", "strain_canonical", "sex", "tissue", "genotype", "treatment",
  "sample_context", "dev_stage", "dataformat", "n_cells_used", "n_classes_used",
  "included_harmonized_cell_classes", "curation_status_summary",
  "age_range_clock_candidate", "control_like_clock_candidate",
  "model_building_eligible", "projection_reference_control", "exclusion_reason",
  "split_role", "primary_training_weight", "inverse_stratum_weight_sensitivity",
  "old_age_stress_test_role", "scrna_atlas_stress_test_role",
  "expression_counts_path", "expression_logcpm_path", "expression_annotation_path",
  "expression_filter_rule", "metadata_source_path",
  "counts_file_exists", "logcpm_file_exists", "annotation_file_exists"
	)
	sample_metadata_publication <- copy(meta[, ..sample_columns])
sample_metadata_publication[, sample_metadata_scope := "frozen_harmonized_metadata"]
sample_metadata_publication[, projection_benchmark_boundary := fifelse(
  sample_key %in% external_samples$sample_key,
  "projection_or_sensitivity_row_available_without_rewriting_frozen_v1_benchmark",
  NA_character_
)]
sample_metadata_publication[, model_building_control_scope := fifelse(
  projection_reference_control,
  "external_projection_reference_control_not_model_building",
  fifelse(
    sample_key %in% benchmark_keys,
    "frozen_v1_benchmark_model_ready",
    fifelse(
      sample_key %in% external_samples$sample_key,
      "harmonized_metadata_sample_with_external_projection_or_sensitivity_role",
      "harmonized_metadata_sample"
    )
  )
)]
sample_metadata_publication[, relative_control_limitation := NA_character_]
	sample_metadata_publication[, frozen_v1_model_ready := sample_key %in% benchmark_keys]
	sample_metadata_publication[, later_external_projection_or_sensitivity := sample_key %in% external_samples$sample_key]

if (nrow(external_samples) > 0) {
  external_only_samples <- copy(external_samples[
    dataset_folder == "PRJNA863443" & !(sample_key %in% sample_metadata_publication$sample_key)
  ])
  if (nrow(external_only_samples) > 0) {
    external_sample_rows <- data.table(
      sample_key = external_only_samples$sample_key,
      modality = external_only_samples$modality,
      dataset_folder = external_only_samples$dataset_folder,
      dataset_name = external_only_samples$dataset_folder,
      dataset_id = external_only_samples$dataset_folder,
      sample_id = external_only_samples$sample_id,
      accession = external_only_samples$accession,
      project_accession = fifelse(external_only_samples$dataset_folder == "PRJNA863443", "PRJNA863443", NA_character_),
      sample_title = external_only_samples$sample_title,
      age_week = external_only_samples$age_week,
      age_month_equivalent = external_only_samples$age_week / 4.345,
      age_raw = fifelse(is.na(external_only_samples$age_week), NA_character_, paste0(external_only_samples$age_week, " weeks")),
      age_curation_confidence = "external_projection_metadata",
      bio_age_bin = external_only_samples$bio_age_bin,
      reproductive_age_bin = NA_character_,
      split_age_stratum = external_only_samples$split_age_stratum,
      strain = NA_character_,
      strain_canonical = NA_character_,
      sex = NA_character_,
      tissue = NA_character_,
      genotype = external_only_samples$genotype,
      treatment = external_only_samples$treatment,
      sample_context = external_only_samples$sample_context,
      dev_stage = NA_character_,
      dataformat = fifelse(external_only_samples$modality == "scrna_pseudobulk", "sample-level pseudobulk RNA-seq", "bulk RNA-seq"),
      n_cells_used = NA_integer_,
      n_classes_used = NA_integer_,
      included_harmonized_cell_classes = NA_character_,
      curation_status_summary = external_only_samples$projection_reason,
      age_range_clock_candidate = FALSE,
      control_like_clock_candidate = FALSE,
      model_building_eligible = FALSE,
      projection_reference_control = FALSE,
      exclusion_reason = external_only_samples$model_building_control_scope,
      split_role = external_only_samples$split_role,
      primary_training_weight = NA_real_,
      inverse_stratum_weight_sensitivity = NA_real_,
      old_age_stress_test_role = NA_character_,
      scrna_atlas_stress_test_role = NA_character_,
      expression_counts_path = NA_character_,
      expression_logcpm_path = NA_character_,
      expression_annotation_path = NA_character_,
      expression_filter_rule = NA_character_,
      metadata_source_path = NA_character_,
      counts_file_exists = NA,
      logcpm_file_exists = NA,
      annotation_file_exists = NA,
      sample_metadata_scope = "external_projection_or_sensitivity_only",
      projection_benchmark_boundary = external_only_samples$projection_benchmark_boundary,
      model_building_control_scope = external_only_samples$model_building_control_scope,
      relative_control_limitation = external_only_samples$relative_control_limitation,
      frozen_v1_model_ready = FALSE,
      later_external_projection_or_sensitivity = TRUE
    )

    if (nrow(prjna_meta) > 0) {
      prjna_lookup <- prjna_meta[, .(
        sample_id,
        prjna_sample_title = sample_title,
        prjna_sex = sex,
        prjna_tissue = tissue,
        prjna_strain = strain,
        prjna_age_raw = age_label,
        prjna_condition = condition,
        prjna_model_detail = model_detail,
        prjna_qc_estimated_cells = qc_estimated_cells
      )]
      external_sample_rows <- merge(external_sample_rows, prjna_lookup, by = "sample_id", all.x = TRUE)
      external_sample_rows[dataset_folder == "PRJNA863443", `:=`(
        sample_title = fcoalesce(sample_title, prjna_sample_title),
        sex = fcoalesce(sex, prjna_sex),
        tissue = fcoalesce(tissue, prjna_tissue),
        strain = fcoalesce(strain, prjna_strain),
        strain_canonical = fcoalesce(strain_canonical, prjna_strain),
        age_raw = fcoalesce(age_raw, prjna_age_raw),
        n_cells_used = as.integer(prjna_qc_estimated_cells),
        curation_status_summary = paste(
          "external projection metadata;",
          fcoalesce(prjna_condition, ""),
          fcoalesce(prjna_model_detail, "")
        ),
        metadata_source_path = prjna_meta_path
      )]
      external_sample_rows[, c(
        "prjna_sample_title", "prjna_sex", "prjna_tissue", "prjna_strain",
        "prjna_age_raw", "prjna_condition", "prjna_model_detail",
        "prjna_qc_estimated_cells"
      ) := NULL]
    }

    if (nrow(clean_pseudobulk_manifest) > 0) {
      prjna_after <- clean_pseudobulk_manifest[dataset_folder == "PRJNA863443" & stage == "after"]
      if (nrow(prjna_after) == 1) {
        external_sample_rows[dataset_folder == "PRJNA863443", `:=`(
          expression_counts_path = prjna_after$counts_tsv,
          expression_logcpm_path = prjna_after$logcpm_tsv,
          expression_annotation_path = file.path(prjna_after$input_dir, "gene_aggregation_map.tsv"),
          expression_filter_rule = "gene retained when count >= 10 in at least 20% of samples within each dataset/modality"
        )]
      }
    }
    external_sample_rows[, `:=`(
      counts_file_exists = file.exists(expression_counts_path),
      logcpm_file_exists = file.exists(expression_logcpm_path),
      annotation_file_exists = file.exists(expression_annotation_path)
	    )]
	    external_sample_rows <- add_missing_cols(external_sample_rows, names(sample_metadata_publication))
	    sample_metadata_publication_cols <- names(sample_metadata_publication)
	    sample_metadata_publication <- rbindlist(
	      list(sample_metadata_publication, external_sample_rows[, ..sample_metadata_publication_cols]),
	      fill = TRUE
	    )
	  }
}
	setorder(sample_metadata_publication, modality, dataset_folder, age_week, sample_id)

benchmark_publication <- merge(
  benchmark_samples,
  sample_metadata_publication,
  by = c("sample_key", "split_role", "modality", "dataset_folder", "bio_age_bin", "split_age_stratum", "age_week"),
  all.x = TRUE
)
if (nrow(alignment) > 0) {
  benchmark_publication <- merge(
    benchmark_publication,
    alignment[, .(sample_key, expression_column, matched_by)],
    by = "sample_key",
    all.x = TRUE
  )
}
setorder(benchmark_publication, split_role, modality, dataset_folder, age_week, sample_key)

external_publication <- copy(external_samples)
if (nrow(external_publication) > 0) {
  setorder(external_publication, dataset_folder, projection_class, age_week, sample_key)
}

external_model_rows <- if (nrow(external_pred) > 0) {
  external_pred[, .(
    external_source = "combined_external_projection_outputs",
    sample_key, dataset_folder, modality, age_week, bio_age_bin,
    projection_class, projection_group, comparison_id, genotype_class,
    treatment, genotype, sample_context, model_family, selector,
    predicted_age_week, delta_age_week, gene_coverage_fraction,
    n_model_genes, n_observed_model_genes, projection_gene_coverage_flag,
    matched_control_method, matched_control_source, n_matched_controls,
    relative_control_reliability,
    primary_external_projection_include,
    primary_effect_interpretation_include
  )]
} else {
  data.table()
}
if (nrow(external_model_rows) > 0) {
  if (nrow(external_samples) > 0) {
    external_model_rows <- merge(
      external_model_rows,
      external_samples[, .(
        sample_key, external_projection_added7,
        external_projection_gse235021_scrna_addition,
        external_projection_prjna1002222_addition,
        frozen_v1_model_ready, projection_benchmark_boundary,
        model_building_control_scope, relative_control_limitation
      )],
      by = "sample_key",
      all.x = TRUE
    )
  }
  setorder(external_model_rows, dataset_folder, sample_key, model_family)
}

expression_paths <- unique(meta[, .(
  dataset_folder, modality,
  counts_path = expression_counts_path,
  logcpm_path = expression_logcpm_path,
  annotation_path = expression_annotation_path,
  expression_filter_rule
)])
if (nrow(clean_pseudobulk_manifest) > 0) {
  external_clean_pseudobulk_paths <- clean_pseudobulk_manifest[
    stage == "after" & !(dataset_folder %in% expression_paths$dataset_folder),
    .(
      dataset_folder,
      modality = "scrna_pseudobulk",
      counts_path = counts_tsv,
      logcpm_path = logcpm_tsv,
      annotation_path = file.path(input_dir, "gene_aggregation_map.tsv"),
      expression_filter_rule = "gene retained when count >= 10 in at least 20% of samples within each dataset/modality"
    )
  ]
  expression_paths <- rbindlist(list(expression_paths, external_clean_pseudobulk_paths), fill = TRUE)
}
expr_inventory <- rbindlist(lapply(seq_len(nrow(expression_paths)), function(i) {
  row <- expression_paths[i]
  counts_shape <- matrix_shape(row$counts_path)
  logcpm_shape <- matrix_shape(row$logcpm_path)
  data.table(
    dataset_folder = row$dataset_folder,
    modality = row$modality,
    expression_filter_rule = row$expression_filter_rule,
    counts_path = row$counts_path,
    counts_file_exists = file.exists(row$counts_path),
    counts_n_rows = counts_shape$n_rows,
    counts_n_sample_columns = counts_shape$n_sample_columns,
    logcpm_path = row$logcpm_path,
    logcpm_file_exists = file.exists(row$logcpm_path),
    logcpm_n_rows = logcpm_shape$n_rows,
    logcpm_n_sample_columns = logcpm_shape$n_sample_columns,
    annotation_path = row$annotation_path,
    annotation_file_exists = file.exists(row$annotation_path)
  )
}), fill = TRUE)
if (nrow(pseudobulk_manifest) > 0) {
  expr_inventory <- merge(
    expr_inventory,
    pseudobulk_manifest[, .(
      dataset_folder,
      pseudobulk_n_cells_used = n_cells_used,
      count_source_seurat_path
    )],
    by = "dataset_folder",
    all.x = TRUE
  )
}
if (nrow(clean_pseudobulk_manifest) > 0) {
  clean_pseudobulk_gene_summary <- dcast(
    clean_pseudobulk_manifest,
    dataset_folder + n_samples ~ stage,
    value.var = "n_gene_rows"
  )
  setnames(
    clean_pseudobulk_gene_summary,
    old = c("n_samples", "before", "after"),
    new = c(
      "pseudobulk_n_samples",
      "pseudobulk_n_genes_before_filter",
      "pseudobulk_n_genes_after_expression_filter"
    ),
    skip_absent = TRUE
  )
  expr_inventory <- merge(
    expr_inventory,
    clean_pseudobulk_gene_summary,
    by = "dataset_folder",
    all.x = TRUE
  )
}
setorder(expr_inventory, modality, dataset_folder)

age_split_definitions <- rbindlist(list(
  data.table(
    dictionary = "biological_age_bin",
    value = c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing"),
    lower_bound_week_inclusive = c(7, 16, 40, 65),
    upper_bound_week_exclusive = c(16, 40, 65, NA_real_),
    description = c(
      "Adult young-cycling/fertile window.",
      "Mature adult before advanced maternal-age bin.",
      "Advanced maternal-age to late reproductive-age window.",
      "Post-reproductive/estropausal-age window."
    )
  ),
  data.table(
    dictionary = "split_age_stratum",
    value = c("early_adult", "reproductive_adult", "aged_adult"),
    lower_bound_week_inclusive = c(7, 13, 39),
    upper_bound_week_exclusive = c(13, 39, NA_real_),
    description = c(
      "Young adult stratum used for split balance.",
      "Reproductive adult stratum used for split balance.",
      "Aged adult stratum used for split balance."
    )
  ),
  data.table(
    dictionary = "analysis_role",
    value = c(
      "age_range_clock_candidate", "control_like_clock_candidate",
      "model_building_eligible", "projection_reference_control",
      "frozen_v1_model_ready", "external_projection_or_sensitivity",
      "primary_effect_interpretation_include"
    ),
    lower_bound_week_inclusive = NA_real_,
    upper_bound_week_exclusive = NA_real_,
    description = c(
      "Adult C57-compatible female ovary sample with age >=7 weeks before perturbation-context filtering.",
      "Adult control or natural-ageing sample used as primary chronological-age clock candidate.",
      "Sample eligible for the prespecified model-development cohort; identical to control_like_clock_candidate after projection-only reference controls are excluded.",
      "Matched control retained only for post-lock external projection and explicitly excluded from model development.",
      "Sample present in the frozen top-level v1 model-ready benchmark prediction table.",
      "Sample projected or evaluated after deployment artifacts were frozen; not used to retroactively change the v1 benchmark.",
      "External projected sample for which matched-control relative interpretation is allowed."
    )
  )
), fill = TRUE)

analysis_use_counts <- rbindlist(list(
  meta[, .(
    source = "harmonized_metadata",
    n_samples = .N,
    n_age_range_clock_candidates = sum(age_range_clock_candidate, na.rm = TRUE),
    n_control_like_clock_candidates = sum(control_like_clock_candidate, na.rm = TRUE),
    n_model_building_eligible = sum(model_building_eligible, na.rm = TRUE),
    n_projection_reference_controls = sum(projection_reference_control, na.rm = TRUE),
    n_frozen_v1_model_ready = sum(sample_key %in% benchmark_keys),
    n_external_projection_or_sensitivity = sum(sample_key %in% external_samples$sample_key)
  ), by = .(modality, dataset_folder, bio_age_bin, split_age_stratum, split_role)],
  if (nrow(external_samples) > 0) {
	    external_samples[, .(
	      source = "external_projection_samples",
	      n_samples = uniqueN(sample_key),
	      n_age_range_clock_candidates = NA_integer_,
	      n_control_like_clock_candidates = NA_integer_,
	      n_model_building_eligible = NA_integer_,
	      n_projection_reference_controls = NA_integer_,
	      n_frozen_v1_model_ready = sum(sample_key %in% benchmark_keys),
	      n_external_projection_or_sensitivity = uniqueN(sample_key)
	    ), by = .(modality, dataset_folder, bio_age_bin, split_age_stratum, split_role)]
	  } else data.table()
	), fill = TRUE)
setorder(analysis_use_counts, source, modality, dataset_folder, bio_age_bin)

source_files <- c(
  metadata_path, dataset_summary_path, split_balance_path, stress_balance_path,
  age_bin_summary_path, pred_path, alignment_path, alignment_issue_path,
  common_gene_path, external_pred_path, bulk_source_path,
  tage_contrast_manifest_path, tage_sample_manifest_path,
  scrna_source_path, source_intake_path, standardized_input_path,
  pseudobulk_manifest_path, clean_pseudobulk_manifest_path, prjna_meta_path,
  list.files(file.path(data_input_dir, "bulkseq_dataset"), pattern = "metadata\\.csv$|gene_annotation_symbol_repair_audit\\.tsv$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(data_input_dir, "aggregated_pseudobulk"), pattern = "sample_pseudobulk_sample_metadata\\.tsv$|sample_gene_count_filter\\.tsv$|gene_aggregation_map\\.tsv$|README\\.md$", recursive = TRUE, full.names = TRUE),
  list.files(file.path(data_input_dir, "standardized_seurat", "00_metadata"), full.names = TRUE)
)
source_file_inventory <- file_info_dt(source_files, "metadata_or_manifest_or_expression_annotation")
source_file_inventory[, relative_path := sub(paste0("^", normalizePath(root_dir, mustWork = FALSE), "/?"), "", normalizePath(path, mustWork = FALSE))]
setcolorder(source_file_inventory, c("source_class", "relative_path", "path", "file_exists", "file_size_bytes", "modified_time", "file_name"))
setorder(source_file_inventory, relative_path)

trig_cd1_generated_rows <- sum(c(
  nrow(dataset_overview[dataset_folder == "Trig_CD1_mice"]),
  nrow(sample_metadata_publication[dataset_folder == "Trig_CD1_mice"]),
  nrow(external_publication[dataset_folder == "Trig_CD1_mice"]),
  nrow(external_model_rows[dataset_folder == "Trig_CD1_mice"]),
  nrow(expr_inventory[dataset_folder == "Trig_CD1_mice"]),
  nrow(analysis_use_counts[dataset_folder == "Trig_CD1_mice"]),
  nrow(source_file_inventory[grepl("Trig_CD1", relative_path, fixed = TRUE)])
))

variable_dictionary <- data.table(
  field = c(
    "sample_key", "dataset_folder", "modality", "age_week", "bio_age_bin",
    "split_age_stratum", "sample_context", "age_range_clock_candidate",
    "control_like_clock_candidate", "model_building_eligible",
    "projection_reference_control", "frozen_v1_model_ready",
    "later_external_projection_or_sensitivity", "split_role",
	    "projection_class", "projection_group", "comparison_id",
	    "matched_control_allowed", "primary_effect_interpretation_include",
	    "projection_benchmark_boundary", "model_building_control_scope",
	    "relative_control_limitation", "external_projection_added7",
      "external_projection_gse235021_scrna_addition",
      "external_projection_prjna1002222_addition",
      "tage_contrast_id", "tage_contrast_type", "tage_placement_only",
	    "expression_filter_rule", "counts_n_rows", "logcpm_n_rows",
	    "n_model_ready_samples", "n_external_projection_samples"
	  ),
  definition = c(
    "Unique harmonized sample identifier used across metadata, expression alignment, predictions, and projection outputs.",
    "Dataset/study folder identifier used by the project.",
    "Expression modality: bulk_rnaseq or scrna_pseudobulk.",
    "Chronological age harmonized to weeks.",
    "Biological reporting age bin used in the manuscript.",
    "Coarser age stratum used for balanced split reporting.",
    "Curated biological context, such as control_or_natural_aging, genetic_model, stress_model, or intervention_or_vehicle.",
    "Adult C57-compatible female ovary sample with age >=7 weeks and basic tissue/strain/sex inclusion.",
    "Age-range sample further restricted to control or natural ageing context for primary clock fitting.",
    "Whether the sample belongs to the prespecified 117-sample model-development cohort eligibility set.",
    "Whether the sample is one of 22 matched controls retained only as an external-projection reference and excluded from model development.",
    "Sample present in the frozen top-level v1 model-ready benchmark prediction table.",
    "Sample used in post hoc external projection or sensitivity analysis after deployment clocks were frozen.",
    "Frozen split assignment: train, validation, test, or not_primary_control_like.",
    "External-projection class used for grouping projected samples.",
	    "External-projection comparison group.",
	    "Matched-control or projection comparison identifier.",
	    "Whether matched-control relative-delta interpretation is allowed for the projection row.",
	    "Whether the external projection row should be interpreted as a primary perturbation/effect comparison.",
	    "Explicit boundary that projection/sensitivity rows do not retroactively modify the frozen v1 benchmark.",
	    "Scope note distinguishing projection/reference context from primary model-building eligibility.",
	    "Dataset-specific limitation note for relative-control interpretation.",
      "Whether a row belongs to one of the seven external cohorts added in the July 2026 projection expansion.",
      "Whether a row is one of the 17 GSE235021 scRNA-derived pseudobulk samples added after the seven-cohort bulk expansion.",
      "Whether a row belongs to the 21-sample PRJNA1002222 whole-ovary bulk external-projection addition.",
      "Unique identifier for a prespecified tAge biological contrast.",
      "Biological class of a tAge contrast, such as intervention, genetic_model, diet, stress_model, natural_age_placement, or developmental_OOD.",
      "Whether a tAge contrast is a chronological-age placement comparison rather than a primary perturbation contrast.",
	    "Expression filter applied before modeling or projection.",
    "Number of rows in the count matrix file, including retained genes/transcripts depending on input level.",
    "Number of rows in the logCPM matrix file, including retained genes/transcripts depending on input level.",
    "Number of samples from the dataset present in the frozen v1 model-ready benchmark matrix.",
    "Number of samples from the dataset present in external projection/sensitivity outputs."
  )
)

added7_expected_counts <- data.table(
  dataset_folder = added7_datasets,
  expected_n = c(18L, 19L, 15L, 8L, 6L, 12L, 12L)
)
added7_observed_counts <- external_samples[
  external_projection_added7 == TRUE,
  .(observed_n = uniqueN(sample_key)),
  by = dataset_folder
]
added7_count_audit <- merge(
  added7_expected_counts, added7_observed_counts,
  by = "dataset_folder", all.x = TRUE
)
added7_count_audit[is.na(observed_n), observed_n := 0L]

integrity_checks <- rbindlist(list(
  data.table(check = "harmonized_sample_metadata_rows", expected = "378", observed = as.character(nrow(meta)), status = ifelse(nrow(meta) == 378, "pass", "review")),
  data.table(check = "frozen_v1_model_ready_unique_samples", expected = "117", observed = as.character(uniqueN(benchmark_samples$sample_key)), status = ifelse(uniqueN(benchmark_samples$sample_key) == 117, "pass", "review")),
  data.table(check = "model_building_eligible_samples", expected = "117", observed = as.character(sum(meta$model_building_eligible, na.rm = TRUE)), status = ifelse(sum(meta$model_building_eligible, na.rm = TRUE) == 117, "pass", "review")),
  data.table(check = "projection_reference_controls_excluded_from_model_building", expected = "22 projection references; 0 model eligible", observed = paste0(sum(meta$projection_reference_control, na.rm = TRUE), " projection references; ", sum(meta$projection_reference_control & meta$model_building_eligible, na.rm = TRUE), " model eligible"), status = ifelse(sum(meta$projection_reference_control, na.rm = TRUE) == 22 && sum(meta$projection_reference_control & meta$model_building_eligible, na.rm = TRUE) == 0, "pass", "review")),
  data.table(check = "dataset_overview_covers_harmonized_datasets", expected = as.character(uniqueN(meta$dataset_folder)), observed = as.character(uniqueN(dataset_overview$dataset_folder[dataset_overview$dataset_folder %in% meta$dataset_folder])), status = ifelse(uniqueN(dataset_overview$dataset_folder[dataset_overview$dataset_folder %in% meta$dataset_folder]) == uniqueN(meta$dataset_folder), "pass", "review")),
  data.table(check = "sample_metadata_publication_rows", expected = "411", observed = as.character(nrow(sample_metadata_publication)), status = ifelse(nrow(sample_metadata_publication) == 411, "pass", "review")),
  data.table(check = "external_projection_samples_present", expected = "300", observed = as.character(if (nrow(external_samples) > 0) uniqueN(external_samples$sample_key) else 0L), status = ifelse(uniqueN(external_samples$sample_key) == 300, "pass", "review")),
  data.table(check = "external_projection_model_rows", expected = "600", observed = as.character(nrow(external_model_rows)), status = ifelse(nrow(external_model_rows) == 600, "pass", "review")),
  data.table(check = "external_projection_datasets", expected = "20", observed = as.character(uniqueN(external_samples$dataset_folder)), status = ifelse(uniqueN(external_samples$dataset_folder) == 20, "pass", "review")),
  data.table(check = "added7_external_projection_samples", expected = "90", observed = as.character(uniqueN(external_samples[external_projection_added7 == TRUE, sample_key])), status = ifelse(uniqueN(external_samples[external_projection_added7 == TRUE, sample_key]) == 90, "pass", "review")),
  data.table(check = "gse235021_scrna_pseudobulk_external_samples", expected = "17", observed = as.character(uniqueN(external_samples[external_projection_gse235021_scrna_addition == TRUE, sample_key])), status = ifelse(uniqueN(external_samples[external_projection_gse235021_scrna_addition == TRUE, sample_key]) == 17, "pass", "review")),
  data.table(check = "added7_dataset_sample_counts", expected = paste(added7_count_audit$expected_n, collapse = ";"), observed = paste(added7_count_audit$observed_n, collapse = ";"), status = ifelse(all(added7_count_audit$expected_n == added7_count_audit$observed_n), "pass", "review")),
  data.table(check = "gse235021_scrna_dataset_summary_row", expected = "1 row; 17 samples", observed = paste0(nrow(dataset_summary_publication[dataset_folder == "GSE235021" & modality == "scrna_pseudobulk"]), " row; ", dataset_summary_publication[dataset_folder == "GSE235021" & modality == "scrna_pseudobulk", n_samples_total], " samples"), status = ifelse(nrow(dataset_summary_publication[dataset_folder == "GSE235021" & modality == "scrna_pseudobulk"]) == 1L && dataset_summary_publication[dataset_folder == "GSE235021" & modality == "scrna_pseudobulk", n_samples_total] == 17L, "pass", "review")),
  data.table(check = "prjna1002222_external_projection_samples", expected = "21", observed = as.character(uniqueN(external_samples[external_projection_prjna1002222_addition == TRUE, sample_key])), status = ifelse(uniqueN(external_samples[external_projection_prjna1002222_addition == TRUE, sample_key]) == 21, "pass", "review")),
  data.table(check = "tage_contrast_definitions", expected = "41", observed = as.character(nrow(tage_contrast_definitions)), status = ifelse(nrow(tage_contrast_definitions) == 41, "pass", "review")),
  data.table(check = "tage_added7_contrasts", expected = "11", observed = as.character(nrow(tage_contrast_definitions[external_projection_added7 == TRUE])), status = ifelse(nrow(tage_contrast_definitions[external_projection_added7 == TRUE]) == 11, "pass", "review")),
  data.table(check = "prjna1002222_contrasts", expected = "2", observed = as.character(nrow(tage_contrast_definitions[external_projection_prjna1002222_addition == TRUE])), status = ifelse(nrow(tage_contrast_definitions[external_projection_prjna1002222_addition == TRUE]) == 2, "pass", "review")),
  data.table(check = "dataset_summary_prjna863443_row", expected = "1 row; 33 external-only samples", observed = paste0(nrow(dataset_summary_publication[dataset_folder == "PRJNA863443"]), " row; ", dataset_summary_publication[dataset_folder == "PRJNA863443", n_samples_total], " samples"), status = ifelse(nrow(dataset_summary_publication[dataset_folder == "PRJNA863443"]) == 1L && dataset_summary_publication[dataset_folder == "PRJNA863443", n_samples_total] == 33L, "pass", "review")),
  data.table(check = "prjna863443_sample_metadata_publication_rows", expected = "33", observed = as.character(nrow(sample_metadata_publication[dataset_folder == "PRJNA863443" & modality == "scrna_pseudobulk"])), status = ifelse(nrow(sample_metadata_publication[dataset_folder == "PRJNA863443" & modality == "scrna_pseudobulk"]) == 33, "pass", "review")),
  data.table(check = "expression_inventory_rows", expected = ">= harmonized dataset count", observed = as.character(nrow(expr_inventory)), status = ifelse(nrow(expr_inventory) >= uniqueN(meta$dataset_folder), "pass", "review")),
  data.table(check = "prjna863443_expression_inventory_present", expected = "PRJNA863443 scrna_pseudobulk row", observed = as.character(nrow(expr_inventory[dataset_folder == "PRJNA863443" & modality == "scrna_pseudobulk"])), status = ifelse(nrow(expr_inventory[dataset_folder == "PRJNA863443" & modality == "scrna_pseudobulk"]) == 1, "pass", "review")),
  data.table(check = "pseudobulk_clean_gene_rows_match_counts_matrix", expected = "all scrna_pseudobulk rows", observed = as.character(nrow(expr_inventory[modality == "scrna_pseudobulk" & counts_n_rows == pseudobulk_n_genes_after_expression_filter])), status = ifelse(nrow(expr_inventory[modality == "scrna_pseudobulk" & counts_n_rows == pseudobulk_n_genes_after_expression_filter]) == nrow(expr_inventory[modality == "scrna_pseudobulk"]), "pass", "review")),
  data.table(check = "trig_cd1_excluded_from_methods_tables", expected = "0", observed = as.character(trig_cd1_generated_rows), status = ifelse(trig_cd1_generated_rows == 0, "pass", "review"))
), fill = TRUE)

manifest <- rbindlist(list(
  write_sheet(dataset_overview, "01_dataset_overview", "Dataset-level publication metadata, source accessions, analysis roles, age ranges, and expression-resource summary."),
  write_sheet(sample_metadata_publication, "02_sample_metadata_publication", "One row per harmonized sample plus PRJNA863443 external-only pseudobulk sample, including the 21 PRJNA1002222 whole-ovary bulk samples, with age, strain, sex, tissue, treatment/context, inclusion flags, split roles, and expression paths."),
  write_sheet(benchmark_publication, "03_frozen_v1_model_ready_samples", "One row per sample present in the frozen v1 model-ready benchmark matrix, including split role and expression alignment."),
  write_sheet(external_publication, "04_external_projection_samples", "One row per sample used in post hoc external projection or sensitivity analysis, including matched-control interpretation fields."),
  write_sheet(external_model_rows, "05_external_projection_model_rows", "One row per projected sample and retained model family, including predicted age, delta age, coverage, and matched-control metadata."),
  write_sheet(expr_inventory, "06_expression_resource_inventory", "Dataset-level expression matrix paths, file presence, matrix dimensions, and pseudobulk cell/gene summaries."),
  write_sheet(age_split_definitions, "07_age_split_definitions", "Definitions for age bins, split strata, and analysis-role flags."),
  write_sheet(analysis_use_counts, "08_analysis_use_counts", "Counts by dataset, modality, age bin, split/projection role, and analysis-use flag."),
  write_sheet(source_file_inventory, "09_source_file_inventory", "Publication-relevant metadata, manifest, annotation, and README source files from bulk-pseudobulk_clock and data_input."),
  write_sheet(variable_dictionary, "10_variable_dictionary", "Definitions for core publication metadata fields."),
  write_sheet(integrity_checks, "11_methods_table_integrity_checks", "Basic internal consistency checks for the generated metadata package."),
  if (nrow(dataset_summary_publication) > 0) write_sheet(dataset_summary_publication, "12_original_dataset_summary", "Project dataset-modality summary including PRJNA1002222, the external-only PRJNA863443 cohort, and the 17-sample GSE235021 scRNA-derived pseudobulk addition.") else NULL,
  if (nrow(split_balance) > 0) write_sheet(split_balance, "13_original_split_balance", "Direct copy of the frozen split balance summary.") else NULL,
  if (nrow(stress_balance) > 0) write_sheet(stress_balance, "14_original_stress_test_balance", "Direct copy of the stress-test balance summary.") else NULL,
  if (nrow(age_bin_summary) > 0) write_sheet(age_bin_summary, "15_original_age_bin_summary", "Direct copy of the age-bin count summary.") else NULL,
  if (nrow(alignment_issues) > 0) write_sheet(alignment_issues, "16_frozen_v1_alignment_issues", "Alignment issues recorded for the frozen v1 benchmark snapshot.") else NULL,
  if (nrow(common_gene_summary) > 0) write_sheet(common_gene_summary, "17_common_gene_summary", "Common-gene counts for bulk, pseudobulk, and integrated model universes.") else NULL
  ,
  if (nrow(tage_contrast_definitions) > 0) write_sheet(tage_contrast_definitions, "18_tage_contrast_definitions", "One row per validated tAge contrast, including matching rules, limitations, and reference/experimental sample counts.") else NULL,
  if (nrow(tage_contrast_manifest) > 0) write_sheet(tage_contrast_manifest, "19_tage_sample_contrast_membership", "Sample-to-contrast membership for all 41 validated tAge and integrated-clock contrasts, including control matching and placement/OOD flags.") else NULL
), fill = TRUE)
fwrite(manifest, file.path(out_dir, "00_table_manifest.tsv"), sep = "\t", na = "")

readme <- c(
  "# Publication Methods Tables",
  "",
  "This directory contains publication-ready metadata sheets for the adult C57 ovary bulk plus pseudobulk clock project.",
  "",
  "The tables are designed for manuscript Methods, reviewer audit, and downstream meta-analysis. They separate:",
  "",
  "- harmonized sample metadata;",
  "- frozen v1 benchmark model-ready samples;",
  "- post hoc external projection and sensitivity-analysis samples;",
  "- 22 matched external-projection reference controls explicitly excluded from the 117-sample model-development cohort;",
  "- the seven added external cohorts and their prespecified tAge contrasts;",
  "- the later 17-sample GSE235021 scRNA-derived pseudobulk addition, kept distinct from the original seven-cohort bulk expansion;",
  "- the later 21-sample PRJNA1002222 whole-ovary bulk addition and its two same-age Cd38-knockout versus WT contrasts;",
  "- expression-resource paths and dimensions;",
  "- age-bin, split-role, and variable definitions;",
  "- source-file provenance.",
  "",
  "Important boundary: the frozen top-level `results/` benchmark is not retroactively modified by later annotation-corrected external projection or sensitivity analyses.",
  "",
  "Recommended entry points:",
  "",
  "1. `00_table_manifest.tsv`: list of all sheets and row/column counts.",
  "2. `12_original_dataset_summary.tsv`: source for Supplementary Table 1.",
  "3. `15_original_age_bin_summary.tsv`: source for Supplementary Table 2.",
  "4. `02_sample_metadata_publication.tsv`: source for Supplementary Tables 3a and 3b (117 frozen model-development samples and 300 external-projection samples, respectively).",
  "5. `03_frozen_v1_model_ready_samples.tsv`: exact frozen benchmark samples.",
  "6. `04_external_projection_samples.tsv` and `05_external_projection_model_rows.tsv`: projection/sensitivity cohorts and model-row metadata.",
  "7. `18_tage_contrast_definitions.tsv` and `19_tage_sample_contrast_membership.tsv`: the same 41 validated contrasts used for integrated-clock and tAge projection, with sample membership.",
  "8. `10_variable_dictionary.tsv`: core field definitions.",
  "",
  paste0("Generated: ", Sys.time()),
  paste0("Generator: ", normalizePath(script_file %||% file.path(script_dir, "12_build_publication_methods_tables.R"), mustWork = FALSE))
)
writeLines(readme, file.path(out_dir, "README.md"))

message("Wrote publication methods tables to: ", out_dir)
message("Manifest: ", file.path(out_dir, "00_table_manifest.tsv"))
