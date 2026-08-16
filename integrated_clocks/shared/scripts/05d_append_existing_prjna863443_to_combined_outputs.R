#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- sub("^--file=", "", script_arg[[1]])
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
projection_dir <- file.path(project_dir, "results", "external_projection")
dataset_dir <- file.path(projection_dir, "PRJNA863443")

na_strings <- c("", "NA", "NaN")
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}
ensure_cols <- function(dt, cols) {
  for (cc in setdiff(cols, names(dt))) dt[, (cc) := NA]
  dt
}
read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing required PRJNA863443 file: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = na_strings)
}
append_replace_dataset <- function(filename, new_dt, dataset = "PRJNA863443") {
  path <- file.path(projection_dir, filename)
  old <- read_tsv(path)
  if ("dataset_folder" %in% names(old)) {
    old <- old[is.na(dataset_folder) | dataset_folder != dataset]
  } else if ("sample_key" %in% names(old) && "sample_key" %in% names(new_dt)) {
    old <- old[!sample_key %chin% new_dt$sample_key]
  }
  old <- ensure_cols(old, names(new_dt))
  new_dt <- ensure_cols(copy(new_dt), names(old))
  out <- rbindlist(
    list(old[, names(old), with = FALSE], new_dt[, names(old), with = FALSE]),
    use.names = TRUE,
    fill = TRUE
  )
  fwrite(out, path, sep = "\t")
  invisible(out)
}

pred <- read_tsv(file.path(dataset_dir, "prjna863443_projection_predictions.tsv"))
relative <- read_tsv(file.path(dataset_dir, "prjna863443_relative_to_controls.tsv"))
meta <- read_tsv(file.path(dataset_dir, "prjna863443_projection_metadata.tsv"))
pred[, primary_external_projection_include := as_bool(primary_external_projection_include)]
pred[, matched_reference_control_flag := as_bool(matched_reference_control_flag)]
relative[, primary_external_projection_include := as_bool(primary_external_projection_include)]

if (uniqueN(pred$sample_key) != 33L || uniqueN(relative$sample_key) != 33L) {
  stop("Existing PRJNA863443 outputs do not contain exactly 33 samples.", call. = FALSE)
}

append_replace_dataset("external_projection_predictions_all.tsv", pred)
append_replace_dataset(
  "external_projection_predictions.tsv",
  pred[primary_external_projection_include == TRUE]
)
append_replace_dataset("external_projection_relative_to_controls_all.tsv", relative)
append_replace_dataset(
  "external_projection_relative_to_controls.tsv",
  relative[primary_external_projection_include == TRUE]
)

reference_controls <- pred[matched_reference_control_flag == TRUE]
append_replace_dataset("external_projection_reference_controls.tsv", reference_controls)

control_pool_candidates <- meta[, .(
  dataset_folder, sample_key, sample_id, accession, modality, age_week, bio_age_bin,
  treatment, sample_context, sample_title, projection_role, projection_class,
  matched_reference_control_flag, matched_reference_control_source,
  reference_control_group_id, reference_control_group_label, control_pool_source,
  control_candidate_status, control_candidate_note
)]
append_replace_dataset(
  "external_projection_control_pool_candidates.tsv",
  control_pool_candidates
)

control_matching_audit <- relative[, .(
  dataset_folder, sample_key, model_family, matched_control_method,
  matched_control_source, n_matched_controls, matched_control_sample_keys,
  relative_control_reliability, relative_control_warning
)]
append_replace_dataset(
  "external_projection_control_matching_audit.tsv",
  control_matching_audit
)

append_replace_dataset(
  "external_projection_group_summary.tsv",
  read_tsv(file.path(dataset_dir, "prjna863443_group_summary.tsv"))
)
append_replace_dataset(
  "external_projection_model_concordance.tsv",
  read_tsv(file.path(dataset_dir, "prjna863443_model_concordance.tsv"))
)
append_replace_dataset(
  "external_projection_gene_coverage.tsv",
  read_tsv(file.path(dataset_dir, "prjna863443_gene_coverage.tsv"))
)

qc_flags <- pred[, .(
  sample_key, dataset_folder, model_family, selector, projection_gene_coverage_flag,
  OOD_flag, projection_distance, projection_distance_percentile,
  primary_external_projection_include, matched_control_allowed
)]
append_replace_dataset("external_projection_qc_flags.tsv", qc_flags)

pca_scores <- read_tsv(file.path(dataset_dir, "prjna863443_pca_scores.tsv"))
if (!"dataset_folder" %in% names(pca_scores)) {
  pca_scores[, dataset_folder := "PRJNA863443"]
}
append_replace_dataset("pca_projection_scores.tsv", pca_scores)

pca_contrib <- read_tsv(file.path(dataset_dir, "prjna863443_pca_pc_contributions.tsv"))
if (!"dataset_folder" %in% names(pca_contrib)) {
  pca_contrib[, dataset_folder := "PRJNA863443"]
}
append_replace_dataset("pca_projection_pc_contributions.tsv", pca_contrib)

availability_path <- file.path(
  projection_dir,
  "external_projection_dataset_control_availability.tsv"
)
availability <- read_tsv(availability_path)
availability <- availability[
  is.na(dataset_folder) | dataset_folder != "PRJNA863443"
]
availability_row <- data.table(
  dataset_folder = "PRJNA863443",
  n_projected_samples_all = 33L,
  n_projected_samples_primary = 33L,
  n_control_meta_candidates = 0L,
  n_projection_internal_control_candidates = uniqueN(
    pred[matched_reference_control_flag == TRUE]$sample_key
  ),
  has_same_dataset_control_meta = FALSE,
  has_projection_internal_control = TRUE,
  recommended_relative_method = paste(
    "PRJNA863443_VCD_same_time_CTL",
    "PRJNA863443_Foxl2_same_age_WT",
    "natural_aging_no_matched_control",
    sep = ";"
  ),
  expected_relative_method = paste(
    "VCD matched to CTL by injection age and time",
    "Foxl2 heterozygous matched to WT by age group",
    "natural ageing has no matched-control effect",
    sep = ";"
  ),
  observed_relative_method = paste(
    sort(unique(relative$matched_control_method)),
    collapse = ";"
  ),
  primary_external_projection_include = TRUE,
  dataset_note = paste(
    "Dedicated PRJNA863443 outputs were reused unchanged;",
    "their 33 samples were appended to all combined outputs."
  ),
  pass_fail_status = "PASS"
)
availability <- ensure_cols(availability, names(availability_row))
availability_row <- ensure_cols(availability_row, names(availability))
availability <- rbindlist(list(
  availability[, names(availability), with = FALSE],
  availability_row[, names(availability), with = FALSE]
), use.names = TRUE, fill = TRUE)
fwrite(availability, availability_path, sep = "\t")

pred_all <- read_tsv(file.path(
  projection_dir,
  "external_projection_predictions_all.tsv"
))
rel_all <- read_tsv(file.path(
  projection_dir,
  "external_projection_relative_to_controls_all.tsv"
))
pred_all[, primary_external_projection_include := as_bool(
  primary_external_projection_include
)]
invalid_relative <- rel_all[
  (matched_control_method == "no_matched_control" &
     is.finite(as.numeric(relative_delta_age_week))) |
  (as.numeric(n_matched_controls) == 0 &
     is.finite(as.numeric(relative_delta_age_week))) |
  (as_bool(matched_control_allowed) == FALSE &
     is.finite(as.numeric(relative_delta_age_week))) |
  (as_bool(primary_external_projection_include) == FALSE &
     is.finite(as.numeric(relative_delta_age_week)))
]
if (nrow(invalid_relative) > 0L) {
  stop("Appending PRJNA863443 introduced invalid relative-control rows.", call. = FALSE)
}

manifest <- data.table(
  field = c(
    "run_label", "output_dir", "prjna863443_appended_by",
    "n_projected_samples_all", "n_projected_samples_primary",
    "n_prediction_rows_all", "n_relative_control_rows_all",
    "n_invalid_relative_control_rows", "created_at"
  ),
  value = c(
    "v1", projection_dir, basename(script_file),
    uniqueN(pred_all$sample_key),
    uniqueN(pred_all[primary_external_projection_include == TRUE]$sample_key),
    nrow(pred_all), nrow(rel_all), nrow(invalid_relative),
    as.character(Sys.time())
  )
)
fwrite(
  manifest,
  file.path(projection_dir, "external_projection_run_manifest.tsv"),
  sep = "\t"
)

summary_path <- file.path(projection_dir, "external_projection_summary.md")
summary_lines <- readLines(summary_path, warn = FALSE)
marker <- "## Combined PRJNA863443 append status"
marker_index <- match(marker, summary_lines)
if (!is.na(marker_index)) summary_lines <- summary_lines[seq_len(marker_index - 1L)]
summary_lines <- c(
  summary_lines,
  "",
  marker,
  "",
  "- Dedicated PRJNA863443 per-dataset outputs were reused without modification.",
  paste0("- Combined projected samples: ", uniqueN(pred_all$sample_key), "."),
  paste0("- Combined prediction rows: ", nrow(pred_all), "."),
  paste0("- Combined relative-control rows: ", nrow(rel_all), "."),
  paste0("- Invalid relative-control rows after append: ", nrow(invalid_relative), ".")
)
writeLines(summary_lines, summary_path)

cat("PRJNA863443 existing-output append PASS\n")
cat("Combined samples:", uniqueN(pred_all$sample_key), "\n")
cat("Combined prediction rows:", nrow(pred_all), "\n")
cat("Combined invalid relative-control rows:", nrow(invalid_relative), "\n")
