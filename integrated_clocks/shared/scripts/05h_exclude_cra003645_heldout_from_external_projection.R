#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- sub("^--file=", "", script_arg[[1]])
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
projection_dir <- file.path(project_dir, "results", "external_projection")
methods_dir <- file.path(project_dir, "methods_tables")

split_path <- file.path(methods_dir, "03_frozen_v1_model_ready_samples.tsv")
split <- fread(split_path, sep = "\t", na.strings = c("", "NA", "NaN"))
heldout_keys <- split[
  dataset_folder == "CRA003645" & split_role == "test",
  unique(sample_key)
]
if (length(heldout_keys) != 16L) {
  stop("Expected exactly 16 CRA003645 held-out test keys; observed ", length(heldout_keys))
}

clear_cra_controls <- function(dt) {
  if (!all(c("dataset_folder", "sample_key") %in% names(dt))) return(dt)
  idx <- !is.na(dt$dataset_folder) & dt$dataset_folder == "CRA003645" &
    !is.na(dt$sample_key) & !dt$sample_key %chin% heldout_keys
  idx[is.na(idx)] <- FALSE
  if (!any(idx)) return(dt)
  values <- list(
    used_in_deployment_training_flag = FALSE,
    matched_reference_control_flag = FALSE,
    matched_reference_control_source = "none",
    matched_control_allowed = FALSE,
    reference_control_group_id = "no_matched_control",
    reference_control_group_label = "No matched control",
    matched_control_method = "no_matched_control",
    matched_control_source = "none",
    n_matched_controls = 0L,
    matched_control_sample_keys = NA_character_,
    relative_control_reliability = "not_available",
    relative_control_warning = paste(
      "Adult CRA003645 references are reserved for the frozen held-out test;",
      "developmental/OOD projections have no matched-control effect."
    ),
    projection_dataset_note = paste(
      "Developmental/OOD placement only; 16 adult CRA003645 samples are excluded",
      "from external projection and retained exclusively as frozen held-out tests."
    )
  )
  for (column in intersect(names(values), names(dt))) {
    set(dt, which(idx), column, values[[column]])
  }
  for (column in intersect(
    c(
      "control_expected_predicted_age_week", "control_expected_delta_age_week",
      "relative_delta_age_week", "control_expected_pred_age",
      "control_expected_delta_age", "relative_delta_age"
    ),
    names(dt)
  )) {
    set(dt, which(idx), column, NA_real_)
  }
  dt
}

filter_file <- function(path) {
  dt <- fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
  if (!"sample_key" %in% names(dt)) return(0L)
  before <- nrow(dt)
  dt <- dt[!sample_key %chin% heldout_keys]
  dt <- clear_cra_controls(dt)
  fwrite(dt, path, sep = "\t", na = "")
  before - nrow(dt)
}

root_files <- list.files(projection_dir, pattern = "[.]tsv$", full.names = TRUE)
simplified_dir <- file.path(projection_dir, "simplified_classes")
simplified_files <- if (dir.exists(simplified_dir)) {
  list.files(simplified_dir, pattern = "[.]tsv$", full.names = TRUE, recursive = TRUE)
} else {
  character()
}
files <- c(root_files, simplified_files)
removed <- vapply(files, filter_file, integer(1))

pred_path <- file.path(projection_dir, "external_projection_predictions_all.tsv")
relative_path <- file.path(projection_dir, "external_projection_relative_to_controls_all.tsv")
pred <- fread(pred_path, sep = "\t", na.strings = c("", "NA", "NaN"))
relative <- fread(relative_path, sep = "\t", na.strings = c("", "NA", "NaN"))
if (uniqueN(pred$sample_key) != 300L || nrow(pred) != 600L) {
  stop("Filtered integrated projection must contain 300 samples and 600 model rows.")
}
if (any(pred$sample_key %chin% heldout_keys)) stop("Held-out CRA003645 key remains in predictions.")

availability_path <- file.path(
  projection_dir, "external_projection_dataset_control_availability.tsv"
)
availability <- fread(availability_path, sep = "\t", na.strings = c("", "NA", "NaN"))
availability[dataset_folder == "CRA003645", `:=`(
  n_projected_samples_all = 6L,
  n_projected_samples_primary = 6L,
  n_control_meta_candidates = 0L,
  n_projection_internal_control_candidates = 0L,
  has_same_dataset_control_meta = FALSE,
  has_projection_internal_control = FALSE,
  recommended_relative_method = "no_matched_control",
  expected_relative_method = "no_matched_control",
  observed_relative_method = "no_matched_control",
  dataset_note = paste(
    "Six developmental/OOD samples retained for descriptive placement only;",
    "16 adult samples reserved exclusively for frozen held-out testing."
  ),
  pass_fail_status = "PASS"
)]
fwrite(availability, availability_path, sep = "\t", na = "")

group_path <- file.path(projection_dir, "external_projection_group_summary.tsv")
group_summary <- fread(group_path, sep = "\t", na.strings = c("", "NA", "NaN"))
group_summary <- group_summary[dataset_folder != "CRA003645"]
cra <- pred[dataset_folder == "CRA003645"]
cra[, group := "developmental/OOD placement"]
cra_summary <- cra[, .(
  n_samples = uniqueN(sample_key),
  mean_age_week = mean(age_week),
  mean_predicted_age_week = mean(predicted_age_week),
  mean_delta_age_week = mean(delta_age_week),
  median_delta_age_week = median(delta_age_week),
  mean_relative_delta_age_week = NA_real_,
  median_relative_delta_age_week = NA_real_
), by = .(dataset_folder, projection_class, group, model_family, selector)]
for (column in setdiff(names(group_summary), names(cra_summary))) cra_summary[, (column) := NA]
setcolorder(cra_summary, names(group_summary))
fwrite(rbindlist(list(group_summary, cra_summary), use.names = TRUE), group_path, sep = "\t", na = "")

manifest_path <- file.path(projection_dir, "external_projection_run_manifest.tsv")
manifest <- fread(manifest_path, sep = "\t")
manifest[field == "n_projected_samples_all", value := "300"]
manifest[field == "n_projected_samples_primary", value := "286"]
manifest[field == "n_prediction_rows_all", value := "600"]
manifest[field == "n_relative_control_rows_all", value := as.character(nrow(relative))]
fwrite(manifest, manifest_path, sep = "\t")

summary_path <- file.path(projection_dir, "external_projection_summary.md")
summary_lines <- readLines(summary_path, warn = FALSE)
summary_lines <- gsub("All projected/reference samples: 316", "All projected/reference samples: 300", summary_lines, fixed = TRUE)
summary_lines <- gsub("Primary projected/reference samples: 302", "Primary projected/reference samples: 286", summary_lines, fixed = TRUE)
summary_lines <- gsub("Prediction rows: 632", "Prediction rows: 600", summary_lines, fixed = TRUE)
summary_lines <- c(
  summary_lines,
  "",
  "## CRA003645 held-out boundary",
  "",
  "- Sixteen adult CRA003645 samples are retained exclusively in the frozen held-out test and are excluded from external projection.",
  "- Six developmental/OOD CRA003645 samples remain as descriptive placements without matched-control relative effects."
)
writeLines(summary_lines, summary_path)

audit <- data.table(
  check = c(
    "heldout_key_count", "projection_sample_count", "prediction_row_count",
    "heldout_keys_absent", "cra_developmental_samples_retained",
    "cra_relative_effects_absent"
  ),
  observed = c(
    length(heldout_keys), uniqueN(pred$sample_key), nrow(pred),
    sum(pred$sample_key %chin% heldout_keys),
    uniqueN(pred[dataset_folder == "CRA003645"]$sample_key),
    sum(is.finite(relative[dataset_folder == "CRA003645"]$relative_delta_age_week))
  ),
  expected = c(16L, 300L, 600L, 0L, 6L, 0L)
)
audit[, status := fifelse(observed == expected, "PASS", "FAIL")]
audit[, detail := paste0("Removed rows across filtered TSVs: ", sum(removed))]
fwrite(
  audit,
  file.path(projection_dir, "cra003645_heldout_exclusion_audit.tsv"),
  sep = "\t"
)
if (any(audit$status != "PASS")) stop("CRA003645 held-out exclusion audit failed.")

cat("PASS: excluded 16 CRA003645 held-out samples; retained 300 projection samples.\n")
