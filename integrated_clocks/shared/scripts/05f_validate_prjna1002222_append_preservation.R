#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("--snapshot", "--validate")) {
  stop("Usage: 05f_validate_prjna1002222_append_preservation.R --snapshot|--validate")
}
mode <- sub("^--", "", args[[1]])
script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][[1]])
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
projection_dir <- file.path(project_dir, "results", "external_projection")
dataset_dir <- file.path(projection_dir, "PRJNA1002222")
snapshot_path <- file.path(dataset_dir, "preappend_preservation_snapshot.rds")
report_path <- file.path(dataset_dir, "preexisting_output_preservation_audit.tsv")
dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

combined_files <- c(
  "external_projection_predictions_all.tsv",
  "external_projection_predictions.tsv",
  "external_projection_relative_to_controls_all.tsv",
  "external_projection_relative_to_controls.tsv",
  "external_projection_reference_controls.tsv",
  "external_projection_control_pool_candidates.tsv",
  "external_projection_control_matching_audit.tsv",
  "external_projection_group_summary.tsv",
  "external_projection_model_concordance.tsv",
  "external_projection_gene_coverage.tsv",
  "external_projection_qc_flags.tsv",
  "pca_projection_scores.tsv",
  "pca_projection_pc_contributions.tsv",
  "external_projection_dataset_control_availability.tsv"
)

canonical_rows <- function(path) {
  x <- fread(path, sep = "\t", na.strings = NULL, colClasses = "character")
  if ("dataset_folder" %in% names(x)) x <- x[dataset_folder != "PRJNA1002222"]
  if (nrow(x) == 0L) return(character())
  sort(do.call(paste, c(x, sep = "\034")), method = "radix")
}

dataset_subdirs <- list.dirs(projection_dir, recursive = FALSE, full.names = TRUE)
dataset_subdirs <- dataset_subdirs[
  basename(dataset_subdirs) != "PRJNA1002222" &
    grepl("^(GSE|PRJNA|CRA|NIHMS)", basename(dataset_subdirs))
]
dataset_files <- sort(unlist(lapply(dataset_subdirs, function(path) {
  list.files(path, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
})))
dataset_files <- dataset_files[file.info(dataset_files)$isdir == FALSE]

if (mode == "snapshot") {
  snapshot <- list(
    combined = setNames(lapply(file.path(projection_dir, combined_files), canonical_rows), combined_files),
    dedicated_paths = dataset_files,
    dedicated_md5 = unname(tools::md5sum(dataset_files)),
    baseline_samples = uniqueN(fread(
      file.path(projection_dir, "external_projection_predictions_all.tsv"), sep = "\t"
    )$sample_key),
    baseline_prediction_rows = nrow(fread(
      file.path(projection_dir, "external_projection_predictions_all.tsv"), sep = "\t"
    ))
  )
  saveRDS(snapshot, snapshot_path)
  message(
    "Snapshot saved: ", snapshot$baseline_samples, " samples; ",
    snapshot$baseline_prediction_rows, " model rows; ", length(dataset_files),
    " files in pre-existing dataset-specific folders."
  )
  quit(save = "no", status = 0)
}

if (!file.exists(snapshot_path)) stop("Missing preservation snapshot: ", snapshot_path)
snapshot <- readRDS(snapshot_path)
combined_audit <- rbindlist(lapply(combined_files, function(filename) {
  observed <- canonical_rows(file.path(projection_dir, filename))
  expected <- snapshot$combined[[filename]]
  data.table(
    scope = "preexisting_combined_rows",
    target = filename,
    expected = length(expected),
    observed = length(observed),
    pass = identical(expected, observed)
  )
}))
current_exists <- file.exists(snapshot$dedicated_paths)
current_md5 <- rep(NA_character_, length(snapshot$dedicated_paths))
current_md5[current_exists] <- unname(tools::md5sum(snapshot$dedicated_paths[current_exists]))
dedicated_audit <- data.table(
  scope = "preexisting_dataset_specific_file",
  target = snapshot$dedicated_paths,
  expected = snapshot$dedicated_md5,
  observed = current_md5,
  pass = current_exists & snapshot$dedicated_md5 == current_md5
)
audit <- rbindlist(list(combined_audit, dedicated_audit), fill = TRUE)
fwrite(audit, report_path, sep = "\t")
if (!all(audit$pass)) {
  stop("Pre-existing external-projection outputs changed; see ", report_path)
}
message("PASS: all pre-existing combined rows and dataset-specific files are unchanged.")
