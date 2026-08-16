#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][[1]])
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
data_dir <- normalizePath(file.path(project_dir, "..", "data_input"), mustWork = TRUE)
projection_dir <- file.path(project_dir, "results", "external_projection")
dataset_result_dir <- file.path(projection_dir, "PRJNA1002222")
dataset_input_dir <- file.path(data_dir, "bulkseq_dataset", "PRJNA1002222")

read_tsv <- function(path) fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

meta <- fread(file.path(dataset_input_dir, "metadata", "PRJNA1002222_metadata.csv"))
counts <- fread(file.path(dataset_input_dir, "input", "PRJNA1002222_gene_counts_filtered.csv"))
logcpm <- fread(file.path(dataset_input_dir, "input", "PRJNA1002222_gene_logcpm_filtered.csv"))
readiness <- read_tsv(file.path(dataset_input_dir, "input", "PRJNA1002222_projection_readiness.tsv"))
source_validation <- read_tsv(file.path(dataset_input_dir, "input", "PRJNA1002222_supplementary_validation_summary.tsv"))
pred <- read_tsv(file.path(dataset_result_dir, "prjna1002222_projection_predictions.tsv"))
relative <- read_tsv(file.path(dataset_result_dir, "prjna1002222_relative_to_controls.tsv"))
pred_all <- read_tsv(file.path(projection_dir, "external_projection_predictions_all.tsv"))
relative_all <- read_tsv(file.path(projection_dir, "external_projection_relative_to_controls_all.tsv"))
preservation <- read_tsv(file.path(dataset_result_dir, "preexisting_output_preservation_audit.tsv"))

pred[, primary_external_projection_include := as_bool(primary_external_projection_include)]
relative[, matched_control_allowed := as_bool(matched_control_allowed)]
expected_groups <- data.table(
  projection_group = c("WT_2M", "WT_8M", "Cd38KO_8M", "WT_12M", "Cd38KO_12M"),
  expected_samples = c(5L, 5L, 5L, 3L, 3L)
)
observed_groups <- unique(pred[, .(sample_key, projection_group)])[, .(
  observed_samples = .N
), by = projection_group]
group_audit <- merge(expected_groups, observed_groups, by = "projection_group", all = TRUE)
group_pass <- all(group_audit$expected_samples == group_audit$observed_samples)

matched <- relative[matched_control_allowed == TRUE]
matched_8m <- matched[comparison_id == "PRJNA1002222_Cd38KO_vs_WT_8M"]
matched_12m <- matched[comparison_id == "PRJNA1002222_Cd38KO_vs_WT_12M"]
placement_2m <- relative[comparison_id == "PRJNA1002222_WT_natural_aging_reference_2M"]

audit <- data.table(
  check = c(
    "metadata_samples_21", "filtered_count_columns_21", "logcpm_columns_21",
    "matrix_columns_match_metadata", "finite_logcpm", "two_models_ready",
    "supplementary_validation_pass", "dedicated_prediction_rows_42",
    "dedicated_samples_21", "all_primary_projection", "group_sizes",
    "8M_controls_n5", "12M_controls_n3", "2M_no_relative_effect",
    "combined_samples_316", "combined_prediction_rows_632",
    "combined_relative_rows_632", "preexisting_outputs_preserved"
  ),
  expected = c(
    "21", "21", "21", "TRUE", "TRUE", "2", "TRUE", "42", "21", "TRUE",
    "TRUE", "5", "3", "TRUE", "316", "632", "632", "TRUE"
  ),
  observed = c(
    nrow(meta), ncol(counts) - 1L, ncol(logcpm) - 1L,
    identical(names(counts)[-1], meta$`Sample name`) && identical(names(logcpm)[-1], meta$`Sample name`),
    all(is.finite(as.matrix(logcpm[, -1, with = FALSE]))),
    sum(readiness$projection_ready == TRUE), all(source_validation$pass), nrow(pred),
    uniqueN(pred$sample_key), all(pred$primary_external_projection_include), group_pass,
    if (nrow(matched_8m)) unique(matched_8m$n_matched_controls) else NA,
    if (nrow(matched_12m)) unique(matched_12m$n_matched_controls) else NA,
    nrow(placement_2m) == 10L && all(!is.finite(placement_2m$relative_delta_age_week)),
    uniqueN(pred_all$sample_key), nrow(pred_all), nrow(relative_all), all(preservation$pass)
  )
)
# Mixed numeric/logical vectors are coerced to numeric by base R; normalize the
# logical audit fields back to their explicit textual representation.
audit[, observed := as.character(observed)]
audit[expected %chin% c("TRUE", "FALSE") & observed == "1", observed := "TRUE"]
audit[expected %chin% c("TRUE", "FALSE") & observed == "0", observed := "FALSE"]
audit[, pass := as.character(observed) == expected]
fwrite(group_audit, file.path(dataset_result_dir, "prjna1002222_group_size_audit.tsv"), sep = "\t")
fwrite(audit, file.path(dataset_result_dir, "prjna1002222_final_validation.tsv"), sep = "\t")
if (!all(audit$pass)) stop("PRJNA1002222 final validation failed.")
message("PASS: all PRJNA1002222 final validation checks succeeded.")
