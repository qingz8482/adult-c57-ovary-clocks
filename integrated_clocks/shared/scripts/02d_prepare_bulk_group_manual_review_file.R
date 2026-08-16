#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
if (!dir.exists(script_dir) || !file.exists(file.path(script_dir, "03_train_bulk_pseudobulk_clocks.R"))) {
  script_dir <- file.path(getwd(), "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
out_dir <- file.path(project_dir, "results", "alignment_audit", "bulk_group_level_alignment")
template_path <- file.path(out_dir, "bulk_group_manual_mapping_template.tsv")
review_path <- file.path(out_dir, "bulk_group_manual_mapping_TO_REVIEW.tsv")

if (!file.exists(template_path)) {
  stop("Manual mapping template is missing. Run scripts/02b_prepare_bulk_group_level_alignment.R first: ", template_path, call. = FALSE)
}

dt <- fread(template_path, sep = "\t")

ensure_col <- function(name, default = "") {
  if (!name %in% names(dt)) dt[, (name) := default]
}

for (col in c(
  "reviewer_initials", "review_date", "review_status",
  "biological_replicate_evidence", "technical_replicate_evidence",
  "gse279446_manual_biological_status"
)) {
  ensure_col(col, "")
}

truthy <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
dt[, manual_approval_TRUE_FALSE := truthy(manual_approval_TRUE_FALSE)]
dt[, include_in_primary_clock_TRUE_FALSE := truthy(include_in_primary_clock_TRUE_FALSE)]

dt[dataset_folder %in% c("GSE273724", "GSE279573") & manual_approval_TRUE_FALSE != TRUE,
   review_status := fifelse(
     include_in_primary_clock_TRUE_FALSE == TRUE,
     "needs_manual_approval_for_primary_clock",
     "excluded_from_primary_clock_review_context"
   )]

dt[dataset_folder == "GSE279446" &
     gse279446_manual_biological_status != "verified_control_like",
   `:=`(
     manual_approval_TRUE_FALSE = FALSE,
     include_in_primary_clock_TRUE_FALSE = FALSE,
     sample_unit_type = "unsafe_unknown",
     aggregation_method = "unsafe_do_not_use",
     review_status = "forced_excluded_until_verified_control_like",
     manual_note = fifelse(
       is.na(manual_note) | trimws(manual_note) == "",
       "excluded_pending_manual_review_gse279446_requires_verified_control_like_status",
       manual_note
     )
   )]

dt[, manual_approval_TRUE_FALSE := ifelse(manual_approval_TRUE_FALSE, "TRUE", "FALSE")]
dt[, include_in_primary_clock_TRUE_FALSE := ifelse(include_in_primary_clock_TRUE_FALSE, "TRUE", "FALSE")]

fwrite(dt, review_path, sep = "\t")

cat("Wrote manual review file:\n")
cat(review_path, "\n\n")
cat("Edit bulk_group_manual_mapping_TO_REVIEW.tsv manually.\n")
cat("Save the approved copy as bulk_group_manual_mapping_APPROVED.tsv.\n")
cat("Then run scripts/02c_validate_bulk_group_manual_mapping.R.\n")
