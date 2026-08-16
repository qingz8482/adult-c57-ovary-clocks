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
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

default_review_path <- file.path(out_dir, "bulk_group_manual_mapping_TO_REVIEW.tsv")
default_approved_path <- file.path(out_dir, "bulk_group_manual_mapping_APPROVED.tsv")
review_path <- Sys.getenv("BULK_GROUP_REVIEW_FILE", unset = default_review_path)
approved_path <- Sys.getenv("BULK_GROUP_APPROVED_FILE", unset = default_approved_path)

summary_tsv <- file.path(out_dir, "bulk_group_manual_approval_edit_summary.tsv")
summary_md <- file.path(out_dir, "bulk_group_manual_approval_edit_summary.md")
warnings_path <- file.path(out_dir, "bulk_group_manual_approval_edit_warnings.txt")

truthy <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")

safe_key <- function(...) {
  x <- paste(..., sep = "__")
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

normalize_label <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

append_note <- function(old, note) {
  old <- as.character(old)
  old[is.na(old)] <- ""
  has_note <- grepl(note, old, fixed = TRUE)
  out <- ifelse(trimws(old) == "", note, paste(old, note, sep = "; "))
  ifelse(has_note, old, out)
}

format_table_md <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return("_None._")
  dt <- as.data.table(dt)
  header <- paste("|", paste(names(dt), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(dt)), collapse = " | "), "|")
  rows <- apply(dt, 1, function(x) paste("|", paste(as.character(x), collapse = " | "), "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

if (!file.exists(review_path)) stop("Review mapping file is missing: ", review_path, call. = FALSE)

required_cols <- c(
  "dataset_folder", "metadata_group_label", "metadata_age_week",
  "metadata_control_like_clock_candidate", "expression_column", "sample_unit_type",
  "aggregation_id", "aggregation_method", "output_sample_key",
  "include_in_primary_clock_TRUE_FALSE", "manual_approval_TRUE_FALSE", "manual_note"
)

dt <- fread(review_path, sep = "\t", colClasses = "character")
missing_cols <- setdiff(required_cols, names(dt))
if (length(missing_cols) > 0) {
  stop("Review mapping file is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

if (file.exists(approved_path)) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  backup_path <- file.path(out_dir, paste0("bulk_group_manual_mapping_APPROVED.backup_", timestamp, ".tsv"))
  if (!file.copy(approved_path, backup_path, overwrite = FALSE)) {
    stop("Failed to back up existing approved mapping to: ", backup_path, call. = FALSE)
  }
} else {
  backup_path <- NA_character_
}

warnings <- character()
add_warning <- function(...) warnings <<- c(warnings, paste0(...))

dt[, metadata_control_like_clock_candidate_bool := truthy(metadata_control_like_clock_candidate)]
dt[, manual_approval_TRUE_FALSE := FALSE]
dt[, include_in_primary_clock_TRUE_FALSE := FALSE]

label_text <- do.call(paste, c(dt[, intersect(c(
  "metadata_row_id", "metadata_sample_key", "metadata_sample_id", "metadata_group_label"
), names(dt)), with = FALSE], sep = " "))
label_norm <- normalize_label(label_text)

is_gse273724_young <- dt$dataset_folder == "GSE273724" &
  dt$metadata_control_like_clock_candidate_bool &
  (grepl("young.*7.*8", tolower(label_text)) | grepl("young78", label_norm) | grepl("young78w", label_norm))

is_gse279573_wt <- dt$dataset_folder == "GSE279573" &
  dt$metadata_control_like_clock_candidate_bool &
  (grepl("wt.*chow.*fertile", tolower(label_text)) | grepl("wtchowfertile", label_norm))

target <- is_gse273724_young | is_gse279573_wt
if (!any(target)) {
  stop("No target approval rows found for GSE273724 young_7_8w or GSE279573 wt_chow_fertile.", call. = FALSE)
}

non_target_gse273724_279573 <- dt$dataset_folder %chin% c("GSE273724", "GSE279573") & !target
dt[non_target_gse273724_279573, manual_note := append_note(manual_note, "not_approved_for_primary_clock")]

dt[dataset_folder == "GSE279446", `:=`(
  manual_approval_TRUE_FALSE = FALSE,
  include_in_primary_clock_TRUE_FALSE = FALSE,
  sample_unit_type = "unsafe_unknown",
  aggregation_method = "unsafe_do_not_use",
  manual_note = append_note(manual_note, "kept_excluded; not manually verified as control-like biological samples")
)]
if ("review_status" %in% names(dt)) {
  dt[dataset_folder == "GSE279446", review_status := "excluded_not_verified_control_like"]
}

dt[target, `:=`(
  manual_approval_TRUE_FALSE = TRUE,
  include_in_primary_clock_TRUE_FALSE = TRUE,
  sample_unit_type = "biological_replicate",
  aggregation_method = "none_expand_replicates",
  aggregation_id = safe_key(metadata_sample_key, expression_column),
  output_sample_key = safe_key(metadata_sample_key, expression_column),
  manual_note = append_note(manual_note, "approved_control_like_biological_replicates")
)]
dt[is_gse273724_young, manual_note := append_note(manual_note, "only young_7_8w control-like rows included")]
dt[is_gse279573_wt, manual_note := append_note(manual_note, "validation split control-like rows included")]
if ("review_status" %in% names(dt)) {
  dt[target, review_status := "approved_for_primary_clock"]
}
if ("biological_replicate_evidence" %in% names(dt)) {
  dt[is_gse273724_young, biological_replicate_evidence := append_note(
    biological_replicate_evidence,
    "SRA title evidence: Ovary young_1 to young_4; control_or_natural_aging young_7_8w"
  )]
  dt[is_gse279573_wt, biological_replicate_evidence := append_note(
    biological_replicate_evidence,
    "SRA title evidence: WT, CD-fed, fertile, rep1 to rep2; control_or_natural_aging wt_chow_fertile"
  )]
}
if ("technical_replicate_evidence" %in% names(dt)) {
  dt[target, technical_replicate_evidence := append_note(
    technical_replicate_evidence,
    "not technical aggregation; each SRR retained as one biological-replicate model-ready sample"
  )]
}

approved <- dt[manual_approval_TRUE_FALSE == TRUE & include_in_primary_clock_TRUE_FALSE == TRUE]
if (nrow(approved) == 0) stop("Zero rows are approved after editing.", call. = FALSE)
if (anyDuplicated(approved$output_sample_key)) stop("Approved rows have duplicated output_sample_key.", call. = FALSE)
if (nrow(approved[is.na(expression_column) | trimws(expression_column) == ""]) > 0) {
  stop("Approved rows have missing expression_column.", call. = FALSE)
}
if (nrow(approved[is.na(metadata_age_week) | trimws(metadata_age_week) == ""]) > 0) {
  stop("Approved rows have missing metadata_age_week.", call. = FALSE)
}
if (nrow(approved[metadata_control_like_clock_candidate_bool != TRUE]) > 0) {
  stop("Approved rows include metadata_control_like_clock_candidate not TRUE.", call. = FALSE)
}
if (nrow(approved[sample_unit_type == "unsafe_unknown" | aggregation_method == "unsafe_do_not_use"]) > 0) {
  stop("Approved rows contain unsafe_unknown or unsafe_do_not_use.", call. = FALSE)
}

allowed_units <- c("biological_replicate", "technical_run", "group_mean")
allowed_methods <- c("none_expand_replicates", "sum_raw_counts", "mean_logcpm")
if (nrow(approved[!sample_unit_type %chin% allowed_units]) > 0) {
  stop("Approved rows contain invalid sample_unit_type.", call. = FALSE)
}
if (nrow(approved[!aggregation_method %chin% allowed_methods]) > 0) {
  stop("Approved rows contain invalid aggregation_method.", call. = FALSE)
}
if (nrow(approved[aggregation_method == "mean_logcpm" & !grepl("approved_group_mean_logcpm", manual_note, ignore.case = TRUE)]) > 0) {
  stop("Approved mean_logcpm rows require manual_note containing approved_group_mean_logcpm.", call. = FALSE)
}

dt[, metadata_control_like_clock_candidate_bool := NULL]
fwrite(dt, approved_path, sep = "\t", quote = TRUE, na = "")

summary <- approved[, .(
  approved_rows = .N,
  output_samples = uniqueN(output_sample_key),
  expression_columns = uniqueN(expression_column)
), by = .(dataset_folder, metadata_group_label, sample_unit_type, aggregation_method)][order(dataset_folder, metadata_group_label)]

fwrite(summary, summary_tsv, sep = "\t")
writeLines(if (length(warnings) > 0) unique(warnings) else "No approval edit warnings.", warnings_path)
writeLines(c(
  "# Bulk Group Manual Approval Edit Summary",
  "",
  paste0("- Review file: `", review_path, "`"),
  paste0("- Approved file: `", approved_path, "`"),
  paste0("- Backup file: `", ifelse(is.na(backup_path), "none", backup_path), "`"),
  paste0("- Approved rows: ", nrow(approved)),
  paste0("- GSE279446 approved rows: ", nrow(approved[dataset_folder == "GSE279446"])),
  "",
  "## Approved Rows",
  format_table_md(summary),
  "",
  "## Warnings",
  if (length(warnings) > 0) paste0("- ", unique(warnings)) else "_None._"
), summary_md)

print(summary)
message("BULK_GROUP_MANUAL_APPROVAL_EDIT_PASS")
