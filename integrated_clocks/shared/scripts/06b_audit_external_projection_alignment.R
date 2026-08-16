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
results_dir <- file.path(project_dir, "results")
projection_dir <- file.path(results_dir, "external_projection")
metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
skipped_path <- file.path(projection_dir, "external_projection_skipped_unaligned_samples.tsv")
out_dir <- file.path(projection_dir, "alignment_audit")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

sanitize <- function(x) {
  x <- tolower(na_chr(x))
  x <- gsub("[^a-z0-9]+", "", x)
  x[x == ""] <- NA_character_
  x
}

collapse_unique <- function(x) {
  x <- unique(na_chr(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

read_required <- function(path, sep = "\t") {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
  fread(path, sep = sep, na.strings = c("", "NA", "NaN"))
}

read_expression_columns <- function(path) {
  if (!file.exists(path)) return(character())
  cols <- names(fread(path, nrows = 0))
  cols[-1]
}

field_values <- function(row, fields) {
  fields <- fields[fields %in% names(row)]
  out <- character()
  for (ff in fields) out <- c(out, as.character(row[[ff]]))
  unique(na_chr(out))
}

candidate_row <- function(row, source_row, expr_cols) {
  sample_key <- as.character(row$sample_key)
  sample_id <- as.character(row$sample_id %||% NA_character_)
  dataset_folder <- as.character(row$dataset_folder)

  accession_fields <- grep("accession|run|srr|srx|srs|sra|gsm", names(row), ignore.case = TRUE, value = TRUE)
  source_accession_fields <- grep("accession|run|srr|srx|srs|sra|gsm", names(source_row), ignore.case = TRUE, value = TRUE)
  label_fields <- intersect(c(
    "sample_key", "sample_id", "sample_title", "sample_name", "sample_context",
    "treatment", "genotype", "condition", "group", "Type", "Sample title", "Sample name",
    "Treatment", "Genotype", "Phenotype"
  ), names(row))
  source_label_fields <- intersect(c(
    "sample_key", "sample_id", "sample_title", "sample_name", "sample_context",
    "treatment", "genotype", "condition", "group", "Type", "Sample title", "Sample name",
    "Treatment", "Genotype", "Phenotype"
  ), names(source_row))

  tier <- "no_match"
  status <- "no_match"
  matched <- character()
  evidence <- character()

  tier1 <- intersect(sample_key, expr_cols)
  if (length(tier1) == 1) {
    tier <- "tier1_exact_sample_key"
    status <- "safe_exact"
    matched <- tier1
    evidence <- paste0("sample_key exactly equals expression column ", tier1)
  } else {
    accession_values <- unique(c(field_values(row, accession_fields), field_values(source_row, source_accession_fields)))
    accession_values <- accession_values[grepl("^(SRR|ERR|DRR|SRX|SRS|SAMN|GSM)", accession_values, ignore.case = TRUE)]
    tier2 <- intersect(accession_values, expr_cols)
    if (length(tier2) == 1) {
      tier <- "tier2_exact_run_accession"
      status <- "safe_exact"
      matched <- tier2
      evidence <- paste0("run/accession exactly equals expression column ", tier2)
    } else if (length(tier2) > 1) {
      tier <- "tier2_exact_run_accession"
      status <- "ambiguous"
      matched <- tier2
      evidence <- "multiple exact run/accession matches"
    } else {
      expr_s <- sanitize(expr_cols)
      names(expr_s) <- expr_cols
      sample_key_s <- sanitize(sample_key)
      tier3 <- names(expr_s)[!is.na(expr_s) & expr_s == sample_key_s]
      if (length(tier3) == 1) {
        tier <- "tier3_sanitized_sample_key"
        status <- "safe_sanitized"
        matched <- tier3
        evidence <- paste0("sanitized sample_key matches expression column ", tier3)
      } else if (length(tier3) > 1) {
        tier <- "tier3_sanitized_sample_key"
        status <- "ambiguous"
        matched <- tier3
        evidence <- "sanitized sample_key matches multiple expression columns"
      } else {
        label_values <- unique(c(field_values(row, label_fields), field_values(source_row, source_label_fields)))
        label_s <- sanitize(label_values)
        label_s <- label_s[!is.na(label_s) & nchar(label_s) >= 4]
        tier4 <- names(expr_s)[!is.na(expr_s) & expr_s %in% label_s]
        if (length(tier4) == 1) {
          tier <- "tier4_sanitized_sample_title_or_name"
          status <- "safe_sanitized"
          matched <- tier4
          evidence <- paste0("sanitized sample label/title matches expression column ", tier4)
        } else if (length(tier4) > 1) {
          tier <- "tier4_sanitized_sample_title_or_name"
          status <- "ambiguous"
          matched <- tier4
          evidence <- "sanitized sample labels/titles match multiple expression columns"
        } else {
          group_tokens <- unique(unlist(strsplit(tolower(paste(label_values, collapse = " ")), "[^a-z0-9]+")))
          group_tokens <- group_tokens[!is.na(group_tokens) & nchar(group_tokens) >= 4]
          group_tokens <- group_tokens[!group_tokens %chin% c("sample", "group", "level", "metadata", "ovary", "week", "mice", "mouse")]
          token_hits <- unique(unlist(lapply(group_tokens, function(tok) {
            expr_cols[grepl(tok, tolower(expr_cols), fixed = TRUE)]
          })))
          if (length(token_hits) == 1) {
            tier <- "tier5_group_label_only"
            status <- "unsafe_group_only"
            matched <- token_hits
            evidence <- "only a group-label token matched an expression column"
          } else if (length(token_hits) > 1) {
            tier <- "tier5_group_label_only"
            status <- "ambiguous"
            matched <- token_hits
            evidence <- "group-label tokens match multiple expression columns"
          } else {
            tier <- "no_match"
            status <- "no_match"
            matched <- character()
            evidence <- "no exact, accession, sanitized, or group-label expression-column match"
          }
        }
      }
    }
  }

  data.table(
    dataset_folder = dataset_folder,
    metadata_sample_key = sample_key,
    metadata_sample_id = sample_id,
    metadata_sample_label = as.character(row$sample_title %||% row$treatment %||% sample_id),
    metadata_age_week = as.numeric(row$age_week),
    metadata_projection_class = as.character(row$projection_class %||% NA_character_),
    metadata_treatment = as.character(row$treatment %||% NA_character_),
    candidate_expression_column = collapse_unique(matched),
    candidate_SRR = collapse_unique(matched[grepl("^SRR", matched, ignore.case = TRUE)]),
    match_tier = tier,
    match_status = status,
    evidence = paste(evidence, collapse = "; "),
    n_expression_columns_in_dataset = length(expr_cols),
    expression_columns_available = paste(expr_cols, collapse = ";"),
    metadata_accession_fields_available = paste(unique(c(accession_fields, source_accession_fields)), collapse = ";"),
    metadata_label_fields_available = paste(unique(c(label_fields, source_label_fields)), collapse = ";")
  )
}

skipped <- read_required(skipped_path)
if (nrow(skipped) == 0) {
  candidate_dt <- data.table(
    dataset_folder = character(),
    metadata_sample_key = character(),
    metadata_sample_id = character(),
    metadata_sample_label = character(),
    metadata_age_week = numeric(),
    metadata_projection_class = character(),
    metadata_treatment = character(),
    candidate_expression_column = character(),
    candidate_SRR = character(),
    match_tier = character(),
    match_status = character(),
    evidence = character(),
    n_expression_columns_in_dataset = integer(),
    expression_columns_available = character(),
    metadata_accession_fields_available = character(),
    metadata_label_fields_available = character()
  )
  fwrite(
    candidate_dt,
    file.path(out_dir, "skipped_alignment_candidate_matches.tsv"),
    sep = "\t",
    na = "NA"
  )

  summary_dt <- data.table(
    dataset_folder = "ALL",
    n_skipped = 0L,
    n_safe_exact = 0L,
    n_safe_sanitized = 0L,
    n_ambiguous = 0L,
    n_unsafe_group_only = 0L,
    n_no_match = 0L,
    can_rerun_projection_automatically = TRUE
  )
  fwrite(summary_dt, file.path(out_dir, "skipped_alignment_summary.tsv"), sep = "\t", na = "NA")

  template <- data.table(
    dataset_folder = character(),
    metadata_sample_key = character(),
    metadata_sample_label = character(),
    metadata_age_week = numeric(),
    metadata_projection_class = character(),
    metadata_treatment = character(),
    candidate_expression_column = character(),
    candidate_SRR = character(),
    match_tier = character(),
    match_status = character(),
    use_for_projection_TRUE_FALSE = logical(),
    manual_note = character()
  )
  fwrite(template, file.path(out_dir, "manual_projection_sample_mapping_template.tsv"), sep = "\t", na = "NA")

  manual_lines <- c(
    "# Skipped External Projection Alignment Manual Review",
    "",
    "No samples are currently skipped by external projection alignment.",
    "",
    "The existing approved manual mapping file is left intact because it is used by the projection script to align group-level metadata rows to SRR expression columns."
  )
  writeLines(manual_lines, file.path(out_dir, "skipped_alignment_manual_review.md"))

  cat("number skipped: 0\n")
  cat("projection can be rerun automatically: TRUE\n")
  cat("output directory:", out_dir, "\n")
  quit(status = 0)
}
meta <- fread(metadata_path, na.strings = c("", "NA", "NaN"))
meta[, age_week := as.numeric(age_week)]
skipped_full <- merge(
  skipped,
  meta,
  by = intersect(names(skipped), names(meta)),
  all.x = TRUE,
  suffixes = c("", ".meta")
)
if (!"sample_key" %in% names(skipped_full)) stop("Skipped table lacks sample_key.", call. = FALSE)

candidate_matches <- list()
dataset_sections <- list()

for (ds in sort(unique(skipped_full$dataset_folder))) {
  ds_rows <- skipped_full[dataset_folder == ds]
  expr_paths <- unique(na_chr(ds_rows$expression_logcpm_path))
  expr_paths <- expr_paths[!is.na(expr_paths)]
  expr_cols <- if (length(expr_paths) == 1) read_expression_columns(expr_paths[[1]]) else character()

  source_paths <- unique(na_chr(ds_rows$metadata_source_path))
  source_paths <- source_paths[!is.na(source_paths)]
  source_dt <- data.table()
  if (length(source_paths) == 1 && file.exists(source_paths[[1]])) {
    source_dt <- fread(source_paths[[1]], na.strings = c("", "NA", "NaN"))
    if ("dataset_folder" %in% names(source_dt)) {
      source_dt <- source_dt[as.character(dataset_folder) == ds]
    } else if ("dataset_id" %in% names(source_dt)) {
      source_dt <- source_dt[as.character(dataset_id) == ds]
    }
  }

  for (i in seq_len(nrow(ds_rows))) {
    row <- ds_rows[i]
    source_row <- data.table()
    if (nrow(source_dt) > 0) {
      source_keys <- intersect(c("Sample name", "sample_id", "sample_key", "sample_title", "Sample title"), names(source_dt))
      row_values <- field_values(row, c("sample_id", "sample_key", "sample_title", "treatment"))
      source_row <- source_dt[Reduce(`|`, lapply(source_keys, function(kk) as.character(get(kk)) %chin% row_values))]
      if (nrow(source_row) == 0 && "Sample name" %in% names(source_dt) && "sample_id" %in% names(row)) {
        source_row <- source_dt[as.character(`Sample name`) == as.character(row$sample_id)]
      }
      if (nrow(source_row) == 0) source_row <- source_dt[0]
      if (nrow(source_row) > 1) source_row <- source_row[1]
    }
    candidate_matches[[length(candidate_matches) + 1L]] <- candidate_row(row, source_row, expr_cols)
  }

  accession_fields <- grep("accession|run|srr|srx|srs|sra|gsm", names(ds_rows), ignore.case = TRUE, value = TRUE)
  source_accession_fields <- if (nrow(source_dt) > 0) grep("accession|run|srr|srx|srs|sra|gsm", names(source_dt), ignore.case = TRUE, value = TRUE) else character()
  label_fields <- grep("sample|title|name|group|condition|treatment|genotype|type|phenotype", names(ds_rows), ignore.case = TRUE, value = TRUE)
  source_label_fields <- if (nrow(source_dt) > 0) grep("sample|title|name|group|condition|treatment|genotype|type|phenotype", names(source_dt), ignore.case = TRUE, value = TRUE) else character()

  dataset_sections[[ds]] <- list(
    skipped = ds_rows,
    source = source_dt,
    expression_columns = expr_cols,
    expression_path = expr_paths,
    metadata_source_path = source_paths,
    accession_fields = unique(c(accession_fields, source_accession_fields)),
    label_fields = unique(c(label_fields, source_label_fields))
  )
}

candidate_dt <- rbindlist(candidate_matches, fill = TRUE)
setorder(candidate_dt, dataset_folder, metadata_sample_key)
fwrite(
  candidate_dt,
  file.path(out_dir, "skipped_alignment_candidate_matches.tsv"),
  sep = "\t",
  na = "NA"
)

summary_dt <- candidate_dt[, .(
  n_skipped = .N,
  n_safe_exact = sum(match_status == "safe_exact"),
  n_safe_sanitized = sum(match_status == "safe_sanitized"),
  n_ambiguous = sum(match_status == "ambiguous"),
  n_unsafe_group_only = sum(match_status == "unsafe_group_only"),
  n_no_match = sum(match_status == "no_match"),
  can_rerun_projection_automatically = all(match_status %chin% c("safe_exact", "safe_sanitized"))
), by = dataset_folder]
summary_dt <- rbind(
  summary_dt,
  candidate_dt[, .(
    n_skipped = .N,
    n_safe_exact = sum(match_status == "safe_exact"),
    n_safe_sanitized = sum(match_status == "safe_sanitized"),
    n_ambiguous = sum(match_status == "ambiguous"),
    n_unsafe_group_only = sum(match_status == "unsafe_group_only"),
    n_no_match = sum(match_status == "no_match"),
    can_rerun_projection_automatically = all(match_status %chin% c("safe_exact", "safe_sanitized"))
  )][, dataset_folder := "ALL"],
  fill = TRUE
)
setcolorder(summary_dt, c("dataset_folder", setdiff(names(summary_dt), "dataset_folder")))
fwrite(
  summary_dt,
  file.path(out_dir, "skipped_alignment_summary.tsv"),
  sep = "\t",
  na = "NA"
)

template <- candidate_dt[, .(
  dataset_folder,
  metadata_sample_key,
  metadata_sample_label,
  metadata_age_week,
  metadata_projection_class,
  metadata_treatment,
  candidate_expression_column,
  candidate_SRR,
  match_tier,
  match_status,
  use_for_projection_TRUE_FALSE = FALSE,
  manual_note = fifelse(
    match_status %chin% c("safe_exact", "safe_sanitized"),
    "review optional; automatic match considered safe",
    "manual SRR/sample mapping required before projection"
  )
)]
fwrite(
  template,
  file.path(out_dir, "manual_projection_sample_mapping_template.tsv"),
  sep = "\t",
  na = "NA"
)

auto_ok <- all(candidate_dt$match_status %chin% c("safe_exact", "safe_sanitized"))
if (auto_ok) {
  approved <- copy(template)
  approved[, use_for_projection_TRUE_FALSE := TRUE]
  approved[, manual_note := "auto-approved because every skipped sample had safe exact/sanitized evidence"]
  fwrite(
    approved,
    file.path(out_dir, "manual_projection_sample_mapping_APPROVED_AUTO.tsv"),
    sep = "\t",
    na = "NA"
  )
}

dataset_md <- unlist(lapply(names(dataset_sections), function(ds) {
  section <- dataset_sections[[ds]]
  skipped_cols <- intersect(c("sample_key", "sample_id", "age_week", "projection_class", "projection_reason", "treatment", "genotype", "sample_title"), names(section$skipped))
  skipped_lines <- apply(section$skipped[, ..skipped_cols], 1, function(x) paste(paste(skipped_cols, x, sep = "="), collapse = "; "))
  source_preview <- character()
  if (nrow(section$source) > 0) {
    source_cols <- intersect(c("Sample name", "Accession", "Project accession", "Sample title", "Treatment", "Genotype", "Type", "dataformat", "Age raw"), names(section$source))
    source_preview <- apply(head(section$source[, ..source_cols], 12), 1, function(x) paste(paste(source_cols, x, sep = "="), collapse = "; "))
  }
  c(
    paste0("## ", ds),
    "",
    paste0("- Expression path: ", paste(section$expression_path, collapse = "; ")),
    paste0("- Metadata source path: ", paste(section$metadata_source_path, collapse = "; ")),
    paste0("- Expression columns available: ", paste(section$expression_columns, collapse = ", ")),
    paste0("- SRR/run/accession fields available: ", ifelse(length(section$accession_fields) == 0, "none", paste(section$accession_fields, collapse = ", "))),
    paste0("- Sample title/name/group/condition/treatment fields available: ", ifelse(length(section$label_fields) == 0, "none", paste(section$label_fields, collapse = ", "))),
    "",
    "Skipped metadata rows:",
    paste0("- ", skipped_lines),
    "",
    "Source metadata preview:",
    if (length(source_preview) == 0) "- Source metadata unavailable or empty." else paste0("- ", source_preview),
    ""
  )
}), use.names = FALSE)

manual_lines <- c(
  "# Skipped External Projection Alignment Manual Review",
  "",
  "This audit does not retrain clocks and does not project samples. It only compares skipped metadata labels against saved expression matrix columns.",
  "",
  "## Decision",
  paste0("- Skipped samples: ", nrow(candidate_dt), "."),
  paste0("- safe_exact: ", sum(candidate_dt$match_status == "safe_exact"), "."),
  paste0("- safe_sanitized: ", sum(candidate_dt$match_status == "safe_sanitized"), "."),
  paste0("- ambiguous: ", sum(candidate_dt$match_status == "ambiguous"), "."),
  paste0("- unsafe_group_only: ", sum(candidate_dt$match_status == "unsafe_group_only"), "."),
  paste0("- no_match: ", sum(candidate_dt$match_status == "no_match"), "."),
  paste0("- Projection can be rerun automatically: ", auto_ok, "."),
  "",
  "Automatic approval is allowed only when every skipped sample has `safe_exact` or `safe_sanitized` evidence. Group-label-only or missing SRR evidence requires manual mapping.",
  "",
  "## Candidate Matches",
  paste0(
    "- ",
    candidate_dt$metadata_sample_key,
    ": ",
    candidate_dt$match_status,
    " via ",
    candidate_dt$match_tier,
    "; candidate=",
    candidate_dt$candidate_expression_column %||% "NA",
    "; evidence=",
    candidate_dt$evidence
  ),
  "",
  dataset_md
)
writeLines(manual_lines, file.path(out_dir, "skipped_alignment_manual_review.md"))

cat("number skipped:", nrow(candidate_dt), "\n")
cat("number safe_exact:", sum(candidate_dt$match_status == "safe_exact"), "\n")
cat("number safe_sanitized:", sum(candidate_dt$match_status == "safe_sanitized"), "\n")
cat("number ambiguous:", sum(candidate_dt$match_status == "ambiguous"), "\n")
cat("number unsafe_group_only:", sum(candidate_dt$match_status == "unsafe_group_only"), "\n")
cat("number no_match:", sum(candidate_dt$match_status == "no_match"), "\n")
cat("projection can be rerun automatically:", auto_ok, "\n")
cat("output directory:", out_dir, "\n")
