#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(glmnet)
})

options(stringsAsFactors = FALSE)
set.seed(20260630)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

ensure_dt <- function(x) {
  if (!data.table::is.data.table(x)) {
    x <- data.table::as.data.table(x)
  }
  x
}

manifest_row <- function(section, key, value) {
  out <- data.table(section = section, manifest_key = key, value = as.character(value))
  setnames(out, "manifest_key", "key")
  out
}

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else NA_character_
if (is.na(script_dir) || !dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
data_input_dir <- normalizePath(file.path(root_dir, "data_input"), mustWork = TRUE)
metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
selector_path <- file.path(project_dir, "scripts", "02_age_dataset_penalized_feature_selection.R")
stopifnot(file.exists(metadata_path), file.exists(selector_path))
source(selector_path)

default_results_dir <- file.path(project_dir, "results")
clock_results_override <- Sys.getenv("CLOCK_RESULTS_DIR", unset = "")
results_dir <- if (nzchar(clock_results_override)) {
  normalizePath(clock_results_override, mustWork = FALSE)
} else {
  default_results_dir
}
run_label <- Sys.getenv("RUN_LABEL", unset = "integrated_bulk_pseudobulk")
pca_dir <- file.path(results_dir, "pca_glmnet_clock")
elasticnet_dir <- file.path(results_dir, "elasticnet_clock")
selected_gene_dir <- file.path(results_dir, "selected_genes")
plot_dir <- file.path(results_dir, "readable_plots")
artifact_root_dir <- file.path(results_dir, "model_artifacts")
benchmark_artifact_dir <- file.path(artifact_root_dir, "benchmark")
deployment_artifact_dir <- file.path(artifact_root_dir, "deployment")
default_bulk_group_mapping_path <- file.path(
  default_results_dir,
  "alignment_audit",
  "bulk_group_level_alignment",
  "bulk_group_manual_mapping_APPROVED.tsv"
)
bulk_group_mapping_path <- Sys.getenv("BULK_GROUP_MAPPING_FILE", unset = default_bulk_group_mapping_path)
bulk_group_mapping_mode <- tolower(Sys.getenv("BULK_GROUP_MAPPING_MODE", unset = "auto"))
allowed_bulk_group_mapping_modes <- c("auto", "disabled", "force")
if (!bulk_group_mapping_mode %chin% allowed_bulk_group_mapping_modes) {
  stop(
    "Invalid BULK_GROUP_MAPPING_MODE: ", bulk_group_mapping_mode,
    ". Allowed values: ", paste(allowed_bulk_group_mapping_modes, collapse = ", "),
    call. = FALSE
  )
}
bulk_group_mapping_loaded <- if (file.exists(bulk_group_mapping_path)) {
  message("  Loading approved bulk group mapping: ", bulk_group_mapping_path)
  fread(bulk_group_mapping_path, sep = "\t")
} else {
  data.table()
}
bulk_group_mapping <- if (bulk_group_mapping_mode == "disabled") {
  data.table()
} else {
  bulk_group_mapping_loaded
}
if (!nzchar(clock_results_override) && Sys.getenv("ALLOW_DELETE_RESULTS", unset = "FALSE") != "TRUE") {
  stop(
    "Refusing to delete/rebuild results/. Set CLOCK_RESULTS_DIR for a versioned run or set ALLOW_DELETE_RESULTS=TRUE.",
    call. = FALSE
  )
}
if (!nzchar(clock_results_override) && dir.exists(results_dir)) {
  unlink(results_dir, recursive = TRUE, force = TRUE)
}
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(elasticnet_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_gene_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(benchmark_artifact_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deployment_artifact_dir, recursive = TRUE, showWarnings = FALSE)

run_status_path <- file.path(results_dir, "run_manifest.tsv")
run_status_md_path <- file.path(results_dir, "run_manifest.md")
current_stage <- "initializing"

write_run_status <- function(status, failure_stage = NA_character_, failure_message = NA_character_) {
  rows <- manifest_row(
    "run",
    c("run_label", "timestamp", "script_path", "results_dir", "status", "failure_stage", "failure_message"),
    c(
      run_label,
      as.character(Sys.time()),
      normalizePath(script_file %||% file.path(script_dir, "03_train_bulk_pseudobulk_clocks.R"), mustWork = FALSE),
      results_dir,
      status,
      failure_stage,
      failure_message
    )
  )
  fwrite(rows, run_status_path, sep = "\t")
  writeLines(c(
    paste0("# Run Manifest: ", run_label),
    "",
    paste0("- Timestamp: ", Sys.time()),
    paste0("- Results dir: `", results_dir, "`"),
    paste0("- Status: ", status),
    if (!is.na(failure_stage)) paste0("- Failure stage: ", failure_stage) else NULL,
    if (!is.na(failure_message)) paste0("- Failure message: ", failure_message) else NULL
  ), run_status_md_path)
  invisible(rows)
}

set_run_stage <- function(stage) {
  current_stage <<- stage
  write_run_status(stage)
  invisible(stage)
}

options(error = function() {
  write_run_status("failed", failure_stage = current_stage, failure_message = geterrmessage())
  q(status = 1, save = "no")
})
set_run_stage("started")

alpha_grid <- c(0, 0.5, 1)
max_pc <- 10L
hvg_n_genes <- 2000L
version_b_n_genes <- age_dataset_penalized_defaults$n_genes
training_scope <- "integrated_bulk_pseudobulk"

bulk_group_audit_dir <- file.path(results_dir, "alignment_audit", "bulk_group_level_alignment")
bulk_group_added_samples <- list()
bulk_group_expression_column_map <- list()
bulk_group_dataset_summary <- list()
bulk_group_gene_qc <- list()
bulk_group_warnings <- character()
bulk_group_mapping_used_datasets <- character()
exact_samplelevel_matching_used_for <- character()
samplelevel_exact_priority_datasets <- c("GSE273724", "GSE279446", "GSE279573")

copy_alignment_audit_inputs <- function() {
  source_dir <- file.path(default_results_dir, "alignment_audit", "bulk_group_level_alignment")
  dest_dir <- file.path(results_dir, "alignment_audit", "bulk_group_level_alignment")
  if (!nzchar(clock_results_override) || !dir.exists(source_dir)) return(invisible(NULL))
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(source_dir, full.names = TRUE)
  if (length(files) > 0) file.copy(files, dest_dir, overwrite = TRUE, copy.date = TRUE)
  invisible(NULL)
}

copy_alignment_audit_inputs()

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("^_|_$", "", x)
}

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

gene_id_symbol_map_cache <- NULL

load_gene_id_symbol_map <- function() {
  if (!is.null(gene_id_symbol_map_cache)) return(gene_id_symbol_map_cache)
  feature_paths <- unique(c(
    Sys.glob(file.path(
      data_input_dir, "standardized_seurat", "01_count_matrices",
      "*", "*", "filtered_feature_bc_matrix", "features.tsv.gz"
    )),
    Sys.glob(file.path(
      data_input_dir, "scRNA_dataset", "*", "cellranger_counts",
      "*", "outs", "filtered_feature_bc_matrix", "features.tsv.gz"
    ))
  ))
  feature_paths <- feature_paths[file.exists(feature_paths)]
  maps <- list()
  for (feature_path in feature_paths) {
    dt <- tryCatch(
      fread(cmd = paste("gzip -cd", shQuote(feature_path)), header = FALSE),
      error = function(e) NULL
    )
    if (is.null(dt) || ncol(dt) < 2) next
    maps[[length(maps) + 1L]] <- dt[, .(gene_id = na_chr(V1), gene_symbol = na_chr(V2))]
  }
  if (length(maps) == 0) {
    gene_id_symbol_map_cache <<- data.table(gene_id = character(), gene_symbol = character())
    return(gene_id_symbol_map_cache)
  }
  mapping <- rbindlist(maps, use.names = TRUE, fill = TRUE)
  mapping <- mapping[
    grepl("^ENSMUSG", gene_id) &
      !is.na(gene_symbol) &
      gene_symbol != "" &
      !grepl("^ENSMUSG", gene_symbol)
  ]
  mapping <- mapping[, .N, by = .(gene_id, gene_symbol)]
  setorder(mapping, gene_id, -N, gene_symbol)
  mapping <- mapping[, .SD[1], by = gene_id][, .(gene_id, gene_symbol)]
  gene_id_symbol_map_cache <<- mapping
  mapping
}

map_ensembl_rows_to_symbols <- function(mat, path) {
  gene_ids <- rownames(mat)
  if (!any(grepl("^ENSMUSG", gene_ids))) return(mat)
  mapping <- load_gene_id_symbol_map()
  if (nrow(mapping) == 0) return(mat)
  symbol_lookup <- mapping$gene_symbol
  names(symbol_lookup) <- mapping$gene_id
  lookup_gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
  mapped_symbols <- unname(symbol_lookup[lookup_gene_ids])
  mapped <- !is.na(mapped_symbols) & mapped_symbols != ""
  if (!any(mapped)) return(mat)
  new_gene_ids <- gene_ids
  new_gene_ids[mapped] <- mapped_symbols[mapped]
  if (!anyDuplicated(new_gene_ids)) {
    rownames(mat) <- new_gene_ids
    message("Mapped ", sum(mapped), " Ensembl gene rows to symbols for ", basename(path), ".")
    return(mat)
  }
  groups <- split(seq_along(new_gene_ids), new_gene_ids)
  collapsed <- t(vapply(
    groups,
    function(idx) colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(mat))
  ))
  colnames(collapsed) <- colnames(mat)
  message(
    "Mapped ", sum(mapped), " Ensembl gene rows to symbols and collapsed ",
    length(new_gene_ids) - length(groups), " duplicate rows for ", basename(path), "."
  )
  collapsed
}

read_expression_matrix <- function(path) {
  dt <- fread(path)
  gene_col <- names(dt)[1]
  genes <- as.character(dt[[gene_col]])
  dt[, (gene_col) := NULL]
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- genes
  keep <- !is.na(rownames(mat)) & rownames(mat) != "" & !duplicated(rownames(mat))
  mat <- mat[keep, , drop = FALSE]
  map_ensembl_rows_to_symbols(mat, path)
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

text_any <- function(dt, cols) {
  cols <- intersect(cols, names(dt))
  if (length(cols) == 0) return(rep("", nrow(dt)))
  do.call(paste, c(lapply(dt[, ..cols], function(x) fifelse(is.na(x), "", as.character(x))), sep = " "))
}

add_dataset_subcohort <- function(meta_dt) {
  meta_dt <- copy(meta_dt)
  if (!"dataset_subcohort" %in% names(meta_dt)) {
    meta_dt[, dataset_subcohort := dataset_folder]
  }
  meta_dt[is.na(dataset_subcohort) | dataset_subcohort == "", dataset_subcohort := dataset_folder]

  gse279446 <- meta_dt$dataset_folder == "GSE279446"
  if (any(gse279446)) {
    row_text <- text_any(meta_dt, c("sample_id", "accession", "sample_key", "sample_title", "genotype", "treatment"))
    wt_ovary <- gse279446 & grepl("(^|[^A-Za-z0-9_])WT_ovary[0-9]", row_text)
    zeb1_wt_young <- gse279446 &
      grepl("Zeb1_WT_ovary", row_text, ignore.case = FALSE) &
      grepl("youngCohort", row_text, ignore.case = FALSE)
    meta_dt[gse279446, dataset_subcohort := "GSE279446__other_review"]
    meta_dt[wt_ovary, dataset_subcohort := "GSE279446__WT_ovary_cohort"]
    meta_dt[zeb1_wt_young, dataset_subcohort := "GSE279446__Zeb1_WT_youngCohort"]
  }
  meta_dt
}

run_severe_metadata_checks <- function(meta_dt) {
  issues <- list()
  add_issue <- function(dataset_folder, severity, issue, n_rows = NA_integer_, sample_keys = character()) {
    issues[[length(issues) + 1L]] <<- data.table(
      dataset_folder = dataset_folder,
      severity = severity,
      issue = issue,
      n_rows = as.integer(n_rows),
      sample_keys = paste(sample_keys, collapse = ";")
    )
  }

  gse279573 <- meta_dt[dataset_folder == "GSE279573"]
  if (nrow(gse279573) > 0) {
    ko_text <- text_any(gse279573, c("sample_key", "sample_id", "accession", "sample_title", "genotype", "treatment", "sample_context"))
    bad <- gse279573[control_like_clock_candidate == TRUE &
      grepl("KO|knockout|rescue|infertile|vegetable-oil-free", ko_text, ignore.case = TRUE)]
    if (nrow(bad) > 0) {
      add_issue("GSE279573", "severe", "KO/rescue/infertile row marked control_like_clock_candidate TRUE.", nrow(bad), bad$sample_key)
    }
    expected_wt <- c("SRR31013339", "SRR31013344")
    actual <- sort(gse279573[control_like_clock_candidate == TRUE, accession])
    if (!identical(actual, expected_wt)) {
      add_issue("GSE279573", "severe", "Primary candidates are not exactly the two WT chow fertile rows.", length(actual), actual)
    }
  }

  gse279446 <- meta_dt[dataset_folder == "GSE279446"]
  if (nrow(gse279446) > 0) {
    primary <- gse279446[control_like_clock_candidate == TRUE]
    row_text <- text_any(primary, c("sample_key", "sample_id", "accession", "sample_title", "genotype", "treatment", "sample_context"))
    suspicious <- primary[grepl("KO|knockout|mutant|perturb|treatment|stress|diet|rescue|infertile", row_text, ignore.case = TRUE)]
    if (nrow(suspicious) > 0) {
      add_issue("GSE279446", "warning", "Control-like row contains suspicious intervention/genetic-model term.", nrow(suspicious), suspicious$sample_key)
    }
    unclear <- primary[!grepl("(^|[^A-Za-z0-9])WT([^A-Za-z0-9]|$)|Zeb1_WT", row_text)]
    if (nrow(unclear) > 0) {
      add_issue("GSE279446", "severe", "Control-like row is not clearly WT/control-like.", nrow(unclear), unclear$sample_key)
    }
    if (nrow(primary) != 8L) {
      add_issue("GSE279446", "warning", "Control-like primary count is not the expected 8 rows.", nrow(primary), primary$sample_key)
    }
  }

  gse273724 <- meta_dt[dataset_folder == "GSE273724"]
  if (nrow(gse273724) > 0) {
    primary <- gse273724[control_like_clock_candidate == TRUE]
    if (nrow(primary) != 4L) {
      add_issue("GSE273724", "warning", "Control-like primary count is not the expected 4 rows.", nrow(primary), primary$sample_key)
    }
    extra <- primary[!grepl("young.*7.*8|young_7_8w", text_any(primary, c("sample_key", "sample_id", "sample_title", "treatment")), ignore.case = TRUE)]
    if (nrow(extra) > 0) {
      add_issue("GSE273724", "warning", "Additional non-young_7_8w rows are marked control-like.", nrow(extra), extra$sample_key)
    }
  }

  out <- if (length(issues) > 0) {
    rbindlist(issues, fill = TRUE)
  } else {
    data.table(dataset_folder = character(), severity = character(), issue = character(), n_rows = integer(), sample_keys = character())
  }
  if (nrow(out[severity == "severe"]) > 0) {
    stop("Severe metadata issues detected before training.", call. = FALSE)
  }
  invisible(out)
}

collapse_duplicate_genes <- function(mat) {
  gene_ids <- trimws(rownames(mat))
  keep <- !is.na(gene_ids) & gene_ids != ""
  mat <- mat[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]
  if (mean(grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", gene_ids), na.rm = TRUE) > 0.5) {
    gene_ids <- sub("\\.[0-9]+$", "", gene_ids)
  }
  rownames(mat) <- gene_ids
  if (!anyDuplicated(gene_ids)) return(list(matrix = mat, duplicated = data.table()))

  row_mean <- rowMeans(mat, na.rm = TRUE)
  keep_idx <- unlist(tapply(seq_along(gene_ids), gene_ids, function(idx) {
    idx[which.max(row_mean[idx])]
  }), use.names = FALSE)
  duplicated_audit <- data.table(
    gene_id = gene_ids[duplicated(gene_ids) | duplicated(gene_ids, fromLast = TRUE)]
  )[, .N, by = gene_id][order(-N)]
  list(matrix = mat[sort(keep_idx), , drop = FALSE], duplicated = duplicated_audit)
}

logcpm_from_counts <- function(count_mat) {
  lib_size <- colSums(count_mat, na.rm = TRUE)
  if (any(!is.finite(lib_size) | lib_size <= 0)) {
    stop("Cannot compute logCPM: one or more mapped count columns has non-positive library size.")
  }
  log1p(t(t(count_mat) / lib_size * 1e6))
}

validate_bulk_group_mapping_columns <- function(mapping) {
  required <- c(
    "dataset_folder", "metadata_sample_key", "metadata_age_week",
    "metadata_control_like_clock_candidate", "expression_column",
    "sample_unit_type", "aggregation_id", "aggregation_method",
    "output_sample_key", "include_in_primary_clock_TRUE_FALSE",
    "manual_approval_TRUE_FALSE", "manual_note"
  )
  missing <- setdiff(required, names(mapping))
  if (length(missing) > 0) {
    stop(
      "Approved bulk group mapping is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

align_dataset_matrix_with_bulk_group_mapping <- function(meta_ds, ds, logcpm_path, counts_path, expr_cols) {
  if (nrow(bulk_group_mapping) == 0) {
    stop("No approved bulk group mapping is available for ", ds, call. = FALSE)
  }
  validate_bulk_group_mapping_columns(bulk_group_mapping)
  map_ds <- copy(bulk_group_mapping[dataset_folder == ds])
  if (nrow(map_ds) == 0) {
    stop("No approved bulk group mapping rows are available for ", ds, call. = FALSE)
  }

  map_ds[, include_in_primary_clock_TRUE_FALSE := truthy(include_in_primary_clock_TRUE_FALSE)]
  map_ds[, manual_approval_TRUE_FALSE := truthy(manual_approval_TRUE_FALSE)]
  map_ds <- map_ds[include_in_primary_clock_TRUE_FALSE == TRUE & manual_approval_TRUE_FALSE == TRUE]
  if (nrow(map_ds) == 0) {
    stop("Approved bulk group mapping for ", ds, " has no included, manually approved rows.", call. = FALSE)
  }

  allowed_sample_units <- c("biological_replicate", "technical_run", "group_mean")
  allowed_methods <- c("none_expand_replicates", "sum_raw_counts", "mean_logcpm")
  if (any(!map_ds$sample_unit_type %chin% allowed_sample_units)) {
    stop("Unsafe sample_unit_type in approved bulk group mapping for ", ds, call. = FALSE)
  }
  if (any(!map_ds$aggregation_method %chin% allowed_methods)) {
    stop("Unsafe aggregation_method in approved bulk group mapping for ", ds, call. = FALSE)
  }
  if (anyNA(map_ds$output_sample_key) || any(map_ds$output_sample_key == "")) {
    stop("Approved bulk group mapping for ", ds, " has missing output_sample_key.", call. = FALSE)
  }
  if (anyNA(map_ds$aggregation_id) || any(map_ds$aggregation_id == "")) {
    stop("Approved bulk group mapping for ", ds, " has missing aggregation_id.", call. = FALSE)
  }
  if (any(map_ds$sample_unit_type == "technical_run" & map_ds$aggregation_method != "sum_raw_counts")) {
    stop(
      "Approved bulk group mapping for ", ds,
      " marks technical_run rows without aggregation_method == sum_raw_counts.",
      call. = FALSE
    )
  }
  expected_cols <- unique(trimws(as.character(map_ds$expression_column)))
  missing_cols <- setdiff(expected_cols, expr_cols)
  if (length(missing_cols) > 0) {
    stop(
      "Approved bulk group mapping for ", ds,
      " references expression columns not present in matrix: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  map_ds[, metadata_age_week_numeric := suppressWarnings(as.numeric(metadata_age_week))]
  if (any(!is.finite(map_ds$metadata_age_week_numeric))) {
    stop("Approved bulk group mapping for ", ds, " has missing or non-numeric metadata_age_week.", call. = FALSE)
  }

  meta_idx <- match(map_ds$metadata_sample_key, meta_ds$sample_key)
  if (anyNA(meta_idx)) {
    stop("Approved bulk group mapping for ", ds, " references metadata_sample_key not in metadata.", call. = FALSE)
  }
  map_ds[, metadata_control_like_clock_candidate := truthy(metadata_control_like_clock_candidate)]
  if (any(map_ds$metadata_control_like_clock_candidate != TRUE)) {
    stop("Approved bulk group mapping for ", ds, " includes non-control-like rows.", call. = FALSE)
  }
  if (any(meta_ds$control_like_clock_candidate[meta_idx] != TRUE)) {
    stop("Approved bulk group mapping for ", ds, " includes rows not marked control_like_clock_candidate in metadata.", call. = FALSE)
  }

  map_ds[, output_per_aggregation := output_sample_key[1], by = aggregation_id]
  output_by_agg <- unique(map_ds[, .(aggregation_id, output_sample_key = output_per_aggregation)])
  if (anyDuplicated(output_by_agg$output_sample_key)) {
    stop("Approved bulk group mapping for ", ds, " has duplicated output sample keys after aggregation.", call. = FALSE)
  }
  agg_output_counts <- map_ds[, .(n_output_sample_keys = uniqueN(output_sample_key)), by = aggregation_id]
  if (nrow(agg_output_counts[n_output_sample_keys > 1]) > 0) {
    stop("Approved bulk group mapping for ", ds, " has aggregation_id values mapping to multiple output_sample_key values.", call. = FALSE)
  }
  technical_groups <- map_ds[sample_unit_type == "technical_run",
    .(
      n_metadata_sample_keys = uniqueN(metadata_sample_key),
      n_ages = uniqueN(metadata_age_week_numeric)
    ),
    by = aggregation_id
  ]
  if (nrow(technical_groups[n_metadata_sample_keys > 1 | n_ages > 1]) > 0) {
    stop("Approved bulk group mapping for ", ds, " has technical_run aggregation_id groups spanning multiple metadata samples or ages.", call. = FALSE)
  }

  methods <- unique(map_ds$aggregation_method)
  if (length(methods) > 1) {
    stop("Approved bulk group mapping for ", ds, " mixes aggregation methods; split the mapping by dataset/sample type.", call. = FALSE)
  }
  method <- methods[[1]]
  manual_notes <- if ("manual_note" %in% names(map_ds)) as.character(map_ds$manual_note) else ""
  if (method == "mean_logcpm" && !any(grepl("approved_group_mean_logcpm", manual_notes, ignore.case = TRUE))) {
    stop(
      "Approved bulk group mapping for ", ds,
      " uses mean_logcpm but manual_note does not contain approved_group_mean_logcpm.",
      call. = FALSE
    )
  }
  input_scale <- if (method == "sum_raw_counts") "raw_counts" else "logcpm"
  source_path <- if (method == "sum_raw_counts") counts_path else logcpm_path
  source_mat <- read_expression_matrix(source_path)
  source_clean <- collapse_duplicate_genes(source_mat)
  source_mat <- source_clean$matrix

  if (method == "none_expand_replicates") {
    out <- source_mat[, map_ds$expression_column, drop = FALSE]
    colnames(out) <- map_ds$output_sample_key
    out_meta_source <- meta_ds[match(map_ds$metadata_sample_key, sample_key)]
    out_meta <- copy(out_meta_source)
    out_meta[, sample_key := map_ds$output_sample_key]
    out_meta[, sample_id := map_ds$output_sample_key]
    out_meta[, age_week := map_ds$metadata_age_week_numeric]
    group_mean_warning <- FALSE
  } else {
    grouped <- split(seq_len(nrow(map_ds)), map_ds$aggregation_id)
    out_list <- lapply(grouped, function(idx) {
      cols <- map_ds$expression_column[idx]
      if (method == "sum_raw_counts") {
        rowSums(source_mat[, cols, drop = FALSE], na.rm = TRUE)
      } else {
        rowMeans(source_mat[, cols, drop = FALSE], na.rm = TRUE)
      }
    })
    out <- do.call(cbind, out_list)
    output_keys <- vapply(grouped, function(idx) map_ds$output_sample_key[idx][[1]], character(1))
    colnames(out) <- output_keys
    if (method == "sum_raw_counts") out <- logcpm_from_counts(out)
    meta_keys <- vapply(grouped, function(idx) map_ds$metadata_sample_key[idx][[1]], character(1))
    out_meta <- copy(meta_ds[match(meta_keys, sample_key)])
    out_meta[, sample_key := output_keys]
    out_meta[, sample_id := output_keys]
    out_meta[, age_week := vapply(grouped, function(idx) map_ds$metadata_age_week_numeric[idx][[1]], numeric(1))]
    group_mean_warning <- method == "mean_logcpm"
  }

  if (anyDuplicated(colnames(out))) {
    stop("Mapped output sample keys are not unique for ", ds, ".", call. = FALSE)
  }
  gene_count <- nrow(out)
  common_universe_path <- file.path(results_dir, "integrated_bulk_pseudobulk_technical_gene_universe.tsv")
  overlap_common <- NA_integer_
  if (file.exists(common_universe_path)) {
    common_dt <- fread(common_universe_path)
    common_col <- if ("gene_id" %in% names(common_dt)) "gene_id" else names(common_dt)[1]
    overlap_common <- length(intersect(rownames(out), as.character(common_dt[[common_col]])))
  }
  if (gene_count < 5000) {
    stop("Mapped dataset ", ds, " has too few genes after cleaning: ", gene_count, call. = FALSE)
  }
  if (is.finite(overlap_common) && overlap_common < 5000) {
    stop("Mapped dataset ", ds, " has too few genes overlapping the integrated gene universe: ", overlap_common, call. = FALSE)
  }
  if (ds == "GSE279446" && is.finite(overlap_common) && overlap_common < 5000) {
    stop("GSE279446 is excluded from primary training because gene overlap remains poor.", call. = FALSE)
  }

  selected_gene_paths <- list.files(selected_gene_dir, pattern = "selected_genes\\.tsv$", full.names = TRUE)
  selected_gene_coverage <- NA_real_
  if (length(selected_gene_paths) > 0) {
    selected <- unique(unlist(lapply(selected_gene_paths, function(path) {
      dt <- fread(path)
      col <- if ("gene_id" %in% names(dt)) "gene_id" else names(dt)[1]
      as.character(dt[[col]])
    }), use.names = FALSE))
    selected_gene_coverage <- length(intersect(rownames(out), selected)) / length(selected)
    if (is.finite(selected_gene_coverage) && selected_gene_coverage < 0.80) {
      bulk_group_warnings <<- c(
        bulk_group_warnings,
        paste0(ds, ": selected-gene coverage below 0.80 after mapping (", round(selected_gene_coverage, 3), ").")
      )
    }
  }
  if (group_mean_warning) {
    bulk_group_warnings <<- c(
      bulk_group_warnings,
      paste0(ds, ": mean_logcpm creates group-level samples; consider dataset-balanced sensitivity analysis.")
    )
  }
  if (ncol(out) > 10) {
    bulk_group_warnings <<- c(
      bulk_group_warnings,
      paste0(ds, ": mapping adds ", ncol(out), " samples; consider later dataset-balanced weighting sensitivity.")
    )
  }

  bulk_group_added_samples[[length(bulk_group_added_samples) + 1L]] <<- out_meta
  bulk_group_mapping_used_datasets <<- unique(c(bulk_group_mapping_used_datasets, ds))
  bulk_group_expression_column_map[[length(bulk_group_expression_column_map) + 1L]] <<- copy(map_ds[, .(
    dataset_folder, metadata_sample_key, expression_column, sample_unit_type,
    aggregation_id, aggregation_method, output_sample_key
  )])
  bulk_group_dataset_summary[[length(bulk_group_dataset_summary) + 1L]] <<- data.table(
    dataset_folder = ds,
    n_output_samples = ncol(out),
    n_mapping_rows = nrow(map_ds),
    aggregation_method = method,
    expression_scale_decision = input_scale,
    group_mean_warning = group_mean_warning
  )
  bulk_group_gene_qc[[length(bulk_group_gene_qc) + 1L]] <<- data.table(
    dataset_folder = ds,
    n_genes_after_cleaning = gene_count,
    n_duplicate_gene_ids_removed = if (nrow(source_clean$duplicated) > 0) sum(source_clean$duplicated$N - 1L) else 0L,
    overlap_with_integrated_gene_universe = overlap_common,
    selected_gene_coverage = selected_gene_coverage
  )

  list(
    matrix = out,
    metadata = out_meta,
    alignment = data.table(
      sample_key = map_ds$output_sample_key,
      dataset_folder = ds,
      expression_column = map_ds$expression_column,
      matched_by = paste0("approved_bulk_group_mapping:", method)
    )
  )
}

metadata_looks_samplelevel <- function(meta_ds, expr_cols, modality) {
  candidates <- if (modality == "bulk_rnaseq") {
    list(sample_id = meta_ds$sample_id, accession = meta_ds$accession, sample_key = meta_ds$sample_key)
  } else {
    list(sample_id = meta_ds$sample_id, sample_key = meta_ds$sample_key)
  }
  any(vapply(candidates, function(x) sum(as.character(x) %chin% expr_cols, na.rm = TRUE), integer(1)) > 0L) ||
    nrow(meta_ds) == length(expr_cols)
}

check_manual_mapping_conflict_with_exact <- function(ds, exact_alignment) {
  if (bulk_group_mapping_mode != "force" || nrow(bulk_group_mapping_loaded) == 0) return(invisible(NULL))
  if (!"dataset_folder" %in% names(bulk_group_mapping_loaded) || !ds %chin% bulk_group_mapping_loaded$dataset_folder) {
    return(invisible(NULL))
  }
  validate_bulk_group_mapping_columns(bulk_group_mapping_loaded)
  map_ds <- copy(bulk_group_mapping_loaded[dataset_folder == ds])
  map_ds[, include_in_primary_clock_TRUE_FALSE := truthy(include_in_primary_clock_TRUE_FALSE)]
  map_ds[, manual_approval_TRUE_FALSE := truthy(manual_approval_TRUE_FALSE)]
  map_ds <- map_ds[include_in_primary_clock_TRUE_FALSE == TRUE & manual_approval_TRUE_FALSE == TRUE]
  if (nrow(map_ds) == 0) return(invisible(NULL))
  exact_ds <- exact_alignment[dataset_folder == ds]
  exact_idx <- match(map_ds$metadata_sample_key, exact_ds$sample_key)
  missing <- map_ds[is.na(exact_idx), metadata_sample_key]
  if (length(missing) > 0) {
    stop(
      "BULK_GROUP_MAPPING_MODE=force conflict for ", ds,
      ": approved mapping references metadata_sample_key values absent from exact alignment: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  exact_matched <- exact_ds$expression_column[exact_idx]
  conflict <- map_ds$expression_column != exact_matched | map_ds$output_sample_key != map_ds$metadata_sample_key
  if (any(conflict)) {
    stop(
      "BULK_GROUP_MAPPING_MODE=force conflict for ", ds,
      ": exact sample-level alignment disagrees with approved manual mapping.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

align_dataset_matrix <- function(mat, meta_ds, modality, allow_position_fallback = TRUE) {
  expr_cols <- colnames(mat)
  candidates <- if (modality == "bulk_rnaseq") {
    list(sample_id = meta_ds$sample_id, accession = meta_ds$accession, sample_key = meta_ds$sample_key)
  } else {
    list(sample_id = meta_ds$sample_id, sample_key = meta_ds$sample_key)
  }

  matched_col <- rep(NA_character_, nrow(meta_ds))
  matched_by <- rep(NA_character_, nrow(meta_ds))
  for (candidate_name in names(candidates)) {
    idx <- match(as.character(candidates[[candidate_name]]), expr_cols)
    use <- is.na(matched_col) & !is.na(idx)
    matched_col[use] <- expr_cols[idx[use]]
    matched_by[use] <- candidate_name
  }

  if (any(is.na(matched_col))) {
    if (allow_position_fallback && ncol(mat) == nrow(meta_ds)) {
      missing_idx <- which(is.na(matched_col))
      matched_col[missing_idx] <- expr_cols[missing_idx]
      matched_by[missing_idx] <- "dataset_position"
    } else {
      missing_samples <- meta_ds[is.na(matched_col), .(sample_id, accession, sample_key)]
      stop(
        "Could not align all samples for ", unique(meta_ds$dataset_folder),
        ". Missing: ", paste(missing_samples$sample_key, collapse = ", ")
      )
    }
  }

  out <- mat[, matched_col, drop = FALSE]
  colnames(out) <- meta_ds$sample_key
  list(matrix = out, alignment = data.table(
    sample_key = meta_ds$sample_key,
    dataset_folder = meta_ds$dataset_folder,
    expression_column = matched_col,
    matched_by = matched_by
  ))
}

metadata_for_keys <- function(meta_dt, keys) {
  idx <- match(keys, meta_dt$sample_key)
  if (anyNA(idx)) {
    missing <- keys[is.na(idx)]
    stop("Missing metadata rows for sample_key values: ", paste(missing, collapse = ", "))
  }
  meta_dt[idx]
}

load_modality_expression <- function(meta, modality_name) {
  modality_filter <- meta[["modality"]] == modality_name
  meta_mod_all <- copy(meta[modality_filter])
  meta_mod_primary <- copy(meta_mod_all[control_like_clock_candidate == TRUE])
  if (nrow(meta_mod_primary) == 0) stop("No control-like samples for modality: ", modality_name)
  pieces <- list()
  metadata_pieces <- list()
  alignment <- list()
  alignment_issues <- list()
  for (ds in unique(meta_mod_primary$dataset_folder)) {
    meta_ds_all <- meta_mod_all[dataset_folder == ds]
    logcpm_path <- unique(meta_ds_all$expression_logcpm_path)[1]
    counts_path <- unique(meta_ds_all$expression_counts_path)[1]
    mat <- read_expression_matrix(logcpm_path)
    samplelevel_candidate <- metadata_looks_samplelevel(meta_ds_all, colnames(mat), modality_name)
    aligned <- tryCatch(
      align_dataset_matrix(
        mat,
        meta_ds_all,
        modality_name,
        allow_position_fallback = modality_name != "bulk_rnaseq"
      ),
      error = function(e) {
        if (
          modality_name == "bulk_rnaseq" &&
            bulk_group_mapping_mode != "disabled" &&
            (bulk_group_mapping_mode == "force" || !samplelevel_candidate)
        ) {
          tryCatch(
            align_dataset_matrix_with_bulk_group_mapping(
              meta_ds = meta_ds_all,
              ds = ds,
              logcpm_path = logcpm_path,
              counts_path = counts_path,
              expr_cols = colnames(mat)
            ),
            error = function(e2) {
              alignment_issues[[length(alignment_issues) + 1L]] <<- data.table(
                modality = modality_name,
                dataset_folder = ds,
                issue = paste(conditionMessage(e), "| bulk_group_mapping:", conditionMessage(e2))
              )
              NULL
            }
          )
        } else {
          issue <- conditionMessage(e)
          if (modality_name == "bulk_rnaseq" && bulk_group_mapping_mode == "disabled") {
            issue <- paste(issue, "| bulk_group_mapping disabled")
          } else if (modality_name == "bulk_rnaseq" && samplelevel_candidate) {
            issue <- paste(issue, "| metadata appears sample-level; manual group mapping not used in auto mode")
          }
          alignment_issues[[length(alignment_issues) + 1L]] <<- data.table(
            modality = modality_name,
            dataset_folder = ds,
            issue = issue
          )
          NULL
        }
      }
    )
    if (is.null(aligned)) next
    if (modality_name == "bulk_rnaseq" && is.null(aligned$metadata)) {
      exact_samplelevel_matching_used_for <<- unique(c(exact_samplelevel_matching_used_for, ds))
      check_manual_mapping_conflict_with_exact(ds, aligned$alignment)
    }
    if (!is.null(aligned$metadata)) {
      pieces[[ds]] <- aligned$matrix
      metadata_pieces[[ds]] <- aligned$metadata
    } else {
      primary_keys <- meta_mod_primary[dataset_folder == ds, sample_key]
      pieces[[ds]] <- aligned$matrix[, primary_keys, drop = FALSE]
      metadata_pieces[[ds]] <- metadata_for_keys(meta_mod_primary, primary_keys)
    }
    alignment[[ds]] <- aligned$alignment
  }
  if (length(pieces) == 0) stop("No alignable datasets for modality: ", modality_name)
  common_genes <- Reduce(intersect, lapply(pieces, rownames))
  common_genes <- sort(common_genes)
  if (length(common_genes) < 100) {
    stop("Too few common genes for ", modality_name, ": ", length(common_genes))
  }
  expr <- do.call(cbind, lapply(pieces, function(x) x[common_genes, , drop = FALSE]))
  meta_mod <- rbindlist(metadata_pieces, fill = TRUE)
  meta_mod <- metadata_for_keys(meta_mod, colnames(expr))
  rownames(meta_mod) <- meta_mod$sample_key
  setkey(meta_mod, sample_key)
  list(
    expression = expr,
    metadata = meta_mod,
    alignment = rbindlist(alignment),
    alignment_issues = if (length(alignment_issues) > 0) {
      rbindlist(alignment_issues, fill = TRUE)
    } else {
      data.table(modality = character(), dataset_folder = character(), issue = character())
    },
    common_genes = common_genes
  )
}

add_integrated_covariates <- function(meta_dt, cv_group_col = "dataset_folder") {
  meta_dt <- copy(meta_dt)
  if (!"dataset_subcohort" %in% names(meta_dt)) {
    meta_dt[, dataset_subcohort := dataset_folder]
  }
  meta_dt[is.na(dataset_subcohort) | dataset_subcohort == "", dataset_subcohort := dataset_folder]
  meta_dt[, source_batch := interaction(modality, dataset_folder, drop = TRUE)]
  meta_dt[, source_batch_for_feature_selection := interaction(
    modality,
    dataset_folder,
    dataset_subcohort,
    drop = TRUE
  )]
  meta_dt[, log_n_cells_sc := fifelse(
    modality == "scrna_pseudobulk",
    log1p(fifelse(is.na(n_cells_used), 0, as.numeric(n_cells_used))),
    0
  )]
  meta_dt[, cv_group := as.character(get(cv_group_col))]
  rownames(meta_dt) <- meta_dt$sample_key
  setkey(meta_dt, sample_key)
  meta_dt
}

load_integrated_expression <- function(meta) {
  bulk_loaded <- load_modality_expression(meta, "bulk_rnaseq")
  scrna_loaded <- load_modality_expression(meta, "scrna_pseudobulk")
  common_genes <- sort(intersect(rownames(bulk_loaded$expression), rownames(scrna_loaded$expression)))
  if (length(common_genes) < 100) {
    stop("Too few common integrated genes: ", length(common_genes))
  }

  expr <- cbind(
    bulk_loaded$expression[common_genes, , drop = FALSE],
    scrna_loaded$expression[common_genes, , drop = FALSE]
  )
  if (anyDuplicated(colnames(expr))) {
    stop("Integrated sample_key values are not unique.")
  }

  meta_integrated_source <- rbindlist(list(bulk_loaded$metadata, scrna_loaded$metadata), fill = TRUE)
  meta_integrated <- metadata_for_keys(meta_integrated_source, colnames(expr))
  meta_integrated <- add_integrated_covariates(meta_integrated, cv_group_col = "dataset_folder")
  if (
    nrow(bulk_group_mapping) > 0 &&
      "dataset_folder" %in% names(bulk_group_mapping) &&
      "GSE279573" %chin% bulk_group_mapping$dataset_folder &&
      meta_integrated[dataset_folder == "GSE279573" & split_role == "validation", .N] == 0
  ) {
    stop("GSE279573 is expected in validation but has zero model-ready samples after mapping.", call. = FALSE)
  }

  list(
    expression = expr,
    metadata = meta_integrated,
    alignment = rbindlist(list(bulk_loaded$alignment, scrna_loaded$alignment), fill = TRUE),
    alignment_issues = rbindlist(list(bulk_loaded$alignment_issues, scrna_loaded$alignment_issues), fill = TRUE),
    common_genes = common_genes,
    modality_common_genes = data.table(
      modality = c("bulk_rnaseq", "scrna_pseudobulk", training_scope),
      n_common_genes = c(
        length(bulk_loaded$common_genes),
        length(scrna_loaded$common_genes),
        length(common_genes)
      )
    )
  )
}

select_hvg <- function(expr_train, n_genes = 2000L) {
  vars <- apply(expr_train, 1, stats::var, na.rm = TRUE)
  vars[!is.finite(vars)] <- NA_real_
  head(names(sort(vars, decreasing = TRUE, na.last = NA)), min(n_genes, sum(is.finite(vars))))
}

selector_label <- function(selector) {
  if (selector == "hvg_baseline") return("HVG baseline")
  if (selector == "version_b_age_dataset_penalized") return("Version B age/source-batch penalized")
  selector
}

selector_file_stub <- function(selector) {
  if (selector == "hvg_baseline") return("hvg")
  if (selector == "version_b_age_dataset_penalized") return("versionB")
  safe_name(selector)
}

plot_top30_versionb_genes <- function(top_dt, family) {
  if (nrow(top_dt) == 0 || !"feature_score" %in% names(top_dt)) return(invisible(NULL))
  top_dt <- copy(top_dt[is.finite(feature_score)])
  if (nrow(top_dt) == 0) return(invisible(NULL))
  top_dt <- top_dt[order(feature_score)]
  prefix <- file.path(plot_dir, paste0("13_top30_versionB_genes_integrated_", family))
  plot_one <- function() {
    old_mar <- par("mar")
    on.exit(par(mar = old_mar), add = TRUE)
    par(mar = c(5, 12, 4, 2))
    barplot(
      height = top_dt$feature_score,
      names.arg = top_dt$gene_id,
      horiz = TRUE,
      las = 1,
      col = "#4C78A8",
      border = NA,
      xlab = "Version B feature score",
      main = paste("Top 30 integrated Version B genes:", family)
    )
  }
  grDevices::png(paste0(prefix, ".png"), width = 2400, height = 1800, res = 220)
  plot_one()
  grDevices::dev.off()
  grDevices::pdf(paste0(prefix, ".pdf"), width = 11, height = 8)
  plot_one()
  grDevices::dev.off()
  invisible(NULL)
}

prepare_version_genes <- function(selector, expr_train, meta_train, out_dir) {
  if (selector == "hvg_baseline") {
    genes <- select_hvg(expr_train, hvg_n_genes)
    table <- data.table(gene_id = genes, rank = seq_along(genes), feature_score = NA_real_)
    return(list(genes = genes, table = table, score_table = table))
  }

  meta_fs <- as.data.frame(meta_train)
  rownames(meta_fs) <- meta_fs$sample_key
  covariates <- character(0)
  dataset_col <- if ("source_batch_for_feature_selection" %in% names(meta_fs)) {
    "source_batch_for_feature_selection"
  } else {
    "source_batch"
  }

  fs <- select_age_dataset_penalized_genes(
    logcpm_train = expr_train,
    meta_train = meta_fs,
    dataset_col = dataset_col,
    covariates = covariates,
    n_genes = version_b_n_genes,
    dataset_penalty = age_dataset_penalized_defaults$dataset_penalty,
    max_dataset_partial_r2 = age_dataset_penalized_defaults$max_dataset_partial_r2,
    balance_direction = age_dataset_penalized_defaults$balance_direction
  )
  fwrite(fs$feature_score_table, file.path(out_dir, "feature_score_table.tsv"), sep = "\t")
  fwrite(fs$selected_gene_table, file.path(out_dir, "selected_gene_table.tsv"), sep = "\t")
  fwrite(head(fs$selected_gene_table[order(-feature_score)], 30), file.path(out_dir, "selected_gene_top30.tsv"), sep = "\t")
  list(genes = fs$selected_genes, table = fs$selected_gene_table, score_table = fs$feature_score_table)
}

scale_by_training <- function(expr, genes, train_keys, apply_keys) {
  x_train <- t(expr[genes, train_keys, drop = FALSE])
  center <- colMeans(x_train, na.rm = TRUE)
  scale <- apply(x_train, 2, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  x_apply <- t(expr[genes, apply_keys, drop = FALSE])
  x_scaled <- sweep(x_apply, 2, center, "-")
  x_scaled <- sweep(x_scaled, 2, scale, "/")
  x_scaled[!is.finite(x_scaled)] <- 0
  list(x = x_scaled, center = center, scale = scale)
}

fit_alpha_grid <- function(x, y, foldid = NULL) {
  if (nrow(x) < 6 || length(unique(y)) < 3) return(NULL)
  records <- list()
  best <- NULL
  nfolds <- if (is.null(foldid)) min(5L, nrow(x)) else length(unique(foldid))
  for (alpha in alpha_grid) {
    fit <- tryCatch(
      cv.glmnet(
        x = x,
        y = y,
        alpha = alpha,
        family = "gaussian",
        foldid = foldid,
        nfolds = nfolds,
        type.measure = "mae",
        standardize = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    rec <- data.table(alpha = alpha, lambda = fit$lambda.min, cv_mae_week = min(fit$cvm, na.rm = TRUE))
    records[[length(records) + 1L]] <- rec
    if (is.null(best) || rec$cv_mae_week < best$cv_mae_week) best <- rec
  }
  if (is.null(best)) return(NULL)
  list(best = best, cv_table = rbindlist(records))
}

calibrate_prediction <- function(y_train, raw_train, raw_apply, calibration_sample_keys) {
  calibration_df <- data.frame(
    age_week = as.numeric(y_train),
    raw_pred_age_week = as.numeric(raw_train)
  )
  if (length(unique(calibration_df$raw_pred_age_week)) < 2 || length(unique(calibration_df$age_week)) < 2) {
    stop("Calibration requires variation in observed age_week and raw_pred_age_week.", call. = FALSE)
  }
  calib_fit <- lm(age_week ~ raw_pred_age_week, data = calibration_df)
  intercept <- unname(coef(calib_fit)[1])
  slope <- unname(coef(calib_fit)[2])
  calibration <- list(
    calibration_type = "linear_lm",
    formula = "age_week ~ raw_pred_age_week",
    intercept = intercept,
    slope = slope,
    raw_prediction_column = "raw_pred_age_week",
    calibrated_prediction_column = "predicted_age_week",
    calibration_sample_keys = as.character(calibration_sample_keys),
    calibration_n = length(calibration_sample_keys),
    calibration_age_range = range(calibration_df$age_week, na.rm = TRUE),
    calibration_raw_prediction_range = range(calibration_df$raw_pred_age_week, na.rm = TRUE),
    clip_predictions = FALSE
  )
  list(
    pred = intercept + slope * as.numeric(raw_apply),
    method = "linear_lm_training_only",
    calibration = calibration
  )
}

prediction_table <- function(meta_apply, raw_pred, pred, model_family, selector, fit_stage,
                             calibration_method) {
  data.table(
    training_scope = training_scope,
    model_family = model_family,
    selector = selector,
    model_version = selector,
    fit_stage = fit_stage,
    calibration_method = calibration_method,
    sample_key = meta_apply$sample_key,
    sample_id = meta_apply$sample_id,
    modality = meta_apply$modality,
    dataset_folder = meta_apply$dataset_folder,
    dataset_subcohort = if ("dataset_subcohort" %in% names(meta_apply)) meta_apply$dataset_subcohort else meta_apply$dataset_folder,
    split_role = meta_apply$split_role,
    biological_age_bin = meta_apply$bio_age_bin,
    bio_age_bin = meta_apply$bio_age_bin,
    split_age_stratum = meta_apply$split_age_stratum,
    age_week = as.numeric(meta_apply$age_week),
    raw_predicted_age_week = raw_pred,
    predicted_age_week = pred,
    residual_week = pred - as.numeric(meta_apply$age_week),
    abs_error_week = abs(pred - as.numeric(meta_apply$age_week))
  )
}

metric_one <- function(dt, by_cols = character()) {
  if (length(by_cols) == 0) {
    groups <- list(overall = dt)
  } else {
    group_id <- do.call(paste, c(as.data.frame(dt[, ..by_cols]), sep = "\r"))
    groups <- split(dt, group_id)
  }
  out <- lapply(groups, function(x) {
    y <- x$age_week
    p <- x$predicted_age_week
    data.table(
      n_samples = nrow(x),
      mae_week = mean(abs(p - y), na.rm = TRUE),
      rmse_week = sqrt(mean((p - y)^2, na.rm = TRUE)),
      median_abs_error_week = stats::median(abs(p - y), na.rm = TRUE),
      mean_residual_week = mean(p - y, na.rm = TRUE),
      residual_sd_week = stats::sd(p - y, na.rm = TRUE),
      pearson = suppressWarnings(stats::cor(y, p, method = "pearson", use = "complete.obs")),
      spearman = suppressWarnings(stats::cor(y, p, method = "spearman", use = "complete.obs"))
    )
  })
  res <- rbindlist(out, idcol = "group_id")
  if (length(by_cols) > 0) {
    keys <- tstrsplit(res$group_id, "\r", fixed = TRUE)
    for (i in seq_along(by_cols)) res[, (by_cols[i]) := keys[[i]]]
    res[, group_id := NULL]
    setcolorder(res, c(by_cols, setdiff(names(res), by_cols)))
  } else {
    res[, group_id := NULL]
  }
  res
}

summarise_predictions <- function(pred) {
  keys <- c("training_scope", "model_family", "selector", "fit_stage")
  out <- list(
    c(keys),
    c(keys, "split_role"),
    c(keys, "split_role", "modality"),
    c(keys, "split_role", "biological_age_bin"),
    c(keys, "split_role", "split_age_stratum"),
    c(keys, "split_role", "dataset_folder"),
    c(keys, "split_role", "dataset_subcohort")
  )
  rbindlist(lapply(out, function(cols) metric_one(pred, cols)), fill = TRUE)
}

fit_elasticnet_clock <- function(expr, meta, train_keys, apply_keys, genes, selector, out_dir,
                                 fit_stage, fold_col = "cv_group") {
  meta_train <- metadata_for_keys(meta, train_keys)
  meta_apply <- metadata_for_keys(meta, apply_keys)
  scaled_train <- scale_by_training(expr, genes, train_keys, train_keys)
  x_train <- scaled_train$x
  keep <- apply(x_train, 2, stats::sd) > 0
  x_train <- x_train[, keep, drop = FALSE]
  genes <- colnames(x_train)
  if (length(genes) < 2) stop("Too few nonzero-variance genes for ElasticNet.")
  fold_source <- if (fold_col %in% names(meta_train)) meta_train[[fold_col]] else meta_train$dataset_folder
  foldid <- as.integer(factor(fold_source))
  tuned <- fit_alpha_grid(x_train, meta_train$age_week, foldid = foldid)
  if (is.null(tuned)) stop("ElasticNet tuning failed.")
  final <- glmnet(
    x = x_train,
    y = meta_train$age_week,
    alpha = tuned$best$alpha,
    lambda = tuned$best$lambda,
    family = "gaussian",
    standardize = FALSE
  )
  x_apply <- scale_by_training(expr, genes, train_keys, apply_keys)$x
  raw_train <- as.numeric(predict(final, newx = x_train))
  raw_apply <- as.numeric(predict(final, newx = x_apply))
  cal <- calibrate_prediction(meta_train$age_week, raw_train, raw_apply, meta_train$sample_key)
  coef_dt <- data.table(gene_id = genes, coefficient = as.numeric(coef(final))[-1])
  nonzero_coef_dt <- coef_dt[coefficient != 0]
  model_component <- list(
    glmnet_fit = final,
    alpha = as.numeric(tuned$best$alpha),
    lambda = as.numeric(tuned$best$lambda),
    glmnet_intercept = unname(as.numeric(coef(final))[1]),
    coefficient_table = coef_dt,
    nonzero_coefficient_table = nonzero_coef_dt,
    raw_prediction_function_inputs = list(
      feature_space = "training_scaled_selected_genes",
      standardize = FALSE
    )
  )
  fwrite(tuned$cv_table, file.path(out_dir, paste0(fit_stage, "_elasticnet_cv_grid.tsv")), sep = "\t")
  fwrite(coef_dt, file.path(out_dir, paste0(fit_stage, "_elasticnet_coefficients.tsv")), sep = "\t")
  saveRDS(
    list(
      model = final,
      genes = genes,
      center = scaled_train$center[genes],
      scale = scaled_train$scale[genes],
      training_scope = training_scope,
      calibration = cal$calibration,
      alpha = as.numeric(tuned$best$alpha),
      lambda = as.numeric(tuned$best$lambda),
      coefficient_table = coef_dt,
      nonzero_coefficient_table = nonzero_coef_dt
    ),
    file.path(out_dir, paste0(fit_stage, "_elasticnet_model.rds"))
  )
  list(
    predictions = prediction_table(meta_apply, raw_apply, cal$pred, "elasticnet", selector, fit_stage, cal$method),
    fit_artifact = list(
      model_family = "elasticnet",
      selector = selector,
      fit_stage = fit_stage,
      selected_genes = genes,
      gene_center = scaled_train$center[genes],
      gene_scale = scaled_train$scale[genes],
      model = model_component,
      calibration = cal$calibration,
      training_sample_keys = as.character(train_keys),
      apply_sample_keys = as.character(apply_keys),
      reference_metadata = copy(meta_train)
    )
  )
}

fit_pca_glmnet_clock <- function(expr, meta, train_keys, apply_keys, genes, selector, out_dir,
                                 fit_stage, fold_col = "cv_group") {
  meta_train <- metadata_for_keys(meta, train_keys)
  meta_apply <- metadata_for_keys(meta, apply_keys)
  scaled_train <- scale_by_training(expr, genes, train_keys, train_keys)
  x_train_genes <- scaled_train$x
  keep <- apply(x_train_genes, 2, stats::sd) > 0
  x_train_genes <- x_train_genes[, keep, drop = FALSE]
  genes <- colnames(x_train_genes)
  if (length(genes) < 2) stop("Too few nonzero-variance genes for PCA/glmnet.")
  pca <- prcomp(x_train_genes, center = FALSE, scale. = FALSE)
  max_k <- min(max_pc, ncol(pca$x), nrow(pca$x) - 2L)
  fold_source <- if (fold_col %in% names(meta_train)) meta_train[[fold_col]] else meta_train$dataset_folder
  foldid <- as.integer(factor(fold_source))
  records <- list()
  best <- NULL
  for (k in seq_len(max_k)) {
    tuned <- fit_alpha_grid(pca$x[, seq_len(k), drop = FALSE], meta_train$age_week, foldid = foldid)
    if (is.null(tuned)) next
    rec <- copy(tuned$best)
    rec[, pc_count := k]
    records[[length(records) + 1L]] <- rec
    if (is.null(best) || rec$cv_mae_week < best$cv_mae_week) best <- rec
  }
  if (is.null(best)) stop("PCA/glmnet tuning failed.")
  pc_train <- pca$x[, seq_len(best$pc_count), drop = FALSE]
  final <- glmnet(
    x = pc_train,
    y = meta_train$age_week,
    alpha = best$alpha,
    lambda = best$lambda,
    family = "gaussian",
    standardize = FALSE
  )
  x_apply_genes <- scale_by_training(expr, genes, train_keys, apply_keys)$x
  pc_apply <- predict(pca, newdata = x_apply_genes)[, seq_len(best$pc_count), drop = FALSE]
  raw_train <- as.numeric(predict(final, newx = pc_train))
  raw_apply <- as.numeric(predict(final, newx = pc_apply))
  cal <- calibrate_prediction(meta_train$age_week, raw_train, raw_apply, meta_train$sample_key)
  variance_dt <- data.table(
    pc = seq_along(pca$sdev),
    variance_explained = pca$sdev^2 / sum(pca$sdev^2)
  )
  selected_pc_names <- paste0("PC", seq_len(best$pc_count))
  pc_beta <- as.numeric(coef(final))[-1][seq_len(best$pc_count)]
  pc_coef_dt <- data.table(
    pc = selected_pc_names,
    raw_pc_coefficient = pc_beta,
    calibrated_pc_coefficient = cal$calibration$slope * pc_beta
  )
  calibrated_intercept <- cal$calibration$intercept + cal$calibration$slope * unname(as.numeric(coef(final))[1])
  model_component <- list(
    glmnet_fit = final,
    pca_fit = pca,
    pca_rotation = pca$rotation,
    pca_center = pca$center,
    selected_pc_names = selected_pc_names,
    selected_pc_count = as.integer(best$pc_count),
    alpha = as.numeric(best$alpha),
    lambda = as.numeric(best$lambda),
    glmnet_intercept = unname(as.numeric(coef(final))[1]),
    pc_coefficient_table = pc_coef_dt,
    pc_variance_explained = variance_dt,
    calibrated_intercept = calibrated_intercept,
    contribution_formula = "calibrated_PCj_contribution = calibrated_PCj_coefficient * PCj_score",
    raw_prediction_function_inputs = list(
      feature_space = "frozen_training_scaled_selected_genes_projected_to_frozen_pca_scores",
      standardize = FALSE
    )
  )
  fwrite(rbindlist(records), file.path(out_dir, paste0(fit_stage, "_pca_glmnet_cv_grid.tsv")), sep = "\t")
  fwrite(variance_dt, file.path(out_dir, paste0(fit_stage, "_pca_variance.tsv")), sep = "\t")
  saveRDS(
    list(
      model = final,
      pca = pca,
      genes = genes,
      center = scaled_train$center[genes],
      scale = scaled_train$scale[genes],
      pc_count = best$pc_count,
      training_scope = training_scope,
      calibration = cal$calibration,
      alpha = as.numeric(best$alpha),
      lambda = as.numeric(best$lambda),
      pc_coefficient_table = pc_coef_dt,
      pc_variance_explained = variance_dt,
      calibrated_intercept = calibrated_intercept
    ),
    file.path(out_dir, paste0(fit_stage, "_pca_glmnet_model.rds"))
  )
  list(
    predictions = prediction_table(meta_apply, raw_apply, cal$pred, "pca_glmnet", selector, fit_stage, cal$method),
    fit_artifact = list(
      model_family = "pca_glmnet",
      selector = selector,
      fit_stage = fit_stage,
      selected_genes = genes,
      gene_center = scaled_train$center[genes],
      gene_scale = scaled_train$scale[genes],
      model = model_component,
      calibration = cal$calibration,
      training_sample_keys = as.character(train_keys),
      apply_sample_keys = as.character(apply_keys),
      reference_metadata = copy(meta_train)
    )
  )
}

write_selected_gene_outputs <- function(fs, family, selector) {
  stub <- paste(training_scope, family, selector_file_stub(selector), "selected_genes", sep = "_")
  selected_path <- file.path(selected_gene_dir, paste0(stub, ".tsv"))
  out <- copy(fs$table)
  if (!"rank" %in% names(out)) out[, rank := seq_len(.N)]
  fwrite(out, selected_path, sep = "\t")
  if (selector == "version_b_age_dataset_penalized") {
    top_path <- file.path(selected_gene_dir, paste0(
      training_scope, "_", family, "_versionB_top30_genes.tsv"
    ))
    top_dt <- head(out[order(-feature_score)], 30)
    fwrite(top_dt, top_path, sep = "\t")
    plot_top30_versionb_genes(top_dt, family)
  }
  selected_path
}

write_bulk_group_alignment_outputs <- function() {
  dir.create(bulk_group_audit_dir, recursive = TRUE, showWarnings = FALSE)
  if (nrow(bulk_group_mapping) > 0) {
    fwrite(
      bulk_group_mapping,
      file.path(bulk_group_audit_dir, "bulk_group_manual_mapping_APPROVED.tsv"),
      sep = "\t"
    )
  }
  write_dt <- function(x, path) {
    if (length(x) > 0) {
      fwrite(rbindlist(x, fill = TRUE), path, sep = "\t")
    } else {
      fwrite(data.table(), path, sep = "\t")
    }
  }
  write_dt(
    bulk_group_added_samples,
    file.path(bulk_group_audit_dir, "model_ready_bulk_group_added_samples.tsv")
  )
  write_dt(
    bulk_group_expression_column_map,
    file.path(bulk_group_audit_dir, "model_ready_bulk_group_expression_column_map.tsv")
  )
  write_dt(
    bulk_group_dataset_summary,
    file.path(bulk_group_audit_dir, "model_ready_bulk_group_dataset_summary.tsv")
  )
  write_dt(
    bulk_group_gene_qc,
    file.path(bulk_group_audit_dir, "model_ready_bulk_group_gene_qc.tsv")
  )
  warnings_path <- file.path(bulk_group_audit_dir, "model_ready_bulk_group_warnings.txt")
  if (length(bulk_group_warnings) > 0) {
    writeLines(unique(bulk_group_warnings), warnings_path)
  } else {
    writeLines("No bulk group-level mapping warnings were emitted.", warnings_path)
  }
}

train_integrated_clock <- function(loaded) {
  expr <- loaded$expression
  meta_mod <- loaded$metadata
  if (!"cv_group" %in% names(meta_mod)) {
    meta_mod <- add_integrated_covariates(meta_mod, cv_group_col = "dataset_folder")
  }
  scope_stub <- safe_name(training_scope)
  message(
    "  Loaded ", training_scope, ": ", ncol(expr), " samples x ",
    nrow(expr), " common genes."
  )
  fwrite(loaded$alignment, file.path(results_dir, paste0(scope_stub, "_sample_expression_alignment.tsv")), sep = "\t")
  fwrite(loaded$alignment_issues, file.path(results_dir, paste0(scope_stub, "_alignment_issues.tsv")), sep = "\t")
  fwrite(data.table(gene_id = loaded$common_genes), file.path(results_dir, paste0(scope_stub, "_technical_gene_universe.tsv")), sep = "\t")
  write_bulk_group_alignment_outputs()
  if (!is.null(loaded$modality_common_genes)) {
    fwrite(loaded$modality_common_genes, file.path(results_dir, paste0(scope_stub, "_common_gene_summary.tsv")), sep = "\t")
  }

  train_keys <- meta_mod[split_role == "train", sample_key]
  validation_keys <- meta_mod[split_role == "validation", sample_key]
  test_keys <- meta_mod[split_role == "test", sample_key]
  train_validation_keys <- meta_mod[split_role %in% c("train", "validation"), sample_key]

  all_predictions <- list()
  all_validation_metrics <- list()
  selected_gene_counts <- list()
  fit_artifacts <- list()

  for (selector in c("hvg_baseline", "version_b_age_dataset_penalized")) {
    for (family in c("pca_glmnet", "elasticnet")) {
      message("  Training ", training_scope, " / ", family, " / ", selector, ".")
      family_dir <- if (family == "pca_glmnet") pca_dir else elasticnet_dir
      out_dir <- file.path(family_dir, scope_stub, selector)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      fs <- prepare_version_genes(selector, expr[, train_keys, drop = FALSE], metadata_for_keys(meta_mod, train_keys), out_dir)
      fwrite(data.table(gene_id = fs$genes, rank = seq_along(fs$genes)), file.path(out_dir, "selected_genes.tsv"), sep = "\t")
      selected_gene_path <- write_selected_gene_outputs(fs, family, selector)
      selected_gene_counts[[length(selected_gene_counts) + 1L]] <- data.table(
        training_scope = training_scope,
        model_family = family,
        selector = selector,
        selected_gene_count = length(fs$genes),
        selected_gene_path = selected_gene_path
      )

      fit_fun <- if (family == "pca_glmnet") fit_pca_glmnet_clock else fit_elasticnet_clock
      pred_val <- fit_fun(expr, meta_mod, train_keys, validation_keys, fs$genes, selector, out_dir, "train_to_validation")
      pred_train <- fit_fun(expr, meta_mod, train_keys, train_keys, fs$genes, selector, out_dir, "train_fit")
      all_predictions[[length(all_predictions) + 1L]] <- pred_val$predictions
      all_predictions[[length(all_predictions) + 1L]] <- pred_train$predictions
      fit_artifacts[[paste(family, selector, "train_to_validation", sep = "\r")]] <- pred_val$fit_artifact
      fit_artifacts[[paste(family, selector, "train_fit", sep = "\r")]] <- pred_train$fit_artifact
      val_metric <- metric_one(pred_val$predictions, c("training_scope", "model_family", "selector", "fit_stage", "split_role"))
      all_validation_metrics[[length(all_validation_metrics) + 1L]] <- val_metric
    }
  }

  validation_metrics <- rbindlist(all_validation_metrics, fill = TRUE)
  locked <- validation_metrics[fit_stage == "train_to_validation"][order(mae_week), .SD[1], by = .(training_scope, model_family)]
  fwrite(locked, file.path(results_dir, paste0(scope_stub, "_locked_model_selection.tsv")), sep = "\t")

  for (i in seq_len(nrow(locked))) {
    family <- locked$model_family[i]
    selector <- locked$selector[i]
    message("  Locked refit ", training_scope, " / ", family, " / ", selector, ".")
    family_dir <- if (family == "pca_glmnet") pca_dir else elasticnet_dir
    out_dir <- file.path(family_dir, scope_stub, selector)
    fs <- prepare_version_genes(selector, expr[, train_validation_keys, drop = FALSE], metadata_for_keys(meta_mod, train_validation_keys), out_dir)
    fwrite(data.table(gene_id = fs$genes, rank = seq_along(fs$genes)), file.path(out_dir, "selected_genes.tsv"), sep = "\t")
    selected_gene_path <- write_selected_gene_outputs(fs, family, selector)
    selected_gene_counts[[length(selected_gene_counts) + 1L]] <- data.table(
      training_scope = training_scope,
      model_family = family,
      selector = paste0(selector, "_locked_train_validation"),
      selected_gene_count = length(fs$genes),
      selected_gene_path = selected_gene_path
    )
    fit_fun <- if (family == "pca_glmnet") fit_pca_glmnet_clock else fit_elasticnet_clock
    pred_test <- fit_fun(expr, meta_mod, train_validation_keys, test_keys, fs$genes, selector, out_dir, "locked_train_validation_to_test")
    all_predictions[[length(all_predictions) + 1L]] <- pred_test$predictions
    fit_artifacts[[paste(family, selector, "locked_train_validation_to_test", sep = "\r")]] <- pred_test$fit_artifact
  }

  pred <- rbindlist(all_predictions, fill = TRUE)
  performance <- summarise_predictions(pred)
  selected_gene_counts <- rbindlist(selected_gene_counts, fill = TRUE)
  fwrite(pred, file.path(results_dir, paste0(scope_stub, "_all_clock_predictions.tsv")), sep = "\t")
  fwrite(performance, file.path(results_dir, paste0(scope_stub, "_all_clock_performance.tsv")), sep = "\t")
  fwrite(selected_gene_counts, file.path(selected_gene_dir, "selected_gene_count_summary.tsv"), sep = "\t")
  message("  Wrote ", training_scope, " predictions and performance tables.")

  for (family in c("pca_glmnet", "elasticnet")) {
    family_pred <- pred[model_family == family]
    family_perf <- performance[model_family == family]
    fwrite(family_pred, file.path(if (family == "pca_glmnet") pca_dir else elasticnet_dir, scope_stub, "predictions.tsv"), sep = "\t")
    fwrite(family_perf, file.path(if (family == "pca_glmnet") pca_dir else elasticnet_dir, scope_stub, "performance.tsv"), sep = "\t")
  }
  list(
    predictions = pred,
    performance = performance,
    locked = locked,
    loaded = loaded,
    selected_gene_counts = selected_gene_counts,
    fit_artifacts = fit_artifacts
  )
}

retained_clock_table <- data.table(model_family = character(), selector = character())

get_retained_clock_table <- function(result = NULL) {
  if (!is.null(result) && "locked" %in% names(result) && nrow(result$locked) > 0) {
    retained <- unique(ensure_dt(result$locked)[, .(model_family, selector)])
    setorder(retained, model_family, selector)
    return(retained)
  }
  retained_clock_table
}

artifact_filename <- function(model_family, selector, artifact_scope) {
  file.path(
    if (artifact_scope == "benchmark") benchmark_artifact_dir else deployment_artifact_dir,
    paste0(paste(training_scope, model_family, selector, artifact_scope, sep = "__"), ".rds")
  )
}

stable_artifact <- function(fit, artifact_scope, control_like_sample_keys, validation_sample_keys = character(),
                            test_sample_keys = character(), benchmark_fits = NULL) {
  artifact <- list(
    schema_version = "integrated_ovary_clock_artifact_v1",
    training_scope = training_scope,
    artifact_scope = artifact_scope,
    model_family = fit$model_family,
    selector = fit$selector,
    expression_scale = "logCPM",
    selected_genes = fit$selected_genes,
    gene_center = fit$gene_center,
    gene_scale = fit$gene_scale,
    missing_gene_policy = "impute_training_center_then_scale_to_zero",
    model = fit$model,
    calibration = fit$calibration,
    reference_metadata = fit$reference_metadata,
    training_sample_keys = fit$training_sample_keys,
    validation_sample_keys = validation_sample_keys,
    test_sample_keys = test_sample_keys,
    control_like_sample_keys = control_like_sample_keys,
    created_at = Sys.time()
  )
  if (!is.null(benchmark_fits)) artifact$benchmark_fits <- benchmark_fits
  artifact
}

write_benchmark_artifacts <- function(result) {
  retained <- get_retained_clock_table(result)
  manifest <- list()
  calibration_rows <- list()
  elastic_rows <- list()
  pca_rows <- list()
  control_like_keys <- result$loaded$metadata$sample_key
  validation_keys <- result$loaded$metadata[split_role == "validation", sample_key]
  test_keys <- result$loaded$metadata[split_role == "test", sample_key]

  for (i in seq_len(nrow(retained))) {
    family <- retained$model_family[i]
    selector <- retained$selector[i]
    stage_keys <- paste(family, selector, c("train_to_validation", "train_fit", "locked_train_validation_to_test"), sep = "\r")
    benchmark_fits <- result$fit_artifacts[stage_keys]
    names(benchmark_fits) <- c("train_to_validation", "train_fit", "locked_train_validation_to_test")
    top_fit <- benchmark_fits[["locked_train_validation_to_test"]]
    artifact <- stable_artifact(top_fit, "benchmark", control_like_keys, validation_keys, test_keys, benchmark_fits)
    path <- artifact_filename(family, selector, "benchmark")
    saveRDS(artifact, path)
    manifest[[length(manifest) + 1L]] <- data.table(
      artifact_scope = "benchmark",
      model_family = family,
      selector = selector,
      artifact_path = path,
      n_selected_genes = length(artifact$selected_genes),
      calibration_n = artifact$calibration$calibration_n
    )
    for (stage in names(benchmark_fits)) {
      cal <- benchmark_fits[[stage]]$calibration
      calibration_rows[[length(calibration_rows) + 1L]] <- data.table(
        artifact_scope = "benchmark",
        fit_stage = stage,
        model_family = family,
        selector = selector,
        calibration_type = cal$calibration_type,
        intercept = cal$intercept,
        slope = cal$slope,
        calibration_n = cal$calibration_n,
        calibration_age_min = cal$calibration_age_range[1],
        calibration_age_max = cal$calibration_age_range[2],
        calibration_raw_prediction_min = cal$calibration_raw_prediction_range[1],
        calibration_raw_prediction_max = cal$calibration_raw_prediction_range[2],
        clip_predictions = cal$clip_predictions
      )
    }
    if (family == "elasticnet") {
      coefficient_table <- ensure_dt(copy(artifact$model$nonzero_coefficient_table))
      coefficient_table[, `:=`(
        artifact_scope = "benchmark",
        model_family = family,
        selector = selector
      )]
      elastic_rows[[length(elastic_rows) + 1L]] <- coefficient_table
    } else {
      coefficient_table <- ensure_dt(copy(artifact$model$pc_coefficient_table))
      coefficient_table[, `:=`(
        artifact_scope = "benchmark",
        model_family = family,
        selector = selector,
        calibrated_intercept = artifact$model$calibrated_intercept
      )]
      pca_rows[[length(pca_rows) + 1L]] <- coefficient_table
    }
  }
  list(
    manifest = ensure_dt(rbindlist(manifest, fill = TRUE)),
    calibration = ensure_dt(rbindlist(calibration_rows, fill = TRUE)),
    elasticnet_nonzero = ensure_dt(rbindlist(elastic_rows, fill = TRUE)),
    pca_coefficients = ensure_dt(rbindlist(pca_rows, fill = TRUE))
  )
}

scale_with_artifact <- function(expr, keys, genes, center, scale) {
  missing <- setdiff(genes, rownames(expr))
  if (length(missing) > 0) stop("Expression matrix is missing artifact genes: ", paste(head(missing, 20), collapse = ", "), call. = FALSE)
  x <- t(expr[genes, keys, drop = FALSE])
  scale[!is.finite(scale) | scale == 0] <- 1
  x <- sweep(x, 2, center[genes], "-")
  x <- sweep(x, 2, scale[genes], "/")
  x[!is.finite(x)] <- 0
  x
}

predict_from_fit_artifact <- function(fit, expr, keys) {
  x <- scale_with_artifact(expr, keys, fit$selected_genes, fit$gene_center, fit$gene_scale)
  if (fit$model_family == "elasticnet") {
    raw <- as.numeric(predict(fit$model$glmnet_fit, newx = x))
  } else {
    pc <- predict(fit$model$pca_fit, newdata = x)[, seq_len(fit$model$selected_pc_count), drop = FALSE]
    raw <- as.numeric(predict(fit$model$glmnet_fit, newx = pc))
  }
  pred <- fit$calibration$intercept + fit$calibration$slope * raw
  data.table(sample_key = keys, raw_predicted_age_week_recomputed = raw, predicted_age_week_recomputed = pred)
}

write_reproducibility_check <- function(result, all_pred) {
  retained <- get_retained_clock_table(result)
  checks <- list()
  expr <- result$loaded$expression
  for (i in seq_len(nrow(retained))) {
    family <- retained$model_family[i]
    selector_value <- retained$selector[i]
    artifact <- readRDS(artifact_filename(family, selector_value, "benchmark"))
    for (stage in names(artifact$benchmark_fits)) {
      pred_rows <- all_pred[model_family == family & selector == selector_value & fit_stage == stage]
      if (nrow(pred_rows) == 0) next
      fit <- artifact$benchmark_fits[[stage]]
      recomputed <- predict_from_fit_artifact(fit, expr, pred_rows$sample_key)
      cmp <- merge(
        pred_rows[, .(sample_key, predicted_age_week_saved = predicted_age_week)],
        recomputed,
        by = "sample_key",
        all = FALSE
      )
      max_abs_diff <- max(abs(cmp$predicted_age_week_saved - cmp$predicted_age_week_recomputed), na.rm = TRUE)
      checks[[length(checks) + 1L]] <- data.table(
        artifact_scope = "benchmark",
        model_family = family,
        selector = selector_value,
        fit_stage = stage,
        n_checked = nrow(cmp),
        max_abs_predicted_age_diff = max_abs_diff,
        reproducibility_status = ifelse(is.finite(max_abs_diff) && max_abs_diff <= 1e-6, "pass", "fail")
      )
    }
  }
  out <- rbindlist(checks, fill = TRUE)
  fwrite(out, file.path(artifact_root_dir, "artifact_reproducibility_check.tsv"), sep = "\t")
  if (any(out$reproducibility_status != "pass")) {
    stop("Benchmark artifact reproducibility check failed.", call. = FALSE)
  }
  out
}

write_deployment_artifacts <- function(result) {
  retained <- get_retained_clock_table(result)
  expr <- result$loaded$expression
  meta_mod <- result$loaded$metadata
  control_like_keys <- meta_mod$sample_key
  manifest <- list()
  calibration_rows <- list()
  elastic_rows <- list()
  pca_rows <- list()
  audits <- list()

  for (i in seq_len(nrow(retained))) {
    family <- retained$model_family[i]
    selector <- retained$selector[i]
    message("  Deployment refit ", training_scope, " / ", family, " / ", selector, ".")
    out_dir <- file.path(deployment_artifact_dir, paste(family, selector, sep = "__"))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fs <- prepare_version_genes(selector, expr[, control_like_keys, drop = FALSE], metadata_for_keys(meta_mod, control_like_keys), out_dir)
    fit_fun <- if (family == "pca_glmnet") fit_pca_glmnet_clock else fit_elasticnet_clock
    dep <- fit_fun(expr, meta_mod, control_like_keys, control_like_keys, fs$genes, selector, out_dir, "deployment_all_control_like")
    artifact <- stable_artifact(dep$fit_artifact, "deployment_all_control_like", control_like_keys)
    path <- artifact_filename(family, selector, "deployment_all_control_like")
    saveRDS(artifact, path)
    manifest[[length(manifest) + 1L]] <- data.table(
      artifact_scope = "deployment_all_control_like",
      model_family = family,
      selector = selector,
      artifact_path = path,
      n_selected_genes = length(artifact$selected_genes),
      calibration_n = artifact$calibration$calibration_n
    )
    cal <- artifact$calibration
    calibration_rows[[length(calibration_rows) + 1L]] <- data.table(
      artifact_scope = "deployment_all_control_like",
      fit_stage = "deployment_all_control_like",
      model_family = family,
      selector = selector,
      calibration_type = cal$calibration_type,
      intercept = cal$intercept,
      slope = cal$slope,
      calibration_n = cal$calibration_n,
      calibration_age_min = cal$calibration_age_range[1],
      calibration_age_max = cal$calibration_age_range[2],
      calibration_raw_prediction_min = cal$calibration_raw_prediction_range[1],
      calibration_raw_prediction_max = cal$calibration_raw_prediction_range[2],
      clip_predictions = cal$clip_predictions
    )
    if (family == "elasticnet") {
      coefficient_table <- ensure_dt(copy(artifact$model$nonzero_coefficient_table))
      coefficient_table[, `:=`(
        artifact_scope = "deployment_all_control_like",
        model_family = family,
        selector = selector
      )]
      elastic_rows[[length(elastic_rows) + 1L]] <- coefficient_table
    } else {
      coefficient_table <- ensure_dt(copy(artifact$model$pc_coefficient_table))
      coefficient_table[, `:=`(
        artifact_scope = "deployment_all_control_like",
        model_family = family,
        selector = selector,
        calibrated_intercept = artifact$model$calibrated_intercept
      )]
      pca_rows[[length(pca_rows) + 1L]] <- coefficient_table
    }
    audit_meta <- meta_mod[sample_key %chin% artifact$training_sample_keys]
    cal_meta <- meta_mod[sample_key %chin% artifact$calibration$calibration_sample_keys]
    audits[[length(audits) + 1L]] <- data.table(
      model_family = family,
      selector = selector,
      artifact_path = path,
      training_sample_count = length(artifact$training_sample_keys),
      calibration_sample_count = length(artifact$calibration$calibration_sample_keys),
      all_training_control_like = all(audit_meta$control_like_clock_candidate == TRUE),
      no_training_exclusion_reason = all(is.na(audit_meta$exclusion_reason) | audit_meta$exclusion_reason == ""),
      all_calibration_control_like = all(cal_meta$control_like_clock_candidate == TRUE),
      selected_genes_present = length(artifact$selected_genes) > 0,
      centers_present = length(artifact$gene_center) == length(artifact$selected_genes),
      scales_present = length(artifact$gene_scale) == length(artifact$selected_genes),
      model_present = !is.null(artifact$model),
      calibration_present = !is.null(artifact$calibration),
      deployment_ready = all(audit_meta$control_like_clock_candidate == TRUE) &&
        all(is.na(audit_meta$exclusion_reason) | audit_meta$exclusion_reason == "") &&
        all(cal_meta$control_like_clock_candidate == TRUE) &&
        length(artifact$selected_genes) > 0 &&
        length(artifact$gene_center) == length(artifact$selected_genes) &&
        length(artifact$gene_scale) == length(artifact$selected_genes) &&
        !is.null(artifact$model) &&
        !is.null(artifact$calibration)
    )
  }
  audit <- rbindlist(audits, fill = TRUE)
  fwrite(audit, file.path(artifact_root_dir, "deployment_artifact_audit.tsv"), sep = "\t")
  if (any(audit$deployment_ready != TRUE)) {
    stop("Deployment artifact audit failed.", call. = FALSE)
  }
  list(
    manifest = ensure_dt(rbindlist(manifest, fill = TRUE)),
    calibration = ensure_dt(rbindlist(calibration_rows, fill = TRUE)),
    elasticnet_nonzero = ensure_dt(rbindlist(elastic_rows, fill = TRUE)),
    pca_coefficients = ensure_dt(rbindlist(pca_rows, fill = TRUE)),
    audit = ensure_dt(audit)
  )
}

field_present <- function(x, field) {
  parts <- strsplit(field, ".", fixed = TRUE)[[1]]
  value <- x
  for (part in parts) {
    if (is.null(value) || is.null(value[[part]])) return(FALSE)
    value <- value[[part]]
  }
  length(value) > 0 && !all(is.na(value))
}

artifact_model_present <- function(artifact) {
  if (is.null(artifact$model)) return(FALSE)
  if (!is.null(artifact$model$glmnet_fit)) return(TRUE)
  if (!is.null(artifact$model$nonzero_coefficient_table) && nrow(ensure_dt(artifact$model$nonzero_coefficient_table)) > 0) return(TRUE)
  if (!is.null(artifact$model$pc_coefficient_table) && nrow(ensure_dt(artifact$model$pc_coefficient_table)) > 0) return(TRUE)
  FALSE
}

write_artifact_schema_check <- function(manifest, metadata) {
  manifest <- ensure_dt(manifest)
  metadata <- ensure_dt(metadata)
  required_fields <- c(
    "schema_version", "training_scope", "artifact_scope", "model_family", "selector",
    "selected_genes", "gene_center", "gene_scale", "calibration.intercept",
    "calibration.slope"
  )
  checks <- list()
  for (i in seq_len(nrow(manifest))) {
    path <- manifest$artifact_path[i]
    artifact <- if (file.exists(path)) readRDS(path) else NULL
    scope <- manifest$artifact_scope[i]
    family <- manifest$model_family[i]
    selector <- manifest$selector[i]
    add_check <- function(check, ok, message) {
      checks[[length(checks) + 1L]] <<- data.table(
        artifact_path = path,
        artifact_scope = scope,
        model_family = family,
        selector = selector,
        check = check,
        status = ifelse(ok, "PASS", "FAIL"),
        message = message
      )
    }
    add_check("artifact_file_exists", file.exists(path), ifelse(file.exists(path), "Artifact file exists.", "Artifact file is missing."))
    if (is.null(artifact)) next
    for (field in required_fields) {
      add_check(
        paste0("required_field:", field),
        field_present(artifact, field),
        ifelse(field_present(artifact, field), "Required field is present.", "Required field is missing.")
      )
    }
    add_check("model_or_coefficients_present", artifact_model_present(artifact), "Model object or coefficient table must be present.")
    add_check("training_scope_matches", identical(artifact$training_scope, training_scope), "Artifact training_scope must match the current run.")
    add_check("artifact_scope_matches", identical(artifact$artifact_scope, scope), "Artifact artifact_scope must match manifest.")
    if (grepl("deployment", scope)) {
      training_keys <- if (is.null(artifact$training_sample_keys)) character() else artifact$training_sample_keys
      training_meta <- metadata[sample_key %chin% training_keys]
      add_check("deployment_scope_contains_deployment", grepl("deployment", artifact$artifact_scope), "Deployment artifact_scope must contain deployment.")
      add_check("deployment_training_keys_exist", length(training_keys) > 0 && nrow(training_meta) == length(unique(training_keys)), "Deployment training sample keys must exist in model-ready metadata.")
      add_check("deployment_training_control_like", nrow(training_meta) > 0 && all(training_meta$control_like_clock_candidate == TRUE), "Deployment training samples must be control-like.")
      add_check(
        "deployment_no_excluded_samples",
        nrow(training_meta) > 0 && all(is.na(training_meta$exclusion_reason) | training_meta$exclusion_reason == ""),
        "Excluded/intervention samples must not be used in deployment artifacts."
      )
    }
  }
  out <- if (length(checks) > 0) {
    ensure_dt(rbindlist(checks, fill = TRUE))
  } else {
    data.table(
      artifact_path = character(), artifact_scope = character(), model_family = character(),
      selector = character(), check = character(), status = character(), message = character()
    )
  }
  fwrite(out, file.path(artifact_root_dir, "artifact_schema_check.tsv"), sep = "\t")
  if (nrow(out[status != "PASS"]) > 0) {
    stop("Artifact schema check failed. See model_artifacts/artifact_schema_check.tsv.", call. = FALSE)
  }
  out
}

write_model_artifact_outputs <- function(result, all_pred) {
  benchmark <- write_benchmark_artifacts(result)
  reproducibility <- write_reproducibility_check(result, all_pred)
  deployment <- write_deployment_artifacts(result)

  manifest <- ensure_dt(rbindlist(list(benchmark$manifest, deployment$manifest), fill = TRUE))
  calibration <- ensure_dt(rbindlist(list(benchmark$calibration, deployment$calibration), fill = TRUE))
  elasticnet_nonzero <- ensure_dt(rbindlist(list(benchmark$elasticnet_nonzero, deployment$elasticnet_nonzero), fill = TRUE))
  pca_coefficients <- ensure_dt(rbindlist(list(benchmark$pca_coefficients, deployment$pca_coefficients), fill = TRUE))

  fwrite(manifest, file.path(artifact_root_dir, "model_artifact_manifest.tsv"), sep = "\t")
  fwrite(calibration, file.path(artifact_root_dir, "calibration_parameters.tsv"), sep = "\t")
  fwrite(pca_coefficients, file.path(artifact_root_dir, "pca_calibrated_pc_coefficients.tsv"), sep = "\t")
  fwrite(elasticnet_nonzero, file.path(artifact_root_dir, "elasticnet_nonzero_coefficients.tsv"), sep = "\t")
  schema_check <- write_artifact_schema_check(manifest, result$loaded$metadata)

  list(
    manifest = manifest,
    calibration = calibration,
    reproducibility = reproducibility,
    deployment_audit = deployment$audit,
    schema_check = schema_check
  )
}

format_table_md <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return("_None._")
  dt <- as.data.table(dt)
  header <- paste("|", paste(names(dt), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(dt)), collapse = " | "), "|")
  rows <- apply(dt, 1, function(x) paste("|", paste(as.character(x), collapse = " | "), "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

write_run_manifest <- function(loaded, result, all_pred, all_perf, all_locked, artifact_outputs, status = "completed") {
  manifest_path <- file.path(results_dir, "run_manifest.tsv")
  manifest_md_path <- file.path(results_dir, "run_manifest.md")
  mapping_used <- length(bulk_group_mapping_used_datasets) > 0
  exact_target_used <- intersect(samplelevel_exact_priority_datasets, sort(unique(exact_samplelevel_matching_used_for)))
  added_samples <- if (length(bulk_group_added_samples) > 0) {
    rbindlist(bulk_group_added_samples, fill = TRUE)
  } else {
    data.table()
  }
  before_counts <- meta[control_like_clock_candidate == TRUE, .N, by = .(modality, split_role)]
  after_counts <- loaded$metadata[, .N, by = .(modality, split_role)]
  split_dataset <- loaded$metadata[, .N, by = .(split_role, dataset_folder)][order(split_role, dataset_folder)]
  split_subcohort <- loaded$metadata[, .N, by = .(split_role, dataset_folder, dataset_subcohort)][order(split_role, dataset_folder, dataset_subcohort)]
  split_modality <- loaded$metadata[, .N, by = .(split_role, modality)][order(split_role, modality)]
  split_stratum <- loaded$metadata[, .N, by = .(split_role, split_age_stratum)][order(split_role, split_age_stratum)]
  locked <- all_locked[, .(model_family, selector, split_role, n_samples, mae_week, rmse_week)]
  retained_eval <- all_perf[
    fit_stage %in% c("train_to_validation", "locked_train_validation_to_test") &
      split_role %in% c("validation", "test") &
      is.na(modality) &
      is.na(dataset_folder) &
      is.na(dataset_subcohort) &
      is.na(biological_age_bin) &
      is.na(split_age_stratum),
    .(model_family, selector, split_role, fit_stage, n_samples, mae_week, rmse_week)
  ][order(model_family, split_role)]
  deploy_paths <- artifact_outputs$manifest[
    artifact_scope == "deployment_all_control_like",
    .(model_family, selector, artifact_path)
  ]

  rows <- rbindlist(list(
    manifest_row("run", "run_label", run_label),
    manifest_row("run", "timestamp", as.character(Sys.time())),
    manifest_row("run", "script_path", normalizePath(script_file %||% file.path(script_dir, "03_train_bulk_pseudobulk_clocks.R"), mustWork = FALSE)),
    manifest_row("run", "results_dir", results_dir),
    manifest_row("run", "status", status),
    manifest_row("run", "failure_stage", NA_character_),
    manifest_row("run", "failure_message", NA_character_),
    manifest_row("mapping", "approved_mapping_path", bulk_group_mapping_path),
    manifest_row("mapping", "bulk_group_mapping_mode", bulk_group_mapping_mode),
    manifest_row("mapping", "group_level_bulk_mapping_used", as.character(mapping_used)),
    manifest_row("mapping", "manual_group_mapping_used", as.character(mapping_used)),
    manifest_row("mapping", "manual_group_mapping_used_for", paste(sort(unique(bulk_group_mapping_used_datasets)), collapse = ",")),
    manifest_row("mapping", "exact_samplelevel_matching_used_for", paste(exact_target_used, collapse = ",")),
    manifest_row("mapping", "added_datasets", paste(sort(unique(added_samples$dataset_folder)), collapse = ";")),
    manifest_row("samples", "control_like_samples_before_mapping", as.character(meta[control_like_clock_candidate == TRUE, .N])),
    manifest_row("samples", "model_ready_samples_after_mapping", as.character(ncol(loaded$expression))),
    manifest_row("genes", "common_gene_count", as.character(nrow(loaded$expression))),
    manifest_row("artifacts", "deployment_artifact_paths", paste(deploy_paths$artifact_path, collapse = ";")),
    manifest_row("artifacts", "artifact_schema_check_path", file.path(artifact_root_dir, "artifact_schema_check.tsv")),
    manifest_row("artifacts", "artifact_reproducibility_check_path", file.path(artifact_root_dir, "artifact_reproducibility_check.tsv"))
  ), fill = TRUE)
  fwrite(rows, manifest_path, sep = "\t")

  md <- c(
    paste0("# Run Manifest: ", run_label),
    "",
    paste0("- Timestamp: ", Sys.time()),
    paste0("- Results dir: `", results_dir, "`"),
    paste0("- Status: ", status),
    paste0("- Approved mapping path: `", bulk_group_mapping_path, "`"),
    paste0("- Bulk group mapping mode: ", bulk_group_mapping_mode),
    paste0("- Manual group mapping used: ", mapping_used),
    paste0("- Exact sample-level matching used for: ", paste(exact_target_used, collapse = ", ")),
    paste0("- Added datasets: ", paste(sort(unique(added_samples$dataset_folder)), collapse = ", ")),
    paste0("- Common gene count: ", nrow(loaded$expression)),
    "",
    "## Sample Counts Before Mapping",
    format_table_md(before_counts[order(modality, split_role)]),
    "",
    "## Sample Counts After Mapping",
    format_table_md(after_counts[order(modality, split_role)]),
    "",
    "## Split Counts By Dataset",
    format_table_md(split_dataset),
    "",
    "## Split Counts By Dataset Subcohort",
    format_table_md(split_subcohort),
    "",
    "## Split Counts By Modality",
    format_table_md(split_modality),
    "",
    "## Split Counts By Split Age Stratum",
    format_table_md(split_stratum),
    "",
    "## Locked Selector Per Model Family",
    format_table_md(locked),
    "",
    "## Validation/Test MAE Per Retained Clock",
    format_table_md(retained_eval),
    "",
    "## Deployment Artifact Paths",
    format_table_md(deploy_paths)
  )
  writeLines(md, manifest_md_path)
  invisible(list(tsv = manifest_path, md = manifest_md_path))
}

meta <- fread(metadata_path)
meta[, age_week := as.numeric(age_week)]
meta[, n_cells_used := suppressWarnings(as.numeric(n_cells_used))]
meta <- add_dataset_subcohort(meta)
metadata_check_issues <- run_severe_metadata_checks(meta)
fwrite(
  metadata_check_issues,
  file.path(results_dir, "samplelevel_metadata_check_issues.tsv"),
  sep = "\t"
)
fwrite(
  meta[, .N, by = .(dataset_folder, dataset_subcohort, control_like_clock_candidate, split_role)][
    order(dataset_folder, dataset_subcohort, -control_like_clock_candidate, split_role)
  ],
  file.path(results_dir, "dataset_subcohort_summary.tsv"),
  sep = "\t"
)
if (nrow(metadata_check_issues[severity == "warning"]) > 0) {
  warning(
    "Non-severe metadata warnings detected: ",
    paste(metadata_check_issues[severity == "warning", issue], collapse = " | "),
    call. = FALSE
  )
}
setkey(meta, sample_key)

message("Training scope: ", training_scope)
current_stage <- "loading_integrated_expression"
integrated_loaded <- load_integrated_expression(meta)
current_stage <- "model_fitting"
result <- train_integrated_clock(integrated_loaded)

all_pred <- result$predictions
all_perf <- result$performance
all_locked <- result$locked
fwrite(all_pred, file.path(results_dir, "all_clock_predictions.tsv"), sep = "\t")
fwrite(all_perf, file.path(results_dir, "all_clock_performance.tsv"), sep = "\t")
fwrite(all_locked, file.path(results_dir, "locked_model_selection.tsv"), sep = "\t")
set_run_stage("model_predictions_written")

if (uniqueN(all_pred$training_scope) != 1L || unique(all_pred$training_scope) != training_scope) {
  stop("Quality check failed: more than one training_scope was produced.")
}
if (uniqueN(all_perf$training_scope) != 1L || unique(all_perf$training_scope) != training_scope) {
  stop("Quality check failed: all_clock_performance.tsv contains an unexpected training_scope.")
}
if (uniqueN(all_locked$training_scope) != 1L || unique(all_locked$training_scope) != training_scope) {
  stop("Quality check failed: locked_model_selection.tsv contains an unexpected training_scope.")
}
if (any(grepl("bulk_rnaseq|scrna_pseudobulk", all_pred$training_scope))) {
  stop("Quality check failed: model scope contains a modality name.")
}
selected_files <- list.files(selected_gene_dir, pattern = "selected_genes\\.tsv$", recursive = TRUE)
if (any(grepl("bulk_rnaseq|scrna_pseudobulk", selected_files))) {
  stop("Quality check failed: selected-gene files are separated by modality.")
}

current_stage <- "artifact_writing"
artifact_outputs <- write_model_artifact_outputs(result, all_pred)
set_run_stage("artifacts_written")
run_manifest <- write_run_manifest(integrated_loaded, result, all_pred, all_perf, all_locked, artifact_outputs, status = "completed")
current_stage <- "completed"

overall_eval <- all_perf[
  fit_stage %in% c("train_to_validation", "locked_train_validation_to_test") &
    split_role %in% c("validation", "test") &
    is.na(modality) &
    is.na(dataset_folder) &
    is.na(dataset_subcohort) &
    is.na(biological_age_bin) &
    is.na(split_age_stratum)
]
modality_eval <- all_perf[
  fit_stage %in% c("train_to_validation", "locked_train_validation_to_test") &
    split_role %in% c("validation", "test") &
    !is.na(modality) &
    is.na(dataset_folder) &
    is.na(dataset_subcohort) &
    is.na(biological_age_bin) &
    is.na(split_age_stratum)
]

cat("\nIntegrated clock run summary\n")
cat("Training scope:", training_scope, "\n")
cat("Integrated samples:", ncol(integrated_loaded$expression), "\n")
cat("Bulk samples:", integrated_loaded$metadata[modality == "bulk_rnaseq", .N], "\n")
cat("Pseudobulk samples:", integrated_loaded$metadata[modality == "scrna_pseudobulk", .N], "\n")
cat("Common genes:", nrow(integrated_loaded$expression), "\n")
cat("\nTrain/validation/test sample counts by modality:\n")
print(integrated_loaded$metadata[, .N, by = .(split_role, modality)][order(split_role, modality)])
cat("\nTrain/validation/test sample counts by dataset_folder:\n")
print(integrated_loaded$metadata[, .N, by = .(split_role, dataset_folder)][order(split_role, dataset_folder)])
cat("\nSelected gene count per model_family + selector:\n")
print(result$selected_gene_counts[!grepl("_locked_train_validation$", selector), .(
  selected_gene_count = max(selected_gene_count)
), by = .(model_family, selector)][order(model_family, selector)])
cat("\nLocked integrated model choices:\n")
print(all_locked[, .(training_scope, model_family, selector, split_role, n_samples, mae_week, rmse_week)])
cat("\nValidation and test MAE overall:\n")
print(overall_eval[, .(model_family, selector, split_role, fit_stage, n_samples, mae_week)][order(model_family, split_role)])
cat("\nValidation and test MAE by modality:\n")
print(modality_eval[, .(model_family, selector, split_role, modality, fit_stage, n_samples, mae_week)][order(model_family, split_role, modality)])

cat("\nProjection-ready artifact paths:\n")
print(artifact_outputs$manifest[, .(artifact_scope, model_family, selector, artifact_path)])
cat("\nCalibration intercept/slope for retained clocks:\n")
print(artifact_outputs$calibration[
  artifact_scope %in% c("benchmark", "deployment_all_control_like") &
    fit_stage %in% c("locked_train_validation_to_test", "deployment_all_control_like"),
  .(artifact_scope, model_family, selector, fit_stage, intercept, slope, calibration_n)
][order(artifact_scope, model_family)])
cat("\nBenchmark reproducibility status:\n")
print(artifact_outputs$reproducibility[, .(model_family, selector, fit_stage, n_checked, max_abs_predicted_age_diff, reproducibility_status)])
cat("\nDeployment artifact readiness status:\n")
print(artifact_outputs$deployment_audit[, .(model_family, selector, training_sample_count, calibration_sample_count, deployment_ready)])

cat("Wrote clock results under:", results_dir, "\n")
cat("PCA/glmnet results:", pca_dir, "\n")
cat("ElasticNet results:", elasticnet_dir, "\n")
cat("Model artifacts:", artifact_root_dir, "\n")
cat("Run manifest TSV:", run_manifest$tsv, "\n")
cat("Run manifest MD:", run_manifest$md, "\n")
