#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
if (!dir.exists(script_dir) || !file.exists(file.path(script_dir, "03_train_bulk_pseudobulk_clocks.R"))) {
  script_dir <- file.path(getwd(), "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_input_dir <- normalizePath(file.path(project_dir, "..", "data_input"), mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
out_dir <- file.path(results_dir, "alignment_audit", "bulk_group_level_alignment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

default_mapping_path <- file.path(out_dir, "bulk_group_manual_mapping_APPROVED.tsv")
mapping_path <- Sys.getenv("BULK_GROUP_MAPPING_FILE", unset = default_mapping_path)
if (!file.exists(mapping_path)) {
  message("APPROVED_MAPPING_VALIDATION_FAIL")
  stop("Approved mapping file is missing: ", mapping_path, call. = FALSE)
}
if (!file.exists(metadata_path)) stop("Missing metadata file: ", metadata_path, call. = FALSE)

required_cols <- c(
  "dataset_folder", "metadata_row_id", "metadata_sample_key", "metadata_sample_id",
  "metadata_group_label", "metadata_age_week", "metadata_condition",
  "metadata_control_like_clock_candidate", "expression_column", "sample_unit_type",
  "aggregation_id", "aggregation_method", "output_sample_key",
  "include_in_primary_clock_TRUE_FALSE", "manual_approval_TRUE_FALSE", "manual_note"
)

truthy <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

format_table_md <- function(dt) {
  if (is.null(dt) || nrow(dt) == 0) return("_None._")
  dt <- as.data.table(dt)
  header <- paste("|", paste(names(dt), collapse = " | "), "|")
  sep <- paste("|", paste(rep("---", ncol(dt)), collapse = " | "), "|")
  rows <- apply(dt, 1, function(x) paste("|", paste(as.character(x), collapse = " | "), "|"))
  paste(c(header, sep, rows), collapse = "\n")
}

read_matrix_header <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  names(fread(path, nrows = 0))
}

read_gene_ids <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  dt <- fread(path, select = 1, showProgress = FALSE)
  unique(trimws(as.character(dt[[1]])))
}

detect_raw_counts <- function(path) {
  if (is.na(path) || !file.exists(path)) return(FALSE)
  dt <- fread(path, nrows = 200, showProgress = FALSE)
  if (ncol(dt) <= 1) return(FALSE)
  vals <- suppressWarnings(as.numeric(unlist(dt[, -1, with = FALSE], use.names = FALSE)))
  vals <- vals[is.finite(vals)]
  length(vals) > 0 && mean(abs(vals - round(vals)) < 1e-8) > 0.99 && max(vals, na.rm = TRUE) > 50
}

gene_id_symbol_map_cache <- NULL
load_gene_id_symbol_map <- function() {
  if (!is.null(gene_id_symbol_map_cache)) return(gene_id_symbol_map_cache)
  feature_paths <- unique(c(
    Sys.glob(file.path(data_input_dir, "standardized_seurat", "01_count_matrices", "*", "*", "filtered_feature_bc_matrix", "features.tsv.gz")),
    Sys.glob(file.path(data_input_dir, "scRNA_dataset", "*", "cellranger_counts", "*", "outs", "filtered_feature_bc_matrix", "features.tsv.gz"))
  ))
  feature_paths <- feature_paths[file.exists(feature_paths)]
  maps <- list()
  for (feature_path in feature_paths) {
    dt <- tryCatch(fread(cmd = paste("gzip -cd", shQuote(feature_path)), header = FALSE), error = function(e) NULL)
    if (is.null(dt) || ncol(dt) < 2) next
    maps[[length(maps) + 1L]] <- dt[, .(gene_id = na_chr(V1), gene_symbol = na_chr(V2))]
  }
  if (length(maps) == 0) {
    gene_id_symbol_map_cache <<- data.table(gene_id = character(), gene_symbol = character())
    return(gene_id_symbol_map_cache)
  }
  mapping <- rbindlist(maps, fill = TRUE)
  mapping <- mapping[grepl("^ENSMUSG", gene_id) & !is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSMUSG", gene_symbol)]
  mapping <- mapping[, .N, by = .(gene_id, gene_symbol)]
  setorder(mapping, gene_id, -N, gene_symbol)
  gene_id_symbol_map_cache <<- mapping[, .SD[1], by = gene_id][, .(gene_id, gene_symbol)]
  gene_id_symbol_map_cache
}

map_gene_ids_for_overlap <- function(genes) {
  genes <- trimws(as.character(genes))
  if (mean(grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", genes), na.rm = TRUE) > 0.5) {
    genes <- sub("\\.[0-9]+$", "", genes)
  }
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (!any(grepl("^ENSMUSG", genes))) return(genes)
  mapping <- load_gene_id_symbol_map()
  if (nrow(mapping) == 0) return(genes)
  lookup <- mapping$gene_symbol
  names(lookup) <- mapping$gene_id
  mapped <- unname(lookup[genes])
  genes[!is.na(mapped) & mapped != ""] <- mapped[!is.na(mapped) & mapped != ""]
  unique(genes[!is.na(genes) & genes != ""])
}

load_optional_genes <- function(paths) {
  genes <- character()
  for (path in paths[file.exists(paths)]) {
    dt <- fread(path, showProgress = FALSE)
    col <- intersect(c("gene_id", "gene", "feature", names(dt)[1]), names(dt))[1]
    genes <- c(genes, as.character(dt[[col]]))
  }
  unique(trimws(genes[!is.na(genes) & genes != ""]))
}

failures <- character()
warnings <- character()
add_failure <- function(...) failures <<- c(failures, paste0(...))
add_warning <- function(...) warnings <<- c(warnings, paste0(...))

mapping <- fread(mapping_path, sep = "\t")
missing_cols <- setdiff(required_cols, names(mapping))
if (length(missing_cols) > 0) add_failure("Missing required mapping columns: ", paste(missing_cols, collapse = ", "))

validation_path <- file.path(out_dir, "bulk_group_manual_mapping_validation.tsv")
warnings_path <- file.path(out_dir, "bulk_group_manual_mapping_validation_warnings.txt")
summary_path <- file.path(out_dir, "bulk_group_manual_mapping_validation_summary.md")

write_zero_approved_failure <- function() {
  msg <- paste(
    "No rows have both manual_approval_TRUE_FALSE == TRUE and",
    "include_in_primary_clock_TRUE_FALSE == TRUE. Edit",
    "bulk_group_manual_mapping_APPROVED.tsv before validation."
  )
  validation_empty <- data.table(
    check = "approved_included_rows",
    status = "FAIL",
    dataset_folder = NA_character_,
    message = msg
  )
  fwrite(validation_empty, validation_path, sep = "\t")
  writeLines("No approved included rows found in manual mapping.", warnings_path)
  writeLines(c(
    "# Bulk Group Manual Mapping Validation",
    "",
    paste0("- Mapping file: `", mapping_path, "`"),
    "- Status: FAIL",
    "",
    "No approved included rows were found.",
    "",
    msg
  ), summary_path)
  message("APPROVED_MAPPING_VALIDATION_FAIL")
  stop("No approved included rows found in manual mapping.", call. = FALSE)
}

if (length(failures) == 0) {
  mapping[, include_in_primary_clock_TRUE_FALSE := truthy(include_in_primary_clock_TRUE_FALSE)]
  mapping[, manual_approval_TRUE_FALSE := truthy(manual_approval_TRUE_FALSE)]
  mapping[, metadata_control_like_clock_candidate := truthy(metadata_control_like_clock_candidate)]
  mapping[, metadata_age_week_numeric := suppressWarnings(as.numeric(metadata_age_week))]
  if (nrow(mapping[include_in_primary_clock_TRUE_FALSE == TRUE & manual_approval_TRUE_FALSE != TRUE]) > 0) {
    add_failure("Rows with include_in_primary_clock_TRUE_FALSE == TRUE must also have manual_approval_TRUE_FALSE == TRUE.")
  }
  included <- mapping[manual_approval_TRUE_FALSE == TRUE & include_in_primary_clock_TRUE_FALSE == TRUE]
  if (nrow(included) == 0) write_zero_approved_failure()
  if (nrow(included[metadata_control_like_clock_candidate != TRUE]) > 0) add_failure("Included rows must have metadata_control_like_clock_candidate == TRUE.")
  if (nrow(included[sample_unit_type %chin% c("unsafe_unknown") | aggregation_method %chin% c("unsafe_do_not_use")]) > 0) {
    add_failure("Included rows contain unsafe sample_unit_type or aggregation_method.")
  }
  allowed_units <- c("biological_replicate", "technical_run", "group_mean")
  allowed_methods <- c("none_expand_replicates", "sum_raw_counts", "mean_logcpm")
  if (nrow(included[!sample_unit_type %chin% allowed_units]) > 0) add_failure("Included rows have invalid sample_unit_type.")
  if (nrow(included[!aggregation_method %chin% allowed_methods]) > 0) add_failure("Included rows have invalid aggregation_method.")
  if (nrow(included[!is.finite(metadata_age_week_numeric)]) > 0) add_failure("Included rows have missing/non-numeric ages.")
  if (nrow(included[aggregation_method == "mean_logcpm" & !grepl("approved_group_mean_logcpm", manual_note, ignore.case = TRUE)]) > 0) {
    add_failure("mean_logcpm rows require manual_note containing approved_group_mean_logcpm.")
  }
  if (nrow(included[sample_unit_type == "technical_run" & aggregation_method != "sum_raw_counts"]) > 0) {
    add_failure("Technical-run rows must use aggregation_method == sum_raw_counts.")
  }
  if (nrow(included[sample_unit_type == "biological_replicate" & aggregation_method != "none_expand_replicates"]) > 0) {
    add_warning("Biological-replicate rows usually should use none_expand_replicates unless group mean is explicitly justified.")
  }
  missing_agg <- included[is.na(aggregation_id) | trimws(as.character(aggregation_id)) == ""]
  if (nrow(missing_agg) > 0) add_failure("Included rows must have non-empty aggregation_id.")
  output_map <- unique(included[, .(output_sample_key, aggregation_id)])
  if (anyDuplicated(output_map$output_sample_key)) {
    add_failure("output_sample_key is duplicated across independent aggregation_id values.")
  }
  aggregation_output_counts <- included[, .(n_output_sample_keys = uniqueN(output_sample_key)), by = aggregation_id]
  if (nrow(aggregation_output_counts[n_output_sample_keys > 1]) > 0) {
    add_failure("Each aggregation_id must map to exactly one output_sample_key.")
  }
  technical_groups <- included[sample_unit_type == "technical_run",
    .(
      n_expression_columns = uniqueN(expression_column),
      n_metadata_sample_keys = uniqueN(metadata_sample_key),
      n_ages = uniqueN(metadata_age_week_numeric)
    ),
    by = aggregation_id
  ]
  if (nrow(technical_groups[n_metadata_sample_keys > 1 | n_ages > 1]) > 0) {
    add_failure("Technical-run aggregation_id groups must represent one shared biological sample metadata row and age.")
  }
  if (nrow(technical_groups[n_expression_columns < 2]) > 0) {
    add_warning("Some technical_run aggregation_id groups contain only one expression column.")
  }
  gse279446 <- included[dataset_folder == "GSE279446"]
  if (nrow(gse279446) > 0) {
    if (!"gse279446_manual_biological_status" %in% names(mapping)) {
      add_failure("GSE279446 included rows require gse279446_manual_biological_status column.")
    } else if (nrow(gse279446[gse279446_manual_biological_status != "verified_control_like"]) > 0) {
      add_failure("GSE279446 included rows require gse279446_manual_biological_status == verified_control_like.")
    }
    if (nrow(gse279446[!grepl("control|WT|non-perturbed|valid", manual_note, ignore.case = TRUE)]) > 0) {
      add_failure("GSE279446 included rows need manual_note explaining WT/control-like, non-perturbed biological validity.")
    }
  }
}

meta <- fread(metadata_path)
common_genes <- load_optional_genes(file.path(results_dir, "integrated_bulk_pseudobulk_technical_gene_universe.tsv"))
selected_genes <- load_optional_genes(list.files(file.path(results_dir, "selected_genes"), pattern = "selected_genes\\.tsv$", full.names = TRUE))

validation_rows <- list()
if (length(failures) == 0) {
  for (ds in sort(unique(mapping$dataset_folder))) {
    meta_ds <- meta[dataset_folder == ds]
    if (nrow(meta_ds) == 0) {
      add_failure("Mapping references dataset not present in metadata: ", ds)
      next
    }
    log_path <- unique(meta_ds$expression_logcpm_path)[1]
    count_path <- unique(meta_ds$expression_counts_path)[1]
    expr_cols <- read_matrix_header(log_path)[-1]
    included_ds <- mapping[dataset_folder == ds & manual_approval_TRUE_FALSE == TRUE & include_in_primary_clock_TRUE_FALSE == TRUE]
    excluded_ds <- mapping[dataset_folder == ds & include_in_primary_clock_TRUE_FALSE != TRUE]
    missing_expr <- setdiff(unique(included_ds$expression_column), expr_cols)
    if (length(missing_expr) > 0) add_failure(ds, " mapping references missing expression columns: ", paste(missing_expr, collapse = ", "))
    raw_counts_detected <- detect_raw_counts(count_path)
    if (raw_counts_detected && nrow(included_ds[aggregation_method == "mean_logcpm" & !grepl("approved_group_mean_logcpm", manual_note, ignore.case = TRUE)]) > 0) {
      add_failure(ds, " uses mean_logcpm despite raw counts being available without approved_group_mean_logcpm.")
    }
    reused <- included_ds[, .(n_output = uniqueN(output_sample_key), aggregations = paste(sort(unique(aggregation_id)), collapse = ";")), by = expression_column][n_output > 1]
    if (nrow(reused) > 0) add_failure(ds, " expression columns reused across multiple output samples: ", paste(reused$expression_column, collapse = ", "))
    unreviewed_excluded <- excluded_ds[is.na(manual_note) | trimws(manual_note) == ""]
    if (nrow(unreviewed_excluded) > 0) add_failure(ds, " has excluded rows without manual_note.")
    genes <- map_gene_ids_for_overlap(read_gene_ids(log_path))
    common_overlap <- if (length(common_genes) > 0) length(intersect(genes, common_genes)) else NA_integer_
    selected_overlap <- if (length(selected_genes) > 0) length(intersect(genes, selected_genes)) else NA_integer_
    if (nrow(included_ds) > 0 && length(genes) < 5000) add_failure(ds, " has too few genes after cleaning: ", length(genes))
    if (nrow(included_ds) > 0 && is.finite(common_overlap) && common_overlap < 5000) {
      add_failure(ds, " has too few genes overlapping integrated universe: ", common_overlap)
    }
    validation_rows[[length(validation_rows) + 1L]] <- data.table(
      dataset_folder = ds,
      n_mapping_rows = nrow(mapping[dataset_folder == ds]),
      n_included_rows = nrow(included_ds),
      n_output_samples = uniqueN(included_ds$output_sample_key),
      n_expression_columns_consumed = uniqueN(included_ds$expression_column),
      raw_counts_detected = raw_counts_detected,
      aggregation_methods = paste(sort(unique(included_ds$aggregation_method)), collapse = ";"),
      n_genes_after_cleaning = length(genes),
      overlap_with_integrated_gene_universe = common_overlap,
      selected_gene_coverage = if (length(selected_genes) > 0) selected_overlap / length(selected_genes) else NA_real_,
      ready_for_training = nrow(included_ds) > 0
    )
  }
}

validation <- if (length(validation_rows) > 0) {
  rbindlist(validation_rows, fill = TRUE)
} else {
  data.table(
    check = character(),
    status = character(),
    dataset_folder = character(),
    message = character()
  )
}
fwrite(validation, validation_path, sep = "\t")
writeLines(if (length(warnings) > 0) unique(warnings) else "No validation warnings.", warnings_path)

approved_counts <- if (length(failures) == 0) {
  mapping[manual_approval_TRUE_FALSE == TRUE & include_in_primary_clock_TRUE_FALSE == TRUE, .(
    approved_rows = .N,
    output_samples = uniqueN(output_sample_key),
    expression_columns = uniqueN(expression_column)
  ), by = dataset_folder][order(dataset_folder)]
} else data.table()
mode_counts <- if (length(failures) == 0) {
  mapping[manual_approval_TRUE_FALSE == TRUE & include_in_primary_clock_TRUE_FALSE == TRUE, .N, by = .(sample_unit_type, aggregation_method)][order(sample_unit_type, aggregation_method)]
} else data.table()
gse273724_ready <- nrow(validation[dataset_folder == "GSE273724" & ready_for_training == TRUE]) > 0
gse279573_ready <- nrow(validation[dataset_folder == "GSE279573" & ready_for_training == TRUE]) > 0
gse279446_excluded <- if (length(failures) == 0) nrow(mapping[dataset_folder == "GSE279446" & manual_approval_TRUE_FALSE == TRUE & include_in_primary_clock_TRUE_FALSE == TRUE]) == 0 else FALSE
can_train <- length(failures) == 0 && gse273724_ready && gse279573_ready && gse279446_excluded

md <- c(
  "# Bulk Group Manual Mapping Validation",
  "",
  paste0("- Mapping file: `", mapping_path, "`"),
  paste0("- GSE273724 ready: ", gse273724_ready),
  paste0("- GSE279573 ready: ", gse279573_ready),
  paste0("- GSE279446 remains excluded: ", gse279446_excluded),
  paste0("- Full training can be safely run: ", can_train),
  "",
  "## Approved Rows By Dataset",
  format_table_md(approved_counts),
  "",
  "## Mapping Mode Counts",
  format_table_md(mode_counts),
  "",
  "## Gene/Readiness QC",
  format_table_md(validation),
  "",
  "## Failures",
  if (length(failures) > 0) paste0("- ", unique(failures)) else "_None._",
  "",
  "## Warnings",
  if (length(warnings) > 0) paste0("- ", unique(warnings)) else "_None._"
)
writeLines(md, summary_path)

if (!can_train) {
  message("APPROVED_MAPPING_VALIDATION_FAIL")
  stop(paste(unique(c(failures, "Required ready states were not all satisfied.")), collapse = "\n"), call. = FALSE)
}

message("APPROVED_MAPPING_VALIDATION_PASS")
