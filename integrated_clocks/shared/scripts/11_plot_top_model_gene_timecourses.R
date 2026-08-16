#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

ensure_dt <- function(x) {
  if (!data.table::is.data.table(x)) x <- data.table::as.data.table(x)
  x
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  gsub("^_|_$", "", x)
}

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
root_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
results_dir <- normalizePath(Sys.getenv("CLOCK_RESULTS_DIR", unset = file.path(project_dir, "results")), mustWork = FALSE)
out_dir <- file.path(results_dir, "top_gene_timecourses")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_sections <- c(
  elasticnet = file.path(out_dir, "elasticnet"),
  pca_glmnet = file.path(out_dir, "pca_glmnet"),
  cross_model = file.path(out_dir, "cross_model"),
  qc = file.path(out_dir, "qc")
)
invisible(lapply(out_sections, dir.create, recursive = TRUE, showWarnings = FALSE))
out_path <- function(section, filename) file.path(out_sections[[section]], filename)

metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
bulk_group_mapping_path <- Sys.getenv(
  "BULK_GROUP_MAPPING_FILE",
  unset = file.path(results_dir, "alignment_audit", "bulk_group_level_alignment", "bulk_group_manual_mapping_APPROVED.tsv")
)

warnings_vec <- character()
add_warning <- function(...) {
  msg <- paste0(...)
  warnings_vec <<- unique(c(warnings_vec, msg))
  warning(msg, call. = FALSE)
}

stop_missing <- function(path, label) {
  if (!file.exists(path)) stop(label, " missing: ", path, call. = FALSE)
}

elasticnet_artifact_path <- file.path(
  results_dir,
  "model_artifacts/deployment/integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds"
)
elasticnet_final_top30_path <- file.path(
  results_dir,
  "integrated_diagnostics/09_selected_gene_diagnostics/elasticnet_hvg_coefficients/top30_final_elasticnet_hvg_coefficients.tsv"
)
pca_version_b_top30_path <- file.path(
  results_dir,
  "selected_genes/integrated_bulk_pseudobulk_pca_glmnet_versionB_top30_genes.tsv"
)
pca_version_b_all_path <- file.path(
  results_dir,
  "model_artifacts/deployment/pca_glmnet__version_b_age_dataset_penalized/selected_gene_table.tsv"
)
pca_artifact_path <- file.path(
  results_dir,
  "model_artifacts/deployment/integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
)

stop_missing(elasticnet_artifact_path, "ElasticNet deployment artifact")
stop_missing(elasticnet_final_top30_path, "ElasticNet final top30 diagnostic table")
stop_missing(pca_version_b_top30_path, "PCA Version B selected top30 table")
stop_missing(pca_version_b_all_path, "PCA Version B selected-gene table")
stop_missing(pca_artifact_path, "PCA/glmnet deployment artifact")
stop_missing(metadata_path, "Metadata table")

elastic_artifact <- readRDS(elasticnet_artifact_path)
pca_artifact <- readRDS(pca_artifact_path)

validate_artifact <- function(artifact, family, selector) {
  if (!identical(artifact$model_family, family)) {
    stop("Artifact model_family mismatch for ", family, ".", call. = FALSE)
  }
  if (!identical(artifact$selector, selector)) {
    stop("Artifact selector mismatch for ", family, ".", call. = FALSE)
  }
  if (!identical(artifact$artifact_scope, "deployment_all_control_like")) {
    stop("Artifact is not a deployment_all_control_like artifact for ", family, ".", call. = FALSE)
  }
  if (!identical(artifact$expression_scale, "logCPM")) {
    stop("Artifact expression scale is not logCPM for ", family, ".", call. = FALSE)
  }
  ref_meta <- ensure_dt(copy(artifact$reference_metadata))
  if (nrow(ref_meta) == 0 || !"sample_key" %in% names(ref_meta)) {
    stop("Artifact reference metadata is missing sample_key for ", family, ".", call. = FALSE)
  }
  if (!"control_like_clock_candidate" %in% names(ref_meta) || any(ref_meta$control_like_clock_candidate != TRUE)) {
    stop("Artifact reference metadata includes non-control-like samples for ", family, ".", call. = FALSE)
  }
  if ("exclusion_reason" %in% names(ref_meta)) {
    bad <- ref_meta[!is.na(exclusion_reason) & exclusion_reason != "", .N]
    if (bad > 0) stop("Artifact reference metadata includes excluded/intervention samples for ", family, ".", call. = FALSE)
  }
  invisible(ref_meta)
}

elastic_ref_meta <- validate_artifact(elastic_artifact, "elasticnet", "hvg_baseline")
pca_ref_meta <- validate_artifact(pca_artifact, "pca_glmnet", "version_b_age_dataset_penalized")
if (!identical(elastic_ref_meta$sample_key, pca_ref_meta$sample_key)) {
  stop("Deployment reference sample keys differ between retained artifacts.", call. = FALSE)
}
ref_meta <- copy(elastic_ref_meta)
if (nrow(ref_meta) != 117L) {
  add_warning("Expected 117 control-like deployment reference samples; found ", nrow(ref_meta), ".")
}

if (!"biological_age_bin" %in% names(ref_meta)) {
  ref_meta[, biological_age_bin := bio_age_bin %||% NA_character_]
}
if (!"source_batch" %in% names(ref_meta)) {
  ref_meta[, source_batch := interaction(modality, dataset_folder, drop = TRUE)]
}
ref_meta[, age_week := as.numeric(age_week)]
ref_meta[, age_point := round(age_week, 1)]
setkey(ref_meta, sample_key)

select_directional_top <- function(dt, score_col, n = 15L) {
  older <- copy(dt[get(score_col) > 0][order(-abs(get(score_col)))])
  younger <- copy(dt[get(score_col) < 0][order(-abs(get(score_col)))])
  if (nrow(older) < n) add_warning("Fewer than ", n, " older-predictive genes for ", unique(dt$model_family), ".")
  if (nrow(younger) < n) add_warning("Fewer than ", n, " younger-predictive genes for ", unique(dt$model_family), ".")
  older <- head(older, n)
  younger <- head(younger, n)
  older[, rank_within_model_direction := seq_len(.N)]
  younger[, rank_within_model_direction := seq_len(.N)]
  rbindlist(list(older, younger), fill = TRUE)
}

elastic_coef <- ensure_dt(copy(elastic_artifact$model$nonzero_coefficient_table))
if (nrow(elastic_coef) == 0) stop("ElasticNet artifact has no nonzero coefficients.", call. = FALSE)
gene_col <- if ("gene_id" %in% names(elastic_coef)) "gene_id" else if ("gene" %in% names(elastic_coef)) "gene" else names(elastic_coef)[1]
raw_col <- if ("raw_glmnet_gene_coefficient" %in% names(elastic_coef)) {
  "raw_glmnet_gene_coefficient"
} else if ("coefficient" %in% names(elastic_coef)) {
  "coefficient"
} else {
  stop("ElasticNet coefficient table has no recognizable coefficient column.", call. = FALSE)
}
if ("calibrated_gene_coefficient" %in% names(elastic_coef)) {
  elastic_coef[, calibrated_gene_coefficient := as.numeric(calibrated_gene_coefficient)]
  elastic_coef[, raw_glmnet_gene_coefficient := as.numeric(get(raw_col))]
} else {
  elastic_coef[, raw_glmnet_gene_coefficient := as.numeric(get(raw_col))]
  elastic_coef[, calibrated_gene_coefficient := elastic_artifact$calibration$slope * raw_glmnet_gene_coefficient]
}
elastic_top_source <- elastic_coef[, .(
  model_family = "elasticnet",
  selector = "hvg_baseline",
  gene = as.character(get(gene_col)),
  model_direction = fifelse(calibrated_gene_coefficient > 0, "older_predictive", "younger_predictive"),
  model_score = calibrated_gene_coefficient,
  abs_model_score = abs(calibrated_gene_coefficient),
  calibrated_coefficient = calibrated_gene_coefficient,
  raw_glmnet_gene_coefficient = raw_glmnet_gene_coefficient
)][is.finite(model_score) & model_score != 0]

elastic_all <- copy(elastic_top_source[order(-abs_model_score, gene)])
if (anyDuplicated(elastic_all$gene)) {
  stop("ElasticNet deployment artifact contains duplicate non-zero gene coefficients.", call. = FALSE)
}
elastic_all[, coefficient_rank := seq_len(.N)]

select_elasticnet_final_top30 <- function(elastic_source, top30_path) {
  final_top30 <- fread(top30_path)
  required_cols <- c("rank", "gene_id", "coefficient", "abs_coefficient")
  missing_cols <- setdiff(required_cols, names(final_top30))
  if (length(missing_cols) > 0) {
    stop("ElasticNet final top30 table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  final_top30 <- final_top30[order(as.integer(rank))]
  final_top30 <- final_top30[seq_len(min(.N, 30L)), .(
    gene = as.character(gene_id),
    rank_within_model_abs = as.integer(rank),
    final_top30_coefficient = as.numeric(coefficient),
    final_top30_abs_coefficient = as.numeric(abs_coefficient)
  )]
  if (anyDuplicated(final_top30$gene)) {
    stop("ElasticNet final top30 table contains duplicate genes.", call. = FALSE)
  }
  merged <- merge(final_top30, elastic_source, by = "gene", all.x = TRUE, sort = FALSE)
  missing_from_artifact <- merged[is.na(model_score), gene]
  if (length(missing_from_artifact) > 0) {
    stop("Final ElasticNet top30 genes are missing from the deployment artifact: ", paste(missing_from_artifact, collapse = ", "), call. = FALSE)
  }
  merged[, `:=`(
    model_score = final_top30_coefficient,
    abs_model_score = final_top30_abs_coefficient,
    model_direction = fifelse(final_top30_coefficient > 0, "older_predictive", "younger_predictive")
  )]
  merged[, coefficient_plot_direction_order := fifelse(model_score >= 0, 1L, 2L)]
  setorder(merged, coefficient_plot_direction_order, -model_score, gene)
  merged[, coefficient_plot_order := seq_len(.N)]
  merged[, rank_within_model_direction := seq_len(.N), by = model_direction]
  merged[order(rank_within_model_abs)]
}

elastic_top <- select_elasticnet_final_top30(elastic_top_source, elasticnet_final_top30_path)

select_pca_version_b_top30 <- function(top30_path) {
  version_b <- fread(top30_path, sep = "\t")
  required_cols <- c("rank", "gene_id", "age_beta", "feature_score")
  missing_cols <- setdiff(required_cols, names(version_b))
  if (length(missing_cols) > 0) {
    stop("PCA Version B selected top30 table is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  version_b <- version_b[order(as.integer(rank))]
  version_b <- version_b[seq_len(min(.N, 30L)), .(
    gene = as.character(gene_id),
    version_b_rank = as.integer(rank),
    model_family = "pca_glmnet",
    selector = "version_b_age_dataset_penalized",
    model_direction = fifelse(as.numeric(age_beta) > 0, "older_predictive", "younger_predictive"),
    model_score = as.numeric(age_beta),
    abs_model_score = as.numeric(feature_score),
    version_b_age_beta = as.numeric(age_beta),
    version_b_feature_score = as.numeric(feature_score),
    dominant_pc = NA_character_,
    dominant_pc_share = NA_real_,
    mixed_axis_flag = NA
  )]
  if (anyDuplicated(version_b$gene)) {
    stop("PCA Version B selected top30 table contains duplicate genes.", call. = FALSE)
  }
  version_b[, coefficient_plot_direction_order := fifelse(model_score >= 0, 1L, 2L)]
  setorder(version_b, coefficient_plot_direction_order, version_b_rank, gene)
  version_b[, coefficient_plot_order := seq_len(.N)]
  version_b[, rank_within_model_direction := seq_len(.N), by = model_direction]
  version_b[order(version_b_rank)]
}

pca_top <- select_pca_version_b_top30(pca_version_b_top30_path)
fwrite(pca_top, out_path("pca_glmnet", "pca_version_b_top30_genes_for_timecourse.tsv"), sep = "\t")

pca_version_b_all <- fread(pca_version_b_all_path, sep = "\t")
required_version_b_cols <- c(
  "gene_id", "age_beta", "age_t", "partial_r2_age", "partial_r2_dataset",
  "n_samples_used", "feature_score"
)
missing_version_b_cols <- setdiff(required_version_b_cols, names(pca_version_b_all))
if (length(missing_version_b_cols) > 0) {
  stop("PCA Version B selected-gene table is missing columns: ", paste(missing_version_b_cols, collapse = ", "), call. = FALSE)
}
if (nrow(pca_version_b_all) != 2000L || uniqueN(pca_version_b_all$gene_id) != 2000L) {
  stop("Expected 2,000 unique PCA Version B genes.", call. = FALSE)
}
if (!setequal(pca_version_b_all$gene_id, pca_artifact$selected_genes)) {
  stop("PCA Version B selected-gene table does not match the deployment artifact.", call. = FALSE)
}
pca_version_b_all[, gene := as.character(gene_id)]
pca_version_b_all[, version_b_rank := seq_len(.N)]

top_genes <- rbindlist(list(elastic_top, pca_top), use.names = TRUE, fill = TRUE)
setorder(top_genes, model_family, model_direction, rank_within_model_direction)
if (top_genes[, .N, by = model_family][N < 25L, .N] > 0) {
  add_warning("Fewer than 25 genes are available for at least one top-30 model plot.")
}
fwrite(top_genes[model_family == "elasticnet"], out_path("elasticnet", "top_model_genes_for_timecourse.tsv"), sep = "\t")
fwrite(top_genes[model_family == "pca_glmnet"], out_path("pca_glmnet", "top_model_genes_for_timecourse.tsv"), sep = "\t")
fwrite(top_genes, out_path("cross_model", "top_model_genes_for_timecourse.tsv"), sep = "\t")

pca_final_heatmap_genes <- copy(pca_top)
pca_final_heatmap_genes[, row_label := as.character(gene)]

selected_genes <- unique(c(
  elastic_all$gene,
  pca_version_b_all$gene,
  top_genes$gene,
  pca_final_heatmap_genes$gene
))

gene_id_symbol_map_cache <- NULL
load_gene_id_symbol_map <- function() {
  if (!is.null(gene_id_symbol_map_cache)) return(gene_id_symbol_map_cache)
  data_input_dir <- file.path(root_dir, "data_input")
  feature_paths <- unique(c(
    Sys.glob(file.path(data_input_dir, "standardized_seurat", "01_count_matrices", "*", "*", "filtered_feature_bc_matrix", "features.tsv.gz")),
    Sys.glob(file.path(data_input_dir, "scRNA_dataset", "*", "cellranger_counts", "*", "outs", "filtered_feature_bc_matrix", "features.tsv.gz"))
  ))
  maps <- list()
  for (feature_path in feature_paths[file.exists(feature_paths)]) {
    dt <- tryCatch(fread(cmd = paste("gzip -cd", shQuote(feature_path)), header = FALSE), error = function(e) NULL)
    if (is.null(dt) || ncol(dt) < 2) next
    maps[[length(maps) + 1L]] <- dt[, .(gene_id = na_chr(V1), gene_symbol = na_chr(V2))]
  }
  if (length(maps) == 0) {
    gene_id_symbol_map_cache <<- data.table(gene_id = character(), gene_symbol = character())
    return(gene_id_symbol_map_cache)
  }
  mapping <- rbindlist(maps, use.names = TRUE, fill = TRUE)
  mapping <- mapping[grepl("^ENSMUSG", gene_id) & !is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSMUSG", gene_symbol)]
  mapping <- mapping[, .N, by = .(gene_id, gene_symbol)]
  setorder(mapping, gene_id, -N, gene_symbol)
  gene_id_symbol_map_cache <<- mapping[, .SD[1], by = gene_id][, .(gene_id, gene_symbol)]
  gene_id_symbol_map_cache
}

map_ensembl_rows_to_symbols <- function(mat) {
  gene_ids <- rownames(mat)
  if (!any(grepl("^ENSMUSG", gene_ids))) return(mat)
  mapping <- load_gene_id_symbol_map()
  if (nrow(mapping) == 0) return(mat)
  symbol_lookup <- mapping$gene_symbol
  names(symbol_lookup) <- mapping$gene_id
  mapped_symbols <- unname(symbol_lookup[gene_ids])
  mapped <- !is.na(mapped_symbols) & mapped_symbols != ""
  if (!any(mapped)) return(mat)
  gene_ids[mapped] <- mapped_symbols[mapped]
  rownames(mat) <- gene_ids
  mat
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
  if (!anyDuplicated(gene_ids)) return(mat)
  row_mean <- rowMeans(mat, na.rm = TRUE)
  keep_idx <- unlist(tapply(seq_along(gene_ids), gene_ids, function(idx) idx[which.max(row_mean[idx])]), use.names = FALSE)
  mat[sort(keep_idx), , drop = FALSE]
}

read_expression_matrix <- function(path, target_genes = NULL) {
  stop_missing(path, "Expression matrix")
  dt <- fread(path)
  if (ncol(dt) < 2) stop("Expression matrix has fewer than two columns: ", path, call. = FALSE)
  gene_col <- names(dt)[1]
  genes <- as.character(dt[[gene_col]])
  dt[, (gene_col) := NULL]
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- genes
  mat <- map_ensembl_rows_to_symbols(collapse_duplicate_genes(mat))
  mat <- collapse_duplicate_genes(mat)
  if (!is.null(target_genes)) {
    mat <- mat[intersect(target_genes, rownames(mat)), , drop = FALSE]
  }
  mat
}

logcpm_from_counts <- function(count_mat) {
  lib_size <- colSums(count_mat, na.rm = TRUE)
  if (any(!is.finite(lib_size) | lib_size <= 0)) {
    stop("Cannot compute logCPM: one or more mapped count columns has non-positive library size.", call. = FALSE)
  }
  log1p(t(t(count_mat) / lib_size * 1e6))
}

align_direct_dataset <- function(mat, meta_ds_all, reference_keys, modality_name) {
  expr_cols <- colnames(mat)
  candidates <- if (modality_name == "bulk_rnaseq") {
    list(sample_id = meta_ds_all$sample_id, accession = meta_ds_all$accession, sample_key = meta_ds_all$sample_key)
  } else {
    list(sample_id = meta_ds_all$sample_id, sample_key = meta_ds_all$sample_key)
  }
  matched_col <- rep(NA_character_, nrow(meta_ds_all))
  for (candidate_name in names(candidates)) {
    idx <- match(as.character(candidates[[candidate_name]]), expr_cols)
    use <- is.na(matched_col) & !is.na(idx)
    matched_col[use] <- expr_cols[idx[use]]
  }
  if (anyNA(matched_col) && ncol(mat) == nrow(meta_ds_all)) {
    missing_idx <- which(is.na(matched_col))
    matched_col[missing_idx] <- expr_cols[missing_idx]
  }
  keep <- !is.na(matched_col) & meta_ds_all$sample_key %chin% reference_keys
  if (!any(keep)) return(NULL)
  out <- mat[, matched_col[keep], drop = FALSE]
  colnames(out) <- meta_ds_all$sample_key[keep]
  out
}

validate_bulk_group_mapping_columns <- function(mapping) {
  required <- c(
    "dataset_folder", "metadata_sample_key", "metadata_age_week",
    "metadata_control_like_clock_candidate", "expression_column",
    "sample_unit_type", "aggregation_id", "aggregation_method",
    "output_sample_key", "include_in_primary_clock_TRUE_FALSE",
    "manual_approval_TRUE_FALSE"
  )
  missing <- setdiff(required, names(mapping))
  if (length(missing) > 0) {
    stop("Approved bulk group mapping is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

align_grouped_bulk_dataset <- function(meta_ds_all, mapping, reference_keys, target_genes) {
  if (nrow(mapping) == 0) return(NULL)
  validate_bulk_group_mapping_columns(mapping)
  ds <- unique(meta_ds_all$dataset_folder)
  map_ds <- copy(mapping[dataset_folder == ds])
  if (nrow(map_ds) == 0) return(NULL)
  map_ds[, include_in_primary_clock_TRUE_FALSE := truthy(include_in_primary_clock_TRUE_FALSE)]
  map_ds[, manual_approval_TRUE_FALSE := truthy(manual_approval_TRUE_FALSE)]
  map_ds[, metadata_control_like_clock_candidate := truthy(metadata_control_like_clock_candidate)]
  map_ds <- map_ds[
    include_in_primary_clock_TRUE_FALSE == TRUE &
      manual_approval_TRUE_FALSE == TRUE &
      metadata_control_like_clock_candidate == TRUE &
      output_sample_key %chin% reference_keys
  ]
  if (nrow(map_ds) == 0) return(NULL)
  methods <- unique(map_ds$aggregation_method)
  if (length(methods) != 1L) stop("Approved bulk group mapping mixes aggregation methods for ", ds, ".", call. = FALSE)
  method <- methods[[1]]
  source_path <- if (method == "sum_raw_counts") unique(meta_ds_all$expression_counts_path)[1] else unique(meta_ds_all$expression_logcpm_path)[1]
  source_mat <- read_expression_matrix(source_path, target_genes = target_genes)
  missing_cols <- setdiff(unique(map_ds$expression_column), colnames(source_mat))
  if (length(missing_cols) > 0) {
    stop("Approved bulk group mapping references missing expression columns for ", ds, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  grouped <- split(seq_len(nrow(map_ds)), map_ds$aggregation_id)
  out_list <- lapply(grouped, function(idx) {
    cols <- map_ds$expression_column[idx]
    if (method == "sum_raw_counts") rowSums(source_mat[, cols, drop = FALSE], na.rm = TRUE) else rowMeans(source_mat[, cols, drop = FALSE], na.rm = TRUE)
  })
  out <- do.call(cbind, out_list)
  output_keys <- vapply(grouped, function(idx) map_ds$output_sample_key[idx][[1]], character(1))
  colnames(out) <- output_keys
  if (method == "sum_raw_counts") out <- logcpm_from_counts(out)
  out
}

find_saved_expression_matrix <- function() {
  candidates <- c(
    "integrated_control_like_logcpm.rds",
    "integrated_control_like_logcpm.tsv",
    "integrated_bulk_pseudobulk_logcpm.rds",
    "integrated_bulk_pseudobulk_logcpm.tsv",
    "deployment_reference_logcpm.rds",
    "deployment_reference_expression_logcpm.rds"
  )
  roots <- c(results_dir, file.path(results_dir, "model_artifacts"), file.path(results_dir, "model_artifacts", "deployment"))
  paths <- unique(file.path(rep(roots, each = length(candidates)), candidates))
  paths[file.exists(paths)][1] %||% NA_character_
}

load_saved_matrix_if_available <- function(path, reference_keys, target_genes) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  message("Loading saved integrated expression matrix: ", path)
  mat <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else read_expression_matrix(path)
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  if (nrow(mat) < ncol(mat) && length(intersect(reference_keys, rownames(mat))) > length(intersect(reference_keys, colnames(mat)))) {
    mat <- t(mat)
  }
  missing_keys <- setdiff(reference_keys, colnames(mat))
  if (length(missing_keys) > 0) {
    stop("Saved expression matrix is missing deployment reference samples: ", paste(head(missing_keys, 20), collapse = ", "), call. = FALSE)
  }
  mat <- mat[intersect(target_genes, rownames(mat)), reference_keys, drop = FALSE]
  if (any(mat < -1e-6, na.rm = TRUE)) stop("Saved expression matrix appears scaled/centered rather than logCPM.", call. = FALSE)
  mat
}

reconstruct_reference_expression <- function(reference_meta, target_genes) {
  saved_path <- find_saved_expression_matrix()
  saved <- load_saved_matrix_if_available(saved_path, reference_meta$sample_key, target_genes)
  if (!is.null(saved)) return(saved)

  message("No saved integrated logCPM matrix found; reconstructing selected top-gene logCPM values.")
  meta_all <- fread(metadata_path)
  meta_all[, control_like_clock_candidate := truthy(control_like_clock_candidate)]
  mapping <- if (file.exists(bulk_group_mapping_path)) fread(bulk_group_mapping_path) else data.table()
  expr <- matrix(NA_real_, nrow = length(target_genes), ncol = nrow(reference_meta), dimnames = list(target_genes, reference_meta$sample_key))

  for (ds in unique(reference_meta$dataset_folder)) {
    meta_ds_all <- copy(meta_all[dataset_folder == ds])
    ref_ds <- reference_meta[dataset_folder == ds]
    if (nrow(meta_ds_all) == 0) {
      stop("No metadata rows found for deployment dataset ", ds, ".", call. = FALSE)
    }
    logcpm_path <- unique(meta_ds_all$expression_logcpm_path)[1]
    mat <- read_expression_matrix(logcpm_path, target_genes = target_genes)
    direct <- align_direct_dataset(mat, meta_ds_all, ref_ds$sample_key, unique(ref_ds$modality))
    if (!is.null(direct)) {
      expr[rownames(direct), colnames(direct)] <- direct
    }
    missing_keys <- ref_ds$sample_key[colSums(!is.na(expr[, ref_ds$sample_key, drop = FALSE])) == 0]
    if (length(missing_keys) > 0 && unique(ref_ds$modality) == "bulk_rnaseq") {
      grouped <- align_grouped_bulk_dataset(meta_ds_all, mapping, missing_keys, target_genes)
      if (!is.null(grouped)) expr[rownames(grouped), colnames(grouped)] <- grouped
    }
  }

  if (any(expr < -1e-6, na.rm = TRUE)) {
    stop("Reconstructed expression matrix contains negative values and appears scaled rather than logCPM.", call. = FALSE)
  }
  expr
}

expr <- reconstruct_reference_expression(ref_meta, selected_genes)
missing_selected_genes <- selected_genes[rowSums(!is.na(expr)) == 0]
if (length(missing_selected_genes) > 0) {
  add_warning("Selected genes missing from expression matrix: ", paste(missing_selected_genes, collapse = ", "))
}
missing_pca_final_heatmap_genes <- intersect(pca_final_heatmap_genes$gene, missing_selected_genes)
if (length(missing_pca_final_heatmap_genes) > 0) {
  stop(
    "Cannot make the requested PCA-final-top30 heatmap because these genes are missing from expression: ",
    paste(missing_pca_final_heatmap_genes, collapse = ", "),
    call. = FALSE
  )
}
missing_elasticnet_genes <- intersect(elastic_all$gene, missing_selected_genes)
if (length(missing_elasticnet_genes) > 0) {
  stop(
    "Cannot calculate deployment ElasticNet gene-age correlations because these genes are missing from expression: ",
    paste(missing_elasticnet_genes, collapse = ", "),
    call. = FALSE
  )
}
missing_version_b_genes <- intersect(pca_version_b_all$gene, missing_selected_genes)
if (length(missing_version_b_genes) > 0) {
  stop(
    "Cannot calculate PCA Version B gene-age correlations because these genes are missing from expression: ",
    paste(head(missing_version_b_genes, 20), collapse = ", "),
    call. = FALSE
  )
}

safe_cor_test <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  n <- sum(keep)
  if (n < 3L || stats::sd(x[keep]) == 0 || stats::sd(y[keep]) == 0) {
    return(list(n = n, estimate = NA_real_, p = NA_real_))
  }
  test <- suppressWarnings(tryCatch(
    stats::cor.test(x[keep], y[keep], method = method, exact = FALSE),
    error = function(e) NULL
  ))
  if (is.null(test)) return(list(n = n, estimate = NA_real_, p = NA_real_))
  list(n = n, estimate = unname(test$estimate), p = test$p.value)
}

elastic_expr_long <- melt(
  as.data.table(expr[elastic_all$gene, , drop = FALSE], keep.rownames = "gene"),
  id.vars = "gene",
  variable.name = "sample_key",
  value.name = "expression_logCPM"
)
elastic_expr_long <- merge(
  elastic_expr_long,
  ref_meta[, .(sample_key, age_week)],
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)
elastic_correlations <- elastic_expr_long[, {
  pearson <- safe_cor_test(expression_logCPM, age_week, "pearson")
  spearman <- safe_cor_test(expression_logCPM, age_week, "spearman")
  .(
    n_samples = pearson$n,
    pearson_r = pearson$estimate,
    pearson_p = pearson$p,
    spearman_rho = spearman$estimate,
    spearman_p = spearman$p
  )
}, by = gene]
elastic_correlations <- merge(
  elastic_all[, .(
    gene,
    coefficient_rank,
    raw_glmnet_gene_coefficient,
    calibrated_gene_coefficient = calibrated_coefficient,
    model_direction
  )],
  elastic_correlations,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
elastic_correlations[, `:=`(
  pearson_fdr_bh = p.adjust(pearson_p, method = "BH"),
  spearman_fdr_bh = p.adjust(spearman_p, method = "BH"),
  pearson_direction_matches_coefficient = sign(pearson_r) == sign(calibrated_gene_coefficient),
  spearman_direction_matches_coefficient = sign(spearman_rho) == sign(calibrated_gene_coefficient),
  artifact_scope = "deployment_all_control_like",
  model_family = "elasticnet",
  selector = "hvg_baseline",
  expression_scale = "logCPM"
)]
setorder(elastic_correlations, coefficient_rank)
setcolorder(elastic_correlations, c(
  "coefficient_rank", "gene", "raw_glmnet_gene_coefficient", "calibrated_gene_coefficient",
  "model_direction", "n_samples", "pearson_r", "pearson_p", "pearson_fdr_bh",
  "spearman_rho", "spearman_p", "spearman_fdr_bh",
  "pearson_direction_matches_coefficient", "spearman_direction_matches_coefficient",
  "artifact_scope", "model_family", "selector", "expression_scale"
))
fwrite(
  elastic_correlations,
  file.path(results_dir, "selected_genes", "integrated_bulk_pseudobulk_elasticnet_hvg_nonzero_gene_age_correlations.tsv"),
  sep = "\t"
)

version_b_expr_long <- melt(
  as.data.table(expr[pca_version_b_all$gene, , drop = FALSE], keep.rownames = "gene"),
  id.vars = "gene",
  variable.name = "sample_key",
  value.name = "expression_logCPM"
)
version_b_expr_long <- merge(
  version_b_expr_long,
  ref_meta[, .(sample_key, age_week)],
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)
version_b_correlations <- version_b_expr_long[, {
  pearson <- safe_cor_test(expression_logCPM, age_week, "pearson")
  spearman <- safe_cor_test(expression_logCPM, age_week, "spearman")
  .(
    n_samples = pearson$n,
    pearson_r = pearson$estimate,
    pearson_p = pearson$p,
    spearman_rho = spearman$estimate,
    spearman_p = spearman$p
  )
}, by = gene]
version_b_correlations <- merge(
  pca_version_b_all[, .(
    gene,
    version_b_rank,
    age_beta,
    age_t,
    partial_r2_age,
    partial_r2_dataset,
    n_samples_used,
    n_nonzero_variance_datasets,
    coverage_weight,
    feature_score
  )],
  version_b_correlations,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
version_b_correlations[, `:=`(
  pearson_fdr_bh = p.adjust(pearson_p, method = "BH"),
  spearman_fdr_bh = p.adjust(spearman_p, method = "BH"),
  pearson_direction_matches_age_beta = sign(pearson_r) == sign(age_beta),
  spearman_direction_matches_age_beta = sign(spearman_rho) == sign(age_beta),
  artifact_scope = "deployment_all_control_like",
  model_family = "pca_glmnet",
  selector = "version_b_age_dataset_penalized",
  expression_scale = "logCPM"
)]
setorder(version_b_correlations, version_b_rank)
setcolorder(version_b_correlations, c(
  "version_b_rank", "gene", "age_beta", "age_t", "partial_r2_age",
  "partial_r2_dataset", "feature_score", "n_samples_used",
  "n_nonzero_variance_datasets", "coverage_weight", "n_samples",
  "pearson_r", "pearson_p", "pearson_fdr_bh", "spearman_rho", "spearman_p",
  "spearman_fdr_bh", "pearson_direction_matches_age_beta",
  "spearman_direction_matches_age_beta", "artifact_scope", "model_family",
  "selector", "expression_scale"
))
fwrite(
  version_b_correlations,
  file.path(results_dir, "selected_genes", "integrated_bulk_pseudobulk_pca_glmnet_versionB_selected_gene_age_correlations.tsv"),
  sep = "\t"
)

top_genes <- top_genes[!gene %chin% missing_selected_genes]
pca_final_heatmap_genes <- pca_final_heatmap_genes[!gene %chin% missing_selected_genes]
expr <- expr[unique(c(top_genes$gene, pca_final_heatmap_genes$gene)), , drop = FALSE]

sample_expr_long <- melt(
  as.data.table(expr, keep.rownames = "gene"),
  id.vars = "gene",
  variable.name = "sample_key",
  value.name = "expression_logCPM"
)
sample_expr_long <- merge(sample_expr_long, ref_meta[, .(
  sample_key, age_week, age_point, modality, dataset_folder,
  biological_age_bin, split_age_stratum, source_batch
)], by = "sample_key", all.x = TRUE, sort = FALSE)
sample_long <- merge(
  sample_expr_long,
  top_genes[, .(
    model_family, selector, gene, model_direction, model_score,
    abs_model_score, dominant_pc, mixed_axis_flag
  )],
  by = "gene", allow.cartesian = TRUE, sort = FALSE
)
setcolorder(sample_long, c(
  "sample_key", "gene", "model_family", "selector", "age_week", "age_point",
  "modality", "dataset_folder", "biological_age_bin", "expression_logCPM",
  "model_direction", "model_score"
))
sample_level_out <- sample_long[, .(
    sample_key, gene, model_family, selector, age_week, age_point, modality,
    dataset_folder, biological_age_bin, expression_logCPM, model_direction, model_score
  )]
fwrite(sample_level_out[model_family == "elasticnet"], out_path("elasticnet", "top_gene_timecourse_sample_level.tsv"), sep = "\t")
fwrite(sample_level_out[model_family == "pca_glmnet"], out_path("pca_glmnet", "top_gene_timecourse_sample_level.tsv"), sep = "\t")
fwrite(sample_level_out, out_path("cross_model", "top_gene_timecourse_sample_level.tsv"), sep = "\t")

age_stats <- unique(ref_meta[, .(sample_key, age_point, modality, dataset_folder)])
age_stats <- age_stats[, .(
  n_samples = .N,
  n_bulk = sum(modality == "bulk_rnaseq"),
  n_pseudobulk = sum(modality == "scrna_pseudobulk"),
  n_datasets = uniqueN(dataset_folder),
  dataset_list = paste(sort(unique(dataset_folder)), collapse = ";")
), by = age_point]
setorder(age_stats, age_point)

eligible <- age_stats[n_samples >= 3]
use_bins <- FALSE
if (nrow(eligible) >= 4) {
  if (nrow(eligible) > 8L) {
    guide <- c(7, 13, 16, 26, 39, 52, 65, 80)
    selected <- numeric()
    pool <- copy(eligible)
    for (g in guide) {
      if (nrow(pool) == 0) break
      pool[, guide_dist := abs(age_point - g)]
      setorder(pool, guide_dist, -n_samples, age_point)
      chosen <- pool$age_point[1]
      selected <- c(selected, chosen)
      pool <- pool[age_point != chosen]
    }
    selected_age_points <- sort(unique(selected))
    age_stats[, selection_reason := fifelse(age_point %in% selected_age_points, "closest_to_guide_grid_with_n_ge_3", "not_selected_representative_limit")]
  } else {
    selected_age_points <- eligible$age_point
    age_stats[, selection_reason := fifelse(age_point %in% selected_age_points, "n_samples_ge_3", "n_samples_lt_3")]
  }
  age_stats[, selected_for_primary_plot := age_point %in% selected_age_points]
  selected_age_table <- copy(age_stats)
  selected_age_table[, age_bin := NA_character_]
} else {
  use_bins <- TRUE
  ref_meta[, age_bin_for_plot := fifelse(
    age_week >= 7 & age_week < 16, "young_adult",
    fifelse(age_week >= 16 & age_week < 40, "mature_adult",
      fifelse(age_week >= 40 & age_week < 65, "advanced_maternal_age", "post_reproductive_ageing")
    )
  )]
  sample_long <- merge(sample_long, ref_meta[, .(sample_key, age_bin_for_plot)], by = "sample_key", all.x = TRUE)
  sample_expr_long <- merge(sample_expr_long, ref_meta[, .(sample_key, age_bin_for_plot)], by = "sample_key", all.x = TRUE)
  selected_age_table <- unique(ref_meta[, .(sample_key, age_bin = age_bin_for_plot, modality, dataset_folder)])[, .(
    n_samples = .N,
    n_bulk = sum(modality == "bulk_rnaseq"),
    n_pseudobulk = sum(modality == "scrna_pseudobulk"),
    n_datasets = uniqueN(dataset_folder),
    dataset_list = paste(sort(unique(dataset_folder)), collapse = ";"),
    selected_for_primary_plot = TRUE,
    selection_reason = "fallback_biological_age_bin"
  ), by = age_bin]
  selected_age_table[, age_point := NA_real_]
}
if (!use_bins && any(selected_age_table[selected_for_primary_plot == TRUE]$n_samples < 3)) {
  add_warning("At least one selected age point has n_samples < 3.")
}
fwrite(
  selected_age_table[, .(
    age_point, age_bin, n_samples, n_bulk, n_pseudobulk, n_datasets,
    dataset_list, selected_for_primary_plot, selection_reason
  )],
  out_path("qc", "selected_age_points_for_timecourse.tsv"),
  sep = "\t"
)

sample_expr_long[, plot_group := if (use_bins) age_bin_for_plot else as.character(age_point)]
sample_long[, plot_group := if (use_bins) age_bin_for_plot else as.character(age_point)]
selected_groups <- if (use_bins) selected_age_table$age_bin else as.character(selected_age_table[selected_for_primary_plot == TRUE]$age_point)
plot_sample_long <- sample_long[plot_group %chin% selected_groups & !is.na(expression_logCPM)]

summarize_group <- function(dt, value_col = "expression_logCPM") {
  dt[, .(
    mean_logCPM = mean(get(value_col), na.rm = TRUE),
    sd_logCPM = sd(get(value_col), na.rm = TRUE),
    se_logCPM = sd(get(value_col), na.rm = TRUE) / sqrt(.N),
    n_samples = .N,
    n_bulk = sum(modality == "bulk_rnaseq"),
    n_pseudobulk = sum(modality == "scrna_pseudobulk"),
    n_datasets = uniqueN(dataset_folder),
    dataset_list = paste(sort(unique(dataset_folder)), collapse = ";")
  ), by = .(model_family, selector, gene, model_direction, model_score, plot_group)]
}

age_mean <- summarize_group(plot_sample_long)
age_mean[, `:=`(
  age_point = suppressWarnings(as.numeric(plot_group)),
  age_bin = if (use_bins) plot_group else NA_character_
)]
fwrite(age_mean[, .(
  model_family, selector, gene, model_direction, model_score, age_point, age_bin,
  mean_logCPM, sd_logCPM, se_logCPM, n_samples, n_bulk, n_pseudobulk, n_datasets, dataset_list
)], out_path("cross_model", "top_gene_timecourse_age_mean.tsv"), sep = "\t")
fwrite(age_mean[model_family == "elasticnet"], out_path("elasticnet", "top_gene_timecourse_age_mean.tsv"), sep = "\t")
fwrite(age_mean[model_family == "pca_glmnet"], out_path("pca_glmnet", "top_gene_timecourse_age_mean.tsv"), sep = "\t")

pca_final_heatmap_sample <- merge(
  sample_expr_long[gene %chin% pca_final_heatmap_genes$gene & plot_group %chin% selected_groups & !is.na(expression_logCPM)],
  pca_final_heatmap_genes[, .(
    gene, version_b_rank, model_family, selector, model_direction, model_score,
    abs_model_score, dominant_pc, dominant_pc_share
  )],
  by = "gene",
  all.x = FALSE,
  all.y = FALSE,
  sort = FALSE
)
pca_final_heatmap_age_mean <- summarize_group(pca_final_heatmap_sample)
pca_final_heatmap_age_mean[, `:=`(
  age_point = suppressWarnings(as.numeric(plot_group)),
  age_bin = if (use_bins) plot_group else NA_character_
)]
pca_final_heatmap_age_mean <- merge(
  pca_final_heatmap_age_mean,
  pca_final_heatmap_genes[, .(gene, version_b_rank, abs_model_score, dominant_pc, dominant_pc_share)],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)
fwrite(
  pca_final_heatmap_age_mean[, .(
    version_b_rank, gene, model_family, selector, model_direction, model_score,
    abs_model_score, dominant_pc, dominant_pc_share, age_point, age_bin,
    mean_logCPM, sd_logCPM, se_logCPM, n_samples, n_bulk, n_pseudobulk,
    n_datasets, dataset_list
  )],
  out_path("pca_glmnet", "pca_version_b_top_gene_timecourse_zscore_heatmap_source.tsv"),
  sep = "\t"
)

dataset_first <- plot_sample_long[, .(
  dataset_mean_logCPM = mean(expression_logCPM, na.rm = TRUE),
  dataset_n_samples = .N
), by = .(model_family, selector, gene, model_direction, model_score, plot_group, dataset_folder)]
dataset_balanced <- dataset_first[, .(
  dataset_balanced_mean_logCPM = mean(dataset_mean_logCPM, na.rm = TRUE),
  sd_dataset_means = sd(dataset_mean_logCPM, na.rm = TRUE),
  se_dataset_means = sd(dataset_mean_logCPM, na.rm = TRUE) / sqrt(.N),
  n_datasets = .N,
  dataset_list = paste(sort(dataset_folder), collapse = ";"),
  n_samples = sum(dataset_n_samples)
), by = .(model_family, selector, gene, model_direction, model_score, plot_group)]
dataset_balanced[, `:=`(age_point = suppressWarnings(as.numeric(plot_group)), age_bin = if (use_bins) plot_group else NA_character_)]
fwrite(dataset_balanced, out_path("cross_model", "top_gene_timecourse_dataset_balanced_mean.tsv"), sep = "\t")
fwrite(dataset_balanced[model_family == "elasticnet"], out_path("elasticnet", "top_gene_timecourse_dataset_balanced_mean.tsv"), sep = "\t")
fwrite(dataset_balanced[model_family == "pca_glmnet"], out_path("pca_glmnet", "top_gene_timecourse_dataset_balanced_mean.tsv"), sep = "\t")

modality_first <- plot_sample_long[, .(
  modality_mean_logCPM = mean(expression_logCPM, na.rm = TRUE),
  modality_n_samples = .N
), by = .(model_family, selector, gene, model_direction, model_score, plot_group, modality)]
modality_balanced <- modality_first[, .(
  modality_balanced_mean_logCPM = mean(modality_mean_logCPM, na.rm = TRUE),
  sd_modality_means = sd(modality_mean_logCPM, na.rm = TRUE),
  se_modality_means = sd(modality_mean_logCPM, na.rm = TRUE) / sqrt(.N),
  n_modalities = .N,
  modality_list = paste(sort(modality), collapse = ";"),
  n_samples = sum(modality_n_samples)
), by = .(model_family, selector, gene, model_direction, model_score, plot_group)]
modality_balanced[, `:=`(age_point = suppressWarnings(as.numeric(plot_group)), age_bin = if (use_bins) plot_group else NA_character_)]
fwrite(modality_balanced, out_path("cross_model", "top_gene_timecourse_modality_balanced_mean.tsv"), sep = "\t")
fwrite(modality_balanced[model_family == "elasticnet"], out_path("elasticnet", "top_gene_timecourse_modality_balanced_mean.tsv"), sep = "\t")
fwrite(modality_balanced[model_family == "pca_glmnet"], out_path("pca_glmnet", "top_gene_timecourse_modality_balanced_mean.tsv"), sep = "\t")

lm_stats <- function(df, formula_obj) {
  fit <- tryCatch(lm(formula_obj, data = df), error = function(e) NULL)
  if (is.null(fit)) return(list(beta = NA_real_, p = NA_real_, r = NA_real_))
  coefs <- summary(fit)$coefficients
  beta <- if ("age_week" %in% rownames(coefs)) coefs["age_week", "Estimate"] else NA_real_
  p <- if ("age_week" %in% rownames(coefs)) coefs["age_week", "Pr(>|t|)"] else NA_real_
  r <- suppressWarnings(cor(df$age_week, df$expression_logCPM, use = "complete.obs", method = "pearson"))
  list(beta = beta, p = p, r = r)
}

trend_summary <- sample_long[!is.na(expression_logCPM), {
  df <- as.data.frame(.SD)
  sample_fit <- lm_stats(df, expression_logCPM ~ age_week)
  dataset_fit <- if (uniqueN(dataset_folder) > 1L) lm_stats(df, expression_logCPM ~ age_week + dataset_folder) else list(beta = NA_real_, p = NA_real_, r = NA_real_)
  batch_fit <- if ("source_batch" %in% names(df) && length(unique(df$source_batch)) > 1L) {
    lm_stats(df, expression_logCPM ~ age_week + source_batch)
  } else {
    list(beta = NA_real_, p = NA_real_, r = NA_real_)
  }
  rho <- suppressWarnings(cor(df$age_week, df$expression_logCPM, use = "complete.obs", method = "spearman"))
  .(
    n_samples = .N,
    sample_level_age_beta = sample_fit$beta,
    sample_level_age_p = sample_fit$p,
    sample_level_age_r = sample_fit$r,
    dataset_adjusted_age_beta = dataset_fit$beta,
    dataset_adjusted_age_p = dataset_fit$p,
    dataset_adjusted_age_r = dataset_fit$r,
    source_batch_adjusted_age_beta = batch_fit$beta,
    source_batch_adjusted_age_p = batch_fit$p,
    source_batch_adjusted_age_r = batch_fit$r,
    spearman_rho = rho
  )
}, by = .(model_family, selector, gene, model_direction, model_score)]
trend_summary[, direction_match_model_score := fifelse(
  model_direction == "older_predictive",
  sample_level_age_beta > 0,
  sample_level_age_beta < 0
)]
fwrite(trend_summary, out_path("cross_model", "top_gene_timecourse_trend_summary.tsv"), sep = "\t")
fwrite(trend_summary[model_family == "elasticnet"], out_path("elasticnet", "top_gene_timecourse_trend_summary.tsv"), sep = "\t")
fwrite(trend_summary[model_family == "pca_glmnet"], out_path("pca_glmnet", "top_gene_timecourse_trend_summary.tsv"), sep = "\t")

theme_set(theme_bw(base_size = 14))
direction_cols <- c(older_predictive = "#8C2D04", younger_predictive = "#225EA8")
plot_manifest <- list()
save_plot <- function(plot, stem, width = 12, height = 8, section = "cross_model") {
  png_path <- out_path(section, paste0(stem, ".png"))
  pdf_path <- out_path(section, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, bg = "white")
  plot_manifest[[length(plot_manifest) + 1L]] <<- data.table(
    section = section,
    plot_id = stem,
    png_path = png_path,
    pdf_path = pdf_path,
    width_in = width,
    height_in = height
  )
}

label_map <- if (use_bins) {
  setNames(selected_age_table$n_samples, selected_age_table$age_bin)
} else {
  setNames(selected_age_table[selected_for_primary_plot == TRUE]$n_samples, as.character(selected_age_table[selected_for_primary_plot == TRUE]$age_point))
}
format_group_label <- function(x) paste0(x, "\nn=", label_map[as.character(x)])

plot_timecourse <- function(model_family, stem, title) {
  mf <- model_family
  dt <- copy(age_mean[model_family == mf])
  if (nrow(dt) == 0) return(invisible(NULL))
  if (model_family == "pca_glmnet") {
    pc_lookup <- top_genes[model_family == "pca_glmnet", .(gene, dominant_pc)]
    dt <- merge(dt, pc_lookup, by = "gene", all.x = TRUE)
    dt[, facet_gene := fifelse(!is.na(dominant_pc) & dominant_pc != "", paste0(gene, " (", dominant_pc, ")"), gene)]
  } else {
    dt[, facet_gene := gene]
  }
  dt[, x_group := if (use_bins) factor(plot_group, levels = selected_groups) else as.numeric(plot_group)]
  p <- ggplot(dt, aes(x = x_group, y = mean_logCPM, group = gene, color = model_direction)) +
    geom_line(linewidth = 0.45) +
    geom_point(size = 1.7) +
    geom_errorbar(aes(ymin = mean_logCPM - se_logCPM, ymax = mean_logCPM + se_logCPM), width = if (use_bins) 0.15 else 0.8, linewidth = 0.35) +
    facet_wrap(~facet_gene, scales = "free_y", ncol = 5) +
    scale_color_manual(values = direction_cols) +
    labs(title = title, x = if (use_bins) "Biological age bin" else "Age (weeks)", y = "Mean logCPM +/- SE", color = "Model direction") +
    theme(strip.text = element_text(size = 8), axis.text.x = element_text(size = 5.5))
  if (!use_bins) {
    p <- p + scale_x_continuous(breaks = as.numeric(selected_groups), labels = format_group_label, guide = guide_axis(n.dodge = 2))
  } else {
    p <- p + scale_x_discrete(labels = format_group_label, guide = guide_axis(n.dodge = 2))
  }
  save_plot(p, stem, width = 15, height = 11, section = if (model_family == "pca_glmnet") "pca_glmnet" else "elasticnet")
}

plot_timecourse("elasticnet", "01_elasticnet_top_gene_timecourse_mean_logcpm", "ElasticNet top model-gene trajectories")
plot_timecourse("pca_glmnet", "01_pca_version_b_top_gene_timecourse_mean_logcpm", "PCA/glmnet Version B selected-gene trajectories")

zscore_dt <- function(dt) {
  out <- copy(dt)
  out[, z_mean_logCPM := {
    z <- as.numeric(scale(mean_logCPM))
    ifelse(is.finite(z), z, 0)
  }, by = .(model_family, gene)]
  out
}

plot_heatmap <- function(model_family, stem, title) {
  mf <- model_family
  dt <- zscore_dt(age_mean[model_family == mf])
  if (nrow(dt) == 0) return(invisible(NULL))
  order_dt <- top_genes[model_family == mf]
  if (model_family == "elasticnet" && "coefficient_plot_order" %in% names(order_dt) && any(!is.na(order_dt$coefficient_plot_order))) {
    setorder(order_dt, coefficient_plot_order)
  } else {
    setorder(order_dt, model_direction, -abs_model_score)
  }
  if (model_family == "pca_glmnet") {
    order_dt[, row_label := fifelse(!is.na(dominant_pc) & dominant_pc != "", paste0(gene, " (", dominant_pc, fifelse(mixed_axis_flag == TRUE, "; mixed", ""), ")"), gene)]
  } else {
    order_dt[, row_label := gene]
  }
  dt <- merge(dt, order_dt[, .(gene, row_label)], by = "gene", all.x = TRUE)
  dt[, row_label := factor(row_label, levels = rev(order_dt$row_label))]
  dt[, x_group := factor(plot_group, levels = selected_groups)]
  p <- ggplot(dt, aes(x = x_group, y = row_label, fill = z_mean_logCPM)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Row z") +
    labs(title = title, x = if (use_bins) "Biological age bin" else "Age point", y = NULL) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
  save_plot(p, stem, width = 9, height = 9, section = if (model_family == "pca_glmnet") "pca_glmnet" else "elasticnet")
}

plot_ordered_heatmap <- function(dt, order_dt, stem, title, subtitle = NULL) {
  dt <- zscore_dt(dt)
  if (nrow(dt) == 0) return(invisible(NULL))
  order_dt <- copy(order_dt)
  setorder(order_dt, version_b_rank)
  dt <- merge(dt, order_dt[, .(gene, row_label)], by = "gene", all.x = FALSE, sort = FALSE)
  dt[, row_label := factor(row_label, levels = rev(order_dt$row_label))]
  dt[, x_group := factor(plot_group, levels = selected_groups)]
  p <- ggplot(dt, aes(x = x_group, y = row_label, fill = z_mean_logCPM)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Row z") +
    labs(title = title, subtitle = subtitle, x = if (use_bins) "Biological age bin" else "Age point", y = NULL) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7),
      plot.subtitle = element_text(size = 9, margin = margin(t = 2, b = 6))
    )
  save_plot(p, stem, width = 9, height = 9, section = "pca_glmnet")
}

plot_heatmap("elasticnet", "02_elasticnet_top_gene_timecourse_zscore_heatmap", "ElasticNet top-gene z-scored mean logCPM")
plot_heatmap("pca_glmnet", "02_pca_version_b_top_gene_timecourse_zscore_heatmap", "PCA/glmnet Version B selected-gene z-scored mean logCPM")

p5 <- ggplot(trend_summary, aes(x = model_direction, fill = direction_match_model_score)) +
  geom_bar(position = "dodge", color = "grey25", linewidth = 0.25) +
  facet_wrap(~model_family) +
  scale_fill_manual(values = c(`TRUE` = "#238B45", `FALSE` = "#CB181D"), na.value = "grey70") +
  labs(title = "Observed age trend direction vs model direction", x = NULL, y = "Selected genes", fill = "Direction match") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_plot(p5, "01_top_gene_timecourse_direction_summary", width = 8, height = 5, section = "cross_model")

comp <- unique(ref_meta[, .(sample_key, plot_group = if (use_bins) age_bin_for_plot else as.character(age_point), modality, dataset_folder)])
comp <- comp[plot_group %chin% selected_groups]
comp[, x_group := factor(plot_group, levels = selected_groups)]
p6 <- ggplot(comp, aes(x = x_group, fill = modality)) +
  geom_bar(color = "grey25", linewidth = 0.2) +
  labs(title = "Selected age-point sample composition", x = if (use_bins) "Biological age bin" else "Age point", y = "Samples", fill = "Modality") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p6, "01_age_point_sample_composition", width = 8, height = 5, section = "qc")

p6b <- ggplot(comp, aes(x = x_group, fill = dataset_folder)) +
  geom_bar(color = "grey25", linewidth = 0.15) +
  labs(title = "Selected age-point dataset composition", x = if (use_bins) "Biological age bin" else "Age point", y = "Samples", fill = "Dataset") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.text = element_text(size = 6), legend.key.size = grid::unit(0.25, "cm"))
save_plot(p6b, "02_age_point_dataset_composition", width = 10, height = 6, section = "qc")

dataset_sens <- merge(
  age_mean[, .(model_family, selector, gene, model_direction, plot_group, sample_mean_logCPM = mean_logCPM)],
  dataset_balanced[, .(model_family, selector, gene, plot_group, dataset_balanced_mean_logCPM)],
  by = c("model_family", "selector", "gene", "plot_group")
)
p7 <- ggplot(dataset_sens, aes(x = sample_mean_logCPM, y = dataset_balanced_mean_logCPM, color = model_direction)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.75, size = 1.5) +
  facet_wrap(~model_family, scales = "free") +
  scale_color_manual(values = direction_cols) +
  labs(title = "Dataset-balanced sensitivity", x = "Sample-weighted mean logCPM", y = "Dataset-balanced mean logCPM", color = "Direction")
save_plot(p7, "02_dataset_balanced_vs_sample_mean_sensitivity", width = 8, height = 5, section = "cross_model")

modality_sens <- merge(
  age_mean[, .(model_family, selector, gene, model_direction, plot_group, sample_mean_logCPM = mean_logCPM)],
  modality_balanced[, .(model_family, selector, gene, plot_group, modality_balanced_mean_logCPM)],
  by = c("model_family", "selector", "gene", "plot_group")
)
p8 <- ggplot(modality_sens, aes(x = sample_mean_logCPM, y = modality_balanced_mean_logCPM, color = model_direction)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.75, size = 1.5) +
  facet_wrap(~model_family, scales = "free") +
  scale_color_manual(values = direction_cols) +
  labs(title = "Modality-balanced sensitivity", x = "Sample-weighted mean logCPM", y = "Modality-balanced mean logCPM", color = "Direction")
save_plot(p8, "03_modality_balanced_vs_sample_mean_sensitivity", width = 8, height = 5, section = "cross_model")

overlap_genes <- top_genes[, uniqueN(model_family), by = gene][V1 > 1, gene]
combined <- zscore_dt(age_mean)
combined_key <- top_genes[, .(model_family, gene, model_direction, abs_model_score)]
combined_key[, row_label := paste0(model_family, " | ", model_direction, " | ", gene, fifelse(gene %chin% overlap_genes, " *", ""))]
setorder(combined_key, model_family, model_direction, -abs_model_score)
combined <- merge(combined, combined_key[, .(model_family, gene, row_label)], by = c("model_family", "gene"))
combined[, row_label := factor(row_label, levels = rev(combined_key$row_label))]
combined[, x_group := factor(plot_group, levels = selected_groups)]
p9 <- ggplot(combined, aes(x = x_group, y = row_label, fill = z_mean_logCPM)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Row z") +
  labs(title = "Combined retained-model top-gene heatmap", subtitle = "* gene selected by both retained models", x = if (use_bins) "Biological age bin" else "Age point", y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 6))
save_plot(p9, "04_combined_model_top_gene_heatmap", width = 11, height = 12, section = "cross_model")

candidate_genes <- top_genes[, head(.SD[order(-abs_model_score)], 6L), by = .(model_family, model_direction)]
candidate_mean <- merge(age_mean, candidate_genes[, .(model_family, gene)], by = c("model_family", "gene"))
candidate_mean[, facet_gene := paste0(model_family, "\n", gene)]
candidate_mean[, x_group := if (use_bins) factor(plot_group, levels = selected_groups) else as.numeric(plot_group)]
p10 <- ggplot(candidate_mean, aes(x = x_group, y = mean_logCPM, color = model_direction, group = gene)) +
  geom_line(linewidth = 0.55) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = mean_logCPM - se_logCPM, ymax = mean_logCPM + se_logCPM), width = if (use_bins) 0.15 else 0.8, linewidth = 0.35) +
  facet_wrap(~facet_gene, scales = "free_y", ncol = 4) +
  scale_color_manual(values = direction_cols) +
  labs(title = "Publication-candidate top model-gene trajectories", x = if (use_bins) "Biological age bin" else "Age (weeks)", y = "Mean logCPM +/- SE", color = "Direction") +
  theme(strip.text = element_text(size = 8), axis.text.x = element_text(size = 6))
if (!use_bins) {
  p10 <- p10 + scale_x_continuous(breaks = as.numeric(selected_groups), labels = format_group_label, guide = guide_axis(n.dodge = 2))
} else {
  p10 <- p10 + scale_x_discrete(labels = format_group_label, guide = guide_axis(n.dodge = 2))
}
save_plot(p10, "05_publication_candidate_top12_timecourses", width = 14, height = 11, section = "cross_model")

manifest_dt <- rbindlist(plot_manifest, fill = TRUE)
fwrite(manifest_dt, out_path("qc", "timecourse_plot_manifest.tsv"), sep = "\t")

dataset_sens_delta <- dataset_sens[, abs(sample_mean_logCPM - dataset_balanced_mean_logCPM)]
modality_sens_delta <- modality_sens[, abs(sample_mean_logCPM - modality_balanced_mean_logCPM)]
trend_counts <- trend_summary[, .(
  n_match = sum(direction_match_model_score == TRUE, na.rm = TRUE),
  n_mismatch = sum(direction_match_model_score == FALSE, na.rm = TRUE),
  n_unknown = sum(is.na(direction_match_model_score))
), by = .(model_family, model_direction)]
age_count_lines <- selected_age_table[selected_for_primary_plot == TRUE, paste0(
  ifelse(is.na(age_point), age_bin, age_point),
  " (n=", n_samples, ")"
)]
top_count_lines <- top_genes[, .N, by = .(model_family, model_direction)][order(model_family, model_direction)]
summary_lines <- c(
  "# Top Model Gene Time-Course Summary",
  "",
  "## Analysis identity",
  "Top model gene expression time courses from the integrated deployment reference, using control-like chronological samples only.",
  "",
  "## Model genes used",
  "- ElasticNet / hvg_baseline",
  "- PCA/glmnet / version_b_age_dataset_penalized",
  "",
  "## Selected top genes",
  paste0("- ", top_count_lines$model_family, " / ", top_count_lines$model_direction, ": ", top_count_lines$N),
  "",
  "## Expression scale",
  "Expression values are logCPM values, not scaled model inputs.",
  "",
  "## Age-point averaging rule",
  "Samples were averaged at selected representative age points; individual sample points are not shown in the main time-course panels.",
  paste0("- Selected ", if (use_bins) "biological bins" else "age points", ": ", paste(age_count_lines, collapse = "; ")),
  "",
  "## Direction consistency",
  paste0("- ", trend_counts$model_family, " / ", trend_counts$model_direction, ": ", trend_counts$n_match, " matched, ", trend_counts$n_mismatch, " mismatched, ", trend_counts$n_unknown, " unknown"),
  "",
  "## Dataset/modality sensitivity",
  paste0("- Median absolute dataset-balanced vs sample-weighted difference: ", signif(stats::median(dataset_sens_delta, na.rm = TRUE), 3), " logCPM."),
  paste0("- Median absolute modality-balanced vs sample-weighted difference: ", signif(stats::median(modality_sens_delta, na.rm = TRUE), 3), " logCPM."),
  "",
  "## Suggested publication use",
  "- Main/supplemental candidate: cross_model/05_publication_candidate_top12_timecourses.",
  "- Supplemental context: model-level top-30 time courses, z-score heatmaps, direction summary, composition, and sensitivity plots.",
  "",
  "## Caution",
  "These plots show descriptive expression trajectories of model-selected genes across the control-like chronological reference cohort. They do not establish that individual genes are causal drivers of ovarian ageing."
)
writeLines(summary_lines, out_path("cross_model", "top_gene_timecourse_summary.md"))

if (length(warnings_vec) == 0) warnings_vec <- "No warnings."
writeLines(warnings_vec, out_path("qc", "timecourse_warnings.txt"))

message("Output directory: ", out_dir)
message("Reference samples used: ", nrow(ref_meta))
message("Selected age groups: ", paste(age_count_lines, collapse = "; "))
message("Top genes per model: ", paste(top_count_lines$model_family, top_count_lines$model_direction, top_count_lines$N, sep = "/", collapse = "; "))
message("Plots written: ", nrow(manifest_dt), " PNG and ", nrow(manifest_dt), " PDF files")
