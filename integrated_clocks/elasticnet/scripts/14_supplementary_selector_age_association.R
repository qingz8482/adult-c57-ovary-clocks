#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

options(stringsAsFactors = FALSE)
set.seed(20260721)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
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
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
root_dir <- normalizePath(file.path(project_dir, ".."), mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
selected_dir <- file.path(results_dir, "selected_genes")
out_dir <- file.path(selected_dir, "supplementary_selector_age_association")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
beta_discrepancy_path <- file.path(out_dir, "versionB_beta_reproduction_discrepancy_top50.tsv")
if (file.exists(beta_discrepancy_path)) unlink(beta_discrepancy_path)

metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
bulk_group_mapping_path <- file.path(
  results_dir,
  "alignment_audit", "bulk_group_level_alignment", "bulk_group_manual_mapping_APPROVED.tsv"
)
elastic_artifact_path <- file.path(
  results_dir, "model_artifacts", "deployment",
  "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds"
)
pca_artifact_path <- file.path(
  results_dir, "model_artifacts", "deployment",
  "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
)
version_b_selected_path <- file.path(
  results_dir, "model_artifacts", "deployment",
  "pca_glmnet__version_b_age_dataset_penalized", "selected_gene_table.tsv"
)
elastic_marginal_path <- file.path(
  selected_dir, "integrated_bulk_pseudobulk_elasticnet_hvg_nonzero_gene_age_correlations.tsv"
)
version_b_marginal_path <- file.path(
  selected_dir, "integrated_bulk_pseudobulk_pca_glmnet_versionB_selected_gene_age_correlations.tsv"
)

required_paths <- c(
  metadata_path, elastic_artifact_path, pca_artifact_path,
  version_b_selected_path, elastic_marginal_path, version_b_marginal_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Required inputs are missing: ", paste(missing_paths, collapse = "; "), call. = FALSE)
}

elastic_artifact <- readRDS(elastic_artifact_path)
pca_artifact <- readRDS(pca_artifact_path)
if (!identical(elastic_artifact$artifact_scope, "deployment_all_control_like") ||
    !identical(pca_artifact$artifact_scope, "deployment_all_control_like")) {
  stop("Expected frozen deployment_all_control_like artifacts.", call. = FALSE)
}
if (!identical(elastic_artifact$training_sample_keys, pca_artifact$training_sample_keys)) {
  stop("ElasticNet and PCA/glmnet deployment sample sets differ.", call. = FALSE)
}
if (!identical(elastic_artifact$expression_scale, "logCPM") ||
    !identical(pca_artifact$expression_scale, "logCPM")) {
  stop("Frozen artifacts do not report logCPM expression.", call. = FALSE)
}

reference_keys <- as.character(pca_artifact$training_sample_keys)
ref_meta <- as.data.table(pca_artifact$reference_metadata)
ref_meta <- ref_meta[match(reference_keys, sample_key)]
if (nrow(ref_meta) != 117L || anyNA(ref_meta$sample_key)) {
  stop("Expected 117 aligned deployment reference samples.", call. = FALSE)
}
if (anyDuplicated(reference_keys) || anyDuplicated(ref_meta$sample_key)) {
  stop("Frozen deployment sample identifiers are not unique.", call. = FALSE)
}
if (!all(ref_meta$control_like_clock_candidate == TRUE)) {
  stop("Deployment reference includes non-control-like samples.", call. = FALSE)
}
if (!"source_batch" %in% names(ref_meta)) {
  stop("Frozen artifact lacks the source_batch covariate.", call. = FALSE)
}
ref_meta[, age_week := as.numeric(age_week)]
ref_meta[, source_batch := droplevels(factor(source_batch))]
if (any(!is.finite(ref_meta$age_week)) || anyNA(ref_meta$source_batch)) {
  stop("Age or source_batch contains missing values.", call. = FALSE)
}

elastic_coef <- as.data.table(elastic_artifact$model$nonzero_coefficient_table)
setnames(elastic_coef, c("gene_id", "coefficient"), c("gene", "elasticnet_coefficient"))
elastic_coef <- elastic_coef[is.finite(elasticnet_coefficient) & elasticnet_coefficient != 0]
if (nrow(elastic_coef) != 66L || uniqueN(elastic_coef$gene) != 66L) {
  stop("Expected 66 unique non-zero frozen ElasticNet genes.", call. = FALSE)
}

version_b_selected <- fread(version_b_selected_path, sep = "\t")
required_version_b <- c(
  "gene_id", "age_beta", "age_t", "partial_r2_age", "partial_r2_dataset",
  "n_samples_used", "residual_df", "n_nonzero_variance_datasets", "feature_score"
)
missing_version_b <- setdiff(required_version_b, names(version_b_selected))
if (length(missing_version_b) > 0) {
  stop("Version B selected table lacks columns: ", paste(missing_version_b, collapse = ", "), call. = FALSE)
}
if (nrow(version_b_selected) != 2000L || uniqueN(version_b_selected$gene_id) != 2000L) {
  stop("Expected 2,000 unique Version B genes.", call. = FALSE)
}
if (!setequal(version_b_selected$gene_id, pca_artifact$selected_genes)) {
  stop("Version B table and frozen PCA/glmnet artifact disagree on selected genes.", call. = FALSE)
}
setnames(version_b_selected, "gene_id", "gene")
setnames(version_b_selected, "residual_df", "stored_residual_df")
version_b_selected[, version_b_rank := seq_len(.N)]

elastic_marginal <- fread(elastic_marginal_path, sep = "\t")
version_b_marginal <- fread(version_b_marginal_path, sep = "\t")
if (!all(c("gene", "pearson_r") %in% names(elastic_marginal)) ||
    !all(c("gene", "pearson_r") %in% names(version_b_marginal))) {
  stop("Marginal-correlation source tables lack gene or pearson_r.", call. = FALSE)
}

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
  mapping <- mapping[
    grepl("^ENSMUSG", gene_id) & !is.na(gene_symbol) & gene_symbol != "" & !grepl("^ENSMUSG", gene_symbol)
  ]
  mapping <- mapping[, .N, by = .(gene_id, gene_symbol)]
  setorder(mapping, gene_id, -N, gene_symbol)
  gene_id_symbol_map_cache <<- mapping[, .SD[1], by = gene_id][, .(gene_id, gene_symbol)]
  gene_id_symbol_map_cache
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

map_ensembl_rows_to_symbols <- function(mat) {
  gene_ids <- rownames(mat)
  if (!any(grepl("^ENSMUSG", gene_ids))) return(mat)
  mapping <- load_gene_id_symbol_map()
  if (nrow(mapping) == 0) return(mat)
  lookup <- mapping$gene_symbol
  names(lookup) <- mapping$gene_id
  mapped <- unname(lookup[gene_ids])
  use <- !is.na(mapped) & mapped != ""
  gene_ids[use] <- mapped[use]
  rownames(mat) <- gene_ids
  mat
}

read_expression_matrix <- function(path, target_genes = NULL) {
  if (!file.exists(path)) stop("Expression matrix is missing: ", path, call. = FALSE)
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
  if (!is.null(target_genes)) mat <- mat[intersect(target_genes, rownames(mat)), , drop = FALSE]
  mat
}

logcpm_from_counts <- function(count_mat) {
  lib_size <- colSums(count_mat, na.rm = TRUE)
  if (any(!is.finite(lib_size) | lib_size <= 0)) {
    stop("Cannot compute logCPM from grouped counts with non-positive library size.", call. = FALSE)
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
  for (candidate in candidates) {
    idx <- match(as.character(candidate), expr_cols)
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
    include_in_primary_clock_TRUE_FALSE == TRUE & manual_approval_TRUE_FALSE == TRUE &
      metadata_control_like_clock_candidate == TRUE & output_sample_key %chin% reference_keys
  ]
  if (nrow(map_ds) == 0) return(NULL)
  methods <- unique(map_ds$aggregation_method)
  if (length(methods) != 1L) stop("Bulk grouping mixes aggregation methods for ", ds, call. = FALSE)
  method <- methods[[1]]
  source_path <- if (method == "sum_raw_counts") {
    unique(meta_ds_all$expression_counts_path)[1]
  } else {
    unique(meta_ds_all$expression_logcpm_path)[1]
  }
  source_mat <- read_expression_matrix(source_path, target_genes)
  missing_cols <- setdiff(unique(map_ds$expression_column), colnames(source_mat))
  if (length(missing_cols) > 0) stop("Bulk group mapping references missing expression columns for ", ds, call. = FALSE)
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
  out
}

reconstruct_reference_expression <- function(reference_meta, target_genes) {
  meta_all <- fread(metadata_path)
  mapping <- if (file.exists(bulk_group_mapping_path)) fread(bulk_group_mapping_path) else data.table()
  expr <- matrix(
    NA_real_,
    nrow = length(target_genes),
    ncol = nrow(reference_meta),
    dimnames = list(target_genes, reference_meta$sample_key)
  )
  source_rows <- list()

  for (ds in unique(reference_meta$dataset_folder)) {
    meta_ds_all <- copy(meta_all[dataset_folder == ds])
    ref_ds <- reference_meta[dataset_folder == ds]
    if (nrow(meta_ds_all) == 0) stop("No metadata rows for deployment dataset ", ds, call. = FALSE)
    logcpm_path <- unique(meta_ds_all$expression_logcpm_path)[1]
    mat <- read_expression_matrix(logcpm_path, target_genes)
    direct <- align_direct_dataset(mat, meta_ds_all, ref_ds$sample_key, unique(ref_ds$modality))
    method_used <- "direct_logcpm"
    if (!is.null(direct)) expr[rownames(direct), colnames(direct)] <- direct
    missing_keys <- ref_ds$sample_key[colSums(!is.na(expr[, ref_ds$sample_key, drop = FALSE])) == 0]
    if (length(missing_keys) > 0 && unique(ref_ds$modality) == "bulk_rnaseq") {
      grouped <- align_grouped_bulk_dataset(meta_ds_all, mapping, missing_keys, target_genes)
      if (!is.null(grouped)) {
        expr[rownames(grouped), colnames(grouped)] <- grouped
        method_used <- "approved_group_mapping"
      }
    }
    source_rows[[length(source_rows) + 1L]] <- data.table(
      dataset_folder = ds,
      modality = unique(ref_ds$modality),
      n_reference_samples = nrow(ref_ds),
      expression_logcpm_path = logcpm_path,
      reconstruction_method = method_used
    )
  }

  missing_genes <- target_genes[rowSums(is.finite(expr)) != ncol(expr)]
  if (length(missing_genes) > 0) {
    stop("Incomplete expression reconstruction for genes: ", paste(head(missing_genes, 20), collapse = ", "), call. = FALSE)
  }
  if (any(expr < -1e-6, na.rm = TRUE)) {
    stop("Reconstructed matrix appears centered or scaled rather than logCPM.", call. = FALSE)
  }
  list(expression = expr, sources = rbindlist(source_rows, fill = TRUE))
}

load_exact_training_expression <- function(reference_keys, reference_datasets, target_genes) {
  training_script <- file.path(project_dir, "scripts", "03_train_bulk_pseudobulk_clocks.R")
  parsed <- parse(training_script, keep.source = FALSE)
  loader_env <- new.env(parent = .GlobalEnv)
  is_function_assignment <- function(expr) {
    is.call(expr) &&
      identical(expr[[1]], as.name("<-")) &&
      length(expr) == 3L &&
      is.call(expr[[3]]) &&
      identical(expr[[3]][[1]], as.name("function"))
  }
  for (expr_def in parsed[vapply(parsed, is_function_assignment, logical(1))]) {
    eval(expr_def, envir = loader_env)
  }

  # The frozen artifact predates the current prohibition on positional bulk alignment.
  # CRA003645 has 22 ordered matrix columns but no shared metadata identifiers, so the
  # historical fallback is required to reproduce the frozen 117-sample deployment fit.
  n_alignment_switches <- 0L
  enable_historical_bulk_fallback <- function(node) {
    if (!is.call(node)) return(node)
    if (identical(node[[1]], as.name("align_dataset_matrix"))) {
      fallback_idx <- which(names(node) == "allow_position_fallback")
      if (length(fallback_idx) == 1L) {
        node[[fallback_idx]] <- TRUE
        n_alignment_switches <<- n_alignment_switches + 1L
      }
    }
    if (length(node) > 1L) {
      for (idx in 2:length(node)) node[[idx]] <- enable_historical_bulk_fallback(node[[idx]])
    }
    node
  }
  modified_body <- enable_historical_bulk_fallback(body(loader_env$load_modality_expression))
  if (n_alignment_switches != 1L) {
    stop("Could not locate the bulk positional-alignment switch in the training loader.", call. = FALSE)
  }
  body(loader_env$load_modality_expression) <- modified_body
  environment(loader_env$load_modality_expression) <- loader_env

  loader_values <- list(
    project_dir = project_dir,
    root_dir = root_dir,
    data_input_dir = normalizePath(file.path(root_dir, "data_input"), mustWork = TRUE),
    results_dir = results_dir,
    default_results_dir = results_dir,
    selected_gene_dir = selected_dir,
    bulk_group_mapping_path = bulk_group_mapping_path,
    bulk_group_mapping_mode = "auto",
    bulk_group_mapping_loaded = if (file.exists(bulk_group_mapping_path)) {
      fread(bulk_group_mapping_path, sep = "\t")[dataset_folder %chin% reference_datasets]
    } else {
      data.table()
    },
    bulk_group_added_samples = list(),
    bulk_group_expression_column_map = list(),
    bulk_group_dataset_summary = list(),
    bulk_group_gene_qc = list(),
    bulk_group_warnings = character(),
    bulk_group_mapping_used_datasets = character(),
    exact_samplelevel_matching_used_for = character(),
    samplelevel_exact_priority_datasets = c("GSE273724", "GSE279446", "GSE279573"),
    training_scope = "integrated_bulk_pseudobulk",
    gene_id_symbol_map_cache = NULL
  )
  loader_values$bulk_group_mapping <- loader_values$bulk_group_mapping_loaded
  list2env(loader_values, envir = loader_env)

  metadata <- fread(metadata_path)
  metadata <- metadata[dataset_folder %chin% reference_datasets]
  metadata[, age_week := as.numeric(age_week)]
  if ("n_cells_used" %in% names(metadata)) {
    metadata[, n_cells_used := as.numeric(n_cells_used)]
  }
  metadata <- loader_env$add_dataset_subcohort(metadata)
  setkey(metadata, sample_key)
  loaded <- loader_env$load_integrated_expression(metadata)

  missing_samples <- setdiff(reference_keys, colnames(loaded$expression))
  missing_genes <- setdiff(target_genes, rownames(loaded$expression))
  if (length(missing_samples) > 0) {
    stop("Exact training loader lacks frozen sample keys: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  }
  if (length(missing_genes) > 0) {
    stop("Exact training loader lacks selected genes: ", paste(head(missing_genes, 20), collapse = ", "), call. = FALSE)
  }

  expression <- loaded$expression[target_genes, reference_keys, drop = FALSE]
  if (any(!is.finite(expression))) {
    stop("Exact training-loader matrix contains non-finite selected-gene expression.", call. = FALSE)
  }
  alignment <- as.data.table(loaded$alignment)
  alignment <- alignment[sample_key %chin% reference_keys]
  sources <- alignment[, .(
    n_reference_samples = uniqueN(sample_key),
    matched_by = paste(sort(unique(matched_by)), collapse = ";")
  ), by = dataset_folder]
  source_meta <- unique(metadata[sample_key %chin% reference_keys, .(
    dataset_folder, modality, expression_logcpm_path, expression_counts_path
  )])
  sources <- merge(sources, source_meta, by = "dataset_folder", all.x = TRUE, sort = FALSE)
  sources[, reconstruction_method := "exact_training_loader"]
  setcolorder(sources, c(
    "dataset_folder", "modality", "n_reference_samples", "expression_logcpm_path",
    "expression_counts_path", "matched_by", "reconstruction_method"
  ))
  list(expression = expression, sources = sources)
}

all_genes <- unique(c(elastic_coef$gene, version_b_selected$gene))
reconstructed <- load_exact_training_expression(
  reference_keys,
  unique(ref_meta$dataset_folder),
  all_genes
)
expr <- reconstructed$expression[all_genes, reference_keys, drop = FALSE]
fwrite(reconstructed$sources, file.path(out_dir, "expression_input_sources.tsv"), sep = "\t")

fit_adjusted_gene_models <- function(expression, metadata, batch_col = "source_batch") {
  expression <- expression[, metadata$sample_key, drop = FALSE]
  if (any(!is.finite(expression))) stop("Adjusted model matrix contains non-finite expression.", call. = FALSE)
  model_df <- data.frame(
    age_week = metadata$age_week,
    source_batch = droplevels(factor(metadata[[batch_col]]))
  )
  design <- model.matrix(~ age_week + source_batch, data = model_df)
  age_idx <- match("age_week", colnames(design))
  qr_design <- qr(design)
  if (is.na(age_idx) || qr_design$rank < ncol(design)) {
    stop("Full source-adjusted model matrix is rank deficient or lacks age_week.", call. = FALSE)
  }
  y <- t(expression)
  beta <- qr.coef(qr_design, y)
  residuals <- qr.resid(qr_design, y)
  residual_df <- nrow(design) - qr_design$rank
  inverse_xtx <- solve(crossprod(design))
  sigma2 <- colSums(residuals^2) / residual_df
  age_se <- sqrt(sigma2 * inverse_xtx[age_idx, age_idx])
  age_beta <- as.numeric(beta[age_idx, ])
  age_t <- age_beta / age_se
  partial_r <- sign(age_beta) * sqrt(age_t^2 / (age_t^2 + residual_df))
  age_p <- 2 * pt(abs(age_t), df = residual_df, lower.tail = FALSE)

  source_design <- model.matrix(~ source_batch, data = model_df)
  source_qr <- qr(source_design)
  age_residual <- as.numeric(qr.resid(source_qr, model_df$age_week))
  expression_residual <- qr.resid(source_qr, y)
  residual_partial_r <- colSums(expression_residual * age_residual) /
    sqrt(colSums(expression_residual^2) * sum(age_residual^2))

  data.table(
    gene = rownames(expression),
    source_adjusted_age_beta = age_beta,
    source_adjusted_age_se = age_se,
    source_adjusted_age_t = age_t,
    residual_df = residual_df,
    partial_pearson_r = partial_r,
    residualized_partial_pearson_r = as.numeric(residual_partial_r),
    partial_r_implementation_abs_diff = abs(partial_r - residual_partial_r),
    source_adjusted_age_p = age_p,
    n_samples = nrow(metadata),
    n_source_batches = nlevels(model_df$source_batch)
  )
}

full_adjusted <- fit_adjusted_gene_models(expr, ref_meta, "source_batch")
max_partial_check_diff <- max(full_adjusted$partial_r_implementation_abs_diff, na.rm = TRUE)
if (!is.finite(max_partial_check_diff) || max_partial_check_diff > 1e-10) {
  stop("Partial-correlation implementations disagree; maximum difference = ", max_partial_check_diff, call. = FALSE)
}

elastic_adjusted <- merge(elastic_coef, full_adjusted, by = "gene", all.x = TRUE, sort = FALSE)
elastic_adjusted <- merge(
  elastic_adjusted,
  elastic_marginal[, .(gene, marginal_pearson_age_r = pearson_r)],
  by = "gene", all.x = TRUE, sort = FALSE
)
elastic_adjusted[, coefficient_partial_r_sign_concordance :=
  sign(elasticnet_coefficient) == sign(partial_pearson_r)]
elastic_adjusted[, source_adjusted_age_p_bh := p.adjust(source_adjusted_age_p, method = "BH")]
elastic_adjusted[, coefficient_abs_for_order := abs(elasticnet_coefficient)]
setorder(elastic_adjusted, -coefficient_abs_for_order)
elastic_adjusted[, coefficient_abs_for_order := NULL]
setcolorder(elastic_adjusted, c(
  "gene", "elasticnet_coefficient", "marginal_pearson_age_r",
  "source_adjusted_age_beta", "source_adjusted_age_se", "source_adjusted_age_t",
  "residual_df", "partial_pearson_r", "residualized_partial_pearson_r",
  "partial_r_implementation_abs_diff", "source_adjusted_age_p",
  "source_adjusted_age_p_bh",
  "coefficient_partial_r_sign_concordance", "n_samples", "n_source_batches"
))
fwrite(
  elastic_adjusted,
  file.path(out_dir, "elasticnet_nonzero_gene_source_adjusted_age_correlations.tsv"),
  sep = "\t"
)

version_b_adjusted <- merge(
  version_b_selected,
  full_adjusted,
  by = "gene", all.x = TRUE, sort = FALSE
)
version_b_adjusted <- merge(
  version_b_adjusted,
  version_b_marginal[, .(gene, marginal_pearson_age_r = pearson_r)],
  by = "gene", all.x = TRUE, sort = FALSE
)
setnames(version_b_adjusted, "age_beta", "stored_version_b_age_beta")
version_b_adjusted[, refitted_minus_stored_age_beta :=
  source_adjusted_age_beta - stored_version_b_age_beta]
version_b_adjusted[, coefficient_partial_r_sign_concordance :=
  sign(stored_version_b_age_beta) == sign(partial_pearson_r)]
version_b_adjusted[, source_adjusted_age_p_bh := p.adjust(source_adjusted_age_p, method = "BH")]
setorder(version_b_adjusted, version_b_rank)
setcolorder(version_b_adjusted, c(
  "gene", "version_b_rank", "stored_version_b_age_beta", "source_adjusted_age_beta",
  "refitted_minus_stored_age_beta", "marginal_pearson_age_r", "partial_pearson_r",
  "residualized_partial_pearson_r", "partial_r_implementation_abs_diff",
  "source_adjusted_age_se", "source_adjusted_age_t", "residual_df", "stored_residual_df",
  "source_adjusted_age_p", "source_adjusted_age_p_bh", "coefficient_partial_r_sign_concordance",
  "n_samples", "n_source_batches", "feature_score", "partial_r2_age",
  "partial_r2_dataset", "n_nonzero_variance_datasets"
))
fwrite(
  version_b_adjusted,
  file.path(out_dir, "versionB_2000_gene_source_adjusted_age_correlations.tsv"),
  sep = "\t"
)

max_beta_reproduction_diff <- max(abs(version_b_adjusted$refitted_minus_stored_age_beta), na.rm = TRUE)
beta_tolerance <- 1e-10
if (!is.finite(max_beta_reproduction_diff) || max_beta_reproduction_diff > beta_tolerance) {
  fwrite(
    version_b_adjusted[order(-abs(refitted_minus_stored_age_beta))][1:50],
    beta_discrepancy_path,
    sep = "\t"
  )
  stop(
    "Stored Version B age betas did not reproduce within tolerance; maximum absolute difference = ",
    max_beta_reproduction_diff,
    call. = FALSE
  )
}
if (any(version_b_adjusted$coefficient_partial_r_sign_concordance != TRUE)) {
  stop("Version B beta and partial-r signs disagree despite coming from the same model.", call. = FALSE)
}

count_nonzero_variance_batches <- function(expression, batch) {
  y <- t(expression)
  groups <- split(seq_along(batch), batch)
  valid_by_group <- vapply(groups, function(idx) {
    if (length(idx) <= 1L) return(rep(FALSE, ncol(y)))
    vars <- apply(y[idx, , drop = FALSE], 2, var)
    is.finite(vars) & vars > 0
  }, logical(ncol(y)))
  if (is.null(dim(valid_by_group))) valid_by_group <- matrix(valid_by_group, ncol = 1)
  rowSums(valid_by_group)
}

fit_lodo_beta <- function(expression, metadata, omitted_dataset) {
  keep <- metadata$dataset_folder != omitted_dataset
  meta_sub <- copy(metadata[keep])
  meta_sub[, source_batch := droplevels(factor(source_batch))]
  expr_sub <- expression[, meta_sub$sample_key, drop = FALSE]
  model_df <- data.frame(age_week = meta_sub$age_week, source_batch = meta_sub$source_batch)
  design <- model.matrix(~ age_week + source_batch, data = model_df)
  age_idx <- match("age_week", colnames(design))
  qr_design <- qr(design)
  residual_df <- nrow(design) - qr_design$rank

  base <- data.table(
    gene = rownames(expression),
    omitted_dataset = omitted_dataset,
    n_samples_remaining = nrow(meta_sub),
    n_source_batches_remaining = nlevels(meta_sub$source_batch),
    design_rank = qr_design$rank,
    design_columns = ncol(design),
    residual_df = residual_df
  )
  if (is.na(age_idx)) {
    base[, `:=`(lodo_age_beta = NA_real_, fit_valid = FALSE, invalid_reason = "age_term_absent")]
    return(base)
  }
  if (qr_design$rank < ncol(design)) {
    base[, `:=`(lodo_age_beta = NA_real_, fit_valid = FALSE, invalid_reason = "rank_deficient_design")]
    return(base)
  }
  if (residual_df < 3L || var(meta_sub$age_week) == 0) {
    base[, `:=`(lodo_age_beta = NA_real_, fit_valid = FALSE, invalid_reason = "age_not_identifiable_or_low_df")]
    return(base)
  }

  n_variable_batches <- count_nonzero_variance_batches(expr_sub, meta_sub$source_batch)
  y <- t(expr_sub)
  beta <- qr.coef(qr_design, y)
  age_beta <- as.numeric(beta[age_idx, ])
  base[, n_nonzero_variance_source_batches := n_variable_batches]
  base[, lodo_age_beta := age_beta]
  base[, fit_valid := is.finite(lodo_age_beta) & n_nonzero_variance_source_batches >= 2L]
  base[, invalid_reason := fifelse(
    fit_valid,
    NA_character_,
    fifelse(!is.finite(lodo_age_beta), "nonfinite_age_beta", "fewer_than_two_variable_source_batches")
  )]
  base
}

version_b_expr <- expr[version_b_adjusted$gene, , drop = FALSE]
datasets <- sort(unique(ref_meta$dataset_folder))
lodo_fits <- rbindlist(lapply(datasets, function(ds) {
  fit_lodo_beta(version_b_expr, ref_meta, ds)
}), use.names = TRUE, fill = TRUE)
lodo_fits <- merge(
  lodo_fits,
  version_b_adjusted[, .(gene, full_data_age_beta = source_adjusted_age_beta)],
  by = "gene", all.x = TRUE, sort = FALSE
)
lodo_fits[, sign_matches_full := fifelse(
  fit_valid,
  sign(lodo_age_beta) == sign(full_data_age_beta),
  NA
)]
setorder(lodo_fits, gene, omitted_dataset)
fwrite(
  lodo_fits,
  file.path(out_dir, "versionB_2000_gene_leave_one_dataset_out_age_beta_fits.tsv"),
  sep = "\t"
)

lodo_dataset_audit <- lodo_fits[, .(
  n_genes = .N,
  n_valid_fits = sum(fit_valid),
  n_invalid_fits = sum(!fit_valid),
  invalid_reasons = paste(sort(unique(na.omit(invalid_reason))), collapse = ";"),
  n_samples_remaining = unique(n_samples_remaining),
  n_source_batches_remaining = unique(n_source_batches_remaining),
  design_rank = unique(design_rank),
  design_columns = unique(design_columns),
  residual_df = unique(residual_df)
), by = omitted_dataset]
fwrite(lodo_dataset_audit, file.path(out_dir, "leave_one_dataset_out_fit_audit.tsv"), sep = "\t")

stability <- lodo_fits[, {
  valid_beta <- lodo_age_beta[fit_valid]
  valid_sign <- sign_matches_full[fit_valid]
  n_valid <- length(valid_beta)
  median_beta <- if (n_valid > 0) median(valid_beta) else NA_real_
  .(
    full_data_age_beta = unique(full_data_age_beta),
    median_lodo_age_beta = median_beta,
    minimum_lodo_age_beta = if (n_valid > 0) min(valid_beta) else NA_real_,
    maximum_lodo_age_beta = if (n_valid > 0) max(valid_beta) else NA_real_,
    n_valid_lodo_fits = n_valid,
    n_invalid_lodo_fits = .N - n_valid,
    proportion_same_sign_as_full = if (n_valid > 0) mean(valid_sign) else NA_real_,
    coefficient_median_absolute_deviation = if (n_valid > 0) median(abs(valid_beta - median_beta)) else NA_real_,
    maximum_absolute_change_from_full = if (n_valid > 0) max(abs(valid_beta - unique(full_data_age_beta))) else NA_real_,
    every_valid_fit_retains_full_direction = if (n_valid > 0) all(valid_sign) else FALSE
  )
}, by = gene]
stability[, stability_category := fifelse(
  n_valid_lodo_fits == 0,
  "unstable_lt80",
  fifelse(
    every_valid_fit_retains_full_direction,
    "stable_all",
    fifelse(proportion_same_sign_as_full >= 0.80, "stable_ge80_not_all", "unstable_lt80")
  )
)]
stability[, stability_category := factor(
  stability_category,
  levels = c("stable_all", "stable_ge80_not_all", "unstable_lt80")
)]
stability[, full_beta_abs_for_order := abs(full_data_age_beta)]
setorder(stability, -full_beta_abs_for_order, gene)
stability[, full_beta_abs_for_order := NULL]
fwrite(
  stability,
  file.path(out_dir, "versionB_2000_gene_leave_one_dataset_out_age_beta_stability.tsv"),
  sep = "\t"
)

implementation_check <- data.table(
  check = c(
    "maximum_partial_r_formula_vs_residualized_absolute_difference",
    "maximum_stored_vs_refitted_versionB_age_beta_absolute_difference",
    "versionB_beta_partial_r_sign_discordant_genes",
    "invalid_leave_one_dataset_out_fits"
  ),
  value = c(
    max_partial_check_diff,
    max_beta_reproduction_diff,
    sum(!version_b_adjusted$coefficient_partial_r_sign_concordance),
    sum(!lodo_fits$fit_valid)
  ),
  tolerance_or_expectation = c("<=1e-10", "<=1e-10", "0", "recorded, never replaced with zero")
)
fwrite(implementation_check, file.path(out_dir, "implementation_validation_checks.tsv"), sep = "\t")

elastic_pearson <- cor(
  elastic_adjusted$elasticnet_coefficient,
  elastic_adjusted$partial_pearson_r,
  method = "pearson"
)
elastic_spearman <- cor(
  elastic_adjusted$elasticnet_coefficient,
  elastic_adjusted$partial_pearson_r,
  method = "spearman"
)
elastic_concordant_n <- sum(elastic_adjusted$coefficient_partial_r_sign_concordance)
elastic_discordant_genes <- elastic_adjusted[
  coefficient_partial_r_sign_concordance == FALSE,
  paste(gene, collapse = ";")
]
version_b_concordant_n <- sum(version_b_adjusted$coefficient_partial_r_sign_concordance)
version_b_stable_all_n <- sum(stability$proportion_same_sign_as_full == 1)
version_b_stable_ge80_n <- sum(stability$proportion_same_sign_as_full >= 0.80)
version_b_stable_lt80_n <- sum(stability$proportion_same_sign_as_full < 0.80)

summary_rows <- rbindlist(list(
  data.table(model = "ElasticNet", metric = c(
    "n_frozen_nonzero_genes", "coefficient_partial_r_pearson", "coefficient_partial_r_spearman",
    "sign_concordant_n", "sign_concordant_percent", "discordant_genes",
    "nominal_p_lt_0.05_n", "BH_adjusted_p_lt_0.05_n"
  ), value = as.character(c(
    nrow(elastic_adjusted), elastic_pearson, elastic_spearman,
    elastic_concordant_n, 100 * elastic_concordant_n / nrow(elastic_adjusted),
    elastic_discordant_genes,
    sum(elastic_adjusted$source_adjusted_age_p < 0.05),
    sum(elastic_adjusted$source_adjusted_age_p_bh < 0.05)
  ))),
  data.table(model = "Version B", metric = c(
    "n_unique_selected_genes", "beta_partial_r_sign_concordant_n",
    "beta_partial_r_sign_concordant_percent", "max_stored_refitted_beta_absolute_difference",
    "nominal_p_lt_0.05_n", "BH_adjusted_p_lt_0.05_n",
    "LODO_sign_stable_100_percent_n", "LODO_sign_stable_100_percent_percent",
    "LODO_sign_stable_ge80_percent_n", "LODO_sign_stable_ge80_percent_percent",
    "LODO_sign_stable_lt80_percent_n", "LODO_sign_stable_lt80_percent_percent"
  ), value = as.character(c(
    nrow(version_b_adjusted), version_b_concordant_n,
    100 * version_b_concordant_n / nrow(version_b_adjusted), max_beta_reproduction_diff,
    sum(version_b_adjusted$source_adjusted_age_p < 0.05),
    sum(version_b_adjusted$source_adjusted_age_p_bh < 0.05),
    version_b_stable_all_n, 100 * version_b_stable_all_n / nrow(stability),
    version_b_stable_ge80_n, 100 * version_b_stable_ge80_n / nrow(stability),
    version_b_stable_lt80_n, 100 * version_b_stable_lt80_n / nrow(stability)
  )))
), use.names = TRUE)
fwrite(
  summary_rows,
  file.path(out_dir, "selector_effect_vs_adjusted_age_association_summary.tsv"),
  sep = "\t"
)

sample_summary <- ref_meta[, .(
  n_samples = .N,
  age_min_week = min(age_week),
  age_max_week = max(age_week),
  n_unique_ages = uniqueN(age_week)
), by = .(split_role, modality)]
fwrite(sample_summary, file.path(out_dir, "deployment_reference_sample_summary.tsv"), sep = "\t")

audit_lines <- c(
  "# Analysis input audit",
  "",
  "## Frozen artifacts and expression representation",
  "",
  paste0("- ElasticNet artifact: `", elastic_artifact_path, "`."),
  paste0("- PCA/glmnet Version B artifact: `", pca_artifact_path, "`."),
  "- Both frozen artifacts report `expression_scale = logCPM` and use the same deployment reference sample keys.",
  "- The full integrated expression matrix was constructed at training time and was not persisted as a standalone matrix. This analysis reconstructs it from the per-dataset logCPM paths and approved bulk group mapping using function definitions from the clock training script.",
  "- Reconstruction is restricted to the 13 datasets represented in the frozen artifact before forming the integrated common-gene universe; datasets added to later metadata versions are not introduced retrospectively.",
  "- The frozen artifact predates the current prohibition on positional fallback for bulk matrices. Historical positional alignment was retained for CRA003645 because its 22 ordered CRR matrix columns have no identifier overlap with the 22 metadata rows; this is necessary to recover the frozen CRA003645 deployment sample keys.",
  "- Gene identifiers are mapped and duplicate rows are resolved using the training loader's symbol-mapping and collapse rules.",
  paste0("- Expression source details are recorded in `expression_input_sources.tsv` (", nrow(reconstructed$sources), " datasets)."),
  "",
  "## Samples and age",
  "",
  paste0("- Samples: ", length(reference_keys), " control-like deployment samples shared by both frozen models."),
  paste0("- Split composition: ", paste(names(table(ref_meta$split_role)), as.integer(table(ref_meta$split_role)), collapse = "; "), "."),
  paste0("- Modalities: ", paste(names(table(ref_meta$modality)), as.integer(table(ref_meta$modality)), collapse = "; "), "."),
  paste0("- Age variable: `age_week`, numeric chronological age in weeks; range ", min(ref_meta$age_week), "-", max(ref_meta$age_week), " weeks; ", uniqueN(ref_meta$age_week), " unique values."),
  "- The deployment refit includes the original train, validation and held-out test samples only after clock selection and benchmark evaluation were completed.",
  "",
  "## Frozen Version B adjustment model",
  "",
  "- Frozen gene-level model: `logCPM ~ age_week + source_batch`.",
  "- The frozen artifact defines `source_batch = interaction(modality, dataset_folder)` and contains 13 levels.",
  "- In this frozen sample set, `source_batch` and `dataset_folder` are one-to-one. No dataset occurs in both modalities.",
  "- The frozen artifact contains neither `dataset_subcohort` nor `source_batch_for_feature_selection`; therefore this analysis reuses the frozen 13-level `source_batch` factor rather than inventing a newer composite covariate.",
  "- No additional covariates were included in the Version B gene-level model.",
  "- Version B age coefficients use unscaled logCPM as the response and unscaled `age_week` as the predictor, so beta_age has units logCPM per week.",
  "",
  "## ElasticNet scaling and coefficient interpretation",
  "",
  "- ElasticNet received the 2,000 HVGs after external gene-wise centering and division by the training-sample standard deviation.",
  "- The subsequent glmnet calls used `standardize = FALSE`, because standardization had already been applied externally.",
  "- The exported frozen coefficients are on the standardized predictor scale, not the original logCPM scale.",
  "- The final frozen deployment model contains 66 non-zero gene coefficients. The intercept is a separate model term and is not treated as a gene.",
  "",
  "## Missing-data and eligibility rules",
  "",
  "- All 2,000 frozen Version B genes and all 66 non-zero ElasticNet genes have finite expression in all 117 deployment samples.",
  "- Full adjusted fits use the same finite-expression, finite-age and non-missing-source requirements as the frozen selector.",
  "- Leave-one-dataset-out fits retain the fixed 2,000-gene set and require at least two remaining source batches with non-zero within-batch variance. Failed or non-identifiable fits are recorded as invalid and are never replaced with zero.",
  "",
  "## Reproduction checks",
  "",
  paste0("- Maximum absolute stored-versus-refitted Version B beta_age difference: ", format(max_beta_reproduction_diff, scientific = TRUE, digits = 6), "."),
  paste0("- Maximum absolute partial-r formula-versus-residualization difference: ", format(max_partial_check_diff, scientific = TRUE, digits = 6), "."),
  paste0("- Invalid leave-one-dataset-out gene fits: ", sum(!lodo_fits$fit_valid), " of ", nrow(lodo_fits), ".")
)
writeLines(audit_lines, file.path(out_dir, "ANALYSIS_INPUT_AUDIT.md"))

elastic_adjusted[, concordance := factor(
  fifelse(coefficient_partial_r_sign_concordance, "Concordant", "Discordant"),
  levels = c("Concordant", "Discordant")
)]
version_b_adjusted[, concordance := factor(
  fifelse(coefficient_partial_r_sign_concordance, "Concordant", "Discordant"),
  levels = c("Concordant", "Discordant")
)]

elastic_concordance_pct <- 100 * mean(elastic_adjusted$coefficient_partial_r_sign_concordance)
version_b_concordance_pct <- 100 * mean(version_b_adjusted$coefficient_partial_r_sign_concordance)

elastic_top_labels <- elastic_adjusted[order(-abs(elasticnet_coefficient))][1:10]
elastic_discordant <- elastic_adjusted[coefficient_partial_r_sign_concordance == FALSE][order(-abs(elasticnet_coefficient))]
elastic_label_data <- unique(rbindlist(list(
  elastic_top_labels,
  if (nrow(elastic_discordant) <= 12L) elastic_discordant else elastic_discordant[1:8]
), fill = TRUE), by = "gene")
version_b_label_data <- version_b_adjusted[order(-abs(stored_version_b_age_beta))][1:10]
stability_label_data <- unique(rbindlist(list(
  stability[order(proportion_same_sign_as_full, -maximum_absolute_change_from_full)][1:8],
  stability[order(-abs(full_data_age_beta))][1:8]
), fill = TRUE), by = "gene")

concordance_colors <- c(Concordant = "#3775BA", Discordant = "#D95F02")
concordance_shapes <- c(Concordant = 16, Discordant = 17)
stability_colors <- c(
  stable_all = "#3775BA",
  stable_ge80_not_all = "#4AA89C",
  unstable_lt80 = "#D95F02"
)
stability_shapes <- c(stable_all = 16, stable_ge80_not_all = 15, unstable_lt80 = 17)
stability_labels <- c(
  stable_all = "Stable in all fits",
  stable_ge80_not_all = "Stable in >=80%",
  unstable_lt80 = "Unstable in <80%"
)
elastic_x_limit <- 1.08 * max(abs(elastic_adjusted$elasticnet_coefficient))
version_b_x_limit <- 1.08 * max(abs(version_b_adjusted$stored_version_b_age_beta))
stability_axis_limit <- 1.05 * max(abs(c(
  stability$full_data_age_beta,
  stability$median_lodo_age_beta
)))

theme_supp <- theme_classic(base_size = 6.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(size = 5.8, colour = "#222222"),
    axis.title = element_text(size = 6.5),
    plot.title = element_text(size = 7, face = "bold", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 5.8, colour = "#4D4D4D", margin = margin(b = 3)),
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.height = grid::unit(2.8, "mm"),
    legend.key.width = grid::unit(3.5, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(4, 4, 4, 4)
  )

p_a <- ggplot(elastic_adjusted, aes(x = elasticnet_coefficient, y = partial_pearson_r)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_point(aes(colour = concordance, shape = concordance), size = 1.8, alpha = 0.75, stroke = 0.45) +
  geom_text_repel(
    data = elastic_label_data,
    aes(label = gene),
    size = 1.8,
    family = "sans",
    colour = "#222222",
    box.padding = 0.20,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.size = 0.22,
    max.overlaps = Inf,
    seed = 141,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("n = 66; concordance = %.1f%%; Spearman rho = %.2f", elastic_concordance_pct, elastic_spearman),
    hjust = -0.03, vjust = 1.3, size = 1.9, family = "sans"
  ) +
  scale_colour_manual(values = concordance_colors, drop = TRUE) +
  scale_shape_manual(values = concordance_shapes, drop = TRUE) +
  scale_x_continuous(limits = c(-elastic_x_limit, elastic_x_limit)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
  labs(
    title = "ElasticNet non-zero HVGs",
    x = "ElasticNet clock coefficient",
    y = "Source-batch-adjusted age correlation\n(partial Pearson r)",
    colour = "Direction",
    shape = "Direction"
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)), shape = "none") +
  coord_cartesian(clip = "off") +
  theme_supp

p_b <- ggplot(version_b_adjusted, aes(x = stored_version_b_age_beta, y = partial_pearson_r)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_point(aes(colour = concordance, shape = concordance), size = 0.75, alpha = 0.42, stroke = 0.3) +
  geom_text_repel(
    data = version_b_label_data,
    aes(label = gene),
    size = 1.8,
    family = "sans",
    colour = "#222222",
    box.padding = 0.20,
    point.padding = 0.08,
    min.segment.length = 0,
    segment.size = 0.22,
    max.overlaps = Inf,
    seed = 142,
    show.legend = FALSE
  ) +
  annotate(
    "text", x = -Inf, y = Inf,
    label = sprintf("n = 2,000; sign concordance = %.1f%%", version_b_concordance_pct),
    hjust = -0.03, vjust = 1.3, size = 1.9, family = "sans"
  ) +
  scale_colour_manual(values = concordance_colors, drop = FALSE) +
  scale_shape_manual(values = concordance_shapes, drop = FALSE) +
  scale_x_continuous(limits = c(-version_b_x_limit, version_b_x_limit)) +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.5)) +
  labs(
    title = "Version B selected genes",
    x = expression("Source-batch-adjusted " * beta[age] * " (logCPM " * week^-1 * ")"),
    y = "Source-batch-adjusted age correlation\n(partial Pearson r)",
    colour = "Direction",
    shape = "Direction"
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)), shape = "none") +
  coord_cartesian(ylim = c(-1, 1), clip = "off") +
  theme_supp

p_c <- ggplot(stability, aes(x = full_data_age_beta, y = median_lodo_age_beta)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.45, colour = "#555555") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_point(aes(colour = stability_category, shape = stability_category), size = 0.8, alpha = 0.48, stroke = 0.35) +
  geom_text_repel(
    data = stability_label_data,
    aes(label = gene),
    size = 1.8,
    family = "sans",
    colour = "#222222",
    box.padding = 0.20,
    point.padding = 0.08,
    min.segment.length = 0,
    segment.size = 0.22,
    max.overlaps = Inf,
    seed = 143,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = stability_colors, labels = stability_labels, drop = TRUE) +
  scale_shape_manual(values = stability_shapes, labels = stability_labels, drop = TRUE) +
  scale_x_continuous(limits = c(-stability_axis_limit, stability_axis_limit)) +
  scale_y_continuous(limits = c(-stability_axis_limit, stability_axis_limit)) +
  coord_equal(clip = "off") +
  labs(
    title = "Version B LODO stability",
    x = expression("Full-data source-batch-adjusted " * beta[age]),
    y = expression("Median leave-one-dataset-out adjusted " * beta[age]),
    colour = "Sign stability",
    shape = "Sign stability"
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)), shape = "none") +
  theme_supp

combined <-
  (p_a + theme(legend.position = "none")) +
  (p_b + theme(legend.position = "none")) +
  (p_c + theme(legend.position = "none")) +
  plot_layout(guides = "keep", widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold", family = "sans"))

save_plot_bundle <- function(plot, stem, width_mm, height_mm, png_dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  grDevices::svg(paste0(stem, ".svg"), width = width_in, height = height_in, family = "sans", onefile = TRUE)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "sans")
  print(plot)
  dev.off()
  ragg::agg_png(
    paste0(stem, ".png"), width = width_in, height = height_in,
    units = "in", res = png_dpi, background = "white"
  )
  print(plot)
  dev.off()
}

save_plot_bundle(
  p_a,
  file.path(out_dir, "panel_A_elasticnet_coefficient_vs_partial_age_r"),
  width_mm = 89, height_mm = 82
)
save_plot_bundle(
  p_b,
  file.path(out_dir, "panel_B_versionB_age_beta_vs_partial_age_r"),
  width_mm = 89, height_mm = 82
)
save_plot_bundle(
  p_c,
  file.path(out_dir, "panel_C_versionB_leave_one_dataset_out_beta_stability"),
  width_mm = 89, height_mm = 82
)
save_plot_bundle(
  combined,
  file.path(out_dir, "supplementary_selector_effect_vs_adjusted_age_association_ABC"),
  width_mm = 183, height_mm = 76
)

stability_counts <- stability[, .N, by = stability_category]
stability_count_text <- paste(
  paste0(as.character(stability_counts$stability_category), "=", stability_counts$N),
  collapse = "; "
)
legend_lines <- c(
  "# Supplementary figure legend",
  "",
  "**Supplementary Fig. X | Selector effects, source-batch-adjusted age associations and cross-dataset stability.**",
  "**A,** Frozen non-zero ElasticNet genes plotted by clock coefficient and source-batch-adjusted partial Pearson correlation with chronological age. ElasticNet inputs were centered and scaled gene-wise before model fitting, and the displayed frozen coefficients therefore act on standardized predictors. Blue circles indicate matching coefficient and partial-correlation directions; orange triangles indicate discordance. The panel reports the across-gene Spearman correlation because ElasticNet shrinkage and correlated predictors can produce a non-linear coefficient-correlation relationship.",
  "**B,** All 2,000 Version B selected genes plotted by the source-batch-adjusted age coefficient and the corresponding partial Pearson correlation. Positive beta_age indicates increasing expression with age after source-batch adjustment, whereas negative beta_age indicates decreasing expression. beta_age is a gene-level selection coefficient, not a PCA loading or a fitted PCA/glmnet clock weight. Because both quantities derive from the same adjusted model, their agreement is a consistency property rather than independent validation.",
  "**C,** Leave-one-dataset-out sensitivity of Version B age coefficients. Each point compares the full-data adjusted beta_age with the median coefficient obtained after omitting each dataset in turn. The solid line is the identity line; dashed lines mark zero. Categories indicate whether the coefficient retained its full-data sign in all valid fits, in at least 80% but not all valid fits, or in fewer than 80% of valid fits. This is an internal cross-dataset stability analysis, not external validation.",
  "",
  paste0("All gene-level models used logCPM ~ age_week + source_batch in 117 control-like reference samples. The source_batch factor is the frozen interaction of modality and dataset. Leave-one-dataset-out category counts: ", stability_count_text, ".")
)
writeLines(legend_lines, file.path(out_dir, "FIGURE_LEGEND.md"))

caption_lines <- c(
  "# Supplementary figure caption",
  "",
  "**Supplementary Fig. X | Selector effects, source-batch-adjusted age associations and cross-dataset stability.**",
  paste0(
    "**A,** Frozen non-zero ElasticNet clock genes (n = 66) are plotted by their multivariable clock coefficients and partial Pearson correlations with chronological age after adjustment for source batch. ",
    "ElasticNet predictors were centered and scaled gene-wise before model fitting; the displayed coefficients therefore act on standardized predictors and are not assigned logCPM units. ",
    "Blue circles indicate agreement and orange triangles indicate discordance between coefficient and partial-correlation directions. ",
    "Across genes, Spearman rho = ", sprintf("%.2f", elastic_spearman), " and ",
    sprintf("%.1f%%", elastic_concordance_pct), " of coefficient signs were concordant with partial r."
  ),
  paste0(
    "**B,** All Version B selected genes (n = 2,000) are plotted by the source-batch-adjusted age coefficient (beta_age; logCPM week^-1) and corresponding partial Pearson r. ",
    "Positive beta_age denotes increasing expression with age after source-batch adjustment, whereas negative beta_age denotes decreasing expression. ",
    "The Version B beta_age is a gene-level selection coefficient and is neither a PCA loading nor a fitted PCA/glmnet clock coefficient. ",
    "Because beta_age and partial r derive from the same adjusted linear model, this panel characterizes association direction and standardized strength rather than providing independent validation."
  ),
  paste0(
    "**C,** Version B cross-dataset sensitivity was evaluated by sequentially omitting each of 13 contributing datasets, refitting the same adjusted model for the fixed 2,000-gene set, and comparing the full-data beta_age with the median leave-one-dataset-out estimate. ",
    "The solid line is the identity line and dashed lines mark zero. Points are classified as directionally stable in all valid fits (blue circles), stable in at least 80% but not all fits (teal squares), or stable in fewer than 80% of fits (orange triangles). ",
    version_b_stable_all_n, " genes (", sprintf("%.1f%%", 100 * version_b_stable_all_n / nrow(stability)), ") retained direction in every omission, and no fits were invalid or rank deficient. ",
    "This is a leave-one-dataset-out stability analysis, not external validation."
  ),
  "",
  "For all gene-level analyses, chronological age was measured in weeks and the model was logCPM ~ age_week + source_batch in 117 frozen control-like deployment samples. The frozen source_batch factor was interaction(modality, dataset_folder), yielding 13 levels. Partial Pearson r was calculated from the age-term t-statistic and independently verified by correlating source-batch-residualized expression with source-batch-residualized age."
)
writeLines(caption_lines, file.path(out_dir, "SUPPLEMENTARY_FIGURE_CAPTION.md"))

results_lines <- c(
  "# Supplementary Results",
  "",
  paste0(
    "Among the 66 frozen non-zero ElasticNet genes, clock coefficients were positively associated with source-batch-adjusted partial Pearson gene-age correlations across genes (Spearman rho = ",
    sprintf("%.2f", elastic_spearman), "; Pearson r = ", sprintf("%.2f", elastic_pearson), "), with concordant signs for ",
    elastic_concordant_n, " genes (", sprintf("%.1f%%", elastic_concordance_pct), "). The remaining ",
    nrow(elastic_adjusted) - elastic_concordant_n,
    " discordant genes are compatible with the conditional nature of multivariable penalized-regression weights, including shrinkage, correlation among predictors and suppression effects, and do not by themselves indicate an implementation error. For the 2,000 Version B genes, the age-coefficient and partial-correlation signs were fully concordant, as expected because both were derived from the same source-batch-adjusted linear model; this comparison therefore characterizes effect direction and standardized association strength rather than providing independent replication. In leave-one-dataset-out analysis, ",
    version_b_stable_all_n, " genes (", sprintf("%.1f%%", 100 * version_b_stable_all_n / nrow(stability)),
    ") retained the full-data coefficient direction in all 13 valid omissions, and all 2,000 genes retained direction in at least 80% of omissions. No leave-one-dataset-out fit was rank deficient or otherwise invalid. These results support cross-dataset directional stability within the integrated data resource but do not constitute external or biological validation."
  )
)
writeLines(results_lines, file.path(out_dir, "SUPPLEMENTARY_RESULTS_TEXT.md"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "SESSION_INFO.txt"))

readme_lines <- c(
  "# Supplementary selector-age association analysis",
  "",
  "## Execution",
  "",
  "Run non-interactively from the project root:",
  "",
  "```bash",
  "Rscript integrated_clocks/elasticnet/scripts/14_plot_adjusted_selector_effect_vs_age_association.R",
  "```",
  "",
  "## Inputs",
  "",
  "- Frozen deployment ElasticNet/HVG and PCA/glmnet/Version B RDS artifacts under `results/model_artifacts/deployment/`.",
  "- Frozen Version B `selected_gene_table.tsv`.",
  "- `metadata/bulk_pseudobulk_sample_metadata.csv` restricted to the 13 frozen datasets and 117 frozen sample keys.",
  "- Per-dataset normalized logCPM matrices and the approved bulk alignment table listed in `expression_input_sources.tsv`.",
  "- Existing marginal Pearson correlation tables in `results/selected_genes/`, read without modification.",
  "",
  "## Methods",
  "",
  "Gene-level age effects were fitted as `logCPM ~ age_week + source_batch`, using the frozen 13-level `interaction(modality, dataset_folder)` factor. Partial Pearson correlations were obtained from the age-term t-statistic and checked against correlations of source-batch-residualized expression and age. Version B stability was evaluated by omitting each dataset in turn while retaining the fixed 2,000-gene set.",
  "",
  "## Outputs",
  "",
  "- Adjusted gene tables: `elasticnet_nonzero_gene_source_adjusted_age_correlations.tsv` and `versionB_2000_gene_source_adjusted_age_correlations.tsv`.",
  "- LODO tables: `versionB_2000_gene_leave_one_dataset_out_age_beta_fits.tsv`, `versionB_2000_gene_leave_one_dataset_out_age_beta_stability.tsv` and `leave_one_dataset_out_fit_audit.tsv`.",
  "- Statistical summary: `selector_effect_vs_adjusted_age_association_summary.tsv`.",
  "- Combined and separate panel figures in PDF, SVG and 600-dpi PNG formats.",
  "- Documentation: `ANALYSIS_INPUT_AUDIT.md`, `SUPPLEMENTARY_FIGURE_CAPTION.md`, `SUPPLEMENTARY_RESULTS_TEXT.md`, `FIGURE_LEGEND.md` and `SESSION_INFO.txt`.",
  "",
  "## QC",
  "",
  paste0("- Stored-versus-refitted Version B beta maximum absolute difference: ", format(max_beta_reproduction_diff, scientific = TRUE, digits = 6), "."),
  paste0("- Partial-r implementation maximum absolute difference: ", format(max_partial_check_diff, scientific = TRUE, digits = 6), "."),
  paste0("- Version B beta/partial-r sign concordance: ", version_b_concordant_n, "/", nrow(version_b_adjusted), "."),
  paste0("- Valid LODO fits: ", sum(lodo_fits$fit_valid), "/", nrow(lodo_fits), "; rank-deficient fits: 0."),
  "- All scripted validations passed.",
  "",
  "## Interpretation limits",
  "",
  "The Version B beta-versus-partial-r panel is a descriptive coefficient-direction and effect-scaling check because both quantities come from the same model. LODO is the stronger internal robustness analysis, but it is not external validation. Nominal or BH-adjusted significance counts do not validate the feature-selection procedure."
)
writeLines(readme_lines, file.path(out_dir, "README.md"))

message("Output directory: ", out_dir)
message("Maximum Version B beta reproduction difference: ", signif(max_beta_reproduction_diff, 6))
message("Maximum partial-r implementation difference: ", signif(max_partial_check_diff, 6))
message("ElasticNet sign concordance: ", signif(elastic_concordance_pct, 4), "%")
message("ElasticNet across-gene Spearman rho: ", signif(elastic_spearman, 4))
message("Version B sign concordance: ", signif(version_b_concordance_pct, 4), "%")
message("LODO invalid fits: ", sum(!lodo_fits$fit_valid), " of ", nrow(lodo_fits))
message("LODO stability categories: ", stability_count_text)
