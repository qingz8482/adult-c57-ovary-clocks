#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}

project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_input_dir <- normalizePath(file.path(project_dir, "..", "data_input"), mustWork = TRUE)
results_dir <- Sys.getenv("RESULTS_ROOT", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = TRUE)
out_dir <- file.path(results_dir, "external_projection")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_external_projection_datasets <- c(
  "CRA003645", "CRA017492", "GSE194196", "GSE252791", "GSE273724",
  "GSE279446", "GSE279573", "GSE329024", "GSE228975", "GSE240547",
  "PRJNA1076509"
)

model_paths <- file.path(
  results_dir,
  "model_artifacts",
  "deployment",
  c(
    "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds",
    "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
  )
)

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

load_gene_id_symbol_map <- function() {
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
    maps[[length(maps) + 1L]] <- dt[, .(
      gene_id_base = sub("[.][0-9]+$", "", na_chr(V1)),
      gene_symbol = na_chr(V2)
    )]
  }

  if (length(maps) == 0) return(data.table(gene_id_base = character(), gene_symbol = character()))
  mapping <- rbindlist(maps, use.names = TRUE, fill = TRUE)
  mapping <- mapping[
    grepl("^ENSMUSG", gene_id_base) &
      !is.na(gene_symbol) &
      gene_symbol != "" &
      !grepl("^ENSMUSG", gene_symbol)
  ]
  mapping <- mapping[, .N, by = .(gene_id_base, gene_symbol)]
  setorder(mapping, gene_id_base, -N, gene_symbol)
  mapping[, .SD[1], by = gene_id_base][, .(gene_id_base, gene_symbol)]
}

if (!all(file.exists(model_paths))) {
  stop("Missing deployment model artifacts: ", paste(model_paths[!file.exists(model_paths)], collapse = ", "), call. = FALSE)
}
models <- lapply(model_paths, readRDS)
model_genes <- lapply(models, function(x) x$selected_genes)
names(model_genes) <- vapply(models, function(x) paste(x$model_family, x$selector, sep = "__"), character(1))

mapping <- load_gene_id_symbol_map()
symbol_lookup <- mapping$gene_symbol
names(symbol_lookup) <- mapping$gene_id_base

matrix_paths <- Sys.glob(file.path(data_input_dir, "bulkseq_dataset", "*", "input", "*_gene_logcpm_filtered.csv"))
if (length(matrix_paths) == 0) stop("No bulk gene_logcpm_filtered matrices found.", call. = FALSE)

audit <- list()
for (matrix_path in matrix_paths) {
  dataset <- basename(dirname(dirname(matrix_path)))
  dt <- fread(matrix_path, select = 1)
  row_ids <- na_chr(dt[[1]])
  row_id_base <- sub("[.][0-9]+$", "", row_ids)
  mapped_symbols <- unname(symbol_lookup[row_id_base])
  mapped_gene_ids <- row_ids
  use_mapped <- grepl("^ENSMUSG", row_ids) & !is.na(mapped_symbols) & mapped_symbols != ""
  mapped_gene_ids[use_mapped] <- mapped_symbols[use_mapped]
  row_id_type <- if (all(grepl("^ENSMUSG", row_ids))) {
    if (any(grepl("^ENSMUSG[0-9]+[.][0-9]+$", row_ids))) "ensembl_versioned" else "ensembl_unversioned"
  } else if (any(grepl("^ENSMUSG", row_ids))) {
    "mixed"
  } else {
    "symbol_like"
  }

  for (model_name in names(model_genes)) {
    genes <- model_genes[[model_name]]
    raw_overlap <- length(intersect(row_ids, genes))
    mapped_overlap <- length(intersect(mapped_gene_ids, genes))
    audit[[length(audit) + 1L]] <- data.table(
      dataset_folder = dataset,
      matrix_path = matrix_path,
      row_id_type = row_id_type,
      model_family = sub("__.*$", "", model_name),
      selector = sub("^[^_]+__", "", model_name),
      n_matrix_rows = length(row_ids),
      n_ensembl_rows = sum(grepl("^ENSMUSG", row_ids)),
      n_ensembl_rows_mapped = sum(use_mapped),
      n_model_genes = length(genes),
      raw_model_gene_overlap = raw_overlap,
      mapped_model_gene_overlap = mapped_overlap,
      raw_gene_coverage_fraction = raw_overlap / length(genes),
      mapped_gene_coverage_fraction = mapped_overlap / length(genes),
      required_external_projection_dataset = dataset %in% required_external_projection_datasets
    )
  }
}

audit <- rbindlist(audit, use.names = TRUE)
audit[, status := fifelse(
  required_external_projection_dataset & mapped_gene_coverage_fraction < 0.80,
  "FAIL_low_mapped_clock_gene_coverage",
  "PASS"
)]
setorder(audit, dataset_folder, model_family, selector)
out_path <- file.path(out_dir, "input_gene_identifier_audit.tsv")
fwrite(audit, out_path, sep = "\t")

failures <- audit[status != "PASS"]
cat("Input gene identifier audit:", out_path, "\n")
cat("Rows:", nrow(audit), "\n")
cat("Failures:", nrow(failures), "\n")
if (nrow(failures) > 0) {
  print(failures)
  stop("Input gene identifier validation failed.", call. = FALSE)
}
cat("INPUT_GENE_IDENTIFIER_VALIDATION_PASS\n")
