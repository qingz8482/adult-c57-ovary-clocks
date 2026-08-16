#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)][1])
repo <- if (length(args)) normalizePath(args[[1]]) else normalizePath(file.path(dirname(file_arg), ".."))
ann <- fread(file.path(repo, "metadata", "sample_annotations.tsv"))
stopifnot(
  nrow(ann) == 411L,
  uniqueN(ann$sample_key) == 411L,
  sum(ann$modality == "bulk_rnaseq") == 294L,
  sum(ann$modality == "scrna_pseudobulk") == 117L
)
splits <- fread(file.path(repo, "metadata", "split_assignments.tsv"))
stopifnot(all(c("sample_key", "split_role") %in% names(splits)), uniqueN(splits$sample_key) == nrow(splits))
folds <- fread(file.path(repo, "metadata", "metacell_cv_assignments.tsv"))
stopifnot(all(c("metacell_id", "fold_id") %in% names(folds)), nrow(folds) > 5000L)
cat("metadata validation passed\n")
