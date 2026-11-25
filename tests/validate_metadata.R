#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
fa <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo <- if (length(args)) normalizePath(args[[1]]) else normalizePath(file.path(dirname(fa), ".."))
pi <- function(f) file.path(repo, "metadata", "pipeline_inputs", f)
mt <- function(f, ...) fread(file.path(repo, "metadata", "manuscript_tables", f), ...)

dev <- fread(pi("integrated_development_samples.csv"))
stopifnot(nrow(dev) == 136L, uniqueN(dev$sample_key) == 136L, uniqueN(dev$study) == 14L,
          sum(dev$modality == "bulk_rnaseq") == 69L, sum(dev$modality == "scrna_pseudobulk") == 67L,
          min(dev$age_week) >= 7, max(dev$age_week) <= 83)
single_age <- dev[, .(n_age = uniqueN(age_week), n = .N), by = study][n_age == 1L]
stopifnot(nrow(single_age) == 5L, sum(single_age$n) == 26L, nrow(dev) - sum(single_age$n) == 110L)

folds <- fread(pi("integrated_fold_assignments.csv"))
stopifnot(setequal(folds$sample_key, dev$sample_key),
          identical(as.integer(table(folds$cv_fold)), c(36L, 48L, 20L, 32L)),
          sum(grepl("^repeat_[0-9]{2}$", names(folds))) == 30L,
          folds[, uniqueN(cv_fold), by = study][, all(V1 == 1L)])  # whole studies held out

st2a <- mt("ST2a_model_development_profiles.tsv")
stopifnot(nrow(st2a) == 136L, setequal(st2a$sample_key, dev$sample_key))
st2b <- mt("ST2b_external_projection_profiles.tsv")
stopifnot(nrow(st2b) == 281L, !anyDuplicated(st2b$sample_key))
st4 <- mt("ST4_prespecified_contrasts.tsv")
stopifnot(nrow(st4) == 36L)

mc <- fread(pi("metacell_fold_assignments.tsv"))
stopifnot(nrow(mc) == 7049L, !anyDuplicated(mc$metacell_id), uniqueN(mc$harmonized_cell_class) == 8L,
          all(mc$min_label_assignment_confidence >= 0.70 | is.na(mc$min_label_assignment_confidence)),
          mc[, uniqueN(fold_id), by = harmonized_cell_class][, all(V1 == 5L)])
cat("metadata validation passed\n")
