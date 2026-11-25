#!/usr/bin/env Rscript

# GSE235021 alias-corrected re-projection through the frozen v2 metacell clocks (diagnostic, option 2).
#
# Problem: training metacell counts (step_02) are not symbol-harmonised across datasets. Some model
# features are rows populated by a single training dataset (2022 Morris: old symbols such as Shfm1,
# Ptrf, Tceb2 and R-sanitised mitochondrial names such as mt.Atp6) and are zero in every other training
# dataset. project_GSE235021_to_metacell_clocks.R resolved these features to the current gene through
# stable Ensembl aliases / R-name reversal and filled them with full GSE235021 expression, which the
# clocks never saw outside Morris.
#
# Correction: an alias-resolved feature (resolution_method != direct_symbol) whose training expression is
# zero in every non-Morris dataset is projected as 0, i.e. the value the clock saw for 7 of 8 training
# datasets. Sensitivity: "direct_only" sets every non-direct feature to 0.
# Outputs go to GSE235021/alias_corrected_reprojection/; release outputs are not modified.

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(glmnet)
})

run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")
gse_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection/GSE235021")
table_dir <- file.path(gse_dir, "tables")
out_dir <- file.path(gse_dir, "alias_corrected_reprojection")
dir.create(out_dir, showWarnings = FALSE)
morris_id <- "2022_eLife_Morris_cycling_murine_ovary_atlas"

models <- readRDS(file.path(run_root, "step_04_model_building/celltype_metacell_elasticnet/celltype_metacell_elasticnet_final_models.rds"))
folds <- fread(file.path(run_root, "step_03_metacell_cv_split/metacell_fold_assignments.tsv"))
train_counts <- readRDS(file.path(run_root, "step_02_preprocessing_filtering/filtered_metacell_counts.rds"))
train_lib <- Matrix::colSums(train_counts)
# Audit only the metacells the eight clocks were trained on (library size from the full matrix, as in training)
train_counts <- train_counts[, intersect(folds$metacell_id, colnames(train_counts))]
train_lib <- train_lib[colnames(train_counts)]
train_ds <- folds$dataset_id[match(colnames(train_counts), folds$metacell_id)]
stopifnot(ncol(train_counts) == nrow(folds), !anyNA(train_ds))

gse_counts <- readRDS(file.path(table_dir, "GSE235021_metacell_raw_counts.rds"))
gse_lib <- Matrix::colSums(gse_counts)
alias <- unique(fread(file.path(table_dir, "GSE235021_model_feature_alias_map.tsv")))
release <- fread(file.path(table_dir, "GSE235021_metacell_clock_predictions.tsv"))

# ---------------------------------------------------------------- feature audit
all_features <- unique(unlist(lapply(models, `[[`, "features")))
present_train <- intersect(all_features, rownames(train_counts))
cpm <- t(t(train_counts[present_train, , drop = FALSE]) / train_lib) * 1e6
by_ds <- sapply(split(seq_along(train_ds), train_ds), function(i) Matrix::rowMeans(cpm[, i, drop = FALSE]))
feature_audit <- data.table(
  model_feature = present_train,
  morris_mean_cpm = by_ds[, morris_id],
  non_morris_max_dataset_mean_cpm = apply(by_ds[, colnames(by_ds) != morris_id, drop = FALSE], 1, max),
  n_datasets_nonzero = rowSums(by_ds > 0)
)
feature_audit[, morris_only := non_morris_max_dataset_mean_cpm == 0 & morris_mean_cpm > 0]
# Same rule as the fixed projection script: any single-training-dataset row (Morris-only or other)
feature_audit[, single_dataset := n_datasets_nonzero <= 1]
feature_audit <- merge(feature_audit, alias[, .(model_feature, gse235021_gene_symbol, resolution_method)],
                       by = "model_feature", all.x = TRUE)
feature_audit[, corrected_action := fcase(
  is.na(resolution_method) | resolution_method == "unresolved_zero_filled", "zero_filled (unchanged)",
  resolution_method == "direct_symbol", "direct match (unchanged)",
  single_dataset, "alias -> set to 0 (single-training-dataset row)",
  default = "alias kept (expressed in >= 2 training datasets)")]

# ---------------------------------------------------------------- projection
project <- function(ct, mode) {
  m <- models[[ct]]
  co <- as.matrix(coef(m$fit)); b <- setNames(co[-1, 1], rownames(co)[-1])[m$features]; b0 <- co[1, 1]
  ids <- release[harmonized_cell_class == ct, metacell_id]
  a <- alias[match(m$features, alias$model_feature)]
  g <- a$gse235021_gene_symbol
  use <- !is.na(g) & g %in% rownames(gse_counts) & a$resolution_method != "unresolved_zero_filled"
  if (mode == "corrected") {
    mo <- feature_audit$single_dataset[match(m$features, feature_audit$model_feature)]
    use <- use & !(a$resolution_method %in% c("stable_ensembl_alias", "r_syntactic_name_reversal") & mo %in% TRUE)
  } else if (mode == "direct_only") {
    use <- use & a$resolution_method == "direct_symbol"
  }
  x <- matrix(0, length(m$features), length(ids), dimnames = list(m$features, ids))
  x[use, ] <- as.matrix(log1p(t(t(gse_counts[g[use], ids, drop = FALSE]) / gse_lib[ids]) * 1e6))
  data.table(metacell_id = ids, harmonized_cell_class = ct, mode = mode,
             predicted_age_weeks = as.numeric(b0 + crossprod(b, x)),
             n_features_used = sum(use), n_features = length(m$features))
}
pred <- rbindlist(lapply(names(models), function(ct) rbindlist(lapply(c("release_replica", "corrected", "direct_only"),
                                                                      function(md) project(ct, md)))))
# Replica must reproduce the release predictions exactly
chk <- merge(pred[mode == "release_replica"], release[, .(metacell_id, rel = predicted_age_weeks)], by = "metacell_id")
stopifnot(nrow(chk) == nrow(release), max(abs(chk$predicted_age_weeks - chk$rel)) < 1e-6)
meta <- unique(release[, .(metacell_id, sample_id, sample_title, diet, age_month, age_weeks)])
pred <- merge(pred, meta, by = "metacell_id")

# ---------------------------------------------------------------- library-level summaries
lib <- pred[, .(median_predicted_age_weeks = median(predicted_age_weeks), n_metacells = .N),
            by = .(mode, harmonized_cell_class, sample_id, sample_title, diet, age_month, age_weeks)]
diets <- c("polished_rice", "white_brown_rice", "black_rice")
posctrl <- lib[, {
  pd <- sapply(diets, function(d) mean(median_predicted_age_weeks[diet == d & age_month == 9]) -
                                    mean(median_predicted_age_weeks[diet == d & age_month == 6]))
  ok <- all(table(factor(diet, levels = diets), factor(age_month, levels = c(6, 9))) > 0)
  if (ok) {
    f <- lm(median_predicted_age_weeks ~ factor(age_month) + diet)
    cf <- summary(f)$coefficients["factor(age_month)9", ]
    pooled <- cf[["Estimate"]]; p1 <- pt(cf[["t value"]], f$df.residual, lower.tail = FALSE)
  } else { pooled <- NA_real_; p1 <- NA_real_ }
  .(group_mean_6m = mean(median_predicted_age_weeks[age_month == 6]),
    group_mean_9m = mean(median_predicted_age_weeks[age_month == 9]),
    diff_polished = pd[["polished_rice"]], diff_white_brown = pd[["white_brown_rice"]], diff_black = pd[["black_rice"]],
    pooled_age_effect = pooled, pooled_p_one_sided = p1,
    positive_control = fifelse(is.na(p1) | anyNA(pd), "NT", fifelse(all(pd > 0) & p1 < 0.05, "PASS", "FAIL")))
}, by = .(mode, harmonized_cell_class)]

diet_tests <- lib[diet != "polished_rice", {
  .SD
}, by = .(mode, harmonized_cell_class)][0]
diet_tests <- rbindlist(lapply(split(lib, by = c("mode", "harmonized_cell_class")), function(d) {
  rbindlist(lapply(c(6, 9), function(a) rbindlist(lapply(c("black_rice", "white_brown_rice"), function(cs) {
    x <- d[age_month == a & diet == cs, median_predicted_age_weeks]
    y <- d[age_month == a & diet == "polished_rice", median_predicted_age_weeks]
    data.table(mode = d$mode[1], harmonized_cell_class = d$harmonized_cell_class[1], age_month = a, case_diet = cs,
               n_case = length(x), n_polished = length(y),
               case_minus_polished_weeks = if (length(x) && length(y)) mean(x) - mean(y) else NA_real_,
               welch_p = if (length(x) >= 2 && length(y) >= 2) t.test(x, y)$p.value else NA_real_)
  }))))
}))

fwrite(feature_audit[order(-morris_only, resolution_method)], file.path(out_dir, "feature_alias_training_audit.tsv"), sep = "\t")
fwrite(pred, file.path(out_dir, "metacell_predictions_by_mode.tsv"), sep = "\t")
fwrite(lib, file.path(out_dir, "library_median_predictions_by_mode.tsv"), sep = "\t")
fwrite(posctrl, file.path(out_dir, "age_positive_control_by_mode.tsv"), sep = "\t")
fwrite(diet_tests, file.path(out_dir, "library_level_diet_contrasts_by_mode.tsv"), sep = "\t")

options(width = 220)
cat("\nAlias-resolved model features by training status:\n")
print(feature_audit[, .N, by = corrected_action])
cat("\nFeatures used per clock:\n")
print(dcast(unique(pred[, .(harmonized_cell_class, mode, n_features_used, n_features)]),
            harmonized_cell_class + n_features ~ mode, value.var = "n_features_used"))
cat("\nGroup mean predicted age (truth 6M = 26.1, 9M = 39.1) and age positive control:\n")
print(posctrl[order(harmonized_cell_class, mode), .(harmonized_cell_class, mode, m6 = round(group_mean_6m, 1), m9 = round(group_mean_9m, 1),
                                                    dP = round(diff_polished, 1), dW = round(diff_white_brown, 1), dB = round(diff_black, 1),
                                                    pooled = round(pooled_age_effect, 1), p = signif(pooled_p_one_sided, 2), positive_control)])
cat("\n9-month diet contrasts (library medians, Welch):\n")
print(diet_tests[age_month == 9][order(harmonized_cell_class, case_diet, mode),
                                 .(harmonized_cell_class, case_diet, mode, n_case, n_polished,
                                   diff = round(case_minus_polished_weeks, 1), p = signif(welch_p, 2))])
message("Wrote ", out_dir)
