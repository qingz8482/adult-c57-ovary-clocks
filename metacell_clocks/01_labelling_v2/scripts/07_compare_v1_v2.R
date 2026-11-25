# Step D3: compare the 20260707 (v1 labels) clocks with v2 main and sensitivity runs.
suppressPackageStartupMessages(library(data.table))
V2 <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
V1 <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/7 scRNAseq_clock_results/clocks/metacell_level_scRNA/scRNA_metacell_elasticnet_metacell_cv_20260707")
out_dir <- file.path(V2, "06_clock_run", "comparison")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- c(v1_20260707 = V1,
          v2_main = file.path(V2, "06_clock_run", "cv"),
          v2_confirmed_only = file.path(V2, "06_clock_run", "cv_confirmed_only"),
          v2_granulosa_no_luteal = file.path(V2, "06_clock_run", "cv_no_luteal"),
          v2_leave_one_dataset_out = file.path(V2, "06_clock_run", "cv_lodo"))

ev <- rbindlist(lapply(names(runs), function(r) {
  f <- file.path(runs[[r]], "step_05_internal_cv_validation", "celltype_clock_cv_evaluation_summary.tsv")
  if (!file.exists(f)) return(NULL)
  d <- fread(f); d[, run := r]; d
}), fill = TRUE)
fwrite(ev, file.path(out_dir, "cv_evaluation_all_runs.tsv"), sep = "\t")
wide <- dcast(ev, harmonized_cell_class ~ run, value.var = c("n_metacells", "mae_weeks", "pearson_r"))
fwrite(wide, file.path(out_dir, "cv_evaluation_wide.tsv"), sep = "\t")

# Final-model features: v1 vs v2 main (union of nonzero genes in the final all-metacell model)
feat <- function(root) {
  f <- file.path(root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_coefficients_by_fold.tsv")
  if (!file.exists(f)) return(NULL)
  d <- fread(f)[nonzero_coefficient == TRUE & split_id == "metacell_internal_cv"]
  # stable features: nonzero in >= 3 outer folds; mean coefficient across folds
  d[, .(n_folds = uniqueN(fold_id), coefficient = mean(coefficient)), by = .(harmonized_cell_class, gene)][n_folds >= 3]
}
f1 <- feat(V1); f2 <- feat(runs[["v2_main"]])
if (!is.null(f1) && !is.null(f2)) {
  cls <- intersect(unique(f1$harmonized_cell_class), unique(f2$harmonized_cell_class))
  ov <- rbindlist(lapply(cls, function(k) {
    a <- f1[harmonized_cell_class == k, gene]; b <- f2[harmonized_cell_class == k, gene]
    m <- merge(f1[harmonized_cell_class == k], f2[harmonized_cell_class == k], by = "gene")
    data.table(harmonized_cell_class = k, n_v1 = length(a), n_v2 = length(b), n_shared = length(intersect(a, b)),
               jaccard = length(intersect(a, b)) / length(union(a, b)),
               sign_agreement_shared = if (nrow(m)) mean(sign(m$coefficient.x) == sign(m$coefficient.y)) else NA_real_)
  }))
  fwrite(ov, file.path(out_dir, "stable_feature_overlap_v1_v2.tsv"), sep = "\t")
  print(ov)
}

# Metacell composition (v2 main): luteal fraction and label sources per class
mm <- file.path(V2, "06_clock_run", "build", "step_02_preprocessing_filtering", "filtered_metacell_metadata.tsv")
if (file.exists(mm)) {
  m <- fread(mm)[model_role == "clock_candidate"]
  comp <- m[, .(n_metacells = .N, n_cells = sum(n_cells), n_datasets = uniqueN(dataset_id),
                mean_frac_luteal = mean(frac_luteal), mean_frac_confirmed = mean(frac_confirmed)), by = harmonized_cell_class]
  fwrite(comp, file.path(out_dir, "v2_metacell_composition_by_class.tsv"), sep = "\t")
  lut <- m[harmonized_cell_class == "Granulosa", .(n_metacells = .N, mean_frac_luteal = mean(frac_luteal)), by = .(dataset_id, reproductive_age_bin)]
  fwrite(lut[order(dataset_id, reproductive_age_bin)], file.path(out_dir, "v2_granulosa_luteal_fraction_by_dataset_agebin.tsv"), sep = "\t")
  print(comp)
}
options(width = 220)
print(wide)
