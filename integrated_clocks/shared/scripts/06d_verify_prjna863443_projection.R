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
results_dir <- Sys.getenv("RESULTS_ROOT", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = TRUE)
if (grepl("results_groupmapped_v2", results_dir, fixed = TRUE) || grepl("_v2", basename(results_dir), fixed = TRUE)) {
  stop("Refusing to verify PRJNA863443 against a v2 results root: ", results_dir, call. = FALSE)
}

projection_dir <- file.path(results_dir, "external_projection")
dataset_dir <- file.path(projection_dir, "PRJNA863443")
out_dir <- file.path(dataset_dir, "verification")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_required <- function(path) {
  if (!file.exists(path)) stop("Required input missing: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

ensure_cols <- function(dt, cols) {
  for (cc in cols) if (!cc %in% names(dt)) dt[, (cc) := NA]
  dt
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

pred <- read_required(file.path(dataset_dir, "prjna863443_projection_predictions.tsv"))
rel <- read_required(file.path(dataset_dir, "prjna863443_relative_to_controls.tsv"))
meta <- read_required(file.path(dataset_dir, "prjna863443_projection_metadata.tsv"))
coverage <- read_required(file.path(dataset_dir, "prjna863443_gene_coverage.tsv"))
main_pred <- read_required(file.path(projection_dir, "external_projection_predictions_all.tsv"))
main_rel <- read_required(file.path(projection_dir, "external_projection_relative_to_controls_all.tsv"))

cols <- c(
  "sample_key", "dataset_folder", "model_family", "selector", "projection_class",
  "projection_role", "projection_group", "comparison_id", "age_group",
  "condition", "treatment", "genotype", "matched_control_method",
  "n_matched_controls", "relative_delta_age_week", "primary_external_projection_include",
  "primary_effect_interpretation_include", "matched_control_allowed",
  "matched_reference_control_flag", "gene_coverage_fraction"
)
pred <- ensure_cols(pred, cols)
rel <- ensure_cols(rel, cols)
main_pred <- ensure_cols(main_pred, cols)
main_rel <- ensure_cols(main_rel, cols)
coverage <- ensure_cols(coverage, c("sample_key", "model_family", "gene_coverage_fraction", "projection_gene_coverage_flag"))

for (dt in list(pred, rel, main_pred, main_rel)) {
  for (cc in c("primary_external_projection_include", "primary_effect_interpretation_include", "matched_control_allowed", "matched_reference_control_flag")) {
    if (cc %in% names(dt)) dt[, (cc) := as_bool(get(cc))]
  }
  if ("n_matched_controls" %in% names(dt)) dt[, n_matched_controls := as.integer(n_matched_controls)]
  if ("relative_delta_age_week" %in% names(dt)) dt[, relative_delta_age_week := safe_num(relative_delta_age_week)]
  if ("gene_coverage_fraction" %in% names(dt)) dt[, gene_coverage_fraction := safe_num(gene_coverage_fraction)]
}
coverage[, gene_coverage_fraction := safe_num(gene_coverage_fraction)]

checks <- list()
add_check <- function(scope, check, passed, detail = "") {
  checks[[length(checks) + 1L]] <<- data.table(
    scope = scope,
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail)
  )
}

add_check("files", "dataset_specific_predictions_exist", nrow(pred) > 0, paste0("rows=", nrow(pred)))
add_check("files", "main_predictions_include_PRJNA863443", nrow(main_pred[dataset_folder == "PRJNA863443"]) > 0, paste0("rows=", nrow(main_pred[dataset_folder == "PRJNA863443"])))
add_check("files", "main_relative_include_PRJNA863443", nrow(main_rel[dataset_folder == "PRJNA863443"]) > 0, paste0("rows=", nrow(main_rel[dataset_folder == "PRJNA863443"])))

add_check("counts", "exactly_33_samples_metadata", nrow(meta) == 33 && uniqueN(meta$sample_key) == 33, paste0("rows=", nrow(meta), ";unique=", uniqueN(meta$sample_key)))
add_check("counts", "exactly_33_samples_predictions", uniqueN(pred$sample_key) == 33, paste0("n=", uniqueN(pred$sample_key)))
add_check("counts", "exactly_66_prediction_rows", nrow(pred) == 66, paste0("rows=", nrow(pred)))
add_check("counts", "both_model_families_present", identical(sort(unique(pred$model_family)), c("elasticnet", "pca_glmnet")), paste(sort(unique(pred$model_family)), collapse = ";"))

natural <- pred[comparison_id == "PRJNA863443_natural_aging_reference"]
natural_rel <- rel[comparison_id == "PRJNA863443_natural_aging_reference"]
add_check("natural_aging", "seven_samples_projected", uniqueN(natural$sample_key) == 7, paste0("n=", uniqueN(natural$sample_key)))
add_check("natural_aging", "raw_predictions_present", all(is.finite(safe_num(natural$predicted_age_week)) & is.finite(safe_num(natural$delta_age_week))), "predicted_age_week/delta_age_week finite")
add_check("natural_aging", "no_primary_effect_interpretation", all(natural$primary_effect_interpretation_include == FALSE), paste(unique(natural$primary_effect_interpretation_include), collapse = ";"))
add_check("natural_aging", "no_matched_control_relative_delta", nrow(natural_rel) == 14 && all(natural_rel$matched_control_method == "no_matched_control") && all(!is.finite(natural_rel$relative_delta_age_week)), paste0("rows=", nrow(natural_rel), ";finite=", sum(is.finite(natural_rel$relative_delta_age_week))))

vcd_expected <- c("3m_30d", "3m_90d", "10m_30d", "10m_90d")
for (grp in vcd_expected) {
  method <- paste0("PRJNA863443_VCD_vs_CTL_", grp)
  pred_grp <- pred[comparison_id == method]
  rel_grp <- rel[comparison_id == method]
  ctl <- unique(pred_grp[grepl(paste0("^CTL_", grp, "$"), projection_group)]$sample_key)
  vcd <- unique(pred_grp[grepl(paste0("^VCD_", grp, "$"), projection_group)]$sample_key)
  add_check("VCD", paste0(grp, "_two_CTL_controls"), length(ctl) == 2, paste0("n=", length(ctl)))
  add_check("VCD", paste0(grp, "_two_VCD_samples"), length(vcd) == 2, paste0("n=", length(vcd)))
  trt_rows <- rel_grp[grepl(paste0("^VCD_", grp, "$"), projection_group)]
  add_check("VCD", paste0(grp, "_VCD_rows_method"), nrow(trt_rows) == 4 && all(trt_rows$matched_control_method == method), paste0("rows=", nrow(trt_rows), ";methods=", paste(unique(trt_rows$matched_control_method), collapse = ";")))
  add_check("VCD", paste0(grp, "_VCD_rows_two_controls"), nrow(trt_rows) == 4 && all(trt_rows$n_matched_controls == 2), paste0("n=", paste(unique(trt_rows$n_matched_controls), collapse = ";")))
  ctl_rows <- rel_grp[grepl(paste0("^CTL_", grp, "$"), projection_group)]
  add_check("VCD", paste0(grp, "_CTL_rows_centered"), nrow(ctl_rows) == 4 && all(ctl_rows$n_matched_controls == 2) && abs(mean(ctl_rows$relative_delta_age_week, na.rm = TRUE)) < 1e-8, paste0("mean=", mean(ctl_rows$relative_delta_age_week, na.rm = TRUE)))
}

fox_expected <- data.table(age_group = c("young", "mid_age"), n_controls = c(2L, 3L), n_het = c(2L, 3L))
for (i in seq_len(nrow(fox_expected))) {
  ag <- fox_expected$age_group[i]
  method <- paste0("PRJNA863443_Foxl2_het_vs_WT_", ag)
  pred_grp <- pred[comparison_id == method]
  rel_grp <- rel[comparison_id == method]
  wt <- unique(pred_grp[projection_group == paste0("Foxl2+/+_", ag)]$sample_key)
  het <- unique(pred_grp[projection_group == paste0("Foxl2+/-_", ag)]$sample_key)
  add_check("Foxl2", paste0(ag, "_WT_controls"), length(wt) == fox_expected$n_controls[i], paste0("n=", length(wt)))
  add_check("Foxl2", paste0(ag, "_het_samples"), length(het) == fox_expected$n_het[i], paste0("n=", length(het)))
  het_rows <- rel_grp[projection_group == paste0("Foxl2+/-_", ag)]
  add_check("Foxl2", paste0(ag, "_het_rows_method"), nrow(het_rows) == fox_expected$n_het[i] * 2 && all(het_rows$matched_control_method == method), paste0("rows=", nrow(het_rows), ";methods=", paste(unique(het_rows$matched_control_method), collapse = ";")))
  add_check("Foxl2", paste0(ag, "_het_rows_expected_controls"), nrow(het_rows) == fox_expected$n_het[i] * 2 && all(het_rows$n_matched_controls == fox_expected$n_controls[i]), paste0("n=", paste(unique(het_rows$n_matched_controls), collapse = ";")))
  wt_rows <- rel_grp[projection_group == paste0("Foxl2+/+_", ag)]
  add_check("Foxl2", paste0(ag, "_WT_rows_centered"), nrow(wt_rows) == fox_expected$n_controls[i] * 2 && abs(mean(wt_rows$relative_delta_age_week, na.rm = TRUE)) < 1e-8, paste0("mean=", mean(wt_rows$relative_delta_age_week, na.rm = TRUE)))
}

invalid <- rel[
  (matched_control_method == "no_matched_control" & is.finite(relative_delta_age_week)) |
    (n_matched_controls == 0 & is.finite(relative_delta_age_week))
]
add_check("relative", "no_invalid_relative_rows", nrow(invalid) == 0, paste0("rows=", nrow(invalid)))

global_invalid <- main_rel[
  (matched_control_method == "no_matched_control" & is.finite(relative_delta_age_week)) |
    (n_matched_controls == 0 & is.finite(relative_delta_age_week)) |
    (matched_control_allowed == FALSE & is.finite(relative_delta_age_week)) |
    (primary_external_projection_include == FALSE & is.finite(relative_delta_age_week))
]
add_check("relative", "global_no_invalid_relative_rows", nrow(global_invalid) == 0, paste0("rows=", nrow(global_invalid)))

cov_summary <- coverage[, .(
  min_gene_coverage = min(gene_coverage_fraction, na.rm = TRUE),
  median_gene_coverage = median(gene_coverage_fraction, na.rm = TRUE),
  n_low_coverage = sum(gene_coverage_fraction < 0.80, na.rm = TRUE)
), by = model_family]
add_check("coverage", "gene_coverage_ge_0_80_all_rows", all(coverage$gene_coverage_fraction >= 0.80), paste(paste0(cov_summary$model_family, " min=", round(cov_summary$min_gene_coverage, 4)), collapse = ";"))

results <- rbindlist(checks, use.names = TRUE, fill = TRUE)
failures <- results[passed != TRUE]
fwrite(results, file.path(out_dir, "prjna863443_projection_verification.tsv"), sep = "\t")
fwrite(failures, file.path(out_dir, "prjna863443_projection_verification_failures.tsv"), sep = "\t")

status <- if (nrow(failures) == 0) "PRJNA863443_PROJECTION_VERIFICATION_PASS" else "PRJNA863443_PROJECTION_VERIFICATION_FAIL"
summary_lines <- c(
  "# PRJNA863443 projection verification summary",
  "",
  paste0("- PRJNA863443 sample count: ", uniqueN(pred$sample_key), " (expected 33)."),
  paste0("- PRJNA863443 prediction rows: ", nrow(pred), " (expected 66)."),
  paste0("- Gene coverage by model: ", paste(paste0(cov_summary$model_family, " min=", round(cov_summary$min_gene_coverage, 4), " median=", round(cov_summary$median_gene_coverage, 4)), collapse = "; "), "."),
  paste0("- VCD matched-control verification status: ", if (all(results[scope == "VCD"]$passed)) "PASS" else "FAIL", "."),
  paste0("- Foxl2 matched-control verification status: ", if (all(results[scope == "Foxl2"]$passed)) "PASS" else "FAIL", "."),
  paste0("- Natural ageing projection status: ", if (all(results[scope == "natural_aging"]$passed)) "PASS" else "FAIL", "."),
  paste0("- Global invalid relative-control row count: ", nrow(global_invalid), " (expected 0)."),
  paste0("- Overall status: ", status, ".")
)
writeLines(summary_lines, file.path(out_dir, "prjna863443_projection_verification_summary.md"))

cat("PRJNA863443 sample count:", uniqueN(pred$sample_key), "\n")
cat("PRJNA863443 prediction rows:", nrow(pred), "\n")
cat("PRJNA863443 gene coverage by model:\n")
print(cov_summary)
cat("VCD matched-control verification status:", if (all(results[scope == "VCD"]$passed)) "PASS" else "FAIL", "\n")
cat("Foxl2 matched-control verification status:", if (all(results[scope == "Foxl2"]$passed)) "PASS" else "FAIL", "\n")
cat("Natural ageing projection status:", if (all(results[scope == "natural_aging"]$passed)) "PASS" else "FAIL", "\n")
cat("Global invalid relative-control row count:", nrow(global_invalid), "\n")
cat(status, "\n")

if (nrow(failures) > 0) quit(status = 1)
