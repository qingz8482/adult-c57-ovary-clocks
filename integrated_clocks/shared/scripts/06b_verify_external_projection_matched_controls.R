#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

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
  stop("Refusing to verify v1 external projection against a v2 results root: ", results_dir, call. = FALSE)
}

projection_dir <- file.path(results_dir, "external_projection")
out_dir <- file.path(projection_dir, "matched_control_verification")
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

pred_all <- read_required(file.path(projection_dir, "external_projection_predictions_all.tsv"))
pred_primary <- read_required(file.path(projection_dir, "external_projection_predictions.tsv"))
rel_all <- read_required(file.path(projection_dir, "external_projection_relative_to_controls_all.tsv"))
rel_primary <- read_required(file.path(projection_dir, "external_projection_relative_to_controls.tsv"))
availability <- read_required(file.path(projection_dir, "external_projection_dataset_control_availability.tsv"))
unreliable <- read_required(file.path(projection_dir, "external_projection_unreliable_relative_controls_audit.tsv"))
reference_controls <- read_required(file.path(projection_dir, "external_projection_reference_controls.tsv"))
simplified_path <- file.path(
  projection_dir, "simplified_classes", "external_projection_predictions_all_simplified.tsv"
)
simplified_all <- if (file.exists(simplified_path)) read_required(simplified_path) else NULL
analysis_domain_count <- if (!is.null(simplified_all) &&
                             all(c("sample_key", "plot_include_primary") %in% names(simplified_all))) {
  uniqueN(simplified_all[as_bool(plot_include_primary), sample_key])
} else {
  NA_integer_
}
developmental_ood_count <- if (is.finite(analysis_domain_count)) {
  uniqueN(pred_all$sample_key) - analysis_domain_count
} else {
  NA_integer_
}

common_cols <- c(
  "dataset_folder", "sample_key", "model_family", "selector", "projection_role",
  "projection_class", "primary_external_projection_include", "primary_effect_interpretation_include",
  "matched_control_allowed", "matched_reference_control_flag",
  "matched_reference_control_source", "used_in_deployment_training_flag",
  "reference_control_group_id", "gse279446_cohort", "treatment", "genotype", "sample_title",
  "projection_group", "comparison_id", "n_observed_model_genes", "gene_coverage_fraction",
  "projection_gene_coverage_flag"
)
pred_all <- ensure_cols(pred_all, common_cols)
pred_primary <- ensure_cols(pred_primary, common_cols)
rel_all <- ensure_cols(rel_all, c(
  common_cols, "matched_control_method", "matched_control_source",
  "n_matched_controls", "matched_control_sample_keys", "relative_delta_age_week",
  "control_expected_delta_age_week", "control_expected_predicted_age_week"
))
rel_primary <- ensure_cols(rel_primary, c(common_cols, "matched_control_method", "n_matched_controls", "relative_delta_age_week"))
reference_controls <- ensure_cols(reference_controls, c(
  "dataset_folder", "sample_key", "model_family", "projection_role",
  "matched_reference_control_source", "used_in_deployment_training_flag",
  "reference_control_group_id"
))

for (dt in list(pred_all, pred_primary, rel_all, rel_primary, reference_controls)) {
  if ("primary_external_projection_include" %in% names(dt)) dt[, primary_external_projection_include := as_bool(primary_external_projection_include)]
  if ("primary_effect_interpretation_include" %in% names(dt)) dt[, primary_effect_interpretation_include := as_bool(primary_effect_interpretation_include)]
  if ("matched_control_allowed" %in% names(dt)) dt[, matched_control_allowed := as_bool(matched_control_allowed)]
  if ("matched_reference_control_flag" %in% names(dt)) dt[, matched_reference_control_flag := as_bool(matched_reference_control_flag)]
  if ("used_in_deployment_training_flag" %in% names(dt)) dt[, used_in_deployment_training_flag := as_bool(used_in_deployment_training_flag)]
}
rel_all[, `:=`(
  n_matched_controls = as.integer(n_matched_controls),
  relative_delta_age_week = suppressWarnings(as.numeric(relative_delta_age_week)),
  control_expected_delta_age_week = suppressWarnings(as.numeric(control_expected_delta_age_week)),
  control_expected_predicted_age_week = suppressWarnings(as.numeric(control_expected_predicted_age_week))
)]
pred_all[, `:=`(
  n_observed_model_genes = suppressWarnings(as.integer(n_observed_model_genes)),
  gene_coverage_fraction = suppressWarnings(as.numeric(gene_coverage_fraction))
)]
pred_primary[, `:=`(
  n_observed_model_genes = suppressWarnings(as.integer(n_observed_model_genes)),
  gene_coverage_fraction = suppressWarnings(as.numeric(gene_coverage_fraction))
)]

checks <- list()
add_check <- function(scope, check, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.table(
    scope = scope,
    check = check,
    passed = isTRUE(passed),
    detail = as.character(detail)
  )
}

required_datasets <- c(
  "CRA003645", "CRA017492", "GSE194196", "GSE252791", "GSE273724",
  "GSE279446", "GSE279573", "GSE329024", "GSE228975", "GSE240547",
  "PRJNA1076509", "PRJNA863443"
)

add_check("global", "all_required_datasets_in_availability", all(required_datasets %in% availability$dataset_folder), paste(setdiff(required_datasets, availability$dataset_folder), collapse = ";"))
add_check("global", "both_model_families_present", all(c("elasticnet", "pca_glmnet") %in% unique(pred_all$model_family)), paste(sort(unique(pred_all$model_family)), collapse = ";"))
add_check("global", "unreliable_relative_control_audit_empty", nrow(unreliable) == 0, paste0("rows=", nrow(unreliable)))

low_coverage_primary <- pred_primary[
  primary_external_projection_include == TRUE &
    (is.na(n_observed_model_genes) | n_observed_model_genes <= 0 |
      is.na(gene_coverage_fraction) | gene_coverage_fraction < 0.80 |
      grepl("low", as.character(projection_gene_coverage_flag), ignore.case = TRUE))
]
add_check(
  "global",
  "primary_projection_gene_coverage_ok",
  nrow(low_coverage_primary) == 0,
  paste0(
    "rows=", nrow(low_coverage_primary),
    if (nrow(low_coverage_primary) > 0) paste0(
      ";datasets=", paste(sort(unique(low_coverage_primary$dataset_folder)), collapse = ";")
    ) else ""
  )
)

invalid_relative <- rel_all[
  (matched_control_method == "no_matched_control" & is.finite(relative_delta_age_week)) |
    (n_matched_controls == 0 & is.finite(relative_delta_age_week)) |
    (matched_control_allowed == FALSE & is.finite(relative_delta_age_week)) |
    (primary_external_projection_include == FALSE & is.finite(relative_delta_age_week))
]
add_check("global", "no_invalid_relative_rows", nrow(invalid_relative) == 0, paste0("rows=", nrow(invalid_relative)))

check_named_controls <- function(ds, method, n_controls) {
  pred_ds <- pred_all[dataset_folder == ds]
  rel_ds <- rel_all[dataset_folder == ds]
  ref_ds <- reference_controls[dataset_folder == ds & reference_control_group_id == method]
  treatment_rel <- rel_ds[projection_role == "external_projected_treatment"]
  add_check(ds, "reference_controls_detected", uniqueN(ref_ds$sample_key) == n_controls, paste0("n=", uniqueN(ref_ds$sample_key)))
  add_check(ds, "reference_controls_in_predictions", uniqueN(pred_ds[matched_reference_control_flag == TRUE & reference_control_group_id == method]$sample_key) == n_controls, paste0("n=", uniqueN(pred_ds[matched_reference_control_flag == TRUE & reference_control_group_id == method]$sample_key)))
  add_check(ds, "treatment_rows_matched_to_named_controls", nrow(treatment_rel) > 0 && all(treatment_rel$matched_control_method == method), paste(unique(treatment_rel$matched_control_method), collapse = ";"))
  add_check(ds, paste0("treatment_rows_have_", n_controls, "_matched_controls"), nrow(treatment_rel) > 0 && all(treatment_rel$n_matched_controls == n_controls), paste(unique(treatment_rel$n_matched_controls), collapse = ";"))
}

check_named_controls("GSE194196", "GSE194196_WT_1_4_7", 3)
check_named_controls("GSE252791", "GSE252791_WT_1_4_7", 3)

gse329024_expected <- data.table(
  design_token = c("2M-8W", "2M-12W", "6M-8W", "6M-12W"),
  reference_control_group_id = c(
    "GSE329024_2M_8W_control_meta",
    "GSE329024_2M_12W_control_meta",
    "GSE329024_6M_8W_control_meta",
    "GSE329024_6M_12W_control_meta"
  ),
  expected_samples = c(8L, 8L, 8L, 7L),
  expected_controls = c(4L, 4L, 4L, 3L)
)
gse329024_pred <- pred_all[dataset_folder == "GSE329024"]
gse329024_rel <- rel_all[dataset_folder == "GSE329024"]
gse329024_ref <- reference_controls[dataset_folder == "GSE329024"]

add_check(
  "GSE329024",
  "four_design_matched_control_groups_present",
  setequal(unique(gse329024_pred$reference_control_group_id), gse329024_expected$reference_control_group_id),
  paste(sort(unique(gse329024_pred$reference_control_group_id)), collapse = ";")
)

for (i in seq_len(nrow(gse329024_expected))) {
  exp_row <- gse329024_expected[i]
  group_id <- exp_row$reference_control_group_id
  design_token <- exp_row$design_token
  pred_group <- gse329024_pred[reference_control_group_id == group_id]
  rel_group <- gse329024_rel[reference_control_group_id == group_id]
  ref_group <- gse329024_ref[reference_control_group_id == group_id]
  sample_text <- paste(pred_group$treatment, pred_group$sample_title, pred_group$sample_key, sep = " ")

  add_check(
    "GSE329024",
    paste0(design_token, "_samples_are_not_mixed"),
    uniqueN(pred_group$sample_key) == exp_row$expected_samples &&
      length(sample_text) > 0 && all(grepl(design_token, sample_text, fixed = TRUE)),
    paste0("samples=", uniqueN(pred_group$sample_key), ";groups=", paste(sort(unique(pred_group$reference_control_group_id)), collapse = ";"))
  )
  add_check(
    "GSE329024",
    paste0(design_token, "_uses_only_its_matched_controls"),
    uniqueN(ref_group$sample_key) == exp_row$expected_controls &&
      nrow(rel_group) > 0 &&
      all(rel_group$matched_control_method == group_id) &&
      all(rel_group$n_matched_controls == exp_row$expected_controls),
    paste0(
      "controls=", uniqueN(ref_group$sample_key),
      ";methods=", paste(sort(unique(rel_group$matched_control_method)), collapse = ";"),
      ";n_matched=", paste(sort(unique(rel_group$n_matched_controls)), collapse = ";")
    )
  )
}

cra_pred <- pred_all[dataset_folder == "CRA017492"]
cra_rel <- rel_all[dataset_folder == "CRA017492"]
cra_young <- cra_pred[grepl("young 2 months control", treatment, ignore.case = TRUE)]
cra_aged_water <- cra_pred[grepl("aged 12 months water control", treatment, ignore.case = TRUE)]
cra_nmn_rel <- cra_rel[grepl("aged 12 months plus NMN", treatment, ignore.case = TRUE)]
add_check("CRA017492", "three_young_reference_controls", uniqueN(cra_young$sample_key) == 3, paste0("n=", uniqueN(cra_young$sample_key)))
add_check("CRA017492", "four_aged_water_controls", uniqueN(cra_aged_water$sample_key) == 4, paste0("n=", uniqueN(cra_aged_water$sample_key)))
add_check("CRA017492", "five_aged_NMN_samples", uniqueN(cra_nmn_rel$sample_key) == 5, paste0("n=", uniqueN(cra_nmn_rel$sample_key)))
add_check("CRA017492", "young_and_aged_controls_have_distinct_groups",
  all(c("CRA017492_young_control_reference", "CRA017492_aged_water_controls") %in% unique(cra_pred$reference_control_group_id)),
  paste(sort(unique(cra_pred$reference_control_group_id)), collapse = ";"))
add_check("CRA017492", "NMN_matches_aged_water_controls_only",
  nrow(cra_nmn_rel) == 10 && all(cra_nmn_rel$matched_control_method == "CRA017492_aged_water_controls"),
  paste0("rows=", nrow(cra_nmn_rel), ";methods=", paste(unique(cra_nmn_rel$matched_control_method), collapse = ";")))
add_check("CRA017492", "NMN_rows_have_four_aged_controls",
  nrow(cra_nmn_rel) == 10 && all(cra_nmn_rel$n_matched_controls == 4),
  paste0("rows=", nrow(cra_nmn_rel), ";n=", paste(unique(cra_nmn_rel$n_matched_controls), collapse = ";")))

gse279446_pred <- pred_all[dataset_folder == "GSE279446"]
gse279446_rel <- rel_all[dataset_folder == "GSE279446"]
add_check("GSE279446", "all_24_samples_projected", uniqueN(gse279446_pred$sample_key) == 24, paste0("n=", uniqueN(gse279446_pred$sample_key)))
add_check("GSE279446", "adult_and_young_cohorts_present", identical(sort(unique(na.omit(gse279446_pred$gse279446_cohort))), c("adult_cohort", "young_cohort")), paste(sort(unique(na.omit(gse279446_pred$gse279446_cohort))), collapse = ";"))
adult_mut <- gse279446_rel[gse279446_cohort == "adult_cohort" & projection_role == "external_projected_treatment"]
young_mut <- gse279446_rel[gse279446_cohort == "young_cohort" & projection_role == "external_projected_treatment"]
add_check("GSE279446", "adult_mutants_match_adult_WT_only", nrow(adult_mut) == 16 && all(adult_mut$matched_control_method == "GSE279446_adult_WT_ovary_controls"), paste0("rows=", nrow(adult_mut), ";methods=", paste(unique(adult_mut$matched_control_method), collapse = ";")))
add_check("GSE279446", "young_mutants_match_young_WT_only", nrow(young_mut) == 16 && all(young_mut$matched_control_method == "GSE279446_young_WT_controls"), paste0("rows=", nrow(young_mut), ";methods=", paste(unique(young_mut$matched_control_method), collapse = ";")))
add_check("GSE279446", "Zeb1_WT_not_heterozygous", nrow(gse279446_pred[grepl("Zeb1_WT|WT_ovary|WT ovary", paste(sample_key, treatment, genotype, sample_title, sep = " "), ignore.case = TRUE) & grepl("heterozygous", projection_class, ignore.case = TRUE)]) == 0, "WT classification checked")

gse279573_pred <- pred_all[dataset_folder == "GSE279573"]
gse279573_rel_treat <- rel_all[dataset_folder == "GSE279573" & projection_role == "external_projected_treatment"]
add_check("GSE279573", "seven_samples_in_predictions_all", uniqueN(gse279573_pred$sample_key) == 7, paste0("n=", uniqueN(gse279573_pred$sample_key)))
add_check("GSE279573", "two_WT_reference_controls", uniqueN(gse279573_pred[matched_reference_control_flag == TRUE]$sample_key) == 2, paste0("n=", uniqueN(gse279573_pred[matched_reference_control_flag == TRUE]$sample_key)))
add_check("GSE279573", "five_treatment_samples", uniqueN(gse279573_pred[projection_role == "external_projected_treatment"]$sample_key) == 5, paste0("n=", uniqueN(gse279573_pred[projection_role == "external_projected_treatment"]$sample_key)))
add_check("GSE279573", "treatment_rows_have_two_controls", nrow(gse279573_rel_treat) == 10 && all(gse279573_rel_treat$n_matched_controls == 2), paste0("rows=", nrow(gse279573_rel_treat), ";n=", paste(unique(gse279573_rel_treat$n_matched_controls), collapse = ";")))

gse240_pred <- pred_all[dataset_folder == "GSE240547"]
gse240_hfhs <- rel_all[dataset_folder == "GSE240547" & projection_class == "diet_intervention"]
add_check("GSE240547", "four_control_diet_samples", uniqueN(gse240_pred[projection_class == "control_diet"]$sample_key) == 4, paste0("n=", uniqueN(gse240_pred[projection_class == "control_diet"]$sample_key)))
add_check("GSE240547", "four_HFHS_samples", uniqueN(gse240_pred[projection_class == "diet_intervention"]$sample_key) == 4, paste0("n=", uniqueN(gse240_pred[projection_class == "diet_intervention"]$sample_key)))
add_check("GSE240547", "HFHS_rows_have_four_controls", nrow(gse240_hfhs) == 8 && all(gse240_hfhs$n_matched_controls == 4), paste0("rows=", nrow(gse240_hfhs), ";n=", paste(unique(gse240_hfhs$n_matched_controls), collapse = ";")))

prjna_pred <- pred_all[dataset_folder == "PRJNA1076509"]
prjna_rel <- rel_all[dataset_folder == "PRJNA1076509"]
add_check("PRJNA1076509", "fourteen_samples_detected", uniqueN(prjna_pred$sample_key) == 14, paste0("n=", uniqueN(prjna_pred$sample_key)))
add_check("PRJNA1076509", "twenty_eight_relative_rows", nrow(prjna_rel) == 28, paste0("rows=", nrow(prjna_rel)))
add_check("PRJNA1076509", "all_rows_use_young_donor_method", nrow(prjna_rel) == 28 && all(prjna_rel$matched_control_method == "PRJNA1076509_young_female_donor_microbiota"), paste(unique(prjna_rel$matched_control_method), collapse = ";"))
add_check("PRJNA1076509", "all_rows_have_seven_controls", nrow(prjna_rel) == 28 && all(prjna_rel$n_matched_controls == 7), paste(unique(prjna_rel$n_matched_controls), collapse = ";"))

prjna863443_pred <- pred_all[dataset_folder == "PRJNA863443"]
prjna863443_rel <- rel_all[dataset_folder == "PRJNA863443"]
add_check("PRJNA863443", "thirty_three_samples_detected", uniqueN(prjna863443_pred$sample_key) == 33, paste0("n=", uniqueN(prjna863443_pred$sample_key)))
add_check("PRJNA863443", "sixty_six_prediction_rows", nrow(prjna863443_pred) == 66, paste0("rows=", nrow(prjna863443_pred)))
vcd_treat <- prjna863443_rel[grepl("^PRJNA863443_VCD_vs_CTL_", comparison_id %||% reference_control_group_id) & projection_class == "intervention"]
fox_het <- prjna863443_rel[grepl("^PRJNA863443_Foxl2_het_vs_WT_", comparison_id %||% reference_control_group_id) & projection_class == "genetic_model"]
natural <- prjna863443_rel[(comparison_id %||% reference_control_group_id) == "PRJNA863443_natural_aging_reference"]
add_check("PRJNA863443", "VCD_treatment_rows_have_two_controls", nrow(vcd_treat) == 16 && all(vcd_treat$n_matched_controls == 2), paste0("rows=", nrow(vcd_treat), ";n=", paste(unique(vcd_treat$n_matched_controls), collapse = ";")))
add_check("PRJNA863443", "Foxl2_het_rows_have_expected_controls", nrow(fox_het) == 10 && all((grepl("_young$", fox_het$comparison_id) & fox_het$n_matched_controls == 2) | (grepl("_mid_age$", fox_het$comparison_id) & fox_het$n_matched_controls == 3)), paste0("rows=", nrow(fox_het), ";n=", paste(unique(fox_het$n_matched_controls), collapse = ";")))
add_check("PRJNA863443", "natural_ageing_no_relative_delta", nrow(natural) == 14 && all(natural$matched_control_method == "no_matched_control") && all(!is.finite(natural$relative_delta_age_week)), paste0("rows=", nrow(natural), ";finite=", sum(is.finite(natural$relative_delta_age_week))))

gse228_pred <- pred_all[dataset_folder == "GSE228975"]
gse228_rel <- rel_all[dataset_folder == "GSE228975"]
add_check("GSE228975", "primary_external_projection_include_false", nrow(gse228_pred) > 0 && all(gse228_pred$primary_external_projection_include == FALSE), paste(unique(gse228_pred$primary_external_projection_include), collapse = ";"))
add_check("GSE228975", "no_relative_delta", nrow(gse228_rel) > 0 && all(!is.finite(gse228_rel$relative_delta_age_week)), paste0("non_na=", sum(is.finite(gse228_rel$relative_delta_age_week))))
add_check("GSE228975", "excluded_from_primary_relative_file", nrow(rel_primary[dataset_folder == "GSE228975"]) == 0, paste0("rows=", nrow(rel_primary[dataset_folder == "GSE228975"])))

results <- rbindlist(checks, use.names = TRUE, fill = TRUE)
failures <- results[passed != TRUE]
fwrite(results, file.path(out_dir, "matched_control_verification_results.tsv"), sep = "\t")
fwrite(failures, file.path(out_dir, "matched_control_verification_failures.tsv"), sep = "\t")

status_for <- function(scope) if (all(results[scope == scope]$passed)) "PASS" else "FAIL"
summary_lines <- c(
  "# Matched-control verification summary",
  "",
  paste0("- All projected/reference sample count: ", uniqueN(pred_all$sample_key), "."),
  paste0("- Analysis-domain sample count: ", analysis_domain_count, "."),
  paste0("- Developmental/OOD placement sample count: ", developmental_ood_count, "."),
  paste0("- Legacy primary projection file sample count: ", uniqueN(pred_primary$sample_key), "."),
  paste0("- Prediction row count: ", nrow(pred_all), "."),
  paste0("- Relative-control row count: ", nrow(rel_all), "."),
  paste0("- Invalid relative-control row count: ", nrow(unreliable) + nrow(invalid_relative), "."),
  paste0("- GSE194196 matched-control status: ", status_for("GSE194196"), "."),
  paste0("- CRA017492 matched-control status: ", status_for("CRA017492"), "."),
  paste0("- GSE252791 matched-control status: ", status_for("GSE252791"), "."),
  paste0("- GSE279446 matched-control status: ", status_for("GSE279446"), "."),
  paste0("- GSE279573 matched-control status: ", status_for("GSE279573"), "."),
  paste0("- GSE329024 four-stratum matched-control status: ", status_for("GSE329024"), "."),
  paste0("- GSE240547 matched-control status: ", status_for("GSE240547"), "."),
  paste0("- PRJNA1076509 matched-control status: ", status_for("PRJNA1076509"), "."),
  paste0("- PRJNA863443 matched-control status: ", status_for("PRJNA863443"), "."),
  paste0("- GSE228975 OOD status: ", status_for("GSE228975"), "."),
  paste0("- Overall status: ", if (nrow(failures) == 0) "MATCHED_CONTROL_VERIFICATION_PASS" else "MATCHED_CONTROL_VERIFICATION_FAIL", ".")
)
writeLines(summary_lines, file.path(out_dir, "matched_control_verification_summary.md"))

cat("All projected/reference sample count:", uniqueN(pred_all$sample_key), "\n")
cat("Analysis-domain sample count:", analysis_domain_count, "\n")
cat("Developmental/OOD placement sample count:", developmental_ood_count, "\n")
cat("Legacy primary projection file sample count:", uniqueN(pred_primary$sample_key), "\n")
cat("Prediction row count:", nrow(pred_all), "\n")
cat("Relative-control row count:", nrow(rel_all), "\n")
cat("Invalid relative-control row count:", nrow(unreliable) + nrow(invalid_relative), "\n")
cat("GSE194196 matched-control status:", status_for("GSE194196"), "\n")
cat("CRA017492 matched-control status:", status_for("CRA017492"), "\n")
cat("GSE252791 matched-control status:", status_for("GSE252791"), "\n")
cat("GSE279446 matched-control status:", status_for("GSE279446"), "\n")
cat("GSE279573 matched-control status:", status_for("GSE279573"), "\n")
cat("GSE329024 four-stratum matched-control status:", status_for("GSE329024"), "\n")
cat("GSE240547 matched-control status:", status_for("GSE240547"), "\n")
cat("PRJNA1076509 matched-control status:", status_for("PRJNA1076509"), "\n")
cat("PRJNA863443 matched-control status:", status_for("PRJNA863443"), "\n")
cat("GSE228975 OOD status:", status_for("GSE228975"), "\n")

if (nrow(failures) == 0) {
  cat("MATCHED_CONTROL_VERIFICATION_PASS\n")
} else {
  cat("MATCHED_CONTROL_VERIFICATION_FAIL\n")
  stop("Matched-control verification failed; see matched_control_verification_failures.tsv", call. = FALSE)
}
