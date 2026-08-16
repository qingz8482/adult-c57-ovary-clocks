#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
})

options(stringsAsFactors = FALSE)

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
data_input_dir <- normalizePath(file.path(project_dir, "..", "data_input"), mustWork = TRUE)
results_dir <- Sys.getenv("RESULTS_ROOT", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = TRUE)
if (grepl("results_groupmapped_v2", results_dir, fixed = TRUE) || grepl("_v2", basename(results_dir), fixed = TRUE)) {
  stop("Refusing to project PRJNA863443 against a v2 results root: ", results_dir, call. = FALSE)
}

projection_dir <- file.path(results_dir, "external_projection")
dataset_dir <- file.path(projection_dir, "PRJNA863443")
plot_dir <- file.path(dataset_dir, "plots")
dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

dataset_input_dir <- file.path(data_input_dir, "aggregated_pseudobulk", "PRJNA863443")
genelevel_dir <- file.path(dataset_input_dir, "input_genelevel")
metadata_path <- file.path(dataset_input_dir, "PRJNA863443_sample_metadata.tsv")
counts_path <- file.path(genelevel_dir, "sample_gene_counts_after_expression_filter.tsv")
log_path <- file.path(genelevel_dir, "sample_gene_logcpm_after_expression_filter.tsv")

retained_clocks <- data.table(
  model_family = c("elasticnet", "pca_glmnet"),
  selector = c("hvg_baseline", "version_b_age_dataset_penalized"),
  model_path = file.path(
    results_dir,
    "model_artifacts",
    "deployment",
    c(
      "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds",
      "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
    )
  )
)

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

ensure_cols <- function(dt, cols) {
  for (cc in cols) if (!cc %in% names(dt)) dt[, (cc) := NA]
  dt
}

read_expression_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing expression matrix: ", path, call. = FALSE)
  dt <- fread(path)
  gene_col <- names(dt)[1]
  genes <- na_chr(dt[[gene_col]])
  dt[, (gene_col) := NULL]
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- genes
  keep <- !is.na(rownames(mat)) & rownames(mat) != "" & !duplicated(rownames(mat))
  mat[keep, , drop = FALSE]
}

counts_to_logcpm <- function(counts) {
  lib_size <- colSums(counts, na.rm = TRUE)
  if (any(!is.finite(lib_size) | lib_size <= 0)) {
    stop("PRJNA863443 raw count matrix contains samples with zero/non-finite library size.", call. = FALSE)
  }
  log1p(t(t(counts) / lib_size * 1e6))
}

extract_calibration <- function(obj) {
  candidates <- c("calibration", "calibration_model", "calibration_fit", "calibration_params", "calibration_parameters")
  hit <- candidates[candidates %in% names(obj)]
  if (length(hit) > 0) {
    cal <- obj[[hit[[1]]]]
    if (inherits(cal, "lm")) return(list(type = "lm", fit = cal))
    if (is.list(cal)) {
      intercept <- cal$intercept %||% cal$calibration_intercept %||% cal$a
      slope <- cal$slope %||% cal$calibration_slope %||% cal$b
      if (length(intercept) == 1 && length(slope) == 1 && is.finite(intercept) && is.finite(slope)) {
        return(list(type = "linear", intercept = as.numeric(intercept), slope = as.numeric(slope)))
      }
    }
    if (is.numeric(cal) && all(c("intercept", "slope") %in% names(cal))) {
      return(list(type = "linear", intercept = unname(cal[["intercept"]]), slope = unname(cal[["slope"]])))
    }
  }
  intercept <- obj$calibration_intercept %||% obj$intercept
  slope <- obj$calibration_slope %||% obj$slope
  if (length(intercept) == 1 && length(slope) == 1 && is.finite(intercept) && is.finite(slope)) {
    return(list(type = "linear", intercept = as.numeric(intercept), slope = as.numeric(slope)))
  }
  NULL
}

apply_calibration <- function(raw_pred, cal) {
  if (is.null(cal)) stop("Cannot apply missing calibration.", call. = FALSE)
  if (identical(cal$type, "lm")) return(as.numeric(predict(cal$fit, newdata = data.frame(raw_train = raw_pred))))
  if (identical(cal$type, "linear")) return(cal$intercept + cal$slope * raw_pred)
  stop("Unsupported calibration object type: ", cal$type, call. = FALSE)
}

normalize_artifact_object <- function(obj, model_family) {
  if (!is.null(obj$schema_version) && obj$schema_version == "integrated_ovary_clock_artifact_v1") {
    out <- list(
      model = obj$model$glmnet_fit,
      genes = obj$selected_genes,
      center = obj$gene_center,
      scale = obj$gene_scale,
      training_scope = obj$training_scope,
      artifact_scope = obj$artifact_scope,
      expression_scale = obj$expression_scale,
      calibration = obj$calibration,
      source_artifact = obj
    )
    if (model_family == "pca_glmnet") {
      out$pca <- obj$model$pca_fit
      out$pc_count <- obj$model$selected_pc_count
      out$pc_coefficient_table <- obj$model$pc_coefficient_table
      out$calibrated_intercept <- obj$model$calibrated_intercept
    }
    return(out)
  }
  obj
}

audit_artifacts <- function(clock_dt) {
  rows <- list()
  objects <- list()
  for (i in seq_len(nrow(clock_dt))) {
    family <- clock_dt$model_family[i]
    selector <- clock_dt$selector[i]
    path <- clock_dt$model_path[i]
    add <- function(requirement, present, detail) {
      rows[[length(rows) + 1L]] <<- data.table(
        model_family = family,
        selector = selector,
        requirement = requirement,
        present = isTRUE(present),
        detail = as.character(detail)
      )
    }
    add("deployment model RDS", file.exists(path), path)
    if (!file.exists(path)) next
    raw_obj <- readRDS(path)
    obj <- normalize_artifact_object(raw_obj, family)
    cal <- extract_calibration(obj)
    objects[[paste(family, selector, sep = "\r")]] <- obj
    add("selected genes", !is.null(obj$genes) && length(obj$genes) > 0, "selected_genes")
    add("gene center", !is.null(obj$center) && length(obj$center) > 0, "gene_center")
    add("gene scale", !is.null(obj$scale) && length(obj$scale) > 0, "gene_scale")
    add("calibration intercept", !is.null(cal) && !is.null(cal$intercept) && is.finite(cal$intercept), "calibration$intercept")
    add("calibration slope", !is.null(cal) && !is.null(cal$slope) && is.finite(cal$slope), "calibration$slope")
    add("model object", !is.null(obj$model), "glmnet model")
    add("artifact scope deployment", !is.null(raw_obj$artifact_scope) && grepl("deployment", raw_obj$artifact_scope), raw_obj$artifact_scope %||% NA_character_)
    add("training scope integrated", identical(raw_obj$training_scope, "integrated_bulk_pseudobulk"), raw_obj$training_scope %||% NA_character_)
    if (family == "pca_glmnet") {
      add("PCA fit", !is.null(obj$pca) && "rotation" %in% names(obj$pca), "pca_fit")
      add("selected PC count", !is.null(obj$pc_count) && is.finite(obj$pc_count) && obj$pc_count > 0, obj$pc_count %||% NA_character_)
    }
  }
  audit <- rbindlist(rows, fill = TRUE)
  fwrite(audit, file.path(dataset_dir, "prjna863443_artifact_audit.tsv"), sep = "\t")
  missing <- audit[present != TRUE]
  if (nrow(missing) > 0) {
    stop("Missing/invalid deployment artifact requirements; see prjna863443_artifact_audit.tsv.", call. = FALSE)
  }
  list(audit = audit, objects = objects)
}

model_matrix_for_genes <- function(expr_pieces, meta_dt, model_genes, center, scale) {
  out <- matrix(NA_real_, nrow = nrow(meta_dt), ncol = length(model_genes))
  rownames(out) <- meta_dt$sample_key
  colnames(out) <- model_genes
  for (piece_name in names(expr_pieces)) {
    mat <- expr_pieces[[piece_name]]
    keys <- intersect(colnames(mat), meta_dt$sample_key)
    if (length(keys) == 0) next
    genes_present <- intersect(model_genes, rownames(mat))
    if (length(genes_present) > 0) out[keys, genes_present] <- t(mat[genes_present, keys, drop = FALSE])
  }
  missing <- is.na(out)
  raw <- out
  for (gene in model_genes) raw[missing[, gene], gene] <- center[[gene]]
  scaled <- sweep(raw, 2, center[model_genes], "-")
  scaled <- sweep(scaled, 2, scale[model_genes], "/")
  scaled[!is.finite(scaled)] <- 0
  coverage <- data.table(
    sample_key = rownames(out),
    n_model_genes = length(model_genes),
    n_observed_model_genes = rowSums(!missing),
    gene_coverage_fraction = rowSums(!missing) / length(model_genes),
    missing_model_genes = apply(missing, 1, function(x) paste(model_genes[x], collapse = ";"))
  )
  coverage[, projection_gene_coverage_flag := fifelse(gene_coverage_fraction < 0.80, "low_coverage", "ok")]
  list(x = scaled, coverage = coverage)
}

robust_mahalanobis <- function(x, ref_x) {
  if (nrow(ref_x) < 3 || ncol(ref_x) < 1) return(rep(NA_real_, nrow(x)))
  center <- colMeans(ref_x, na.rm = TRUE)
  cov_mat <- stats::cov(ref_x, use = "pairwise.complete.obs")
  cov_mat <- as.matrix(cov_mat)
  cov_mat[!is.finite(cov_mat)] <- 0
  diag(cov_mat) <- diag(cov_mat) + 1e-6
  as.numeric(stats::mahalanobis(x, center = center, cov = cov_mat))
}

distance_percentile <- function(distance, ref_distance) {
  vapply(distance, function(d) {
    if (!is.finite(d) || all(!is.finite(ref_distance))) return(NA_real_)
    mean(ref_distance <= d, na.rm = TRUE)
  }, numeric(1))
}

clock_predict <- function(clock, obj, cal, expr_pieces, meta_apply, meta_reference) {
  genes <- as.character(obj$genes)
  center <- obj$center[genes]
  scale <- obj$scale[genes]
  scale[!is.finite(scale) | scale == 0] <- 1
  apply_x <- model_matrix_for_genes(expr_pieces, meta_apply, genes, center, scale)
  ref_x <- model_matrix_for_genes(expr_pieces, meta_reference, genes, center, scale)
  if (clock$model_family == "elasticnet") {
    raw_apply <- as.numeric(predict(obj$model, newx = apply_x$x))
    raw_ref <- as.numeric(predict(obj$model, newx = ref_x$x))
    pred <- apply_calibration(raw_apply, cal)
    pred_ref <- apply_calibration(raw_ref, cal)
    centroid <- colMeans(ref_x$x, na.rm = TRUE)
    dist_apply <- sqrt(rowSums(sweep(apply_x$x, 2, centroid, "-")^2))
    dist_ref <- sqrt(rowSums(sweep(ref_x$x, 2, centroid, "-")^2))
    pctl <- distance_percentile(dist_apply, dist_ref)
    pred_dt <- data.table(
      sample_key = meta_apply$sample_key,
      model_family = clock$model_family,
      selector = clock$selector,
      raw_predicted_age_week = raw_apply,
      predicted_age_week = pred,
      projection_distance = dist_apply,
      projection_distance_percentile = pctl,
      OOD_flag = fifelse(pctl > 0.99, "OOD_gt_99th", fifelse(pctl > 0.95, "OOD_gt_95th", "in_reference_range"))
    )
    ref_dt <- data.table(
      sample_key = meta_reference$sample_key,
      model_family = clock$model_family,
      selector = clock$selector,
      predicted_age_week = pred_ref,
      delta_age_week = pred_ref - safe_num(meta_reference$age_week),
      projection_distance = dist_ref
    )
    return(list(predictions = merge(pred_dt, apply_x$coverage, by = "sample_key"), reference_predictions = ref_dt, pca_scores = NULL))
  }
  pc_count <- as.integer(obj$pc_count)
  pc_apply_all <- predict(obj$pca, newdata = apply_x$x)
  pc_ref_all <- predict(obj$pca, newdata = ref_x$x)
  pc_cols <- paste0("PC", seq_len(pc_count))
  pc_apply <- pc_apply_all[, seq_len(pc_count), drop = FALSE]
  pc_ref <- pc_ref_all[, seq_len(pc_count), drop = FALSE]
  colnames(pc_apply) <- pc_cols
  colnames(pc_ref) <- pc_cols
  raw_apply <- as.numeric(predict(obj$model, newx = pc_apply))
  raw_ref <- as.numeric(predict(obj$model, newx = pc_ref))
  pred <- apply_calibration(raw_apply, cal)
  pred_ref <- apply_calibration(raw_ref, cal)
  dist_apply <- robust_mahalanobis(pc_apply, pc_ref)
  dist_ref <- robust_mahalanobis(pc_ref, pc_ref)
  pctl <- distance_percentile(dist_apply, dist_ref)
  pred_dt <- data.table(
    sample_key = meta_apply$sample_key,
    model_family = clock$model_family,
    selector = clock$selector,
    raw_predicted_age_week = raw_apply,
    predicted_age_week = pred,
    projection_distance = dist_apply,
    projection_distance_percentile = pctl,
    OOD_flag = fifelse(pctl > 0.99, "OOD_gt_99th", fifelse(pctl > 0.95, "OOD_gt_95th", "in_reference_range"))
  )
  beta <- as.numeric(coef(obj$model))[-1]
  beta <- beta[seq_len(pc_count)]
  score_dt <- as.data.table(pc_apply)
  score_dt[, sample_key := meta_apply$sample_key]
  setcolorder(score_dt, "sample_key")
  for (j in seq_len(pc_count)) score_dt[, paste0(pc_cols[[j]], "_contribution_week") := pc_apply[, j] * beta[[j]]]
  score_dt[, `:=`(
    model_family = clock$model_family,
    selector = clock$selector,
    total_calibrated_predicted_age_week = pred,
    delta_age_week = pred - safe_num(meta_apply$age_week)
  )]
  ref_dt <- data.table(
    sample_key = meta_reference$sample_key,
    model_family = clock$model_family,
    selector = clock$selector,
    predicted_age_week = pred_ref,
    delta_age_week = pred_ref - safe_num(meta_reference$age_week),
    projection_distance = dist_ref
  )
  list(predictions = merge(pred_dt, apply_x$coverage, by = "sample_key"), reference_predictions = ref_dt, pca_scores = score_dt)
}

age_group_from_week <- function(age_week) {
  fifelse(age_week < 25, "young", fifelse(age_week < 65, "mid_age", "old"))
}

natural_age_reference_label <- function(age_group) {
  fifelse(
    age_group == "young", "natural ageing young adult reference",
    fifelse(age_group == "mid_age", "natural ageing mid-age reference", "natural ageing old reference")
  )
}

natural_aging_plot_group <- function(age_week) {
  fifelse(
    age_week < 25,
    "young_13w-17.4w",
    fifelse(age_week < 86, "mid_76w", "old_86w-87w")
  )
}

make_prjna_metadata <- function(meta) {
  meta <- copy(meta)
  meta[, age_week := safe_num(age_weeks)]
  meta[, age_group := age_group_from_week(age_week)]
  meta[, natural_age_reference_group := natural_age_reference_label(age_group)]
  meta[, natural_aging_plot_group := natural_aging_plot_group(age_week)]
  meta[, dataset_folder := "PRJNA863443"]
  meta[, modality := "scrna_pseudobulk"]
  meta[, sample_key := paste("scrna_pseudobulk", dataset_folder, sample_id, sep = "__")]
  meta[, accession := run_accession]
  meta[, bio_age_bin := fifelse(age_week < 16, "young_adult",
    fifelse(age_week < 40, "mature_adult",
      fifelse(age_week < 65, "advanced_maternal_age", "post_reproductive_ageing")
    )
  )]
  meta[, sample_context := paste(model, condition, genotype, treatment, sep = " | ")]
  meta[, projection_class := "other"]
  meta[, projection_role := "external_treatment"]
  meta[, projection_reason := "PRJNA863443_external_pseudobulk_projection"]
  meta[, projection_group := NA_character_]
  meta[, comparison_id := NA_character_]
  meta[, genotype_class := NA_character_]
  meta[, cohort := age_group]
  meta[, control_group_label := NA_character_]
  meta[, matched_reference_control_flag := FALSE]
  meta[, matched_reference_control_source := "none"]
  meta[, primary_external_projection_include := TRUE]
  meta[, primary_effect_interpretation_include := TRUE]
  meta[, matched_control_allowed := TRUE]
  meta[, reference_control_group_id := NA_character_]
  meta[, reference_control_group_label := NA_character_]
  meta[, gse279446_cohort := NA_character_]
  meta[, gse279446_platform := NA_character_]
  meta[, projection_dataset_note := "PRJNA863443 scRNA-seq pseudobulk external projection; no fitting/refitting."]
  meta[, control_pool_source := "none"]
  meta[, control_candidate_status := "rejected_or_projection_sample"]
  meta[, control_candidate_note := projection_reason]
  meta[, time_id := paste0(injection_age, "_", post_injection_days, "d")]

  is_natural <- grepl("Aging", meta$model, ignore.case = TRUE) & meta$treatment == "none" & meta$genotype == "WT"
  is_vcd_ctl <- meta$model == "VCD" & (meta$condition == "CTL" | meta$treatment == "vehicle")
  is_vcd <- meta$model == "VCD" & (meta$condition == "VCD" | meta$treatment == "VCD")
  is_fox_wt <- meta$genotype == "Foxl2+/+" | meta$condition == "Foxl2_wt"
  is_fox_het <- meta$genotype == "Foxl2+/-" | meta$condition == "Foxl2_het"
  meta[is_natural, `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = natural_age_reference_group,
    comparison_id = "PRJNA863443_natural_aging_reference",
    primary_effect_interpretation_include = FALSE,
    matched_control_allowed = FALSE,
    reference_control_group_id = "no_matched_control",
    reference_control_group_label = "No matched control",
    genotype_class = "WT",
    control_group_label = natural_age_reference_group,
    projection_reason = "PRJNA863443_natural_ageing_reference_placement"
  )]
  meta[is_vcd_ctl, `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = paste0("CTL_", time_id),
    comparison_id = paste0("PRJNA863443_VCD_vs_CTL_", time_id),
    reference_control_group_id = paste0("PRJNA863443_VCD_vs_CTL_", time_id),
    reference_control_group_label = paste0("PRJNA863443 CTL ", time_id, " controls"),
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    genotype_class = "WT",
    control_group_label = paste0("CTL_", time_id),
    control_pool_source = "projection_internal_control",
    control_candidate_status = "accepted_projection_internal_control",
    control_candidate_note = "PRJNA863443 VCD same injection-age/post-injection controls"
  )]
  meta[is_vcd, `:=`(
    projection_class = "intervention",
    projection_role = "external_treatment",
    projection_group = paste0("VCD_", time_id),
    comparison_id = paste0("PRJNA863443_VCD_vs_CTL_", time_id),
    reference_control_group_id = paste0("PRJNA863443_VCD_vs_CTL_", time_id),
    reference_control_group_label = paste0("PRJNA863443 CTL ", time_id, " controls"),
    genotype_class = "WT",
    control_group_label = paste0("CTL_", time_id),
    projection_reason = "PRJNA863443_VCD_treatment"
  )]
  meta[is_fox_wt, `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = paste0("Foxl2+/+_", age_group),
    comparison_id = paste0("PRJNA863443_Foxl2_het_vs_WT_", age_group),
    reference_control_group_id = paste0("PRJNA863443_Foxl2_het_vs_WT_", age_group),
    reference_control_group_label = paste0("PRJNA863443 Foxl2+/+ ", age_group, " controls"),
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    genotype_class = "WT",
    control_group_label = paste0("Foxl2+/+_", age_group),
    control_pool_source = "projection_internal_control",
    control_candidate_status = "accepted_projection_internal_control",
    control_candidate_note = "PRJNA863443 Foxl2 same age-group WT controls"
  )]
  meta[is_fox_het, `:=`(
    projection_class = "genetic_model",
    projection_role = "external_treatment",
    projection_group = paste0("Foxl2+/-_", age_group),
    comparison_id = paste0("PRJNA863443_Foxl2_het_vs_WT_", age_group),
    reference_control_group_id = paste0("PRJNA863443_Foxl2_het_vs_WT_", age_group),
    reference_control_group_label = paste0("PRJNA863443 Foxl2+/+ ", age_group, " controls"),
    genotype_class = "heterozygous",
    control_group_label = paste0("Foxl2+/+_", age_group),
    projection_reason = "PRJNA863443_Foxl2_heterozygous"
  )]
  meta[]
}

compute_prjna_relative <- function(pred) {
  rel <- copy(pred)
  rel[, `:=`(
    control_expected_predicted_age_week = NA_real_,
    control_expected_delta_age_week = NA_real_,
    relative_delta_age_week = NA_real_,
    matched_control_method = "no_matched_control",
    matched_control_source = "none",
    n_matched_controls = 0L,
    matched_control_sample_keys = NA_character_,
    relative_control_reliability = "not_available",
    relative_control_warning = "no matched control; placement only"
  )]
  for (key in unique(rel[matched_control_allowed == TRUE]$reference_control_group_id)) {
    if (is.na(key) || key == "" || key == "no_matched_control") next
    idx <- rel$matched_control_allowed == TRUE & rel$reference_control_group_id == key
    ctrl <- rel[idx & matched_reference_control_flag == TRUE]
    if (nrow(ctrl) == 0) next
    ctrl_mean <- ctrl[, .(
      control_expected_predicted_age_week = mean(predicted_age_week, na.rm = TRUE),
      control_expected_delta_age_week = mean(delta_age_week, na.rm = TRUE),
      n_matched_controls = uniqueN(sample_key),
      matched_control_sample_keys = paste(sort(unique(sample_key)), collapse = ";")
    ), by = .(model_family, selector)]
    for (j in seq_len(nrow(ctrl_mean))) {
      idx_model <- idx &
        rel$model_family == ctrl_mean$model_family[j] &
        rel$selector == ctrl_mean$selector[j]
      rel[idx_model, `:=`(
        control_expected_predicted_age_week = ctrl_mean$control_expected_predicted_age_week[j],
        control_expected_delta_age_week = ctrl_mean$control_expected_delta_age_week[j],
        n_matched_controls = ctrl_mean$n_matched_controls[j],
        matched_control_sample_keys = ctrl_mean$matched_control_sample_keys[j]
      )]
    }
    rel[idx, `:=`(
      relative_delta_age_week = delta_age_week - control_expected_delta_age_week,
      matched_control_method = key,
      matched_control_source = "projection_internal_control",
      relative_control_reliability = fifelse(n_matched_controls >= 3, "valid_same_study_control_pool", "limited_control_pool"),
      relative_control_warning = NA_character_
    )]
  }
  rel[, `:=`(
    control_expected_pred_age = control_expected_predicted_age_week,
    control_expected_delta_age = control_expected_delta_age_week,
    relative_delta_age = relative_delta_age_week
  )]
  rel[matched_control_allowed == FALSE, `:=`(
    control_expected_predicted_age_week = NA_real_,
    control_expected_delta_age_week = NA_real_,
    relative_delta_age_week = NA_real_,
    control_expected_pred_age = NA_real_,
    control_expected_delta_age = NA_real_,
    relative_delta_age = NA_real_,
    matched_control_method = "no_matched_control",
    matched_control_source = "none",
    n_matched_controls = 0L,
    matched_control_sample_keys = NA_character_,
    relative_control_reliability = "not_available",
    relative_control_warning = "natural ageing reference; no matched-control relative effect"
  )]
  rel[]
}

append_replace_dataset <- function(filename, new_dt, dataset = "PRJNA863443") {
  path <- file.path(projection_dir, filename)
  if (file.exists(path)) {
    old <- fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
    old <- ensure_cols(old, names(new_dt))
    new_dt <- ensure_cols(copy(new_dt), names(old))
    old <- old[is.na(dataset_folder) | dataset_folder != dataset]
    out <- rbindlist(list(old[, names(old), with = FALSE], new_dt[, names(old), with = FALSE]), use.names = TRUE, fill = TRUE)
  } else {
    out <- copy(new_dt)
  }
  fwrite(out, path, sep = "\t")
  invisible(out)
}

update_dataset_availability <- function(prjna_pred, prjna_rel) {
  path <- file.path(projection_dir, "external_projection_dataset_control_availability.tsv")
  dt <- if (file.exists(path)) fread(path, sep = "\t", na.strings = c("", "NA", "NaN")) else data.table()
  dt <- dt[is.na(dataset_folder) | dataset_folder != "PRJNA863443"]
  methods <- sort(unique(prjna_rel$matched_control_method))
  row <- data.table(
    dataset_folder = "PRJNA863443",
    n_projected_samples_all = uniqueN(prjna_pred$sample_key),
    n_projected_samples_primary = uniqueN(prjna_pred[primary_external_projection_include == TRUE]$sample_key),
    n_control_meta_candidates = 0L,
    n_projection_internal_control_candidates = uniqueN(prjna_pred[matched_reference_control_flag == TRUE]$sample_key),
    has_same_dataset_control_meta = FALSE,
    has_projection_internal_control = TRUE,
    recommended_relative_method = "PRJNA863443_VCD_same_time_CTL;PRJNA863443_Foxl2_same_age_WT;natural_aging_no_matched_control",
    expected_relative_method = "VCD matched to CTL by injection_age/post_injection_days; Foxl2+/- matched to Foxl2+/+ by age_group; natural ageing no_matched_control",
    observed_relative_method = paste(methods, collapse = ";"),
    primary_external_projection_include = TRUE,
    dataset_note = "PRJNA863443 scRNA-seq pseudobulk: VCD and Foxl2 use projection-internal controls; natural ageing rows are placement-only references.",
    pass_fail_status = "PASS"
  )
  dt <- ensure_cols(dt, names(row))
  row <- ensure_cols(row, names(dt))
  out <- rbindlist(list(dt[, names(dt), with = FALSE], row[, names(dt), with = FALSE]), use.names = TRUE, fill = TRUE)
  fwrite(out, path, sep = "\t")
  invisible(out)
}

update_run_manifest_and_summary <- function() {
  pred_all <- fread(file.path(projection_dir, "external_projection_predictions_all.tsv"), sep = "\t", na.strings = c("", "NA", "NaN"))
  rel_all <- fread(file.path(projection_dir, "external_projection_relative_to_controls_all.tsv"), sep = "\t", na.strings = c("", "NA", "NaN"))
  pred_all[, primary_external_projection_include := as_bool(primary_external_projection_include)]
  invalid <- rel_all[
    (matched_control_method == "no_matched_control" & is.finite(safe_num(relative_delta_age_week))) |
      (safe_num(n_matched_controls) == 0 & is.finite(safe_num(relative_delta_age_week))) |
      (matched_control_allowed == FALSE & is.finite(safe_num(relative_delta_age_week))) |
      (primary_external_projection_include == FALSE & is.finite(safe_num(relative_delta_age_week)))
  ]
  manifest <- data.table(
    field = c(
      "run_label", "results_root", "output_dir", "prjna863443_appended_by",
      "n_projected_samples_all", "n_projected_samples_primary",
      "n_prediction_rows_all", "n_relative_control_rows_all",
      "n_invalid_relative_control_rows", "created_at"
    ),
    value = c(
      Sys.getenv("RUN_LABEL", unset = "v1"), results_dir, projection_dir,
      basename(script_file %||% "05c_project_prjna863443_pseudobulk_external.R"),
      uniqueN(pred_all$sample_key),
      uniqueN(pred_all[primary_external_projection_include == TRUE]$sample_key),
      nrow(pred_all), nrow(rel_all), nrow(invalid), as.character(Sys.time())
    )
  )
  fwrite(manifest, file.path(projection_dir, "external_projection_run_manifest.tsv"), sep = "\t")
  summary_path <- file.path(projection_dir, "external_projection_summary.md")
  lines <- if (file.exists(summary_path)) readLines(summary_path, warn = FALSE) else c("# External projection summary")
  lines <- lines[!grepl("^## PRJNA863443$", lines) & !grepl("^[-] PRJNA863443", lines)]
  prjna <- pred_all[dataset_folder == "PRJNA863443"]
  cov <- prjna[, .(
    min_gene_coverage = round(min(gene_coverage_fraction, na.rm = TRUE), 4),
    median_gene_coverage = round(median(gene_coverage_fraction, na.rm = TRUE), 4)
  ), by = model_family]
  add <- c(
    "",
    "## PRJNA863443",
    paste0("- PRJNA863443 sample count: ", uniqueN(prjna$sample_key), "."),
    paste0("- PRJNA863443 prediction rows: ", nrow(prjna), "."),
    paste0("- Gene coverage by model: ", paste(paste0(cov$model_family, " min=", cov$min_gene_coverage, " median=", cov$median_gene_coverage), collapse = "; "), "."),
    "- VCD relative deltas use CTL controls matched by injection_age and post_injection_days.",
    "- Foxl2 relative deltas use Foxl2+/+ controls matched by age_group.",
    "- Natural ageing rows are projected as reference-control placement rows with no matched-control relative effect."
  )
  writeLines(c(lines, add), summary_path)
  invisible(invalid)
}

theme_prjna <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

save_prjna_plot <- function(p, stem, width = 10, height = 7) {
  ggsave(file.path(plot_dir, paste0(stem, ".pdf")), p, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(plot_dir, paste0(stem, ".png")), p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

load_deployment_training_expression <- function(meta_dt) {
  alignment <- fread(
    file.path(results_dir, "integrated_bulk_pseudobulk_sample_expression_alignment.tsv"),
    sep = "\t",
    na.strings = c("", "NA", "NaN")
  )
  pieces <- list()
  for (path in unique(meta_dt$expression_logcpm_path)) {
    path_meta <- meta_dt[expression_logcpm_path == path]
    mat <- read_expression_matrix(path)
    expression_cols <- alignment$expression_column[match(path_meta$sample_key, alignment$sample_key)]
    missing_alignment <- is.na(expression_cols) | !expression_cols %in% colnames(mat)
    if (any(missing_alignment)) {
      stop(
        "Deployment training expression alignment is missing samples in ", path, ": ",
        paste(path_meta$sample_key[missing_alignment], collapse = ", "),
        call. = FALSE
      )
    }
    mat <- mat[, expression_cols, drop = FALSE]
    colnames(mat) <- path_meta$sample_key
    pieces[[path]] <- mat
  }
  pieces
}

write_dataset_plots <- function(pred, rel, concordance, coverage, pca_scores) {
  deployment_meta <- as.data.table(readRDS(retained_clocks$model_path[[1]])$reference_metadata)
  deployment_expr <- load_deployment_training_expression(deployment_meta)
  deployment_fit_list <- vector("list", nrow(retained_clocks))
  for (i in seq_len(nrow(retained_clocks))) {
    clock <- retained_clocks[i]
    obj <- normalize_artifact_object(readRDS(clock$model_path), clock$model_family)
    deployment_fit <- clock_predict(
      clock,
      obj,
      extract_calibration(obj),
      deployment_expr,
      deployment_meta,
      deployment_meta
    )$predictions
    deployment_fit <- merge(
      deployment_fit,
      deployment_meta[, .(sample_key, age_week)],
      by = "sample_key"
    )
    deployment_fit_list[[i]] <- deployment_fit
  }
  locked_line <- rbindlist(deployment_fit_list, use.names = TRUE, fill = TRUE)
  locked_line_stats <- locked_line[, .(
    locked_mae_week = mean(abs(predicted_age_week - age_week), na.rm = TRUE),
    locked_r2 = summary(lm(predicted_age_week ~ age_week))$r.squared
  ), by = .(model_family, selector)]
  locked_line_stats[, locked_metric_label := sprintf(
    "Deployment training fit (in-sample)\nMAE = %.1f weeks\nR² = %.2f",
    locked_mae_week,
    locked_r2
  )]

  pred_plot <- copy(pred)
  pred_plot[, predicted_vs_age_color_group := fifelse(
    comparison_id == "PRJNA863443_natural_aging_reference",
    natural_aging_plot_group,
    projection_group
  )]
  pred_plot[, predicted_vs_age_color_group := factor(
    predicted_vs_age_color_group,
    levels = unique(c(
      "young_13w-17.4w", "mid_76w", "old_86w-87w",
      sort(setdiff(unique(predicted_vs_age_color_group), c("young_13w-17.4w", "mid_76w", "old_86w-87w")))
    ))
  )]
  p <- ggplot(pred_plot, aes(age_week, predicted_age_week, color = predicted_vs_age_color_group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_smooth(
      data = locked_line,
      aes(x = age_week, y = predicted_age_week, group = 1),
      inherit.aes = FALSE,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      linewidth = 0.9,
      linetype = "solid",
      color = "black"
    ) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_text(
      data = locked_line_stats,
      aes(x = -Inf, y = Inf, label = locked_metric_label),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.08,
      size = 3.8,
      lineheight = 1.05
    ) +
    facet_wrap(~model_family) +
    labs(x = "Chronological age (weeks)", y = "Predicted age (weeks)", color = "Projection group") +
    theme_prjna()
  save_prjna_plot(p, "01_prjna863443_predicted_vs_age", width = 12, height = 7)

  p <- ggplot(pred, aes(projection_group, delta_age_week, color = projection_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, color = "grey55") +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.85) +
    facet_wrap(~model_family, ncol = 1) +
    labs(x = "Projection group", y = "Raw delta age (weeks)", color = "Projection class") +
    theme_prjna()
  save_prjna_plot(p, "02_prjna863443_raw_delta_by_group", width = 14, height = 9)

  natural <- pred[comparison_id == "PRJNA863443_natural_aging_reference"]
  natural[, natural_aging_plot_group := factor(
    natural_aging_plot_group,
    levels = c("young_13w-17.4w", "mid_76w", "old_86w-87w")
  )]
  natural[, natural_aging_panel := "PRJNA863443_Natural aging"]
  p <- ggplot(natural, aes(natural_aging_plot_group, delta_age_week, color = natural_aging_plot_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, color = "grey55") +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.4, alpha = 0.9) +
    facet_grid(model_family ~ natural_aging_panel, scales = "free_x", space = "free_x") +
    labs(
      x = "Natural aging group",
      y = "Raw delta age (weeks)",
      color = "Natural aging group"
    ) +
    theme_prjna()
  save_prjna_plot(p, "08_prjna863443_natural_aging_raw_delta_by_group", width = 10, height = 8)

  p <- ggplot(natural, aes(age_week, predicted_age_week, color = natural_aging_plot_group)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(size = 2.6, alpha = 0.9) +
    facet_wrap(~model_family) +
    labs(
      x = "Chronological age (weeks)",
      y = "Predicted age (weeks)",
      color = "Natural aging group"
    ) +
    theme_prjna()
  save_prjna_plot(p, "09_prjna863443_predicted_vs_age_natural_aging_bins", width = 11, height = 7)

  vcd <- rel[grepl("PRJNA863443_VCD_vs_CTL", comparison_id) & matched_control_allowed == TRUE]
  p <- ggplot(vcd, aes(projection_group, relative_delta_age_week, color = projection_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, color = "grey55") +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.85) +
    facet_grid(model_family ~ comparison_id, scales = "free_x", space = "free_x") +
    labs(x = "Projection group", y = "Relative delta age (weeks)", color = "Projection class") +
    theme_prjna() +
    theme(strip.text.x = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 35, hjust = 1))
  save_prjna_plot(p, "03_prjna863443_vcd_relative_delta", width = 16, height = 8)

  fox <- rel[grepl("PRJNA863443_Foxl2_het_vs_WT", comparison_id) & matched_control_allowed == TRUE]
  p <- ggplot(fox, aes(projection_group, relative_delta_age_week, color = projection_class)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, color = "grey55") +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.85) +
    facet_grid(model_family ~ comparison_id, scales = "free_x", space = "free_x") +
    labs(x = "Projection group", y = "Relative delta age (weeks)", color = "Projection class") +
    theme_prjna()
  save_prjna_plot(p, "04_prjna863443_foxl2_relative_delta", width = 11, height = 8)

  p <- ggplot(concordance, aes(elasticnet, pca_glmnet, color = projection_group)) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(size = 2.5, alpha = 0.9) +
    labs(x = "ElasticNet delta age (weeks)", y = "PCA/glmnet delta age (weeks)", color = "Projection group") +
    theme_prjna()
  save_prjna_plot(p, "05_prjna863443_clock_concordance", width = 9, height = 7)

  p <- ggplot(coverage, aes(projection_group, gene_coverage_fraction, color = model_family)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12, color = "grey55") +
    geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.85) +
    labs(x = "Projection group", y = "Model gene coverage fraction", color = "Model family") +
    theme_prjna()
  save_prjna_plot(p, "06_prjna863443_gene_coverage", width = 14, height = 7)

  contrib_cols <- grep("_contribution_week$", names(pca_scores), value = TRUE)
  if (length(contrib_cols) > 0) {
    heat <- merge(pca_scores, unique(pred[, .(sample_key, projection_group)]), by = "sample_key", all.x = TRUE)
    heat <- melt(
      heat,
      id.vars = c("sample_key", "projection_group", "model_family", "selector"),
      measure.vars = contrib_cols,
      variable.name = "pc",
      value.name = "pc_contribution_week"
    )
    heat[, pc := gsub("_contribution_week$", "", pc)]
    group_heat <- heat[, .(pc_contribution_week = mean(pc_contribution_week, na.rm = TRUE)), by = .(projection_group, pc)]
    p <- ggplot(group_heat, aes(pc, projection_group, fill = pc_contribution_week)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44E52", midpoint = 0) +
      labs(x = "Frozen selected PC", y = "Projection group", fill = "Weeks") +
      theme_prjna()
    save_prjna_plot(p, "07_prjna863443_pca_pc_contribution_heatmap", width = 9, height = 7)
  }
}

message("Projecting PRJNA863443 with frozen v1 deployment artifacts.")
for (path in c(metadata_path, counts_path, log_path)) {
  if (!file.exists(path)) stop("Missing required PRJNA863443 input: ", path, call. = FALSE)
}

artifact_audit <- audit_artifacts(retained_clocks)
meta_raw <- fread(metadata_path, sep = "\t", na.strings = c("", "NA", "NaN"))
if (nrow(meta_raw) != 33 || uniqueN(meta_raw$sample_id) != 33) {
  stop("Expected exactly 33 PRJNA863443 samples; observed rows=", nrow(meta_raw), ", unique sample_id=", uniqueN(meta_raw$sample_id), call. = FALSE)
}
meta <- make_prjna_metadata(meta_raw)

counts <- read_expression_matrix(counts_path)
log_from_counts <- counts_to_logcpm(counts)
col_map <- meta$sample_key
names(col_map) <- meta$sample_id
missing_cols <- setdiff(meta$sample_id, colnames(log_from_counts))
if (length(missing_cols) > 0) stop("PRJNA863443 expression matrix missing samples: ", paste(missing_cols, collapse = ", "), call. = FALSE)
log_from_counts <- log_from_counts[, meta$sample_id, drop = FALSE]
colnames(log_from_counts) <- col_map[colnames(log_from_counts)]
expr_pieces <- list(PRJNA863443 = log_from_counts)

artifact_scales <- unique(vapply(artifact_audit$objects, function(obj) obj$expression_scale %||% NA_character_, character(1)))
scale_audit <- data.table(
  input_count_file = counts_path,
  input_log_file = log_path,
  expression_scale_detected = "raw counts transformed to log1p(CPM)",
  expression_scale_required_by_artifact = paste(artifact_scales, collapse = ";"),
  transformation_used = "log1p(raw_counts / column_sum(raw_counts) * 1e6)",
  n_genes_raw = nrow(counts),
  n_genes_after_filter = nrow(log_from_counts),
  n_samples = ncol(log_from_counts),
  decision = "used raw gene-level counts and transformed to logCPM because deployment artifacts require logCPM",
  warning = fifelse(all(grepl("logCPM", artifact_scales, ignore.case = TRUE)), NA_character_, "At least one artifact did not declare logCPM")
)
fwrite(scale_audit, file.path(dataset_dir, "prjna863443_expression_scale_audit.tsv"), sep = "\t")

reference_meta <- meta[projection_role == "reference_control"]
if (nrow(reference_meta) < 3) stop("Too few PRJNA863443 reference rows for OOD distance calculations.", call. = FALSE)

prediction_list <- list()
pca_score_list <- list()
for (i in seq_len(nrow(retained_clocks))) {
  clock <- retained_clocks[i]
  key <- paste(clock$model_family, clock$selector, sep = "\r")
  obj <- artifact_audit$objects[[key]]
  cal <- extract_calibration(obj)
  message("Projecting PRJNA863443 ", clock$model_family, " / ", clock$selector, ".")
  res <- clock_predict(clock, obj, cal, expr_pieces, meta, reference_meta)
  prediction_list[[key]] <- res$predictions
  if (!is.null(res$pca_scores)) pca_score_list[[key]] <- res$pca_scores
}

pred <- rbindlist(prediction_list, use.names = TRUE, fill = TRUE)
pred <- merge(pred, meta, by = "sample_key", all.x = TRUE)
pred[, delta_age_week := predicted_age_week - age_week]
pred[, used_in_deployment_training_flag := FALSE]
pred[, same_dataset_control_available := TRUE]
pred[, matched_control_available := matched_control_allowed]
pred[, primary_external_projection_include := as_bool(primary_external_projection_include)]
pred[, primary_effect_interpretation_include := as_bool(primary_effect_interpretation_include)]
pred[, matched_control_allowed := as_bool(matched_control_allowed)]

relative <- compute_prjna_relative(pred)
pred_match <- unique(relative[, .(
  sample_key, model_family, selector, matched_control_method, matched_control_source,
  n_matched_controls, relative_control_reliability
)])
pred <- merge(pred, pred_match, by = c("sample_key", "model_family", "selector"), all.x = TRUE)

treatment_cols <- c("treatment", "genotype", "sample_context", "sample_title", "model", "condition", "injection_age", "post_injection_days", "age_group", "natural_aging_plot_group")
pred_cols <- c(
  "sample_key", "dataset_folder", "modality", "age_week", "bio_age_bin",
  "projection_role", "projection_class", "projection_reason",
  "matched_reference_control_flag", "matched_reference_control_source",
  "used_in_deployment_training_flag", "primary_external_projection_include",
  "primary_effect_interpretation_include", "matched_control_allowed",
  "reference_control_group_id", "reference_control_group_label",
  "gse279446_cohort", "gse279446_platform", "projection_dataset_note",
  treatment_cols, "sample_id", "accession", "projection_group", "comparison_id",
  "genotype_class", "cohort", "control_group_label",
  "model_family", "selector", "raw_predicted_age_week", "predicted_age_week", "delta_age_week",
  "gene_coverage_fraction", "n_model_genes", "n_observed_model_genes", "missing_model_genes",
  "projection_gene_coverage_flag", "projection_distance", "projection_distance_percentile", "OOD_flag",
  "matched_control_method", "matched_control_source", "n_matched_controls", "relative_control_reliability"
)
pred_cols <- pred_cols[pred_cols %in% names(pred)]
pred_out <- pred[, ..pred_cols]

rel_cols <- c(
  "sample_key", "dataset_folder", "modality", "age_week", "bio_age_bin",
  "projection_role", "projection_class", "projection_reason",
  "matched_reference_control_flag", "matched_reference_control_source",
  "used_in_deployment_training_flag", "primary_external_projection_include",
  "primary_effect_interpretation_include", "matched_control_allowed",
  "reference_control_group_id", "reference_control_group_label",
  "gse279446_cohort", "gse279446_platform",
  treatment_cols, "sample_id", "accession", "projection_group", "comparison_id",
  "genotype_class", "cohort", "control_group_label", "model_family", "selector",
  "control_expected_predicted_age_week", "control_expected_delta_age_week",
  "relative_delta_age_week", "control_expected_pred_age", "control_expected_delta_age",
  "relative_delta_age", "matched_control_method", "matched_control_source",
  "n_matched_controls", "matched_control_sample_keys", "relative_control_reliability",
  "relative_control_warning", "projection_dataset_note"
)
rel_cols <- rel_cols[rel_cols %in% names(relative)]
relative_out <- relative[, ..rel_cols]

fwrite(meta, file.path(dataset_dir, "prjna863443_projection_metadata.tsv"), sep = "\t")
fwrite(pred_out, file.path(dataset_dir, "prjna863443_projection_predictions.tsv"), sep = "\t")
fwrite(relative_out, file.path(dataset_dir, "prjna863443_relative_to_controls.tsv"), sep = "\t")
coverage <- pred[, .(
  sample_key, dataset_folder, projection_group, model_family, selector,
  n_model_genes, n_observed_model_genes, gene_coverage_fraction,
  projection_gene_coverage_flag, missing_model_genes
)]
fwrite(coverage, file.path(dataset_dir, "prjna863443_gene_coverage.tsv"), sep = "\t")

group_summary <- merge(
  pred[, .(
    n_samples = uniqueN(sample_key),
    mean_age_week = mean(age_week, na.rm = TRUE),
    mean_predicted_age_week = mean(predicted_age_week, na.rm = TRUE),
    mean_delta_age_week = mean(delta_age_week, na.rm = TRUE),
    median_delta_age_week = median(delta_age_week, na.rm = TRUE)
  ), by = .(dataset_folder, projection_class, group = projection_group, comparison_id, model_family, selector)],
  relative[, .(
    mean_relative_delta_age_week = mean(relative_delta_age_week, na.rm = TRUE),
    median_relative_delta_age_week = median(relative_delta_age_week, na.rm = TRUE),
    mean_n_matched_controls = mean(n_matched_controls, na.rm = TRUE)
  ), by = .(dataset_folder, projection_class, group = projection_group, comparison_id, model_family, selector)],
  by = c("dataset_folder", "projection_class", "group", "comparison_id", "model_family", "selector"),
  all.x = TRUE
)
fwrite(group_summary, file.path(dataset_dir, "prjna863443_group_summary.tsv"), sep = "\t")

concordance <- dcast(
  pred[, .(sample_key, dataset_folder, projection_class, projection_group, delta_age_week, model_family)],
  sample_key + dataset_folder + projection_class + projection_group ~ model_family,
  value.var = "delta_age_week"
)
if (all(c("elasticnet", "pca_glmnet") %in% names(concordance))) {
  concordance[, `:=`(
    delta_age_difference_pca_minus_elasticnet = pca_glmnet - elasticnet,
    absolute_delta_age_disagreement = abs(pca_glmnet - elasticnet)
  )]
}
fwrite(concordance, file.path(dataset_dir, "prjna863443_model_concordance.tsv"), sep = "\t")

pca_scores <- rbindlist(pca_score_list, use.names = TRUE, fill = TRUE)
if (nrow(pca_scores) > 0) {
  fwrite(pca_scores, file.path(dataset_dir, "prjna863443_pca_scores.tsv"), sep = "\t")
  contrib_cols <- grep("_contribution_week$", names(pca_scores), value = TRUE)
  if (length(contrib_cols) > 0) {
    contrib <- melt(
      pca_scores,
      id.vars = intersect(c("sample_key", "model_family", "selector"), names(pca_scores)),
      measure.vars = contrib_cols,
      variable.name = "pc",
      value.name = "pc_contribution_week"
    )
    contrib[, pc := gsub("_contribution_week$", "", pc)]
    fwrite(contrib, file.path(dataset_dir, "prjna863443_pca_pc_contributions.tsv"), sep = "\t")
  }
}

warnings <- character()
low_cov <- coverage[gene_coverage_fraction < 0.80]
if (nrow(low_cov) > 0) warnings <- c(warnings, paste0("Low coverage rows: ", nrow(low_cov)))
if (nrow(relative[matched_control_method == "no_matched_control" & is.finite(relative_delta_age_week)]) > 0) {
  warnings <- c(warnings, "Invalid no_matched_control relative delta rows detected.")
}
if (length(warnings) == 0) warnings <- "No PRJNA863443 projection warnings."
writeLines(warnings, file.path(dataset_dir, "prjna863443_projection_warnings.txt"))

summary_lines <- c(
  "# PRJNA863443 external projection summary",
  "",
  "- v1 deployment artifacts were used; no fitting, scaling, PCA, glmnet, or calibration was refit.",
  paste0("- Samples: ", uniqueN(pred$sample_key), "."),
  paste0("- Prediction rows: ", nrow(pred), "."),
  paste0("- Relative-control rows: ", nrow(relative), "."),
  paste0("- Natural ageing samples: ", uniqueN(pred[comparison_id == "PRJNA863443_natural_aging_reference"]$sample_key), " placement-only reference samples."),
  paste0("- VCD samples: ", uniqueN(pred[grepl("PRJNA863443_VCD_vs_CTL", comparison_id)]$sample_key), " across four matched CTL/VCD comparisons."),
  paste0("- Foxl2 samples: ", uniqueN(pred[grepl("PRJNA863443_Foxl2_het_vs_WT", comparison_id)]$sample_key), " across young and mid_age comparisons."),
  paste0("- Minimum gene coverage: ", paste(coverage[, paste0(model_family, "=", round(min(gene_coverage_fraction), 4)), by = model_family]$V1, collapse = "; "), ".")
)
writeLines(summary_lines, file.path(dataset_dir, "prjna863443_projection_summary.md"))

reference_controls <- pred[matched_reference_control_flag == TRUE, .(
  dataset_folder, sample_key, age_week, model_family, selector,
  projection_role, matched_reference_control_source,
  used_in_deployment_training_flag, predicted_age_week, delta_age_week,
  reference_control_group_id, reference_control_group_label
)]
control_pool_candidates <- meta[, .(
  dataset_folder, sample_key, sample_id, accession, modality, age_week, bio_age_bin,
  treatment, sample_context, sample_title, projection_role, projection_class,
  matched_reference_control_flag, matched_reference_control_source,
  reference_control_group_id, reference_control_group_label, control_pool_source,
  control_candidate_status, control_candidate_note
)]
control_matching_audit <- relative[, .(
  dataset_folder, sample_key, model_family, matched_control_method, matched_control_source,
  n_matched_controls, matched_control_sample_keys, relative_control_reliability,
  relative_control_warning
)]
qc_flags <- pred[, .(
  sample_key, dataset_folder, model_family, selector, projection_gene_coverage_flag,
  OOD_flag, projection_distance, projection_distance_percentile,
  primary_external_projection_include, matched_control_allowed
)]

append_replace_dataset("external_projection_predictions_all.tsv", pred_out)
append_replace_dataset("external_projection_predictions.tsv", pred_out[primary_external_projection_include == TRUE])
append_replace_dataset("external_projection_relative_to_controls_all.tsv", relative_out)
append_replace_dataset("external_projection_relative_to_controls.tsv", relative_out[primary_external_projection_include == TRUE])
append_replace_dataset("external_projection_reference_controls.tsv", reference_controls)
append_replace_dataset("external_projection_control_pool_candidates.tsv", control_pool_candidates)
append_replace_dataset("external_projection_control_matching_audit.tsv", control_matching_audit)
append_replace_dataset("external_projection_group_summary.tsv", group_summary)
append_replace_dataset("external_projection_model_concordance.tsv", concordance)
append_replace_dataset("external_projection_gene_coverage.tsv", coverage)
append_replace_dataset("external_projection_qc_flags.tsv", qc_flags)
if (nrow(pca_scores) > 0) {
  sample_dataset_lookup <- unique(pred[, .(sample_key, dataset_folder)])
  append_replace_dataset("pca_projection_scores.tsv", merge(pca_scores, sample_dataset_lookup, by = "sample_key", all.x = TRUE))
  contrib_path <- file.path(dataset_dir, "prjna863443_pca_pc_contributions.tsv")
  if (file.exists(contrib_path)) {
    append_replace_dataset("pca_projection_pc_contributions.tsv", merge(fread(contrib_path, sep = "\t"), sample_dataset_lookup, by = "sample_key", all.x = TRUE))
  }
}
update_dataset_availability(pred, relative)
invalid_global <- update_run_manifest_and_summary()

write_dataset_plots(pred, relative, concordance, coverage, pca_scores)

cat("PRJNA863443 sample count:", uniqueN(pred$sample_key), "\n")
cat("PRJNA863443 prediction rows:", nrow(pred), "\n")
cat("PRJNA863443 gene coverage by model:\n")
print(coverage[, .(
  min_gene_coverage = min(gene_coverage_fraction, na.rm = TRUE),
  median_gene_coverage = median(gene_coverage_fraction, na.rm = TRUE),
  n_low_coverage = sum(gene_coverage_fraction < 0.80, na.rm = TRUE)
), by = model_family])
cat("Global invalid relative-control rows:", nrow(invalid_global), "\n")
message("PRJNA863443 projection appended to: ", projection_dir)
