#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x
}

safe_chr <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_dir <- normalizePath(file.path(project_dir, "..", "data_input", "bulkseq_dataset", "CD1mice"), mustWork = TRUE)
out_dir <- file.path(project_dir, "results", "external_projection_CD1mice")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

warnings_vec <- character()
add_warning <- function(message) warnings_vec <<- unique(c(warnings_vec, message))

datasets <- c("GSE179888", "GSE264327", "GSE299728", "GSE310764", "GSE325001", "CRA027217")

model_table <- data.table(
  model_family = c("elasticnet", "pca_glmnet"),
  selector = c("hvg_baseline", "version_b_age_dataset_penalized"),
  model_path = file.path(
    project_dir, "results", "model_artifacts", "deployment",
    c(
      "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds",
      "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
    )
  )
)

age_to_week <- function(age, unit) {
  unit <- tolower(safe_chr(unit))
  fifelse(
    unit %chin% c("day", "days"), age / 7,
    fifelse(
      unit %chin% c("week", "weeks"), age,
      fifelse(
        unit %chin% c("month", "months"), age * 52.1775 / 12,
        fifelse(unit %chin% c("year", "years"), age * 52.1775, NA_real_)
      )
    )
  )
}

bio_age_bin <- function(age_week) {
  fifelse(
    age_week < 7, "developmental_or_prepubertal",
    fifelse(
      age_week < 16, "young_adult",
      fifelse(age_week < 40, "mature_adult", fifelse(age_week < 65, "advanced_maternal_age", "post_reproductive_ageing"))
    )
  )
}

read_dataset_metadata <- function(dataset) {
  path <- file.path(data_dir, dataset, "metadata", paste0(dataset, "_metadata.csv"))
  metadata <- fread(path, check.names = FALSE)
  required <- c(
    "Sample name", "sample_column", "Chronological age", "Chronological age unit",
    "Breed", "Strain", "Sex", "Tissue", "Treatment", "Phenotype",
    "Genotype", "Sample title", "Age raw", "Age curation confidence"
  )
  missing <- setdiff(required, names(metadata))
  if (length(missing)) stop(dataset, " metadata is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(metadata$sample_column)) stop(dataset, " has duplicated sample_column values.", call. = FALSE)

  out <- data.table(
    sample_key = paste("CD1", dataset, metadata$sample_column, sep = "__"),
    sample_id = safe_chr(metadata[["Sample name"]]),
    sample_column = safe_chr(metadata$sample_column),
    accession = safe_chr(metadata$Accession),
    dataset_folder = dataset,
    modality = "bulk_rnaseq",
    strain = safe_chr(metadata$Strain),
    breed = safe_chr(metadata$Breed),
    sex = tolower(safe_chr(metadata$Sex)),
    tissue = safe_chr(metadata$Tissue),
    age_value = safe_num(metadata[["Chronological age"]]),
    age_unit = safe_chr(metadata[["Chronological age unit"]]),
    age_raw = safe_chr(metadata[["Age raw"]]),
    age_curation_confidence = safe_chr(metadata[["Age curation confidence"]]),
    treatment = safe_chr(metadata$Treatment),
    phenotype = safe_chr(metadata$Phenotype),
    health_state = safe_chr(metadata[["Health state"]]),
    genotype = safe_chr(metadata$Genotype),
    sample_title = safe_chr(metadata[["Sample title"]]),
    expression_logcpm_path = file.path(data_dir, dataset, "input", paste0(dataset, "_gene_logcpm_filtered.csv"))
  )
  out[, age_week := age_to_week(age_value, age_unit)]
  out[, age_month_equivalent := age_week * 12 / 52.1775]
  out[, bio_age_bin := bio_age_bin(age_week)]
  out[, `:=`(
    projection_class = "reference_control",
    projection_role = "reference_control",
    projection_group = treatment,
    comparison_id = paste0(dataset, "_natural_age_reference"),
    reference_control_group_id = paste0(dataset, "_age_", sprintf("%.6f", age_week)),
    projection_reason = "untreated_natural_age_reference",
    matched_reference_control_flag = TRUE,
    primary_external_projection_include = age_week >= 7,
    matched_control_allowed = age_week >= 7
  )]

  if (dataset == "GSE310764") {
    out[, projection_group := fcase(
      phenotype == "young_control", "Young adult reference (2 months)",
      phenotype == "reproductively_old", "Naturally aged, reproductively old (18 months)",
      phenotype == "post_reproductive", "Naturally aged, post-reproductive (24 months)",
      default = phenotype
    )]
  }

  if (dataset == "GSE264327") {
    is_oocyte <- grepl("oocyte", out$sample_id, ignore.case = TRUE)
    is_control <- grepl("Ctrl", out$sample_id, ignore.case = TRUE)
    is_cy <- grepl("_Cy_", out$sample_id, ignore.case = TRUE)
    is_dq <- grepl("_DQ_", out$sample_id, ignore.case = TRUE)
    out[, `:=`(
      projection_class = fifelse(is_control, "reference_control", "intervention"),
      projection_role = fifelse(is_control, "reference_control", "external_treatment"),
      projection_group = fifelse(is_control, "untreated control", fifelse(is_cy, "cyclophosphamide POF", "dasatinib + quercetin rescue")),
      comparison_id = "GSE264327_whole_ovary_DQ_rescue",
      reference_control_group_id = "GSE264327_whole_ovary_controls",
      projection_reason = fifelse(is_control, "same_study_WT_reference_control", fifelse(is_cy, "cyclophosphamide_POI_model", "senolytic_rescue_intervention")),
      matched_reference_control_flag = is_control
    )]
    out[is_oocyte, `:=`(
      projection_class = "other",
      projection_role = "OOD_supplement",
      projection_reason = "isolated_oocyte_composition_OOD_for_whole_ovary_clock",
      comparison_id = "GSE264327_isolated_oocyte_OOD",
      reference_control_group_id = "no_matched_control",
      matched_reference_control_flag = FALSE,
      primary_external_projection_include = FALSE,
      matched_control_allowed = FALSE
    )]
  }

  if (dataset == "GSE299728") {
    is_control <- grepl("^NC_", out$sample_id)
    group <- fifelse(
      is_control, "normal control",
      fifelse(grepl("^Cy_", out$sample_id), "cyclophosphamide", fifelse(grepl("^UA_", out$sample_id), "urolithin A", "cyclophosphamide + urolithin A"))
    )
    out[, `:=`(
      projection_class = "developmental_OOD",
      projection_role = "OOD_supplement",
      projection_group = group,
      comparison_id = "GSE299728_neonatal_chemotherapy_OOD",
      reference_control_group_id = "no_matched_control",
      projection_reason = "six_day_postnatal_outside_adult_clock_age_frame",
      matched_reference_control_flag = FALSE,
      primary_external_projection_include = FALSE,
      matched_control_allowed = FALSE
    )]
  }

  if (dataset == "GSE325001") {
    is_a <- grepl("^A[0-9]+$", out$sample_id)
    is_g <- grepl("^G[0-9]+$", out$sample_id)
    is_rya <- grepl("^RYA", out$sample_id)
    is_raa <- grepl("^RAA", out$sample_id)
    is_ryy <- grepl("^RYY", out$sample_id)
    out[, `:=`(
      projection_class = fifelse(is_a | is_raa | is_ryy, "reference_control", fifelse(is_rya, "microbiome", "intervention")),
      projection_role = fifelse(is_a | is_raa | is_ryy, "reference_control", "external_treatment"),
      projection_group = fifelse(is_a, "aged untreated control", fifelse(is_g, "aged + glutamic acid", fifelse(is_rya, "aged + young-donor FMT", fifelse(is_raa, "aged + aged-donor FMT control", "young + young-donor FMT control")))),
      comparison_id = fifelse(is_a | is_g, "GSE325001_glutamic_acid_vs_aged", fifelse(is_rya | is_raa, "GSE325001_young_FMT_vs_aged_FMT", "GSE325001_young_FMT_reference")),
      reference_control_group_id = fifelse(is_a | is_g, "GSE325001_aged_untreated_controls", fifelse(is_rya | is_raa, "GSE325001_aged_aged_FMT_controls", "GSE325001_young_young_FMT_controls")),
      projection_reason = fifelse(is_a, "same_study_aged_reference_control", fifelse(is_g, "glutamic_acid_intervention", fifelse(is_rya, "young_donor_microbiome_rescue", fifelse(is_raa, "age_matched_aged_donor_microbiome_control", "young_recipient_young_donor_microbiome_reference")))),
      matched_reference_control_flag = is_a | is_raa | is_ryy
    )]
  }

  if (dataset == "CRA027217") {
    is_control <- grepl("_Con", out$sample_id)
    out[, `:=`(
      projection_class = fifelse(is_control, "reference_control", "intervention"),
      projection_role = fifelse(is_control, "reference_control", "external_treatment"),
      projection_group = fifelse(is_control, "F3 vehicle control", "F3 transgenerational propylparaben DOR"),
      comparison_id = "CRA027217_F3_PrP_DOR_vs_control",
      reference_control_group_id = "CRA027217_F3_vehicle_controls",
      projection_reason = fifelse(is_control, "same_study_vehicle_reference_control", "transgenerational_prenatal_propylparaben_exposure"),
      matched_reference_control_flag = is_control
    )]
  }

  out[, primary_effect_interpretation_include := primary_external_projection_include & projection_role != "OOD_supplement"]
  out
}

metadata <- rbindlist(lapply(datasets, read_dataset_metadata), use.names = TRUE, fill = TRUE)
if (nrow(metadata) != 75L) stop("Expected 75 CD1 samples; found ", nrow(metadata), call. = FALSE)
if (any(!is.finite(metadata$age_week))) stop("Non-finite age_week found.", call. = FALSE)
if (any(metadata$breed != "ICR/CD1")) stop("Non-ICR/CD1 metadata row found.", call. = FALSE)
if (any(metadata$sex != "female")) stop("Non-female metadata row found.", call. = FALSE)
if (anyDuplicated(metadata$sample_key)) stop("Duplicated sample_key found.", call. = FALSE)

read_logcpm <- function(dataset, metadata_subset) {
  path <- unique(metadata_subset$expression_logcpm_path)
  if (length(path) != 1L || !file.exists(path)) stop("Missing logCPM matrix for ", dataset, call. = FALSE)
  input <- fread(path, check.names = FALSE)
  gene_col <- names(input)[1]
  genes <- safe_chr(input[[gene_col]])
  sample_columns <- metadata_subset$sample_column
  missing_samples <- setdiff(sample_columns, names(input))
  if (length(missing_samples)) stop(dataset, " logCPM matrix lacks: ", paste(missing_samples, collapse = ", "), call. = FALSE)
  matrix_values <- as.matrix(input[, ..sample_columns])
  mode(matrix_values) <- "numeric"
  rownames(matrix_values) <- genes
  if (any(!is.finite(matrix_values))) stop(dataset, " logCPM matrix contains non-finite values.", call. = FALSE)
  if (anyDuplicated(rownames(matrix_values))) stop(dataset, " logCPM matrix has duplicated gene identifiers.", call. = FALSE)
  matrix_values
}

expression_by_dataset <- setNames(lapply(datasets, function(dataset) {
  read_logcpm(dataset, metadata[dataset_folder == dataset])
}), datasets)

extract_calibration <- function(object) {
  calibration <- object$calibration
  intercept <- calibration$intercept %||% calibration$calibration_intercept %||% calibration$a
  slope <- calibration$slope %||% calibration$calibration_slope %||% calibration$b
  if (length(intercept) != 1L || length(slope) != 1L || !is.finite(intercept) || !is.finite(slope)) {
    stop("Invalid deployment calibration.", call. = FALSE)
  }
  list(intercept = as.numeric(intercept), slope = as.numeric(slope))
}

validate_artifact <- function(object, clock) {
  required <- c("selected_genes", "gene_center", "gene_scale", "model", "calibration")
  missing <- required[vapply(required, function(name) is.null(object[[name]]), logical(1))]
  if (length(missing)) stop(clock$model_family, " artifact lacks: ", paste(missing, collapse = ", "), call. = FALSE)
  stopifnot(
    identical(object$training_scope, "integrated_bulk_pseudobulk"),
    identical(object$artifact_scope, "deployment_all_control_like"),
    identical(object$model_family, clock$model_family),
    identical(object$selector, clock$selector),
    identical(object$expression_scale, "logCPM")
  )
  if (clock$model_family == "pca_glmnet") {
    stopifnot(identical(as.character(object$model$selected_pc_names), paste0("PC", 1:7)))
  }
  TRUE
}

build_scaled_matrix <- function(logcpm, metadata_subset, object) {
  genes <- as.character(object$selected_genes)
  raw <- matrix(
    NA_real_, nrow = nrow(metadata_subset), ncol = length(genes),
    dimnames = list(metadata_subset$sample_key, genes)
  )
  observed_genes <- intersect(genes, rownames(logcpm))
  raw[, observed_genes] <- t(logcpm[observed_genes, metadata_subset$sample_column, drop = FALSE])
  missing <- is.na(raw)
  for (gene in genes) raw[missing[, gene], gene] <- object$gene_center[[gene]]
  scale_use <- object$gene_scale[genes]
  scale_use[!is.finite(scale_use) | scale_use == 0] <- 1
  scaled <- sweep(raw, 2, object$gene_center[genes], "-")
  scaled <- sweep(scaled, 2, scale_use, "/")
  scaled[!is.finite(scaled)] <- 0
  coverage <- data.table(
    sample_key = metadata_subset$sample_key,
    n_model_genes = length(genes),
    n_observed_model_genes = rowSums(!missing),
    gene_coverage_fraction = rowSums(!missing) / length(genes),
    missing_model_genes = apply(missing, 1, function(flag) paste(genes[flag], collapse = ";"))
  )
  coverage[, low_coverage_flag := gene_coverage_fraction < 0.80]
  list(x = scaled, coverage = coverage)
}

robust_mahalanobis <- function(x, reference) {
  center <- colMeans(reference, na.rm = TRUE)
  covariance <- as.matrix(stats::cov(reference, use = "pairwise.complete.obs"))
  covariance[!is.finite(covariance)] <- 0
  diag(covariance) <- diag(covariance) + 1e-6
  as.numeric(stats::mahalanobis(x, center = center, cov = covariance))
}

distance_percentile <- function(distance, reference_distance) {
  vapply(distance, function(value) mean(reference_distance <= value, na.rm = TRUE), numeric(1))
}

projection_results <- list()
pca_scores_list <- list()
pca_contribution_list <- list()
artifact_audit <- list()
expression_audit <- list()

for (model_index in seq_len(nrow(model_table))) {
  clock <- model_table[model_index]
  if (!file.exists(clock$model_path)) stop("Missing deployment artifact: ", clock$model_path, call. = FALSE)
  object <- readRDS(clock$model_path)
  validate_artifact(object, clock)
  calibration <- extract_calibration(object)

  artifact_audit[[model_index]] <- data.table(
    model_family = clock$model_family,
    selector = clock$selector,
    artifact_path = clock$model_path,
    training_scope = object$training_scope,
    artifact_scope = object$artifact_scope,
    expression_scale = object$expression_scale,
    n_selected_genes = length(object$selected_genes),
    n_training_samples = length(object$training_sample_keys),
    calibration_intercept = calibration$intercept,
    calibration_slope = calibration$slope,
    model_retrained = FALSE
  )

  for (dataset in datasets) {
    metadata_subset <- metadata[dataset_folder == dataset]
    scaled <- build_scaled_matrix(expression_by_dataset[[dataset]], metadata_subset, object)

    if (clock$model_family == "elasticnet") {
      raw_prediction <- as.numeric(predict(object$model$glmnet_fit, newx = scaled$x))
      predicted <- calibration$intercept + calibration$slope * raw_prediction
      result <- cbind(
        metadata_subset,
        data.table(
          model_family = clock$model_family,
          selector = clock$selector,
          raw_predicted_age_week = raw_prediction,
          predicted_age_week = predicted,
          delta_age_week = predicted - metadata_subset$age_week
        )
      )
    } else {
      pc_names <- as.character(object$model$selected_pc_names)
      pc_all <- predict(object$model$pca_fit, newdata = scaled$x)
      pc_matrix <- pc_all[, pc_names, drop = FALSE]
      raw_prediction <- as.numeric(predict(object$model$glmnet_fit, newx = pc_matrix))
      predicted <- calibration$intercept + calibration$slope * raw_prediction
      training_pc <- object$model$pca_fit$x[, pc_names, drop = FALSE]
      distance <- robust_mahalanobis(pc_matrix, training_pc)
      training_distance <- robust_mahalanobis(training_pc, training_pc)
      percentile <- distance_percentile(distance, training_distance)
      result <- cbind(
        metadata_subset,
        data.table(
          model_family = clock$model_family,
          selector = clock$selector,
          raw_predicted_age_week = raw_prediction,
          predicted_age_week = predicted,
          delta_age_week = predicted - metadata_subset$age_week,
          pca_projection_distance = distance,
          pca_projection_distance_percentile = percentile,
          OOD_flag = fifelse(percentile > 0.99, "OOD_gt_99th", fifelse(percentile > 0.95, "OOD_gt_95th", "in_reference_range"))
        )
      )

      score_table <- as.data.table(pc_matrix)
      score_table[, `:=`(
        sample_key = metadata_subset$sample_key,
        dataset_folder = dataset,
        sample_id = metadata_subset$sample_id,
        projection_class = metadata_subset$projection_class,
        projection_group = metadata_subset$projection_group,
        age_week = metadata_subset$age_week,
        predicted_age_week = predicted,
        delta_age_week = predicted - metadata_subset$age_week,
        pca_projection_distance = distance,
        pca_projection_distance_percentile = percentile,
        OOD_flag = result$OOD_flag
      )]
      setcolorder(score_table, c("sample_key", "dataset_folder", "sample_id", "projection_class", "projection_group", "age_week", pc_names))
      pca_scores_list[[length(pca_scores_list) + 1L]] <- score_table

      coefficients <- as.data.table(object$model$pc_coefficient_table)
      contribution <- melt(
        score_table[, c("sample_key", "dataset_folder", "sample_id", "projection_class", "projection_group", "age_week", pc_names), with = FALSE],
        id.vars = c("sample_key", "dataset_folder", "sample_id", "projection_class", "projection_group", "age_week"),
        variable.name = "pc", value.name = "pc_score"
      )
      contribution <- merge(contribution, coefficients, by = "pc", all.x = TRUE)
      contribution[, calibrated_PC_contribution_week := pc_score * calibrated_pc_coefficient]
      pca_contribution_list[[length(pca_contribution_list) + 1L]] <- contribution
    }

    result <- merge(result, scaled$coverage, by = "sample_key", all.x = TRUE, sort = FALSE)
    projection_results[[length(projection_results) + 1L]] <- result
    expression_audit[[length(expression_audit) + 1L]] <- data.table(
      dataset_folder = dataset,
      model_family = clock$model_family,
      selector = clock$selector,
      expression_path = unique(metadata_subset$expression_logcpm_path),
      n_expression_genes = nrow(expression_by_dataset[[dataset]]),
      n_samples = nrow(metadata_subset),
      n_model_genes = unique(scaled$coverage$n_model_genes),
      n_observed_model_genes = unique(scaled$coverage$n_observed_model_genes),
      gene_coverage_fraction = unique(scaled$coverage$gene_coverage_fraction),
      low_coverage_flag = unique(scaled$coverage$low_coverage_flag)
    )
  }
}

predictions <- rbindlist(projection_results, use.names = TRUE, fill = TRUE)
pca_scores <- rbindlist(pca_scores_list, use.names = TRUE, fill = TRUE)
pca_contributions <- rbindlist(pca_contribution_list, use.names = TRUE, fill = TRUE)
artifact_audit <- rbindlist(artifact_audit)
expression_audit <- rbindlist(expression_audit)

pca_qc <- unique(predictions[model_family == "pca_glmnet", .(
  sample_key, pca_projection_distance, pca_projection_distance_percentile, OOD_flag
)])
predictions <- merge(predictions[, !c("pca_projection_distance", "pca_projection_distance_percentile", "OOD_flag"), with = FALSE], pca_qc, by = "sample_key", all.x = TRUE)

reference_controls <- predictions[
  matched_reference_control_flag == TRUE & primary_external_projection_include == TRUE,
  .(
    reference_delta_age_week = mean(delta_age_week, na.rm = TRUE),
    n_matched_controls = uniqueN(sample_key),
    matched_control_sample_keys = paste(sort(unique(sample_key)), collapse = ";")
  ),
  by = .(model_family, selector, reference_control_group_id)
]

relative <- merge(
  predictions,
  reference_controls,
  by = c("model_family", "selector", "reference_control_group_id"),
  all.x = TRUE,
  sort = FALSE
)
relative[, relative_delta_age_week := delta_age_week - reference_delta_age_week]
relative[matched_control_allowed == FALSE | primary_external_projection_include == FALSE, `:=`(
  relative_delta_age_week = NA_real_,
  reference_delta_age_week = NA_real_,
  n_matched_controls = 0L,
  matched_control_sample_keys = NA_character_
)]
relative[, matched_control_method := fifelse(
  is.finite(relative_delta_age_week), reference_control_group_id, "no_matched_control"
)]
relative[, relative_control_reliability := fifelse(
  matched_control_method == "no_matched_control", "unavailable_OOD_or_unmatched",
  fifelse(n_matched_controls >= 3, "same_study_control_n_ge_3", "same_study_control_n_lt_3")
)]

predictions <- relative[, setdiff(names(predictions), character()), with = FALSE]
predictions[, relative_delta_age_week := relative$relative_delta_age_week]
predictions[, matched_control_method := relative$matched_control_method]
predictions[, n_matched_controls := relative$n_matched_controls]
predictions[, relative_control_reliability := relative$relative_control_reliability]

group_summary <- predictions[, .(
  n_samples = uniqueN(sample_key),
  mean_age_week = mean(age_week),
  mean_predicted_age_week = mean(predicted_age_week),
  sd_predicted_age_week = sd(predicted_age_week),
  mean_delta_age_week = mean(delta_age_week),
  sd_delta_age_week = sd(delta_age_week),
  mean_relative_delta_age_week = mean(relative_delta_age_week, na.rm = TRUE),
  sd_relative_delta_age_week = sd(relative_delta_age_week, na.rm = TRUE),
  mean_gene_coverage_fraction = mean(gene_coverage_fraction),
  n_OOD_gt_95th = sum(OOD_flag != "in_reference_range", na.rm = TRUE)
), by = .(
  dataset_folder, projection_class, projection_role, projection_group,
  comparison_id, model_family, selector, primary_external_projection_include
)]
for (column in c("mean_relative_delta_age_week", "sd_relative_delta_age_week")) {
  set(group_summary, which(is.nan(group_summary[[column]])), column, NA_real_)
}

concordance <- dcast(
  predictions[, .(sample_key, dataset_folder, projection_class, projection_group, age_week, model_family, predicted_age_week, delta_age_week, relative_delta_age_week)],
  sample_key + dataset_folder + projection_class + projection_group + age_week ~ model_family,
  value.var = c("predicted_age_week", "delta_age_week", "relative_delta_age_week")
)
concordance[, `:=`(
  absolute_delta_disagreement = abs(delta_age_week_pca_glmnet - delta_age_week_elasticnet),
  pca_minus_elasticnet_delta = delta_age_week_pca_glmnet - delta_age_week_elasticnet,
  same_delta_direction_flag = sign(delta_age_week_pca_glmnet) == sign(delta_age_week_elasticnet)
)]

contrast_map <- unique(metadata[
  primary_external_projection_include == TRUE & projection_role == "external_treatment",
  .(dataset_folder, comparison_id, projection_group, reference_control_group_id)
])
condition_tests <- rbindlist(lapply(seq_len(nrow(contrast_map)), function(index) {
  contrast <- contrast_map[index]
  rbindlist(lapply(unique(predictions$model_family), function(model_name) {
    condition <- predictions[
      dataset_folder == contrast$dataset_folder & model_family == model_name &
        comparison_id == contrast$comparison_id & projection_group == contrast$projection_group,
      delta_age_week
    ]
    controls <- predictions[
      dataset_folder == contrast$dataset_folder & model_family == model_name &
        matched_reference_control_flag == TRUE & reference_control_group_id == contrast$reference_control_group_id,
      delta_age_week
    ]
    if (length(condition) < 2L || length(controls) < 2L) return(NULL)
    test <- t.test(condition, controls)
    data.table(
      dataset_folder = contrast$dataset_folder,
      comparison_id = contrast$comparison_id,
      condition_group = contrast$projection_group,
      reference_control_group_id = contrast$reference_control_group_id,
      model_family = model_name,
      test = "Welch two-sample t-test",
      n_condition = length(condition),
      n_control = length(controls),
      mean_condition_delta_age_week = mean(condition),
      mean_control_delta_age_week = mean(controls),
      mean_difference_week = mean(condition) - mean(controls),
      t_statistic = unname(test$statistic),
      degrees_freedom = unname(test$parameter),
      p_value = test$p.value,
      confidence_low_week = test$conf.int[1],
      confidence_high_week = test$conf.int[2]
    )
  }), fill = TRUE)
}), fill = TRUE)
if (nrow(condition_tests)) condition_tests[, p_adjust_BH := p.adjust(p_value, method = "BH"), by = model_family]

dataset_availability <- metadata[, .(
  n_samples_all = .N,
  n_samples_primary = sum(primary_external_projection_include),
  n_reference_controls = sum(matched_reference_control_flag & primary_external_projection_include),
  n_external_treatments = sum(projection_role == "external_treatment" & primary_external_projection_include),
  n_OOD_supplement = sum(projection_role == "OOD_supplement"),
  projection_classes = paste(sort(unique(projection_class)), collapse = ";"),
  control_group_ids = paste(sort(unique(reference_control_group_id)), collapse = ";")
), by = dataset_folder]

class_dictionary <- data.table(
  projection_class = c("reference_control", "intervention", "microbiome", "developmental_OOD", "other"),
  definition = c(
    "Same-study untreated, vehicle, natural-age, or design-matched control",
    "Drug, chemotherapy, rescue compound, or inherited exposure intervention",
    "FMT or donor-microbiome manipulation",
    "Chronological age below the adult clock frame (<7 weeks)",
    "Composition or design OOD not represented by another established class"
  ),
  primary_role = c("reference_control", "external_treatment", "external_treatment", "OOD_supplement", "OOD_supplement")
)

validation <- data.table(
  check = c(
    "metadata_sample_count_75", "prediction_rows_150", "unique_sample_keys",
    "all_six_datasets_present", "all_female_ICR_CD1", "all_expression_matrices_aligned",
    "all_model_artifacts_frozen", "no_PC8", "all_GSE299728_developmental_OOD",
    "all_GSE264327_oocytes_OOD", "primary_rows_have_controls",
    "relative_delta_missing_for_OOD", "finite_predictions", "gene_coverage_at_least_80_percent"
  ),
  pass = c(
    nrow(metadata) == 75L,
    nrow(predictions) == 150L,
    !anyDuplicated(metadata$sample_key),
    setequal(unique(metadata$dataset_folder), datasets),
    all(metadata$sex == "female" & metadata$breed == "ICR/CD1"),
    nrow(expression_audit) == 12L,
    all(artifact_audit$model_retrained == FALSE),
    !"PC8" %in% unique(pca_contributions$pc),
    all(metadata[dataset_folder == "GSE299728", projection_class == "developmental_OOD" & projection_role == "OOD_supplement"]),
    all(metadata[dataset_folder == "GSE264327" & grepl("oocyte", sample_id, ignore.case = TRUE), projection_role == "OOD_supplement"]),
    all(predictions[primary_external_projection_include == TRUE, n_matched_controls > 0]),
    all(!is.finite(relative[primary_external_projection_include == FALSE, relative_delta_age_week])),
    all(is.finite(predictions$predicted_age_week)),
    all(predictions$gene_coverage_fraction >= 0.80)
  )
)
validation[, status := fifelse(pass, "PASS", "FAIL")]

fwrite(metadata, file.path(out_dir, "external_projection_CD1mice_metadata.tsv"), sep = "\t", na = "NA")
fwrite(predictions, file.path(out_dir, "external_projection_predictions_all.tsv"), sep = "\t", na = "NA")
fwrite(predictions[primary_external_projection_include == TRUE], file.path(out_dir, "external_projection_predictions.tsv"), sep = "\t", na = "NA")
fwrite(relative, file.path(out_dir, "external_projection_relative_to_controls_all.tsv"), sep = "\t", na = "NA")
fwrite(relative[primary_external_projection_include == TRUE], file.path(out_dir, "external_projection_relative_to_controls.tsv"), sep = "\t", na = "NA")
fwrite(reference_controls, file.path(out_dir, "external_projection_reference_controls.tsv"), sep = "\t", na = "NA")
fwrite(group_summary, file.path(out_dir, "external_projection_group_summary.tsv"), sep = "\t", na = "NA")
fwrite(concordance, file.path(out_dir, "external_projection_model_concordance.tsv"), sep = "\t", na = "NA")
fwrite(condition_tests, file.path(out_dir, "external_projection_condition_vs_control_tests.tsv"), sep = "\t", na = "NA")
fwrite(pca_scores, file.path(out_dir, "pca_projection_scores.tsv"), sep = "\t", na = "NA")
fwrite(pca_contributions, file.path(out_dir, "pca_projection_pc_contributions.tsv"), sep = "\t", na = "NA")
fwrite(expression_audit, file.path(out_dir, "external_projection_gene_coverage.tsv"), sep = "\t", na = "NA")
fwrite(artifact_audit, file.path(out_dir, "projection_artifact_audit.tsv"), sep = "\t", na = "NA")
fwrite(dataset_availability, file.path(out_dir, "external_projection_dataset_control_availability.tsv"), sep = "\t", na = "NA")
fwrite(class_dictionary, file.path(out_dir, "external_projection_class_dictionary.tsv"), sep = "\t", na = "NA")
fwrite(validation, file.path(out_dir, "external_projection_validation_checks.tsv"), sep = "\t", na = "NA")

theme_projection <- function() {
  theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(plot, name, width = 11, height = 7) {
  ggsave(file.path(plot_dir, paste0(name, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(plot_dir, paste0(name, ".png")), plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
}

primary <- predictions[primary_external_projection_include == TRUE]
plot <- ggplot(primary, aes(age_week, predicted_age_week, color = projection_class, shape = projection_role)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(size = 2.5, alpha = 0.85) +
  facet_wrap(~model_family) +
  labs(title = "CD1 external projection", subtitle = "Primary adult whole-ovary projection set", x = "Chronological age (weeks)", y = "Predicted age (weeks)", color = "Projection class", shape = "Projection role") +
  theme_projection()
save_plot(plot, "01_external_projection_predicted_vs_chronological_age")

plot <- ggplot(predictions, aes(dataset_folder, delta_age_week, color = projection_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.16, size = 1.5, alpha = 0.75) +
  facet_wrap(~model_family, scales = "free_y") +
  labs(title = "Raw delta age across all CD1 projections", subtitle = "Includes developmental and isolated-oocyte OOD supplements", x = NULL, y = "Delta age (weeks)", color = "Projection class") +
  theme_projection()
save_plot(plot, "02_raw_delta_age_by_dataset_including_OOD", width = 12, height = 7)

plot <- ggplot(primary, aes(projection_group, delta_age_week, color = projection_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  facet_grid(model_family ~ dataset_folder, scales = "free_x", space = "free_x") +
  labs(title = "Primary CD1 raw delta age by study group", x = NULL, y = "Delta age (weeks)", color = "Projection class") +
  theme_projection()
save_plot(plot, "03_raw_delta_age_by_dataset_treatment", width = 18, height = 8)

relative_primary <- relative[primary_external_projection_include == TRUE & is.finite(relative_delta_age_week)]
plot <- ggplot(relative_primary, aes(projection_group, relative_delta_age_week, color = projection_class)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.8) +
  facet_grid(model_family ~ dataset_folder, scales = "free_x", space = "free_x") +
  labs(title = "CD1 delta age relative to same-study controls", x = NULL, y = "Relative delta age (weeks)", color = "Projection class") +
  theme_projection()
save_plot(plot, "04_relative_delta_age_matched_controls", width = 18, height = 8)

plot <- ggplot(concordance, aes(delta_age_week_elasticnet, delta_age_week_pca_glmnet, color = projection_class)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey55") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.2, alpha = 0.8) +
  labs(title = "ElasticNet versus PCA/glmnet delta-age concordance", x = "ElasticNet delta age (weeks)", y = "PCA/glmnet delta age (weeks)", color = "Projection class") +
  theme_projection()
save_plot(plot, "05_clock_concordance_delta_age", width = 9, height = 7)

plot <- ggplot(predictions, aes(gene_coverage_fraction, pca_projection_distance, color = OOD_flag, shape = projection_role)) +
  geom_vline(xintercept = 0.80, linetype = "dashed", color = "grey45") +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~dataset_folder, scales = "free_y") +
  labs(title = "CD1 projection coverage and frozen-PCA OOD distance", x = "Model gene coverage fraction", y = "PCA Mahalanobis distance", color = "OOD flag", shape = "Projection role") +
  theme_projection()
save_plot(plot, "06_projection_qc_coverage_ood", width = 13, height = 8)

plot <- ggplot(pca_scores, aes(PC1, PC2, color = projection_class, shape = dataset_folder)) +
  geom_point(size = 2.4, alpha = 0.85) +
  labs(title = "CD1 samples in the frozen PCA space", x = "Frozen PC1", y = "Frozen PC2", color = "Projection class", shape = "Dataset") +
  theme_projection()
save_plot(plot, "07_pca_pc1_pc2_projection_map", width = 10, height = 7)

heat <- copy(pca_contributions)
heat[, sample_key := factor(sample_key, levels = rev(unique(metadata$sample_key)))]
heat[, pc := factor(pc, levels = paste0("PC", 1:7))]
plot <- ggplot(heat, aes(pc, sample_key, fill = calibrated_PC_contribution_week)) +
  geom_tile(color = "white", linewidth = 0.1) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b8322a", midpoint = 0) +
  facet_grid(dataset_folder ~ ., scales = "free_y", space = "free_y") +
  labs(title = "PCA/glmnet calibrated PC contributions", x = NULL, y = NULL, fill = "Contribution\n(weeks)") +
  theme_projection() + theme(axis.text.y = element_text(size = 5))
save_plot(plot, "08_pca_selected_pc_contribution_heatmap", width = 10, height = 15)

flow <- melt(
  dataset_availability,
  id.vars = "dataset_folder",
  measure.vars = c("n_samples_primary", "n_OOD_supplement"),
  variable.name = "projection_set", value.name = "n_samples"
)
plot <- ggplot(flow, aes(dataset_folder, n_samples, fill = projection_set)) +
  geom_col(position = "stack") +
  labs(title = "CD1 external-projection cohort flow", x = NULL, y = "Samples", fill = "Projection set") +
  theme_projection()
save_plot(plot, "09_projection_cohort_flow_dataset", width = 10, height = 6)

overall_correlation <- cor(concordance$delta_age_week_elasticnet, concordance$delta_age_week_pca_glmnet, use = "complete.obs")
mean_disagreement <- mean(concordance$absolute_delta_disagreement, na.rm = TRUE)
summary_lines <- c(
  "# CD1 Mice External Projection Summary",
  "",
  "## Method",
  "",
  "- Frozen C57-compatible integrated bulk+pseudobulk deployment clocks were applied without retraining.",
  "- Input matrices are the prepared per-dataset `log(1 + CPM)` matrices used by the CD1 preprocessing workflow.",
  "- Missing model genes were imputed at the frozen training center and transformed with frozen gene center/scale parameters.",
  "- ElasticNet uses the frozen HVG baseline artifact; PCA/glmnet uses the frozen age/dataset-penalized 2,000-gene artifact and PC1-PC7 only.",
  "- Predictions and delta ages are in weeks. Relative delta age uses explicitly assigned same-study controls.",
  "",
  "## Projection Classes",
  "",
  "- `reference_control`: untreated, vehicle, natural-age, or design-matched control.",
  "- `intervention`: chemotherapy, rescue compound, glutamic acid, or inherited propylparaben exposure.",
  "- `microbiome`: young-donor FMT manipulation.",
  "- `developmental_OOD`: six-day postnatal GSE299728 samples, retained only in all-projection supplements.",
  "- `other`: isolated-oocyte GSE264327 profiles, retained as composition-OOD supplements.",
  "",
  "## Cohort",
  "",
  paste0("- All samples: ", uniqueN(predictions$sample_key), "."),
  paste0("- Primary adult whole-ovary samples: ", uniqueN(primary$sample_key), "."),
  paste0("- OOD supplements: ", uniqueN(predictions[primary_external_projection_include == FALSE]$sample_key), "."),
  "",
  "## QC",
  "",
  paste0("- Validation checks passed: ", sum(validation$pass), "/", nrow(validation), "."),
  paste0("- Minimum model-gene coverage: ", round(100 * min(predictions$gene_coverage_fraction), 1), "%."),
  paste0("- ElasticNet/PCA delta-age correlation: ", round(overall_correlation, 3), "."),
  paste0("- Mean absolute model disagreement: ", round(mean_disagreement, 3), " weeks."),
  "",
  "## Interpretation Boundary",
  "",
  "These are cross-strain projections of ICR/CD1 bulk RNA-seq into C57-compatible clocks. They position samples on trained transcriptomic ageing axes but do not establish rejuvenation, accelerated ageing, mechanism, or causality. OOD supplements must not be used for primary biological effect claims."
)
writeLines(summary_lines, file.path(out_dir, "external_projection_summary.md"))

manifest <- data.table(
  item = c("script", "output_directory", "models", "model_retrained", "n_samples_all", "n_samples_primary", "n_prediction_rows", "validation_status"),
  value = c(
    normalizePath(script_file), normalizePath(out_dir),
    paste(model_table$model_family, model_table$selector, sep = "/", collapse = ";"),
    "FALSE", uniqueN(predictions$sample_key), uniqueN(primary$sample_key), nrow(predictions),
    ifelse(all(validation$pass), "PASS", "FAIL")
  )
)
fwrite(manifest, file.path(out_dir, "external_projection_run_manifest.tsv"), sep = "\t")

if (any(predictions$OOD_flag != "in_reference_range")) {
  add_warning(paste0(uniqueN(predictions[OOD_flag != "in_reference_range"]$sample_key), " samples exceed the 95th percentile of the frozen PCA training-distance distribution."))
}
if (any(predictions$age_week > 83)) {
  add_warning(paste0(uniqueN(predictions[age_week > 83]$sample_key), " samples exceed the deployment calibration age maximum of 83 weeks; treat extrapolated ages cautiously."))
}
if (length(warnings_vec) == 0L) warnings_vec <- "No projection warnings."
writeLines(warnings_vec, file.path(out_dir, "external_projection_warnings.txt"))

if (!all(validation$pass)) {
  stop("CD1 external projection completed with failed validation checks: ", paste(validation[pass == FALSE]$check, collapse = ", "), call. = FALSE)
}

message("CD1 samples projected: ", uniqueN(predictions$sample_key))
message("Primary samples: ", uniqueN(primary$sample_key))
message("Prediction rows: ", nrow(predictions))
message("Validation checks: ", sum(validation$pass), "/", nrow(validation), " PASS")
message("Output directory: ", out_dir)
