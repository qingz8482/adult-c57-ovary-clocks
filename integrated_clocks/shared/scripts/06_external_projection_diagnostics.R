#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

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
  stop("Refusing to run v1 external projection diagnostics against a v2 results root: ", results_dir, call. = FALSE)
}
artifact_dir <- file.path(results_dir, "model_artifacts")
projection_dir <- file.path(results_dir, "external_projection")
simplified_dir <- file.path(projection_dir, "simplified_classes")
out_dir <- file.path(projection_dir, "diagnostics")
root_plot_dir <- file.path(projection_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(root_plot_dir, recursive = TRUE, showWarnings = FALSE)

warnings_log <- character()
tables_written <- character()
plot_prefixes_written <- character()

warn <- function(...) {
  msg <- paste0(...)
  warnings_log <<- unique(c(warnings_log, msg))
  message("WARNING: ", msg)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

read_required <- function(path) {
  if (!file.exists(path)) stop("Required input missing: ", path, call. = FALSE)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

read_optional <- function(path) {
  if (!file.exists(path)) {
    warn("Optional input not found: ", basename(path))
    return(NULL)
  }
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

read_optional_quiet <- function(path) {
  if (!file.exists(path)) return(NULL)
  fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
}

prefer_simplified <- function(original_filename, simplified_filename) {
  simplified_path <- file.path(simplified_dir, simplified_filename)
  if (file.exists(simplified_path)) return(simplified_path)
  file.path(projection_dir, original_filename)
}

ensure_cols <- function(dt, cols, value = NA_character_) {
  for (cc in cols) {
    if (!cc %in% names(dt)) dt[, (cc) := value]
  }
  dt
}

coerce_num <- function(dt, cols) {
  for (cc in intersect(cols, names(dt))) dt[, (cc) := as.numeric(get(cc))]
  dt
}

first_present <- function(candidates, names_vec) {
  hit <- intersect(candidates, names_vec)
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "missing", x)
}

label_projection_class <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & x == "stress"] <- "stress_model"
  x
}

add_projection_class_label <- function(dt) {
  if (!is.null(dt) && "projection_class" %in% names(dt)) {
    dt[, projection_class_label := label_projection_class(projection_class)]
  }
  dt
}

resolve_merged_col <- function(dt, cc) {
  if (cc %in% names(dt)) return(dt)
  x <- paste0(cc, ".x")
  y <- paste0(cc, ".y")
  if (x %in% names(dt) && y %in% names(dt)) {
    dt[, (cc) := get(x)]
    dt[is.na(get(cc)) | get(cc) == "", (cc) := get(y)]
  } else if (x %in% names(dt)) {
    dt[, (cc) := get(x)]
  } else if (y %in% names(dt)) {
    dt[, (cc) := get(y)]
  }
  dt
}

merge_sample_lookup <- function(dt, lookup) {
  if (is.null(dt) || nrow(dt) == 0 || !"sample_key" %in% names(dt)) return(dt)
  add_cols <- setdiff(c("dataset_folder", "modality"), names(dt))
  if (length(add_cols) == 0) return(dt)
  merge(dt, unique(lookup[, c("sample_key", add_cols), with = FALSE]), by = "sample_key", all.x = TRUE)
}

augment_pca_scores <- function(pca_scores, pred_lookup, pca_scores_path) {
  score_list <- list()
  if (!is.null(pca_scores) && nrow(pca_scores) > 0) {
    dt <- copy(pca_scores)
    dt <- merge_sample_lookup(dt, pred_lookup)
    dt[, pca_score_source := basename(pca_scores_path)]
    score_list[[length(score_list) + 1L]] <- dt
  }

  sidecar_paths <- c(
    file.path(projection_dir, "PRJNA863443", "prjna863443_pca_scores.tsv")
  )
  for (path in sidecar_paths) {
    sidecar <- read_optional_quiet(path)
    if (is.null(sidecar) || nrow(sidecar) == 0 || !"sample_key" %in% names(sidecar)) next
    sidecar <- merge_sample_lookup(sidecar, pred_lookup)
    sidecar[, pca_score_source := file.path(basename(dirname(path)), basename(path))]
    score_list[[length(score_list) + 1L]] <- sidecar
  }

  if (length(score_list) == 0) return(NULL)
  out <- rbindlist(score_list, use.names = TRUE, fill = TRUE)
  if (!"dataset_folder" %in% names(out)) out[, dataset_folder := NA_character_]
  if (!"modality" %in% names(out)) out[, modality := NA_character_]
  out[, pca_score_source_priority := fifelse(grepl("PRJNA863443/", pca_score_source, fixed = TRUE), 1L, 0L)]
  key_cols <- intersect(c("sample_key", "model_family", "selector"), names(out))
  if (length(key_cols) > 0) {
    setorderv(out, c(key_cols, "pca_score_source_priority"))
    out <- unique(out, by = key_cols, fromLast = TRUE)
  }
  out[, pca_score_source_priority := NULL]

  pca_source_audit <- out[, .(
    n_rows = .N,
    n_samples = uniqueN(sample_key)
  ), by = .(pca_score_source, dataset_folder, modality)][order(pca_score_source, dataset_folder, modality)]
  write_tsv(pca_source_audit, "pca_projection_score_source_audit.tsv")
  out
}

load_reference_pca_scores <- function() {
  saved <- read_optional_quiet(file.path(projection_dir, "reference_pca_projection_scores.tsv"))
  if (!is.null(saved) && nrow(saved) > 0 && all(c("PC1", "PC2", "sample_key") %in% names(saved))) {
    saved[, reference_pca_score_source := "reference_pca_projection_scores.tsv"]
    return(saved)
  }

  pca_row <- deployment_audit[model_family == "pca_glmnet" & selector == "version_b_age_dataset_penalized"]
  if (nrow(pca_row) == 0 || !"artifact_path" %in% names(pca_row) || !file.exists(pca_row$artifact_path[[1]])) {
    warn("Reference PCA scores unavailable: deployment PCA artifact path not found.")
    return(NULL)
  }
  obj <- readRDS(pca_row$artifact_path[[1]])
  pca_fit <- obj$model$pca_fit %||% obj$pca
  if (is.null(pca_fit) || is.null(pca_fit$x) || !all(c("PC1", "PC2") %in% colnames(pca_fit$x))) {
    warn("Reference PCA scores unavailable: deployment PCA artifact lacks pca$x PC1/PC2.")
    return(NULL)
  }
  ref <- as.data.table(pca_fit$x[, c("PC1", "PC2"), drop = FALSE])
  ref[, sample_key := rownames(pca_fit$x)]
  setcolorder(ref, "sample_key")
  ref[, `:=`(
    model_family = "pca_glmnet",
    selector = "version_b_age_dataset_penalized",
    reference_pca_score_source = "deployment_artifact_pca_fit"
  )]

  metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
  if (file.exists(metadata_path)) {
    meta <- fread(metadata_path, na.strings = c("", "NA", "NaN"))
    keep <- intersect(c("sample_key", "dataset_folder", "modality", "age_week"), names(meta))
    ref <- merge(ref, unique(meta[, ..keep]), by = "sample_key", all.x = TRUE)
  } else {
    warn("Reference PCA metadata labels unavailable: ", metadata_path)
  }
  write_tsv(ref, "reference_pca_projection_scores.tsv")
  ref
}

write_tsv <- function(dt, filename) {
  path <- file.path(out_dir, filename)
  fwrite(dt, path, sep = "\t", na = "NA")
  tables_written <<- c(tables_written, path)
  invisible(path)
}

theme_projection <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text.x = element_text(angle = 35, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

save_plot <- function(p, prefix, width = 10, height = 7) {
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  ggsave(pdf_path, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  plot_prefixes_written <<- c(plot_prefixes_written, prefix)
  invisible(c(pdf_path, png_path))
}

extract_ggplot_legend <- function(p) {
  g <- ggplotGrob(p)
  idx <- which(vapply(g$grobs, function(x) x$name, character(1)) == "guide-box")
  if (length(idx) == 0) return(grid::nullGrob())
  g$grobs[[idx[[1]]]]
}

save_plot_fixed_body_wide_legend <- function(p, prefix, body_width = 11, canvas_width = 14, height = 7, legend_height = 1.15) {
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  p_body <- p + theme(legend.position = "none")
  legend <- extract_ggplot_legend(p)
  draw_composite <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      x = 0,
      y = legend_height / height,
      width = body_width / canvas_width,
      height = (height - legend_height) / height,
      just = c("left", "bottom")
    ))
    grid::grid.draw(ggplotGrob(p_body))
    grid::popViewport()
    grid::pushViewport(grid::viewport(
      x = 0.5,
      y = (legend_height / 2) / height,
      width = 1,
      height = legend_height / height,
      just = c("center", "center")
    ))
    grid::grid.draw(legend)
    grid::popViewport()
  }
  grDevices::pdf(pdf_path, width = canvas_width, height = height)
  draw_composite()
  grDevices::dev.off()
  grDevices::png(png_path, width = canvas_width, height = height, units = "in", res = 300, bg = "white")
  draw_composite()
  grDevices::dev.off()
  plot_prefixes_written <<- c(plot_prefixes_written, prefix)
  invisible(c(pdf_path, png_path))
}

save_root_plot <- function(p, prefix, width = 10, height = 7) {
  pdf_path <- file.path(root_plot_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(root_plot_dir, paste0(prefix, ".png"))
  ggsave(pdf_path, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  invisible(c(pdf_path, png_path))
}

save_root_plot_fixed_body_wide_legend <- function(p, prefix, body_width = 10.5, canvas_width = 14, height = 8, legend_height = 1.15) {
  pdf_path <- file.path(root_plot_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(root_plot_dir, paste0(prefix, ".png"))
  p_body <- p + theme(legend.position = "none")
  legend <- extract_ggplot_legend(p)
  draw_composite <- function() {
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      x = 0,
      y = legend_height / height,
      width = body_width / canvas_width,
      height = (height - legend_height) / height,
      just = c("left", "bottom")
    ))
    grid::grid.draw(ggplotGrob(p_body))
    grid::popViewport()
    grid::pushViewport(grid::viewport(
      x = 0.5,
      y = (legend_height / 2) / height,
      width = 1,
      height = legend_height / height,
      just = c("center", "center")
    ))
    grid::grid.draw(legend)
    grid::popViewport()
  }
  grDevices::pdf(pdf_path, width = canvas_width, height = height)
  draw_composite()
  grDevices::dev.off()
  grDevices::png(png_path, width = canvas_width, height = height, units = "in", res = 300, bg = "white")
  draw_composite()
  grDevices::dev.off()
  invisible(c(pdf_path, png_path))
}

is_bad_ood <- function(x) {
  y <- tolower(as.character(x))
  !is.na(y) & y != "" & !y %chin% c("ok", "false", "in_distribution", "not_ood", "none", "0")
}

model_artifact_manifest <- read_required(file.path(artifact_dir, "model_artifact_manifest.tsv"))
calibration_parameters <- read_required(file.path(artifact_dir, "calibration_parameters.tsv"))
artifact_repro_check <- read_required(file.path(artifact_dir, "artifact_reproducibility_check.tsv"))
deployment_audit <- read_required(file.path(artifact_dir, "deployment_artifact_audit.tsv"))
pred_all <- read_required(prefer_simplified(
  "external_projection_predictions_all.tsv",
  "external_projection_predictions_all_simplified.tsv"
))
pred <- read_required(prefer_simplified(
  "external_projection_predictions.tsv",
  "external_projection_predictions_simplified.tsv"
))
skipped <- read_required(file.path(projection_dir, "external_projection_skipped_unaligned_samples.tsv"))
projection_summary_path <- file.path(projection_dir, "external_projection_summary.md")
if (!file.exists(projection_summary_path)) warn("Required markdown input not found: external_projection_summary.md")

relative_all <- read_optional(prefer_simplified(
  "external_projection_relative_to_controls_all.tsv",
  "external_projection_relative_to_controls_all_simplified.tsv"
))
relative <- read_optional(prefer_simplified(
  "external_projection_relative_to_controls.tsv",
  "external_projection_relative_to_controls_simplified.tsv"
))
dataset_control_availability <- read_optional(prefer_simplified(
  "external_projection_dataset_control_availability.tsv",
  "external_projection_dataset_control_availability_simplified.tsv"
))
pca_scores_path <- prefer_simplified(
  "pca_projection_scores.tsv",
  "pca_projection_scores_simplified.tsv"
)
pca_scores <- read_optional(pca_scores_path)
pca_scores <- augment_pca_scores(pca_scores, pred_all, pca_scores_path)
projection_artifact_audit <- read_optional(file.path(projection_dir, "projection_artifact_audit.tsv"))

required_external_projection_datasets <- c(
  "CRA003645", "CRA017492", "GSE194196", "GSE252791", "GSE273724",
  "GSE279446", "GSE279573", "GSE329024", "GSE228975", "GSE240547",
  "PRJNA1076509", "PRJNA863443"
)

optional_found <- data.table(
  optional_file = list.files(
    projection_dir,
    pattern = "(coverage|ood|pc|projection_artifact_audit|external_projection_relative_to_controls|pca_projection_scores).*\\.tsv$",
    full.names = FALSE,
    ignore.case = TRUE
  )
)
write_tsv(optional_found, "optional_projection_files_found.tsv")

needed_cols <- c(
  "sample_key", "dataset_folder", "modality", "age_week", "model_family", "selector",
  "predicted_age_week", "delta_age_week", "projection_class", "projection_reason",
  "treatment", "intervention", "genotype", "diet", "stress", "microbiome", "postpartum",
  "lactation", "gene_coverage_fraction", "projection_distance", "projection_distance_percentile",
  "OOD_flag", "projection_gene_coverage_flag"
)
pred_all <- ensure_cols(pred_all, c(needed_cols, "primary_external_projection_include", "matched_control_allowed"))
pred_all <- coerce_num(pred_all, c(
  "age_week", "predicted_age_week", "delta_age_week", "gene_coverage_fraction",
  "projection_distance", "projection_distance_percentile"
))
if (all(is.na(pred_all$delta_age_week))) {
  pred_all[, delta_age_week := predicted_age_week - age_week]
}
if (any(is.na(pred_all$delta_age_week) & is.finite(pred_all$predicted_age_week) & is.finite(pred_all$age_week))) {
  pred_all[is.na(delta_age_week) & is.finite(predicted_age_week) & is.finite(age_week),
    delta_age_week := predicted_age_week - age_week
  ]
}
pred_all[, primary_external_projection_include := as.logical(primary_external_projection_include)]
pred_all[is.na(primary_external_projection_include), primary_external_projection_include := TRUE]
if ("plot_include_primary" %in% names(pred_all)) {
  pred_all[, primary_external_projection_include := as.logical(plot_include_primary)]
  pred_all[is.na(primary_external_projection_include), primary_external_projection_include := FALSE]
}
pred_all[, matched_control_allowed := as.logical(matched_control_allowed)]
pred_all[is.na(matched_control_allowed), matched_control_allowed := TRUE]

pred <- ensure_cols(pred, needed_cols)
pred <- coerce_num(pred, c(
  "age_week", "predicted_age_week", "delta_age_week", "gene_coverage_fraction",
  "projection_distance", "projection_distance_percentile"
))
if (all(is.na(pred$delta_age_week))) {
  pred[, delta_age_week := predicted_age_week - age_week]
}
if (any(is.na(pred$delta_age_week) & is.finite(pred$predicted_age_week) & is.finite(pred$age_week))) {
  pred[is.na(delta_age_week) & is.finite(predicted_age_week) & is.finite(age_week),
    delta_age_week := predicted_age_week - age_week
  ]
}
if ("plot_include_primary" %in% names(pred)) {
  pred <- pred[as.logical(plot_include_primary) == TRUE]
}

expected_prediction_rows <- nrow(pred)
expected_unique_projected_samples <- uniqueN(pred$sample_key)
n_projected_all <- uniqueN(pred_all$sample_key)
n_prediction_rows_all <- nrow(pred_all)
n_projected <- uniqueN(pred$sample_key)
n_model_families <- uniqueN(pred$model_family)
if (n_model_families != 2) warn("Expected exactly two retained model families but found ", n_model_families, ".")

expected_cal <- data.table(
  model_family = c("elasticnet", "pca_glmnet"),
  selector = c("hvg_baseline", "version_b_age_dataset_penalized"),
  expected_intercept = c(-2.4891763, -0.3013835),
  expected_slope = c(1.087985, 1.010653)
)
deploy_cal <- calibration_parameters[artifact_scope == "deployment_all_control_like"]
deploy_cal <- merge(deploy_cal, expected_cal, by = c("model_family", "selector"), all.x = TRUE)
deploy_cal[, intercept_matches_expected := is.na(expected_intercept) | abs(intercept - expected_intercept) < 1e-4]
deploy_cal[, slope_matches_expected := is.na(expected_slope) | abs(slope - expected_slope) < 1e-4]
if (any(deploy_cal$intercept_matches_expected == FALSE | deploy_cal$slope_matches_expected == FALSE, na.rm = TRUE)) {
  warn("Deployment calibration parameters differ from expected values.")
}

deployment_paths_ok <- TRUE
if (!is.null(projection_artifact_audit) && all(c("requirement", "detail") %in% names(projection_artifact_audit))) {
  artifact_rows <- projection_artifact_audit[grepl("model RDS", requirement)]
  deployment_paths_ok <- all(grepl("deployment", artifact_rows$detail), na.rm = TRUE)
}
if (!deployment_paths_ok) warn("Could not confirm that external projection used deployment artifacts.")

training_overlap <- NA_integer_
all_clock_predictions_path <- file.path(results_dir, "all_clock_predictions.tsv")
if (file.exists(all_clock_predictions_path)) {
  training_keys <- fread(all_clock_predictions_path, sep = "\t", select = "sample_key")
  pred_for_overlap <- pred
  if ("projection_role" %in% names(pred_for_overlap)) {
    pred_for_overlap <- pred_for_overlap[projection_role == "external_projected_treatment"]
  }
  training_overlap <- length(intersect(unique(pred_for_overlap$sample_key), unique(training_keys$sample_key)))
  if (training_overlap > 0) warn(training_overlap, " external treatment projection sample keys overlap all_clock_predictions.tsv sample keys.")
} else {
  warn("Cannot audit projection/training sample overlap because all_clock_predictions.tsv is absent.")
}

deployment_control_like_ok <- NA
if ("all_training_control_like" %in% names(deployment_audit)) {
  deployment_control_like_ok <- all(as.logical(deployment_audit$all_training_control_like), na.rm = TRUE)
}

validation_checks <- data.table(
  check = c(
    "prediction_rows_primary_inferred",
    "unique_projected_samples_primary_inferred",
    "prediction_rows_all_inferred",
    "unique_projected_samples_all_inferred",
    "two_retained_model_families",
    "deployment_calibration_matches_expected",
    "projection_artifacts_are_deployment",
    "no_projected_samples_overlap_clock_training_predictions",
    "deployment_training_set_all_control_like"
  ),
  passed = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    n_model_families == 2,
    all(deploy_cal$intercept_matches_expected & deploy_cal$slope_matches_expected, na.rm = TRUE),
    deployment_paths_ok,
    ifelse(is.na(training_overlap), NA, training_overlap == 0),
    deployment_control_like_ok
  ),
  detail = c(
    paste0("observed=", nrow(pred)),
    paste0("observed=", n_projected),
    paste0("observed=", n_prediction_rows_all),
    paste0("observed=", n_projected_all),
    paste0("observed=", paste(sort(unique(pred$model_family)), collapse = ", ")),
    paste0("deployment rows=", nrow(deploy_cal)),
    if (!is.null(projection_artifact_audit)) "projection_artifact_audit.tsv inspected" else "audit unavailable",
    ifelse(is.na(training_overlap), "all_clock_predictions.tsv unavailable", paste0("overlap=", training_overlap)),
    paste0("all_training_control_like=", deployment_control_like_ok)
  )
)
write_tsv(validation_checks, "projection_diagnostic_validation_checks.tsv")

reference_ranges <- calibration_parameters[
  artifact_scope == "deployment_all_control_like",
  .(model_family, selector, reference_age_min = calibration_age_min, reference_age_max = calibration_age_max)
]

pred <- merge(pred, reference_ranges, by = c("model_family", "selector"), all.x = TRUE)
pred[, low_coverage_flag := fifelse(
  (!is.na(gene_coverage_fraction) & gene_coverage_fraction < 0.80) |
    grepl("low", projection_gene_coverage_flag %||% "", ignore.case = TRUE),
  TRUE, FALSE
)]
pred[, OOD_flag_bool := is_bad_ood(OOD_flag)]
pred[, outside_reference_age_range_flag := fifelse(
  is.finite(age_week) & is.finite(reference_age_min) & is.finite(reference_age_max) &
    (age_week < reference_age_min | age_week > reference_age_max),
  TRUE, FALSE
)]
pred[, claim_eligibility_flag := fifelse(
  low_coverage_flag | OOD_flag_bool | outside_reference_age_range_flag,
  "qc_caution",
  "eligible_for_projection_description"
)]

metadata_cols <- c(
  "treatment", "intervention", "genotype", "diet", "stress", "microbiome",
  "postpartum", "lactation", "sample_context", "sample_title"
)
pred <- ensure_cols(pred, metadata_cols)
group_metadata_cols <- intersect(c(
  "treatment", "intervention", "genotype", "diet", "stress", "microbiome", "postpartum", "lactation"
), names(pred))
active_group_metadata_cols <- group_metadata_cols[
  vapply(group_metadata_cols, function(cc) any(!is.na(pred[[cc]]) & pred[[cc]] != ""), logical(1))
]
if (length(active_group_metadata_cols) == 0) active_group_metadata_cols <- "projection_reason"

qc_cols <- unique(c(
  "sample_key", "dataset_folder", "modality", "age_week", "projection_class", "projection_reason",
  metadata_cols, "model_family", "selector", "predicted_age_week", "delta_age_week",
  "gene_coverage_fraction", "projection_distance", "projection_distance_percentile", "OOD_flag",
  "low_coverage_flag", "outside_reference_age_range_flag", "claim_eligibility_flag"
))
qc_cols <- qc_cols[qc_cols %in% names(pred)]
projection_qc <- pred[, ..qc_cols]
setorder(projection_qc, dataset_folder, sample_key, model_family)
write_tsv(projection_qc, "projection_qc_sample_level.tsv")

group_cols <- unique(c(
  "model_family", "selector", "dataset_folder", "modality", "projection_class",
  active_group_metadata_cols
))
projection_group_summary <- pred[, .(
  n = uniqueN(sample_key),
  n_prediction_rows = .N,
  mean_age_week = mean(age_week, na.rm = TRUE),
  mean_predicted_age_week = mean(predicted_age_week, na.rm = TRUE),
  mean_delta_age_week = mean(delta_age_week, na.rm = TRUE),
  median_delta_age_week = median(delta_age_week, na.rm = TRUE),
  sd_delta_age_week = sd(delta_age_week, na.rm = TRUE),
  mean_gene_coverage = mean(gene_coverage_fraction, na.rm = TRUE),
  n_low_coverage = sum(low_coverage_flag, na.rm = TRUE),
  n_OOD = sum(OOD_flag_bool, na.rm = TRUE)
), by = group_cols]
setorder(projection_group_summary, model_family, dataset_folder, projection_class)
write_tsv(projection_group_summary, "projection_group_summary.tsv")

sample_meta <- unique(pred[, c(
  "sample_key", "dataset_folder", "modality", "age_week", "projection_class",
  "projection_reason", metadata_cols
), with = FALSE])
delta_wide <- dcast(
  pred,
  sample_key ~ model_family,
  value.var = "delta_age_week",
  fun.aggregate = function(x) if (length(x) == 0) NA_real_ else x[[1]]
)
setnames(delta_wide, intersect(c("elasticnet", "pca_glmnet"), names(delta_wide)),
  paste0("delta_age_", intersect(c("elasticnet", "pca_glmnet"), names(delta_wide)))
)
model_concordance <- merge(sample_meta, delta_wide, by = "sample_key", all.x = TRUE)
model_concordance <- ensure_cols(model_concordance, c("delta_age_elasticnet", "delta_age_pca_glmnet"), NA_real_)
model_concordance[, delta_difference := delta_age_pca_glmnet - delta_age_elasticnet]
model_concordance[, same_direction_flag := fifelse(
  is.finite(delta_age_elasticnet) & is.finite(delta_age_pca_glmnet) &
    sign(delta_age_elasticnet) == sign(delta_age_pca_glmnet),
  TRUE, FALSE
)]
model_concordance[, absolute_disagreement := abs(delta_difference)]
setorder(model_concordance, -absolute_disagreement)
write_tsv(model_concordance, "model_concordance_sample_level.tsv")

concordance_group_cols <- unique(c("dataset_folder", "modality", "projection_class", active_group_metadata_cols))
model_concordance_group <- model_concordance[, .(
  n = .N,
  mean_delta_elasticnet = mean(delta_age_elasticnet, na.rm = TRUE),
  mean_delta_pca_glmnet = mean(delta_age_pca_glmnet, na.rm = TRUE),
  correlation = if (.N >= 3 && uniqueN(delta_age_elasticnet[is.finite(delta_age_elasticnet)]) > 1 &&
    uniqueN(delta_age_pca_glmnet[is.finite(delta_age_pca_glmnet)]) > 1) {
    suppressWarnings(cor(delta_age_elasticnet, delta_age_pca_glmnet, use = "complete.obs"))
  } else NA_real_,
  fraction_same_direction = mean(same_direction_flag, na.rm = TRUE),
  mean_absolute_disagreement = mean(absolute_disagreement, na.rm = TRUE)
), by = concordance_group_cols]
setorder(model_concordance_group, dataset_folder, projection_class)
write_tsv(model_concordance_group, "model_concordance_group_summary.tsv")

if (!is.null(relative)) {
  before_relative <- nrow(relative)
  relative <- unique(relative)
  if (nrow(relative) < before_relative) {
    warn("Removed ", before_relative - nrow(relative), " exact duplicate rows from external_projection_relative_to_controls.tsv for diagnostics.")
  }
  if ("plot_include_primary" %in% names(relative)) {
    relative <- relative[as.logical(plot_include_primary) == TRUE]
  }
  relative <- coerce_num(relative, c(
    "age_week", "delta_age_week", "control_expected_delta_age", "relative_delta_age", "n_matched_controls"
  ))
  rel_key <- intersect(c("sample_key", "model_family", "selector"), names(relative))
  if (length(rel_key) == 3) {
    relative <- relative[, .SD[1], by = rel_key]
  }
  rel_meta <- unique(pred[, c(
    "sample_key", "model_family", "selector", "dataset_folder", "modality",
    "projection_class", "delta_age_week", metadata_cols
  ), with = FALSE])
  relative <- merge(relative, rel_meta, by = c("sample_key", "model_family", "selector"), all.x = TRUE, suffixes = c("", ".pred"))
  for (cc in c("dataset_folder", "modality", "projection_class", "delta_age_week", metadata_cols)) {
    alt <- paste0(cc, ".pred")
    if (alt %in% names(relative)) {
      if (!cc %in% names(relative)) relative[, (cc) := get(alt)]
      relative[is.na(get(cc)) | get(cc) == "", (cc) := get(alt)]
      relative[, (alt) := NULL]
    }
  }
  relative_group_cols <- unique(c("model_family", "selector", "dataset_folder", "modality", "projection_class", active_group_metadata_cols))
  matched_control_relative_delta_summary <- relative[, .(
    n = uniqueN(sample_key),
    mean_raw_delta_age_week = mean(delta_age_week, na.rm = TRUE),
    median_raw_delta_age_week = median(delta_age_week, na.rm = TRUE),
    mean_control_expected_delta_age = mean(control_expected_delta_age, na.rm = TRUE),
    median_control_expected_delta_age = median(control_expected_delta_age, na.rm = TRUE),
    mean_relative_delta_age = mean(relative_delta_age, na.rm = TRUE),
    median_relative_delta_age = median(relative_delta_age, na.rm = TRUE),
    sd_relative_delta_age = sd(relative_delta_age, na.rm = TRUE),
    mean_n_matched_controls = mean(n_matched_controls, na.rm = TRUE)
  ), by = relative_group_cols]
} else {
  matched_control_relative_delta_summary <- data.table(
    note = "external_projection_relative_to_controls.tsv was not available; matched-control summaries and plots were skipped."
  )
}
write_tsv(matched_control_relative_delta_summary, "matched_control_relative_delta_summary.tsv")

audit_relative <- if (!is.null(relative_all)) copy(relative_all) else if (!is.null(relative)) copy(relative) else data.table()
audit_relative <- ensure_cols(audit_relative, c(
  "dataset_folder", "model_family", "matched_control_method", "primary_external_projection_include",
  "matched_control_allowed", "relative_delta_age_week", "relative_delta_age", "n_matched_controls",
  "treatment", "sample_title", "sample_key"
))
if (nrow(audit_relative) > 0) {
  audit_relative <- coerce_num(audit_relative, c("relative_delta_age_week", "relative_delta_age", "n_matched_controls"))
  if (all(is.na(audit_relative$relative_delta_age_week)) && any(is.finite(audit_relative$relative_delta_age))) {
    audit_relative[, relative_delta_age_week := relative_delta_age]
  }
  audit_relative[, primary_external_projection_include := as.logical(primary_external_projection_include)]
  audit_relative[, matched_control_allowed := as.logical(matched_control_allowed)]
  audit_relative[, invalid_no_matched_control := matched_control_method == "no_matched_control" & (
    is.finite(relative_delta_age_week) | n_matched_controls != 0
  )]
  audit_relative[, primary_excluded_with_relative_delta := primary_external_projection_include == FALSE & is.finite(relative_delta_age_week)]
  audit_relative[, matched_disallowed_with_relative_delta := matched_control_allowed == FALSE & is.finite(relative_delta_age_week)]
  audit_relative[, invalid_any := invalid_no_matched_control | primary_excluded_with_relative_delta | matched_disallowed_with_relative_delta]

  matched_control_integrity_audit <- audit_relative[, .(
    n_rows = .N,
    n_no_matched_control = sum(matched_control_method == "no_matched_control", na.rm = TRUE),
    n_invalid_no_matched_control = sum(invalid_no_matched_control, na.rm = TRUE),
    n_non_na_relative_delta = sum(is.finite(relative_delta_age_week), na.rm = TRUE),
    n_primary_excluded_with_relative_delta = sum(primary_excluded_with_relative_delta, na.rm = TRUE),
    pass_fail = fifelse(any(invalid_any, na.rm = TRUE), "FAIL", "PASS")
  ), by = .(dataset_folder, model_family)]
} else {
  audit_relative[, `:=`(
    invalid_no_matched_control = logical(),
    primary_excluded_with_relative_delta = logical(),
    matched_disallowed_with_relative_delta = logical(),
    invalid_any = logical()
  )]
  matched_control_integrity_audit <- data.table(
    dataset_folder = character(), model_family = character(), n_rows = integer(),
    n_no_matched_control = integer(), n_invalid_no_matched_control = integer(),
    n_non_na_relative_delta = integer(), n_primary_excluded_with_relative_delta = integer(),
    pass_fail = character()
  )
}
write_tsv(matched_control_integrity_audit, "matched_control_integrity_audit.tsv")

pred_all_audit <- ensure_cols(copy(pred_all), c(
  "dataset_folder", "sample_key", "primary_external_projection_include",
  "treatment", "genotype", "sample_title", "projection_class", "projection_role",
  "matched_reference_control_flag", "reference_control_group_id", "gse279446_cohort",
  "projection_group", "comparison_id"
))
if (!"gse279446_cohort" %in% names(audit_relative)) audit_relative[, gse279446_cohort := NA_character_]
if (!"projection_role" %in% names(audit_relative)) audit_relative[, projection_role := NA_character_]
if (!"projection_group" %in% names(audit_relative)) audit_relative[, projection_group := NA_character_]
if (!"comparison_id" %in% names(audit_relative)) audit_relative[, comparison_id := reference_control_group_id]

audit_rows <- list()
add_audit <- function(dataset_folder, check_name, observed, expected) {
  audit_rows[[length(audit_rows) + 1L]] <<- data.table(
    dataset_folder = dataset_folder,
    check_name = check_name,
    observed = as.character(observed),
    expected = as.character(expected)
  )
}
eq_count <- function(x) as.integer(x)

for (ds in c("GSE194196", "GSE252791")) {
  method <- if (ds == "GSE194196") "GSE194196_WT_1_4_7" else "GSE252791_WT_1_4_7"
  add_audit(ds, "three_WT_reference_controls_detected",
    uniqueN(pred_all_audit[dataset_folder == ds & matched_reference_control_flag == TRUE]$sample_key), 3)
  add_audit(ds, "treatment_rows_matched_to_named_WT_controls",
    sum(audit_relative[dataset_folder == ds & projection_role == "external_projected_treatment", matched_control_method == method], na.rm = TRUE),
    nrow(audit_relative[dataset_folder == ds & projection_role == "external_projected_treatment"]))
  add_audit(ds, "treatment_rows_n_matched_controls_3",
    sum(audit_relative[dataset_folder == ds & projection_role == "external_projected_treatment", n_matched_controls == 3], na.rm = TRUE),
    nrow(audit_relative[dataset_folder == ds & projection_role == "external_projected_treatment"]))
}

add_audit("GSE279446", "all_projected_samples_detected", uniqueN(pred_all_audit[dataset_folder == "GSE279446"]$sample_key), 24)
add_audit("GSE279446", "adult_and_young_cohorts_detected",
  paste(sort(unique(na.omit(pred_all_audit[dataset_folder == "GSE279446"]$gse279446_cohort))), collapse = ";"),
  "adult_cohort;young_cohort")
add_audit("GSE279446", "adult_mutants_matched_to_adult_WT",
  sum(audit_relative[dataset_folder == "GSE279446" & gse279446_cohort == "adult_cohort" & grepl("heterozygous|homozygous", projection_class), matched_control_method == "GSE279446_adult_WT_ovary_controls"], na.rm = TRUE),
  nrow(audit_relative[dataset_folder == "GSE279446" & gse279446_cohort == "adult_cohort" & grepl("heterozygous|homozygous", projection_class)]))
add_audit("GSE279446", "young_mutants_matched_to_young_WT",
  sum(audit_relative[dataset_folder == "GSE279446" & gse279446_cohort == "young_cohort" & grepl("heterozygous|homozygous", projection_class), matched_control_method == "GSE279446_young_WT_controls"], na.rm = TRUE),
  nrow(audit_relative[dataset_folder == "GSE279446" & gse279446_cohort == "young_cohort" & grepl("heterozygous|homozygous", projection_class)]))
add_audit("GSE279446", "Zeb1_WT_never_heterozygous",
  sum(pred_all_audit[dataset_folder == "GSE279446" & grepl("Zeb1_WT|WT_ovary|WT ovary", paste(sample_key, treatment, genotype, sample_title, sep = " "), ignore.case = TRUE), grepl("heterozygous", projection_class, ignore.case = TRUE)], na.rm = TRUE),
  0)

add_audit("GSE279573", "seven_samples_in_predictions_all", uniqueN(pred_all_audit[dataset_folder == "GSE279573"]$sample_key), 7)
add_audit("GSE279573", "two_WT_reference_controls", uniqueN(pred_all_audit[dataset_folder == "GSE279573" & matched_reference_control_flag == TRUE]$sample_key), 2)
add_audit("GSE279573", "five_treatment_samples", uniqueN(pred_all_audit[dataset_folder == "GSE279573" & projection_role == "external_projected_treatment"]$sample_key), 5)
add_audit("GSE279573", "treatment_rows_n_matched_controls_2",
  sum(audit_relative[dataset_folder == "GSE279573" & projection_role == "external_projected_treatment", n_matched_controls == 2], na.rm = TRUE),
  nrow(audit_relative[dataset_folder == "GSE279573" & projection_role == "external_projected_treatment"]))

add_audit("GSE240547", "four_control_diet_samples", uniqueN(pred_all_audit[dataset_folder == "GSE240547" & projection_class == "control_diet"]$sample_key), 4)
add_audit("GSE240547", "four_HFHS_samples", uniqueN(pred_all_audit[dataset_folder == "GSE240547" & projection_class == "diet_intervention"]$sample_key), 4)
add_audit("GSE240547", "HFHS_n_matched_controls_4",
  sum(audit_relative[dataset_folder == "GSE240547" & projection_class == "diet_intervention", n_matched_controls == 4], na.rm = TRUE),
  nrow(audit_relative[dataset_folder == "GSE240547" & projection_class == "diet_intervention"]))

add_audit("PRJNA1076509", "fourteen_samples_detected", uniqueN(pred_all_audit[dataset_folder == "PRJNA1076509"]$sample_key), 14)
add_audit("PRJNA1076509", "twenty_eight_relative_rows", nrow(audit_relative[dataset_folder == "PRJNA1076509"]), 28)
add_audit("PRJNA1076509", "all_method_young_female_donor_microbiota",
  sum(audit_relative[dataset_folder == "PRJNA1076509", matched_control_method == "PRJNA1076509_young_female_donor_microbiota"], na.rm = TRUE), 28)
add_audit("PRJNA1076509", "all_n_matched_controls_7",
  sum(audit_relative[dataset_folder == "PRJNA1076509", n_matched_controls == 7], na.rm = TRUE), 28)

add_audit("PRJNA863443", "thirty_three_samples_detected", uniqueN(pred_all_audit[dataset_folder == "PRJNA863443"]$sample_key), 33)
add_audit("PRJNA863443", "sixty_six_prediction_rows", nrow(pred_all_audit[dataset_folder == "PRJNA863443"]), 66)
add_audit("PRJNA863443", "four_VCD_comparisons_detected",
  uniqueN(audit_relative[dataset_folder == "PRJNA863443" & grepl("^PRJNA863443_VCD_vs_CTL_", comparison_id)]$comparison_id), 4)
add_audit("PRJNA863443", "VCD_treatment_rows_have_two_controls",
  sum(audit_relative[
    dataset_folder == "PRJNA863443" &
      grepl("^PRJNA863443_VCD_vs_CTL_", comparison_id) &
      projection_class == "intervention",
    n_matched_controls == 2
  ], na.rm = TRUE),
  nrow(audit_relative[
    dataset_folder == "PRJNA863443" &
      grepl("^PRJNA863443_VCD_vs_CTL_", comparison_id) &
      projection_class == "intervention"
  ]))
add_audit("PRJNA863443", "two_Foxl2_comparisons_detected",
  uniqueN(audit_relative[dataset_folder == "PRJNA863443" & grepl("^PRJNA863443_Foxl2_het_vs_WT_", comparison_id)]$comparison_id), 2)
add_audit("PRJNA863443", "Foxl2_het_rows_have_expected_controls",
  nrow(audit_relative[
    dataset_folder == "PRJNA863443" &
      grepl("^PRJNA863443_Foxl2_het_vs_WT_", comparison_id) &
      projection_class == "genetic_model" &
      ((grepl("_young$", comparison_id) & n_matched_controls == 2) |
         (grepl("_mid_age$", comparison_id) & n_matched_controls == 3))
  ]),
  nrow(audit_relative[
    dataset_folder == "PRJNA863443" &
      grepl("^PRJNA863443_Foxl2_het_vs_WT_", comparison_id) &
      projection_class == "genetic_model"
  ]))
add_audit("PRJNA863443", "natural_ageing_projected_no_relative_delta",
  nrow(audit_relative[
    dataset_folder == "PRJNA863443" &
      comparison_id == "PRJNA863443_natural_aging_reference" &
      matched_control_method == "no_matched_control" &
      !is.finite(relative_delta_age_week)
  ]),
  nrow(audit_relative[
    dataset_folder == "PRJNA863443" &
      comparison_id == "PRJNA863443_natural_aging_reference"
  ]))

add_audit("GSE228975", "primary_external_projection_include_FALSE",
  sum(pred_all_audit[dataset_folder == "GSE228975", primary_external_projection_include == FALSE], na.rm = TRUE),
  nrow(pred_all_audit[dataset_folder == "GSE228975"]))
add_audit("GSE228975", "relative_delta_age_week_NA",
  sum(audit_relative[dataset_folder == "GSE228975", !is.finite(relative_delta_age_week)], na.rm = TRUE),
  nrow(audit_relative[dataset_folder == "GSE228975"]))

special_dataset_control_audit <- rbindlist(audit_rows, use.names = TRUE, fill = TRUE)
special_dataset_control_audit[, pass_fail := fifelse(observed == expected, "PASS", "FAIL")]
write_tsv(special_dataset_control_audit, "special_dataset_control_audit.tsv")

if (!is.null(dataset_control_availability)) {
  dataset_relative_control_status <- copy(dataset_control_availability)
  dataset_relative_control_status <- ensure_cols(dataset_relative_control_status, c("dataset_folder", "observed_relative_method", "n_projected_samples_all", "n_projected_samples_primary", "pass_fail_status", "dataset_note"))
  rel_counts <- audit_relative[, .(
    n_matched_control_rows = sum(matched_control_method != "no_matched_control", na.rm = TRUE),
    n_unmatched_rows = sum(matched_control_method == "no_matched_control", na.rm = TRUE),
    n_invalid_rows = sum(invalid_any, na.rm = TRUE)
  ), by = dataset_folder]
  dataset_relative_control_status <- merge(
    data.table(dataset_folder = required_external_projection_datasets),
    dataset_relative_control_status,
    by = "dataset_folder",
    all.x = TRUE
  )
  dataset_relative_control_status <- merge(dataset_relative_control_status, rel_counts, by = "dataset_folder", all.x = TRUE)
	  dataset_relative_control_status[, expected_control_strategy := fifelse(
	    dataset_folder == "GSE228975", "no_matched_control",
	    fifelse(dataset_folder == "GSE240547", "GSE240547_control_diet",
	      fifelse(dataset_folder == "PRJNA1076509", "PRJNA1076509_young_female_donor_microbiota",
	        fifelse(dataset_folder == "GSE329024", "GSE329024_design_matched_controls",
	          fifelse(dataset_folder == "PRJNA863443", "VCD_same_time_CTL;Foxl2_same_age_WT;natural_aging_no_matched_control", "same_dataset_control_meta_when_available")
	        )
	      )
	    )
	  )]
  dataset_relative_control_status[, observed_control_strategy := observed_relative_method]
  dataset_relative_control_status[, pass_fail := pass_fail_status]
  dataset_relative_control_status[, note := dataset_note]
  for (cc in c("n_matched_control_rows", "n_unmatched_rows", "n_invalid_rows")) {
    dataset_relative_control_status[is.na(get(cc)), (cc) := 0L]
  }
  status_cols <- c(
    "dataset_folder", "expected_control_strategy", "observed_control_strategy",
    "n_projected_samples_all", "n_projected_samples_primary",
    "n_matched_control_rows", "n_unmatched_rows", "n_invalid_rows",
    "pass_fail", "note"
  )
  dataset_relative_control_status <- dataset_relative_control_status[, ..status_cols]
} else {
  dataset_relative_control_status <- data.table(dataset_folder = required_external_projection_datasets)
}
write_tsv(dataset_relative_control_status, "dataset_relative_control_status.tsv")

matched_pairs <- if (!is.null(relative)) {
  unique(relative[is.finite(relative_delta_age), .(sample_key, model_family, selector, matched_control_available = TRUE)])
} else {
  data.table(sample_key = character(), model_family = character(), selector = character(), matched_control_available = logical())
}
qc_for_grade <- merge(
  pred,
  matched_pairs,
  by = c("sample_key", "model_family", "selector"),
  all.x = TRUE
)
qc_for_grade[is.na(matched_control_available), matched_control_available := FALSE]
paired_group <- model_concordance_group[, c(concordance_group_cols, "fraction_same_direction", "mean_delta_elasticnet", "mean_delta_pca_glmnet"), with = FALSE]
claim_group_cols <- unique(c("dataset_folder", "modality", "projection_class", active_group_metadata_cols))
projection_claim_grading <- qc_for_grade[, .(
  n = uniqueN(sample_key),
  n_prediction_rows = .N,
  matched_controls_available = any(matched_control_available, na.rm = TRUE),
  mean_gene_coverage = mean(gene_coverage_fraction, na.rm = TRUE),
  low_coverage_fraction = mean(low_coverage_flag, na.rm = TRUE),
  OOD_fraction = mean(OOD_flag_bool, na.rm = TRUE),
  missing_key_metadata = any(is.na(dataset_folder) | is.na(projection_class) | is.na(modality))
), by = claim_group_cols]
projection_claim_grading <- merge(projection_claim_grading, paired_group, by = claim_group_cols, all.x = TRUE)
projection_claim_grading[, both_clocks_agree_in_direction := is.finite(mean_delta_elasticnet) &
  is.finite(mean_delta_pca_glmnet) &
  sign(mean_delta_elasticnet) == sign(mean_delta_pca_glmnet) &
  fraction_same_direction >= 0.75]
projection_claim_grading[, claim_grade := fifelse(
  missing_key_metadata | low_coverage_fraction > 0.5 | OOD_fraction > 0.5,
  "exclude_from_claims",
  fifelse(
    n >= 3 & matched_controls_available & (is.na(mean_gene_coverage) | mean_gene_coverage >= 0.8) &
      OOD_fraction <= 0.25 & both_clocks_agree_in_direction,
    "strong_projection_signal",
    fifelse(
      n >= 2 & (matched_controls_available | both_clocks_agree_in_direction) &
        low_coverage_fraction <= 0.25 & OOD_fraction <= 0.5,
      "moderate_projection_signal",
      "descriptive_only"
    )
  )
)]
projection_claim_grading[, grading_reason := fifelse(
  claim_grade == "exclude_from_claims",
  "major QC issue, high OOD fraction, low coverage, or missing metadata",
  fifelse(
    claim_grade == "strong_projection_signal",
    "n>=3, matched controls available, acceptable coverage/OOD, and both clocks agree",
    fifelse(
      claim_grade == "moderate_projection_signal",
      "n>=2 with matched controls or clock agreement and no major QC issue",
      "descriptive or supplementary interpretation only"
    )
  )
)]
setorder(projection_claim_grading, claim_grade, dataset_folder, projection_class)
write_tsv(projection_claim_grading, "projection_claim_grading.tsv")

skipped <- ensure_cols(skipped, c("dataset_folder", "projection_reason", "projection_class", "modality", "sample_key"))
skipped_projection_alignment_summary <- skipped[, .(
  n_skipped = .N,
  skipped_samples = paste(sample_key, collapse = "; ")
), by = .(dataset_folder, modality, projection_class, reason = projection_reason)]
setorder(skipped_projection_alignment_summary, dataset_folder, reason)
write_tsv(skipped_projection_alignment_summary, "skipped_projection_alignment_summary.tsv")

largest_abs_delta <- pred[order(-abs(delta_age_week))][seq_len(min(.N, 10))]
largest_abs_delta[, discrepancy_category := "largest_absolute_delta_age"]
largest_abs_delta[, ranking_value := abs(delta_age_week)]
largest_disagreement <- model_concordance[order(-absolute_disagreement)][seq_len(min(.N, 10))]
largest_disagreement[, discrepancy_category := "largest_elasticnet_pca_disagreement"]
largest_disagreement[, ranking_value := absolute_disagreement]
lowest_coverage <- pred[is.finite(gene_coverage_fraction)][order(gene_coverage_fraction)][seq_len(min(.N, 10))]
lowest_coverage[, discrepancy_category := "lowest_gene_coverage"]
lowest_coverage[, ranking_value := gene_coverage_fraction]
highest_ood <- pred[is.finite(projection_distance)][order(-projection_distance)][seq_len(min(.N, 10))]
highest_ood[, discrepancy_category := "highest_projection_distance"]
highest_ood[, ranking_value := projection_distance]
skipped_worst <- copy(skipped)
skipped_worst[, `:=`(
  discrepancy_category = "skipped_unaligned_sample",
  ranking_value = NA_real_,
  model_family = NA_character_,
  selector = NA_character_,
  predicted_age_week = NA_real_,
  delta_age_week = NA_real_,
  gene_coverage_fraction = NA_real_,
  projection_distance = NA_real_,
  OOD_flag = NA_character_,
  delta_age_elasticnet = NA_real_,
  delta_age_pca_glmnet = NA_real_,
  absolute_disagreement = NA_real_
)]
worst_cols <- unique(c(
  "discrepancy_category", "ranking_value", "sample_key", "dataset_folder", "modality",
  "age_week", "projection_class", "projection_reason", "model_family", "selector",
  "predicted_age_week", "delta_age_week", "delta_age_elasticnet", "delta_age_pca_glmnet",
  "absolute_disagreement", "gene_coverage_fraction", "projection_distance", "OOD_flag"
))
worst_discrepant_samples <- rbindlist(lapply(list(
  largest_abs_delta, largest_disagreement, lowest_coverage, highest_ood, skipped_worst
), function(x) {
  x <- ensure_cols(x, worst_cols)
  x[, ..worst_cols]
}), fill = TRUE)
write_tsv(worst_discrepant_samples, "worst_discrepant_samples.tsv")

pred <- add_projection_class_label(pred)
pred_all <- add_projection_class_label(pred_all)
relative <- add_projection_class_label(relative)
model_concordance <- add_projection_class_label(model_concordance)

projected_counts <- unique(pred_all[, .(sample_key, dataset_folder, primary_external_projection_include)])
projected_summary <- projected_counts[, .(
  n = .N,
  n_primary = sum(primary_external_projection_include == TRUE, na.rm = TRUE),
  n_non_primary = sum(primary_external_projection_include == FALSE, na.rm = TRUE)
), by = dataset_folder]
projected_bar_counts <- projected_summary[, .(
  dataset_folder,
  n,
  status = fifelse(n_non_primary > 0 & n_primary == 0, "OOD/not primary", "Dataset")
)]
skipped_counts <- unique(skipped[, .(sample_key, dataset_folder)])
skipped_counts <- skipped_counts[, .(n = .N, status = "skipped_unaligned"), by = dataset_folder]
flow_counts <- rbindlist(list(projected_bar_counts, skipped_counts), fill = TRUE)
flow_labels <- projected_summary[n_non_primary > 0]
flow_labels[, label := paste0("OOD/not primary: ", n_non_primary)]
p <- ggplot(flow_counts, aes(dataset_folder, n, fill = status)) +
  geom_col(position = "dodge", color = "grey25") +
  geom_text(
    data = flow_labels,
    aes(dataset_folder, n, label = label),
    inherit.aes = FALSE,
    vjust = -0.35,
    size = 3.2
  ) +
  expand_limits(y = max(flow_counts$n, na.rm = TRUE) * 1.12) +
  labs(
    title = "External projection cohort flow by dataset",
    subtitle = paste0(n_projected_all, " all projected samples; OOD/not-primary datasets are colored and labeled."),
    x = "Dataset",
    y = "Sample count",
    fill = "Status"
  ) +
  scale_fill_manual(values = c(
    "Dataset" = "#F8766D",
    "OOD/not primary" = "#4E79A7",
    "skipped_unaligned" = "grey70"
  )) +
  theme_projection()
save_plot(p, "01_projection_cohort_flow_dataset", width = 11, height = 6)
save_root_plot(p, "08_all_projection_cohort_flow_dataset", width = 11, height = 6)

class_counts <- unique(pred_all[, .(sample_key, projection_class, projection_class_label, modality, primary_external_projection_include)])
class_counts <- class_counts[, .(
  n = .N,
  n_non_primary = sum(primary_external_projection_include == FALSE, na.rm = TRUE)
), by = .(projection_class, projection_class_label, modality)]
class_grid <- CJ(
  projection_class = sort(unique(class_counts$projection_class)),
  modality = sort(unique(class_counts$modality))
)
class_grid[, projection_class_label := label_projection_class(projection_class)]
class_counts <- merge(class_grid, class_counts, by = c("projection_class", "projection_class_label", "modality"), all.x = TRUE)
class_counts[is.na(n), n := 0L]
class_counts[is.na(n_non_primary), n_non_primary := 0L]
class_labels <- class_counts[, .(
  n = max(n, na.rm = TRUE),
  n_non_primary = sum(n_non_primary, na.rm = TRUE)
), by = projection_class_label][n_non_primary > 0]
class_labels[, label := paste0("OOD/not primary: ", n_non_primary)]
p <- ggplot(class_counts, aes(projection_class_label, n, fill = modality)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.68, color = "grey25") +
  geom_text(
    data = class_labels,
    aes(projection_class_label, n, label = label),
    inherit.aes = FALSE,
    hjust = -0.05,
    size = 3.2
  ) +
  expand_limits(y = max(class_counts$n, na.rm = TRUE) * 1.15) +
  coord_flip() +
  labs(
    title = "Projected samples by projection class and modality",
    subtitle = paste0(n_projected_all, " all projected samples; non-primary/OOD classes are labeled."),
    x = "Projection class",
    y = "Projected sample count",
    fill = "Modality"
  ) +
  theme_projection() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )
save_plot(p, "01_projection_cohort_flow_class_modality", width = 10, height = 7)
save_root_plot(p, "09_all_projection_class_modality_primary_vs_ood", width = 10, height = 7)

has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
all_plot <- copy(pred_all)
all_plot[, interpretation_frame := fifelse(primary_external_projection_include == TRUE, "primary", "OOD/supplementary")]
all_plot_labels <- all_plot[
  primary_external_projection_include == FALSE |
    dataset_folder %chin% c("GSE240547", "PRJNA1076509"),
  .SD[order(dataset_folder, sample_key)][seq_len(min(.N, 4))],
  by = model_family
]
p <- ggplot(all_plot, aes(age_week, predicted_age_week, color = dataset_folder, shape = interpretation_frame)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(size = 2.4, alpha = 0.85) +
  facet_wrap(~model_family) +
  labs(
    title = "All external projections including OOD/supplementary samples",
    subtitle = "Primary and non-primary samples are both shown; GSE228975 is intentionally OOD/supplementary.",
    x = "Chronological age (weeks)",
    y = "Projected transcriptomic age (weeks)",
    color = "Dataset",
    shape = "Interpretation frame"
  ) +
  theme_projection()
if (nrow(all_plot_labels) > 0) {
  if (has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(data = all_plot_labels, aes(label = dataset_folder), size = 2.4, max.overlaps = 20, show.legend = FALSE)
  } else {
    p <- p + geom_text(data = all_plot_labels, aes(label = dataset_folder), size = 2.2, vjust = -0.7, show.legend = FALSE)
  }
}
  save_plot_fixed_body_wide_legend(p, "13_all_projection_predicted_vs_age_including_ood", body_width = 11, canvas_width = 16, height = 7)
  save_root_plot(p, "10_all_projection_predicted_vs_age_including_ood", width = 11, height = 7)

p <- ggplot(all_plot, aes(dataset_folder, delta_age_week, color = interpretation_frame)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, alpha = 0.12) +
  geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
  facet_wrap(~model_family, scales = "free_x") +
  labs(
    title = "All external raw delta age by dataset",
    subtitle = "Includes unmatched and OOD/supplementary datasets; use primary relative-delta plots only for matched-control interpretation.",
    x = "Dataset",
    y = "Delta age (projected - chronological weeks)",
    color = "Interpretation frame"
  ) +
  theme_projection()
save_plot(p, "14_all_projection_raw_delta_by_dataset_including_unmatched", width = 13, height = 7)
save_root_plot(p, "11_all_projection_raw_delta_by_dataset_including_unmatched", width = 13, height = 7)

if (nrow(dataset_relative_control_status) > 0 && all(c("n_matched_control_rows", "n_unmatched_rows") %in% names(dataset_relative_control_status))) {
  status_plot <- melt(
    dataset_relative_control_status,
    id.vars = intersect(c("dataset_folder", "pass_fail"), names(dataset_relative_control_status)),
    measure.vars = c("n_matched_control_rows", "n_unmatched_rows"),
    variable.name = "relative_control_status",
    value.name = "n_rows"
  )
  status_plot[, relative_control_status := fifelse(
    relative_control_status == "n_matched_control_rows",
    "matched_control_rows",
    "unmatched_or_no_relative_rows"
  )]
  p <- ggplot(status_plot, aes(dataset_folder, n_rows, fill = relative_control_status)) +
    geom_col(position = "stack", color = "grey25") +
    labs(
      title = "Matched-control availability by dataset",
      subtitle = "Unmatched/OOD rows are shown instead of disappearing from relative-delta plots.",
      x = "Dataset",
      y = "Sample x model rows",
      fill = "Relative-control status"
    ) +
    theme_projection()
  save_plot(p, "15_relative_control_status_by_dataset", width = 12, height = 6)
  save_root_plot(p, "12_relative_control_status_by_dataset", width = 12, height = 6)
}

outlier_labels <- pred[, .SD[order(-abs(delta_age_week))][seq_len(min(.N, 4))], by = model_family]
outlier_labels[, compact_sample_id := sub("^.*__", "", sample_key)]
outlier_labels[, outlier_label := paste(dataset_folder, compact_sample_id, sep = "-")]
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
p <- ggplot(pred, aes(age_week, predicted_age_week, color = projection_class_label, shape = modality)) +
  annotate(
    "rect",
    xmin = min(reference_ranges$reference_age_min, na.rm = TRUE),
    xmax = max(reference_ranges$reference_age_max, na.rm = TRUE),
    ymin = -Inf,
    ymax = Inf,
    alpha = 0.06,
    fill = "grey50"
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
  geom_point(size = 2.4, alpha = 0.85) +
  facet_wrap(~model_family) +
  labs(
    title = "External predicted versus chronological age",
    subtitle = "Grey band marks deployment calibration/reference age range.",
    x = "Chronological age (weeks)",
    y = "Projected transcriptomic age (weeks)",
    color = "Projection class",
    shape = "Modality"
  ) +
  theme_projection()
if (nrow(outlier_labels) > 0) {
	  if (has_ggrepel) {
	    p <- p + ggrepel::geom_text_repel(
	      data = outlier_labels,
	      aes(label = outlier_label),
	      size = 2.4,
	      max.overlaps = 20,
	      show.legend = FALSE
	    )
	  } else {
	    p <- p + geom_text(
	      data = outlier_labels,
	      aes(label = outlier_label),
	      size = 2.2,
	      vjust = -0.7,
	      show.legend = FALSE
	    )
	  }
	}
save_plot_fixed_body_wide_legend(p, "02_external_predicted_vs_age", body_width = 11, canvas_width = 14, height = 7)
save_root_plot(p, "01_external_projection_predicted_vs_chronological_age", width = 14, height = 7)

p <- ggplot(pred, aes(projection_class_label, delta_age_week, color = projection_class_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.18, height = 0, size = 1.8, alpha = 0.8) +
  facet_wrap(~model_family, scales = "free_x") +
  labs(
    title = "Raw delta age by projection class",
    x = "Projection class",
    y = "Delta age (projected - chronological weeks)",
    color = "Projection class"
  ) +
  theme_projection()
save_plot_fixed_body_wide_legend(p, "03_raw_delta_age_by_projection_class", body_width = 11, canvas_width = 14, height = 7)
save_root_plot(p, "02_delta_age_by_projection_class", width = 14, height = 7)

treatment_axis <- first_present(c("projection_group", "treatment", "intervention", "genotype", "diet", "stress", "microbiome", "projection_class"), names(pred))
plot_dt <- copy(pred)
plot_dt[, plot_group := get(treatment_axis)]
if (identical(treatment_axis, "projection_class")) plot_dt[, plot_group := projection_class_label]
plot_dt[is.na(plot_group) | plot_group == "", plot_group := projection_class_label]
for (ds in sort(unique(plot_dt$dataset_folder))) {
  dt <- plot_dt[dataset_folder == ds]
  if (nrow(dt) == 0) next
  if (ds == "GSE329024") {
    dt[, plot_group := factor(plot_group, levels = c(
      "2M-8W control", "2M-8W CUMS stress",
      "2M-12W control", "2M-12W CUMS stress",
      "6M-8W control", "6M-8W CUMS stress",
      "6M-12W control", "6M-12W CUMS stress"
    ))]
  }
  p <- ggplot(dt, aes(plot_group, delta_age_week, color = projection_class_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12) +
    geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~model_family, scales = "free_x") +
    labs(
      title = paste0("Raw delta age by ", treatment_axis, ": ", ds),
      x = treatment_axis,
      y = "Delta age (weeks)",
      color = "Projection class"
    ) +
    theme_projection()
  save_plot(p, paste0("04_raw_delta_age_by_dataset_treatment_", safe_name(ds)), width = 10, height = 6.5)
}

if (!is.null(relative) && "relative_delta_age" %in% names(relative) && any(is.finite(relative$relative_delta_age))) {
  rel_plot <- copy(relative)
  rel_plot[, plot_group := get(treatment_axis)]
  if (identical(treatment_axis, "projection_class")) rel_plot[, plot_group := projection_class_label]
  rel_plot[is.na(plot_group) | plot_group == "", plot_group := projection_class_label]
  for (ds in sort(unique(rel_plot$dataset_folder))) {
    dt <- rel_plot[dataset_folder == ds & is.finite(relative_delta_age)]
    if (nrow(dt) == 0) next
    if (ds == "GSE329024") {
      dt[, plot_group := factor(plot_group, levels = c(
        "2M-8W control", "2M-8W CUMS stress",
        "2M-12W control", "2M-12W CUMS stress",
        "6M-8W control", "6M-8W CUMS stress",
        "6M-12W control", "6M-12W CUMS stress"
      ))]
    }
    p <- ggplot(dt, aes(plot_group, relative_delta_age, color = projection_class_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_boxplot(outlier.shape = NA, alpha = 0.12) +
      geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
      facet_wrap(~model_family, scales = "free_x") +
      labs(
        title = paste0("Matched-control relative delta age: ", ds),
        x = treatment_axis,
        y = "Relative delta age (weeks)",
        color = "Projection class"
      ) +
      theme_projection()
    save_plot(p, paste0("05_relative_delta_age_matched_controls_", safe_name(ds)), width = 10, height = 6.5)
  }
} else {
  warn("Matched-control relative delta plot skipped because usable relative_delta_age values were absent.")
}

if (nrow(pred) > 0) {
  pred_root <- copy(pred)
  pred_root[, plot_group := get(treatment_axis)]
  if (identical(treatment_axis, "projection_class")) pred_root[, plot_group := projection_class_label]
  pred_root[is.na(plot_group) | plot_group == "", plot_group := projection_class_label]
  role_shape <- if ("projection_role" %in% names(pred_root)) "projection_role" else "modality"

  p <- ggplot(pred_root, aes(age_week, predicted_age_week, color = projection_class_label, shape = get(role_shape))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(size = 2.3, alpha = 0.86) +
    facet_wrap(~model_family) +
    labs(x = "Chronological age (weeks)", y = "Projected transcriptomic age (weeks)", color = "Projection class", shape = role_shape) +
    theme_projection()
  save_root_plot(p, "external_projection_predicted_vs_age", width = 12, height = 7)

  p <- ggplot(pred_root, aes(plot_group, delta_age_week, color = projection_class_label, shape = get(role_shape))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(position = position_jitter(width = 0.16, height = 0), size = 1.8, alpha = 0.82) +
    facet_grid(model_family ~ dataset_folder, scales = "free_x", space = "free_x") +
    labs(x = treatment_axis, y = "Raw delta age (weeks)", color = "Projection class", shape = role_shape) +
    theme_projection()
  save_root_plot(p, "external_projection_raw_delta_by_dataset_treatment", width = 16, height = 8.5)
  save_root_plot(p, "03_delta_age_by_treatment_within_dataset", width = 16, height = 8.5)
}

if (!is.null(relative) && "relative_delta_age" %in% names(relative) && any(is.finite(relative$relative_delta_age))) {
  rel_root <- copy(relative[is.finite(relative_delta_age)])
  rel_root[, plot_group := get(treatment_axis)]
  if (identical(treatment_axis, "projection_class")) rel_root[, plot_group := projection_class_label]
  rel_root[is.na(plot_group) | plot_group == "", plot_group := projection_class_label]
  role_shape <- if ("projection_role" %in% names(rel_root)) "projection_role" else "modality"

  p <- ggplot(rel_root, aes(dataset_folder, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(position = position_jitter(width = 0.16, height = 0), size = 1.8, alpha = 0.82) +
    facet_wrap(~model_family, scales = "free_x") +
    labs(x = "Dataset", y = "Matched-control relative delta age (weeks)", color = "Projection class", shape = role_shape) +
    theme_projection()
  save_root_plot(p, "external_projection_relative_delta_by_dataset", width = 13, height = 7)

  p <- ggplot(rel_root, aes(projection_class_label, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(position = position_jitter(width = 0.16, height = 0), size = 1.8, alpha = 0.82) +
    facet_wrap(~model_family, scales = "free_x") +
    labs(x = "Projection class", y = "Matched-control relative delta age (weeks)", color = "Projection class", shape = role_shape) +
    theme_projection()
  save_root_plot(p, "external_projection_relative_delta_by_projection_class", width = 13, height = 7)
  save_root_plot(p, "07_matched_control_relative_delta", width = 13, height = 7)

  ds <- rel_root[dataset_folder == "GSE279446"]
  if (nrow(ds) > 0) {
    p <- ggplot(ds, aes(projection_class_label, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_point(position = position_jitter(width = 0.16, height = 0), size = 2, alpha = 0.86) +
      facet_grid(model_family ~ gse279446_cohort, scales = "free_x", space = "free_x") +
      labs(x = NULL, y = "Relative delta age (weeks)", color = "Projection class", shape = role_shape) +
      theme_projection()
    save_root_plot(p, "GSE279446_cohort_relative_delta", width = 11, height = 7)
  }

  ds <- rel_root[dataset_folder == "GSE279573"]
  if (nrow(ds) > 0) {
    p <- ggplot(ds, aes(plot_group, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.2, alpha = 0.9) +
      facet_wrap(~model_family) +
      labs(title = "GSE279573 relative delta age (n_controls = 2)", x = NULL, y = "Relative delta age (weeks)", color = "Projection class", shape = role_shape) +
      theme_projection()
    save_root_plot(p, "GSE279573_relative_delta", width = 10, height = 6.5)
  }

  ds <- rel_root[dataset_folder %chin% c("GSE194196", "GSE252791")]
  if (nrow(ds) > 0) {
    p <- ggplot(ds, aes(plot_group, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.88) +
      facet_grid(model_family ~ dataset_folder, scales = "free_x", space = "free_x") +
      labs(x = NULL, y = "Relative delta age (weeks)", color = "Projection class", shape = role_shape) +
      theme_projection()
    save_root_plot(p, "GSE194196_GSE252791_Bnc1_relative_delta", width = 11, height = 7)
  }

  ds <- rel_root[dataset_folder == "GSE240547"]
  if (nrow(ds) > 0) {
    p <- ggplot(ds, aes(plot_group, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.88) +
      facet_wrap(~model_family) +
      labs(x = NULL, y = "Relative delta age (weeks)", color = "Projection class", shape = role_shape) +
      theme_projection()
    save_root_plot(p, "GSE240547_HFHS_relative_delta", width = 9, height = 6)
  }

  ds <- rel_root[dataset_folder == "PRJNA1076509"]
  if (nrow(ds) > 0) {
    p <- ggplot(ds, aes(plot_group, relative_delta_age, color = projection_class_label, shape = get(role_shape))) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_point(position = position_jitter(width = 0.12, height = 0), size = 2.1, alpha = 0.88) +
      facet_wrap(~model_family) +
      labs(x = NULL, y = "Relative delta age (weeks)", color = "Projection class", shape = role_shape) +
      theme_projection()
    save_root_plot(p, "PRJNA1076509_FMT_relative_delta", width = 10, height = 6.5)
  }
}

discordant_labels <- model_concordance[order(-absolute_disagreement)][seq_len(min(.N, 8))]
p <- ggplot(model_concordance, aes(delta_age_elasticnet, delta_age_pca_glmnet, color = projection_class_label, shape = modality)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey35") +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(
    title = "ElasticNet versus PCA/glmnet delta-age concordance",
    x = "ElasticNet delta age (weeks)",
    y = "PCA/glmnet delta age (weeks)",
    color = "Projection class",
    shape = "Modality"
  ) +
  theme_projection()
if (nrow(discordant_labels) > 0) {
  if (has_ggrepel) {
    p <- p + ggrepel::geom_text_repel(data = discordant_labels, aes(label = sample_key), size = 2.4, max.overlaps = 20, show.legend = FALSE)
  } else {
    p <- p + geom_text(data = discordant_labels, aes(label = sample_key), size = 2.2, vjust = -0.7, show.legend = FALSE)
  }
}
save_plot(p, "06_clock_concordance_delta_age", width = 9, height = 7)
save_root_plot(p, "external_projection_model_concordance", width = 9, height = 7)

p <- ggplot(model_concordance, aes(projection_class_label, absolute_disagreement, color = projection_class_label)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.8) +
  facet_wrap(~dataset_folder, scales = "free_x") +
  labs(
    title = "Model disagreement by group",
    x = "Projection class",
    y = "Absolute PCA - ElasticNet delta-age difference (weeks)",
    color = "Projection class"
  ) +
  theme_projection()
save_plot(p, "07_model_disagreement_by_group", width = 13, height = 10)

if ("gene_coverage_fraction" %in% names(pred) && any(is.finite(pred$gene_coverage_fraction))) {
  p <- ggplot(pred, aes(projection_class_label, gene_coverage_fraction, color = model_family)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12) +
    geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~dataset_folder, scales = "free_x") +
    labs(
      title = "Projection gene coverage",
      x = "Projection class",
      y = "Gene coverage fraction",
      color = "Model family"
  ) +
    theme_projection()
  save_plot(p, "08_projection_qc_coverage_ood_gene_coverage", width = 13, height = 8)
  save_root_plot(p, "external_projection_gene_coverage", width = 13, height = 8)
}
if ("projection_distance" %in% names(pred) && any(is.finite(pred$projection_distance))) {
  p <- ggplot(pred, aes(projection_class_label, projection_distance, color = OOD_flag)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.12) +
    geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~model_family, scales = "free_y") +
    labs(
      title = "Projection OOD distance",
      x = "Projection class",
      y = "Projection distance",
      color = "OOD flag"
    ) +
    theme_projection()
  save_plot(p, "08_projection_qc_coverage_ood_distance", width = 11, height = 7)
}
if (all(c("gene_coverage_fraction", "projection_distance") %in% names(pred)) &&
  any(is.finite(pred$gene_coverage_fraction)) && any(is.finite(pred$projection_distance))) {
  p <- ggplot(pred, aes(gene_coverage_fraction, projection_distance, color = OOD_flag, shape = low_coverage_flag)) +
    geom_vline(xintercept = 0.80, linetype = "dashed", color = "grey45") +
    geom_point(alpha = 0.85, size = 2.2) +
    facet_wrap(~model_family, scales = "free_y") +
    labs(
      title = "Projection coverage and OOD diagnostic scatter",
      x = "Gene coverage fraction",
      y = "Projection distance",
      color = "OOD flag",
      shape = "Low coverage"
    ) +
    theme_projection()
  save_plot(p, "08_projection_qc_coverage_ood_scatter", width = 10, height = 7)
  save_root_plot(p, "06_ood_gene_coverage", width = 10, height = 7)
}

if (!is.null(pca_scores) && all(c("PC1", "PC2") %in% names(pca_scores))) {
  pca_plot <- merge(
    pca_scores,
    unique(pred[, .(sample_key, projection_class, projection_class_label, modality)]),
    by = "sample_key",
    all.x = TRUE
  )
  for (cc in c("projection_class", "projection_class_label", "modality")) pca_plot <- resolve_merged_col(pca_plot, cc)
  pca_plot <- add_projection_class_label(pca_plot)
  p <- ggplot(pca_plot, aes(PC1, PC2, color = projection_class_label, shape = modality)) +
    geom_hline(yintercept = 0, color = "grey85") +
    geom_vline(xintercept = 0, color = "grey85") +
    geom_point(size = 2.5, alpha = 0.88) +
    labs(
      title = "PCA PC1/PC2 projection map",
      subtitle = "External samples only.",
      x = "Frozen deployment PC1 score",
      y = "Frozen deployment PC2 score",
      color = "Projection class",
      shape = "Modality"
    ) +
    theme_projection()
  save_plot_fixed_body_wide_legend(p, "09_pca_projection_map", body_width = 10.5, canvas_width = 14, height = 8)
  save_root_plot_fixed_body_wide_legend(p, "04_pca_pc1_pc2_projection_map", body_width = 10.5, canvas_width = 14, height = 8)
} else {
  warn("PCA PC1/PC2 projection map skipped because pca_projection_scores.tsv is absent or lacks PC1/PC2.")
}

if (!is.null(pca_scores)) {
  contrib_cols <- grep("_contribution_week$", names(pca_scores), value = TRUE)
  if (length(contrib_cols) > 0) {
    pca_meta <- merge(
      pca_scores,
      unique(pred[, .(sample_key, dataset_folder, projection_class, projection_class_label, modality)]),
      by = "sample_key",
      all.x = TRUE
    )
    for (cc in c("dataset_folder", "projection_class", "projection_class_label", "modality")) pca_meta <- resolve_merged_col(pca_meta, cc)
    pca_meta <- add_projection_class_label(pca_meta)
    pca_meta[, row_group := paste(dataset_folder, projection_class_label, sep = " | ")]
    heat <- melt(
      pca_meta,
      id.vars = c("row_group", "dataset_folder", "projection_class"),
      measure.vars = contrib_cols,
      variable.name = "PC",
      value.name = "contribution_week"
    )
    heat[, PC := gsub("_contribution_week$", "", PC)]
    heat_summary <- heat[, .(mean_contribution_week = mean(contribution_week, na.rm = TRUE)), by = .(row_group, dataset_folder, projection_class, PC)]
    p <- ggplot(heat_summary, aes(PC, row_group, fill = mean_contribution_week)) +
      geom_tile(color = "white", linewidth = 0.2) +
      scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44E52", midpoint = 0) +
      labs(
        title = "PCA selected-PC contribution heatmap",
        subtitle = "Rows are dataset/projection-class means; values are calibrated week contributions.",
        x = "Selected PC",
        y = NULL,
        fill = "Weeks"
      ) +
      theme_projection() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text.y = element_text(size = 9))
    save_plot(p, "10_pca_pc_contribution_heatmap", width = 9, height = 8)
    save_root_plot(p, "pca_projection_pc_contribution_heatmap", width = 9, height = 8)
    save_root_plot(p, "05_pca_selected_pc_contribution_heatmap", width = 9, height = 8)

    pca_contribution_summary <- dcast(
      heat_summary,
      row_group + dataset_folder + projection_class ~ PC,
      value.var = "mean_contribution_week"
    )
    write_tsv(pca_contribution_summary, "pca_pc_contribution_group_summary.tsv")
  } else {
    warn("PCA contribution heatmap skipped because no *_contribution_week columns were found.")
    write_tsv(data.table(note = "No *_contribution_week columns found in pca_projection_scores.tsv."), "pca_pc_contribution_group_summary.tsv")
  }
} else {
  warn("PCA contribution heatmap skipped because pca_projection_scores.tsv is absent.")
  write_tsv(data.table(note = "pca_projection_scores.tsv absent."), "pca_pc_contribution_group_summary.tsv")
}

grade_counts <- projection_claim_grading[, .(n_groups = .N), by = claim_grade]
p <- ggplot(grade_counts, aes(claim_grade, n_groups, fill = claim_grade)) +
  geom_col(color = "grey25") +
  geom_text(aes(label = n_groups), vjust = -0.25, size = 4) +
  labs(
    title = "Projection claim grading",
    x = "Claim grade",
    y = "Number of groups",
    fill = "Claim grade"
  ) +
  theme_projection() +
  theme(legend.position = "none")
save_plot(p, "11_projection_claim_grading", width = 10, height = 6)

worst_plot <- worst_discrepant_samples[
  discrepancy_category != "skipped_unaligned_sample"
][order(discrepancy_category, -ranking_value)][seq_len(min(.N, 30))]
if (nrow(worst_plot) > 0) {
  worst_plot[, row_label := paste(discrepancy_category, sample_key, sep = " | ")]
  worst_plot[, row_label := make.unique(row_label, sep = " #")]
  worst_plot[, row_label := factor(row_label, levels = rev(row_label))]
  p <- ggplot(worst_plot, aes(ranking_value, row_label, color = discrepancy_category)) +
    geom_point(size = 2.5) +
    labs(
      title = "Worst/discrepant projected samples",
      subtitle = "Full table, including skipped samples, is written to worst_discrepant_samples.tsv.",
      x = "Ranking value",
      y = NULL,
      color = "Category"
    ) +
    theme_projection() +
    theme(axis.text.y = element_text(size = 7))
  save_plot(p, "12_worst_discrepant_samples", width = 12, height = 8)
}

coverage_summary <- pred[, .(
  n_prediction_rows = .N,
  mean_gene_coverage = mean(gene_coverage_fraction, na.rm = TRUE),
  min_gene_coverage = min(gene_coverage_fraction, na.rm = TRUE),
  n_low_coverage = sum(low_coverage_flag, na.rm = TRUE),
  n_OOD = sum(OOD_flag_bool, na.rm = TRUE),
  OOD_fraction = mean(OOD_flag_bool, na.rm = TRUE)
), by = model_family]
write_tsv(coverage_summary, "projection_qc_coverage_ood_summary.tsv")

raw_delta_by_clock <- pred[, .(
  n = .N,
  n_samples = uniqueN(sample_key),
  mean_delta_age_week = mean(delta_age_week, na.rm = TRUE),
  median_delta_age_week = median(delta_age_week, na.rm = TRUE),
  sd_delta_age_week = sd(delta_age_week, na.rm = TRUE),
  min_delta_age_week = min(delta_age_week, na.rm = TRUE),
  max_delta_age_week = max(delta_age_week, na.rm = TRUE)
), by = .(model_family, selector)]
write_tsv(raw_delta_by_clock, "raw_delta_age_summary_by_clock.tsv")

concordance_overall <- model_concordance[, .(
  n = .N,
  correlation = if (.N >= 3) suppressWarnings(cor(delta_age_elasticnet, delta_age_pca_glmnet, use = "complete.obs")) else NA_real_,
  fraction_same_direction = mean(same_direction_flag, na.rm = TRUE),
  mean_absolute_disagreement = mean(absolute_disagreement, na.rm = TRUE),
  median_absolute_disagreement = median(absolute_disagreement, na.rm = TRUE)
)]
write_tsv(concordance_overall, "model_concordance_overall_summary.tsv")

main_text_groups <- projection_claim_grading[claim_grade %chin% c("strong_projection_signal", "moderate_projection_signal")]
supplementary_groups <- projection_claim_grading[claim_grade %chin% c("descriptive_only", "exclude_from_claims")]

fmt_table <- function(dt, cols, n = Inf) {
  if (nrow(dt) == 0) return("_None._")
  cols <- cols[cols %in% names(dt)]
  x <- head(dt[, ..cols], n)
  header <- paste(cols, collapse = " | ")
  sep <- paste(rep("---", length(cols)), collapse = " | ")
  rows <- apply(x, 1, function(r) paste(r, collapse = " | "))
  paste(c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |")), collapse = "\n")
}

summary_lines <- c(
  "# External projection diagnostic summary",
  "",
  "## Run summary",
  paste0("- Prediction rows, all: ", n_prediction_rows_all, "."),
  paste0("- Prediction rows, primary: ", nrow(pred), "."),
  paste0("- Unique projected samples, all: ", n_projected_all, "."),
  paste0("- Unique projected samples, primary: ", n_projected, "."),
  paste0("- Retained clocks: ", paste(unique(paste(pred$model_family, pred$selector, sep = " / ")), collapse = "; "), "."),
  paste0("- Skipped unaligned samples: ", uniqueN(skipped$sample_key), "."),
  paste0("- Diagnostic warnings: ", length(warnings_log), "."),
  "",
  "## Artifact readiness and calibration",
  fmt_table(deploy_cal, c("model_family", "selector", "intercept", "slope", "calibration_n", "intercept_matches_expected", "slope_matches_expected")),
  "",
  "## Projected sample counts: all vs primary",
  fmt_table(validation_checks, c("check", "passed", "detail")),
  "",
  "## Control matching logic",
  "- Standard datasets use same-dataset control_meta pools when available.",
  "- GSE240547 HFHS diet samples use same-study control_diet projection-internal controls.",
  "- PRJNA1076509 FMT samples use young female donor microbiota samples as the explicit control pool.",
  "- PRJNA863443 VCD samples use CTL controls matched by injection_age/post_injection_days; Foxl2+/- samples use Foxl2+/+ controls matched by age_group; natural ageing rows are placement-only references.",
  "- GSE228975 is retained only as OOD/supplementary all-projection output with no matched-control relative delta.",
  "",
  "## Dataset control-status table",
  fmt_table(dataset_relative_control_status, c("dataset_folder", "expected_control_strategy", "observed_control_strategy", "n_projected_samples_all", "n_projected_samples_primary", "n_matched_control_rows", "n_unmatched_rows", "n_invalid_rows", "pass_fail"), 20),
  "",
  "## Special section: GSE240547",
  fmt_table(special_dataset_control_audit[dataset_folder == "GSE240547"], c("check_name", "observed", "expected", "pass_fail")),
  "",
  "## Special section: PRJNA1076509",
  fmt_table(special_dataset_control_audit[dataset_folder == "PRJNA1076509"], c("check_name", "observed", "expected", "pass_fail")),
  "",
  "## Special section: PRJNA863443",
  fmt_table(special_dataset_control_audit[dataset_folder == "PRJNA863443"], c("check_name", "observed", "expected", "pass_fail")),
  "",
  "## Special section: GSE228975",
  fmt_table(special_dataset_control_audit[dataset_folder == "GSE228975"], c("check_name", "observed", "expected", "pass_fail")),
  "",
  "## Matched-control integrity audit",
  fmt_table(matched_control_integrity_audit, c("dataset_folder", "model_family", "n_rows", "n_invalid_no_matched_control", "n_non_na_relative_delta", "n_primary_excluded_with_relative_delta", "pass_fail"), 30),
  "",
  "## Projection classes represented",
  paste0("- ", paste(sort(unique(pred$projection_class)), collapse = ", "), "."),
  "",
  "## Gene coverage and OOD summary",
  fmt_table(coverage_summary, c("model_family", "n_prediction_rows", "mean_gene_coverage", "min_gene_coverage", "n_low_coverage", "n_OOD", "OOD_fraction")),
  "",
  "## Raw delta-age summary by clock",
  fmt_table(raw_delta_by_clock, c("model_family", "selector", "n_samples", "mean_delta_age_week", "median_delta_age_week", "sd_delta_age_week")),
  "",
  "## Matched-control relative delta-age summary",
  if (!is.null(relative)) {
    fmt_table(matched_control_relative_delta_summary, c("model_family", "selector", "dataset_folder", "projection_class", "n", "mean_relative_delta_age", "median_relative_delta_age"), 12)
  } else {
    "- Matched-control file was absent; relative summaries were skipped."
  },
  "",
  "## ElasticNet/PCA concordance",
  fmt_table(concordance_overall, c("n", "correlation", "fraction_same_direction", "mean_absolute_disagreement", "median_absolute_disagreement")),
  "",
  "## PCA-axis contribution summary",
  if (!is.null(pca_scores) && length(grep("_contribution_week$", names(pca_scores), value = TRUE)) > 0) {
    "- PCA selected-PC contribution summaries were generated as group means in `pca_pc_contribution_group_summary.tsv` and `10_pca_pc_contribution_heatmap.*`."
  } else {
    "- PCA selected-PC contribution columns were unavailable."
  },
  "",
  "## Claim grading table",
  fmt_table(projection_claim_grading, c("dataset_folder", "modality", "projection_class", active_group_metadata_cols, "n", "OOD_fraction", "fraction_same_direction", "claim_grade"), 30),
  "",
  "## Groups safe for main-text interpretation",
  fmt_table(main_text_groups, c("dataset_folder", "modality", "projection_class", active_group_metadata_cols, "n", "claim_grade"), 30),
  "",
  "## Supplementary/descriptive or excluded groups",
  fmt_table(supplementary_groups, c("dataset_folder", "modality", "projection_class", active_group_metadata_cols, "n", "claim_grade", "grading_reason"), 30),
  "",
  "## Explicit caution",
  "External projections position samples along the trained transcriptomic ageing reference. They do not prove rejuvenation, accelerated ageing, or causality without independent functional evidence.",
  "",
  "## Warnings",
  if (length(warnings_log) > 0) paste0("- ", warnings_log) else "- None."
)
summary_path <- file.path(out_dir, "external_projection_diagnostic_summary.md")
writeLines(summary_lines, summary_path)

if (length(warnings_log) > 0) {
  write_tsv(data.table(warning = warnings_log), "external_projection_diagnostic_warnings.tsv")
} else {
  write_tsv(data.table(warning = character()), "external_projection_diagnostic_warnings.tsv")
}

cat("Output directory:", out_dir, "\n")
cat("TSV files written:", length(unique(tables_written)), "\n")
cat("Plots written:", length(unique(plot_prefixes_written)), "\n")
cat("Summary markdown:", summary_path, "\n")
