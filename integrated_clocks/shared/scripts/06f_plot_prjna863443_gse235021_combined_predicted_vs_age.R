#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0L) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
clock_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

results_dir <- file.path(clock_dir, "results")
projection_dir <- file.path(results_dir, "external_projection")
output_dir <- file.path(projection_dir, "GSE235021")
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

prediction_path <- file.path(projection_dir, "external_projection_predictions_all.tsv")
alignment_path <- file.path(results_dir, "integrated_bulk_pseudobulk_sample_expression_alignment.tsv")
artifact_paths <- c(
  elasticnet = file.path(
    results_dir, "model_artifacts", "deployment",
    "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds"
  ),
  pca_glmnet = file.path(
    results_dir, "model_artifacts", "deployment",
    "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
  )
)
for (path in c(prediction_path, alignment_path, artifact_paths)) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}

log_path <- file.path(log_dir, "PRJNA863443_GSE235021_combined_predicted_vs_age.log")
log_lines <- character()
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "  ", paste0(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
}

read_expression_matrix <- function(path) {
  dt <- fread(path)
  gene_col <- names(dt)[1]
  genes <- as.character(dt[[gene_col]])
  dt[, (gene_col) := NULL]
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- genes
  keep <- !is.na(genes) & genes != "" & !duplicated(genes)
  mat[keep, , drop = FALSE]
}

load_reference_expression <- function(reference_meta) {
  alignment <- fread(alignment_path, sep = "\t")
  pieces <- list()
  for (path in unique(reference_meta$expression_logcpm_path)) {
    path_meta <- reference_meta[expression_logcpm_path == path]
    mat <- read_expression_matrix(path)
    expression_cols <- alignment$expression_column[match(path_meta$sample_key, alignment$sample_key)]
    missing <- is.na(expression_cols) | !expression_cols %in% colnames(mat)
    if (any(missing)) {
      stop(
        "Reference-expression alignment is missing samples in ", path, ": ",
        paste(path_meta$sample_key[missing], collapse = ", "),
        call. = FALSE
      )
    }
    mat <- mat[, expression_cols, drop = FALSE]
    colnames(mat) <- path_meta$sample_key
    pieces[[path]] <- mat
  }
  pieces
}

build_scaled_matrix <- function(expr_pieces, reference_meta, genes, center, scale) {
  out <- matrix(
    NA_real_,
    nrow = nrow(reference_meta),
    ncol = length(genes),
    dimnames = list(reference_meta$sample_key, genes)
  )
  for (mat in expr_pieces) {
    keys <- intersect(colnames(mat), reference_meta$sample_key)
    present_genes <- intersect(genes, rownames(mat))
    if (length(keys) > 0L && length(present_genes) > 0L) {
      out[keys, present_genes] <- t(mat[present_genes, keys, drop = FALSE])
    }
  }
  for (gene in genes) out[is.na(out[, gene]), gene] <- center[[gene]]
  scaled <- sweep(out, 2, center[genes], "-")
  scale <- scale[genes]
  scale[!is.finite(scale) | scale == 0] <- 1
  scaled <- sweep(scaled, 2, scale, "/")
  scaled[!is.finite(scaled)] <- 0
  scaled
}

predict_reference <- function(artifact, model_family, expr_pieces) {
  reference_meta <- as.data.table(artifact$reference_metadata)
  genes <- as.character(artifact$selected_genes)
  x <- build_scaled_matrix(
    expr_pieces,
    reference_meta,
    genes,
    artifact$gene_center,
    artifact$gene_scale
  )
  if (model_family == "elasticnet") {
    raw <- as.numeric(predict(artifact$model$glmnet_fit, newx = x))
    selector <- "hvg_baseline"
  } else {
    pc_count <- as.integer(artifact$model$selected_pc_count)
    pcs <- predict(artifact$model$pca_fit, newdata = x)[, seq_len(pc_count), drop = FALSE]
    raw <- as.numeric(predict(artifact$model$glmnet_fit, newx = pcs))
    selector <- "version_b_age_dataset_penalized"
  }
  calibrated <- artifact$calibration$intercept + artifact$calibration$slope * raw
  data.table(
    sample_key = reference_meta$sample_key,
    age_week = as.numeric(reference_meta$age_week),
    predicted_age_week = calibrated,
    model_family = model_family,
    selector = selector
  )
}

theme_combined <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

log_msg("Reading the two external cohorts.")
pred_all <- fread(prediction_path, sep = "\t", na.strings = c("", "NA", "NaN"))
combined <- pred_all[dataset_folder %chin% c("PRJNA863443", "GSE235021")]
combined <- combined[, .(
  sample_key,
  dataset_folder,
  modality,
  age_week,
  treatment,
  genotype,
  sample_title,
  projection_class,
  model_family,
  selector,
  predicted_age_week,
  delta_age_week
)]

if (nrow(combined) != 136L ||
    uniqueN(combined$sample_key) != 68L ||
    anyDuplicated(combined[, .(sample_key, model_family)])) {
  stop(
    "Expected 68 samples × two clocks = 136 unique sample-clock rows; observed ",
    uniqueN(combined$sample_key), " samples and ", nrow(combined), " rows.",
    call. = FALSE
  )
}
if (any(!is.finite(combined$age_week)) || any(!is.finite(combined$predicted_age_week))) {
  stop("Chronological or predicted ages contain non-finite values.", call. = FALSE)
}

combined[, dataset_label := factor(
  dataset_folder,
  levels = c("PRJNA863443", "GSE235021"),
  labels = c("PRJNA863443", "GSE235021")
)]
combined[, modality_label := factor(
  modality,
  levels = c("bulk_rnaseq", "scrna_pseudobulk"),
  labels = c("Bulk RNA-seq", "scRNA pseudobulk")
)]

sample_audit <- unique(combined[, .(sample_key, dataset_folder, modality)])
sample_counts <- sample_audit[, .(n_samples = .N), by = .(dataset_folder, modality)]

expected_counts <- data.table(
  dataset_folder = c("PRJNA863443", "GSE235021", "GSE235021"),
  modality = c("scrna_pseudobulk", "bulk_rnaseq", "scrna_pseudobulk"),
  expected_n = c(33L, 18L, 17L)
)
count_audit <- merge(expected_counts, sample_counts, by = c("dataset_folder", "modality"), all.x = TRUE)
count_audit[is.na(n_samples), n_samples := 0L]
if (any(count_audit$n_samples != count_audit$expected_n)) {
  stop("Combined cohort modality counts do not match 33 + 18 + 17 samples.", call. = FALSE)
}

log_msg("Reconstructing the unchanged 117-sample frozen deployment reference fits.")
artifact_list <- lapply(artifact_paths, readRDS)
reference_meta <- as.data.table(artifact_list[[1]]$reference_metadata)
if (nrow(reference_meta) != 117L || uniqueN(reference_meta$sample_key) != 117L) {
  stop("Expected 117 unique frozen deployment reference samples.", call. = FALSE)
}
expr_pieces <- load_reference_expression(reference_meta)
locked_line <- rbindlist(list(
  predict_reference(artifact_list$elasticnet, "elasticnet", expr_pieces),
  predict_reference(artifact_list$pca_glmnet, "pca_glmnet", expr_pieces)
))
locked_stats <- locked_line[, .(
  n_reference_samples = .N,
  locked_mae_week = mean(abs(predicted_age_week - age_week)),
  locked_r2 = summary(lm(predicted_age_week ~ age_week))$r.squared
), by = .(model_family, selector)]
locked_stats[, locked_metric_label := sprintf(
  "Deployment training fit (in-sample)\nMAE = %.1f weeks\nR² = %.2f",
  locked_mae_week,
  locked_r2
)]

fwrite(
  combined[order(model_family, dataset_folder, modality, age_week, sample_key)],
  file.path(table_dir, "PRJNA863443_GSE235021_combined_predicted_vs_age_source.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  locked_line,
  file.path(table_dir, "PRJNA863443_GSE235021_frozen_deployment_reference_fit.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  locked_stats,
  file.path(table_dir, "PRJNA863443_GSE235021_frozen_deployment_reference_fit_stats.tsv"),
  sep = "\t",
  na = "NA"
)

p <- ggplot(
  combined,
  aes(
    age_week,
    predicted_age_week,
    color = dataset_label,
    shape = modality_label
  )
) +
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
  geom_point(size = 2.5, alpha = 0.82, stroke = 0.8) +
  geom_text(
    data = locked_stats,
    aes(x = -Inf, y = Inf, label = locked_metric_label),
    inherit.aes = FALSE,
    hjust = -0.08,
    vjust = 1.08,
    size = 3.8,
    lineheight = 1.05
  ) +
  facet_wrap(~model_family) +
  scale_color_manual(values = c(PRJNA863443 = "#0072B2", GSE235021 = "#D55E00")) +
  scale_shape_manual(values = c("Bulk RNA-seq" = 17, "scRNA pseudobulk" = 16)) +
  labs(
    x = "Chronological age (weeks)",
    y = "Predicted age (weeks)",
    color = "Dataset",
    shape = "Modality"
  ) +
  theme_combined()

stem <- "01_prjna863443_gse235021_combined_predicted_vs_age"
pdf_path <- file.path(plot_dir, paste0(stem, ".pdf"))
png_path <- file.path(plot_dir, paste0(stem, ".png"))
ggsave(pdf_path, p, width = 12, height = 7, limitsize = FALSE)
ggsave(png_path, p, width = 12, height = 7, dpi = 300, limitsize = FALSE)

validation <- rbindlist(list(
  data.table(
    check = "combined_unique_samples",
    expected = "68",
    observed = as.character(uniqueN(combined$sample_key)),
    status = ifelse(uniqueN(combined$sample_key) == 68L, "pass", "fail")
  ),
  data.table(
    check = "combined_sample_clock_rows",
    expected = "136",
    observed = as.character(nrow(combined)),
    status = ifelse(nrow(combined) == 136L, "pass", "fail")
  ),
  data.table(
    check = "PRJNA863443_scRNA_pseudobulk_samples",
    expected = "33",
    observed = as.character(sample_audit[dataset_folder == "PRJNA863443", .N]),
    status = ifelse(sample_audit[dataset_folder == "PRJNA863443", .N] == 33L, "pass", "fail")
  ),
  data.table(
    check = "GSE235021_bulk_samples",
    expected = "18",
    observed = as.character(sample_audit[dataset_folder == "GSE235021" & modality == "bulk_rnaseq", .N]),
    status = ifelse(sample_audit[dataset_folder == "GSE235021" & modality == "bulk_rnaseq", .N] == 18L, "pass", "fail")
  ),
  data.table(
    check = "GSE235021_scRNA_pseudobulk_samples",
    expected = "17",
    observed = as.character(sample_audit[dataset_folder == "GSE235021" & modality == "scrna_pseudobulk", .N]),
    status = ifelse(sample_audit[dataset_folder == "GSE235021" & modality == "scrna_pseudobulk", .N] == 17L, "pass", "fail")
  ),
  data.table(
    check = "finite_external_predictions",
    expected = "136",
    observed = as.character(sum(is.finite(combined$predicted_age_week))),
    status = ifelse(all(is.finite(combined$predicted_age_week)), "pass", "fail")
  ),
  data.table(
    check = "frozen_reference_samples_per_clock",
    expected = "117;117",
    observed = paste(locked_stats$n_reference_samples, collapse = ";"),
    status = ifelse(all(locked_stats$n_reference_samples == 117L), "pass", "fail")
  ),
  data.table(
    check = "reference_metrics_match_original_plot",
    expected = "elasticnet MAE=4.3 R2=0.94; pca_glmnet MAE=9.8 R2=0.64",
    observed = paste(
      sprintf(
        "%s MAE=%.1f R2=%.2f",
        locked_stats$model_family,
        locked_stats$locked_mae_week,
        locked_stats$locked_r2
      ),
      collapse = "; "
    ),
    status = ifelse(
      all(round(locked_stats$locked_mae_week, 1) == c(4.3, 9.8)) &&
        all(round(locked_stats$locked_r2, 2) == c(0.94, 0.64)),
      "pass",
      "fail"
    )
  ),
  data.table(
    check = "nonempty_PDF_PNG_exports",
    expected = "2",
    observed = as.character(sum(file.exists(c(pdf_path, png_path)) & file.info(c(pdf_path, png_path))$size > 0)),
    status = ifelse(
      all(file.exists(c(pdf_path, png_path))) &&
        all(file.info(c(pdf_path, png_path))$size > 0),
      "pass",
      "fail"
    )
  )
))
fwrite(
  validation,
  file.path(table_dir, "PRJNA863443_GSE235021_combined_predicted_vs_age_validation.tsv"),
  sep = "\t"
)
if (any(validation$status != "pass")) {
  stop("One or more combined predicted-versus-age validation checks failed.", call. = FALSE)
}

readme <- c(
  "# PRJNA863443 + GSE235021 combined predicted-versus-age projection",
  "",
  "This plot follows the PRJNA863443 predicted-versus-age format and combines",
  "33 PRJNA863443 scRNA-derived pseudobulk samples with 18 GSE235021 bulk RNA-seq",
  "and 17 GSE235021 scRNA-derived pseudobulk samples.",
  "",
  "The dashed line is identity. The solid black line and MAE/R² labels are reconstructed",
  "from the unchanged frozen 117-sample deployment reference for each clock; no external",
  "sample is used to fit that line. Dataset is encoded by color and modality by shape.",
  "",
  "The figure is an external projection diagnostic, not an external accuracy estimate,",
  "because diet, genotype, VCD treatment, and modality can shift predicted age."
)
writeLines(
  readme,
  file.path(output_dir, "README_PRJNA863443_GSE235021_combined_predicted_vs_age.md")
)

log_msg(
  "Completed combined predicted-versus-age figure: 68 samples × two clocks; ",
  "33 PRJNA863443 pseudobulk, 18 GSE235021 bulk, 17 GSE235021 pseudobulk; ",
  sum(validation$status == "pass"), "/", nrow(validation), " validation checks passed."
)
writeLines(log_lines, log_path)
