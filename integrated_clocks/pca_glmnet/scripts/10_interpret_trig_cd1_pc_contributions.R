#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
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
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
external_dir <- file.path(project_dir, "results", "external_test", "Trig_CD1_mice")
out_dir <- file.path(external_dir, "pc_interpretation")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

warnings_vec <- character()
add_warning <- function(x) warnings_vec <<- unique(c(warnings_vec, x))

selected_pcs <- paste0("PC", 1:7)
contrib_path <- file.path(external_dir, "trig_cd1_pca_pc_contributions.tsv")
score_path <- file.path(external_dir, "trig_cd1_pca_projection_scores.tsv")
pred_path <- file.path(external_dir, "trig_cd1_external_projection_predictions.tsv")
group_summary_path <- file.path(external_dir, "trig_cd1_external_projection_group_summary.tsv")
concordance_path <- file.path(external_dir, "trig_cd1_model_concordance.tsv")
artifact_path <- file.path(
  project_dir,
  "results", "model_artifacts", "deployment",
  "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
)

required_inputs <- c(contrib_path, score_path, pred_path, group_summary_path, concordance_path, artifact_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop("Missing required input files: ", paste(missing_inputs, collapse = "; "), call. = FALSE)
}

contrib_raw <- fread(contrib_path)
scores_raw <- fread(score_path)
pred <- fread(pred_path)
group_summary <- fread(group_summary_path)
concordance <- fread(concordance_path)
artifact <- readRDS(artifact_path)

artifact_pcs <- artifact$model$selected_pc_names %||% paste0("PC", seq_len(as.integer(artifact$model$selected_pc_count)))
if (!identical(as.character(artifact_pcs), selected_pcs)) {
  stop("Deployment artifact selected PCs are not exactly PC1-PC7: ", paste(artifact_pcs, collapse = ", "), call. = FALSE)
}
if (!identical(artifact$training_scope, "integrated_bulk_pseudobulk")) {
  stop("Unexpected deployment artifact training_scope: ", artifact$training_scope, call. = FALSE)
}
if (!identical(artifact$model_family, "pca_glmnet") || !identical(artifact$selector, "version_b_age_dataset_penalized")) {
  stop("Deployment artifact is not the final PCA/glmnet Version B clock.", call. = FALSE)
}
if ("PC8" %in% contrib_raw$pc || "PC8" %in% names(scores_raw)) {
  stop("PC8 appears in the Trig CD1 external PC projection tables.", call. = FALSE)
}
if (!setequal(unique(contrib_raw$pc), selected_pcs)) {
  stop("PCA contribution table must contain exactly PC1-PC7.", call. = FALSE)
}
if (!all(selected_pcs %in% names(scores_raw))) {
  stop("PCA score table must contain PC1-PC7 score columns.", call. = FALSE)
}
if (any(is.na(contrib_raw$treatment_group))) stop("A PCA contribution row is missing treatment_group.", call. = FALSE)

group_counts <- unique(contrib_raw[, .(sample_key, treatment_group)])[, .N, by = treatment_group]
expected_counts <- data.table(
  treatment_group = c("young_3mo_control", "ama_9mo_control", "trg_9mo"),
  expected_n = c(4L, 4L, 4L)
)
group_counts <- merge(expected_counts, group_counts, by = "treatment_group", all.x = TRUE)
group_counts[is.na(N), N := 0L]
bad_counts <- group_counts[N != expected_n]
if (nrow(bad_counts) > 0) {
  stop(
    "Unexpected group sample counts: ",
    paste(paste0(bad_counts$treatment_group, "=", bad_counts$N, " expected ", bad_counts$expected_n), collapse = "; "),
    call. = FALSE
  )
}

pca_pred <- pred[model_family == "pca_glmnet" & selector == "version_b_age_dataset_penalized"]
if (nrow(pca_pred) == 0) stop("PCA/glmnet projection predictions are missing.", call. = FALSE)

meta_cols <- c(
  "sample_key", "sample_id", "projection_class", "age_week", "model_family", "selector",
  "predicted_age_week", "delta_age_week", "relative_delta_age_week"
)
missing_pred_cols <- setdiff(meta_cols, names(pca_pred))
if (length(missing_pred_cols) > 0) {
  stop("PCA prediction table is missing columns: ", paste(missing_pred_cols, collapse = ", "), call. = FALSE)
}

pc_long <- merge(
  contrib_raw[, .(
    sample_key, sample_column, treatment_group, age_month, pc, pc_score,
    raw_pc_coefficient, calibrated_pc_coefficient,
    calibrated_pc_contribution_week = calibrated_PC_contribution_week
  )],
  pca_pred[, ..meta_cols],
  by = "sample_key",
  all.x = TRUE,
  sort = FALSE
)
if (any(is.na(pc_long$projection_class))) stop("Projection metadata failed to join onto PC contributions.", call. = FALSE)
pc_long <- pc_long[, .(
  sample_key, sample_id, treatment_group, projection_class, age_week,
  model_family, selector, pc, pc_score, calibrated_pc_contribution_week,
  predicted_age_week, delta_age_week, relative_delta_age_week
)]
pc_long[, pc := factor(pc, levels = selected_pcs)]
setorder(pc_long, sample_key, pc)
fwrite(pc_long, file.path(out_dir, "trig_cd1_pc_contributions_long.tsv"), sep = "\t")

pc_group_means <- pc_long[, .(
  n = .N,
  mean_pc_score = mean(pc_score, na.rm = TRUE),
  sd_pc_score = sd(pc_score, na.rm = TRUE),
  mean_pc_contribution_week = mean(calibrated_pc_contribution_week, na.rm = TRUE),
  sd_pc_contribution_week = sd(calibrated_pc_contribution_week, na.rm = TRUE),
  se_pc_contribution_week = sd(calibrated_pc_contribution_week, na.rm = TRUE) / sqrt(.N),
  mean_predicted_age_week = mean(predicted_age_week, na.rm = TRUE),
  mean_delta_age_week = mean(delta_age_week, na.rm = TRUE)
), by = .(pc, treatment_group)]
pc_group_means[, pc := as.character(pc)]
setorder(pc_group_means, pc, treatment_group)
fwrite(pc_group_means, file.path(out_dir, "trig_cd1_pc_group_means.tsv"), sep = "\t")

exact_perm_p <- function(values, groups, group_a, group_b) {
  use <- groups %in% c(group_a, group_b) & is.finite(values)
  values <- values[use]
  groups <- groups[use]
  if (length(values) != 8 || sum(groups == group_a) != 4 || sum(groups == group_b) != 4) return(NA_real_)
  obs <- mean(values[groups == group_a]) - mean(values[groups == group_b])
  idx <- seq_along(values)
  cmb <- combn(idx, 4)
  perm <- apply(cmb, 2, function(a_idx) {
    mean(values[a_idx]) - mean(values[setdiff(idx, a_idx)])
  })
  mean(abs(perm) >= abs(obs) - 1e-12)
}

smd_two_groups <- function(values, groups, group_a, group_b) {
  x <- values[groups == group_a]
  y <- values[groups == group_b]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  sx <- stats::var(x, na.rm = TRUE)
  sy <- stats::var(y, na.rm = TRUE)
  pooled <- sqrt(((length(x) - 1) * sx + (length(y) - 1) * sy) / (length(x) + length(y) - 2))
  if (!is.finite(pooled) || pooled == 0) return(NA_real_)
  (mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / pooled
}

pc_summary_path <- file.path(project_dir, "results", "integrated_diagnostics", "10_pca_diagnostics", "final_pca_clock_selected_pc_summary.tsv")
pc_summary <- if (file.exists(pc_summary_path)) fread(pc_summary_path) else data.table()
pc_meta <- data.table(pc = selected_pcs)
if (nrow(pc_summary) > 0) {
  deployment_rows <- pc_summary[artifact_scope == "deployment_all_control_like_external_projection" & pc %in% selected_pcs]
  if (nrow(deployment_rows) == length(selected_pcs)) {
    pc_meta <- deployment_rows[, .(
      pc,
      pc_variance_explained = variance_explained,
      pc_cumulative_variance = cumulative_variance_through_pc,
      pc_pearson_age_cor = pc_age_pearson,
      pc_spearman_age_cor = pc_age_spearman,
      clock_direction = clock_coefficient_direction,
      calibrated_pc_coefficient
    )]
  }
}
if (!"pc_variance_explained" %in% names(pc_meta)) {
  variance <- as.data.table(artifact$model$pc_variance_explained)
  if (!all(c("pc", "variance_explained") %in% names(variance))) {
    variance <- data.table(pc = seq_len(nrow(variance)), variance_explained = variance[[1]])
  }
  variance[, pc := paste0("PC", pc)]
  variance <- variance[pc %in% selected_pcs]
  variance[, pc_cumulative_variance := cumsum(variance_explained)]
  coef_dt <- as.data.table(artifact$model$pc_coefficient_table)
  pc_meta <- merge(
    variance[, .(pc, pc_variance_explained = variance_explained, pc_cumulative_variance)],
    coef_dt[, .(
      pc,
      calibrated_pc_coefficient,
      clock_direction = ifelse(calibrated_pc_coefficient >= 0, "higher_PC_score_increases_predicted_age", "higher_PC_score_decreases_predicted_age")
    )],
    by = "pc",
    all.x = TRUE
  )
  pc_meta[, `:=`(pc_pearson_age_cor = NA_real_, pc_spearman_age_cor = NA_real_)]
}

mean_wide <- dcast(
  pc_group_means,
  pc ~ treatment_group,
  value.var = "mean_pc_contribution_week"
)
setnames(mean_wide, c("young_3mo_control", "ama_9mo_control", "trg_9mo"), c(
  "mean_contribution_young_3mo",
  "mean_contribution_ama_9mo",
  "mean_contribution_trg_9mo"
))
pc_diff <- merge(pc_meta, mean_wide, by = "pc", all.x = TRUE, sort = FALSE)
pc_diff[, `:=`(
  age_shift_pc_ama_minus_young = mean_contribution_ama_9mo - mean_contribution_young_3mo,
  trg_effect_pc_trg_minus_ama = mean_contribution_trg_9mo - mean_contribution_ama_9mo
)]
pc_diff[, reversal_fraction_pc := fifelse(
  abs(age_shift_pc_ama_minus_young) < 0.25,
  NA_real_,
  -trg_effect_pc_trg_minus_ama / age_shift_pc_ama_minus_young
)]
pc_diff[, absolute_trg_effect_share_pc := abs(trg_effect_pc_trg_minus_ama) / sum(abs(trg_effect_pc_trg_minus_ama), na.rm = TRUE)]
pc_diff[, signed_fraction_of_total_pca_relative_delta := trg_effect_pc_trg_minus_ama / sum(trg_effect_pc_trg_minus_ama, na.rm = TRUE)]
pc_diff[, permutation_p_trg_vs_ama := vapply(as.character(pc), function(pcx) {
  dt <- pc_long[pc == pcx]
  exact_perm_p(dt$calibrated_pc_contribution_week, dt$treatment_group, "trg_9mo", "ama_9mo_control")
}, numeric(1))]
pc_diff[, permutation_p_ama_vs_young := vapply(as.character(pc), function(pcx) {
  dt <- pc_long[pc == pcx]
  exact_perm_p(dt$calibrated_pc_contribution_week, dt$treatment_group, "ama_9mo_control", "young_3mo_control")
}, numeric(1))]
pc_diff[, standardized_mean_difference_trg_vs_ama := vapply(as.character(pc), function(pcx) {
  dt <- pc_long[pc == pcx]
  smd_two_groups(dt$calibrated_pc_contribution_week, dt$treatment_group, "trg_9mo", "ama_9mo_control")
}, numeric(1))]
pc_diff[, pc_interpretation_direction := fifelse(
  trg_effect_pc_trg_minus_ama < -0.25,
  "TRG_younger_shift",
  fifelse(trg_effect_pc_trg_minus_ama > 0.25, "TRG_older_shift", "minimal_shift")
)]
pc_diff[, pc := factor(pc, levels = selected_pcs)]
setorder(pc_diff, pc)
pc_diff[, pc := as.character(pc)]
diff_cols <- c(
  "pc", "pc_variance_explained", "pc_cumulative_variance", "pc_pearson_age_cor",
  "pc_spearman_age_cor", "clock_direction", "calibrated_pc_coefficient",
  "mean_contribution_young_3mo", "mean_contribution_ama_9mo", "mean_contribution_trg_9mo",
  "age_shift_pc_ama_minus_young", "trg_effect_pc_trg_minus_ama", "reversal_fraction_pc",
  "absolute_trg_effect_share_pc", "signed_fraction_of_total_pca_relative_delta",
  "permutation_p_trg_vs_ama", "permutation_p_ama_vs_young",
  "standardized_mean_difference_trg_vs_ama", "pc_interpretation_direction"
)
fwrite(pc_diff[, ..diff_cols], file.path(out_dir, "trig_cd1_pc_trg_vs_ama_differences.tsv"), sep = "\t")

gsea_dirs <- c(
  file.path(project_dir, "results", "integrated_diagnostics", "11_final_pca_pc_gsea"),
  file.path(project_dir, "results", "integrated_diagnostics", "12_PCA deployment model GSEA")
)
gsea_dir <- gsea_dirs[dir.exists(gsea_dirs)][1] %||% NA_character_
gsea_all_path <- if (!is.na(gsea_dir)) file.path(gsea_dir, "final_pca_pc_gsea_all_results.tsv") else NA_character_
gsea_joined <- !is.na(gsea_all_path) && file.exists(gsea_all_path)
gsea_all <- if (gsea_joined) fread(gsea_all_path) else data.table()
if (!gsea_joined) add_warning("GSEA pathway labels were unavailable; falling back to top gene loading annotations.")

top_gene_path <- file.path(project_dir, "results", "integrated_diagnostics", "10_pca_diagnostics", "final_pca_clock_top50_gene_pc_loadings.tsv")
top_genes <- if (file.exists(top_gene_path)) fread(top_gene_path) else data.table()

if (gsea_joined) {
  gsea_all <- gsea_all[pc %in% selected_pcs]
  gsea_all[, abs_NES := abs(NES)]
  setorder(gsea_all, pc, direction, fdr_within_pc, empirical_p, -abs_NES)
  pathway_rows <- rbindlist(lapply(selected_pcs, function(pcx) {
    direction <- pc_diff[pc == pcx, pc_interpretation_direction]
    side <- if (identical(direction, "TRG_younger_shift")) {
      "younger_aligned"
    } else if (identical(direction, "TRG_older_shift")) {
      "older_aligned"
    } else {
      "minimal_shift"
    }
    if (identical(side, "minimal_shift")) {
      candidates <- gsea_all[pc == pcx][order(fdr_within_pc, empirical_p, -abs_NES)][seq_len(min(.N, 3))]
    } else {
      candidates <- gsea_all[pc == pcx & direction == side][order(fdr_within_pc, empirical_p, -abs_NES)][seq_len(min(.N, 5))]
    }
    if (nrow(candidates) == 0) return(data.table())
    candidates[, `:=`(
      pc_interpretation_direction = pc_diff[pc == pcx, pc_interpretation_direction],
      pathway_side_used_for_trg_interpretation = side,
      annotation_rank = seq_len(.N)
    )]
    candidates
  }), fill = TRUE)
  pathway_cols <- c(
    "pc", "pc_interpretation_direction", "pathway_side_used_for_trg_interpretation",
    "collection", "gene_set_name", "NES", "empirical_p", "fdr_within_pc",
    "direction", "leading_edge_genes", "annotation_rank"
  )
  pathway_annotation <- pathway_rows[, ..pathway_cols]
} else if (nrow(top_genes) > 0) {
  gene_rows <- rbindlist(lapply(selected_pcs, function(pcx) {
    genes <- top_genes[pc == pcx][seq_len(min(.N, 20))]
    if (nrow(genes) == 0) return(data.table())
    data.table(
      pc = pcx,
      pc_interpretation_direction = pc_diff[pc == pcx, pc_interpretation_direction],
      pathway_side_used_for_trg_interpretation = "top_gene_loading_fallback",
      collection = "top_gene_loadings",
      gene_set_name = paste(genes$gene_id, collapse = ";"),
      NES = NA_real_,
      empirical_p = NA_real_,
      fdr_within_pc = NA_real_,
      direction = "gene_loading_fallback",
      leading_edge_genes = paste(genes$gene_id, collapse = ";"),
      annotation_rank = 1L
    )
  }), fill = TRUE)
  pathway_annotation <- gene_rows
} else {
  pathway_annotation <- data.table(
    pc = selected_pcs,
    pc_interpretation_direction = pc_diff$pc_interpretation_direction,
    pathway_side_used_for_trg_interpretation = "unavailable",
    collection = NA_character_,
    gene_set_name = NA_character_,
    NES = NA_real_,
    empirical_p = NA_real_,
    fdr_within_pc = NA_real_,
    direction = NA_character_,
    leading_edge_genes = NA_character_,
    annotation_rank = NA_integer_
  )
}
fwrite(pathway_annotation, file.path(out_dir, "trig_cd1_pc_pathway_annotation.tsv"), sep = "\t")

top_terms_by_pc <- pathway_annotation[!is.na(gene_set_name), .(
  top_pathways_for_shift = paste(head(gene_set_name, 3), collapse = "; "),
  top_leading_edge_genes_for_shift = paste(unique(unlist(strsplit(paste(head(leading_edge_genes, 3), collapse = ";"), ";", fixed = TRUE)))[1:min(12, length(unique(unlist(strsplit(paste(head(leading_edge_genes, 3), collapse = ";"), ";", fixed = TRUE)))))], collapse = ";")
), by = pc]

pc_interpretation_summary <- merge(pc_diff, top_terms_by_pc, by = "pc", all.x = TRUE, sort = FALSE)
pc_interpretation_summary[, effect_rank_by_abs_trg_shift := frank(-abs(trg_effect_pc_trg_minus_ama), ties.method = "first")]
pc_interpretation_summary[, contribution_share := absolute_trg_effect_share_pc]
pc_interpretation_summary[, caution_flag := fifelse(
  abs(trg_effect_pc_trg_minus_ama) < 0.25,
  "minimal_contribution",
  fifelse(!is.na(pc_pearson_age_cor) & abs(pc_pearson_age_cor) < 0.15, "weak_age_correlation", "candidate_axis")
)]
pc_interpretation_summary[, main_interpretation_sentence := paste0(
  pc, " shifts TRG ",
  fifelse(trg_effect_pc_trg_minus_ama < -0.25, "younger", fifelse(trg_effect_pc_trg_minus_ama > 0.25, "older", "minimally")),
  " relative to 9-month AMA controls by ",
  sprintf("%.3f", trg_effect_pc_trg_minus_ama), " weeks."
)]
pc_interpretation_summary[, age_corr := pc_pearson_age_cor]
summary_cols <- c(
  "pc", "pc_variance_explained", "clock_direction", "age_corr",
  "trg_effect_pc_trg_minus_ama", "effect_rank_by_abs_trg_shift",
  "contribution_share", "main_interpretation_sentence",
  "top_pathways_for_shift", "top_leading_edge_genes_for_shift", "caution_flag"
)
summary_table <- pc_interpretation_summary[, ..summary_cols]
setnames(summary_table, c("pc_variance_explained", "trg_effect_pc_trg_minus_ama"), c("variance_explained", "trg_effect_pc"))
fwrite(summary_table, file.path(out_dir, "trig_cd1_pc_interpretation_summary_table.tsv"), sep = "\t")

reversal <- pc_diff[, .(
  pc,
  age_shift_pc_ama_minus_young,
  trg_effect_pc_trg_minus_ama,
  reversal_fraction_pc
)]
reversal[, reversal_class := fifelse(
  abs(age_shift_pc_ama_minus_young) < 0.25 | abs(trg_effect_pc_trg_minus_ama) < 0.25 | is.na(reversal_fraction_pc),
  "unstable_or_minimal",
  fifelse(reversal_fraction_pc > 0, "opposes_external_aging_shift", "same_direction_as_external_aging_shift")
)]
reversal[, interpretation := fifelse(
  reversal_class == "opposes_external_aging_shift",
  paste0(pc, " TRG effect opposes the external 3mo-to-9mo AMA ageing shift."),
  fifelse(
    reversal_class == "same_direction_as_external_aging_shift",
    paste0(pc, " TRG effect moves in the same direction as the external 3mo-to-9mo AMA ageing shift."),
    paste0(pc, " has a small TRG effect or unstable AMA-vs-young denominator.")
  )
)]
fwrite(reversal, file.path(out_dir, "trig_cd1_pc_reversal_against_aging_axis.tsv"), sep = "\t")

total_pca_relative_delta <- group_summary[
  model_family == "pca_glmnet" &
    selector == "version_b_age_dataset_penalized" &
    treatment_group == "trg_9mo",
  mean_relative_delta_age_week
][1]
sum_trg_effect_pc <- sum(pc_diff$trg_effect_pc_trg_minus_ama, na.rm = TRUE)
qc <- data.table(
  total_pca_relative_delta_from_prediction_table = total_pca_relative_delta,
  sum_trg_effect_pc = sum_trg_effect_pc,
  difference = sum_trg_effect_pc - total_pca_relative_delta,
  tolerance_pass = abs(sum_trg_effect_pc - total_pca_relative_delta) <= 0.1
)
if (!qc$tolerance_pass) {
  add_warning(
    paste0(
      "Sum of PC TRG-vs-AMA effects differs from PCA relative delta by ",
      round(qc$difference, 4), " weeks."
    )
  )
}
fwrite(qc, file.path(out_dir, "trig_cd1_pca_relative_delta_decomposition_qc.tsv"), sep = "\t")

theme_set(theme_bw(base_size = 14))
pc_long[, pc := factor(pc, levels = selected_pcs)]
pc_diff[, pc := factor(pc, levels = selected_pcs)]
pc_group_means[, pc := factor(pc, levels = selected_pcs)]

save_plot <- function(p, filename, width = 10, height = 7) {
  ggsave(file.path(plot_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(plot_dir, paste0(filename, ".pdf")), p, width = width, height = height, bg = "white", limitsize = FALSE)
}

p <- ggplot(pc_group_means, aes(pc, mean_pc_contribution_week, color = treatment_group, group = treatment_group)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_pointrange(aes(ymin = mean_pc_contribution_week - se_pc_contribution_week, ymax = mean_pc_contribution_week + se_pc_contribution_week),
    position = position_dodge(width = 0.45), linewidth = 0.5
  ) +
  labs(
    title = "Trig CD1 external projection: PC-level age contributions",
    x = NULL, y = "Mean calibrated PC contribution (weeks)", color = NULL
  )
save_plot(p, "01_trig_cd1_pc_group_contributions")

p <- ggplot(pc_diff, aes(pc, trg_effect_pc_trg_minus_ama, fill = trg_effect_pc_trg_minus_ama < 0)) +
  geom_hline(yintercept = 0, color = "grey45") +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", trg_effect_pc_trg_minus_ama)), vjust = ifelse(pc_diff$trg_effect_pc_trg_minus_ama < 0, 1.3, -0.4), size = 3.5) +
  scale_fill_manual(values = c("TRUE" = "#2b6cb0", "FALSE" = "#b8322a"), labels = c("TRUE" = "Younger vs AMA", "FALSE" = "Older vs AMA")) +
  labs(x = NULL, y = "TRG - AMA contribution (weeks)", fill = NULL)
save_plot(p, "02_trig_cd1_trg_vs_ama_pc_shift")

p <- ggplot(pc_diff, aes(pc, absolute_trg_effect_share_pc)) +
  geom_col(fill = "#3f7f93", width = 0.65) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "Share of absolute TRG-vs-AMA PC shift")
save_plot(p, "03_trig_cd1_pc_shift_share")

age_trg_long <- melt(
  pc_diff[, .(pc, age_shift_pc_ama_minus_young, trg_effect_pc_trg_minus_ama)],
  id.vars = "pc",
  variable.name = "comparison",
  value.name = "shift_week"
)
age_trg_long[, comparison := factor(comparison,
  levels = c("age_shift_pc_ama_minus_young", "trg_effect_pc_trg_minus_ama"),
  labels = c("AMA 9mo - Young 3mo", "TRG 9mo - AMA 9mo")
)]
p <- ggplot(age_trg_long, aes(pc, shift_week, fill = comparison)) +
  geom_hline(yintercept = 0, color = "grey45") +
  geom_col(position = position_dodge(width = 0.72), width = 0.65) +
  labs(x = NULL, y = "PC contribution shift (weeks)", fill = NULL)
save_plot(p, "04_trig_cd1_ama_vs_young_age_axis_and_trg_reversal")

p <- ggplot(reversal, aes(factor(pc, levels = selected_pcs), reversal_fraction_pc, color = reversal_class)) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", color = "grey50") +
  geom_point(size = 3, na.rm = TRUE) +
  labs(x = NULL, y = "Reversal fraction (-TRG effect / AMA-young shift)", color = NULL)
save_plot(p, "05_trig_cd1_pc_reversal_fraction")

heat_dt <- merge(
  pc_long,
  pca_pred[, .(sample_key, sample_column, predicted_age_week)],
  by = "sample_key",
  all.x = TRUE,
  suffixes = c("", "_pred")
)
sample_order <- heat_dt[, .(treatment_group = treatment_group[1], predicted_age_week = predicted_age_week[1]), by = sample_key][
  order(factor(treatment_group, levels = c("young_3mo_control", "ama_9mo_control", "trg_9mo")), predicted_age_week),
  sample_key
]
heat_dt[, sample_key := factor(sample_key, levels = sample_order)]
p <- ggplot(heat_dt, aes(pc, sample_key, fill = calibrated_pc_contribution_week)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b8322a", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Contribution\n(weeks)") +
  theme(panel.grid = element_blank(), axis.text.y = element_text(size = 8))
save_plot(p, "06_trig_cd1_sample_pc_contribution_heatmap", width = 8.5, height = 7.5)

group_heat <- copy(pc_group_means)
group_heat[, treatment_group := factor(treatment_group, levels = c("young_3mo_control", "ama_9mo_control", "trg_9mo"))]
p <- ggplot(group_heat, aes(pc, treatment_group, fill = mean_pc_contribution_week)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(low = "#2b6cb0", mid = "white", high = "#b8322a", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "Mean contribution\n(weeks)") +
  theme(panel.grid = element_blank())
save_plot(p, "07_trig_cd1_group_mean_pc_contribution_heatmap", width = 8.5, height = 4.8)

plot_pathways <- pathway_annotation[!is.na(gene_set_name)]
if (nrow(plot_pathways) > 0 && any(is.finite(plot_pathways$NES))) {
  plot_pathways[, neglog := -log10(pmax(fdr_within_pc %||% empirical_p, 1e-300))]
  plot_pathways[, short_name := substr(gene_set_name, 1, 70)]
  p <- ggplot(plot_pathways, aes(factor(pc, levels = selected_pcs), reorder(short_name, NES), color = NES, size = neglog)) +
    geom_point(alpha = 0.85) +
    scale_color_gradient2(low = "#2b6cb0", mid = "grey85", high = "#b8322a", midpoint = 0) +
    labs(x = NULL, y = NULL, color = "NES", size = "-log10(FDR)") +
    theme(axis.text.y = element_text(size = 7))
  save_plot(p, "08_trig_cd1_pc_pathway_annotation_dotplot", width = 10, height = 8)
}

score_long <- melt(
  scores_raw[, c("sample_key", "sample_column", "treatment_group", selected_pcs), with = FALSE],
  id.vars = c("sample_key", "sample_column", "treatment_group"),
  variable.name = "pc",
  value.name = "pc_score"
)
score_long[, pc := factor(pc, levels = selected_pcs)]
p <- ggplot(score_long, aes(pc, pc_score, color = treatment_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.75) +
  labs(x = NULL, y = "Frozen PCA score", color = NULL)
save_plot(p, "09_trig_cd1_pc_score_shift")

waterfall <- pc_diff[, .(label = as.character(pc), value = trg_effect_pc_trg_minus_ama)]
waterfall[, start := shift(cumsum(value), fill = 0)]
waterfall[, end := start + value]
total_row <- data.table(label = "Total PCA\nrelative delta", value = total_pca_relative_delta, start = 0, end = total_pca_relative_delta)
waterfall_plot <- rbind(waterfall, total_row, fill = TRUE)
waterfall_plot[, idx := seq_len(.N)]
p <- ggplot(waterfall_plot, aes(x = idx, fill = value < 0)) +
  geom_rect(aes(xmin = idx - 0.35, xmax = idx + 0.35, ymin = pmin(start, end), ymax = pmax(start, end))) +
  geom_hline(yintercept = 0, color = "grey45") +
  scale_x_continuous(breaks = waterfall_plot$idx, labels = waterfall_plot$label) +
  scale_fill_manual(values = c("TRUE" = "#2b6cb0", "FALSE" = "#b8322a"), guide = "none") +
  labs(x = NULL, y = "Cumulative TRG-vs-AMA contribution (weeks)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p, "10_trig_cd1_waterfall_trg_relative_delta", width = 9, height = 6)

elasticnet_rel <- group_summary[model_family == "elasticnet" & treatment_group == "trg_9mo", mean_relative_delta_age_week][1]
concordance_cor <- if (nrow(concordance) > 1) {
  round(cor(concordance$delta_age_elasticnet, concordance$delta_age_pca_glmnet, use = "complete.obs"), 3)
} else {
  NA_real_
}

ranked <- pc_diff[order(-abs(trg_effect_pc_trg_minus_ama))]
younger_pcs <- pc_diff[trg_effect_pc_trg_minus_ama < -0.25, paste0(pc, " (", round(trg_effect_pc_trg_minus_ama, 3), " weeks)")]
older_pcs <- pc_diff[trg_effect_pc_trg_minus_ama > 0.25, paste0(pc, " (+", round(trg_effect_pc_trg_minus_ama, 3), " weeks)")]

pc_detail_lines <- unlist(lapply(selected_pcs, function(pcx) {
  row <- pc_interpretation_summary[pc == pcx]
  rev_row <- reversal[pc == pcx]
  terms <- pathway_annotation[pc == pcx & !is.na(gene_set_name), head(gene_set_name, 3)]
  genes <- row$top_leading_edge_genes_for_shift %||% NA_character_
  c(
    paste0("### ", pcx),
    "",
    paste0(
      "- Variance explained: ", round(row$pc_variance_explained, 4),
      "; age correlation: ", ifelse(is.na(row$pc_pearson_age_cor), "NA", round(row$pc_pearson_age_cor, 3)),
      "; clock direction: ", row$clock_direction, "."
    ),
    paste0(
      "- TRG-vs-AMA contribution difference: ",
      round(row$trg_effect_pc_trg_minus_ama, 3), " weeks; direction: ",
      row$pc_interpretation_direction, "."
    ),
    paste0(
      "- External ageing-axis comparison: AMA-young shift = ",
      round(row$age_shift_pc_ama_minus_young, 3),
      " weeks; reversal class = ", rev_row$reversal_class, "."
    ),
    paste0("- Top relevant pathway terms: ", if (length(terms) > 0) paste(terms, collapse = "; ") else "unavailable", "."),
    paste0("- Leading-edge genes: ", if (!is.na(genes) && nzchar(genes)) genes else "unavailable", "."),
    paste0("- Claim grade: ", row$caution_flag, "."),
    ""
  )
}))

md <- c(
  "# Trig CD1 PC-Level Interpretation",
  "",
  "## Analysis Identity",
  "",
  "- Dataset: Trig CD1 external bulk RNA-seq projection.",
  "- Model: final integrated PCA/glmnet deployment clock with Version B selector.",
  "- PCs interpreted: PC1-PC7 only.",
  "- No model was retrained. No scaling, PCA, glmnet, calibration, gene ranking, or GSEA was refit.",
  "",
  "## Main PCA Result",
  "",
  paste0("- PCA/glmnet TRG mean relative delta vs 9-month AMA controls: ", round(total_pca_relative_delta, 3), " weeks."),
  "",
  "## PC-Level Decomposition",
  "",
  paste0("- Ranked PC effects by absolute TRG-vs-AMA contribution: ", paste(paste0(ranked$pc, "=", round(ranked$trg_effect_pc_trg_minus_ama, 3)), collapse = ", "), "."),
  paste0("- PCs pushing TRG younger: ", if (length(younger_pcs) > 0) paste(younger_pcs, collapse = ", ") else "none", "."),
  paste0("- PCs pushing TRG older: ", if (length(older_pcs) > 0) paste(older_pcs, collapse = ", ") else "none", "."),
  paste0("- Sum of PC TRG-vs-AMA effects: ", round(sum_trg_effect_pc, 3), " weeks; prediction-table relative delta: ", round(total_pca_relative_delta, 3), " weeks; difference: ", round(qc$difference, 4), " weeks."),
  "",
  "## PC-by-PC Interpretation",
  "",
  pc_detail_lines,
  "## External Ageing-Axis Comparison",
  "",
  paste0(
    "- PCs opposing the external AMA-vs-young ageing shift: ",
    paste(reversal[reversal_class == "opposes_external_aging_shift", pc], collapse = ", ") %||% "none",
    "."
  ),
  paste0(
    "- PCs moving in the same direction as the external AMA-vs-young ageing shift: ",
    paste(reversal[reversal_class == "same_direction_as_external_aging_shift", pc], collapse = ", ") %||% "none",
    "."
  ),
  paste0(
    "- PCs with unstable or minimal reversal estimates: ",
    paste(reversal[reversal_class == "unstable_or_minimal", pc], collapse = ", ") %||% "none",
    "."
  ),
  "",
  "## ElasticNet Context",
  "",
  paste0("- ElasticNet TRG mean relative delta vs AMA: ", round(elasticnet_rel, 3), " weeks."),
  paste0("- PCA/glmnet TRG mean relative delta vs AMA: ", round(total_pca_relative_delta, 3), " weeks."),
  paste0("- ElasticNet/PCA delta-age correlation: ", concordance_cor, "."),
  "- ElasticNet has no PC axes, so PC interpretation applies only to PCA/glmnet.",
  "",
  "## Interpretation Caution",
  "",
  "This is an external projection of CD1 bulk RNA-seq into a C57-compatible integrated bulk+pseudobulk ovarian ageing reference. PC contributions position samples along trained transcriptomic ageing axes, but they do not by themselves prove rejuvenation, accelerated ageing, or causality.",
  "",
  "## Claim Grading",
  "",
  "- PCs with strong TRG shift, good pathway annotation, and direction opposing AMA-vs-young shift can be described as candidate axes.",
  "- PCs with small shifts or weak age correlation should be described as predictive modifiers.",
  "- No PC should be described as causal without independent validation."
)
writeLines(md, file.path(out_dir, "trig_cd1_pc_interpretation_summary.md"))

if (length(warnings_vec) == 0) warnings_vec <- "No PC interpretation warnings."
writeLines(warnings_vec, file.path(out_dir, "trig_cd1_pc_interpretation_warnings.txt"))

tables_written <- c(
  "trig_cd1_pc_contributions_long.tsv",
  "trig_cd1_pc_group_means.tsv",
  "trig_cd1_pc_trg_vs_ama_differences.tsv",
  "trig_cd1_pc_pathway_annotation.tsv",
  "trig_cd1_pc_interpretation_summary_table.tsv",
  "trig_cd1_pc_reversal_against_aging_axis.tsv",
  "trig_cd1_pca_relative_delta_decomposition_qc.tsv"
)
plots_written <- list.files(plot_dir, pattern = "\\.(png|pdf)$")

message("Output directory: ", out_dir)
message("Total PCA TRG relative delta: ", round(total_pca_relative_delta, 3), " weeks")
younger_status <- pc_diff[trg_effect_pc_trg_minus_ama < 0][order(trg_effect_pc_trg_minus_ama)]
older_status <- pc_diff[trg_effect_pc_trg_minus_ama > 0][order(-trg_effect_pc_trg_minus_ama)]
message("Top 3 younger-shifting PCs: ", paste(head(younger_status$pc, 3), collapse = ", "))
message("Top 3 older-shifting PCs: ", paste(head(older_status$pc, 3), collapse = ", "))
message("GSEA annotation joined: ", gsea_joined)
message("Tables written: ", length(tables_written), "; plots written: ", length(plots_written))
