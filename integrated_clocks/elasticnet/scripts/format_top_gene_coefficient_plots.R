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
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
results_dir <- Sys.getenv("CLOCK_RESULTS_DIR", unset = file.path(project_dir, "results"))
gene_dir <- file.path(results_dir, "integrated_diagnostics", "09_selected_gene_diagnostics")
stopifnot(dir.exists(gene_dir))
elasticnet_dir <- file.path(gene_dir, "elasticnet_hvg_coefficients")
pca_effect_dir <- file.path(gene_dir, "pca_versionB_gene_effects")
pca_pc_dir <- file.path(gene_dir, "pca_versionB_pc_contributions")
invisible(lapply(c(elasticnet_dir, pca_effect_dir, pca_pc_dir), dir.create, recursive = TRUE, showWarnings = FALSE))
run_mode <- Sys.getenv("TOP_GENE_PLOT_MODE", unset = "all")
if (!run_mode %chin% c("all", "elasticnet", "pca")) {
  stop("TOP_GENE_PLOT_MODE must be one of: all, elasticnet, pca")
}

save_pdf_png <- function(plot, stem, width = 5.0, height = 7.6) {
  ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, bg = "white", limitsize = FALSE)
  ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
}

select_top_n_by_rank <- function(dt, rank_col, n = 20L) {
  out <- copy(dt)
  setorderv(out, rank_col)
  out[seq_len(min(.N, n))]
}

add_age_association_q <- function(dt, feature_file) {
  feature_dt <- fread(feature_file, sep = "\t")
  needed <- c("gene_id", "age_t", "residual_df")
  missing <- setdiff(needed, names(feature_dt))
  if (length(missing) > 0) {
    stop("Feature table lacks required columns: ", paste(missing, collapse = ", "))
  }
  feature_dt[, age_association_p := 2 * pt(-abs(age_t), df = residual_df)]
  feature_dt[, age_association_q := p.adjust(age_association_p, method = "BH")]
  merge(
    dt,
    feature_dt[, .(gene_id, age_t, residual_df, age_association_p, age_association_q)],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
  )
}

format_q_label <- function(q) {
  fifelse(
    is.na(q), "q = NA",
    fifelse(q < 1e-4, "q < 1e-4",
            fifelse(q < 0.01, sprintf("q = %.1e", q),
                    sprintf("q = %.2f", q)))
  )
}

plot_top_gene_coefficients <- function(dt, gene_col, value_col, x_label, title, subtitle) {
  plot_dt <- copy(dt)
  setnames(plot_dt, c(gene_col, value_col), c("gene_id", "signed_value"))
  plot_dt <- plot_dt[is.finite(signed_value)]
  plot_dt[, direction := fifelse(signed_value >= 0, "Positive", "Negative")]
  plot_dt[, direction_order := fifelse(direction == "Positive", 1L, 2L)]
  setorder(plot_dt, direction_order, -signed_value, gene_id)
  plot_dt[, gene_id := factor(gene_id, levels = rev(gene_id))]
  plot_dt[, q_label := format_q_label(age_association_q)]
  n_genes <- nrow(plot_dt)
  y_text_size <- if (n_genes > 25L) 6.1 else 7.4
  q_text_size <- if (n_genes > 25L) 1.75 else 2.05
  title_size <- if (n_genes > 25L) 9.3 else 10.2
  subtitle_size <- if (n_genes > 25L) 6.6 else 7.1

  bar_limit <- max(abs(plot_dt$signed_value), na.rm = TRUE)
  x_limit <- bar_limit * 1.48
  label_nudge <- bar_limit * 0.045
  if (!is.finite(x_limit) || x_limit <= 0) x_limit <- 1
  if (!is.finite(label_nudge) || label_nudge <= 0) label_nudge <- 0.04
  plot_dt[, q_x := fifelse(signed_value >= 0, signed_value + label_nudge, signed_value - label_nudge)]
  plot_dt[, q_hjust := fifelse(signed_value >= 0, 0, 1)]

  ggplot(plot_dt, aes(x = signed_value, y = gene_id, fill = direction)) +
    geom_col(width = 0.78, color = "black", linewidth = 0.28) +
    geom_text(
      aes(x = q_x, label = q_label, hjust = q_hjust),
      size = q_text_size,
      color = "black",
      family = "sans"
    ) +
    geom_vline(xintercept = 0, linewidth = 0.42, color = "black") +
    scale_fill_manual(
      values = c("Positive" = "#F04E37", "Negative" = "#4DB6CC"),
      guide = "none"
    ) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = pretty(c(-x_limit, x_limit), n = 5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      axis.line = element_line(linewidth = 0.42, color = "black"),
      axis.ticks = element_line(linewidth = 0.38, color = "black"),
      axis.ticks.length = unit(2.4, "pt"),
      axis.text.x = element_text(color = "black", size = 8.5),
      axis.text.y = element_text(color = "black", size = y_text_size, face = "italic"),
      axis.title.x = element_text(color = "black", size = 9.5, margin = margin(t = 8)),
      plot.title = element_text(color = "black", size = title_size, face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "black", size = subtitle_size, hjust = 0, margin = margin(t = 2, b = 5)),
      plot.margin = margin(t = 7, r = 8, b = 8, l = 8),
      panel.grid = element_blank()
    )
}

load_pca_gene_pc_contributions <- function() {
  pca_gene_file <- file.path(
    results_dir, "integrated_diagnostics", "10_pca_diagnostics",
    "final_pca_clock_gene_pc_loadings_with_feature_scores.tsv"
  )
  pca_gene <- fread(pca_gene_file, sep = "\t")
  pca_gene <- pca_gene[artifact_scope == "deployment_all_control_like_external_projection"]
  needed <- c("gene_id", "pc", "derived_gene_pc_coefficient", "abs_derived_gene_pc_coefficient")
  missing <- setdiff(needed, names(pca_gene))
  if (length(missing) > 0) {
    stop("PCA loading table lacks required columns: ", paste(missing, collapse = ", "))
  }
  for (optional_col in c("feature_score", "age_beta", "partial_r2_age", "partial_r2_dataset")) {
    if (!optional_col %in% names(pca_gene)) pca_gene[, (optional_col) := NA_real_]
  }
  pca_gene
}

read_pca_version_b_selected_top30 <- function(n = 30L) {
  selected_file <- file.path(
    results_dir,
    "selected_genes",
    "integrated_bulk_pseudobulk_pca_glmnet_versionB_top30_genes.tsv"
  )
  selected <- fread(selected_file, sep = "\t")
  required <- c(
    "gene_id", "age_beta", "age_t", "residual_df", "partial_r2_age",
    "partial_r2_dataset", "n_samples_used", "feature_score", "rank"
  )
  missing <- setdiff(required, names(selected))
  if (length(missing) > 0) {
    stop("PCA Version B selected top30 table lacks required columns: ", paste(missing, collapse = ", "))
  }
  selected <- selected[seq_len(min(.N, n))]
  selected[, `:=`(
    version_b_rank = seq_len(.N),
    selected_gene_rank = as.integer(rank),
    pc_gene_effect = as.numeric(age_beta),
    abs_pca_gene_effect = abs(as.numeric(age_beta)),
    pca_gene_effect_direction = fifelse(as.numeric(age_beta) >= 0, "older-predictor", "younger-predictor"),
    n_pcs_used = NA_integer_,
    abs_sum_pca_gene_effect = NA_real_,
    dominant_pc = NA_character_,
    dominant_pc_contribution = NA_real_,
    dominant_pc_abs_contribution = NA_real_,
    dominant_pc_share = NA_real_,
    selection_source = "selected_genes/integrated_bulk_pseudobulk_pca_glmnet_versionB_top30_genes.tsv",
    age_association_p = 2 * pt(-abs(as.numeric(age_t)), df = as.numeric(residual_df)),
    age_association_q = p.adjust(2 * pt(-abs(as.numeric(age_t)), df = as.numeric(residual_df)), method = "BH")
  )]
  setcolorder(selected, c(
    "version_b_rank", "selected_gene_rank", "gene_id", "feature_score",
    "age_beta", "age_t", "residual_df", "age_association_p",
    "age_association_q", "partial_r2_age", "partial_r2_dataset",
    "n_samples_used", "pc_gene_effect", "abs_pca_gene_effect",
    "pca_gene_effect_direction", "n_pcs_used", "abs_sum_pca_gene_effect",
    "dominant_pc", "dominant_pc_contribution",
    "dominant_pc_abs_contribution", "dominant_pc_share",
    "selection_source"
  ))
  selected[]
}

derive_pca_per_pc_top_genes <- function(pca_gene, n_per_direction = 10L) {
  pca_pc <- pca_gene[
    is.finite(derived_gene_pc_coefficient) & derived_gene_pc_coefficient != 0,
    .(
      pc,
      pc_number,
      gene_id,
      pc_loading,
      calibrated_pc_coefficient,
      pc_gene_contribution = derived_gene_pc_coefficient,
      abs_pc_gene_contribution = abs_derived_gene_pc_coefficient,
      pc_gene_effect_direction = fifelse(derived_gene_pc_coefficient >= 0, "older-predictor", "younger-predictor"),
      variance_explained = pc_age_pearson * NA_real_,
      feature_score,
      age_beta,
      age_t,
      partial_r2_age,
      partial_r2_dataset,
      n_samples_used,
      residual_df
    )
  ]
  variance_file <- file.path(
    results_dir, "integrated_diagnostics", "10_pca_diagnostics",
    "pca_variance_summary.tsv"
  )
  if (file.exists(variance_file)) {
    variance_dt <- fread(variance_file, sep = "\t")
    variance_dt <- variance_dt[
      selector == "version_b_age_dataset_penalized" & fit_stage == "train_fit",
      .(pc = paste0("PC", pc), variance_explained)
    ]
    variance_dt <- unique(variance_dt, by = "pc")
    pca_pc[, variance_explained := NULL]
    pca_pc <- merge(pca_pc, variance_dt, by = "pc", all.x = TRUE, sort = FALSE)
  }

  pca_pc <- pca_pc[, {
    older <- head(.SD[pc_gene_contribution > 0][order(-abs_pc_gene_contribution, gene_id)], n_per_direction)
    younger <- head(.SD[pc_gene_contribution < 0][order(-abs_pc_gene_contribution, gene_id)], n_per_direction)
    rbind(older, younger, fill = TRUE)
  }, by = pc]
  pca_pc[, rank_within_pc_direction := seq_len(.N), by = .(pc, pc_gene_effect_direction)]
  pca_pc[, age_association_p := 2 * pt(-abs(age_t), df = residual_df)]
  pca_pc[, age_association_q := p.adjust(age_association_p, method = "BH")]
  setorder(pca_pc, pc_number, pc_gene_effect_direction, rank_within_pc_direction)
  setcolorder(pca_pc, c(
    "pc", "pc_number", "variance_explained", "gene_id", "pc_gene_effect_direction",
    "rank_within_pc_direction", "pc_gene_contribution", "abs_pc_gene_contribution",
    "pc_loading", "calibrated_pc_coefficient", "feature_score", "age_beta",
    "age_t", "residual_df", "age_association_p", "age_association_q",
    "partial_r2_age", "partial_r2_dataset", "n_samples_used"
  ))
  pca_pc
}

plot_pca_per_pc_top_genes <- function(pca_per_pc, pc_name) {
  pc_dt <- copy(pca_per_pc[pc == pc_name])
  if (nrow(pc_dt) == 0) stop("No PCA per-PC rows found for ", pc_name)
  variance_label <- if (all(is.na(pc_dt$variance_explained))) {
    "variance explained unavailable"
  } else {
    paste0(round(unique(pc_dt$variance_explained)[1] * 100, 1), "% variance explained")
  }
  coefficient_label <- sprintf("calibrated PC coefficient %.3f", unique(pc_dt$calibrated_pc_coefficient)[1])
  subtitle <- paste0(
    "Top 10 per direction by |PC loading x calibrated PC coefficient|; ",
    variance_label, "; ", coefficient_label, "; q = BH-adjusted age association"
  )
  plot_top_gene_coefficients(
    pc_dt,
    gene_col = "gene_id",
    value_col = "pc_gene_contribution",
    x_label = paste0(pc_name, " gene contribution"),
    title = paste0("Top PCA Version B ", pc_name, " genes"),
    subtitle = paste(strwrap(subtitle, width = 76), collapse = "\n")
  )
}

if (run_mode %chin% c("all", "elasticnet")) {
  elastic_top <- fread(file.path(elasticnet_dir, "top30_final_elasticnet_hvg_coefficients.tsv"), sep = "\t")
  elastic_top <- add_age_association_q(
    elastic_top,
    file.path(results_dir, "elasticnet_clock", "integrated_bulk_pseudobulk",
              "version_b_age_dataset_penalized", "feature_score_table.tsv")
  )
  fwrite(
    elastic_top,
    file.path(elasticnet_dir, "top30_final_elasticnet_hvg_coefficients_plot_source.tsv"),
    sep = "\t"
  )
  p_elastic <- plot_top_gene_coefficients(
    elastic_top,
    gene_col = "gene_id",
    value_col = "coefficient",
    x_label = "Gene coefficient",
    title = "Top ElasticNet HVG genes",
    subtitle = "Deployment coefficients; q values are BH-adjusted age-association tests"
  )
  save_pdf_png(p_elastic, file.path(elasticnet_dir, "top30_final_elasticnet_hvg_coefficients"), height = 5.0)

  elastic_top20 <- select_top_n_by_rank(elastic_top, "rank", 20L)
  fwrite(
    elastic_top20[, names(fread(file.path(elasticnet_dir, "top30_final_elasticnet_hvg_coefficients.tsv"), nrows = 0)), with = FALSE],
    file.path(elasticnet_dir, "top20_final_elasticnet_hvg_coefficients.tsv"),
    sep = "\t"
  )
  fwrite(
    elastic_top20,
    file.path(elasticnet_dir, "top20_final_elasticnet_hvg_coefficients_plot_source.tsv"),
    sep = "\t"
  )
  p_elastic20 <- plot_top_gene_coefficients(
    elastic_top20,
    gene_col = "gene_id",
    value_col = "coefficient",
    x_label = "Gene coefficient",
    title = "Top 20 ElasticNet HVG genes",
    subtitle = "Top 20 by global |coefficient| rank; q values are BH-adjusted age-association tests"
  )
  save_pdf_png(p_elastic20, file.path(elasticnet_dir, "top20_final_elasticnet_hvg_coefficients"), width = 4.8, height = 5.0)
}

if (run_mode %chin% c("all", "pca")) {
  pca_gene <- load_pca_gene_pc_contributions()
  pca_top <- read_pca_version_b_selected_top30(30L)
  fwrite(pca_top, file.path(pca_effect_dir, "top30_final_pca_versionB_gene_effects.tsv"), sep = "\t")
  pca_per_pc <- derive_pca_per_pc_top_genes(pca_gene)
  fwrite(
    pca_per_pc,
    file.path(pca_pc_dir, "top_pca_versionB_genes_per_pc_final_contribution.tsv"),
    sep = "\t"
  )
  for (pc_name in c("PC1", "PC2")) {
    pc_source <- copy(pca_per_pc[pc == pc_name])
    fwrite(
      pc_source,
      file.path(pca_pc_dir, paste0("top_pca_versionB_", pc_name, "_genes_final_contribution_plot_source.tsv")),
      sep = "\t"
    )
    p_pc <- plot_pca_per_pc_top_genes(pca_per_pc, pc_name)
    save_pdf_png(
      p_pc,
      file.path(pca_pc_dir, paste0("top_pca_versionB_", pc_name, "_genes_final_contribution")),
      width = 5.8,
      height = 6.1
    )
  }
  fwrite(
    pca_top,
    file.path(pca_effect_dir, "top30_final_pca_versionB_gene_effects_plot_source.tsv"),
    sep = "\t"
  )
  pca_top20 <- pca_top[seq_len(min(.N, 20L))]
  fwrite(
    pca_top20[, names(fread(file.path(pca_effect_dir, "top30_final_pca_versionB_gene_effects.tsv"), nrows = 0)), with = FALSE],
    file.path(pca_effect_dir, "top20_final_pca_versionB_gene_effects.tsv"),
    sep = "\t"
  )
  fwrite(
    pca_top20,
    file.path(pca_effect_dir, "top20_final_pca_versionB_gene_effects_plot_source.tsv"),
    sep = "\t"
  )
  p_pca20 <- plot_top_gene_coefficients(
    pca_top20,
    gene_col = "gene_id",
    value_col = "age_beta",
    x_label = "Version B selected-gene age beta",
    title = "Top 20 PCA Version B genes",
    subtitle = "Top 20 in the Version B selected-gene table; q = BH-adjusted age association"
  )
  save_pdf_png(p_pca20, file.path(pca_effect_dir, "top20_final_pca_versionB_gene_effects"), width = 4.8, height = 5.0)
  p_pca <- plot_top_gene_coefficients(
    pca_top,
    gene_col = "gene_id",
    value_col = "age_beta",
    x_label = "Version B selected-gene age beta",
    title = "Top PCA Version B genes",
    subtitle = "Top 30 directly from the Version B selected-gene table; q = BH-adjusted age association"
  )
  save_pdf_png(p_pca, file.path(pca_effect_dir, "top30_final_pca_versionB_gene_effects"), height = 5.0)
}

cat("Formatted top-gene coefficient plots under:\n")
cat(gene_dir, "\n")
