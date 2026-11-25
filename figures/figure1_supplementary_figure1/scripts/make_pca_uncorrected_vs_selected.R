#!/usr/bin/env Rscript
# Uncorrected PCA of all 136 profiles versus PCA restricted to the genes chosen by the rebutal4 selectors.
#   a  all shared genes: saved fit of the rebutal3 current-input batch assessment
#      (uncorrected_batch_diagnostics/03_PCA_centered_logCPM_PC1_PC2, age panel)
#   b  HVG-EN locked input genes (release hvg_exclusion_vbpca_20261012)
#   c  Version B-PCA locked input genes
#   d  HVG-EN genes with nonzero locked coefficients
# All panels: gene-centered lnCPM, no scaling, no batch correction, 136 profiles, sign convention of the
# assessment (largest |loading| positive). b and c are refitted here with the assessment's prcomp call
# (descriptive only; the locked Version B-PCA rotation uses standardized genes on its own training fit).
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
assess <- file.path(dirname(r4), 'rebutal3_batcheffects removal_CV_LOO/batch-effect assessment/current_input_assessment_2026-10-10')
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
out_dir <- file.path(fig_dir, 'diagnostics', '13_pca_scores')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- inputs
obj <- readRDS(file.path(assess, 'provenance/reconstructed_inputs/verified_expression.rds'))
meta <- as.data.table(obj$metadata)[, .(sample_key, study_id, modality, age_week)]
saved <- fread(file.path(assess, 'uncorrected_batch_diagnostics/centered_logCPM__all_filtered__scores.tsv'))
saved_var <- fread(file.path(assess, 'uncorrected_batch_diagnostics/pca_variance_explained.tsv'))[
  scaling == 'centered_logCPM' & gene_set == 'all_filtered']
stopifnot(nrow(meta) == 136, setequal(saved$sample_key, meta$sample_key))
dat <- readRDS(file.path(dev4, 'shared_inputs/prepared_expression.rds'))
hvg_model <- readRDS(file.path(lock, 'HVG-ElasticNet-exclude-single-age/model.rds'))
vb_model <- readRDS(file.path(lock, 'VersionB-PCA-glmnet/model.rds'))
keys <- meta$sample_key
stopifnot(setequal(keys, colnames(dat$expression)))

# Same PCA as 00_build_current_inputs.R of the assessment.
assessment_pca <- function(E) {
  E <- E[apply(E, 1, var) > 0, , drop = FALSE]
  fit <- prcomp(t(E), center = TRUE, scale. = FALSE)
  for (j in seq_len(ncol(fit$rotation))) {
    sgn <- sign(fit$rotation[which.max(abs(fit$rotation[, j])), j])
    fit$rotation[, j] <- fit$rotation[, j] * sgn; fit$x[, j] <- fit$x[, j] * sgn
  }
  list(x = fit$x, pve = fit$sdev^2 / sum(fit$sdev^2), n_genes = nrow(E))
}
# Sanity check: refitting panel a with this function reproduces the saved coordinates.
re_a <- assessment_pca(obj$expression[, saved$sample_key])
stopifnot(max(abs(re_a$x[, 1:2] - as.matrix(saved[, .(PC1, PC2)]))) < 1e-6)

# Locked gene sets; genes not measured in every one of the 136 profiles are dropped (HVG-EN only:
# it was trained without the single-age studies, some of which do not measure all of its genes).
gene_set_pca <- function(genes) {
  E <- dat$expression[genes, keys, drop = FALSE]
  keep <- rowSums(!is.finite(E)) == 0
  c(assessment_pca(E[keep, , drop = FALSE]), list(n_input = length(genes), n_dropped = sum(!keep)))
}
fit_b <- gene_set_pca(hvg_model$representation$genes)
fit_c <- gene_set_pca(vb_model$representation$genes)
stopifnot(fit_c$n_dropped == 0)
hvg_nonzero <- hvg_model$coefficients[coefficient != 0, gene]
stopifnot(all(hvg_nonzero %in% hvg_model$representation$genes))
fit_d <- gene_set_pca(hvg_nonzero)

panel_info <- data.table(
  panel = c('a', 'b', 'c', 'd'),
  title = c(sprintf('All %s shared genes', format(nrow(obj$expression), big.mark = ',')),
            sprintf('HVG–EN genes (n = %s)', format(fit_b$n_genes, big.mark = ',')),
            sprintf('Version B–PCA genes (n = %s)', format(fit_c$n_genes, big.mark = ',')),
            sprintf('HVG–EN nonzero-coefficient genes (n = %d)', fit_d$n_genes)),
  gene_set = c('all_shared_genes_rebutal3_assessment', 'HVG-ElasticNet-exclude-single-age', 'VersionB-PCA-glmnet',
               'HVG-ElasticNet-exclude-single-age_nonzero_coefficients'),
  n_input_genes = c(nrow(obj$expression), fit_b$n_input, fit_c$n_input, fit_d$n_input),
  n_genes_used = c(nrow(obj$expression), fit_b$n_genes, fit_c$n_genes, fit_d$n_genes),
  n_dropped_incomplete = c(0L, fit_b$n_dropped, fit_c$n_dropped, fit_d$n_dropped),
  PC1_pve = c(saved_var[component == 'PC1', fraction], fit_b$pve[1], fit_c$pve[1], fit_d$pve[1]),
  PC2_pve = c(saved_var[component == 'PC2', fraction], fit_b$pve[2], fit_c$pve[2], fit_d$pve[2]))
fwrite(panel_info, file.path(out_dir, '13_pca_uncorrected_vs_selected_panel_info.tsv'), sep = '\t')

scores <- rbind(
  saved[, .(sample_key, PC1, PC2, panel = 'a')],
  data.table(sample_key = rownames(fit_b$x), PC1 = fit_b$x[, 1], PC2 = fit_b$x[, 2], panel = 'b'),
  data.table(sample_key = rownames(fit_c$x), PC1 = fit_c$x[, 1], PC2 = fit_c$x[, 2], panel = 'c'),
  data.table(sample_key = rownames(fit_d$x), PC1 = fit_d$x[, 1], PC2 = fit_d$x[, 2], panel = 'd'))
scores <- merge(scores, meta, by = 'sample_key', sort = FALSE)
stopifnot(nrow(scores) == 4 * 136)
fwrite(scores, file.path(out_dir, '13_pca_uncorrected_vs_selected_plot_source.tsv'), sep = '\t')

# ---------------------------------------------------------------- figure
theme_set(theme_classic(base_size = 10) + theme(strip.background = element_blank(), plot.title = element_text(face = 'bold')))
age_lim <- range(meta$age_week)
make_panel <- function(p_id) {
  inf <- panel_info[panel == p_id]
  ggplot(scores[panel == p_id], aes(PC1, PC2, color = age_week)) +
    geom_point(size = 1.7) +
    scale_color_viridis_c(option = 'C', end = 0.9, limits = age_lim, name = 'Age (weeks)') +
    labs(title = inf$title, x = sprintf('PC1 (%.1f%%)', 100 * inf$PC1_pve), y = sprintf('PC2 (%.1f%%)', 100 * inf$PC2_pve))
}
p <- wrap_plots(lapply(c('a', 'b', 'c', 'd'), make_panel), ncol = 2) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a',
                  title = 'Uncorrected PCA: all genes versus genes chosen by the rebutal4 selectors',
                  subtitle = paste0('Gene-centered lnCPM, no batch correction, 136 profiles; PCA refitted for each gene set.',
                                    if (fit_b$n_dropped > 0) sprintf('\nHVG–EN: %d of %d genes not measured in all profiles were omitted.',
                                                                     fit_b$n_dropped, fit_b$n_input) else '',
                                    if (fit_d$n_dropped > 0) sprintf(' Nonzero-coefficient set: %d of %d omitted.',
                                                                     fit_d$n_dropped, fit_d$n_input) else '')) &
  theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold'))
stem <- file.path(out_dir, '13_pca_uncorrected_vs_selected_genes')
ggsave(paste0(stem, '.pdf'), p, width = 9.5, height = 9, device = cairo_pdf, bg = 'white')
ggsave(paste0(stem, '.png'), p, width = 9.5, height = 9, dpi = 300, bg = 'white')

# Three-panel row: all genes, Version B-PCA genes, HVG-EN nonzero-coefficient genes (relettered a-c).
p3 <- wrap_plots(lapply(c('a', 'c', 'd'), make_panel), nrow = 1) +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a',
                  title = 'Uncorrected PCA: all genes, Version B–PCA genes and HVG–EN nonzero-coefficient genes',
                  subtitle = paste0('Gene-centered lnCPM, no batch correction, 136 profiles; PCA refitted for each gene set.',
                                    if (fit_d$n_dropped > 0) sprintf(' HVG–EN: %d of %d nonzero-coefficient genes not measured in all profiles were omitted.',
                                                                     fit_d$n_dropped, fit_d$n_input) else '')) &
  theme(legend.position = 'bottom', plot.tag = element_text(face = 'bold'))
stem3 <- file.path(out_dir, '13_pca_all_vs_versionB_vs_hvgEN_nonzero')
ggsave(paste0(stem3, '.pdf'), p3, width = 13, height = 4.8, device = cairo_pdf, bg = 'white')
ggsave(paste0(stem3, '.png'), p3, width = 13, height = 4.8, dpi = 300, bg = 'white')
print(panel_info)
cat('done\n')
