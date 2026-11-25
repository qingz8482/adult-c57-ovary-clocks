#!/usr/bin/env Rscript
# Two diagnostics of bulk-pseudobulk_clock/results/integrated_diagnostics, reproduced for the locked
# release hvg_exclusion_vbpca_20261012:
#   12  age vs source-batch partial R2 (12_versionB_partial_r2_integrated_pca_glmnet), HVG vs Version B
#   13  PC1/PC2 sample scores (13_pca_diagnostics_pc1_pc2_scores)
# Partial R2 uses the Version B age + source-batch regression of fit_ranking(), computed for every gene
# of each locked input set on that model's own training profiles (Version B-PCA 136, HVG-EN 110).
# PCA scores are the locked Version B-PCA rotation applied to its own training profiles. HVG-EN has no
# PCA step; its panel is an unsupervised PCA of the locked HVG-EN standardized inputs (diagnostic only).
# Nothing is refitted.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
out_dir <- file.path(fig_dir, 'diagnostics')
r2_dir <- file.path(out_dir, '12_partial_r2_age_vs_source'); pca_dir <- file.path(out_dir, '13_pca_scores')
for (d in c(r2_dir, pca_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
source(file.path(lock, 'scripts/core.R'))  # fit_ranking(), choose_genes(), project_representation()

save_plot <- function(p, stem, width, height) {
  ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.pdf'), p, width = width, height = height, bg = 'white', limitsize = FALSE, device = cairo_pdf)
}
theme_pub <- function() {
  theme_classic(base_size = 11) +
    theme(axis.line = element_line(linewidth = 0.42, color = 'black'),
          axis.ticks = element_line(linewidth = 0.38, color = 'black'),
          axis.text = element_text(color = 'black'), axis.title = element_text(color = 'black', face = 'bold'),
          plot.title = element_text(face = 'bold', size = 11), plot.subtitle = element_text(color = 'grey30', size = 8),
          strip.background = element_blank(), strip.text = element_text(face = 'bold', size = 10))
}

# ---------------------------------------------------------------- inputs
hvg_model <- readRDS(file.path(lock, 'HVG-ElasticNet-exclude-single-age/model.rds'))
vb_model <- readRDS(file.path(lock, 'VersionB-PCA-glmnet/model.rds'))
dat <- readRDS(file.path(dev4, 'shared_inputs/prepared_expression.rds'))
meta <- as.data.table(dat$metadata)
vb_keys <- vb_model$training_keys; hvg_keys <- hvg_model$training_keys
stopifnot(length(vb_keys) == 136, length(hvg_keys) == 110, all(hvg_keys %in% vb_keys))

# ================================================================= 12 | partial R2
# The regenerated ranking must reproduce the locked Version B gene list exactly.
m_vb <- meta[match(vb_keys, sample_key)]
vb_rank <- fit_ranking(dat$expression[, vb_keys, drop = FALSE], m_vb, 'VersionB', vb_model$config)
stopifnot(identical(choose_genes(vb_rank, vb_model$hyperparameters$N), vb_model$representation$genes))

# Same formulas as fit_ranking(), for every gene of a set (fit_ranking drops NA-score genes). Each set is
# evaluated on its own model's training profiles (HVG-EN excludes the single-age studies, some of which
# do not measure all HVG genes).
partial_r2 <- function(genes, keys) {
  m <- meta[match(keys, sample_key)]
  y <- t(dat$expression[genes, keys, drop = FALSE]); stopifnot(!anyNA(y))
  d <- data.frame(age = m$age_week, batch = factor(m$source_batch_for_feature_selection))
  rss <- colSums(qr.resid(qr(model.matrix(~age + batch, d)), y)^2)
  ra <- colSums(qr.resid(qr(model.matrix(~batch, d)), y)^2)
  rd <- colSums(qr.resid(qr(model.matrix(~age, d)), y)^2)
  data.table(gene = genes, partial_r2_age = pmax(0, pmin(1, (ra - rss) / ra)),
             partial_r2_dataset = pmax(0, pmin(1, (rd - rss) / rd)), n_training_profiles = length(keys))
}
vb_r2 <- partial_r2(vb_model$representation$genes, vb_keys)
chk <- merge(vb_r2, vb_rank$table[, .(gene, pa = partial_r2_age, pd = partial_r2_dataset)], by = 'gene')
stopifnot(nrow(chk) == nrow(vb_r2), isTRUE(all.equal(chk$partial_r2_age, chk$pa)), isTRUE(all.equal(chk$partial_r2_dataset, chk$pd)))
vb_r2[, `:=`(set = 'Version B', version_b_score = vb_rank$table$score[match(gene, vb_rank$table$gene)],
             version_b_rank = match(gene, vb_rank$table$gene))]
hvg_r2 <- partial_r2(hvg_model$representation$genes, hvg_keys)
hvg_r2[, `:=`(set = 'HVG', hvg_en_coefficient = hvg_model$coefficients$coefficient[match(gene, hvg_model$coefficients$gene)])]
plot_dt <- rbind(hvg_r2, vb_r2, fill = TRUE)
plot_dt[, set := factor(set, levels = c('HVG', 'Version B'))]
plot_dt[, in_both_sets := gene %in% intersect(hvg_r2$gene, vb_r2$gene)]
setorder(plot_dt, set, gene)
fwrite(plot_dt, file.path(r2_dir, '12_partial_r2_age_vs_source_plot_source.tsv'), sep = '\t')

label_dt <- head(vb_r2[order(-version_b_score)], 12)
summ <- plot_dt[, .(n_genes = .N, median_partial_r2_age = median(partial_r2_age),
                    median_partial_r2_dataset = median(partial_r2_dataset),
                    frac_dataset_gt_0.9 = mean(partial_r2_dataset > 0.9), n_training_profiles = n_training_profiles[1]), by = set]
summ[, n_shared := sum(hvg_r2$gene %in% vb_r2$gene)]
fwrite(summ, file.path(r2_dir, '12_partial_r2_age_vs_source_summary.tsv'), sep = '\t')
print(summ)

set_col <- c('HVG' = 'grey60', 'Version B' = '#D2691E')
p_r2 <- ggplot(plot_dt, aes(partial_r2_age, partial_r2_dataset, color = set)) +
  geom_point(size = 1.1, alpha = 0.45, stroke = 0) +
  geom_hline(yintercept = 0.9, linetype = 'dashed', color = '#A0413C', linewidth = 0.45) +
  geom_text_repel(data = label_dt, aes(label = gene), color = 'black', fontface = 'italic', size = 3,
                  min.segment.length = 0.2, segment.size = 0.25, segment.color = 'grey50', max.overlaps = Inf, seed = 1) +
  scale_color_manual(values = set_col, name = NULL,
                     labels = c(HVG = sprintf('HVG (HVG–EN input, n = %d)', nrow(hvg_r2)),
                                `Version B` = sprintf('Version B (Version B–PCA input, n = %d)', nrow(vb_r2)))) +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.04))) +
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(x = 'Partial R2 age', y = 'Partial R2 dataset/source') +
  theme_pub() + theme(legend.position = c(0.02, 1.0), legend.justification = c(0, 1),
                      legend.background = element_blank(), legend.key = element_blank(), legend.key.size = unit(10, 'pt'))
save_plot(p_r2, file.path(r2_dir, '12_partial_r2_age_vs_source'), 5.2, 4.2)

# ================================================================= 13 | PC1/PC2 scores
sample_info <- meta[, .(sample_key, study, modality, age_week)]
# Version B-PCA: locked rotation applied to its training profiles.
vb_proj <- project_representation(vb_model$representation, dat$expression[, vb_keys, drop = FALSE])
stopifnot(all(vb_proj$missing_fraction == 0))
vb_var <- apply(vb_proj$x, 2, var)
# Locked rotation is fit on standardized training genes, so total variance = number of genes.
vb_pve <- vb_var / length(vb_model$representation$genes)
# HVG-EN: unsupervised PCA of its standardized training inputs (not part of the model).
hvg_z <- project_representation(hvg_model$representation, dat$expression[, hvg_keys, drop = FALSE])$x
hvg_pca <- prcomp(hvg_z, center = FALSE, scale. = FALSE, rank. = 10)
hvg_pve <- hvg_pca$sdev^2 / sum(hvg_pca$sdev^2)

panel_lab <- c(
  hvg = sprintf('HVG–EN inputs (unsupervised PCA, %d profiles)', length(hvg_keys)),
  vb = sprintf('Version B–PCA locked model (%d profiles)', length(vb_keys)))
scores <- rbind(
  data.table(sample_key = rownames(hvg_z), PC1 = hvg_pca$x[, 1], PC2 = hvg_pca$x[, 2], model = 'hvg',
             pc1_pve = hvg_pve[1], pc2_pve = hvg_pve[2]),
  data.table(sample_key = vb_keys, PC1 = vb_proj$x[, 1], PC2 = vb_proj$x[, 2], model = 'vb',
             pc1_pve = vb_pve[1], pc2_pve = vb_pve[2]))
scores <- merge(scores, sample_info, by = 'sample_key', all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(scores$age_week))
scores[, panel := factor(sprintf('%s\nPC1 %.1f%%, PC2 %.1f%%', panel_lab[model], 100 * pc1_pve, 100 * pc2_pve),
                         levels = unique(sprintf('%s\nPC1 %.1f%%, PC2 %.1f%%', panel_lab[model], 100 * pc1_pve, 100 * pc2_pve)))]
fwrite(scores, file.path(pca_dir, '13_pca_pc1_pc2_scores_plot_source.tsv'), sep = '\t')
fwrite(data.table(model = rep(c('HVG-EN inputs (unsupervised)', 'VersionB-PCA locked'), each = 10), pc = rep(1:10, 2),
                  variance_explained = c(hvg_pve[1:10], vb_pve[1:10]),
                  pc_coefficient = c(rep(NA_real_, 10), vb_model$pc_coefficients)),
       file.path(pca_dir, '13_pca_variance_explained.tsv'), sep = '\t')

p_pca <- ggplot(scores, aes(PC1, PC2, color = age_week, shape = modality)) +
  geom_point(size = 2.3, alpha = 0.9) +
  facet_wrap(~panel, scales = 'free', ncol = 2) +
  scale_color_viridis_c(option = 'C', name = 'Age (weeks)') +
  scale_shape_manual(values = c(bulk_rnaseq = 16, scrna_pseudobulk = 17), name = 'Modality') +
  labs(title = 'PCA score diagnostics of the locked clocks', x = 'PC1', y = 'PC2') +
  theme_pub() + theme(legend.position = 'bottom', panel.border = element_rect(fill = NA, color = 'black', linewidth = 0.4))
save_plot(p_pca, file.path(pca_dir, '13_pca_pc1_pc2_scores'), 9, 5)
cat('done\n')
