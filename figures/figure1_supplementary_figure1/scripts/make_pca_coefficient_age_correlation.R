#!/usr/bin/env Rscript
# Correlation between the locked Version B-PCA-glmnet gene-level weighted coefficients (PC loadings x PC
# coefficients, back-transformed to genes) and each gene's signed age association (age t from the
# age + source-batch model on the model's own 136 training profiles), for all selected genes.
# Age statistics follow version_b_rank()/age_q() in make_diagnostics.R; nothing is refitted.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
pca_dir <- file.path(fig_dir, 'diagnostics/09_selected_gene_diagnostics/pca_versionB_gene_effects')
source(file.path(lock, 'scripts/core.R'))  # fit_ranking(), choose_genes()

save_plot <- function(p, stem, width, height) {
  ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.pdf'), p, width = width, height = height, bg = 'white', limitsize = FALSE, device = cairo_pdf)
}

en_model <- readRDS(file.path(lock, 'HVG-ElasticNet-exclude-single-age/model.rds'))
pca_model <- readRDS(file.path(lock, 'VersionB-PCA-glmnet/model.rds'))
pca_coef <- fread(file.path(lock, 'VersionB-PCA-glmnet/coefficients.csv'))
stopifnot(all.equal(pca_coef$coefficient, pca_model$coefficients$coefficient))
dat <- readRDS(file.path(dev4, 'shared_inputs/prepared_expression.rds'))

# Same as make_diagnostics.R: Version B ranking uses the HVG-EN model config (models[[1]]$config).
version_b_rank <- function(keys) {
  m <- dat$metadata[match(keys, sample_key)]
  stopifnot(identical(m$sample_key, keys))
  fit_ranking(dat$expression[, keys, drop = FALSE], m, 'VersionB', en_model$config)
}
age_q <- function(rank) {
  m <- dat$metadata[match(rank$training_keys, sample_key)]
  q <- qr(model.matrix(~age + batch, data.frame(age = m$age_week, batch = factor(m$source_batch_for_feature_selection))))
  y <- t(dat$expression[rank$filtered_genes, rank$training_keys, drop = FALSE]); rdf <- nrow(y) - q$rank
  beta <- qr.coef(q, y)[2, ]; se <- sqrt(colSums(qr.resid(q, y)^2) / rdf * chol2inv(qr.R(q))[2, 2])
  tab <- data.table(gene = rank$filtered_genes, age_beta = beta, age_t = beta / se, residual_df = rdf)
  tab[, age_association_p := 2 * pt(-abs(age_t), df = rdf)][, age_association_q := p.adjust(age_association_p, method = 'BH')]
  stopifnot(all.equal(tab[match(rank$table$gene, gene), age_t], rank$table$age_t))
  tab
}
pca_rank <- version_b_rank(pca_model$training_keys)
stopifnot(identical(choose_genes(pca_rank, pca_model$hyperparameters$N), pca_model$representation$genes))
pca_q <- age_q(pca_rank)

pc <- merge(pca_coef, pca_q, by = 'gene', all.x = TRUE, sort = FALSE)
stopifnot(nrow(pc) == pca_model$hyperparameters$N, !anyNA(pc$age_t))
pc[, `:=`(version_b_rank = match(gene, pca_rank$table[order(-score, gene, na.last = TRUE)]$gene),
          significant = fifelse(age_association_q < 0.05, 'q < 0.05', 'q ≥ 0.05'))]
setnames(pc, 'gene', 'gene_id'); setorder(pc, -coefficient)
fwrite(pc, file.path(pca_dir, 'selected_gene_coefficient_vs_age_t_correlation_plot_source.tsv'), sep = '\t')

ct_p <- cor.test(pc$coefficient, pc$age_t, method = 'pearson')
ct_s <- suppressWarnings(cor.test(pc$coefficient, pc$age_t, method = 'spearman', exact = FALSE))
stats <- data.table(term = c('pearson', 'spearman'), estimate = c(ct_p$estimate, ct_s$estimate),
                    p_value = c(ct_p$p.value, ct_s$p.value), n = nrow(pc))
fwrite(stats, file.path(pca_dir, 'selected_gene_coefficient_vs_age_t_correlation_stats.tsv'), sep = '\t')
fmt_p <- function(p) if (p < 1e-15) 'p < 1e-15' else sprintf('p = %s', format(signif(p, 2)))
lab <- sprintf('Pearson r = %.2f, %s\nSpearman ρ = %.2f, %s',
               ct_p$estimate, fmt_p(ct_p$p.value), ct_s$estimate, fmt_p(ct_s$p.value))
# Label the 10 largest positive and 10 largest negative gene coefficients.
lab_dt <- rbind(pc[order(-coefficient)][1:10], pc[order(coefficient)][1:10])

p_cor <- ggplot(pc, aes(age_t, coefficient)) +
  geom_hline(yintercept = 0, color = 'grey70', linewidth = 0.3) +
  geom_vline(xintercept = 0, color = 'grey70', linewidth = 0.3) +
  geom_point(aes(color = significant), size = 0.9, alpha = 0.55, stroke = 0) +
  geom_smooth(method = 'lm', formula = y ~ x, color = 'grey20', fill = 'grey70', linewidth = 0.5) +
  ggrepel::geom_text_repel(data = lab_dt, aes(label = gene_id, color = significant), size = 2.2, max.overlaps = Inf,
                           segment.size = 0.2, min.segment.length = 0, seed = 1, show.legend = FALSE) +
  annotate('label', x = -Inf, y = Inf, label = lab, hjust = 0, vjust = 1, size = 2.6, fill = alpha('white', 0.8)) +
  scale_color_manual(values = c('q < 0.05' = '#B2182B', 'q ≥ 0.05' = 'grey55'), name = NULL) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(title = 'Version B–PCA gene coefficient vs gene–age association',
       subtitle = sprintf('All %d selected genes; gene coefficient = PC loadings × glmnet PC coefficients (%d PCs)\nx = age t statistic (age + source-batch model, %d training profiles); line = least-squares fit with 95%% CI\nlabelled: 10 largest positive and negative gene coefficients',
                          nrow(pc), pca_model$hyperparameters$pc, length(pca_model$training_keys)),
       x = expression('Age-association'~italic(t)~'statistic'), y = 'Gene coefficient (weeks per ln(1+CPM) unit)') +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = 'bold', size = 10), plot.subtitle = element_text(size = 7, color = 'grey30'),
        panel.grid.minor = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7.5))
save_plot(p_cor, file.path(pca_dir, 'selected_gene_coefficient_vs_age_t_correlation'), 6.5, 5.8)
print(stats)
cat('done\n')
