#!/usr/bin/env Rscript
# Source table for all 3000 HVG-ElasticNet input genes (79 nonzero, the rest exactly 0): locked coefficient and
# each gene's signed age association (age t from the age + source-batch model on the model's own 110 training
# profiles). Used by make_coefficient_age_correlation_panels.R.
# Age statistics follow version_b_rank()/age_q() in make_diagnostics.R; nothing is refitted.
suppressPackageStartupMessages({
  library(data.table)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
en_dir <- file.path(fig_dir, 'diagnostics/09_selected_gene_diagnostics/elasticnet_hvg_coefficients')
source(file.path(lock, 'scripts/core.R'))  # fit_ranking(), choose_genes()

en_model <- readRDS(file.path(lock, 'HVG-ElasticNet-exclude-single-age/model.rds'))
en_coef <- fread(file.path(lock, 'HVG-ElasticNet-exclude-single-age/coefficients.csv'))
stopifnot(all.equal(en_coef$coefficient, en_model$coefficients$coefficient), length(en_model$training_keys) == 110)
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
en_q <- age_q(version_b_rank(en_model$training_keys))
hv <- merge(en_coef, en_q, by = 'gene', all.x = TRUE, sort = FALSE)
stopifnot(nrow(hv) == 3000, !anyNA(hv$age_t))
hv[, `:=`(nonzero = coefficient != 0,
          significant = fifelse(age_association_q < 0.05, 'q < 0.05', 'q ≥ 0.05'))]
# Must agree with the nonzero-gene table written by make_diagnostics.R.
nz <- fread(file.path(en_dir, 'nonzero_coefficient_vs_age_association_p_plot_source.tsv'))
chk <- merge(nz, hv, by.x = 'gene_id', by.y = 'gene')
stopifnot(nrow(chk) == 79, sum(hv$nonzero) == 79, isTRUE(all.equal(chk$age_t.x, chk$age_t.y)))
setnames(hv, 'gene', 'gene_id'); setorder(hv, -nonzero, -coefficient)
fwrite(hv, file.path(en_dir, 'all_hvg_coefficient_vs_age_t_plot_source.tsv'), sep = '\t')
cat(sprintf('%d HVG genes, %d nonzero, %d with q < 0.05\n', nrow(hv), sum(hv$nonzero), sum(hv$age_association_q < 0.05)))
