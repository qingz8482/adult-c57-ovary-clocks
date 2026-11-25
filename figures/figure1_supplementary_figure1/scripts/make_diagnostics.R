#!/usr/bin/env Rscript
# Diagnostics for the locked release hvg_exclusion_vbpca_20261012, reproducing three plots of
# bulk-pseudobulk_clock/results/integrated_diagnostics for the rebutal4 clocks:
#   07  absolute held-out error by biological age bin (04_integrated_clock_diagnostics.R)
#   09  top-30 HVG-EN coefficients and top-30 Version B genes (format_top_gene_coefficient_plots.R)
# Errors are out-of-study predictions from stages A/B/C; in-sample final-fit error is not shown.
# Version B age-association statistics are recomputed with the pipeline's own fit_ranking() on
# each locked model's training set and checked against the locked gene list. Nothing is refitted.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
r3 <- file.path(dirname(r4), 'rebutal3_batcheffects removal_CV_LOO/model_development')
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
out_dir <- file.path(fig_dir, 'diagnostics')
gene_dir <- file.path(out_dir, '09_selected_gene_diagnostics')
en_dir <- file.path(gene_dir, 'elasticnet_hvg_coefficients'); pca_dir <- file.path(gene_dir, 'pca_versionB_gene_effects')
for (d in c(out_dir, en_dir, pca_dir)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
source(file.path(lock, 'scripts/core.R'))  # fit_ranking(), choose_genes()

# ---------------------------------------------------------------- constants
recipes <- data.table(
  recipe = c('HVG-ElasticNet', 'VersionB-PCA'),
  dev    = c(dev4, r3),
  dir    = c('HVG-ElasticNet', 'VersionB-PCA/glmnet'),
  locked = c('HVG-ElasticNet-exclude-single-age', 'VersionB-PCA-glmnet'),
  label  = c('HVG–EN\n(single-age excl.)', 'Version B–PCA'))
rec_lab <- setNames(recipes$label, recipes$recipe)
stage_lab <- c(A = 'Stage A · 4 folds (OOF)', B = 'Stage B · NIHMS', C = 'Stage C · Winkler')
bio_levels <- c('young_adult', 'mature_adult', 'advanced_maternal_age', 'post_reproductive_ageing')
bio_bin_from_age <- function(age) fcase(age >= 7 & age < 16, 'young_adult', age >= 16 & age < 40, 'mature_adult',
                                        age >= 40 & age < 65, 'advanced_maternal_age', age >= 65, 'post_reproductive_ageing')

save_plot <- function(p, stem, width, height) {
  ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.pdf'), p, width = width, height = height, bg = 'white', limitsize = FALSE, device = cairo_pdf)
}
theme_diag <- function() {
  theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = 'bold'), plot.subtitle = element_text(color = 'grey30'),
          axis.text.x = element_text(angle = 30, hjust = 1), strip.text = element_text(face = 'bold'),
          legend.position = 'bottom')
}

# ---------------------------------------------------------------- inputs
rpath <- function(r, ...) file.path(recipes[recipe == r, dev], recipes[recipe == r, dir], ...)
pred <- rbindlist(lapply(recipes$recipe, function(r) rbindlist(lapply(c('A', 'B', 'C'), function(s)
  fread(rpath(r, 'count_tuning/results', paste0('predictions_', s, '.csv')))))))
stopifnot(pred[analysis == 'A', .N, by = recipe]$N == 136, uniqueN(pred[analysis == 'A']$sample_key) == 136)
models <- setNames(lapply(recipes$locked, function(l) readRDS(file.path(lock, l, 'model.rds'))), recipes$recipe)
# Locked artifacts must be the exact final fits of the development runs.
for (r in recipes$recipe) stopifnot(identical(models[[r]], readRDS(rpath(r, 'count_tuning/results/final_model/model.rds'))))
dat <- readRDS(file.path(dev4, 'shared_inputs/prepared_expression.rds'))
stopifnot(identical(unname(tools::md5sum(file.path(dev4, 'shared_inputs/prepared_expression.rds'))),
                    unname(tools::md5sum(file.path(r3, 'shared_inputs/prepared_expression.rds')))))

# ================================================================= 07 | error by age bin
ev <- copy(pred)
ev[, `:=`(residual_week = pred - age_week, abs_error_week = abs(pred - age_week),
          biological_age_bin = factor(bio_bin_from_age(age_week), levels = bio_levels),
          clock_label = factor(rec_lab[recipe], levels = recipes$label),
          stage_label = factor(stage_lab[analysis], levels = stage_lab))]
stopifnot(!anyNA(ev$biological_age_bin))
perf <- ev[, .(n = .N, n_studies = uniqueN(study), mae_week = mean(abs_error_week), median_abs_error_week = median(abs_error_week),
               rmse_week = sqrt(mean(residual_week^2)), mean_residual_week = mean(residual_week)),
           by = .(recipe, analysis, stage_label, biological_age_bin)]
setorder(perf, recipe, analysis, biological_age_bin)
fwrite(perf, file.path(out_dir, 'performance_by_biological_age_bin.tsv'), sep = '\t')
fwrite(ev[, .(recipe, analysis, fold, sample_key, study, modality, age_week, biological_age_bin, pred, residual_week, abs_error_week)],
       file.path(out_dir, '07_heldout_error_by_biological_age_bin_plot_source.tsv'), sep = '\t')
p <- ggplot(ev, aes(x = biological_age_bin, y = abs_error_week, color = biological_age_bin)) +
  geom_boxplot(outlier.shape = NA, color = 'grey30') +
  geom_point(alpha = 0.8, size = 2.3, position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  scale_x_discrete(drop = FALSE) + scale_color_discrete(drop = FALSE) +
  facet_grid(clock_label ~ stage_label) +
  labs(title = 'Held-out absolute error by biological age bin',
       subtitle = 'Out-of-study predictions of the locked recipes (stages A, B, C); no calibration; in-sample final-fit error not shown',
       x = 'Biological age bin', y = 'Absolute error (weeks)', color = 'Age bin') +
  theme_diag() + theme(legend.position = 'right')
save_plot(p, file.path(out_dir, '07_heldout_error_by_biological_age_bin'), 15, 7)

# ================================================================= 09 | top-30 genes
# Version B age-plus-source-batch regression on each locked model's own training samples.
version_b_rank <- function(keys) {
  m <- dat$metadata[match(keys, sample_key)]
  stopifnot(identical(m$sample_key, keys))
  fit_ranking(dat$expression[, keys, drop = FALSE], m, 'VersionB', models[[1]]$config)
}
# Age t-test from the same age + batch design, for every filtered gene (Version B scores are NA
# for some genes); BH is applied across all filtered genes.
age_q <- function(rank) {
  m <- dat$metadata[match(rank$training_keys, sample_key)]
  q <- qr(model.matrix(~age + batch, data.frame(age = m$age_week, batch = factor(m$source_batch_for_feature_selection))))
  y <- t(dat$expression[rank$filtered_genes, rank$training_keys, drop = FALSE]); rdf <- nrow(y) - q$rank
  beta <- qr.coef(q, y)[2, ]; se <- sqrt(colSums(qr.resid(q, y)^2) / rdf * chol2inv(qr.R(q))[2, 2])
  tab <- data.table(gene = rank$filtered_genes, age_t = beta / se, residual_df = rdf)
  tab[, age_association_p := 2 * pt(-abs(age_t), df = rdf)][, age_association_q := p.adjust(age_association_p, method = 'BH')]
  stopifnot(all.equal(tab[match(rank$table$gene, gene), age_t], rank$table$age_t))
  merge(tab, rank$table[, !'age_t'], by = 'gene', all.x = TRUE, sort = FALSE)[order(-score, gene, na.last = TRUE)]
}
format_q_label <- function(q) {
  fifelse(is.na(q), 'q = NA', fifelse(q < 1e-4, 'q < 1e-4', fifelse(q < 0.01, sprintf('q = %.1e', q), sprintf('q = %.2f', q))))
}
plot_top_gene_coefficients <- function(dt, gene_col, value_col, x_label, title, subtitle, symmetric = TRUE,
                                       x_break = NULL, pad_rows = 0L, y_text_size = NULL, q_text_size = NULL,
                                       x_n_breaks = 5, by_magnitude = FALSE) {
  plot_dt <- copy(dt)
  setnames(plot_dt, c(gene_col, value_col), c('gene_id', 'signed_value'))
  plot_dt <- plot_dt[is.finite(signed_value)]
  plot_dt[, direction := fifelse(signed_value >= 0, 'Positive', 'Negative')]
  plot_dt[, direction_order := fifelse(direction == 'Positive', 1L, 2L)]
  if (by_magnitude) plot_dt[, signed_value_order := -abs(signed_value)] else plot_dt[, signed_value_order := -signed_value]
  setorder(plot_dt, direction_order, signed_value_order, gene_id)
  plot_dt[, gene_id := factor(gene_id, levels = rev(gene_id))]
  plot_dt[, q_label := format_q_label(age_association_q)]
  n_genes <- nrow(plot_dt)
  if (is.null(y_text_size)) y_text_size <- if (n_genes > 25L) 6.1 else 7.4
  if (is.null(q_text_size)) q_text_size <- if (n_genes > 25L) 1.75 else 2.05
  title_size <- if (n_genes > 25L) 9.3 else 10.2
  subtitle_size <- if (n_genes > 25L) 6.6 else 7.1
  # Optional positive-side axis break: values in x_break = c(lo, hi) are squeezed into gap_w plot units.
  gap_w <- 0.6
  to_plot <- if (is.null(x_break)) identity else function(x)
    fifelse(x <= x_break[1], x, fifelse(x >= x_break[2], x - diff(x_break) + gap_w,
                                        x_break[1] + (x - x_break[1]) * gap_w / diff(x_break)))
  plot_dt[, plot_value := to_plot(signed_value)]
  bar_limit <- max(abs(if (is.null(x_break)) plot_dt$signed_value else plot_dt[signed_value <= x_break[1], signed_value]), na.rm = TRUE)
  x_limit <- bar_limit * 1.48; label_nudge <- bar_limit * 0.045
  # Asymmetric: each side spans its own longest bar plus room for the q labels.
  x_lims <- if (symmetric) c(-x_limit, x_limit) else
    c(min(0, plot_dt$plot_value) - 0.3 * bar_limit * any(plot_dt$plot_value < 0),
      max(0, plot_dt$plot_value) + 0.3 * bar_limit * any(plot_dt$plot_value > 0))
  x_brks <- pretty(c(x_lims[1], if (is.null(x_break)) x_lims[2] else x_break[1]), n = x_n_breaks)
  if (!is.null(x_break)) x_brks <- c(x_brks[x_brks < x_break[1]],
                                     floor(max(plot_dt$signed_value)))  # one tick past the break
  x_brks <- x_brks[to_plot(x_brks) >= x_lims[1] & to_plot(x_brks) <= x_lims[2]]
  plot_dt[, q_x := fifelse(signed_value >= 0, plot_value + label_nudge, plot_value - label_nudge)]
  plot_dt[, q_hjust := fifelse(signed_value >= 0, 0, 1)]
  break_layers <- if (is.null(x_break)) list() else {
    # Draw the x axis line by hand so it can be interrupted at the break, then add // marks.
    gx <- x_break[1] + c(0.12, gap_w - 0.12); y0 <- 0.4 - pad_rows; tick <- 0.32
    x_ends <- x_lims + c(-1, 1) * 0.02 * diff(x_lims)
    list(annotate('segment', x = c(x_ends[1], gx[2]), xend = c(gx[1], x_ends[2]), y = y0, yend = y0,
                  linewidth = 0.42, color = 'black'),
         annotate('segment', x = gx - 0.1, xend = gx + 0.1, y = y0 - tick, yend = y0 + tick,
                  linewidth = 0.42, color = 'black'),
         # Fixed ranges so the // marks below the panel do not stretch the y range.
         coord_cartesian(xlim = x_ends, ylim = c(0.4 - pad_rows, n_genes + 0.6), expand = FALSE, clip = 'off'),
         theme(axis.line.x = element_blank()))
  }
  p <- ggplot(plot_dt, aes(x = plot_value, y = gene_id, fill = direction)) +
    geom_col(width = 0.78, color = 'black', linewidth = 0.28) +
    geom_text(aes(x = q_x, label = q_label, hjust = q_hjust), size = q_text_size, color = 'black', family = 'sans') +
    geom_vline(xintercept = 0, linewidth = 0.42, color = 'black') +
    scale_fill_manual(values = c('Positive' = '#F04E37', 'Negative' = '#4DB6CC'), guide = 'none') +
    scale_x_continuous(limits = x_lims, breaks = to_plot(x_brks), labels = x_brks, oob = scales::oob_keep,
                       expand = expansion(mult = c(0.02, 0.02))) +
    # pad_rows adds empty rows at the bottom so panels with fewer genes stay top-aligned.
    scale_y_discrete(expand = expansion(add = c(0.6 + pad_rows, 0.6))) +
    labs(title = title, subtitle = subtitle, x = x_label, y = NULL) +
    theme_classic(base_size = 9) +
    theme(axis.line = element_line(linewidth = 0.42, color = 'black'),
          axis.ticks = element_line(linewidth = 0.38, color = 'black'), axis.ticks.length = unit(2.4, 'pt'),
          axis.text.x = element_text(color = 'black', size = 8.5),
          axis.text.y = element_text(color = 'black', size = y_text_size, face = 'italic'),
          axis.title.x = element_text(color = 'black', size = 9.5, margin = margin(t = 8)),
          plot.title = element_text(color = 'black', size = title_size, face = 'bold', hjust = 0),
          plot.subtitle = element_text(color = 'black', size = subtitle_size, hjust = 0, margin = margin(t = 2, b = 5)),
          plot.margin = margin(t = 7, r = 8, b = 8, l = 8), panel.grid = element_blank()) +
    break_layers
  attr(p, 'x_span') <- diff(x_lims)
  p
}

## HVG-EN: 30 largest |coefficient| of the locked model (110 training profiles) ----------
en_model <- models[['HVG-ElasticNet']]
en_coef <- fread(file.path(lock, 'HVG-ElasticNet-exclude-single-age/coefficients.csv'))
stopifnot(all.equal(en_coef$coefficient, en_model$coefficients$coefficient), length(en_model$training_keys) == 110)
en_q <- age_q(version_b_rank(en_model$training_keys))
en_top <- en_coef[coefficient != 0][order(-abs(coefficient), gene)][1:30]
en_top[, `:=`(rank = .I, abs_coefficient = abs(coefficient))]
en_top <- merge(en_top, en_q[, .(gene, age_beta = beta, age_t, residual_df, partial_r2_age, partial_r2_dataset,
                                age_association_p, age_association_q)], by = 'gene', all.x = TRUE, sort = FALSE)
setorder(en_top, rank); setnames(en_top, 'gene', 'gene_id'); setcolorder(en_top, c('rank', 'gene_id'))
fwrite(en_top, file.path(en_dir, 'top30_final_elasticnet_hvg_coefficients_plot_source.tsv'), sep = '\t')
p_en <- plot_top_gene_coefficients(
  en_top, gene_col = 'gene_id', value_col = 'coefficient', x_label = 'Gene coefficient (weeks per ln(1+CPM) unit)',
  title = 'Top HVG–EN genes (single-age excl.)',
  subtitle = sprintf('Locked coefficients, 30 of %d nonzero; q = BH-adjusted age association\n(age + source-batch model, %d training profiles)',
                     sum(en_coef$coefficient != 0), length(en_model$training_keys)),
  symmetric = FALSE, x_break = c(4.8, 11.5))  # Pcdhb9 (~12) is far beyond the next-largest |coefficient| (~4.4)
save_plot(p_en, file.path(en_dir, 'top30_final_elasticnet_hvg_coefficients'), 4.3, 5.0)

## HVG-EN: all nonzero coefficients with q < 0.05, split by direction --------------------
en_sig <- merge(en_coef[coefficient != 0], en_q[, .(gene, age_beta = beta, age_t, residual_df, age_association_p,
                                                     age_association_q)], by = 'gene', all.x = TRUE, sort = FALSE)
en_sig <- en_sig[age_association_q < 0.05][order(-abs(coefficient), gene)]
en_sig[, `:=`(rank = .I, direction = fifelse(coefficient > 0, 'older', 'younger'))]
setnames(en_sig, 'gene', 'gene_id'); setcolorder(en_sig, c('rank', 'gene_id', 'direction'))
fwrite(en_sig, file.path(en_dir, 'significant_final_elasticnet_hvg_coefficients_plot_source.tsv'), sep = '\t')
n_dir <- en_sig[, .N, by = direction][, setNames(N, direction)]
sig_panel <- function(dir, title, ...) {
  plot_top_gene_coefficients(en_sig[direction == dir], gene_col = 'gene_id', value_col = 'coefficient',
                             x_label = NULL, title = title, subtitle = NULL, symmetric = FALSE,
                             pad_rows = max(n_dir) - n_dir[[dir]], y_text_size = 6.6, q_text_size = 1.85,
                             x_n_breaks = 4, by_magnitude = TRUE, ...) +
    theme(plot.title = element_text(size = 8.5, face = 'bold', hjust = 0.5))
}
p_young <- sig_panel('younger', sprintf('Toward younger predicted age (n = %d)', n_dir[['younger']]))
p_old <- sig_panel('older', sprintf('Toward older predicted age (n = %d)', n_dir[['older']]), x_break = c(3.4, 11.5))
p_sig <- patchwork::wrap_plots(p_young, p_old, nrow = 1, widths = c(attr(p_young, 'x_span'), attr(p_old, 'x_span'))) +
  patchwork::plot_annotation(
    title = 'Age-associated HVG–EN genes (single-age excl.)',
    subtitle = sprintf('%d of %d nonzero locked coefficients with q < 0.05 (BH-adjusted age association,\nage + source-batch model, %d training profiles); x = gene coefficient (weeks per ln(1+CPM) unit)',
                       nrow(en_sig), sum(en_coef$coefficient != 0), length(en_model$training_keys)),
    theme = theme(plot.title = element_text(size = 10, face = 'bold'), plot.subtitle = element_text(size = 7)))
save_plot(p_sig, file.path(en_dir, 'significant_final_elasticnet_hvg_coefficients'), 6.2, 5.6)

## HVG-EN: coefficient vs age-association p for all nonzero coefficients -----------------
en_nz <- merge(en_coef[coefficient != 0], en_q[, .(gene, age_t, age_association_p, age_association_q)],
               by = 'gene', all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(en_nz$age_association_p))
en_nz[, `:=`(significant = fifelse(age_association_q < 0.05, 'q < 0.05', 'q ≥ 0.05'),
             sign_agreement = fifelse(sign(coefficient) == sign(age_t), 'Same sign as age t', 'Opposite sign to age t'))]
setorder(en_nz, age_association_p); setnames(en_nz, 'gene', 'gene_id')
fwrite(en_nz, file.path(en_dir, 'nonzero_coefficient_vs_age_association_p_plot_source.tsv'), sep = '\t')
p_cut <- en_q[age_association_q < 0.05, max(age_association_p)]  # BH q = 0.05 boundary across all filtered genes
p_assoc <- ggplot(en_nz, aes(coefficient, -log10(age_association_p))) +
  geom_vline(xintercept = 0, color = 'grey70', linewidth = 0.3) +
  geom_hline(yintercept = -log10(p_cut), linetype = 'dashed', color = 'grey40', linewidth = 0.3) +
  geom_point(aes(color = significant, shape = sign_agreement), size = 2, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = gene_id, color = significant), size = 2.2, max.overlaps = 30,
                           segment.size = 0.2, min.segment.length = 0, seed = 1, show.legend = FALSE) +
  scale_color_manual(values = c('q < 0.05' = '#B2182B', 'q ≥ 0.05' = 'grey55'), name = NULL) +
  scale_shape_manual(values = c('Same sign as age t' = 16, 'Opposite sign to age t' = 17), name = NULL) +
  labs(title = 'HVG–EN coefficient vs age association (single-age excl.)',
       subtitle = sprintf('All %d nonzero locked coefficients; dashed line = BH q = 0.05 (across %d filtered genes)\nage + source-batch model, %d training profiles',
                          nrow(en_nz), nrow(en_q), length(en_model$training_keys)),
       x = 'Gene coefficient (weeks per ln(1+CPM) unit)', y = expression(-log[10]~'age-association'~italic(p))) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = 'bold', size = 10), plot.subtitle = element_text(size = 7, color = 'grey30'),
        panel.grid.minor = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7.5))
save_plot(p_assoc, file.path(en_dir, 'nonzero_coefficient_vs_age_association_p'), 6.5, 5.8)

## Version B-PCA: first 30 genes of the Version B ranking (136 training profiles) --------
pca_model <- models[['VersionB-PCA']]
pca_rank <- version_b_rank(pca_model$training_keys)
# The regenerated ranking must reproduce the locked gene list exactly.
stopifnot(identical(choose_genes(pca_rank, pca_model$hyperparameters$N), pca_model$representation$genes),
          identical(pca_rank$filtered_genes, pca_model$representation$filtered_genes))
pca_q <- age_q(pca_rank)
pca_top <- pca_q[gene %in% pca_model$representation$genes][1:30]
pca_top[, `:=`(version_b_rank = match(gene, pca_q$gene), gene_coefficient = pca_model$coefficients$coefficient[match(gene, pca_model$coefficients$gene)])]
setnames(pca_top, c('gene', 'score', 'beta'), c('gene_id', 'feature_score', 'age_beta'))
setcolorder(pca_top, c('version_b_rank', 'gene_id', 'feature_score', 'age_beta', 'age_t', 'residual_df',
                       'age_association_p', 'age_association_q'))
fwrite(pca_top, file.path(pca_dir, 'top30_final_pca_versionB_gene_effects_plot_source.tsv'), sep = '\t')
p_pca <- plot_top_gene_coefficients(
  pca_top, gene_col = 'gene_id', value_col = 'age_beta', x_label = 'Version B selected-gene age beta (ln(1+CPM) per week)',
  title = 'Top Version B–PCA genes',
  subtitle = sprintf('Top 30 of the Version B ranking (%d training profiles); q = BH-adjusted age association',
                     length(pca_model$training_keys)))
save_plot(p_pca, file.path(pca_dir, 'top30_final_pca_versionB_gene_effects'), 5.0, 5.0)
cat('done\n')
