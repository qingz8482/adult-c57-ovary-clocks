#!/usr/bin/env Rscript
# Two-panel figure: model gene coefficient vs gene-age association (age t, age + source-batch model on each
# model's own training profiles) for (A) HVG-EN nonzero coefficients and (B) Version B-PCA-glmnet gene-level
# weighted coefficients. A second, three-panel figure adds, between them, the distribution of the age t statistic
# across all 3000 HVG-EN input genes with the 79 nonzero genes marked (most HVG coefficients are shrunk to 0, so
# coefficient vs age t over all HVGs is not informative). Reads the plot-source tables written by make_diagnostics.R, make_hvg_coefficient_age_table.R and
# make_pca_coefficient_age_correlation.R; nothing is refitted.
# Format follows 11_shared_nonzero_gene_coefficients (make_shared_gene_coefficients.R): each panel 110 x 115 mm,
# classic axes, base 7 pt, italic gene labels, broken y axis for the Pcdhb9 outlier; font Arial.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
gene_dir <- file.path(fig_dir, 'diagnostics/09_selected_gene_diagnostics')
stopifnot(basename(systemfonts::match_fonts('Arial')$path) == 'arial.ttf')

en <- fread(file.path(gene_dir, 'elasticnet_hvg_coefficients/nonzero_coefficient_vs_age_association_p_plot_source.tsv'), encoding = 'UTF-8')
hv <- fread(file.path(gene_dir, 'elasticnet_hvg_coefficients/all_hvg_coefficient_vs_age_t_plot_source.tsv'), encoding = 'UTF-8')
pc <- fread(file.path(gene_dir, 'pca_versionB_gene_effects/selected_gene_coefficient_vs_age_t_correlation_plot_source.tsv'), encoding = 'UTF-8')
sig_levels <- c('q < 0.05', 'q ≥ 0.05')
for (d in list(en, hv, pc)) d[, significant := factor(significant, levels = sig_levels)]

# ---------------------------------------------------------------- statistics
cor_stats <- function(d, model) {
  sides <- list(all = rep(TRUE, nrow(d)), age_up = d$age_t > 0, age_down = d$age_t < 0)
  rbindlist(lapply(names(sides), function(s) {
    x <- d[sides[[s]]]
    rbindlist(lapply(c('pearson', 'spearman'), function(m) {
      ct <- suppressWarnings(cor.test(x$coefficient, x$age_t, method = m, exact = if (m == 'spearman') FALSE))
      data.table(model = model, genes = s, method = m, estimate = unname(ct$estimate), p_value = ct$p.value, n = nrow(x))
    }))
  }))
}
stats <- rbind(cor_stats(en, 'HVG-EN'), cor_stats(pc, 'VersionB-PCA'))
fwrite(stats, file.path(gene_dir, 'coefficient_vs_age_t_correlation_panels_stats.tsv'), sep = '\t')
st <- function(m, g, k) stats[model == m & genes == g & method == k]
fmt_p <- function(p) if (p < 1e-15) 'P < 1e-15' else sprintf('P = %.1e', p)
cor_line <- function(m) sprintf('Pearson r = %.2f (%s); Spearman ρ = %.2f (%s)',
                                st(m, 'all', 'pearson')$estimate, fmt_p(st(m, 'all', 'pearson')$p_value),
                                st(m, 'all', 'spearman')$estimate, fmt_p(st(m, 'all', 'spearman')$p_value))

# ---------------------------------------------------------------- style (as make_shared_gene_coefficients.R)
font <- 'Arial'
ink <- '#0b0b0b'; ink2 <- '#52514e'; sig_cols <- c('q < 0.05' = '#B2182B', 'q ≥ 0.05' = 'grey55')
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = font) +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = base - 1),
          axis.line = element_line(linewidth = .3, colour = ink2), axis.ticks = element_line(linewidth = .3, colour = ink2),
          legend.key.size = unit(3, 'mm'), legend.text = element_text(size = base - 1), legend.title = element_text(size = base - 1),
          legend.position = 'bottom', plot.title = element_text(size = base, face = 'bold'),
          plot.subtitle = element_text(size = base - 1.5, colour = ink2), plot.margin = margin(3, 4, 3, 3),
          plot.tag = element_text(size = base + 2, face = 'bold'))
}
x_title <- expression(paste('Age-association ', italic(t), ' statistic (age + source-batch model)'))  # paste(): ~ falls back to a non-Arial font
lab <- function(d) ggrepel::geom_text_repel(data = d, aes(label = gene_id, colour = significant), family = font, size = 1.7,
                                            fontface = 'italic', segment.size = .2, segment.colour = ink2, max.overlaps = Inf,
                                            seed = 20261013, show.legend = FALSE)
# Least-squares fit to all genes of a model, drawn from predictions limited to each part of a broken axis.
fit_layers <- function(d, x_lim, ylim) {
  g <- data.table(age_t = seq(x_lim[1], x_lim[2], length.out = 400))
  g <- cbind(g, as.data.table(predict(lm(coefficient ~ age_t, d), g, interval = 'confidence')))
  g[, `:=`(lwr = pmin(pmax(lwr, ylim[1]), ylim[2]), upr = pmin(pmax(upr, ylim[1]), ylim[2]),
           fit = fifelse(fit >= ylim[1] & fit <= ylim[2], fit, NA_real_))]
  if (all(g$lwr == g$upr)) return(NULL)
  list(geom_ribbon(data = g, aes(age_t, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = 'grey85', alpha = .6),
       geom_line(data = g[!is.na(fit)], aes(age_t, fit), inherit.aes = FALSE, colour = ink2, linewidth = .4))
}
# Diagonal break mark on the y axis line at the bottom (0) or top (1) of a broken-axis part.
break_mark <- function(edge) annotation_custom(grid::segmentsGrob(
  x0 = unit(0, 'npc') + unit(-1.2, 'mm'), x1 = unit(0, 'npc') + unit(1.2, 'mm'),
  y0 = unit(edge, 'npc') + unit(-.8, 'mm'), y1 = unit(edge, 'npc') + unit(.8, 'mm'), gp = grid::gpar(col = ink2, lwd = .6)))
base <- function(d_all, d, lab_dt, x_lim, ylim, point_size, point_alpha) ggplot(d, aes(age_t, coefficient)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = ink2) + geom_vline(xintercept = 0, linewidth = .3, colour = ink2) +
  fit_layers(d_all, x_lim, ylim) +
  geom_point(aes(colour = significant), size = point_size, alpha = point_alpha, stroke = 0) + lab(lab_dt) +
  scale_colour_manual(values = sig_cols, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(size = 1.8, alpha = 1))) +
  scale_x_continuous(limits = x_lim, expand = c(0, 0)) + coord_cartesian(ylim = ylim, expand = FALSE, clip = 'off') + theme_pub()
pad_range <- function(v, f) range(v) + c(-1, 1) * f * diff(range(v))

# ---------------------------------------------------------------- HVG-EN panels, broken y axis
# Main range for all genes but Pcdhb9, a short upper strip for Pcdhb9; same units per mm in both parts and one
# fit to all genes. pt_size/pt_alpha may be per-gene vectors (e.g. small points for zero coefficients).
broken_panel <- function(d, lab_dt, tag, title, subtitle, pt_size, pt_alpha) {
  d <- copy(d)[, `:=`(pt_size = pt_size, pt_alpha = pt_alpha)]
  hi <- d[coefficient > 6]; lo <- d[coefficient <= 6]; stopifnot(identical(hi$gene_id, 'Pcdhb9'))
  x_lim <- pad_range(d$age_t, .04)
  y_lo <- range(lo$coefficient) + c(-.5, .5); y_hi <- hi$coefficient + c(-.6, .6)
  pts <- function(x) geom_point(data = x, aes(colour = significant), size = x$pt_size, alpha = x$pt_alpha, stroke = 0)
  p_hi <- base(d, hi[0], lab_dt[gene_id %in% hi$gene_id], x_lim, y_hi, 0, 0) + pts(hi) +
    scale_y_continuous(breaks = 12) + break_mark(0) +
    labs(x = NULL, y = NULL, colour = NULL, tag = tag, title = title,
         subtitle = sub('<BREAK>', ceiling(y_lo[2]), subtitle, fixed = TRUE)) +
    theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = 'none')
  p_lo <- base(d, lo[0], lab_dt[gene_id %in% lo$gene_id], x_lim, y_lo, 0, 0) + pts(lo) + break_mark(1) +
    labs(x = x_title, y = 'HVG–EN coefficient (weeks per ln(1+CPM))', colour = NULL)
  p_hi / p_lo + plot_layout(heights = c(diff(y_hi), diff(y_lo)))
}

# (A) HVG-EN nonzero coefficients, all 79 labelled.
p_en <- broken_panel(en, en, 'A', sprintf('HVG–EN (single-age excl.): %d nonzero genes', nrow(en)),
                     paste0(cor_line('HVG-EN'), sprintf('\nLine: least-squares fit to all %d genes with 95%% CI; y axis broken above <BREAK>', nrow(en))),
                     1.4, 1)
# Age t across all 3000 HVG-EN input genes (same age + source-batch model, 110 training profiles); rug = the 79
# nonzero genes; dashed = smallest |t| reaching BH q < 0.05 (BH across all filtered genes, as in the other panels).
stopifnot(sum(hv$nonzero) == nrow(en), setequal(hv[nonzero == TRUE, gene_id], en$gene_id))
t_cut <- hv[age_association_q < 0.05, min(abs(age_t))]
p_hv <- ggplot(hv, aes(age_t)) +
  geom_histogram(binwidth = .5, boundary = 0, fill = 'grey80', colour = 'white', linewidth = .15) +
  geom_vline(xintercept = 0, linewidth = .3, colour = ink2) +
  geom_vline(xintercept = c(-1, 1) * t_cut, linewidth = .3, colour = ink2, linetype = 'dashed') +
  geom_rug(data = hv[nonzero == TRUE], aes(colour = significant), sides = 'b', length = unit(2.5, 'mm'),
           linewidth = .35, alpha = .9) +
  scale_colour_manual(values = sig_cols, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) +
  scale_x_continuous(limits = pad_range(hv$age_t, .04), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  coord_cartesian(clip = 'off') + theme_pub() +
  labs(x = x_title, y = 'Number of HVG input genes', colour = NULL, tag = 'B',
       title = sprintf('HVG–EN input: age association of all %d HVGs', nrow(hv)),
       subtitle = sprintf('%d of %d HVGs with q < 0.05 (dashed: |t| = %.2f); rug: %d nonzero HVG–EN genes (%d with q < 0.05)\n%d HVGs (%.0f%%) shrunk to coefficient 0 by the elastic net',
                          hv[age_association_q < 0.05, .N], nrow(hv), t_cut, sum(hv$nonzero),
                          hv[nonzero == TRUE & age_association_q < 0.05, .N], sum(!hv$nonzero), 100 * mean(!hv$nonzero)))

# ---------------------------------------------------------------- (B) Version B-PCA
# Label the 10 largest positive and 10 largest negative gene coefficients.
pc_y <- pad_range(pc$coefficient, .04)
p_pc <- base(pc, pc, rbind(pc[order(-coefficient)][1:10], pc[order(coefficient)][1:10]),
             pad_range(pc$age_t, .04), pc_y, .5, .5) +
  labs(x = x_title, y = 'Version B–PCA weighted gene coefficient\n(rotation × PC coefficients; weeks per ln(1+CPM))',
       colour = NULL, tag = 'B', title = sprintf('Version B–PCA: %d selected genes', nrow(pc)),
       subtitle = sprintf('%s; within age-up / age-down ρ = %.2f / %.2f\nLine: least-squares fit to all %d genes with 95%% CI; genes pre-selected by age association',
                          cor_line('VersionB-PCA'), st('VersionB-PCA', 'age_up', 'spearman')$estimate,
                          st('VersionB-PCA', 'age_down', 'spearman')$estimate, nrow(pc)))

# Each panel 110 x 115 mm, as in 11_shared_nonzero_gene_coefficients; each keeps its own bottom legend.
save_fig <- function(fig, name, n_panels) {
  stem <- file.path(gene_dir, name)
  ggsave(paste0(stem, '.pdf'), fig, width = 110 * n_panels, height = 115, units = 'mm', device = cairo_pdf)
  ggsave(paste0(stem, '.png'), fig, width = 110 * n_panels, height = 115, units = 'mm', dpi = 300, bg = 'white', device = ragg::agg_png)
}
save_fig(p_en | (p_pc + labs(tag = 'B')), 'coefficient_vs_age_t_correlation_panels', 2)
save_fig(p_en | p_hv | (p_pc + labs(tag = 'C')), 'coefficient_vs_age_t_correlation_panels_with_all_hvg', 3)
print(stats)
cat('done\n')
