#!/usr/bin/env Rscript
# Coefficient agreement for the HVG-ElasticNet nonzero genes that are also Version B-PCA input genes
# (36 of 79). x: Version B-PCA effective gene coefficient (rotation %*% PC coefficients, weeks per SD);
# y: HVG-EN coefficient (weeks per SD). The y axis is broken so the Pcdhb9 outlier does not compress the
# other 35 genes; both panels share one y scale and the fit uses all 36 genes. Reads 10_model_gene_overlap outputs only; nothing is refitted.
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
diag_dir <- file.path(fig_dir, 'diagnostics')
stem <- file.path(diag_dir, '11_shared_nonzero_gene_coefficients')

ink <- '#0b0b0b'; ink2 <- '#52514e'; en_col <- '#2a78d6'; off_col <- '#d6542a'
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = 'sans') +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = base - 1),
          axis.line = element_line(linewidth = .3, colour = ink2), axis.ticks = element_line(linewidth = .3, colour = ink2),
          legend.key.size = unit(3, 'mm'), legend.text = element_text(size = base - 1), legend.title = element_text(size = base - 1),
          legend.position = 'bottom', plot.title = element_text(size = base, face = 'bold'),
          plot.subtitle = element_text(size = base - 1.5, colour = ink2), plot.margin = margin(3, 4, 3, 3))
}

loc <- fread(file.path(diag_dir, '10_model_gene_overlap', 'c_en_nonzero_genes_in_versionB.tsv'))
sh <- loc[in_vb_pool == TRUE]
stopifnot(nrow(loc) == 79, nrow(sh) == 36, !anyNA(sh$vb_coef), all(sh$en_coef != 0))
sh[, Sign := factor(fifelse(same_sign, 'Same sign', 'Opposite sign'), levels = c('Same sign', 'Opposite sign'))]

pe <- cor.test(sh$vb_coef, sh$en_coef, method = 'pearson')
sp <- cor.test(sh$vb_coef, sh$en_coef, method = 'spearman', exact = FALSE)
sg <- binom.test(sum(sh$same_sign), nrow(sh), .5)
fwrite(sh[, .(gene, en_coef, en_rank, vb_coef, vb_rank, same_sign)][order(en_rank)], paste0(stem, '.tsv'), sep = '\t')
fwrite(data.table(metric = c('n_genes', 'pearson_r', 'pearson_P', 'spearman_rho', 'spearman_P', 'same_sign', 'same_sign_binomial_P'),
                  value = c(nrow(sh), pe$estimate, pe$p.value, sp$estimate, sp$p.value, sum(sh$same_sign), sg$p.value)),
       paste0(stem, '_stats.tsv'), sep = '\t')

lab <- function(d) if (requireNamespace('ggrepel', quietly = TRUE))
  ggrepel::geom_text_repel(data = d, aes(label = gene, colour = Sign), size = 1.7, fontface = 'italic', segment.size = .2,
                           segment.colour = ink2, max.overlaps = Inf, seed = 20261013, show.legend = FALSE) else
  geom_text(data = d, aes(label = gene, colour = Sign), size = 1.6, fontface = 'italic', vjust = -.7, show.legend = FALSE)
# Broken y axis: main range for 35 genes, a short upper strip for the outlier, same units per mm in both.
hi <- sh[en_coef > 6]; lo <- sh[en_coef <= 6]; stopifnot(identical(hi$gene, 'Pcdhb9'))
y_lo <- range(lo$en_coef) + c(-.5, .5); y_hi <- hi$en_coef + c(-.6, .6)
x_lim <- range(sh$vb_coef) + c(-1, 1) * .04 * diff(range(sh$vb_coef))
# Least-squares fit to all 36 genes, drawn from predictions limited to each panel (clip is off for the break marks).
fit <- lm(en_coef ~ vb_coef, sh); fg <- data.table(vb_coef = seq(x_lim[1], x_lim[2], length.out = 400))
fg <- cbind(fg, as.data.table(predict(fit, fg, interval = 'confidence')))
fit_layers <- function(ylim) {
  g <- copy(fg)[, `:=`(lwr = pmin(pmax(lwr, ylim[1]), ylim[2]), upr = pmin(pmax(upr, ylim[1]), ylim[2]),
                       fit = fifelse(fit >= ylim[1] & fit <= ylim[2], fit, NA_real_))]
  if (all(g$lwr == g$upr)) return(NULL)
  list(geom_ribbon(data = g, aes(vb_coef, ymin = lwr, ymax = upr), inherit.aes = FALSE, fill = 'grey85', alpha = .6),
       geom_line(data = g[!is.na(fit)], aes(vb_coef, fit), inherit.aes = FALSE, colour = ink2, linewidth = .4))
}
base <- function(d, ylim) ggplot(sh, aes(vb_coef, en_coef)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = ink2) + geom_vline(xintercept = 0, linewidth = .3, colour = ink2) +
  fit_layers(ylim) +
  geom_point(data = d, aes(colour = Sign), size = 1.4) + lab(d) +
  scale_colour_manual(values = c('Same sign' = en_col, 'Opposite sign' = off_col), drop = FALSE) +
  scale_x_continuous(limits = x_lim, expand = c(0, 0)) + coord_cartesian(ylim = ylim, expand = FALSE, clip = 'off') + theme_pub()
p_hi <- base(hi, y_hi) + scale_y_continuous(breaks = 10) +
  annotation_custom(grid::segmentsGrob(x0 = unit(0, 'npc') + unit(-1.2, 'mm'), x1 = unit(0, 'npc') + unit(1.2, 'mm'),
                                       y0 = unit(0, 'npc') + unit(-.8, 'mm'), y1 = unit(0, 'npc') + unit(.8, 'mm'), gp = grid::gpar(col = ink2, lwd = .6))) +
  labs(x = NULL, y = NULL, colour = NULL,
       title = sprintf('HVG–EN nonzero genes in the Version B–PCA input (%d of 79)', nrow(sh)),
       subtitle = sprintf('Pearson r = %.2f (P = %.1e); Spearman ρ = %.2f (P = %.1e); same sign %d/%d (binomial P = %.1e)\nLine: least-squares fit to all %d genes with 95%% CI; y axis broken above %d',
                          pe$estimate, pe$p.value, sp$estimate, sp$p.value, sum(sh$same_sign), nrow(sh), sg$p.value, nrow(sh), ceiling(y_lo[2]))) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = 'none')
p_lo <- base(lo, y_lo) +
  annotation_custom(grid::segmentsGrob(x0 = unit(0, 'npc') + unit(-1.2, 'mm'), x1 = unit(0, 'npc') + unit(1.2, 'mm'),
                                       y0 = unit(1, 'npc') + unit(-.8, 'mm'), y1 = unit(1, 'npc') + unit(.8, 'mm'), gp = grid::gpar(col = ink2, lwd = .6))) +
  labs(x = 'Version B–PCA weighted gene coefficient\n(rotation × PC coefficients; weeks per SD)',
       y = 'HVG–EN coefficient (weeks per SD)', colour = NULL)
p <- p_hi / p_lo + plot_layout(heights = c(diff(y_hi), diff(y_lo)))
ggsave(paste0(stem, '.pdf'), p, width = 110, height = 115, units = 'mm', device = cairo_pdf)
ggsave(paste0(stem, '.png'), p, width = 110, height = 115, units = 'mm', dpi = 300, bg = 'white')
cat(sprintf('n = %d; Pearson r = %.3f (P = %.2e); Spearman rho = %.3f (P = %.2e); same sign %d/%d\n',
            nrow(sh), pe$estimate, pe$p.value, sp$estimate, sp$p.value, sum(sh$same_sign), nrow(sh)))
