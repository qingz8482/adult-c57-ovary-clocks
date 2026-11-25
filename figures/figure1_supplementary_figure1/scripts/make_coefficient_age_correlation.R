#!/usr/bin/env Rscript
# Correlation between the locked HVG-EN nonzero coefficients and each gene's signed age association
# (age t from the age + source-batch model). Reads the plot-source table written by make_diagnostics.R;
# nothing is refitted.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
en_dir <- file.path(fig_dir, 'diagnostics/09_selected_gene_diagnostics/elasticnet_hvg_coefficients')
en_nz <- fread(file.path(en_dir, 'nonzero_coefficient_vs_age_association_p_plot_source.tsv'), encoding = 'UTF-8')

save_plot <- function(p, stem, width, height) {
  ggsave(paste0(stem, '.png'), p, width = width, height = height, dpi = 300, bg = 'white', limitsize = FALSE)
  ggsave(paste0(stem, '.pdf'), p, width = width, height = height, bg = 'white', limitsize = FALSE, device = cairo_pdf)
}

ct_p <- cor.test(en_nz$coefficient, en_nz$age_t, method = 'pearson')
ct_s <- suppressWarnings(cor.test(en_nz$coefficient, en_nz$age_t, method = 'spearman', exact = FALSE))
n_same <- en_nz[sign_agreement == 'Same sign as age t', .N]
fmt_p <- function(p) if (p < 1e-15) 'p < 1e-15' else sprintf('p = %s', format(signif(p, 2)))
stats <- data.table(term = c('pearson', 'spearman', 'sign_agreement'),
                    estimate = c(ct_p$estimate, ct_s$estimate, n_same / nrow(en_nz)),
                    p_value = c(ct_p$p.value, ct_s$p.value, binom.test(n_same, nrow(en_nz))$p.value),
                    n = nrow(en_nz))
fwrite(stats, file.path(en_dir, 'nonzero_coefficient_vs_age_t_correlation_stats.tsv'), sep = '\t')
lab <- sprintf('Pearson r = %.2f, %s\nSpearman ρ = %.2f, %s',
               ct_p$estimate, fmt_p(ct_p$p.value), ct_s$estimate, fmt_p(ct_s$p.value))

p_cor <- ggplot(en_nz, aes(age_t, coefficient)) +
  geom_hline(yintercept = 0, color = 'grey70', linewidth = 0.3) +
  geom_vline(xintercept = 0, color = 'grey70', linewidth = 0.3) +
  geom_smooth(method = 'lm', formula = y ~ x, color = 'grey30', fill = 'grey80', linewidth = 0.5) +
  geom_point(aes(color = significant), size = 2, alpha = 0.85) +
  ggrepel::geom_text_repel(aes(label = gene_id, color = significant), size = 2.2, max.overlaps = 30,
                           segment.size = 0.2, min.segment.length = 0, seed = 1, show.legend = FALSE) +
  annotate('label', x = -Inf, y = Inf, label = lab, hjust = 0, vjust = 1, size = 2.6,
           fill = alpha('white', 0.8)) +
  scale_color_manual(values = c('q < 0.05' = '#B2182B', 'q ≥ 0.05' = 'grey55'), name = NULL) +
  labs(title = 'HVG–EN coefficient vs gene–age association (single-age excl.)',
       subtitle = sprintf('All %d nonzero locked coefficients; x = age t statistic (age + source-batch model)\nline = least-squares fit with 95%% CI', nrow(en_nz)),
       x = expression('Age-association'~italic(t)~'statistic'), y = 'Gene coefficient (weeks per ln(1+CPM) unit)') +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = 'bold', size = 10), plot.subtitle = element_text(size = 7, color = 'grey30'),
        panel.grid.minor = element_blank(), legend.position = 'bottom', legend.text = element_text(size = 7.5))
save_plot(p_cor, file.path(en_dir, 'nonzero_coefficient_vs_age_t_correlation'), 6.5, 5.8)
print(stats)
cat('done\n')
