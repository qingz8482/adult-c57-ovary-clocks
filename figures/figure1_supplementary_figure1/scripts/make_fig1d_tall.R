#!/usr/bin/env Rscript
# Taller-tile layout variant of Figure 1d (inner-CV count selection per stage) for manual figure assembly.
# Same encoding as panels/Fig1d_inner_count.pdf (make_figures.R, panel d), re-drawn from its saved source
# table only (selected-N boxes drawn thicker); the original panel and the Figure 1 composite are not modified.
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
stem <- file.path(fig_dir, 'panels', 'Fig1d_inner_count_tall')

ink <- '#0b0b0b'; ink2 <- '#52514e'
stage_levels <- c('A·f1', 'A·f2', 'A·f3', 'A·f4', 'B', 'C', 'Final')
theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = 'sans') +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = base - 1),
          axis.line = element_line(linewidth = .3, colour = ink2),
          axis.ticks = element_line(linewidth = .3, colour = ink2),
          strip.background = element_blank(), strip.text = element_text(face = 'bold', size = base),
          legend.key.size = unit(3, 'mm'), legend.text = element_text(size = base - 1),
          legend.title = element_text(size = base - 1), plot.title = element_text(size = base + 1, face = 'bold'),
          plot.margin = margin(3, 4, 3, 3))
}

best_n <- fread(file.path(fig_dir, 'source_data', 'Fig1d_inner_cv_by_count.csv'), encoding = 'UTF-8')
stopifnot(nrow(best_n) == 70, best_n[selected == TRUE, .N] == 14)
best_n[, `:=`(Recipe = factor(Recipe, levels = unique(Recipe)), stage = factor(stage, levels = rev(stage_levels)),
              Nf = factor(Nf, levels = c('500', '1k', '2k', '3k', '5k')))]
stopifnot(!anyNA(best_n$stage), !anyNA(best_n$Nf))

pc <- ggplot(best_n, aes(Nf, stage)) +
  geom_tile(aes(fill = pmin(delta, 4)), colour = 'white', linewidth = .5) +
  geom_tile(data = best_n[selected == TRUE], fill = NA, colour = ink, linewidth = 1.2) +
  geom_text(aes(label = sprintf('%.1f', inner_MAE), colour = delta < 1.6), size = 1.7) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = ink), guide = 'none') +
  scale_fill_gradient('Δ inner MAE vs best\n(weeks, capped at 4)', low = '#173f73', high = '#eef3fa', limits = c(0, 4)) +
  facet_wrap(~Recipe, nrow = 1) + labs(x = 'Input gene count N', y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = 'right',
                      legend.key.height = unit(10, 'mm'))

# Original: 183 x 36 mm. Same width, taller canvas so each tile row is about 2.4x as high.
ggsave(paste0(stem, '.pdf'), pc, width = 183, height = 80, units = 'mm', device = cairo_pdf)
ggsave(paste0(stem, '.png'), pc, width = 183, height = 80, units = 'mm', dpi = 300, bg = 'white')
