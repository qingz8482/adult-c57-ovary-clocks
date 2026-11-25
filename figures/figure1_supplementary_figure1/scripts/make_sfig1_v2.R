#!/usr/bin/env Rscript
# Supplementary Figure 1, version 2: the final-fit lambda path (old panel a) and the per-stage
# lambda paths (old panel c) are merged into one panel with Final as the last column.
# Runs make_figures.R up to its Supplementary Figure assembly with all file writes disabled,
# so the existing Figure 1, Supplementary Figure 1, panels and source data are left untouched.
# Writes only *_v2 / SFig1v2_* files. Nothing is refitted.

# Mask writers before make_figures.R runs; global definitions take precedence over attached packages.
ggsave <- function(...) invisible(NULL)
fwrite <- function(...) invisible(NULL)
dir.create <- function(...) invisible(TRUE)

a <- commandArgs(FALSE)
this <- normalizePath(gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE))
src <- readLines(file.path(dirname(this), 'make_figures.R'))
stop_at <- grep('^sfig <- ', src)
stopifnot(length(stop_at) == 1)
eval(parse(text = src[seq_len(stop_at - 1)]), envir = globalenv())
rm(ggsave, fwrite, dir.create)

## a | lambda paths of the selected candidate: every outer stage plus the final fit
fin <- copy(lam_final)[, `:=`(stage = 'Final', sel = sel_lambda)][, sel_lambda := NULL]
lam_all <- rbind(lam_st, fin, use.names = TRUE)
lam_all[, stage := factor(stage, levels = stage_levels)]
lm_all <- lam_all[, .(sel = sel[1], minl = lambda[which.min(mean_mae)], thr = min(mean_mae) + se[which.min(mean_mae)]), by = .(Recipe, stage)]
stopifnot(nrow(lm_all) == 14, all(abs(lm_all[stage == 'Final', thr] - lam_mark[order(Recipe), thr]) < 1e-9))
fin_bg <- unique(lam_all[stage == 'Final', .(Recipe, stage)])
data.table::fwrite(lam_all[, .(Recipe, recipe, stage, N, pc, alpha, lambda, mean_mae, se, n_studies, selected_lambda = sel)],
                   file.path(out_src, 'SFig1v2a_lambda_paths_all_fits.csv'))
sl2 <- ggplot(lam_all, aes(log10(lambda), mean_mae)) +
  geom_rect(data = fin_bg, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), inherit.aes = FALSE, fill = '#f1f0ea') +
  geom_ribbon(aes(ymin = mean_mae - se, ymax = mean_mae + se, fill = Recipe), alpha = .18) +
  geom_line(aes(colour = Recipe), linewidth = .5) +
  geom_hline(data = lm_all, aes(yintercept = thr), linetype = 'dotted', linewidth = .3, colour = ink2) +
  geom_vline(data = lm_all, aes(xintercept = log10(minl)), linetype = 'dashed', linewidth = .3, colour = ink2) +
  geom_vline(data = lm_all, aes(xintercept = log10(sel)), linewidth = .4, colour = '#b2182b') +
  scale_colour_manual(values = rec_col, guide = 'none') + scale_fill_manual(values = rec_col, guide = 'none') +
  scale_x_continuous(expression(log[10]~lambda), breaks = c(-4, -2, 0, 2)) +
  facet_grid(Recipe ~ stage, scales = 'free_y', labeller = labeller(Recipe = function(x) sub(' \\(', '\n(', x))) +
  labs(y = 'Inner LOSO MAE (weeks)') +
  theme_pub(6.5) + theme(panel.border = element_rect(fill = NA, colour = grid_col, linewidth = .3), strip.text.y = element_text(size = 5.5))

sfig2 <- tag(sl2, 'a') / tag(sa, 'b') / tag(sc, 'c') / tag(se, 'd') / (tag(sd, 'e') | tag(sr, 'f')) /
  (tag(sf, 'g') | tag(sg, 'h')) / tag(sh, 'i') +
  plot_layout(heights = c(.95, .85, .52, .72, .85, .55, .45))
ggplot2::ggsave(file.path(out_panels, 'SFig1v2a_lambda_all_fits.pdf'), sl2, width = 183, height = 60, units = 'mm', device = cairo_pdf)
W2 <- 183; H2 <- 268
ggplot2::ggsave(file.path(out_final, 'Supplementary_Figure1_model_development_v2.pdf'), sfig2, width = W2, height = H2, units = 'mm', device = cairo_pdf)
ggplot2::ggsave(file.path(out_final, 'Supplementary_Figure1_model_development_v2.png'), sfig2, width = W2, height = H2, units = 'mm', dpi = 300)
cat('done\n')
