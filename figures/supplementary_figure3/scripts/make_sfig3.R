#!/usr/bin/env Rscript
# Supplementary Figure 3 (plus two main-figure panels): complete comparison of the locked ovarian clocks (rebutal4 release
# hvg_exclusion_vbpca_20261012) with the two tAge variants. Reads the frozen outputs of
# ../compare_tAge only (sample predictions, contrast definitions, gene weights, model-resolved
# pathway statistics); no clock or tAge model is refitted. New statistics computed here:
# per-clock Welch tests per contrast, agreement sensitivities (bulk only, training-study
# exclusion, leave-one-dataset-out) and sample-level agreement of target shifts.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(ggrepel)
})
setDTthreads(1)

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
cmp <- file.path(dirname(fig_dir), 'compare_tAge')
tab <- function(f) fread(file.path(cmp, 'tables', f), na.strings = c('NA', ''))
out_panels <- file.path(fig_dir, 'panels'); out_final <- file.path(fig_dir, 'final')
out_src <- file.path(fig_dir, 'source_data')
for (d in c(out_panels, out_final, out_src)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- constants
models <- c('HVG-ElasticNet', 'VersionB-PCA/glmnet', 'tAge scaled-diff', 'tAge YuGene-diff')
mlab <- setNames(c('HVG–EN', 'Version B–PCA', 'tAge scaled-diff', 'tAge YuGene-diff'), models)
ovary <- models[1:2]; tage <- models[3:4]
col_lab <- setNames(c('HVG\nEN', 'VB\nPCA', 'tAge\nscaled', 'tAge\nYuGene'), mlab)
type_col <- c(diet = '#F28E2B', genetic_model = '#E15759', intervention = '#4E79A7',
              microbiome = '#59A14F', stress_model = '#9C755F')
type_lab <- c(diet = 'Diet', genetic_model = 'Genetic model', intervention = 'Intervention',
              microbiome = 'Microbiome', stress_model = 'Stress model')
ink <- '#0b0b0b'; ink2 <- '#52514e'; grid_col <- '#e6e5e0'
div_lo <- '#2166ac'; div_hi <- '#b2182b'

theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = 'sans') +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = base - 1),
          axis.line = element_line(linewidth = .3, colour = ink2),
          axis.ticks = element_line(linewidth = .3, colour = ink2),
          strip.background = element_blank(), strip.text = element_text(face = 'bold', size = base),
          legend.key.size = unit(3, 'mm'), legend.text = element_text(size = base - 1),
          legend.title = element_text(size = base - 1), plot.title = element_text(size = base + 1, face = 'bold'),
          plot.subtitle = element_text(size = base - 1, colour = ink2), plot.margin = margin(3, 4, 3, 3))
}
save_panel <- function(p, name, w, h) {
  ggsave(file.path(out_panels, paste0(name, '.pdf')), p, width = w, height = h, units = 'mm', device = cairo_pdf)
}
tag <- function(p, t) wrap_elements(full = p) + labs(tag = t) +
  theme(plot.tag = element_text(size = 9, face = 'bold'), plot.tag.position = c(0, 1), plot.margin = margin(2, 1, 1, 9))
stars <- function(p) fcase(is.na(p), '', p < 1e-4, '****', p < 1e-3, '***', p < .01, '**', p < .05, '*', default = '')
cor_row <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  data.table(n = length(x), pearson = cor(x, y), spearman = cor(x, y, method = 'spearman'),
             pearson_p = cor.test(x, y)$p.value,
             spearman_p = suppressWarnings(cor.test(x, y, method = 'spearman', exact = FALSE))$p.value,
             direction_agreement = mean(sign(x) == sign(y)))
}
pairs <- as.data.table(t(combn(models, 2))); setnames(pairs, c('model_1', 'model_2'))
pairs[, pair := paste(mlab[model_1], 'vs', mlab[model_2])]

# ---------------------------------------------------------------- inputs
defs <- tab('contrast_definitions.tsv')
pred <- tab('sample_contrast_predictions.tsv')
eff <- tab('four_clock_contrast_effects.tsv')
agree <- tab('clock_pairwise_agreement.tsv')
cover <- tab('tage_gene_coverage.tsv')
prim <- defs[primary == TRUE, analysis_contrast_id]
stopifnot(length(prim) == 38, setequal(eff[primary == TRUE, analysis_contrast_id], prim),
          pred[, .N, by = model][, all(N == N[1])])

# Short, unique contrast labels; modality is appended where a study has both.
defs[, label := gsub(' versus ', ' vs ', contrast_name)]
defs[, label := sub('available sesame-oil vehicle|sesame-oil vehicle', 'vehicle', label)]
defs[, label := sub('oocyte-conditional knockout', 'oocyte cKO', label)]
defs[, label := sub('Cd38 knockout', 'Cd38 KO', label)]
defs[, label := sub('Rnf111 heterozygous KO', 'Rnf111 het KO', label)]
defs[, label := sub('eight vs two mice per cage', '8 vs 2 mice per cage', label)]
defs[, label := sub(' at (\\d+) months', ' \\1m', label)]
defs[, label := sub('_age', ' age', label)]
defs[, label := sub(' rice vs polished rice', ' vs polished', label)]
defs[, label := sub('white brown', 'brown', label)]
defs[, label := sub('estropausal-donor vs young-donor FMT', 'estropausal vs young FMT', label)]
defs[dataset_folder == 'GSE235021', label := paste0(label, fifelse(modality == 'bulk RNA-seq', ' (bulk)', ' (pb)'))]
stopifnot(!anyDuplicated(defs[primary == TRUE, label]))

# Training overlap: any sample of a study used by either locked clock (controls only in practice).
ovl <- pred[model == 'HVG-ElasticNet', .(training_study = any(locked_training_study),
                                         training_profiles = sum(locked_training_sample)), by = analysis_contrast_id]
defs <- merge(defs, ovl, by = 'analysis_contrast_id', all.x = TRUE)
cov_c <- cover[, .(tage_coverage = min(min_sample_coverage)), by = analysis_contrast_id]
defs <- merge(defs, cov_c, by = 'analysis_contrast_id', all.x = TRUE)

# Effects (target mean − matched-control mean, weeks) and per-clock Welch tests on the same values.
welch <- pred[primary == TRUE, {
  x <- value_weeks[role == 'experimental']; y <- value_weeks[role == 'control']
  .(effect_weeks = mean(x) - mean(y), welch_p = t.test(x, y)$p.value)
}, by = .(analysis_contrast_id, model)]
welch[, welch_bh_q := p.adjust(welch_p, 'BH'), by = model]
chk <- melt(eff[primary == TRUE], id.vars = 'analysis_contrast_id', measure.vars = models,
            variable.name = 'model', value.name = 'saved_effect')
chk <- merge(chk, welch, by = c('analysis_contrast_id', 'model'))
stopifnot(max(abs(chk$saved_effect - chk$effect_weeks)) < 1e-8)
wide <- dcast(welch, analysis_contrast_id ~ model, value.var = 'effect_weeks')
wide <- merge(wide, defs[, .(analysis_contrast_id, dataset_folder, modality, contrast_type, label,
                             training_study, training_profiles, tage_coverage)], by = 'analysis_contrast_id')

# ---------------------------------------------------------------- a | contrast-effect scatter
pa_dat <- rbindlist(lapply(ovary, function(o) rbindlist(lapply(tage, function(t)
  wide[, .(analysis_contrast_id, label, contrast_type, modality, ovary_model = o, tage_model = t,
           ovary_effect = get(o), tage_effect = get(t))]))))
pa_stats <- pa_dat[, cor_row(tage_effect, ovary_effect), by = .(ovary_model, tage_model)]
saved <- agree[analysis_set == 'primary_adult_prespecified']
for (i in seq_len(nrow(pa_stats))) {
  s <- saved[(model_1 == pa_stats$ovary_model[i] & model_2 == pa_stats$tage_model[i]) |
               (model_2 == pa_stats$ovary_model[i] & model_1 == pa_stats$tage_model[i])]
  stopifnot(nrow(s) == 1, abs(s$pearson - pa_stats$pearson[i]) < 1e-10, abs(s$spearman - pa_stats$spearman[i]) < 1e-10)
}
slope <- pa_dat[, .(ols_slope = coef(lm(ovary_effect ~ tage_effect))[2]), by = .(ovary_model, tage_model)]
pa_stats <- merge(pa_stats, slope, by = c('ovary_model', 'tage_model'))
pa_stats[, lab := sprintf('r = %.2f%s\nρ = %.2f%s\nsign = %d/%d', pearson, stars(pearson_p), spearman,
                          stars(spearman_p), round(direction_agreement * n), n)]
for (d in list(pa_dat, pa_stats)) d[, `:=`(ovary_model = factor(mlab[ovary_model], mlab[ovary]),
                                           tage_model = factor(mlab[tage_model], mlab[tage]))]
pa_dat[, contrast_type := factor(contrast_type, names(type_col))]
pa <- ggplot(pa_dat, aes(tage_effect, ovary_effect)) +
  geom_hline(yintercept = 0, colour = grid_col, linewidth = .3) + geom_vline(xintercept = 0, colour = grid_col, linewidth = .3) +
  geom_smooth(method = 'lm', formula = y ~ x, se = TRUE, colour = ink2, fill = grid_col, linewidth = .35) +
  geom_point(aes(fill = contrast_type, shape = modality), size = 1.5, stroke = .25, colour = 'white') +
  geom_text(data = pa_stats, aes(x = -Inf, y = Inf, label = lab), hjust = -.05, vjust = 1.1, size = 1.9,
            lineheight = .9, colour = ink) +
  facet_grid(ovary_model ~ tage_model, scales = 'free', switch = 'y') +
  scale_fill_manual(values = type_col, labels = type_lab, name = NULL) +
  scale_shape_manual(values = c('bulk RNA-seq' = 21, 'scRNA-derived pseudobulk' = 24),
                     labels = c('Bulk', 'Pseudobulk'), name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 2), nrow = 2),
         shape = guide_legend(override.aes = list(fill = ink2), nrow = 2)) +
  labs(x = 'tAge contrast effect (weeks)', y = 'Ovary clock contrast effect (weeks)',
       title = 'Contrast effects, 38 primary contrasts') +
  theme_pub() + theme(strip.placement = 'outside', legend.position = 'bottom', legend.box = 'horizontal',
                      legend.margin = margin(0, 0, 0, 0), panel.spacing = unit(2.5, 'mm'))
fwrite(pa_dat, file.path(out_src, 'Main_contrast_effects.csv'))
fwrite(pa_stats[, !'lab'], file.path(out_src, 'Main_contrast_effect_agreement.csv'))

# ---------------------------------------------------------------- b | robustness of agreement
agree_set <- function(ids, set) {
  w <- wide[analysis_contrast_id %in% ids]
  pairs[, cor_row(w[[model_1]], w[[model_2]]), by = .(model_1, model_2, pair)][, set := set]
}
all_eff <- dcast(pred[analysis_family == 'prespecified', .(e = mean(value_weeks[role == 'experimental']) -
                     mean(value_weeks[role == 'control'])), by = .(analysis_contrast_id, model)],
                 analysis_contrast_id ~ model, value.var = 'e')
agree_ext <- function(ids, set) {
  w <- all_eff[analysis_contrast_id %in% ids]
  pairs[, cor_row(w[[model_1]], w[[model_2]]), by = .(model_1, model_2, pair)][, set := set]
}
nonplace <- defs[analysis_family == 'prespecified' & placement_only == FALSE, analysis_contrast_id]
allpre <- defs[analysis_family == 'prespecified', analysis_contrast_id]
sets <- rbind(
  agree_set(prim, 'Primary'),
  agree_ext(nonplace, 'All non-placement'),
  agree_ext(allpre, 'All prespecified'),
  agree_set(wide[modality == 'bulk RNA-seq', analysis_contrast_id], 'Bulk only'),
  agree_set(wide[training_study == FALSE, analysis_contrast_id], 'No training-study overlap'))
# Saved sensitivity sets must reproduce.
for (s in list(c('All non-placement', 'all_prespecified_nonplacement'), c('All prespecified', 'all_prespecified'))) {
  m <- merge(sets[set == s[1]], agree[analysis_set == s[2]], by = c('model_1', 'model_2'))
  stopifnot(nrow(m) == 6, max(abs(m$pearson.x - m$pearson.y)) < 1e-10)
}
# Sample-level: control-centred shift of each target profile within its contrast.
tg <- dcast(pred[primary == TRUE & role == 'experimental'], analysis_contrast_id + sample_key ~ model,
            value.var = 'relative_shift_weeks')
samp <- pairs[, cor_row(tg[[model_1]], tg[[model_2]]), by = .(model_1, model_2, pair)][, set := 'Sample level (targets)']
sets <- rbind(sets, samp)
# Leave-one-dataset-out on the primary set.
loo <- rbindlist(lapply(unique(wide$dataset_folder), function(ds)
  agree_set(wide[dataset_folder != ds, analysis_contrast_id], ds)))
loo_rng <- loo[, .(pearson_min = min(pearson), pearson_max = max(pearson), spearman_min = min(spearman),
                   spearman_max = max(spearman), dataset_min_pearson = set[which.min(pearson)],
                   dataset_min_spearman = set[which.min(spearman)]), by = .(model_1, model_2, pair)]
set_levels <- c('Primary', 'All non-placement', 'All prespecified', 'Bulk only', 'No training-study overlap', 'Sample level (targets)')
set_n <- sets[, .(n = n[1]), by = set]
set_lab <- setNames(sprintf('%s (%d)', set_n$set, set_n$n), set_n$set)
set_col <- setNames(c(ink, '#8c8c8c', '#c2c2c2', '#e69f00', '#009e73', '#cc79a7'), set_levels)
set_shape <- setNames(c(16, 1, 1, 17, 15, 18), set_levels)
pb_long <- melt(sets, id.vars = c('pair', 'set'), measure.vars = c('pearson', 'spearman'), variable.name = 'metric')
lr_long <- rbind(loo_rng[, .(pair, metric = 'pearson', lo = pearson_min, hi = pearson_max)],
                 loo_rng[, .(pair, metric = 'spearman', lo = spearman_min, hi = spearman_max)])
pair_levels <- rev(pairs$pair)
for (d in list(pb_long, lr_long)) d[, `:=`(y = match(pair, pair_levels),
                                           metric = factor(metric, c('pearson', 'spearman'), c('Pearson r', 'Spearman ρ')))]
pb_long[, `:=`(set = factor(set, set_levels), y = y + (match(set, rev(set_levels)) - 3.5) * .12)]
pb <- ggplot(pb_long, aes(value, y)) +
  geom_hline(yintercept = seq(1.5, 5.5), colour = grid_col, linewidth = .2) +
  geom_linerange(data = lr_long, aes(xmin = lo, xmax = hi, y = y + .3), inherit.aes = FALSE,
                 colour = '#56B4E9', linewidth = .9) +
  geom_point(aes(colour = set, shape = set), size = 1.1, stroke = .4) +
  scale_y_continuous(breaks = 1:6, labels = pair_levels, expand = expansion(add = .45)) +
  facet_wrap(~metric, nrow = 1) +
  scale_colour_manual(values = set_col, labels = set_lab, name = NULL) +
  scale_shape_manual(values = set_shape, labels = set_lab, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .25), labels = c('0', '.25', '.5', '.75', '1')) +
  guides(colour = guide_legend(ncol = 2), shape = guide_legend(ncol = 2)) +
  labs(x = 'Agreement', y = NULL, title = 'Robustness of clock agreement',
       subtitle = 'Points: analysis set (n contrasts or target profiles)\nBlue bar: leave-one-dataset-out range, primary set') +
  theme_pub() + theme(legend.position = 'bottom', panel.grid.major.x = element_line(colour = grid_col, linewidth = .2),
                      legend.margin = margin(0, 0, 0, 0), panel.spacing = unit(3, 'mm'), legend.text = element_text(size = 5.5),
                      axis.text.y = element_text(size = 5.8), legend.key.size = unit(2.5, 'mm'))
fwrite(sets[, set := factor(set, set_levels)][order(set, pair)], file.path(out_src, 'Main_agreement_sensitivity.csv'))
fwrite(loo, file.path(out_src, 'Main_leave_one_dataset_out.csv'))
fwrite(loo_rng, file.path(out_src, 'Main_leave_one_dataset_out_range.csv'))

# ---------------------------------------------------------------- c | per-contrast effect map
wide[, mean_eff := rowMeans(.SD), .SDcols = models]
wide[, contrast_type := factor(contrast_type, names(type_col))]
setorder(wide, contrast_type, -mean_eff)
wide[, row := factor(label, rev(label))]
wide[, n_pos := rowSums(.SD > 0), .SDcols = models]
wide[, concord := fifelse(n_pos %in% c(0, 4), '4/4', fifelse(n_pos %in% c(1, 3), '3/4', '2/2'))]
pc_dat <- merge(welch, wide[, .(analysis_contrast_id, row, contrast_type)], by = 'analysis_contrast_id')
pc_dat[, `:=`(model = factor(mlab[model], mlab), txt = paste0(sprintf('%.1f', effect_weeks), fifelse(welch_p < .05, '*', '')))]
type_short <- c(diet = 'Diet', genetic_model = 'Genetic', intervention = 'Intervention', microbiome = 'M', stress_model = 'Stress')
lim_c <- 40
pc_main <- ggplot(pc_dat, aes(model, row)) +
  geom_tile(aes(fill = pmax(pmin(effect_weeks, lim_c), -lim_c)), colour = 'white', linewidth = .3) +
  geom_text(aes(label = txt, colour = abs(effect_weeks) > 25), size = 1.55) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = ink), guide = 'none') +
  scale_fill_gradient2(low = div_lo, mid = 'white', high = div_hi, limits = c(-lim_c, lim_c),
                       breaks = c(-40, -20, 0, 20, 40), labels = c('≤−40', '−20', '0', '20', '≥40'),
                       name = 'Effect (weeks)') +
  scale_x_discrete(position = 'top', labels = col_lab) +
  facet_grid(contrast_type ~ ., scales = 'free_y', space = 'free_y', switch = 'y',
             labeller = labeller(contrast_type = type_short)) +
  labs(x = NULL, y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.2),
                      axis.text.x.top = element_text(size = 5.5, colour = ink, lineheight = .85),
                      strip.placement = 'outside', strip.text.y.left = element_text(angle = 90, size = 5.5, face = 'bold'),
                      panel.spacing = unit(.8, 'mm'),
                      legend.position = 'bottom', legend.key.width = unit(6, 'mm'), legend.key.height = unit(2, 'mm'),
                      legend.title = element_text(vjust = .8))
ann <- melt(wide[, .(row, contrast_type, Mod. = fifelse(modality == 'bulk RNA-seq', 'B', 'PB'),
                     Train = fifelse(training_study, 'C', ''), `Cov.` = sprintf('%.0f', 100 * tage_coverage),
                     Sign = concord)], id.vars = c('row', 'contrast_type'), variable.name = 'col')
ann[, col := factor(col, c('Mod.', 'Train', 'Cov.', 'Sign'))]
ann[, hl := col == 'Sign' & value == '4/4']
pc_ann <- ggplot(ann, aes(col, row)) +
  geom_tile(aes(fill = hl), colour = 'white', linewidth = .3) +
  geom_text(aes(label = value), size = 1.55, colour = ink) +
  scale_fill_manual(values = c(`TRUE` = '#d9ebe0', `FALSE` = '#f4f3ef'), guide = 'none') +
  scale_x_discrete(position = 'top') +
  facet_grid(contrast_type ~ ., scales = 'free_y', space = 'free_y') +
  labs(x = NULL, y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
                      strip.text = element_blank(), panel.spacing = unit(.8, 'mm'),
                      axis.text.x.top = element_text(size = 5.5, colour = ink))
pc <- (pc_main | pc_ann) + plot_layout(widths = c(4, 3)) +
  plot_annotation(title = 'Per-contrast effects across the four clocks',
                  theme = theme(plot.title = element_text(size = 8, face = 'bold', colour = ink)))
fwrite(merge(welch, wide[, .(analysis_contrast_id, label, contrast_type, modality, training_study,
                             training_profiles, tage_coverage, concord)], by = 'analysis_contrast_id'),
       file.path(out_src, 'SFig3a_per_contrast_effects.csv'))

# ---------------------------------------------------------------- d | gene-weight concordance
gw <- tab('feature_coefficient_concordance_all_pca_source.tsv')
gs <- tab('feature_direction_statistics_all_pca.tsv')
gw[, `:=`(ovary = fifelse(grepl('^HVG', comparison), 'HVG–EN', 'Version B–PCA'),
          tage = fifelse(grepl('scaled', comparison), 'tAge scaled-diff', 'tAge YuGene-diff'))]
gs[, `:=`(ovary = fifelse(grepl('^HVG', comparison), 'HVG–EN', 'Version B–PCA'),
          tage = fifelse(grepl('scaled', comparison), 'tAge scaled-diff', 'tAge YuGene-diff'))]
gs[, lab := sprintf('n = %s\nr = %.2f%s\nρ = %.2f%s\nsign = %.0f%%', format(n, big.mark = ','), pearson,
                    stars(pearson_p), spearman, stars(spearman_p), 100 * direction_agreement)]
gene_panel <- function(ov, tg) {
  d <- gw[ovary == ov & tage == tg]; s <- gs[ovary == ov & tage == tg]
  p <- ggplot(d, aes(tage_weight, ovary_weight)) +
    geom_hline(yintercept = 0, colour = grid_col, linewidth = .3) + geom_vline(xintercept = 0, colour = grid_col, linewidth = .3)
  if (ov == 'HVG–EN') {
    d[, agree := fifelse(sign(tage_weight) == sign(ovary_weight), 'Same sign', 'Opposite sign')]
    top <- d[order(-abs(ovary_weight))][1:8]
    p <- p + geom_point(data = d, aes(colour = agree), size = 1, alpha = .9) +
      geom_text_repel(data = top, aes(label = gene), size = 1.6, fontface = 'italic', colour = ink,
                      segment.size = .2, segment.colour = ink2, min.segment.length = 0, box.padding = .2,
                      max.overlaps = Inf, seed = 1) +
      scale_colour_manual(values = c('Same sign' = '#087E8B', 'Opposite sign' = '#BB6A57'), name = 'HVG–EN gene')
  } else {
    p <- p + geom_hex(bins = 35) + scale_fill_gradient(low = '#e3e8ec', high = '#30343B', trans = 'log10', name = 'Genes')
  }
  p + geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, colour = ink, linewidth = .35) +
    annotate('text', x = -Inf, y = Inf, label = s$lab, hjust = -.05, vjust = 1.08, size = 1.75, lineheight = .9) +
    labs(x = sprintf('%s coefficient', tg), y = if (ov == 'HVG–EN') 'HVG–EN std. coefficient' else 'Version B–PCA backprojected weight',
         title = sprintf('%s vs %s', sub('Version B–PCA', 'VB–PCA', ov), sub('-diff', '', tg))) +
    theme_pub(6.5) + theme(plot.title = element_text(size = 6.5), legend.position = 'bottom')
}
pd_list <- list(gene_panel('HVG–EN', 'tAge scaled-diff'), gene_panel('HVG–EN', 'tAge YuGene-diff'),
                gene_panel('Version B–PCA', 'tAge scaled-diff'), gene_panel('Version B–PCA', 'tAge YuGene-diff'))
pd <- wrap_plots(pd_list, ncol = 2) + plot_layout(guides = 'collect') +
  plot_annotation(title = 'Gene-weight concordance (shared genes)',
                  theme = theme(plot.title = element_text(size = 8, face = 'bold', colour = ink), legend.position = 'bottom'))
fwrite(gw, file.path(out_src, 'SFig3b_gene_weights.csv'))
fwrite(gs[, !'lab'], file.path(out_src, 'SFig3b_gene_weight_agreement.csv'))

# ---------------------------------------------------------------- e | model-resolved pathway weighting
sel <- tab('model_resolved_selected_pathways.tsv')
pst <- tab('model_resolved_pathway_statistics.tsv')[pathway_id %in% sel$pathway_id]
stopifnot(nrow(sel) == 15, nrow(pst) == 15 * 4 * 2)
sel[, plab := sprintf('%d. %s', row_index, sub(' \\(GO:\\d+\\)$', '', pathway))]
sel[, plab := ifelse(nchar(plab) > 48, paste0(substr(plab, 1, 46), '…'), plab)]
grp_lab <- c('Shared pathway support' = 'Shared', 'Ovary-preferential candidates' = 'Ovary-\npreferential',
             'tAge-preferential candidates' = 'tAge-\npreferential')
pst <- merge(pst, sel[, .(pathway_id, plab, display_group, row_index)], by = 'pathway_id')
pst[, `:=`(plab = factor(plab, rev(sel[order(row_index), plab])), model = factor(mlab[model], mlab),
           analysis = factor(analysis, c('absolute', 'signed'), c('Importance (|weight| rank)', 'Direction (signed rank)')),
           display_group = factor(display_group, names(grp_lab)), star = fifelse(FDR < .05, '*', ''))]
lim_e <- max(abs(pst$rank_biserial))
pe_heat <- ggplot(pst, aes(model, plab)) +
  geom_tile(aes(fill = rank_biserial), colour = 'white', linewidth = .3) +
  geom_text(aes(label = star), size = 2.6, vjust = .75, colour = ink) +
  scale_fill_gradient2(low = div_lo, mid = 'white', high = div_hi, limits = c(-lim_e, lim_e),
                       name = 'Rank-biserial\neffect') +
  scale_x_discrete(labels = col_lab) +
  facet_grid(display_group ~ analysis, scales = 'free_y', space = 'free_y',
             labeller = labeller(display_group = grp_lab)) +
  labs(x = NULL, y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(size = 5.5),
                      axis.text.x = element_text(size = 5.3), strip.text.y = element_text(angle = 0, hjust = 0, size = 5.5),
                      panel.spacing = unit(1.2, 'mm'), legend.position = 'left', legend.key.height = unit(5, 'mm'))
sel[, `:=`(plab = factor(plab, levels(pst$plab)), display_group = factor(display_group, names(grp_lab)))]
pe_diff <- ggplot(sel, aes(difference, plab)) +
  geom_vline(xintercept = 0, linetype = 'dashed', colour = ink2, linewidth = .3) +
  geom_segment(aes(x = 0, xend = difference, yend = plab), colour = '#b9b8b3', linewidth = .5) +
  geom_point(aes(colour = difference > 0), size = 1.4) +
  scale_colour_manual(values = c(`TRUE` = div_hi, `FALSE` = div_lo), guide = 'none') +
  facet_grid(display_group ~ 'Ovary−tAge', scales = 'free_y', space = 'free_y') +
  scale_x_continuous(breaks = c(-.4, 0, .4), labels = c('−.4', '0', '.4')) +
  labs(x = 'Mean importance\ndifference', y = NULL) +
  theme_pub() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
                      strip.text.y = element_blank(), strip.text.x = element_text(size = 5.5), panel.spacing = unit(1.2, 'mm'),
                      axis.title.x = element_text(size = 5.5))
cnt <- tab('model_resolved_support_counts.tsv')
cnt_txt <- paste(sprintf('%s %s', formatC(cnt[[2]], big.mark = ',', format = 'd'), sub('^(.)', '\\L\\1', cnt[[1]], perl = TRUE)), collapse = ' · ')
pe <- (pe_heat | pe_diff) + plot_layout(widths = c(8, 1.5)) +
  plot_annotation(title = 'Pathway weighting of the four clocks (15 representative pathways)',
                  subtitle = sprintf('Of 6,637 pathways: %s; 0 supported in all four clocks', cnt_txt),
                  theme = theme(plot.title = element_text(size = 8, face = 'bold', colour = ink),
                                plot.subtitle = element_text(size = 6, colour = ink2)))
fwrite(pst, file.path(out_src, 'SFig3c_pathway_weighting.csv'))

# ---------------------------------------------------------------- assemble
# a and b (contrast scatter, agreement robustness) are used in the main figure; kept as panels only.
save_panel(pa, 'Main_contrast_scatter', 95, 85)
save_panel(pb, 'Main_agreement_robustness', 88, 85)
save_panel(pc, 'SFig3a_per_contrast_effects', 95, 120)
save_panel(pd, 'SFig3b_gene_weights', 88, 120)
save_panel(pe, 'SFig3c_pathway_weighting', 183, 75)
W <- 183; H <- 205
row1 <- (tag(pc, 'a') | tag(pd, 'b')) + plot_layout(widths = c(1.3, 1))
sfig <- (row1 / tag(pe, 'c')) + plot_layout(heights = c(1.25, .78))
ggsave(file.path(out_final, 'Supplementary_Figure3_tAge_comparison.pdf'), sfig, width = W, height = H, units = 'mm', device = cairo_pdf)
ggsave(file.path(out_final, 'Supplementary_Figure3_tAge_comparison.png'), sfig, width = W, height = H, units = 'mm', dpi = 300, bg = 'white')

# Numbers quoted in the legend.
key <- list(
  scatter = pa_stats[, .(ovary_model, tage_model, n, pearson, spearman, direction_agreement, ols_slope)],
  loo = loo_rng,
  concord = wide[, .N, by = concord],
  welch_sig = welch[, .(nominal_p05 = sum(welch_p < .05), bh_q05 = sum(welch_bh_q < .05)), by = model],
  median_abs_effect = welch[, .(median_abs_effect = median(abs(effect_weeks))), by = model],
  training_overlap = wide[, .(n_overlap = sum(training_study), n_no_overlap = sum(!training_study))],
  coverage = wide[, .(min = min(tage_coverage), median = median(tage_coverage), max = max(tage_coverage))])
sink(file.path(out_src, 'legend_numbers.txt')); print(key); sink()
print(key)
writeLines(capture.output(sessionInfo()), file.path(out_src, 'sessionInfo.txt'))
message('PASS: Supplementary Figure 3 written.')
