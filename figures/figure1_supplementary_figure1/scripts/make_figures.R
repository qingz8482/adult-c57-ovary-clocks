#!/usr/bin/env Rscript
# Figure 1 and Supplementary Figure 1 for the locked release hvg_exclusion_vbpca_20261012:
# HVG-ElasticNet trained without single-age studies (rebutal4) and Version B-PCA/glmnet
# trained on all studies (rebutal3). Reads saved pipeline outputs only; nothing is refitted.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(jsonlite)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
r3 <- file.path(dirname(r4), 'rebutal3_batcheffects removal_CV_LOO/model_development')
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
out_panels <- file.path(fig_dir, 'panels'); out_final <- file.path(fig_dir, 'final')
out_src <- file.path(fig_dir, 'source_data')
for (d in c(out_panels, out_final, out_src)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------- constants
recipes <- data.table(
  recipe = c('HVG-ElasticNet', 'VersionB-PCA'),
  dev    = c(dev4, r3),
  dir    = c('HVG-ElasticNet', 'VersionB-PCA/glmnet'),
  locked = c('HVG-ElasticNet-exclude-single-age', 'VersionB-PCA-glmnet'),
  label  = c('HVG–EN (single-age excl.)', 'Version B–PCA'))
rec_lab <- setNames(recipes$label, recipes$recipe)
rec_col <- setNames(c('#2a78d6', '#1baf7a'), recipes$label)
fold_col <- c('A·f1' = '#6250d6', 'A·f2' = '#eda100', 'A·f3' = '#008300', 'A·f4' = '#e87ba4')
fold_shape <- c('A·f1' = 16, 'A·f2' = 17, 'A·f3' = 15, 'A·f4' = 18)
sign_col <- c('Older (+)' = '#b2182b', 'Younger (−)' = '#2166ac')
ink <- '#0b0b0b'; ink2 <- '#52514e'; grid_col <- '#e6e5e0'
stage_levels <- c('A·f1', 'A·f2', 'A·f3', 'A·f4', 'B', 'C', 'Final')
stage_dirs <- c('A_r00_f1', 'A_r00_f2', 'A_r00_f3', 'A_r00_f4', 'B_r00_f1', 'C_r00_f1', 'final_r00_f0')
single_age <- c('2022_eLife_Morris_cycling_murine_ovary_atlas', 'CRA017492', 'GSE252791',
                '2025_CommBio_Zhao_developing_ovary_scRNA_spatial', 'GSE194196')
short <- c(
  '2022_eLife_Morris_cycling_murine_ovary_atlas' = 'Morris',
  '2024_Cell_Winkler_FRT_scRNA_spatial' = 'Winkler',
  '2024_LifeMedicine_Jiang_whole_life_mouse_ovary' = 'Jiang',
  '2024_NatureAging_Isola_aging_mouse_ovary_atlas' = 'Isola',
  '2025_CommBio_Zhao_developing_ovary_scRNA_spatial' = 'Zhao',
  '2025_GeroScience_Zhang_multiomic_aging_mouse_ovary' = 'Zhang',
  '2026_CommBio_Galligos_ovarian_inflammaging' = 'Galligos',
  'NIHMS1554107' = 'NIHMS')
short_name <- function(s) ifelse(s %in% names(short), short[s], s)

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
save_panel <- function(p, name, w, h) {
  ggsave(file.path(out_panels, paste0(name, '.pdf')), p, width = w, height = h, units = 'mm', device = cairo_pdf)
}
tag <- function(p, t) wrap_elements(full = p) + labs(tag = t) +
  theme(plot.tag = element_text(size = 9, face = 'bold'), plot.tag.position = c(0, 1), plot.margin = margin(2, 1, 1, 9))

# ---------------------------------------------------------------- inputs
rpath <- function(r, ...) file.path(recipes[recipe == r, dev], recipes[recipe == r, dir], ...)
rd <- function(r, f, mode = 'count_tuning') fread(rpath(r, mode, 'results', f))
fa <- rbindlist(lapply(recipes$recipe, function(r) rd(r, 'fit_audit.csv')[analysis %in% c('A', 'B', 'C', 'final') & repeat_id == 0][, recipe := r]))
fa[, stage := fifelse(analysis == 'A', paste0('A·f', fold), fifelse(analysis == 'final', 'Final', analysis))]
pred <- rbindlist(lapply(recipes$recipe, function(r) rbindlist(lapply(c('A', 'B', 'C'), function(s) rd(r, paste0('predictions_', s, '.csv'))))))
stopifnot(pred[analysis == 'A', .N, by = recipe]$N == 136, uniqueN(pred[analysis == 'A']$sample_key) == 136)
met <- rbindlist(lapply(recipes$recipe, function(r) fread(file.path(recipes[recipe == r, dev], 'comparison/all_metrics.csv'))[recipe == r]))
meta <- fread(file.path(dev4, 'shared_inputs/sample_metadata.csv'), select = c('sample_key', 'study', 'modality', 'age_week'))
stopifnot(nrow(meta) == 136)
manifest <- fromJSON(file.path(lock, 'manifest.json'))
hp <- rbindlist(lapply(recipes$recipe, function(r) as.data.table(fromJSON(file.path(lock, recipes[recipe == r, locked], 'hyperparameters.json')))[, recipe := r]))
# Locked artifacts must be the exact final fits whose development outputs are plotted.
for (r in recipes$recipe) stopifnot(identical(readRDS(file.path(lock, recipes[recipe == r, locked], 'model.rds')),
                                              readRDS(rpath(r, 'count_tuning/results/final_model/model.rds'))))
stopifnot(fa[recipe == 'HVG-ElasticNet' & stage == 'Final', n_train] == 110, fa[recipe == 'VersionB-PCA' & stage == 'Final', n_train] == 136)
stopifnot(!any(single_age %in% unlist(strsplit(fa[recipe == 'HVG-ElasticNet', train_studies], ';'))))

# ================================================================= FIGURE 1
## a | nested tuning schematic -------------------------------------------
boxes <- data.table(
  id = 1:6, x = 1:6,
  head = c('1  Outer split', '2  Training rule', '3  Fold-local genes', '4  Regression', '5  Select (inner LOSO)', '6  Refit + predict'),
  body = c('whole studies held out\nA: 4 folds (36/48/20/32)\nB: NIHMS  ·  C: Winkler\nFinal: no held-out set',
           'HVG–EN: drop single-age\nstudies from training\n(5 studies, 26 profiles)\nVersion B–PCA: all studies',
           'ln(1+CPM); CPM≥1 in ≥50%\nN ∈ 500/1k/2k/3k/5k\nHVG: variance rank\nVersion B: age partial R²',
           'scale → [PCA 1–10 PCs]\n→ study-weighted\nelastic net, 6 α ×\n100 fixed λ',
           'each training study held\nout once; mean study MAE;\nλ: largest within 1 SE;\nthen best (N, α, PCs)',
           'refit on outer-training\nset; frozen prediction\nof held-out studies;\nno calibration'),
  fill = c('#eef1f5', '#fff6e0', '#ffffff', '#ffffff', '#ffffff', '#eef1f5'))
bw <- .44
pb <- ggplot(boxes) +
  annotate('rect', xmin = 2 - bw - .08, xmax = 5 + bw + .08, ymin = -.62, ymax = .78, fill = NA, colour = ink2, linetype = 'dashed', linewidth = .3) +
  annotate('text', x = 3.5, y = .74, label = 'applied inside every outer training set', vjust = 1, size = 1.7, colour = ink2, fontface = 'italic') +
  geom_rect(aes(xmin = x - bw, xmax = x + bw, ymin = -.5, ymax = .5, fill = fill), colour = ink2, linewidth = .3) +
  geom_text(aes(x - bw + .04, .42, label = head), hjust = 0, vjust = 1, size = 2, fontface = 'bold', colour = ink) +
  geom_text(aes(x - bw + .04, .2, label = body), hjust = 0, vjust = 1, size = 1.65, colour = ink2, lineheight = .95) +
  geom_segment(data = boxes[id < 6], aes(x = x + bw, xend = x + 1 - bw, y = 0, yend = 0),
               arrow = arrow(length = unit(.9, 'mm'), type = 'closed'), linewidth = .3, colour = ink2) +
  scale_fill_identity() + scale_x_continuous(limits = c(.5, 6.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-.66, .82), expand = c(0, 0)) + theme_void()

## b | cohort and validation design --------------------------------------
test_map <- fa[recipe == 'VersionB-PCA' & stage != 'Final', .(study = unlist(strsplit(test_studies, ';'))), by = stage]
st <- meta[, .(n = .N, modality = modality[1], n_ages = uniqueN(age_week)), by = study]
st[, fold := test_map[grepl('^A', stage)][match(st$study, study), stage]]
setorder(st, fold, -n)
st[, lab := sprintf('%s%s (n=%d)', short_name(study), ifelse(study %in% single_age, '†', ''), n)]
st[, lab := factor(lab, levels = rev(lab))]
design <- CJ(study = st$study, stage = stage_levels)
design[, role := fifelse(study %in% single_age, 'Train: Version B–PCA only', 'Train: both clocks')]
design[test_map, on = .(study, stage), role := 'Test (held out)']
design[, `:=`(lab = st$lab[match(study, st$study)], stage = factor(stage, levels = stage_levels),
              role = factor(role, levels = c('Test (held out)', 'Train: both clocks', 'Train: Version B–PCA only')))]
# Training membership in the plotted design must match the fit audits of both clocks.
for (r in recipes$recipe) for (s in stage_levels) {
  used <- unlist(strsplit(fa[recipe == r & stage == s, train_studies], ';'))
  shown <- design[stage == s & (role == 'Train: both clocks' | (r == 'VersionB-PCA' & role == 'Train: Version B–PCA only')), study]
  stopifnot(setequal(used, shown))
}
fwrite(design, file.path(out_src, 'Fig1b_design.csv'))
ages <- merge(meta, st[, .(study, lab)], by = 'study')
ages[, Modality := fifelse(modality == 'bulk_rnaseq', 'Bulk RNA-seq', 'scRNA-seq pseudobulk')]
pa1 <- ggplot(ages, aes(age_week, lab)) +
  geom_point(aes(shape = Modality), size = 1.1, colour = ink2, alpha = .8, position = position_jitter(height = .15, width = 0, seed = 1)) +
  scale_shape_manual(values = c('Bulk RNA-seq' = 16, 'scRNA-seq pseudobulk' = 1)) +
  scale_x_continuous('Age (weeks)', breaks = c(0, 26, 52, 78)) + labs(y = NULL) +
  theme_pub() + theme(legend.position = 'bottom', legend.direction = 'vertical', legend.title = element_blank(),
                      panel.grid.major.y = element_line(colour = grid_col, linewidth = .25))
pa2 <- ggplot(design, aes(stage, lab)) +
  geom_tile(aes(fill = role), colour = 'white', linewidth = .6) +
  scale_fill_manual(NULL, values = c('Test (held out)' = '#34373d', 'Train: both clocks' = '#d7dde5', 'Train: Version B–PCA only' = '#f2d79b')) +
  scale_x_discrete(position = 'top', labels = c('f1', 'f2', 'f3', 'f4', 'B', 'C', 'Final')) +
  labs(x = 'Stage A folds  ·  B  ·  C  ·  final fit', y = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  theme_pub() + theme(axis.text.y = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
                      legend.position = 'bottom', axis.title.x.top = element_text(size = 5.5, colour = ink2))
pa <- (pa1 | pa2) + plot_layout(widths = c(1, 1.15))

## c | per-study held-out MAE --------------------------------------------
ps <- met[mode == 'count_tuning' & repeat_id == 0 & ((analysis == 'A' & scope == 'study') | (analysis %in% c('B', 'C') & scope == 'pooled'))]
ps[analysis %in% c('B', 'C'), label := c(B = 'NIHMS1554107', C = '2024_Cell_Winkler_FRT_scRNA_spatial')[analysis]]
ps[, Recipe := factor(rec_lab[recipe], levels = recipes$label)]
ps[, Block := factor(fifelse(analysis == 'A', 'A', 'B, C'), levels = c('A', 'B, C'))]
ps[, fold := st$fold[match(label, st$study)]]
ps[, ylab := fifelse(analysis == 'A', sprintf('%s%s · %s', short_name(label), ifelse(label %in% single_age, '†', ''), sub('A·', '', fold)),
                     sprintf('%s · %s', short_name(label), analysis))]
ord <- ps[analysis == 'A', .(m = mean(MAE)), by = ylab][order(m), ylab]
ps[, ylab := factor(ylab, levels = c(rev(sort(unique(ps[analysis != 'A', ylab]))), ord))]
fwrite(ps[, .(Recipe, analysis, study = label, fold, n, MAE, bias)], file.path(out_src, 'Fig1c_per_study_MAE.csv'))
pe <- ggplot(ps, aes(MAE, ylab)) +
  geom_line(aes(group = ylab), colour = grid_col, linewidth = .8) +
  geom_point(aes(colour = Recipe, shape = Recipe), size = 1.5) +
  scale_colour_manual(NULL, values = rec_col) + scale_shape_manual(NULL, values = c(16, 15)) +
  facet_grid(Block ~ ., scales = 'free_y', space = 'free_y') +
  labs(x = 'Held-out MAE (weeks)', y = NULL) +
  theme_pub() + theme(legend.position = 'bottom', legend.direction = 'horizontal',
                      panel.grid.major.y = element_line(colour = grid_col, linewidth = .2), strip.text.y = element_text(size = 6),
                      plot.margin = margin(8, 4, 3, 8))

## d | inner-CV count selection per stage ---------------------------------
cand <- rbindlist(lapply(recipes$recipe, function(r) rbindlist(lapply(seq_along(stage_dirs), function(i) {
  fread(rpath(r, 'tuning_audit', stage_dirs[i], 'count_tuning_candidate_summary.csv'))[, `:=`(recipe = r, stage = stage_levels[i])]
}))))
best_n <- cand[, .SD[which.min(inner_MAE)], by = .(recipe, stage, N)]
best_n[, delta := inner_MAE - min(inner_MAE), by = .(recipe, stage)]
sel <- fa[, .(recipe, stage, N = requested_N)]
best_n[, selected := FALSE][sel, on = .(recipe, stage, N), selected := TRUE]
stopifnot(best_n[selected == TRUE, .N] == 14, all(best_n[selected == TRUE, delta] < 1e-9))
best_n[, `:=`(Recipe = factor(rec_lab[recipe], levels = recipes$label), stage = factor(stage, levels = rev(stage_levels)),
              Nf = factor(N, levels = c(500, 1000, 2000, 3000, 5000), labels = c('500', '1k', '2k', '3k', '5k')))]
fwrite(best_n, file.path(out_src, 'Fig1d_inner_cv_by_count.csv'))
pc <- ggplot(best_n, aes(Nf, stage)) +
  geom_tile(aes(fill = pmin(delta, 4)), colour = 'white', linewidth = .5) +
  geom_tile(data = best_n[selected == TRUE], fill = NA, colour = ink, linewidth = .55) +
  geom_text(aes(label = sprintf('%.1f', inner_MAE), colour = delta < 1.6), size = 1.7) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = ink), guide = 'none') +
  scale_fill_gradient('Δ inner MAE vs best\n(weeks, capped at 4)', low = '#173f73', high = '#eef3fa', limits = c(0, 4)) +
  facet_wrap(~Recipe, nrow = 1) + labs(x = 'Input gene count N', y = NULL) +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = 'right',
                      legend.key.height = unit(5, 'mm'))

## e | held-out predictions, stages A/B/C ---------------------------------
pr <- copy(pred)
pr[, Recipe := factor(rec_lab[recipe], levels = recipes$label)]
pr[, Stage := factor(c(A = 'A · 4 folds', B = 'B · NIHMS', C = 'C · Winkler')[analysis], levels = c('A · 4 folds', 'B · NIHMS', 'C · Winkler'))]
pr[, Fold := fifelse(analysis == 'A', paste0('A·f', fold), 'single-study holdout')]
pr[, `:=`(err = pred - age_week, single = study %in% single_age)]
fwrite(pr[, .(sample_key, study, modality, age_week, pred, analysis, fold, recipe, selected_N, selected_pc)], file.path(out_src, 'Fig1e_heldout_predictions.csv'))
ann <- pr[, .(MAE = mean(abs(err)), r = cor(pred, age_week), n = .N), by = .(Recipe, Stage, analysis, recipe)]
ann <- merge(ann, met[scope == 'macro_study' & mode == 'count_tuning' & repeat_id == 0, .(recipe, analysis, macro = MAE)], by = c('recipe', 'analysis'))
ann[, txt := ifelse(analysis == 'A', sprintf('MAE %.1f (study %.1f)\nr = %.2f', MAE, macro, r), sprintf('MAE %.1f\nr = %.2f', MAE, r))]
lim <- c(min(0, floor(min(pr$pred))) - 1, max(c(pr$pred, pr$age_week)) * 1.03)  # uncalibrated predictions can be < 0
pd <- ggplot(pr, aes(age_week, pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = .3, colour = ink2) +
  geom_point(data = pr[analysis == 'A'], aes(colour = Fold, shape = Fold), size = 1.1, alpha = .85) +
  geom_point(data = pr[analysis != 'A'], colour = '#34373d', shape = 1, size = 1.1, stroke = .4) +
  geom_text(data = ann, aes(x = lim[1] + 1, y = lim[2], label = txt), hjust = 0, vjust = 1, size = 1.8, colour = ink, lineheight = .9) +
  scale_colour_manual(NULL, values = fold_col) + scale_shape_manual(NULL, values = fold_shape) +
  scale_x_continuous('Chronological age (weeks)', limits = lim, breaks = c(0, 26, 52, 78)) +
  scale_y_continuous('Predicted age (weeks)', limits = lim, breaks = c(0, 26, 52, 78)) +
  coord_equal() + facet_grid(Stage ~ Recipe) +
  theme_pub() + theme(legend.position = 'top', legend.margin = margin(0, 0, -4, 0),
                      panel.border = element_rect(fill = NA, colour = grid_col, linewidth = .3), strip.text = element_text(size = 6.2))

## f | locked final models -----------------------------------------------
coefs <- fread(file.path(lock, 'HVG-ElasticNet-exclude-single-age/coefficients.csv'))
en <- coefs[coefficient != 0]
stopifnot(nrow(en) == manifest$models$nonzero_effective_gene_coefficients[1])
en[, Direction := fifelse(standardized_coefficient > 0, 'Older (+)', 'Younger (−)')]
fwrite(en, file.path(out_src, 'Fig1f_final_nonzero_gene_coefficients.csv'))
hp_txt <- function(r) { h <- hp[recipe == r]; sprintf('%s  ·  N = %s, α = %g, λ = %.3g%s', rec_lab[r], format(h$N, big.mark = ','), h$alpha, h$lambda,
                                                     ifelse(h$pc > 0, sprintf(', %d PCs', h$pc), '')) }
top <- en[order(-abs(standardized_coefficient))][1:20]
top[, gene := factor(gene, levels = gene[order(standardized_coefficient)])]
pg1 <- ggplot(top, aes(standardized_coefficient, gene)) +
  geom_vline(xintercept = 0, linewidth = .3, colour = ink2) +
  geom_segment(aes(x = 0, xend = standardized_coefficient, yend = gene, colour = Direction), linewidth = .5) +
  geom_point(aes(colour = Direction), size = 1.1) +
  scale_colour_manual(NULL, values = sign_col) +
  labs(x = 'Coefficient (weeks per SD)', y = NULL, title = hp_txt('HVG-ElasticNet'),
       subtitle = sprintf('top 20 of %d nonzero genes; trained on 110 profiles / 9 studies', nrow(en))) +
  theme_pub(6.5) + theme(axis.text.y = element_text(size = 5, face = 'italic'), plot.title = element_text(size = 6.2),
                         plot.subtitle = element_text(size = 5, colour = ink2), legend.position = 'bottom')
pcc <- fread(file.path(lock, 'VersionB-PCA-glmnet/pc_coefficients.csv'))
pcc[, `:=`(PCf = factor(paste0('PC', PC), levels = paste0('PC', 1:10)),
           Direction = fifelse(coefficient > 0, 'Older (+)', fifelse(coefficient < 0, 'Younger (−)', 'Zero')))]
fwrite(pcc, file.path(out_src, 'Fig1f_final_pc_coefficients.csv'))
pg2 <- ggplot(pcc, aes(PCf, coefficient)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = ink2) +
  geom_col(aes(fill = Direction), width = .65) +
  geom_point(data = pcc[coefficient == 0], shape = 4, size = .8, colour = ink2) +
  scale_fill_manual(values = c(sign_col, Zero = 'grey80'), guide = 'none') +
  labs(x = NULL, y = 'Weeks per PC unit', title = hp_txt('VersionB-PCA'),
       subtitle = '136 profiles / 14 studies; 4 of 10 PCs nonzero (× = 0)\ngene loadings in Supplementary Fig. 2') +
  theme_pub(6.5) + theme(plot.title = element_text(size = 6.2), plot.subtitle = element_text(size = 5, colour = ink2))
pg <- pg1 / free(pg2) + plot_layout(heights = c(1.9, .8))

## assemble Figure 1 -------------------------------------------------------
fig1 <- tag(pb, 'a') / (tag(pa, 'b') | tag(pe, 'c')) / tag(pc, 'd') / (tag(pd, 'e') | tag(pg, 'f')) +
  plot_layout(heights = c(.26, .8, .42, 1.25))
for (x in list(list(pb, 'Fig1a_schematic', 183, 30), list(pa, 'Fig1b_design', 92, 64), list(pe, 'Fig1c_per_study', 91, 64),
               list(pc, 'Fig1d_inner_count', 183, 36), list(pd, 'Fig1e_predictions', 95, 105), list(pg, 'Fig1f_final_models', 88, 105))) save_panel(x[[1]], x[[2]], x[[3]], x[[4]])
W1 <- 183; H1 <- 240
ggsave(file.path(out_final, 'Figure1_model_development.pdf'), fig1, width = W1, height = H1, units = 'mm', device = cairo_pdf)
ggsave(file.path(out_final, 'Figure1_model_development.png'), fig1, width = W1, height = H1, units = 'mm', dpi = 300)

# ======================================================== SUPPLEMENTARY FIG 1
## a | final-fit inner-CV lambda path ------------------------------------
lam_final <- rbindlist(lapply(recipes$recipe, function(r) {
  h <- fa[recipe == r & stage == 'Final']
  x <- fread(rpath(r, 'tuning_audit/final_r00_f0/count_tuning_lambda_summary.csv'))
  x[N == h$requested_N & pc == h$PCs & abs(alpha - h$alpha) < 1e-9][, `:=`(recipe = r, sel_lambda = h$lambda)]
}))
lam_final[, Recipe := factor(rec_lab[recipe], levels = recipes$label)]
nz <- c('HVG-ElasticNet' = nrow(en), 'VersionB-PCA' = sum(pcc$coefficient != 0))
lam_mark <- lam_final[, .(sel = sel_lambda[1], minl = lambda[which.min(mean_mae)], thr = min(mean_mae) + se[which.min(mean_mae)]), by = .(recipe, Recipe)]
lam_mark <- merge(lam_mark, hp[, .(recipe, N, pc, alpha, threshold)], by = 'recipe')
stopifnot(all(abs(lam_mark$thr - lam_mark$threshold) < 1e-3))
lam_mark[, txt := sprintf('N = %d%s · α = %g\nλ = %.3g · %d %s · %d inner studies', N, ifelse(pc > 0, sprintf(' · %d PCs', pc), ''), alpha, sel,
                          nz[recipe], ifelse(recipe == 'VersionB-PCA', 'active PCs', 'nonzero genes'), ifelse(recipe == 'VersionB-PCA', 14L, 9L))]
fwrite(lam_final, file.path(out_src, 'SFig1a_final_lambda_path.csv'))
sa0 <- ggplot(lam_final, aes(log10(lambda), mean_mae)) +
  geom_ribbon(aes(ymin = mean_mae - se, ymax = mean_mae + se, fill = Recipe), alpha = .18) +
  geom_line(aes(colour = Recipe), linewidth = .6) +
  geom_hline(data = lam_mark, aes(yintercept = thr), linetype = 'dotted', linewidth = .3, colour = ink2) +
  geom_vline(data = lam_mark, aes(xintercept = log10(minl)), linetype = 'dashed', linewidth = .3, colour = ink2) +
  geom_vline(data = lam_mark, aes(xintercept = log10(sel)), linewidth = .45, colour = '#b2182b') +
  geom_label(data = lam_mark, aes(x = -4, y = Inf, label = txt), hjust = 0, vjust = 1.1, size = 1.75, colour = ink, lineheight = .9,
             linewidth = 0, fill = 'white', alpha = .85) +
  scale_colour_manual(values = rec_col, guide = 'none') + scale_fill_manual(values = rec_col, guide = 'none') +
  scale_x_continuous(expression(log[10]~lambda), breaks = seq(-4, 3, 1)) +
  facet_wrap(~Recipe, nrow = 1) + labs(y = 'Inner MAE\n(weeks ± SE)') + theme_pub()

## b | final-fit inner-CV landscape --------------------------------------
fc <- cand[stage == 'Final']
land <- fc[, .SD[which.min(inner_MAE)], by = .(recipe, N, alpha)]
land[, `:=`(Nf = factor(N), af = factor(alpha))]
land[, sel := FALSE][hp, on = .(recipe, N, alpha), sel := TRUE]
pcl <- fc[recipe == 'VersionB-PCA', .SD[which.min(inner_MAE)], by = .(N, pc)]
pcl[, `:=`(sel = N == hp[recipe == 'VersionB-PCA', N] & pc == hp[recipe == 'VersionB-PCA', pc], Nf = factor(N), pcf = factor(pc))]
fwrite(rbind(land[, .(recipe, N, alpha, pc, lambda, inner_MAE, inner_SE)], pcl[, .(recipe = 'VersionB-PCA (N x PC)', N, alpha, pc, lambda, inner_MAE, inner_SE)]), file.path(out_src, 'SFig1b_final_landscape.csv'))
heat <- function(d, x, y, xl, yl, ttl, lims) ggplot(d, aes({{x}}, {{y}})) +
  geom_tile(aes(fill = inner_MAE), colour = 'white', linewidth = .4) +
  geom_tile(data = d[sel == TRUE], fill = NA, colour = ink, linewidth = .55) +
  geom_text(aes(label = sprintf('%.1f', inner_MAE), colour = inner_MAE < mean(lims)), size = 1.5) +
  scale_colour_manual(values = c(`TRUE` = 'white', `FALSE` = ink), guide = 'none') +
  scale_fill_gradient('Inner MAE\n(weeks)', low = '#173f73', high = '#eef3fa', limits = lims, oob = scales::squish) +
  labs(x = xl, y = yl, title = ttl) + theme_pub() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), plot.title = element_text(size = 6.5),
        legend.position = 'bottom', legend.key.width = unit(6, 'mm'), legend.key.height = unit(2, 'mm'), legend.title = element_text(size = 5.5, vjust = .8))
sa <- heat(land[recipe == 'HVG-ElasticNet'], Nf, af, 'N', 'α', paste(rec_lab[1], '(9 inner studies)'), c(6.5, 10)) +
  heat(land[recipe == 'VersionB-PCA'], Nf, af, 'N', 'α', paste(rec_lab[2], '(best PC; 14 studies)'), c(13, 17.5)) +
  heat(pcl, Nf, pcf, 'N', 'PCs', paste(rec_lab[2], '(best α)'), c(13, 17.5)) +
  plot_layout(nrow = 1)

## c | lambda paths of the selected candidate, each outer stage ----------
lam_st <- rbindlist(lapply(recipes$recipe, function(r) rbindlist(lapply(1:6, function(i) {
  h <- fa[recipe == r & stage == stage_levels[i]]
  x <- fread(rpath(r, 'tuning_audit', stage_dirs[i], 'count_tuning_lambda_summary.csv'))
  x[N == h$requested_N & pc == h$PCs & abs(alpha - h$alpha) < 1e-9][, `:=`(recipe = r, stage = stage_levels[i], sel = h$lambda)]
}))))
lam_st[, `:=`(Recipe = factor(rec_lab[recipe], levels = recipes$label), stage = factor(stage, levels = stage_levels))]
lsm <- lam_st[, .(sel = sel[1], thr = min(mean_mae) + se[which.min(mean_mae)]), by = .(Recipe, stage)]
fwrite(lam_st, file.path(out_src, 'SFig1c_stage_lambda_paths.csv'))
sb <- ggplot(lam_st, aes(log10(lambda), mean_mae)) +
  geom_ribbon(aes(ymin = mean_mae - se, ymax = mean_mae + se, fill = Recipe), alpha = .18) +
  geom_line(aes(colour = Recipe), linewidth = .5) +
  geom_hline(data = lsm, aes(yintercept = thr), linetype = 'dotted', linewidth = .3, colour = ink2) +
  geom_vline(data = lsm, aes(xintercept = log10(sel)), linewidth = .4, colour = '#b2182b') +
  scale_colour_manual(values = rec_col, guide = 'none') + scale_fill_manual(values = rec_col, guide = 'none') +
  scale_x_continuous(expression(log[10]~lambda), breaks = c(-4, -2, 0, 2)) +
  facet_grid(Recipe ~ stage, scales = 'free_y', labeller = labeller(Recipe = function(x) sub(' \\(', '\n(', x))) +
  labs(y = 'Inner LOSO MAE (weeks)') +
  theme_pub(6.5) + theme(panel.border = element_rect(fill = NA, colour = grid_col, linewidth = .3), strip.text.y = element_text(size = 5.5))

## d | selected hyperparameters table ------------------------------------
tb <- fa[, .(Recipe = factor(rec_lab[recipe], levels = recipes$label), stage = factor(stage, levels = rev(stage_levels)),
             `Train/\ntest` = ifelse(n_test > 0, sprintf('%d/%d', n_train, n_test), sprintf('%d/–', n_train)),
             `Train\nstudies` = as.character(lengths(strsplit(train_studies, ';'))),
             N = format(requested_N, big.mark = ','), `α` = format(alpha), `λ` = as.character(signif(lambda, 3)),
             PCs = ifelse(PCs > 0, as.character(PCs), '–'),
             `Nonzero\ngenes` = ifelse(recipe == 'VersionB-PCA', 'dense', as.character(n_nonzero_effective_genes)),
             `Inner\nMAE` = sprintf('%.2f', inner_macro_MAE))]
fwrite(tb, file.path(out_src, 'SFig1d_selected_hyperparameters.csv'))
tl <- melt(tb, id.vars = c('Recipe', 'stage'), variable.name = 'col', value.name = 'v', variable.factor = TRUE)
sc <- ggplot(tl, aes(col, stage)) +
  geom_tile(data = tl[stage == 'Final'], fill = '#eef1f5', colour = NA) +
  geom_text(aes(label = v), size = 1.75, colour = ink) +
  scale_x_discrete(position = 'top') + facet_wrap(~Recipe, nrow = 1) + labs(x = NULL, y = NULL) +
  theme_pub(6.5) + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(size = 5.2, face = 'bold', colour = ink))

## e | HVG-EN training rule: all studies vs single-age exclusion ---------
hv <- rbind(fread(file.path(r3, 'comparison/all_metrics.csv'))[recipe == 'HVG-ElasticNet' & mode == 'count_tuning'][, Rule := 'All studies (previous lock)'],
            met[recipe == 'HVG-ElasticNet' & mode == 'count_tuning'][, Rule := 'Single-age excluded (locked)'])
hv[, Rule := factor(Rule, levels = c('All studies (previous lock)', 'Single-age excluded (locked)'))]
h1 <- hv[repeat_id == 0 & ((analysis == 'A' & scope %in% c('pooled', 'macro_study')) | (analysis %in% c('B', 'C') & scope == 'pooled'))]
h1[, Metric := factor(fifelse(analysis == 'A', fifelse(scope == 'pooled', 'A pooled', 'A study-mean'), paste(analysis, c(B = 'NIHMS', C = 'Winkler')[analysis])),
                      levels = rev(c('A pooled', 'A study-mean', 'B NIHMS', 'C Winkler')))]
h2 <- hv[analysis == 'A_repeated' & scope == 'pooled']
stopifnot(h2[, .N, by = Rule]$N == 30)
fwrite(rbind(h1[, .(Rule, Metric = as.character(Metric), repeat_id, MAE)], h2[, .(Rule, Metric = 'A_repeated pooled', repeat_id, MAE)]), file.path(out_src, 'SFig1e_hvg_training_rule.csv'))
rule_shape <- c('All studies (previous lock)' = 1, 'Single-age excluded (locked)' = 16)
se1 <- ggplot(h1, aes(MAE, Metric)) +
  geom_line(aes(group = Metric), colour = grid_col, linewidth = .8) +
  geom_point(aes(shape = Rule), colour = rec_col[1], size = 1.6, stroke = .5) +
  scale_shape_manual(NULL, values = rule_shape) + labs(x = 'Held-out MAE (weeks)', y = NULL, title = 'HVG–EN: main split and single-study holdouts') +
  theme_pub(6.5) + theme(legend.position = 'bottom', plot.title = element_text(size = 6.5),
                         panel.grid.major.y = element_line(colour = grid_col, linewidth = .2))
se2 <- ggplot(h2, aes(Rule, MAE)) +
  geom_line(aes(group = repeat_id), colour = grid_col, linewidth = .3) +
  geom_boxplot(width = .4, outlier.shape = NA, fill = NA, colour = rec_col[1], linewidth = .35) +
  geom_point(aes(shape = Rule), colour = rec_col[1], size = .7, position = position_jitter(width = .05, seed = 3)) +
  scale_shape_manual(values = rule_shape, guide = 'none') +
  scale_x_discrete(labels = c('All studies', 'Single-age\nexcluded')) +
  labs(x = NULL, y = 'Pooled MAE (weeks)', title = sprintf('30 repeated splits (excluded better in %d/30)',
       sum(h2[Rule == levels(Rule)[2]][order(repeat_id), MAE] < h2[Rule == levels(Rule)[1]][order(repeat_id), MAE]))) +
  theme_pub(6.5) + theme(plot.title = element_text(size = 6.5))
se <- (se1 | se2) + plot_layout(widths = c(1.4, 1))

## f | count tuning vs fixed 2,000 ---------------------------------------
cmp <- met[repeat_id == 0 & ((analysis == 'A' & scope %in% c('pooled', 'macro_study')) | (analysis %in% c('B', 'C') & scope == 'pooled'))]
cmp[, Metric := factor(fifelse(analysis == 'A', fifelse(scope == 'pooled', 'A pooled', 'A study-mean'), paste(analysis, c(B = 'NIHMS', C = 'Winkler')[analysis])),
                       levels = rev(c('A pooled', 'A study-mean', 'B NIHMS', 'C Winkler')))]
cmp[, `:=`(Recipe = factor(rec_lab[recipe], levels = recipes$label), Procedure = factor(fifelse(mode == 'fixed_2000', 'Fixed N = 2,000', 'Count-tuned (locked)'), levels = c('Fixed N = 2,000', 'Count-tuned (locked)')))]
fwrite(cmp[, .(Recipe, Metric, Procedure, MAE)], file.path(out_src, 'SFig1f_fixed_vs_tuned.csv'))
sd <- ggplot(cmp, aes(MAE, Metric)) +
  geom_line(aes(group = Metric), colour = grid_col, linewidth = .8) +
  geom_point(aes(shape = Procedure, colour = Recipe), size = 1.5, stroke = .5) +
  scale_shape_manual(NULL, values = c(1, 16)) + scale_colour_manual(values = rec_col, guide = 'none') +
  facet_wrap(~Recipe, nrow = 1, labeller = labeller(Recipe = function(x) sub(' \\(', '\n(', x))) + labs(x = 'Held-out MAE (weeks)', y = NULL) +
  theme_pub(6.5) + theme(legend.position = 'bottom', panel.grid.major.y = element_line(colour = grid_col, linewidth = .2))

## g | 30 repeated study partitions, fixed vs tuned ----------------------
rep <- met[analysis == 'A_repeated' & scope %in% c('pooled', 'macro_study')]
stopifnot(rep[, .N, by = .(recipe, mode, scope)]$N == 30)
rep[, `:=`(Recipe = factor(rec_lab[recipe], levels = recipes$label),
           Procedure = factor(fifelse(mode == 'fixed_2000', 'Fixed 2k', 'Tuned'), levels = c('Fixed 2k', 'Tuned')),
           Metric = factor(fifelse(scope == 'pooled', 'Pooled', 'Study-mean'), levels = c('Pooled', 'Study-mean')))]
fwrite(rep[, .(Recipe, Procedure, Metric, repeat_id, MAE)], file.path(out_src, 'SFig1g_repeated_partitions.csv'))
sr <- ggplot(rep, aes(Procedure, MAE)) +
  geom_line(aes(group = repeat_id), colour = grid_col, linewidth = .25) +
  geom_boxplot(aes(colour = Recipe), width = .45, outlier.shape = NA, fill = NA, linewidth = .35) +
  geom_point(aes(colour = Recipe), size = .5, alpha = .6, position = position_jitter(width = .06, seed = 2)) +
  scale_colour_manual(values = rec_col, guide = 'none') +
  facet_grid(Metric ~ Recipe, scales = 'free_y', labeller = labeller(Recipe = function(x) sub(' \\(', '\n(', x))) +
  labs(x = NULL, y = 'MAE (weeks), 30 repeats') +
  theme_pub(6.5) + theme(panel.border = element_rect(fill = NA, colour = grid_col, linewidth = .3), strip.text.y = element_text(size = 5.5))

## h | selected count frequency across repeated outer fits ----------------
fr <- rbindlist(lapply(recipes$recipe, function(r) rd(r, 'selected_count_frequencies.csv')[analysis == 'A_repeated'][, recipe := r]))
fr <- fr[, .(n_models = sum(n_models)), by = .(recipe, requested_N)]
fr <- merge(CJ(recipe = recipes$recipe, requested_N = c(500, 1000, 2000, 3000, 5000)), fr, all.x = TRUE)[is.na(n_models), n_models := 0]
stopifnot(fr[, sum(n_models), by = recipe]$V1 == 120)
fr[, `:=`(Recipe = factor(rec_lab[recipe], levels = recipes$label), Nf = factor(requested_N, labels = c('500', '1k', '2k', '3k', '5k')))]
fr[, final := FALSE][hp, on = .(recipe, requested_N = N), final := TRUE]
fwrite(fr, file.path(out_src, 'SFig1h_selected_N_frequency.csv'))
sf <- ggplot(fr, aes(Nf, n_models, fill = Recipe)) +
  geom_col(width = .7) +
  geom_text(data = fr[final == TRUE], aes(label = 'final'), vjust = -.4, size = 1.6, colour = ink) +
  scale_fill_manual(values = rec_col, guide = 'none') + scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  facet_wrap(~Recipe, nrow = 1) + labs(x = 'Selected N (inner CV)', y = 'Outer fits\n(of 120)') + theme_pub(6.5)

## i | stability of final HVG genes across repeated outer fits ------------
sfreq <- rd('HVG-ElasticNet', 'feature_stability_A_repeated.csv')
stab <- copy(en)
stab[, `:=`(nonzero_frequency = sfreq$nonzero_frequency[match(gene, sfreq$gene)], selection_frequency = sfreq$selection_frequency[match(gene, sfreq$gene)])]
stab[is.na(nonzero_frequency), nonzero_frequency := 0][is.na(selection_frequency), selection_frequency := 0]
stab[, rank := frank(-nonzero_frequency, ties.method = 'first')]
fwrite(stab[, .(gene, standardized_coefficient, nonzero_frequency, selection_frequency)], file.path(out_src, 'SFig1i_final_gene_stability.csv'))
sl <- melt(stab[, .(gene, rank, Direction, `In candidate pool` = selection_frequency, `Nonzero coefficient` = nonzero_frequency)],
           id.vars = c('gene', 'rank', 'Direction'), variable.name = 'Measure', value.name = 'frac')
sg <- ggplot(sl, aes(rank, frac)) +
  geom_line(aes(group = gene), colour = grid_col, linewidth = .3) +
  geom_point(data = sl[Measure == 'In candidate pool'], shape = 3, size = .8, colour = ink2) +
  geom_point(data = sl[Measure == 'Nonzero coefficient'], aes(colour = Direction), size = 1.1) +
  scale_colour_manual(NULL, values = sign_col) +
  scale_y_continuous('Fraction of 120\nrepeated outer fits', limits = c(0, 1.05), breaks = c(0, .5, 1)) +
  labs(x = 'Final HVG–EN nonzero genes (ranked)', subtitle = '+ in candidate gene pool  ·  dot nonzero coefficient') +
  theme_pub(6.5) + theme(legend.position = 'right', plot.subtitle = element_text(size = 5, colour = ink2))

## j | signed error vs age, stage A --------------------------------------
sh <- ggplot(pr[analysis == 'A'], aes(age_week, err)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = ink2) +
  geom_point(aes(colour = Fold, shape = Fold), size = .9, alpha = .85) +
  geom_smooth(method = 'loess', formula = y ~ x, se = FALSE, colour = ink, linewidth = .45, span = .9) +
  scale_colour_manual(NULL, values = fold_col) + scale_shape_manual(NULL, values = fold_shape) +
  scale_x_continuous('Chronological age (weeks)', breaks = c(0, 26, 52, 78)) +
  facet_wrap(~Recipe, nrow = 1) + labs(y = 'Pred. − chron.\n(weeks)') +
  theme_pub(6.5) + theme(legend.position = 'right')

sfig <- tag(sa0, 'a') / tag(sa, 'b') / tag(sb, 'c') / tag(sc, 'd') / tag(se, 'e') / (tag(sd, 'f') | tag(sr, 'g')) /
  (tag(sf, 'h') | tag(sg, 'i')) / tag(sh, 'j') +
  plot_layout(heights = c(.5, .85, .8, .52, .72, .85, .55, .45))
for (x in list(list(sa0, 'SFig1a_final_lambda', 183, 38), list(sa, 'SFig1b_final_landscape', 183, 40), list(sb, 'SFig1c_stage_lambda', 183, 55),
               list(sc, 'SFig1d_hyperparameters', 183, 35), list(se, 'SFig1e_hvg_training_rule', 183, 42), list(sd, 'SFig1f_fixed_vs_tuned', 90, 55),
               list(sr, 'SFig1g_repeats', 93, 55), list(sf, 'SFig1h_N_frequency', 80, 35), list(sg, 'SFig1i_gene_stability', 103, 35),
               list(sh, 'SFig1j_error_vs_age', 183, 32))) save_panel(x[[1]], x[[2]], x[[3]], x[[4]])
W2 <- 183; H2 <- 278
ggsave(file.path(out_final, 'Supplementary_Figure1_model_development.pdf'), sfig, width = W2, height = H2, units = 'mm', device = cairo_pdf)
ggsave(file.path(out_final, 'Supplementary_Figure1_model_development.png'), sfig, width = W2, height = H2, units = 'mm', dpi = 300)
cat('done\n')
