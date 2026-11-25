#!/usr/bin/env Rscript
# Overlap of highly ranked genes between the two locked clocks of release hvg_exclusion_vbpca_20261012:
# HVG-ElasticNet (single-age studies excluded; 79 nonzero genes) and Version B-PCA/glmnet (dense in
# gene space through 10 PCs). Compares the clocks at four levels: input gene pools, top-ranked genes,
# coefficient signs, and correlated expression signal; adds out-of-fold prediction concordance and
# stability across the 30 repeated study partitions. Reads saved outputs only; nothing is refitted.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

a <- commandArgs(FALSE)
script <- gsub('~+~', ' ', sub('^--file=', '', a[grepl('^--file=', a)][1]), fixed = TRUE)
fig_dir <- dirname(dirname(normalizePath(script)))
r4 <- dirname(fig_dir)
r3 <- file.path(dirname(r4), 'rebutal3_batcheffects removal_CV_LOO/model_development')
dev4 <- file.path(r4, 'model_development'); lock <- file.path(r4, 'lockedmodel')
out_dir <- file.path(fig_dir, 'diagnostics', '10_model_gene_overlap')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
source(file.path(lock, 'scripts/core.R'))  # fit_ranking()

# ---------------------------------------------------------------- constants
recipes <- data.table(
  recipe = c('HVG-ElasticNet', 'VersionB-PCA'),
  dev    = c(dev4, r3),
  dir    = c('HVG-ElasticNet', 'VersionB-PCA/glmnet'),
  locked = c('HVG-ElasticNet-exclude-single-age', 'VersionB-PCA-glmnet'),
  label  = c('HVG–EN (single-age excl.)', 'Version B–PCA'))
rec_col <- setNames(c('#2a78d6', '#1baf7a'), recipes$recipe)
ink <- '#0b0b0b'; ink2 <- '#52514e'; grid_col <- '#e6e5e0'; null_col <- '#9a9893'
top_k <- c(20L, 30L, 50L, 79L, 100L, 200L, 500L)
n_null <- 2000L; set.seed(20261013)

theme_pub <- function(base = 7) {
  theme_classic(base_size = base, base_family = 'sans') +
    theme(text = element_text(colour = ink), axis.text = element_text(colour = ink2, size = base - 1),
          axis.line = element_line(linewidth = .3, colour = ink2), axis.ticks = element_line(linewidth = .3, colour = ink2),
          strip.background = element_blank(), strip.text = element_text(face = 'bold', size = base),
          legend.key.size = unit(3, 'mm'), legend.text = element_text(size = base - 1), legend.title = element_text(size = base - 1),
          plot.title = element_text(size = base, face = 'bold'), plot.subtitle = element_text(size = base - 1.5, colour = ink2),
          plot.margin = margin(3, 4, 3, 3))
}
save_plot <- function(p, stem, w, h) {
  ggsave(paste0(stem, '.pdf'), p, width = w, height = h, units = 'mm', device = cairo_pdf)
  ggsave(paste0(stem, '.png'), p, width = w, height = h, units = 'mm', dpi = 300, bg = 'white')
}
tag <- function(p, t) wrap_elements(full = p) + labs(tag = t) +
  theme(plot.tag = element_text(size = 9, face = 'bold'), plot.tag.position = c(0, 1), plot.margin = margin(2, 1, 1, 9))

# ---------------------------------------------------------------- inputs
rpath <- function(r, ...) file.path(recipes[recipe == r, dev], recipes[recipe == r, dir], ...)
models <- setNames(lapply(recipes$locked, function(l) readRDS(file.path(lock, l, 'model.rds'))), recipes$recipe)
for (r in recipes$recipe) stopifnot(identical(models[[r]], readRDS(rpath(r, 'count_tuning/results/final_model/model.rds'))))
en <- models[['HVG-ElasticNet']]; vb <- models[['VersionB-PCA']]
dat <- readRDS(file.path(dev4, 'shared_inputs/prepared_expression.rds'))
stopifnot(identical(unname(tools::md5sum(file.path(dev4, 'shared_inputs/prepared_expression.rds'))),
                    unname(tools::md5sum(file.path(r3, 'shared_inputs/prepared_expression.rds')))))

# Gene universe: genes passing the expression filter in both training sets, so either clock could use them.
U <- intersect(en$representation$filtered_genes, vb$representation$filtered_genes)
en_pool <- intersect(en$representation$genes, U); vb_pool <- intersect(vb$representation$genes, U)
cat(sprintf('input genes outside the shared universe: HVG-EN %d, Version B-PCA %d\n', 3000 - length(en_pool), 3000 - length(vb_pool)))
# Standardized coefficients are in weeks per gene SD for both clocks; Version B-PCA genes carry the
# effective coefficient rotation %*% PC coefficients.
ec <- copy(en$coefficients)[, .(gene, en_coef = standardized_coefficient)]
vc <- copy(vb$coefficients)[, .(gene, vb_coef = standardized_coefficient)]
stopifnot(all.equal(vc$vb_coef, as.vector(vb$representation$rotation %*% vb$pc_coefficients)))
en_nz <- ec[en_coef != 0][order(-abs(en_coef), gene)][, en_rank := .I]
vc[, vb_rank := frank(-abs(vb_coef), ties.method = 'first')]
# Version B feature-score ranking (the selector behind Version B-PCA) on its 136 training profiles.
vbr <- fit_ranking(dat$expression[, vb$training_keys], dat$metadata[match(vb$training_keys, sample_key)], 'VersionB', vb$config)
stopifnot(identical(head(vbr$table$gene, 30), head(vbr$table[gene %in% vb_pool, gene], 30)))
vs <- vbr$table[, .(gene, vb_score_rank = .I, vb_age_t = age_t)]
# Some HVG-EN genes fail the expression filter on the 136 Version B-PCA training profiles (the five
# single-age studies are added), so Version B-PCA could never use them.
en_nz[, in_U := gene %in% U]
nU <- length(U); nE <- nrow(en_nz); nE_U <- en_nz[in_U == TRUE, .N]

# ================================================================= a | input pools
n_shared <- length(intersect(en_pool, vb_pool)); exp_shared <- length(en_pool) * length(vb_pool) / length(U)
p_pool_depl <- phyper(n_shared, length(vb_pool), nU - length(vb_pool), length(en_pool))
pool <- data.table(set = c('HVG–EN input', 'Version B–PCA input', 'Shared, observed', 'Shared, expected by chance'),
                   n = c(length(en_pool), length(vb_pool), n_shared, exp_shared))
pool[, set := factor(set, levels = rev(set))]
fwrite(pool, file.path(out_dir, 'a_input_pool_overlap.tsv'), sep = '\t')
pa <- ggplot(pool, aes(n, set)) +
  geom_col(aes(fill = set), width = .65) +
  geom_text(aes(label = format(round(n), big.mark = ',')), hjust = -.15, size = 2, colour = ink) +
  scale_fill_manual(values = setNames(c(rec_col[[1]], rec_col[[2]], '#34373d', null_col), levels(pool$set)[4:1]), guide = 'none') +
  scale_x_continuous(expand = expansion(mult = c(0, .2))) +
  labs(x = sprintf('Genes (of %s passing both expression filters)', format(nU, big.mark = ',')), y = NULL,
       title = 'Input gene pools',
       subtitle = sprintf('Shared %d vs %.0f expected (%.2f×); depletion P = %.1e', n_shared, exp_shared, n_shared / exp_shared, p_pool_depl)) +
  theme_pub()

# ================================================================= b | top-k overlap
# EN genes are ranked by |coefficient| (79 nonzero); Version B-PCA genes by |effective coefficient| (3,000)
# or by Version B feature score (all scored genes). Overlap of the top k of each is tested against
# a hypergeometric draw from U.
ov <- rbindlist(lapply(top_k, function(k) {
  ke <- min(k, nE); eg <- intersect(en_nz[en_rank <= ke, gene], U)
  rbind(data.table(k, vb_ranking = 'Effective coefficient', vbg = list(vc[vb_rank <= k, gene])),
        data.table(k, vb_ranking = 'Version B feature score', vbg = list(vs[vb_score_rank <= k, gene])))[, {
    g <- vbg[[1]]; o <- intersect(eg, g)
    .(k_en = ke, k_en_in_universe = length(eg), observed = length(o), expected = length(eg) * k / nU,
      p_enrichment = phyper(length(o) - 1, k, nU - k, length(eg), lower.tail = FALSE), genes = paste(o, collapse = ';'))
  }, by = .(k, vb_ranking)]
}))
fwrite(ov, file.path(out_dir, 'b_topk_overlap.tsv'), sep = '\t')
pb <- ggplot(ov, aes(k, observed)) +
  geom_line(aes(y = expected), linetype = 'dashed', colour = null_col, linewidth = .4) +
  geom_line(aes(colour = vb_ranking), linewidth = .5) + geom_point(aes(colour = vb_ranking), size = 1.2) +
  geom_text(data = ov[observed > 0], aes(label = observed, colour = vb_ranking), vjust = -.8, size = 1.8, show.legend = FALSE) +
  scale_colour_manual('Version B–PCA ranked by', values = c('Effective coefficient' = rec_col[[2]], 'Version B feature score' = '#8a5a00')) +
  scale_x_log10(breaks = top_k) + scale_y_continuous(expand = expansion(mult = c(.02, .25))) +
  labs(x = 'Top k genes of each clock (HVG–EN capped at its 79 nonzero)', y = 'Shared genes',
       title = 'Overlap of top-ranked genes', subtitle = 'dashed: expected by chance') +
  theme_pub() + theme(legend.position = 'bottom', legend.direction = 'vertical')

# ================================================================= c | where EN genes sit in Version B rankings
loc <- merge(en_nz, vc, by = 'gene', all.x = TRUE)
loc <- merge(loc, vs, by = 'gene', all.x = TRUE)
loc[, `:=`(in_vb_pool = !is.na(vb_coef), same_sign = sign(en_coef) == sign(vb_coef), sign_vs_age_t = sign(en_coef) == sign(vb_age_t))]
setorder(loc, en_rank)
fwrite(loc, file.path(out_dir, 'c_en_nonzero_genes_in_versionB.tsv'), sep = '\t')
n_scored <- nrow(vs)
# Rank-based test: are EN genes higher in the Version B score ranking than random genes of the scored set?
pw_score <- wilcox.test(loc[!is.na(vb_score_rank), vb_score_rank], setdiff(vs$vb_score_rank, loc$vb_score_rank), alternative = 'less')$p.value
pw_coef <- wilcox.test(loc[in_vb_pool == TRUE, vb_rank], setdiff(vc$vb_rank, loc$vb_rank), alternative = 'less')$p.value
n_in <- loc[in_vb_pool == TRUE, .N]; n_sign <- loc[in_vb_pool == TRUE, sum(same_sign)]
p_sign <- binom.test(n_sign, n_in, .5, alternative = 'greater')$p.value
pc_dt <- loc[!is.na(vb_score_rank)]
pc_dt[, Status := factor(fifelse(!in_vb_pool, 'not in Version B input',
                                 fifelse(same_sign, 'in input, same sign', 'in input, opposite sign')),
                         levels = c('in input, same sign', 'in input, opposite sign', 'not in Version B input'))]
pc_ <- ggplot(pc_dt, aes(vb_score_rank, en_coef)) +
  geom_vline(xintercept = length(vb_pool), linetype = 'dashed', linewidth = .3, colour = ink2) +
  geom_hline(yintercept = 0, linewidth = .3, colour = ink2) +
  geom_point(aes(colour = Status, shape = Status), size = 1.2) +
  geom_text(data = pc_dt[en_rank <= 12 | vb_score_rank <= 300], aes(label = gene), size = 1.6, fontface = 'italic',
            vjust = -.7, colour = ink2, check_overlap = TRUE) +
  annotate('text', x = length(vb_pool), y = Inf, label = 'Version B–PCA\ninput cut-off (3,000)', hjust = 1.05, vjust = 1.1, size = 1.7, colour = ink2) +
  scale_colour_manual(NULL, values = c('in input, same sign' = '#34373d', 'in input, opposite sign' = '#b2182b', 'not in Version B input' = null_col)) +
  scale_shape_manual(NULL, values = c(16, 17, 1)) +
  scale_x_log10(sprintf('Rank in Version B feature score (of %s scored genes)', format(n_scored, big.mark = ','))) +
  labs(y = 'HVG–EN coefficient\n(weeks per SD)', title = 'Where the 79 HVG–EN genes sit in the Version B ranking',
       subtitle = sprintf('%d/79 fail the Version B expression filter (not shown); %d/79 in Version B–PCA input; same sign %d/%d (P = %.1e); rank shift P = %.1e',
                          nE - nE_U, n_in, n_sign, n_in, p_sign, pw_score)) +
  theme_pub() + theme(legend.position = 'bottom')

# ================================================================= d | coefficient agreement, shared pool
sh <- merge(ec, vc, by = 'gene')
sh[, Type := factor(fifelse(en_coef != 0, 'HVG–EN nonzero', 'HVG–EN zero'), levels = c('HVG–EN nonzero', 'HVG–EN zero'))]
rho_all <- cor(sh$en_coef, sh$vb_coef, method = 'spearman')
fwrite(sh, file.path(out_dir, 'd_shared_pool_coefficients.tsv'), sep = '\t')
pd <- ggplot(sh, aes(vb_coef, en_coef)) +
  geom_hline(yintercept = 0, linewidth = .3, colour = ink2) + geom_vline(xintercept = 0, linewidth = .3, colour = ink2) +
  geom_point(data = sh[Type == 'HVG–EN zero'], colour = null_col, size = .5, alpha = .5) +
  geom_point(data = sh[Type == 'HVG–EN nonzero'], colour = rec_col[[1]], size = 1.3) +
  geom_text(data = sh[Type == 'HVG–EN nonzero'][order(-abs(en_coef))][1:10], aes(label = gene), size = 1.6,
            fontface = 'italic', vjust = -.7, colour = ink2, check_overlap = TRUE) +
  labs(x = 'Version B–PCA effective coefficient (weeks per SD)', y = 'HVG–EN coefficient\n(weeks per SD)',
       title = sprintf('Coefficients of the %d genes in both inputs', nrow(sh)),
       subtitle = sprintf('blue: HVG–EN nonzero (%d); grey: zero in HVG–EN', sh[Type == 'HVG–EN nonzero', .N])) +
  theme_pub()

# ================================================================= e | correlated signal (proxy genes)
# Elastic net keeps one gene of a correlated group, so distinct genes can carry the same signal. For each
# HVG-EN nonzero gene: max |Pearson r| with the top 100 Version B-PCA genes across all 136 profiles,
# against random 79-gene draws from the HVG-EN input pool.
X <- t(dat$expression[U, vb$training_keys]); stopifnot(!anyNA(X))
vb_top100 <- vc[vb_rank <= 100, gene]
max_r <- function(g) apply(abs(cor(X[, g, drop = FALSE], X[, vb_top100])), 1, max)
en_u <- en_nz[in_U == TRUE]
obs_r <- max_r(en_u$gene)
null_pool <- setdiff(en_pool, en_nz$gene)
null_med <- replicate(n_null, median(max_r(sample(null_pool, nE_U))))
p_proxy <- (1 + sum(null_med >= median(obs_r))) / (1 + n_null)
prox <- data.table(gene = en_u$gene, en_coef = en_u$en_coef, max_abs_r_vb_top100 = obs_r,
                   partner = vb_top100[apply(abs(cor(X[, en_u$gene], X[, vb_top100])), 1, which.max)])
fwrite(prox, file.path(out_dir, 'e_en_genes_max_correlation_with_vb_top100.tsv'), sep = '\t')
fwrite(data.table(draw = seq_len(n_null), median_max_abs_r = null_med), file.path(out_dir, 'e_null_median_max_correlation.tsv'), sep = '\t')
pe <- ggplot() +
  geom_histogram(data = data.table(x = null_med), aes(x), bins = 40, fill = null_col, colour = NA) +
  geom_vline(xintercept = median(obs_r), colour = rec_col[[1]], linewidth = .6) +
  annotate('text', x = median(obs_r), y = Inf, label = sprintf(' observed\n median %.2f', median(obs_r)), hjust = 0, vjust = 1.2, size = 1.8, colour = rec_col[[1]]) +
  labs(x = sprintf('Median over %d genes of max |r| with Version B–PCA top 100', nE_U),y = sprintf('Random draws\n(n = %d)', n_null),
       title = 'Correlated signal beyond shared genes',
       subtitle = sprintf('vs random HVG-input genes: P = %.2f\n%d/%d HVG–EN genes have a partner with |r| ≥ 0.7',
                          p_proxy, sum(obs_r >= .7), nE_U)) +
  theme_pub()

# ================================================================= f | prediction concordance, Stage A OOF
pr <- rbindlist(lapply(recipes$recipe, function(r) fread(rpath(r, 'count_tuning/results/predictions_A.csv'))[, .(sample_key, study, age_week, pred, recipe)]))
pw <- dcast(pr, sample_key + study + age_week ~ recipe, value.var = 'pred')
setnames(pw, recipes$recipe, c('en', 'vb'))
stopifnot(nrow(pw) == 136, !anyNA(pw))
r_pred <- cor(pw$en, pw$vb); r_res <- cor(pw$en - pw$age_week, pw$vb - pw$age_week)
r_part <- cor(resid(lm(en ~ age_week, pw)), resid(lm(vb ~ age_week, pw)))
fwrite(pw, file.path(out_dir, 'f_stageA_oof_prediction_pairs.tsv'), sep = '\t')
lim <- range(c(pw$en, pw$vb, pw$age_week))
pf <- ggplot(pw, aes(vb, en)) +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', linewidth = .3, colour = ink2) +
  geom_point(aes(colour = age_week), size = 1.1) +
  scale_colour_viridis_c('Age (wk)', option = 'D', end = .92) +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = 'Version B–PCA prediction (weeks)', y = 'HVG–EN prediction\n(weeks)', title = 'Stage A out-of-fold predictions (136)',
       subtitle = sprintf('r = %.2f; residual r = %.2f; age-adjusted r = %.2f', r_pred, r_res, r_part)) +
  theme_pub() + theme(legend.position = 'right', legend.key.height = unit(4, 'mm'))

# ================================================================= g | robustness: repeated-partition gene sets
fs_en <- fread(rpath('HVG-ElasticNet', 'count_tuning/results/feature_stability_A_repeated.csv'))
fs_vb <- fread(rpath('VersionB-PCA', 'count_tuning/results/feature_stability_A_repeated.csv'))
stopifnot(all(fs_en$n_completed_models == 120), all(fs_vb$n_completed_models == 120))
rob <- rbindlist(lapply(c(.1, .25, .5), function(th) {
  ge <- intersect(fs_en[nonzero_frequency >= th, gene], U); gv <- intersect(fs_vb[selection_frequency >= .9, gene], U)
  o <- intersect(ge, gv)
  data.table(en_nonzero_frequency_min = th, n_en = length(ge), n_vb_stable_input = length(gv), observed = length(o),
             expected = length(ge) * length(gv) / nU, p_enrichment = phyper(length(o) - 1, length(gv), nU - length(gv), length(ge), lower.tail = FALSE),
             genes = paste(sort(o), collapse = ';'))
}))
fwrite(rob, file.path(out_dir, 'g_repeated_partition_overlap.tsv'), sep = '\t')
rl <- melt(rob[, .(th = factor(sprintf('≥ %.2f\n(%d genes)', en_nonzero_frequency_min, n_en), levels = sprintf('≥ %.2f\n(%d genes)', en_nonzero_frequency_min, n_en)),
                   Observed = observed, Expected = expected, p_enrichment)], id.vars = c('th', 'p_enrichment'), variable.name = 'What', value.name = 'n')
pg <- ggplot(rl, aes(th, n, fill = What)) +
  geom_col(position = position_dodge(width = .75), width = .7) +
  geom_text(data = rl[What == 'Observed'], aes(label = sprintf('%d\nP = %.0e', n, p_enrichment)), vjust = -.2, size = 1.7, lineheight = .9, colour = ink) +
  scale_fill_manual(NULL, values = c(Observed = '#34373d', Expected = null_col)) +
  scale_y_continuous(expand = expansion(mult = c(0, .35))) +
  labs(x = 'HVG–EN nonzero frequency across 120 repeated outer fits', y = 'Shared with stable\nVersion B–PCA input',
       title = 'Robustness: repeated study partitions',
       subtitle = sprintf('Version B–PCA input genes selected in ≥ 90%% of 120 fits (%d genes)', rob$n_vb_stable_input[1])) +
  theme_pub() + theme(legend.position = 'bottom')

# ---------------------------------------------------------------- assemble
fig <- (tag(pa, 'a') | tag(pb, 'b')) / tag(pc_, 'c') / (tag(pd, 'd') | tag(pe, 'e')) / (tag(pf, 'f') | tag(pg, 'g')) +
  plot_layout(heights = c(.8, 1, .9, .9))
for (x in list(list(pa, 'a_input_pool_overlap', 90, 50), list(pb, 'b_topk_overlap', 90, 60), list(pc_, 'c_en_genes_in_versionB_ranking', 183, 70),
               list(pd, 'd_shared_pool_coefficients', 90, 65), list(pe, 'e_correlated_signal', 90, 55), list(pf, 'f_prediction_concordance', 90, 75), list(pg, 'g_repeated_partition_overlap', 90, 60)))
  save_plot(x[[1]], file.path(out_dir, x[[2]]), x[[3]], x[[4]])
save_plot(fig, file.path(out_dir, 'model_gene_overlap'), 183, 250)

# ---------------------------------------------------------------- summary numbers
summ <- data.table(metric = c(
  'universe_genes', 'en_nonzero', 'en_nonzero_failing_vb_filter', 'input_pool_shared', 'input_pool_expected', 'input_pool_depletion_P',
  'en_nonzero_in_vb_input', 'same_sign_in_vb_input', 'same_sign_P', 'en_vb_score_rank_shift_P', 'en_vb_coef_rank_shift_P',
  'en_nonzero_in_vb_top100_coef', 'en_nonzero_in_vb_top100_score', 'spearman_coef_shared_pool',
  'median_max_r_vb_top100', 'null_median_max_r', 'proxy_P', 'en_genes_partner_r_ge_0.7',
  'oof_pred_r', 'oof_residual_r', 'oof_age_adjusted_r'),
  value = c(nU, nE, nE - nE_U, n_shared, exp_shared, p_pool_depl, n_in, n_sign, p_sign, pw_score, pw_coef,
            ov[k == 100 & vb_ranking == 'Effective coefficient', observed], ov[k == 100 & vb_ranking == 'Version B feature score', observed],
            rho_all, median(obs_r), median(null_med), p_proxy, sum(obs_r >= .7), r_pred, r_res, r_part))
fwrite(summ, file.path(out_dir, 'summary_metrics.tsv'), sep = '\t')
print(summ); print(ov[, !'genes']); print(rob); print(prox[order(-max_abs_r_vb_top100)][1:15])
cat('done\n')
