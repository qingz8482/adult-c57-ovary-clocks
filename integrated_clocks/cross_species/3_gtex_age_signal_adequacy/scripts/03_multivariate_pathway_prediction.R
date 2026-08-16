#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common.R"))

message("[3/4] Centroid, multivariate, pathway, feature-set, histology, and prediction audits")
set.seed(SEED)
x <- load_core_inputs()
gtex <- x$gtex; meta <- x$meta; stagea <- x$stagea
qc_obj <- readRDS(paudit("derived_data/gtex_age_signal_qc_objects.rds"))
de_obj <- readRDS(paudit("derived_data/gtex_age_de_objects.rds"))
audit <- qc_obj$audit[match(meta$sample_id, sample_id)]
ad <- de_obj$analysis_data
logcpm <- gtex$logcpm
ordered_primary <- de_obj$ordered_primary$table
ordered_qc <- de_obj$ordered_qc$table
symbols <- gtex$age_effects$human_gene_symbol[match(rownames(logcpm), gtex$age_effects$human_gene_id_stable)]

## Technically residualized expression for visualization/distance only --------
C <- model.matrix(~ rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center, ad)
C <- C[, sort(qr(C)$pivot[seq_len(qr(C)$rank)]), drop = FALSE]
H <- C %*% solve(crossprod(C), t(C))
E <- logcpm %*% (diag(nrow(C)) - H)
resid_expr <- E + rowMeans(logcpm)

## Age-bin centroid correlations/distances with donor bootstrap ----------------
message("Computing age-bin centroid correlations and bootstrap intervals")
centroid_metric <- function(A, B) {
  ca <- rowMeans(A); cb <- rowMeans(B)
  zca <- zsafe(ca); zcb <- zsafe(cb)
  c(pearson = cor(ca, cb), spearman = cor(ca, cb, method = "spearman"),
    one_minus_pearson = 1 - cor(ca, cb), euclidean_standardized = sqrt(mean((zca - zcb)^2)))
}

raw_centroids <- sapply(age_levels, function(a) rowMeans(logcpm[, audit$age_bracket == a, drop = FALSE]))
res_centroids <- sapply(age_levels, function(a) rowMeans(resid_expr[, audit$age_bracket == a, drop = FALSE]))

B_BOOT <- 200L
centroid_rows <- list(); boot_rows <- list()
for (i in seq_len(length(age_levels) - 1L)) for (j in (i + 1L):length(age_levels)) {
  a <- age_levels[i]; b <- age_levels[j]
  ia <- which(audit$age_bracket == a); ib <- which(audit$age_bracket == b)
  obs_raw <- centroid_metric(logcpm[, ia, drop = FALSE], logcpm[, ib, drop = FALSE])
  obs_res <- centroid_metric(resid_expr[, ia, drop = FALSE], resid_expr[, ib, drop = FALSE])
  bt <- matrix(NA_real_, B_BOOT, 8)
  colnames(bt) <- c(paste0("raw_", names(obs_raw)), paste0("residualized_", names(obs_res)))
  for (k in seq_len(B_BOOT)) {
    xa <- sample(ia, length(ia), TRUE); xb <- sample(ib, length(ib), TRUE)
    bt[k, ] <- c(centroid_metric(logcpm[, xa, drop = FALSE], logcpm[, xb, drop = FALSE]),
      centroid_metric(resid_expr[, xa, drop = FALSE], resid_expr[, xb, drop = FALSE]))
  }
  obs <- c(obs_raw, obs_res); names(obs) <- colnames(bt)
  centroid_rows[[length(centroid_rows) + 1L]] <- rbindlist(lapply(seq_along(obs), function(q) data.table(
    age_bin_a = a, age_bin_b = b, age_bin_separation = j - i, metric = names(obs)[q], observed = obs[q],
    bootstrap_ci_low = quantile(bt[, q], .025, na.rm = TRUE), bootstrap_ci_high = quantile(bt[, q], .975, na.rm = TRUE),
    bootstrap_replicates = B_BOOT, n_a = length(ia), n_b = length(ib)
  )))
  boot_rows[[length(boot_rows) + 1L]] <- data.table(age_bin_a = a, age_bin_b = b, iteration = seq_len(B_BOOT), bt)
}
centroid_tab <- rbindlist(centroid_rows)
write_tsv(centroid_tab, paudit("results/table_age_bin_centroid_correlations_distances.tsv"))
write_tsv(rbindlist(boot_rows), paudit("results/table_age_bin_centroid_bootstrap_draws.tsv"))

## Distance versus age-bin separation; label permutation among six centroids ---
dist_test <- centroid_tab[grepl("one_minus|euclidean", metric), {
  rho <- cor(age_bin_separation, observed, method = "spearman")
  labs <- seq_along(age_levels); names(labs) <- age_levels
  null <- replicate(10000L, {
    lp <- setNames(sample(unname(labs)), names(labs))
    sep <- abs(lp[age_bin_a] - lp[age_bin_b])
    cor(sep, observed, method = "spearman")
  })
  .(spearman_rho_distance_vs_separation = rho,
    empirical_p = (1 + sum(null >= rho)) / (length(null) + 1),
    null_q025 = quantile(null, .025), null_q975 = quantile(null, .975), permutations = 10000L, seed = SEED)
}, by = metric]
write_tsv(dist_test, paudit("results/table_centroid_distance_age_separation_test.tsv"))

## Sample-level within/adjacent/distant similarity -----------------------------
sample_cor_p <- cor(logcpm, method = "pearson")
sample_cor_s <- cor(logcpm, method = "spearman")
sample_cor_pr <- cor(resid_expr, method = "pearson")
sample_cor_sr <- cor(resid_expr, method = "spearman")
ij <- which(upper.tri(sample_cor_p), arr.ind = TRUE)
sep <- abs(as.integer(audit$age_bracket[ij[, 1]]) - as.integer(audit$age_bracket[ij[, 2]]))
sim <- data.table(sample_a = audit$sample_id[ij[, 1]], sample_b = audit$sample_id[ij[, 2]], age_bin_separation = sep,
  similarity_class = factor(fifelse(sep == 0, "within_bin", fifelse(sep == 1, "adjacent_bin", "distant_bin")),
    levels = c("within_bin", "adjacent_bin", "distant_bin")),
  raw_pearson = sample_cor_p[ij], raw_spearman = sample_cor_s[ij],
  residualized_pearson = sample_cor_pr[ij], residualized_spearman = sample_cor_sr[ij])
write_tsv(sim, paudit("results/table_sample_level_age_bin_similarity.tsv"))
sim_summary <- melt(sim, id.vars = c("sample_a", "sample_b", "age_bin_separation", "similarity_class"),
  measure.vars = c("raw_pearson", "raw_spearman", "residualized_pearson", "residualized_spearman"),
  variable.name = "metric", value.name = "similarity")[, .(
    n_pairs = .N, mean = mean(similarity), sd = sd(similarity), median = median(similarity),
    q25 = quantile(similarity, .25), q75 = quantile(similarity, .75)
  ), by = .(metric, similarity_class)]
write_tsv(sim_summary, paudit("results/table_sample_level_age_bin_similarity_summary.tsv"))

## Multivariate variance partition, PERMANOVA, dispersion ----------------------
message("Running multivariate age-structure tests")
Y <- qc_obj$pca$x[, 1:20, drop = FALSE]
dat <- data.frame(
  age = factor(as.character(audit$age_bracket), levels = age_levels), rin = zsafe(audit$rin),
  ischemia = zsafe(audit$ischemic_time_minutes), mapping = zsafe(audit$mapping_rate),
  exonic = zsafe(audit$exonic_read_rate), duplication = zsafe(audit$duplicate_mapped_read_rate),
  center = factor(audit$sequencing_center), seq_batch = factor(audit$sequencing_batch),
  lib_batch = factor(audit$library_extraction_batch)
)
sse_mv <- function(formula, data = dat, response = Y) sum(residuals(lm(formula, data = data))^2)
sse0 <- sum(scale(Y, center = TRUE, scale = FALSE)^2)
full_formula <- Y ~ age + rin + ischemia + mapping + exonic + duplication + center
full_sse <- sse_mv(full_formula)
blocks <- list(age = "age", RIN = "rin", ischemic_time = "ischemia", mapping_rate = "mapping",
  exonic_rate = "exonic", duplication_rate = "duplication", sequencing_center = "center")
vp <- rbindlist(lapply(names(blocks), function(nm) {
  reduced <- update(full_formula, paste(". ~ . -", blocks[[nm]]))
  ss <- sse_mv(reduced) - full_sse
  data.table(component = nm, unique_partial_variance_fraction = ss / sse0)
}))
unique_sum <- sum(vp$unique_partial_variance_fraction)
vp <- rbind(vp,
  data.table(component = "shared_explained_overlap", unique_partial_variance_fraction = max(0, (sse0 - full_sse) / sse0 - unique_sum)),
  data.table(component = "residual", unique_partial_variance_fraction = full_sse / sse0))

## Sparse batches are reported as marginal adjusted R2, not added to primary partition.
marginal_adj_r2 <- function(v) {
  fit <- lm(Y ~ v)
  rank_x <- qr(model.matrix(fit))$rank
  1 - (sum(residuals(fit)^2) / (nrow(Y) - rank_x)) / (sse0 / (nrow(Y) - 1))
}
batch_rows <- data.table(component = c("sequencing_batch_marginal", "library_batch_marginal"),
  unique_partial_variance_fraction = c(marginal_adj_r2(dat$seq_batch), marginal_adj_r2(dat$lib_batch)),
  note = c("118 sparse levels; marginal adjusted R2 only", "153 sparse levels; marginal adjusted R2 only"))
vp[, note := "unique partial fraction in prespecified multivariate model"]
vp <- rbind(vp, batch_rows, fill = TRUE)
write_tsv(vp, paudit("results/table_multivariate_variance_partition.tsv"))

df_age <- nlevels(dat$age) - 1L
df_res <- nrow(dat) - qr(model.matrix(full_formula, dat))$rank
reduced_age_sse <- sse_mv(update(full_formula, . ~ . - age))
age_ss <- reduced_age_sse - full_sse
F_obs <- (age_ss / df_age) / (full_sse / df_res)
R2_partial <- age_ss / reduced_age_sse
strata <- dat$center
F_null <- numeric(999L)
for (b in seq_along(F_null)) {
  dp <- dat
  for (s in levels(strata)) {
    ii <- which(strata == s); dp$age[ii] <- sample(dp$age[ii], length(ii), FALSE)
  }
  fs <- sse_mv(full_formula, data = dp)
  rs <- sse_mv(update(full_formula, . ~ . - age), data = dp)
  F_null[b] <- ((rs - fs) / df_age) / (fs / df_res)
}
permanova <- data.table(test = "partial multivariate age effect on expression PCs 1-20",
  pseudo_F = F_obs, partial_R2 = R2_partial, empirical_p = (1 + sum(F_null >= F_obs)) / (length(F_null) + 1),
  null_median = median(F_null), null_q025 = quantile(F_null, .025), null_q975 = quantile(F_null, .975),
  permutations = length(F_null), permutation_strata = "sequencing_center", seed = SEED)
write_tsv(permanova, paudit("results/table_permanova_age_structure.tsv"))

cent_y <- rowsum(Y, dat$age) / as.numeric(table(dat$age))
disp <- sqrt(rowSums((Y - cent_y[dat$age, , drop = FALSE])^2))
disp_fit <- anova(lm(disp ~ age, dat))
disp_F <- disp_fit[1, "F value"]
disp_null <- numeric(999L)
for (b in seq_along(disp_null)) {
  gp <- dat$age
  for (s in levels(strata)) { ii <- which(strata == s); gp[ii] <- sample(gp[ii]) }
  cc <- rowsum(Y, gp) / as.numeric(table(gp))
  ddv <- sqrt(rowSums((Y - cc[gp, , drop = FALSE])^2))
  disp_null[b] <- anova(lm(ddv ~ gp))[1, "F value"]
}
dispersion <- data.table(test = "homogeneity of multivariate dispersion", F = disp_F,
  empirical_p = (1 + sum(disp_null >= disp_F)) / (length(disp_null) + 1),
  null_q025 = quantile(disp_null, .025), null_q975 = quantile(disp_null, .975), permutations = length(disp_null), seed = SEED)
write_tsv(dispersion, paudit("results/table_multivariate_dispersion_test.tsv"))

## Leave one age bin out -------------------------------------------------------
loo <- rbindlist(lapply(age_levels, function(drop_bin) {
  keep <- dat$age != drop_bin; d <- droplevels(dat[keep, ]); yy <- Y[keep, , drop = FALSE]
  ff <- yy ~ age + rin + ischemia + mapping + exonic + duplication + center
  s0 <- sum(scale(yy, center = TRUE, scale = FALSE)^2)
  sf <- sum(residuals(lm(ff, d))^2)
  sr <- sum(residuals(lm(update(ff, . ~ . - age), d))^2)
  dfa <- nlevels(d$age) - 1; dfr <- nrow(d) - qr(model.matrix(ff, d))$rank
  data.table(omitted_age_bracket = drop_bin, n_remaining = sum(keep), age_partial_R2 = (sr - sf) / sr,
    pseudo_F = ((sr - sf) / dfa) / (sf / dfr), parametric_p = pf(((sr - sf) / dfa) / (sf / dfr), dfa, dfr, lower.tail = FALSE))
}))
write_tsv(loo, paudit("results/table_leave_one_age_bin_out_multivariate.tsv"))

## Reactome preranked enrichment across specifications -------------------------
message("Running preranked Reactome pathway enrichment")
gmt <- read_gmt(paudit("external_data/ReactomePathways.gmt"), min_size = 15L, max_size = 500L)
indices <- lapply(gmt, function(gs) which(symbols %in% gs))
indices <- indices[lengths(indices) >= 10L & lengths(indices) <= 500L]
stat_primary <- ordered_primary$signed_quasi_likelihood_t[match(rownames(logcpm), ordered_primary$human_gene_id_stable)]
stat_qc <- ordered_qc$signed_quasi_likelihood_t[match(rownames(logcpm), ordered_qc$human_gene_id_stable)]
stat_ageonly <- apply(logcpm, 1, function(z) {
  r <- cor(z, ad$age_order); r * sqrt((nrow(ad) - 2) / pmax(1e-12, 1 - r^2))
})
path_specs <- list(primary_edgeR_technical_adjusted = stat_primary,
  qc_pass_edgeR_technical_adjusted = stat_qc, age_only_logcpm = stat_ageonly)
path_res <- list()
for (nm in names(path_specs)) {
  st <- path_specs[[nm]]; st[!is.finite(st)] <- 0
  cam <- as.data.table(cameraPR(st, indices, use.ranks = TRUE, sort = FALSE), keep.rownames = "pathway")
  cam[, `:=`(specification = nm,
    normalized_enrichment_score = qnorm(pmax(PValue / 2, 1e-300), lower.tail = FALSE) * fifelse(Direction == "Up", 1, -1))]
  cam[, leading_edge_genes := vapply(pathway, function(pw) {
    ii <- indices[[pw]]; zz <- st[ii]; ord <- if (Direction[pathway == pw][1] == "Up") order(zz, decreasing = TRUE) else order(zz)
    paste(head(unique(symbols[ii][ord]), 20), collapse = ";")
  }, character(1))]
  path_res[[nm]] <- cam
}
path_res <- rbindlist(path_res, fill = TRUE)
setnames(path_res, c("NGenes", "Direction", "PValue", "FDR"), c("n_genes", "direction", "p", "fdr"))
write_tsv(path_res, paudit("results/table_reactome_preranked_pathway_enrichment.tsv"))
path_stability <- path_res[, .(
  specifications_tested = .N, specifications_fdr_lt_0_05 = sum(fdr < .05),
  direction_consistent = uniqueN(direction) == 1L,
  min_fdr = min(fdr), median_normalized_enrichment_score = median(normalized_enrichment_score),
  directions = paste(direction, collapse = ";")
), by = pathway]
write_tsv(path_stability, paudit("results/table_reactome_pathway_stability.tsv"))

## Clock-feature-specific age information and matched random sets --------------
message("Comparing frozen clock feature sets with matched random genes")
universe <- data.table(human_gene_id_stable = rownames(logcpm), symbol = symbols,
  mean_expression = rowMeans(logcpm), detection = rowMeans(gtex$counts_filtered >= 10),
  age_t = stat_primary,
  age_fdr = ordered_primary$bh_fdr[match(rownames(logcpm), ordered_primary$human_gene_id_stable)])
universe[, `:=`(
  expr_bin = pmax(1L, pmin(10L, ceiling(frank(mean_expression) / .N * 10))),
  detect_bin = pmax(1L, pmin(5L, ceiling(frank(detection) / .N * 5)))
)]
universe[, stratum := paste(expr_bin, detect_bin, sep = "_")]

flag_sets <- list(
  all_expressed_genes = universe$human_gene_id_stable,
  elasticnet_hvg_2000 = unique(stagea[elasticnet_hvg_2000 == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable]),
  elasticnet_nonzero = unique(stagea[elasticnet_nonzero == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable]),
  pca_versionB_2000 = unique(stagea[pca_versionB_2000 == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable]),
  artifact_shared_297 = unique(stagea[deployment_shared_artifact == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable])
)
flag_sets <- lapply(flag_sets, intersect, y = universe$human_gene_id_stable)
pools <- split(seq_len(nrow(universe)), universe$stratum)
feature_summary <- list(); feature_null <- list()
for (nm in names(flag_sets)) {
  ids <- flag_sets[[nm]]; d <- universe[human_gene_id_stable %in% ids]
  obs <- mean(abs(d$age_t), na.rm = TRUE)
  if (nm == "all_expressed_genes") {
    emp <- NA_real_; nq <- NA_real_; null <- numeric()
  } else {
    strata_target <- d$stratum
    null <- numeric(1000L)
    for (b in seq_along(null)) {
      sel <- vapply(strata_target, function(s) sample(pools[[s]], 1L), integer(1))
      null[b] <- mean(abs(universe$age_t[sel]), na.rm = TRUE)
    }
    emp <- (1 + sum(null >= obs)) / (length(null) + 1)
    nq <- mean(null)
    feature_null[[nm]] <- data.table(feature_set = nm, iteration = seq_along(null), matched_random_mean_absolute_age_t = null)
  }
  feature_summary[[nm]] <- data.table(feature_set = nm, n_genes = nrow(d), mean_absolute_age_t = obs,
    median_absolute_age_t = median(abs(d$age_t), na.rm = TRUE), proportion_fdr_lt_0_05 = mean(d$age_fdr < .05, na.rm = TRUE),
    sum_squared_age_t = sum(d$age_t^2, na.rm = TRUE), matched_random_null_mean = nq,
    matched_random_empirical_p_greater = emp,
    age_informativeness = ifelse(is.na(emp), "reference_universe", ifelse(emp < .05, "more_age_informative_than_matched_random", ifelse(obs < nq, "less_or_equal_age_informative", "equal_or_marginal"))))
}
feature_summary <- rbindlist(feature_summary)
write_tsv(feature_summary, paudit("results/table_clock_feature_set_age_information.tsv"))
write_tsv(rbindlist(feature_null, fill = TRUE), paudit("results/table_clock_feature_set_matched_random_nulls.tsv"))

## Direct histology/pathology balance and exploratory clock adjustment ----------
hist_flags <- c("histology_atrophy_or_postmenopausal", "histology_follicular_material",
  "histology_stromal_predominance", "histology_adipose_predominance", "histology_cystic_or_pathological",
  "histology_cortex_mentioned", "histology_medulla_mentioned", "histology_possible_adjacent_tissue")
hist_balance <- rbindlist(lapply(hist_flags, function(v) {
  tab <- table(audit$age_bracket, audit[[v]])
  tst <- suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000))
  data.table(annotation = v, n_positive = sum(audit[[v]], na.rm = TRUE), cramers_v = cramers_v(tab),
    age_association_p = tst$p.value, test = "Monte Carlo chi-square on direct GTEx pathology-note keyword classification")
}))
hist_balance[, age_association_fdr := p.adjust(age_association_p, "BH")]
write_tsv(hist_balance, paudit("results/table_histology_pathology_age_balance.tsv"))

clock_rows <- list()
for (fam in c("elasticnet", "pca_glmnet")) {
  score_col <- paste0("existing_score_", fam)
  d <- data.frame(score = audit[[score_col]], age = ad$age_order, rin = ad$rin_z, ischemia = ad$ischemia_z,
    mapping = ad$mapping_z, exonic = ad$exonic_z, duplication = ad$duplication_z, center = ad$sequencing_center,
    audit[, ..hist_flags])
  f0 <- lm(score ~ age + rin + ischemia + mapping + exonic + duplication + center, d)
  usable_hist <- hist_flags[vapply(d[, hist_flags, drop = FALSE], function(z) sum(z, na.rm = TRUE) >= 5 && sum(!z, na.rm = TRUE) >= 5, logical(1))]
  f1 <- lm(as.formula(paste("score ~ age + rin + ischemia + mapping + exonic + duplication + center +", paste(usable_hist, collapse = "+"))), d)
  for (nm in c("technical_only", "technical_plus_direct_pathology_flags")) {
    f <- if (nm == "technical_only") f0 else f1
    co <- summary(f)$coef["age", ]
    clock_rows[[length(clock_rows) + 1L]] <- data.table(model_family = fam, model = nm, age_slope_per_bracket = co[1], se = co[2], p = co[4], n = nobs(f), pathology_flags = paste(usable_hist, collapse = ";"))
  }
}
write_tsv(rbindlist(clock_rows), paudit("results/table_clock_slopes_histology_sensitivity.tsv"))
write_tsv(data.table(
  analysis = c("validated_cell_composition_deconvolution", "composition_adjusted_age_DE", "composition_adjusted_clock_slopes"),
  status = "NOT_RUN_NO_VALIDATED_REFERENCE_SUPPLIED",
  rationale = "No validated ovarian bulk-RNA reference signature was present in the completed official input module; expression-derived composition was not inferred and reused as an unlabelled primary covariate."
), paudit("results/table_cell_composition_sensitivity_status.tsv"))

## Optional fully nested ridge/PCA-style age-information prediction audit -------
message("Running optional fully nested ridge prediction adequacy audit")
make_stratified_folds <- function(g, k = 5L, seed = SEED) {
  set.seed(seed); out <- integer(length(g))
  for (lev in unique(g)) { ii <- which(g == lev); out[ii] <- sample(rep(seq_len(k), length.out = length(ii))) }
  out
}
outer_fold <- make_stratified_folds(audit$age_bracket, 5L)
X <- logcpm
y_mid <- audit$age_bracket_midpoint_sensitivity_only

nested_ridge <- function(y, seed_offset = 0L) {
  pred <- rep(NA_real_, length(y))
  for (fo in 1:5) {
    tr <- outer_fold != fo; te <- !tr
    vv <- rowVars(X[, tr, drop = FALSE])
    sel <- order(vv, decreasing = TRUE)[seq_len(min(2000L, length(vv)))]
    inner <- make_stratified_folds(audit$age_bracket[tr], 5L, SEED + seed_offset + fo)
    cv <- cv.glmnet(t(X[sel, tr, drop = FALSE]), y[tr], alpha = 0, family = "gaussian",
      foldid = inner, standardize = TRUE, nlambda = 40, type.measure = "mse")
    pred[te] <- as.numeric(predict(cv, newx = t(X[sel, te, drop = FALSE]), s = "lambda.min"))
  }
  pred
}

pred_obs <- nested_ridge(y_mid)
nearest_bin <- function(z) age_levels[max.col(-abs(outer(z, c(24.5, 34.5, 44.5, 54.5, 64.5, 74.5), "-")), ties.method = "first")]
metrics_pred <- function(pred, y) {
  true_bin <- nearest_bin(y)
  c(correlation = cor(pred, y), mae_years = mean(abs(pred - y)),
    ordinal_exact_accuracy = mean(nearest_bin(pred) == true_bin),
    ordinal_within_one_bin_accuracy = mean(abs(match(nearest_bin(pred), age_levels) - match(true_bin, age_levels)) <= 1))
}
obs_metrics <- metrics_pred(pred_obs, y_mid)

N_PRED_PERM <- 100L
pred_null <- matrix(NA_real_, N_PRED_PERM, length(obs_metrics), dimnames = list(NULL, names(obs_metrics)))
for (b in seq_len(N_PRED_PERM)) {
  yp <- y_mid
  for (s in levels(ad$sequencing_center)) { ii <- which(ad$sequencing_center == s); yp[ii] <- sample(yp[ii]) }
  pred_null[b, ] <- metrics_pred(nested_ridge(yp, 1000L + b * 10L), yp)
  if (b %% 10L == 0L) message("  nested prediction permutation ", b, "/", N_PRED_PERM)
}
pred_summary <- rbindlist(lapply(names(obs_metrics), function(nm) {
  higher <- grepl("correlation|accuracy", nm)
  empirical_summary(obs_metrics[nm], pred_null[, nm], if (higher) "greater" else "less")[,
    `:=`(metric = nm, n_outer_folds = 5L, n_permutations = N_PRED_PERM, seed = SEED,
      leakage_control = "variance selection, scaling, ridge tuning and prediction repeated inside each outer training fold")]
}))
boot_idx <- replicate(2000L, sample.int(length(y_mid), length(y_mid), TRUE))
boot_cor <- apply(boot_idx, 2, function(ii) cor(pred_obs[ii], y_mid[ii]))
boot_mae <- apply(boot_idx, 2, function(ii) mean(abs(pred_obs[ii] - y_mid[ii])))
pred_summary[metric == "correlation", `:=`(bootstrap_ci_low = quantile(boot_cor, .025), bootstrap_ci_high = quantile(boot_cor, .975))]
pred_summary[metric == "mae_years", `:=`(bootstrap_ci_low = quantile(boot_mae, .025), bootstrap_ci_high = quantile(boot_mae, .975))]
write_tsv(data.table(sample_id = audit$sample_id, donor_id = audit$donor_id, age_bracket = audit$age_bracket,
  midpoint_descriptive_only = y_mid, nested_cv_prediction = pred_obs, outer_fold = outer_fold),
  paudit("results/table_nested_ridge_age_prediction.tsv"))
write_tsv(pred_summary, paudit("results/table_nested_ridge_age_prediction_summary.tsv"))
write_tsv(data.table(iteration = seq_len(N_PRED_PERM), pred_null), paudit("results/table_nested_ridge_age_prediction_permutations.tsv"))

## Figures ---------------------------------------------------------------------
cent_plot <- centroid_tab[metric %in% c("raw_one_minus_pearson", "residualized_one_minus_pearson", "raw_euclidean_standardized", "residualized_euclidean_standardized")]
p_cent <- ggplot(cent_plot, aes(age_bin_separation, observed, colour = metric)) + geom_point(size = 1.5) + geom_smooth(method = "lm", se = TRUE, linewidth = .5) +
  scale_colour_manual(values = c("raw_one_minus_pearson" = "#9ECAE1", "residualized_one_minus_pearson" = "#3182BD",
    "raw_euclidean_standardized" = "#F6AE6B", "residualized_euclidean_standardized" = "#D24B40")) +
  labs(x = "Age-bin separation (decades)", y = "Centroid distance", colour = NULL, title = "Centroid distance increases with age separation") + theme(legend.position = "top")
save_pub_r(p_cent, paudit("figures/age_bin_centroid_distance"), 120, 82)

vp_plot <- vp[!grepl("marginal", component)]
p_vp <- ggplot(vp_plot, aes(reorder(component, unique_partial_variance_fraction), unique_partial_variance_fraction, fill = component == "age")) +
  geom_col() + coord_flip() + scale_fill_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = palette_contract["accent_red"]), guide = "none") +
  labs(x = NULL, y = "Variance fraction", title = "Expression variance partition (PC1-20)")
save_pub_r(p_vp, paudit("figures/variance_explained_panel"), 105, 80)

ptop <- path_res[specification == "primary_edgeR_technical_adjusted"][order(fdr, -abs(normalized_enrichment_score))][1:min(.N, 25)]
p_path <- ggplot(ptop, aes(reorder(pathway, normalized_enrichment_score), normalized_enrichment_score, fill = direction)) + geom_col() + coord_flip() +
  scale_fill_manual(values = c(Up = unname(palette_contract["accent_red"]), Down = unname(palette_contract["signal_blue"]))) +
  labs(x = NULL, y = "Normalized enrichment score", fill = NULL, title = "Preranked Reactome age pathways")
save_pub_r(p_path, paudit("figures/reactome_preranked_age_pathways"), 145, 115)

pred_df <- data.table(observed = y_mid, predicted = pred_obs, age_bracket = audit$age_bracket)
p_pred <- ggplot(pred_df, aes(observed, predicted, colour = age_bracket)) + geom_abline(slope = 1, intercept = 0, colour = "#BDBDBD") + geom_point(size = 1.3, alpha = .8) +
  scale_colour_manual(values = palette_age) + labs(x = "Age-bracket midpoint (descriptive years)", y = "Nested-CV prediction", colour = "Age", title = "Dataset-level age-information prediction audit")
save_pub_r(p_pred, paudit("figures/nested_ridge_age_prediction"), 100, 82)

saveRDS(list(residualized_expression = resid_expr, raw_centroids = raw_centroids, residualized_centroids = res_centroids,
  centroid_table = centroid_tab, variance_partition = vp, permanova = permanova, dispersion = dispersion,
  pathway = path_res, feature_summary = feature_summary, prediction_summary = pred_summary),
  paudit("derived_data/gtex_age_multivariate_objects.rds"))
message("Multivariate/pathway/prediction audit complete. PERMANOVA p=", permanova$empirical_p,
  "; primary Reactome FDR<0.05 pathways=", path_res[specification == "primary_edgeR_technical_adjusted", sum(fdr < .05)],
  "; nested-CV age correlation=", fmt(obs_metrics["correlation"]))
