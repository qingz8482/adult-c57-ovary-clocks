#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common.R"))

message("[2/4] Count-aware ordered, omnibus, nonlinear, and pairwise age analyses")
x <- load_core_inputs()
gtex <- x$gtex; meta <- x$meta
qc_obj <- readRDS(paudit("derived_data/gtex_age_signal_qc_objects.rds"))
audit <- qc_obj$audit[match(meta$sample_id, sample_id)]
counts <- gtex$counts_filtered
logcpm <- gtex$logcpm
gene_map <- copy(gtex$age_effects[, .(human_gene_id_stable, human_gene_symbol, human_mean_expression,
  human_detection_fraction)])
gene_map <- gene_map[match(rownames(counts), human_gene_id_stable)]
stopifnot(identical(gene_map$human_gene_id_stable, rownames(counts)))

## Covariates fixed before inspecting gene-level age results -------------------
analysis_data <- data.frame(
  sample_id = audit$sample_id,
  age_bracket = factor(as.character(audit$age_bracket), levels = age_levels),
  age_order = as.integer(audit$age_bracket) - 1L,
  rin_z = zsafe(audit$rin), ischemia_z = zsafe(audit$ischemic_time_minutes),
  mapping_z = zsafe(audit$mapping_rate), exonic_z = zsafe(audit$exonic_read_rate),
  duplication_z = zsafe(audit$duplicate_mapped_read_rate),
  sequencing_center = factor(audit$sequencing_center),
  technical_qc_pass = audit$sensitivity_technical_qc_pass
)
covariate_terms <- c("rin_z", "ischemia_z", "mapping_z", "exonic_z", "duplication_z", "sequencing_center")
covariate_selection <- data.table(
  covariate = c("RIN", "sample-specific ischemic time", "mapping rate", "exonic-read rate",
    "mapped-read duplication rate", "sequencing center", "extraction batch", "sequencing batch", "library size"),
  primary_status = c(rep("included", 6), "audited_not_adjusted", "audited_not_adjusted", "offset_not_covariate"),
  rationale = c(
    "complete continuous RNA-quality metric", "complete continuous pre-analytic metric",
    "complete continuous mapping-QC metric", "complete continuous RNA-seq composition metric",
    "complete continuous duplication metric", "three-level major technical stratum; included despite one rare level",
    "153 levels among 193 donors; no stable common support for fixed-effect adjustment",
    "118 levels among 193 donors; no stable common support for fixed-effect adjustment",
    "edgeR library-size/TMM offset; not duplicated as a model covariate"
  ), selected_before_gene_results = TRUE
)
write_tsv(covariate_selection, paudit("results/table_primary_technical_covariate_selection.tsv"))

full_rank <- function(design) {
  q <- qr(design)
  keep <- sort(q$pivot[seq_len(q$rank)])
  dropped <- setdiff(seq_len(ncol(design)), keep)
  if (length(dropped)) message("Dropping aliased design columns: ", paste(colnames(design)[dropped], collapse = ", "))
  design[, keep, drop = FALSE]
}

fit_ordered <- function(keep, label) {
  dd <- droplevels(analysis_data[keep, , drop = FALSE])
  dge <- DGEList(counts[, keep, drop = FALSE])
  dge$samples$norm.factors <- gtex$normalization_factors[keep]
  design <- full_rank(model.matrix(~ age_order + rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center, dd))
  stopifnot("age_order" %in% colnames(design))
  fit <- glmQLFit(dge, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = match("age_order", colnames(design)))
  tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "human_gene_id_stable")
  signed_t <- sign(tab$logFC) * sqrt(pmax(tab$F, 0))
  gene_sd <- rowSds(logcpm[, keep, drop = FALSE])
  out <- merge(gene_map, tab, by = "human_gene_id_stable", sort = FALSE)
  out[, `:=`(
    analysis_set = label, ordered_age_coefficient_log2fc_per_decade = logFC,
    standardized_coefficient_per_decade = logFC / gene_sd[match(human_gene_id_stable, rownames(logcpm))],
    empirical_bayes_se = abs(logFC / signed_t), signed_quasi_likelihood_t = signed_t,
    nominal_p = PValue, bh_fdr = FDR, fitted_direction = fifelse(logFC > 0, "increases_with_age", "decreases_with_age"),
    expression_level_logcpm = logCPM, detection_prevalence = human_detection_fraction,
    n_donors = sum(keep), model = "edgeR quasi-likelihood: raw counts ~ ordered age bracket + prespecified technical covariates"
  )]
  keep_cols <- c("human_gene_id_stable", "human_gene_symbol", "analysis_set",
    "ordered_age_coefficient_log2fc_per_decade", "standardized_coefficient_per_decade", "empirical_bayes_se",
    "signed_quasi_likelihood_t", "nominal_p", "bh_fdr", "fitted_direction", "expression_level_logcpm",
    "detection_prevalence", "n_donors", "model")
  list(table = out[, ..keep_cols], fit = fit, design = design, keep = keep)
}

ord_all <- fit_ordered(rep(TRUE, nrow(analysis_data)), "primary_all_193")
ord_qc <- fit_ordered(analysis_data$technical_qc_pass, "sensitivity_predefined_technical_qc_pass")
ordered <- rbindlist(list(ord_all$table, ord_qc$table), fill = TRUE)
write_tsv(ordered, paudit("results/table_ordered_age_trend_full.tsv"))

ordered_counts <- ordered[, .(
  total_genes = .N,
  fdr_lt_0_05 = sum(bh_fdr < .05), fdr_lt_0_10 = sum(bh_fdr < .10),
  fdr_lt_0_05_abs_standardized_gt_0_10 = sum(bh_fdr < .05 & abs(standardized_coefficient_per_decade) > .10),
  fdr_lt_0_05_abs_log2fc_gt_0_10 = sum(bh_fdr < .05 & abs(ordered_age_coefficient_log2fc_per_decade) > .10),
  fdr_lt_0_05_abs_log2fc_gt_0_25 = sum(bh_fdr < .05 & abs(ordered_age_coefficient_log2fc_per_decade) > .25),
  fdr_lt_0_05_abs_log2fc_gt_0_50 = sum(bh_fdr < .05 & abs(ordered_age_coefficient_log2fc_per_decade) > .50)
), by = analysis_set]
write_tsv(ordered_counts, paudit("results/table_ordered_age_trend_summary.tsv"))

## Categorical omnibus and polynomial nonlinear tests --------------------------
message("Fitting categorical omnibus and nonlinear age tests")
dd <- analysis_data
dge_all <- DGEList(counts)
dge_all$samples$norm.factors <- gtex$normalization_factors
design_cat <- full_rank(model.matrix(~ age_bracket + rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center, dd))
fit_cat <- glmQLFit(dge_all, design_cat, robust = TRUE)
age_cat_cols <- grep("^age_bracket", colnames(design_cat))
qlf_cat <- glmQLFTest(fit_cat, coef = age_cat_cols)
cat_tab <- as.data.table(topTags(qlf_cat, n = Inf, sort.by = "none")$table, keep.rownames = "human_gene_id_stable")

dd$age_poly <- ordered(dd$age_bracket, levels = age_levels)
contrasts(dd$age_poly) <- contr.poly(length(age_levels))
design_poly <- full_rank(model.matrix(~ age_poly + rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center, dd))
fit_poly <- glmQLFit(dge_all, design_poly, robust = TRUE)
nonlinear_cols <- grep("^age_poly\\.(Q|C|\\^4|\\^5)", colnames(design_poly))
qlf_nonlin <- glmQLFTest(fit_poly, coef = nonlinear_cols)
nonlin_tab <- as.data.table(topTags(qlf_nonlin, n = Inf, sort.by = "none")$table, keep.rownames = "human_gene_id_stable")
setnames(nonlin_tab, c("F", "PValue", "FDR"), c("nonlinear_F", "nonlinear_p", "nonlinear_fdr"))

group_means <- sapply(age_levels, function(a) rowMeans(logcpm[, dd$age_bracket == a, drop = FALSE]))
diffs <- t(apply(group_means, 1, diff))
mono_up <- rowSums(diffs >= -0.02) == ncol(diffs) & rowSums(diffs > 0.02) >= 3
mono_dn <- rowSums(diffs <= 0.02) == ncol(diffs) & rowSums(diffs < -0.02) >= 3
mixed <- rowSums(diffs > 0.02) > 0 & rowSums(diffs < -0.02) > 0
oldest_jump <- abs(diffs[, 5]) > 2 * (rowMedians(abs(diffs[, 1:4, drop = FALSE])) + 0.02)

omnibus <- merge(gene_map, cat_tab, by = "human_gene_id_stable", sort = FALSE)
omnibus <- merge(omnibus, nonlin_tab[, .(human_gene_id_stable, nonlinear_F, nonlinear_p, nonlinear_fdr)], by = "human_gene_id_stable", sort = FALSE)
prim <- ord_all$table[match(omnibus$human_gene_id_stable, human_gene_id_stable)]
omnibus[, `:=`(
  omnibus_F = F, omnibus_p = PValue, omnibus_fdr = FDR,
  ordered_log2fc_per_decade = prim$ordered_age_coefficient_log2fc_per_decade,
  ordered_p = prim$nominal_p, ordered_fdr = prim$bh_fdr,
  monotonic_effect = mono_up | mono_dn,
  monotonic_direction = fifelse(mono_up, "monotonic_increase", fifelse(mono_dn, "monotonic_decrease", "not_monotonic")),
  nonlinear_effect = nonlinear_fdr < .05,
  inconsistent_directions_across_brackets = mixed,
  age_bin_means_logcpm = apply(group_means, 1, function(z) paste(sprintf("%.3f", z), collapse = ";"))
)]
omnibus[, oldest_bin_specific_effect := omnibus_fdr < .05 & oldest_jump]
write_tsv(omnibus[, .(human_gene_id_stable, human_gene_symbol, omnibus_F, omnibus_p, omnibus_fdr,
  ordered_log2fc_per_decade, ordered_p, ordered_fdr, nonlinear_F, nonlinear_p, nonlinear_fdr,
  monotonic_effect, monotonic_direction, nonlinear_effect, oldest_bin_specific_effect,
  inconsistent_directions_across_brackets, age_bin_means_logcpm)],
  paudit("results/table_age_omnibus_nonlinear_full.tsv"))
write_tsv(omnibus[, .(
  total_genes = .N, omnibus_fdr_lt_0_05 = sum(omnibus_fdr < .05), nonlinear_fdr_lt_0_05 = sum(nonlinear_fdr < .05),
  monotonic_genes = sum(monotonic_effect), oldest_bin_specific_genes = sum(oldest_bin_specific_effect),
  inconsistent_direction_genes = sum(inconsistent_directions_across_brackets)
)], paudit("results/table_age_omnibus_nonlinear_summary.tsv"))

## Pairwise age-bin contrasts --------------------------------------------------
message("Running prespecified pairwise age-bin contrasts")
design_group <- full_rank(model.matrix(~ 0 + age_bracket + rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center, dd))
fit_group <- glmQLFit(dge_all, design_group, robust = TRUE)
age_cols <- setNames(match(paste0("age_bracket", age_levels), colnames(design_group)), age_levels)
if (anyNA(age_cols)) stop("Age-group design columns missing")
n_by <- table(dd$age_bracket)
center_overlap <- sapply(age_levels, function(a) all(c("B1", "C1") %in% unique(as.character(dd$sequencing_center[dd$age_bracket == a]))))
young_primary <- if (n_by["20-29"] >= 20 && center_overlap["20-29"]) "20-29" else "20-39"
write_tsv(data.table(age_bracket = age_levels, n = as.integer(n_by[age_levels]), major_center_common_support = center_overlap,
  selected_primary_young_reference = age_levels == young_primary,
  decision = ifelse(young_primary == "20-29", "20-29 retained: n>=20 and both major sequencing centers represented", "20-39 pooled primary; 20-29 sensitivity")),
  paudit("results/table_pairwise_reference_support.tsv"))

contrast_specs <- list()
add_contrast <- function(name, num, den, family, primary = FALSE) {
  v <- rep(0, ncol(design_group)); names(v) <- colnames(design_group)
  for (a in names(num)) v[age_cols[a]] <- v[age_cols[a]] + num[a]
  for (a in names(den)) v[age_cols[a]] <- v[age_cols[a]] - den[a]
  contrast_specs[[name]] <<- list(vector = v, numerator = paste(names(num), collapse = "+"),
    denominator = paste(names(den), collapse = "+"), family = family, primary = primary,
    n_num = sum(n_by[names(num)]), n_den = sum(n_by[names(den)]))
}

for (old in age_levels[-1]) add_contrast(paste0(old, "_vs_20-29"), setNames(1, old), c("20-29" = 1), "older_vs_20-29", young_primary == "20-29")
for (i in 2:length(age_levels)) add_contrast(paste0(age_levels[i], "_vs_", age_levels[i-1]), setNames(1, age_levels[i]), setNames(1, age_levels[i-1]), "adjacent", TRUE)
w2039 <- as.numeric(n_by[c("20-29", "30-39")]); w2039 <- w2039 / sum(w2039); names(w2039) <- c("20-29", "30-39")
for (old in age_levels[3:6]) add_contrast(paste0(old, "_vs_20-39"), setNames(1, old), w2039, "older_vs_20-39", young_primary == "20-39")
if (!"70-79_vs_20-29" %in% names(contrast_specs)) add_contrast("70-79_vs_20-29", c("70-79" = 1), c("20-29" = 1), "prespecified_extreme", FALSE)
if (!"70-79_vs_20-39" %in% names(contrast_specs)) add_contrast("70-79_vs_20-39", c("70-79" = 1), w2039, "prespecified_extreme", FALSE)

pair_full <- list(); pair_counts <- list()
for (nm in names(contrast_specs)) {
  s <- contrast_specs[[nm]]
  q <- glmQLFTest(fit_group, contrast = s$vector)
  tt <- as.data.table(topTags(q, n = Inf, sort.by = "none")$table, keep.rownames = "human_gene_id_stable")
  tt <- merge(gene_map[, .(human_gene_id_stable, human_gene_symbol)], tt, by = "human_gene_id_stable", sort = FALSE)
  tt[, `:=`(contrast = nm, contrast_family = s$family, primary_reference_contrast = s$primary,
    numerator = s$numerator, denominator = s$denominator, n_numerator = s$n_num, n_denominator = s$n_den)]
  pair_full[[nm]] <- tt
  pair_counts[[nm]] <- rbindlist(lapply(c(0, .10, .25, .50), function(th) tt[, .(
    contrast = nm, contrast_family = s$family, primary_reference_contrast = s$primary,
    numerator = s$numerator, denominator = s$denominator, n_numerator = s$n_num, n_denominator = s$n_den,
    fdr_threshold = .05, abs_log2fc_threshold = th,
    upregulated_deg_count = sum(FDR < .05 & logFC > th), downregulated_deg_count = sum(FDR < .05 & logFC < -th),
    median_absolute_effect_all_genes = median(abs(logFC)), maximum_absolute_effect_all_genes = max(abs(logFC)),
    median_absolute_effect_significant = median(abs(logFC)[FDR < .05 & abs(logFC) > th], na.rm = TRUE),
    maximum_absolute_effect_significant = max(abs(logFC)[FDR < .05 & abs(logFC) > th], na.rm = TRUE)
  )]))
}
pair_full <- rbindlist(pair_full, fill = TRUE)
pair_counts <- rbindlist(pair_counts, fill = TRUE)
pair_counts[!is.finite(median_absolute_effect_significant), median_absolute_effect_significant := NA_real_]
pair_counts[!is.finite(maximum_absolute_effect_significant), maximum_absolute_effect_significant := NA_real_]
write_tsv(pair_full, paudit("results/table_pairwise_DEG_full.tsv"))
write_tsv(pair_counts, paudit("results/table_pairwise_DEG_counts.tsv"))

## Pairwise overlap and directional consistency --------------------------------
sig_sets <- split(pair_full[FDR < .05 & abs(logFC) > .25, human_gene_id_stable], pair_full[FDR < .05 & abs(logFC) > .25, contrast])
overlap_rows <- list()
cn <- names(contrast_specs)
for (i in seq_along(cn)) for (j in i:length(cn)) {
  a <- pair_full[contrast == cn[i], .(human_gene_id_stable, logFC_a = logFC)]
  b <- pair_full[contrast == cn[j], .(human_gene_id_stable, logFC_b = logFC)]
  m <- merge(a, b, by = "human_gene_id_stable")
  sa <- sig_sets[[cn[i]]] %||% character(); sb <- sig_sets[[cn[j]]] %||% character()
  overlap_rows[[length(overlap_rows) + 1L]] <- data.table(contrast_a = cn[i], contrast_b = cn[j],
    deg_overlap = length(intersect(sa, sb)), deg_union = length(union(sa, sb)),
    jaccard = ifelse(length(union(sa, sb)), length(intersect(sa, sb)) / length(union(sa, sb)), NA_real_),
    directional_consistency_all_genes = mean(sign(m$logFC_a) == sign(m$logFC_b)),
    logfc_spearman = cor(m$logFC_a, m$logFC_b, method = "spearman"))
}
write_tsv(rbindlist(overlap_rows), paudit("results/table_pairwise_DEG_overlap_direction.tsv"))

## Pairwise figures ------------------------------------------------------------
count_plot_dat <- pair_counts[abs_log2fc_threshold == .25]
count_plot_dat <- melt(count_plot_dat, id.vars = c("contrast", "contrast_family", "primary_reference_contrast"),
  measure.vars = c("upregulated_deg_count", "downregulated_deg_count"), variable.name = "direction", value.name = "count")
count_plot_dat[direction == "downregulated_deg_count", count := -count]
p_count <- ggplot(count_plot_dat, aes(contrast, count, fill = direction)) + geom_col(width = .72) + coord_flip() +
  scale_fill_manual(values = c(upregulated_deg_count = unname(palette_contract["accent_red"]), downregulated_deg_count = unname(palette_contract["signal_blue"])),
    labels = c("Down", "Up")) + geom_hline(yintercept = 0, linewidth = .25) +
  labs(x = NULL, y = "DEGs (FDR <0.05, |log2FC| >0.25)", fill = NULL, title = "Age-bin differential-expression counts")
save_pub_r(p_count, paudit("figures/age_bin_DEG_count_matrix"), 150, 110)
write_tsv(count_plot_dat, paudit("figures/source_data/age_bin_DEG_count_matrix.tsv"))

primary_contrasts <- names(contrast_specs)[vapply(contrast_specs, function(z) isTRUE(z$primary) && grepl("older_vs", z$family), logical(1))]
mem <- unique(pair_full[contrast %in% primary_contrasts & FDR < .05 & abs(logFC) > .25, .(human_gene_id_stable, contrast)])
if (nrow(mem)) {
  mw <- dcast(mem[, present := TRUE], human_gene_id_stable ~ contrast, value.var = "present", fill = FALSE)
  bool_cols <- setdiff(names(mw), "human_gene_id_stable")
  mw[, intersection := apply(.SD, 1, function(z) paste(bool_cols[as.logical(z)], collapse = " & ")), .SDcols = bool_cols]
  upset <- mw[, .N, by = intersection][order(-N)][1:min(.N, 20)]
  p_upset <- ggplot(upset, aes(reorder(intersection, N), N)) + geom_col(fill = palette_contract["signal_teal"]) +
    coord_flip() + labs(x = "Contrast intersection", y = "Genes", title = "Overlap of primary age contrasts")
  save_pub_r(p_upset, paudit("figures/age_bin_DEG_upset"), 140, 90)
  write_tsv(upset, paudit("figures/source_data/age_bin_DEG_upset.tsv"))
}

pf <- pair_full[contrast %in% primary_contrasts]
topgenes <- pf[, .(max_abs = max(abs(logFC)), min_fdr = min(FDR)), by = .(human_gene_id_stable, human_gene_symbol)][order(min_fdr, -max_abs)][1:min(.N, 60)]
heat <- pf[human_gene_id_stable %in% topgenes$human_gene_id_stable]
heat[, gene_label := fifelse(!is.na(human_gene_symbol) & human_gene_symbol != "", human_gene_symbol, human_gene_id_stable)]
heat[, gene_label := factor(gene_label, levels = rev(unique(topgenes$human_gene_symbol)))]
p_heat <- ggplot(heat, aes(contrast, gene_label, fill = pmax(-1.5, pmin(1.5, logFC)))) + geom_tile() +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-1.5, 1.5), name = "log2FC") +
  labs(x = NULL, y = NULL, title = "Top age-bin expression effects") + theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(size = 5))
save_pub_r(p_heat, paudit("figures/age_bin_logFC_heatmap"), 150, 150)
write_tsv(heat, paudit("figures/source_data/age_bin_logFC_heatmap.tsv"))

## Ordered-trend figures -------------------------------------------------------
vdat <- ord_all$table
vdat[, significant := bh_fdr < .05]
p_vol <- ggplot(vdat, aes(ordered_age_coefficient_log2fc_per_decade, -log10(pmax(nominal_p, 1e-300)), colour = significant)) +
  geom_point(size = .45, alpha = .45) + scale_colour_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = palette_contract["accent_red"])) +
  labs(x = "Ordered-age log2FC per decade", y = "-log10(p)", colour = "FDR <0.05", title = "Transcriptome-wide ordered age trend") +
  theme(legend.position = "top")
p_coef <- ggplot(vdat, aes(ordered_age_coefficient_log2fc_per_decade)) + geom_histogram(bins = 80, fill = palette_contract["signal_blue"], colour = "white", linewidth = .1) +
  geom_vline(xintercept = 0, linewidth = .3) + labs(x = "log2FC per decade", y = "Genes", title = "Coefficient distribution")
p_p <- ggplot(vdat, aes(nominal_p)) + geom_histogram(bins = 50, fill = palette_contract["signal_teal"], colour = "white", linewidth = .1) +
  labs(x = "Nominal p", y = "Genes", title = "P-value distribution")
save_pub_r(p_vol, paudit("figures/ordered_age_volcano"), 95, 80)
save_pub_r(p_coef, paudit("figures/ordered_age_coefficient_distribution"), 90, 70)
save_pub_r(p_p, paudit("figures/ordered_age_pvalue_histogram"), 90, 70)
save_pub_r(p_vol | (p_coef / p_p), paudit("figures/ordered_age_primary_panel"), 183, 105)

## Technically stratified age-label permutation adequacy test ------------------
message("Running 1,000 technically stratified age-label permutations")
N_PERM <- 1000L
set.seed(SEED)
C <- model.matrix(~ rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center, analysis_data)
C <- full_rank(C)
H <- C %*% solve(crossprod(C), t(C))
M <- diag(nrow(C)) - H
E <- logcpm %*% M
gene_var <- rowSums(E^2) / (nrow(C) - ncol(C))
gene_sd <- sqrt(gene_var)
df_res <- nrow(C) - ncol(C) - 1L
strata <- interaction(analysis_data$sequencing_center, drop = TRUE)

gmt <- read_gmt(paudit("external_data/ReactomePathways.gmt"), min_size = 15L, max_size = 500L)
symbols <- gene_map$human_gene_symbol
path_indices <- lapply(gmt, function(gs) which(symbols %in% gs))
path_indices <- path_indices[lengths(path_indices) >= 10L & lengths(path_indices) <= 500L]
pi <- rep(seq_along(path_indices), lengths(path_indices))
pj <- unlist(path_indices, use.names = FALSE)
px <- 1 / sqrt(lengths(path_indices))[pi]
Pmat <- sparseMatrix(i = pi, j = pj, x = px, dims = c(length(path_indices), nrow(E)))

compute_perm_stats <- function(age) {
  ar <- as.numeric(M %*% age)
  denom <- sum(ar^2)
  beta <- as.numeric(E %*% ar) / denom
  r <- as.numeric(E %*% ar) / sqrt(rowSums(E^2) * denom)
  r <- pmax(-.999999, pmin(.999999, r))
  tt <- r * sqrt(df_res / pmax(1e-12, 1 - r^2))
  pp <- 2 * pt(abs(tt), df = df_res, lower.tail = FALSE)
  ff <- p.adjust(pp, "BH")
  std_beta <- beta / gene_sd
  pz <- as.numeric(Pmat %*% tt)
  path_p <- 2 * pnorm(abs(pz), lower.tail = FALSE)
  c(
    n_fdr_significant_genes = sum(ff < .05),
    n_fdr_and_abs_standardized_gt_0_10 = sum(ff < .05 & abs(std_beta) > .10),
    n_abs_standardized_gt_0_10 = sum(abs(std_beta) > .10),
    n_fdr_and_abs_log2fc_gt_0_10 = sum(ff < .05 & abs(beta) > .10),
    n_fdr_and_abs_log2fc_gt_0_25 = sum(ff < .05 & abs(beta) > .25),
    n_fdr_and_abs_log2fc_gt_0_50 = sum(ff < .05 & abs(beta) > .50),
    sum_squared_age_t = sum(tt^2), mean_absolute_age_t = mean(abs(tt)),
    n_significant_pathways_approx = sum(p.adjust(path_p, "BH") < .05)
  )
}

observed_perm_stats <- compute_perm_stats(analysis_data$age_order)
perm_mat <- matrix(NA_real_, nrow = N_PERM, ncol = length(observed_perm_stats),
  dimnames = list(NULL, names(observed_perm_stats)))
for (b in seq_len(N_PERM)) {
  ap <- analysis_data$age_order
  for (s in levels(strata)) {
    ii <- which(strata == s)
    ap[ii] <- sample(ap[ii], length(ii), replace = FALSE)
  }
  perm_mat[b, ] <- compute_perm_stats(ap)
  if (b %% 100L == 0L) message("  permutation ", b, "/", N_PERM)
}
perm_draws <- as.data.table(perm_mat)
perm_draws[, `:=`(iteration = seq_len(N_PERM), seed = SEED, permutation_strata = "sequencing_center")]
write_tsv(perm_draws, paudit("results/table_age_label_permutation_draws.tsv"))

perm_summary <- rbindlist(lapply(names(observed_perm_stats), function(nm) {
  z <- empirical_summary(observed_perm_stats[nm], perm_mat[, nm], tail = "greater")
  z[, `:=`(statistic = nm, n_permutations = N_PERM, seed = SEED,
    permutation_method = "age labels permuted within GTEx sequencing center; all expression and technical covariates fixed",
    pathway_statistic_note = ifelse(grepl("pathway", nm), "approximate Reactome competitive z-score count; definitive cameraPR results reported separately", NA_character_))]
  z
}))
setcolorder(perm_summary, c("statistic", "observed", "null_median", "null_q025", "null_q975", "empirical_p", "monte_carlo_se", "n_permutations", "seed", "permutation_method", "pathway_statistic_note"))
write_tsv(perm_summary, paudit("results/table_age_label_permutation_adequacy_summary.tsv"))

ps <- copy(perm_draws[, .(sum_squared_age_t)])
p_perm <- ggplot(ps, aes(sum_squared_age_t)) + geom_histogram(bins = 45, fill = palette_contract["neutral_light"], colour = "white") +
  geom_vline(xintercept = observed_perm_stats["sum_squared_age_t"], colour = palette_contract["accent_red"], linewidth = .7) +
  labs(x = "Sum of squared age t statistics", y = "Permutations", title = "Technically stratified age-label null")
save_pub_r(p_perm, paudit("figures/age_label_permutation_adequacy"), 105, 75)

saveRDS(list(ordered_primary = ord_all, ordered_qc = ord_qc, omnibus = omnibus,
  pairwise = pair_full, analysis_data = analysis_data, covariate_terms = covariate_terms,
  permutation_observed = observed_perm_stats, permutation_draws = perm_draws,
  young_primary = young_primary), paudit("derived_data/gtex_age_de_objects.rds"))
message("Count-aware age analysis complete. Primary FDR<0.05 genes: ", ordered_counts[analysis_set == "primary_all_193", fdr_lt_0_05],
  "; permutation empirical p(sum t^2)=", perm_summary[statistic == "sum_squared_age_t", empirical_p])
