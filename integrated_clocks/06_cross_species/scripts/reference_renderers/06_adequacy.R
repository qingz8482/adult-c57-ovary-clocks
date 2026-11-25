palette_age <- c(
  "20-29" = "#4C78A8", "30-39" = "#72A0C1", "40-49" = "#9ECAC3",
  "50-59" = "#E3B65F", "60-69" = "#E07A5F", "70-79" = "#A8434C"
)
palette_contract <- c(
  neutral_dark = "#272727", neutral_mid = "#767676", neutral_light = "#D8D8D8",
  signal_blue = "#3182BD", signal_teal = "#33B5A5", accent_red = "#D24B40",
  accent_orange = "#E28E2C", elasticnet = "#2878B5", pca_glmnet = "#D65F5F"
)

theme_nature_contract <- function(base_size = 7, base_family = "sans") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.7, face = "bold"),
      panel.grid = element_blank()
    )
}
theme_set(theme_nature_contract())

qc<-readRDS(file.path(old,'3_gtex_age_signal_adequacy/derived_data/gtex_age_signal_qc_objects.rds'));audit<-copy(qc$audit)
for(ff in c('elasticnet','pca_glmnet')){d<-cp[cohort=='GTEx'&model_family==ff];ii<-match(audit$sample_id,d$sample_id);stopifnot(!anyNA(ii));audit[[paste0('existing_score_',ff)]]<-d$cross_species_mouse_ovarian_clock_score[ii];audit[[paste0('projection_ood_distance_',ff)]]<-d$projection_ood_distance[ii]}
pcs<-data.table(sample_id=audit$sample_id,qc$pca$x);pve<-qc$pve
p_string<-function(x)ifelse(x<.001,format(x,scientific=TRUE,digits=2),sprintf('%.3f',x))
pdat <- merge(pcs[, 1:3], audit[, .(sample_id, age_bracket, rin, ischemic_time_minutes, sequencing_center,
  sequencing_batch, expressed_genes_count_ge10, projection_ood_distance_elasticnet,
  projection_ood_distance_pca_glmnet, sensitivity_technical_qc_pass)], by = "sample_id")
axis_labs <- c(x = sprintf("PC1 (%.1f%%)", 100 * pve[1]), y = sprintf("PC2 (%.1f%%)", 100 * pve[2]))
plot_cont <- function(v, title, low = "#D8E9F3", high = "#B2182B") {
  ggplot(pdat, aes(PC1, PC2, colour = .data[[v]])) + geom_point(size = 1.4, alpha = .9) +
    scale_colour_gradient(low = low, high = high, na.value = "#BDBDBD") +
    labs(x = axis_labs["x"], y = axis_labs["y"], colour = title, title = title)
}
p_age <- ggplot(pdat, aes(PC1, PC2, colour = age_bracket)) + geom_point(size = 1.4, alpha = .9) +
  scale_colour_manual(values = palette_age, drop = FALSE) + labs(x = axis_labs["x"], y = axis_labs["y"], colour = "Age", title = "Age bracket")
p_rin <- plot_cont("rin", "RIN")
p_isch <- plot_cont("ischemic_time_minutes", "Ischemic time")
p_det <- plot_cont("expressed_genes_count_ge10", "Detected genes")
p_ood_el <- plot_cont("projection_ood_distance_elasticnet", "ElasticNet OOD")
p_ood_pc <- plot_cont("projection_ood_distance_pca_glmnet", "PCA/glmnet OOD")
p_center <- ggplot(pdat, aes(PC1, PC2, colour = sequencing_center)) + geom_point(size = 1.4, alpha = .9) +
  scale_colour_brewer(palette = "Dark2") + labs(x = axis_labs["x"], y = axis_labs["y"], colour = "Center", title = "Sequencing center")

save_pub_r(p_ood_el, paudit("figures/pca_by_elasticnet_ood_distance"), 90, 75)
save_pub_r(p_ood_pc, paudit("figures/pca_by_pca_glmnet_ood_distance"), 90, 75)
p_grid <- (p_age | p_rin | p_isch) / (p_center | p_det | p_ood_el) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8), legend.position = "right")
save_pub_r(p_grid, paudit("figures/expression_pca_qc_panel"), 183, 118)

gtex<-readRDS(file.path(old,'derived_data/gtex/gtex_ovary_prepared.rds'));logcpm<-gtex$logcpm
ordered_primary<-readRDS(file.path(old,'3_gtex_age_signal_adequacy/derived_data/gtex_age_de_objects.rds'))$ordered_primary$table
stat_primary<-ordered_primary$signed_quasi_likelihood_t[match(rownames(logcpm),ordered_primary$human_gene_id_stable)]
symbols<-gtex$age_effects$human_gene_symbol[match(rownames(logcpm),gtex$age_effects$human_gene_id_stable)]
stagea<-copy(master0);stagea[,deployment_shared_artifact:=hvg_selected&pca_selected]
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
  elasticnet_hvg_3000 = unique(stagea[elasticnet_hvg_3000 == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable]),
  elasticnet_nonzero = unique(stagea[elasticnet_nonzero == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable]),
  pca_versionB_3000 = unique(stagea[pca_versionB_3000 == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable]),
  artifact_shared_current = unique(stagea[deployment_shared_artifact == TRUE & !is.na(human_gene_id_stable), human_gene_id_stable])
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
      sel <- vapply(strata_target, function(s) pools[[s]][sample.int(length(pools[[s]]), 1L)], integer(1))
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
write_tsv(feature_summary, file.path(ref_tables,"table_clock_feature_set_age_information.tsv"))
write_tsv(rbindlist(feature_null, fill = TRUE), file.path(ref_tables,"table_clock_feature_set_matched_random_nulls.tsv"))

age_counts <- fread(paudit("results/table_gtex_ovary_age_bin_counts.tsv"))
age_counts[, age_bracket := factor(age_bracket, levels = age_levels, ordered = TRUE)]
setorder(age_counts, age_bracket)
ord <- fread(paudit("results/table_ordered_age_trend_summary.tsv"))
perm <- fread(paudit("results/table_age_label_permutation_adequacy_summary.tsv"))
perm_draw <- fread(paudit("results/table_age_label_permutation_draws.tsv"))
omni <- fread(paudit("results/table_age_omnibus_nonlinear_summary.tsv"))
pair_counts <- fread(paudit("results/table_pairwise_DEG_counts.tsv"))
pair_support <- fread(paudit("results/table_pairwise_reference_support.tsv"))
vp <- fread(paudit("results/table_multivariate_variance_partition.tsv"))
permanova <- fread(paudit("results/table_permanova_age_structure.tsv"))
disp <- fread(paudit("results/table_multivariate_dispersion_test.tsv"))
cent_test <- fread(paudit("results/table_centroid_distance_age_separation_test.tsv"))
path <- fread(paudit("results/table_reactome_preranked_pathway_enrichment.tsv"))
path_stab <- fread(paudit("results/table_reactome_pathway_stability.tsv"))
features <- copy(feature_summary)
pred <- fread(paudit("results/table_nested_ridge_age_prediction_summary.tsv"))
pred_sample <- fread(paudit("results/table_nested_ridge_age_prediction.tsv"))
hist <- fread(paudit("results/table_histology_pathology_age_balance.tsv"))
hist_clock <- fread(paudit("results/table_clock_slopes_histology_sensitivity.tsv"))
balance <- fread(paudit("results/table_gtex_ovary_metadata_balance.tsv"))
outliers <- fread(paudit("results/table_predefined_technical_qc_outliers.tsv"))
agonal <- fread(paudit("results/table_ordered_age_trend_donor_death_sensitivity_summary.tsv"))
feature_null <- rbindlist(feature_null)

primary <- ord[analysis_set == "primary_all_193"]
qcpass <- ord[analysis_set == "sensitivity_predefined_technical_qc_pass"]
perm_p <- perm[statistic == "sum_squared_age_t", empirical_p]
perm_obs <- perm[statistic == "sum_squared_age_t", observed]
perm_q975 <- perm[statistic == "sum_squared_age_t", null_q975]
qc_retention <- qcpass$fdr_lt_0_05 / primary$fdr_lt_0_05
agonal_retention <- agonal$fdr_lt_0_05 / primary$fdr_lt_0_05
max_exclusion_fraction <- max(outliers$n_any_qc_exclusion / outliers$n)
path_primary_n <- path[specification == "primary_edgeR_technical_adjusted", sum(fdr < .05)]
path_stable_n <- path_stab[specifications_fdr_lt_0_05 == 3 & direction_consistent == TRUE, .N]
pred_cor <- pred[metric == "correlation"]
pred_mae <- pred[metric == "mae_years"]

## Integrated quantitative-grid audit figure ----------------------------------
p_perm <- ggplot(perm_draw, aes(sum_squared_age_t)) + geom_histogram(bins = 45, fill = "#D8D8D8", colour = "white") +
  geom_vline(xintercept = perm_obs, colour = palette_contract["accent_red"], linewidth = .8) +
  annotate("text", x = perm_obs, y = Inf, label = paste0("Observed\np=", p_string(perm_p)), vjust = 1.2, hjust = 1.05, size = 2.4, colour = palette_contract["accent_red"]) +
  labs(x = "Sum of squared age t statistics", y = "Permutations", title = "Age structure exceeds the stratified null")

p_n <- ggplot(age_counts, aes(age_bracket, n_samples, fill = age_bracket)) + geom_col(width = .72) +
  geom_text(aes(label = n_samples), vjust = -.25, size = 2.2) + scale_fill_manual(values = palette_age, guide = "none") +
  coord_cartesian(ylim = c(0, max(age_counts$n_samples) * 1.12)) + labs(x = "Age bracket", y = "Donors", title = "Age-bin counts") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))

ordlong <- melt(ord, id.vars = "analysis_set", measure.vars = c("fdr_lt_0_05", "fdr_lt_0_05_abs_log2fc_gt_0_25", "fdr_lt_0_05_abs_log2fc_gt_0_50"),
  variable.name = "threshold", value.name = "genes")
p_de <- ggplot(ordlong, aes(threshold, genes, fill = analysis_set)) + geom_col(position = position_dodge(width = .78), width = .7) +
  scale_fill_manual(values = c(primary_all_193 = unname(palette_contract["signal_blue"]), sensitivity_predefined_technical_qc_pass = unname(palette_contract["signal_teal"])),
    labels = c("All 193", "QC-pass")) +
  scale_x_discrete(labels = c(fdr_lt_0_05 = "FDR<.05", fdr_lt_0_05_abs_log2fc_gt_0_25 = "+ |LFC|>.25", fdr_lt_0_05_abs_log2fc_gt_0_50 = "+ |LFC|>.50")) +
  labs(x = NULL, y = "Genes", fill = NULL, title = "Ordered-age DE") + theme(legend.position = "bottom", legend.box.margin = margin(0, 0, 0, 0), axis.text.x = element_text(angle = 25, hjust = 1))

vpplot <- vp[!grepl("marginal", component)]
p_vp <- ggplot(vpplot, aes(reorder(component, unique_partial_variance_fraction), unique_partial_variance_fraction, fill = component == "age")) +
  geom_col() + coord_flip() + scale_fill_manual(values = c(`FALSE` = "#BDBDBD", `TRUE` = unname(palette_contract["accent_red"])), guide = "none") +
  labs(x = NULL, y = "Variance fraction", title = "PC variance partition")

p_feat <- ggplot(features[feature_set != "all_expressed_genes"], aes(reorder(feature_set, mean_absolute_age_t), mean_absolute_age_t)) +
  geom_col(fill = palette_contract["signal_teal"]) + geom_hline(yintercept = features[feature_set == "all_expressed_genes", mean_absolute_age_t], linetype = 2, linewidth = .4) +
  coord_flip() + labs(x = NULL, y = "Mean |age t|", title = "Frozen feature-set age signal")

pred_sample[, age_bracket := factor(age_bracket, levels = age_levels)]
p_pred <- ggplot(pred_sample, aes(midpoint_descriptive_only, nested_cv_prediction, colour = age_bracket)) + geom_abline(slope = 1, intercept = 0, colour = "#BDBDBD") +
  geom_point(size = 1.1, alpha = .75) + scale_colour_manual(values = palette_age, guide = "none") +
  labs(x = "Bracket midpoint (descriptive)", y = "Nested-CV prediction", title = sprintf("Nested prediction: r=%.2f", pred_cor$observed))

ptop <- path[specification == "primary_edgeR_technical_adjusted"][order(fdr, -abs(normalized_enrichment_score))][1:min(.N, 6)]
ptop[, pathway_short := fifelse(nchar(pathway) > 46, paste0(substr(pathway, 1, 43), "..."), pathway)]
p_path <- ggplot(ptop, aes(reorder(pathway_short, normalized_enrichment_score), normalized_enrichment_score, fill = direction)) + geom_col() + coord_flip() +
  scale_fill_manual(values = c(Up = unname(palette_contract["accent_red"]), Down = unname(palette_contract["signal_blue"]))) +
  labs(x = NULL, y = "NES", fill = NULL, title = "Top Reactome pathways") + theme(axis.text.y = element_text(size = 5.2), legend.position = "none")

design <- "
AAAA
BBCC
DDEE
FFGG
"
main_fig <- p_perm + p_n + p_de + p_vp + p_feat + p_pred + p_path +
  plot_layout(design = design, heights = c(1.15, 1, 1, 1.15)) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8), plot.title = element_text(size = 7, face = "bold"), plot.margin = margin(3, 4, 3, 3))
save_pub_r(main_fig, paudit("figures/gtex_age_signal_adequacy_integrated"), 183, 190)

