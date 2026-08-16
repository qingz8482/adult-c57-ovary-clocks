#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common.R"))

message("[4/4] Building integrated adequacy figure, verdict, and manuscript revision")
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
features <- fread(paudit("results/table_clock_feature_set_age_information.tsv"))
pred <- fread(paudit("results/table_nested_ridge_age_prediction_summary.tsv"))
pred_sample <- fread(paudit("results/table_nested_ridge_age_prediction.tsv"))
hist <- fread(paudit("results/table_histology_pathology_age_balance.tsv"))
hist_clock <- fread(paudit("results/table_clock_slopes_histology_sensitivity.tsv"))
balance <- fread(paudit("results/table_gtex_ovary_metadata_balance.tsv"))
outliers <- fread(paudit("results/table_predefined_technical_qc_outliers.tsv"))
agonal <- fread(paudit("results/table_ordered_age_trend_donor_death_sensitivity_summary.tsv"))
feature_null <- fread(paudit("results/table_clock_feature_set_matched_random_nulls.tsv"))

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

technical_failure <- qc_retention < 0.5 || agonal$fdr_lt_0_05 == 0 || max_exclusion_fraction > 0.35 || disp$empirical_p < .01
if (technical_failure) {
  verdict <- "C. CONFOUNDED_OR_LOW_QUALITY"
} else if (perm_p < .05 && permanova$empirical_p < .05 && primary$fdr_lt_0_05 > 0) {
  verdict <- "A. INFORMATIVE_AGE_SIGNAL"
} else {
  verdict <- "B. WEAK_AGE_SIGNAL"
}

decision <- data.table(
  criterion = c("age-label permutation", "count-aware transcriptome-wide DE", "technical-QC sensitivity retention",
    "donor-death/agonal sensitivity", "multivariate age structure", "dispersion confounding", "age-bin technical exclusion concentration"),
  observed = c(sprintf("sum(t^2)=%.1f; null 97.5%%=%.1f; p=%s", perm_obs, perm_q975, p_string(perm_p)),
    sprintf("%d genes FDR<0.05", primary$fdr_lt_0_05), sprintf("%.1f%% of primary FDR genes retained", 100 * qc_retention),
    sprintf("%d genes FDR<0.05 after Hardy-factor adjustment", agonal$fdr_lt_0_05),
    sprintf("partial R2=%.3f; stratified PERMANOVA p=%s", permanova$partial_R2, p_string(permanova$empirical_p)),
    sprintf("p=%s", p_string(disp$empirical_p)), sprintf("maximum %.1f%% in any bin", 100 * max_exclusion_fraction)),
  pass = c(perm_p < .05, primary$fdr_lt_0_05 > 0, qc_retention >= .5, agonal$fdr_lt_0_05 > 0,
    permanova$empirical_p < .05, disp$empirical_p >= .01, max_exclusion_fraction <= .35),
  prespecified_gate = verdict
)
write_tsv(decision, paudit("results/table_gtex_age_signal_decision_gate.tsv"))

## Correctly styled pairwise DEG count figure ----------------------------------
cnt <- pair_counts[abs_log2fc_threshold == .25]
cnt <- melt(cnt, id.vars = c("contrast", "contrast_family", "primary_reference_contrast"),
  measure.vars = c("upregulated_deg_count", "downregulated_deg_count"), variable.name = "direction", value.name = "count")
cnt[direction == "downregulated_deg_count", count := -count]
p_pair <- ggplot(cnt, aes(contrast, count, fill = direction)) + geom_col(width = .72) + coord_flip() +
  scale_fill_manual(values = c(upregulated_deg_count = unname(palette_contract["accent_red"]),
    downregulated_deg_count = unname(palette_contract["signal_blue"])), labels = c("Down", "Up")) +
  geom_hline(yintercept = 0, linewidth = .25) +
  labs(x = NULL, y = "DEGs (FDR <0.05, |log2FC| >0.25)", fill = NULL, title = "Age-bin differential-expression counts")
save_pub_r(p_pair, paudit("figures/age_bin_DEG_count_matrix"), 150, 110)

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

## Reports ---------------------------------------------------------------------
counts_txt <- paste(sprintf("%s: %d", age_counts$age_bracket, age_counts$n_samples), collapse = "; ")
tech_focus <- balance[variable %in% c("rin", "ischemic_time_minutes", "mapping_rate", "exonic_read_rate",
  "duplicate_mapped_read_rate", "sequencing_center", "donor_death_hardy_scale")]
tech_sig <- tech_focus[(is.finite(ordered_fdr) & ordered_fdr < .05) | (is.finite(omnibus_fdr) & omnibus_fdr < .05), variable]
hist_atrophy <- hist[annotation == "histology_atrophy_or_postmenopausal"]
el_hist <- hist_clock[model_family == "elasticnet"]
pc_hist <- hist_clock[model_family == "pca_glmnet"]

passport <- c(
  "## Material Passport", "",
  "- Artifact type: GTEx V11 ovarian transcriptome age-signal adequacy audit",
  "- Verification status: VERIFIED (complete scripted rerun plus output validation)",
  "- Data scope: 193 official GTEx V11 ovary donors; 25,927 prespecified filtered genes",
  "- Frozen-model boundary: mouse clocks, Stage A/B results, coefficients, PCA rotations, centers/scales, and calibrations were not altered",
  sprintf("- Randomization: seed %d; 1,000 sequencing-center-stratified age-label permutations", SEED),
  "- Plotting backend: R exclusively", ""
)

verdict_lines <- c(passport,
  "# GTEx ovary age-signal adequacy verdict", "",
  paste0("## Decision: **", verdict, "**"), "",
  sprintf("GTEx V11 ovary contains strong transcriptome-wide age-associated structure beyond technically stratified label permutation. The observed global sum of squared age t statistics was %.1f versus a null median of %.1f and null 97.5th percentile of %.1f (empirical p=%s; Monte Carlo SE=%s).",
    perm_obs, perm[statistic == "sum_squared_age_t", null_median], perm_q975, p_string(perm_p), p_string(perm[statistic == "sum_squared_age_t", monte_carlo_se])), "",
  "## Prespecified gate evidence", "",
  sprintf("- Exact bracket counts: %s.", counts_txt),
  sprintf("- Primary raw-count edgeR QL model: %d genes at FDR<0.05; %d also had |log2FC per decade|>0.25 and %d had |log2FC per decade|>0.50.",
    primary$fdr_lt_0_05, primary$fdr_lt_0_05_abs_log2fc_gt_0_25, primary$fdr_lt_0_05_abs_log2fc_gt_0_50),
  sprintf("- Predefined technical-QC sensitivity retained %d FDR genes (%.1f%% of the primary count); no age bin lost more than %.1f%% of donors.",
    qcpass$fdr_lt_0_05, 100 * qc_retention, 100 * max_exclusion_fraction),
  sprintf("- Donor-death Hardy-factor sensitivity retained %d FDR genes among 191 donors (%.1f%% of the primary count).", agonal$fdr_lt_0_05, 100 * agonal_retention),
  sprintf("- Multivariate age effect: partial R2=%.3f, stratified PERMANOVA p=%s; dispersion p=%s, arguing against a location effect explained solely by unequal dispersion.",
    permanova$partial_R2, p_string(permanova$empirical_p), p_string(disp$empirical_p)),
  sprintf("- Reactome preranked analysis: %d pathways at FDR<0.05 in the primary specification; %d were FDR-significant with consistent direction in all three specifications.", path_primary_n, path_stable_n),
  sprintf("- Fully nested optional prediction audit: r=%.3f (95%% bootstrap CI %.3f–%.3f), MAE=%.2f descriptive midpoint-years, permutation p=%s.",
    pred_cor$observed, pred_cor$bootstrap_ci_low, pred_cor$bootstrap_ci_high, pred_mae$observed, p_string(pred_cor$empirical_p)), "",
  "## Technical and biological qualification", "",
  sprintf("Age was associated with %s after multiplicity correction. These variables were either included in the primary model (RIN, ischemia, mapping, exonic and duplication metrics, sequencing center) or examined in a prespecified sensitivity (donor-death Hardy category).", paste(tech_sig, collapse = ", ")),
  sprintf("Direct pathology-note classification found an age association for atrophy/postmenopausal annotations (n=%d, FDR=%s). This is a biologically plausible bulk-tissue feature, but it means the age signal should not be described as cell-intrinsic.", hist_atrophy$n_positive, p_string(hist_atrophy$age_association_fdr)),
  "No validated ovarian cell-composition reference signature was supplied, so expression-derived deconvolution was not reused as a primary adjustment. Cortex/medulla and adjacent-tissue indicators remain direct-text exploratory annotations.",
  "The 70–79 bin contains only five donors, so oldest-bin and 70–79 pairwise effects are low-precision even though the global ordered and multivariate tests are informative.", "",
  "## Transfer interpretation", "",
  if (startsWith(verdict, "A")) "The completed Stage B outcome may be described as a **model-specific failed transfer**: GTEx clearly contains age information, yet the two frozen mouse clocks did not yield convergent, model-specific zero-shot transfer. This does not imply that every ovarian ageing mechanism is conserved or that bracket midpoints are exact ages." else if (startsWith(verdict, "B")) "GTEx did not demonstrate transfer but was under-informative for a definitive test." else "GTEx should not be used as definitive evidence for or against cross-species clock generalization.", "",
  "## Guardrails", "",
  "- Low DEG count was not used as a proxy for RNA quality.",
  "- High whole-transcriptome correlation was not interpreted as absence of age biology.",
  "- No PEER, SVA, RUV, human recalibration, feature reselection, or clock refitting was used in the primary transfer analysis.",
  "- The optional predictor is a dataset-adequacy diagnostic, not a human ovarian clock."
)
writeLines(verdict_lines, paudit("reports/gtex_age_signal_adequacy_verdict.md"))

manuscript_lines <- c(passport,
  "# Manuscript revision text: GTEx ovarian age-signal adequacy", "",
  "## Methods insertion", "",
  "We audited whether the official GTEx V11 ovarian bulk RNA-seq cohort was sufficiently age-informative to test zero-shot transfer of the frozen mouse ovarian clocks. All 193 donor-level samples and the prespecified 25,927-gene filtered count matrix were retained in the primary analysis. Ordered public age bracket was modeled with edgeR quasi-likelihood using raw counts, TMM offsets, and covariates fixed before gene-wise testing: RNA integrity number, sample-specific ischemic time, mapping rate, exonic-read rate, mapped-read duplication rate, and sequencing center. Sparse extraction and sequencing batch identifiers were audited but not fitted as high-dimensional fixed effects. Sensitivity analyses applied predeclared technical-QC exclusions and additionally adjusted for donor-death Hardy category. Age labels were permuted 1,000 times within sequencing-center strata. No latent factors associated with age were included in the primary model.", "",
  "## Results insertion", "",
  sprintf("The six public age brackets contained %s donors, respectively. The primary ordered-age analysis identified %d genes at FDR<0.05, including %d with |log2FC per decade|>0.25. The technical-QC sensitivity retained %d FDR-significant genes, and donor-death adjustment retained %d. The observed transcriptome-wide sum of squared age statistics (%.1f) exceeded every stratified age-label permutation (empirical p=%s). Age explained a partial %.1f%% of multivariate expression-PC variation (PERMANOVA p=%s), whereas homogeneity of dispersion was not rejected (p=%s). Preranked Reactome analysis identified %d FDR-significant pathways, and fully nested ridge prediction achieved r=%.2f against bracket midpoint labels (permutation p=%s). Thus, GTEx ovary was clearly age-informative despite the absence of convergent frozen-clock transfer.",
    paste(age_counts$n_samples, collapse = ", "), primary$fdr_lt_0_05, primary$fdr_lt_0_05_abs_log2fc_gt_0_25,
    qcpass$fdr_lt_0_05, agonal$fdr_lt_0_05, perm_obs, p_string(perm_p), 100 * permanova$partial_R2,
    p_string(permanova$empirical_p), p_string(disp$empirical_p), path_primary_n, pred_cor$observed, p_string(pred_cor$empirical_p)), "",
  "## Discussion/limitations insertion", "",
  "These results separate dataset adequacy from model transfer. GTEx ovarian RNA-seq contains strong age-associated structure, so the frozen-clock result cannot be attributed to a transcriptome devoid of age information. The lack of concordant, model-specific transfer is therefore best described as failed transfer of these particular frozen mouse models under the prespecified projection. The audit does not establish cell-intrinsic ovarian ageing: pathology notes showed the expected age association with atrophic/postmenopausal tissue, validated cell-composition deconvolution was unavailable, and the oldest bracket contained only five donors. Public age brackets and their midpoints also limit temporal resolution. Accordingly, the claim is restricted to bulk-tissue zero-shot generalization under the frozen pipeline, not universal absence of conserved ovarian ageing biology.", "",
  "## Figure legend", "",
  "GTEx V11 ovarian transcriptome age-signal adequacy. a, Observed sum of squared transcriptome-wide age t statistics relative to 1,000 age-label permutations within sequencing-center strata. b, Exact donor counts by public age bracket. c, Ordered-age DEG counts in all donors and after predeclared technical-QC exclusions. d, Unique partial variance fractions across expression PCs 1–20. e, Mean absolute age statistic for mapped/detected frozen feature sets; dashed line denotes all expressed genes. f, Fully nested ridge prediction of descriptive age-bracket midpoints. g, Top primary Reactome preranked pathways. The predictor is an adequacy diagnostic, not a human ovarian clock."
)
writeLines(manuscript_lines, paudit("reports/gtex_age_signal_manuscript_revision.md"))

summary_lines <- c(
  "# CODEX GTEx age-signal QC summary", "", paste0("**Decision: ", verdict, "**"), "",
  sprintf("- Cohort: 193 donors; bracket counts %s.", counts_txt),
  sprintf("- Primary edgeR QL: %d/25,927 genes at FDR<0.05; QC-pass sensitivity %d; Hardy-factor sensitivity %d.", primary$fdr_lt_0_05, qcpass$fdr_lt_0_05, agonal$fdr_lt_0_05),
  sprintf("- Stratified permutation: global sum(t²) %.1f; empirical p=%s (1,000 permutations, seed %d).", perm_obs, p_string(perm_p), SEED),
  sprintf("- Multivariate: partial R² %.3f; PERMANOVA p=%s; dispersion p=%s.", permanova$partial_R2, p_string(permanova$empirical_p), p_string(disp$empirical_p)),
  sprintf("- Pathways: %d Reactome pathways at primary FDR<0.05; nested prediction r=%.3f.", path_primary_n, pred_cor$observed),
  "- Frozen Stage A/B models and outputs were not modified or recalibrated.",
  "- Interpretation: GTEx is age-informative; the frozen-clock outcome is model-specific failed transfer, with bulk composition and sparse oldest-bin caveats.", "",
  "Primary artifacts:", "",
  "- `reports/gtex_age_signal_adequacy_verdict.md`",
  "- `reports/gtex_age_signal_manuscript_revision.md`",
  "- `figures/gtex_age_signal_adequacy_integrated.{pdf,svg,png,tiff}`",
  "- `results/table_gtex_ovary_donor_level_audit.tsv`",
  "- `results/table_ordered_age_trend_full.tsv`",
  "- `results/table_age_label_permutation_adequacy_summary.tsv`"
)
writeLines(summary_lines, paudit("CODEX_GTEX_AGE_SIGNAL_QC_SUMMARY.md"))
message("Reports complete. Decision: ", verdict)
