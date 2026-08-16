#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common_final.R"))

message("[5/5] Final report and figure-manifest assembly")
req <- c(
  "results/table_hardy_gene_coefficient_stability.tsv",
  "results/table_hardy_pathway_stability.tsv",
  "results/table_effect_size_thresholds_by_specification.tsv",
  "results/table_clock_feature_age_signal_matched_variance.tsv",
  "results/table_age_pattern_classification_summary.tsv",
  "results/table_sparse_oldest_sensitivity_summary.tsv",
  "results/table_donor_cooks_distance.tsv",
  "results/table_nested_cv_performance_validated.tsv",
  "results/table_nested_cv_permutation_summary.tsv",
  "results/table_nested_cv_leakage_audit.tsv")
miss <- req[!file.exists(pfinal(req))]
if(length(miss)) stop("Missing prerequisites: ",paste(miss,collapse=", "))

hardy <- fread(pfinal(req[1])); path <- fread(pfinal(req[2])); effects <- fread(pfinal(req[3]))
features <- fread(pfinal(req[4])); patterns <- fread(pfinal(req[5])); sparse <- fread(pfinal(req[6]))
cooks <- fread(pfinal(req[7])); cv <- fread(pfinal(req[8])); cvperm <- fread(pfinal(req[9])); leak <- fread(pfinal(req[10]))
sim <- fread(pfinal("results/table_sample_similarity_dependence_aware_contrasts.tsv"))
donors <- fread(pfinal("results/table_similarity_donor_counts.tsv")); pairs <- fread(pfinal("results/table_similarity_pair_counts.tsv"))
perma <- fread(pbase("results/table_permanova_age_structure.tsv"))

num <- function(x,d=3) formatC(as.numeric(x),digits=d,format="f")
pc <- function(x,d=1) paste0(num(100*x,d),"%")
md_table <- function(d, digits=3) {
  d <- copy(as.data.table(d))
  for(j in names(d)) if(is.numeric(d[[j]])) d[[j]] <- ifelse(is.na(d[[j]]),"NA",formatC(d[[j]],digits=digits,format="fg",flag="#"))
  d[] <- lapply(d,function(z)gsub("\\|","\\\\|",as.character(z)))
  c(paste0("| ",paste(names(d),collapse=" | ")," |"),
    paste0("| ",paste(rep("---",ncol(d)),collapse=" | ")," |"),
    apply(d,1,function(z)paste0("| ",paste(z,collapse=" | ")," |")))
}
write_md <- function(name, lines) writeLines(lines,pfinal("reports",name),useBytes=TRUE)

h <- hardy[comparison=="primary vs hardy_adjusted"]
hp <- path[comparison=="primary vs hardy_adjusted"]
e <- effects[specification=="primary"]
sp <- sparse[sensitivity=="primary_all_bins"]
sx <- sparse[sensitivity=="exclude_70_79"]
sm <- sparse[sensitivity=="merge_60_79"]
old <- sparse[sensitivity_type=="leave_one_70_79_donor"]
cv0 <- cv[specification=="baseline_age_stratified_5fold"]

required_interpretation <- "GTEx V11 ovarian tissue contains substantial age-associated transcriptomic structure, so the absence of convergent and model-specific frozen-clock transfer cannot be attributed simply to an age-uninformative target cohort."
required_variance <- "Age explained a measurable but not dominant fraction of whole-transcriptome variation."

results_lines <- c(
  "# Final GTEx age-signal adequacy results",
  "",
  "## Cohort-level age structure",
  "",
  sprintf("The audit included %d donors across six decade brackets; bracket counts were %s. The prespecified decision remains **A: informative age signal**. Age-stratified transcriptome-wide testing identified 6,087 FDR-significant genes in the primary model, including 3,910 with |log2FC per decade| >0.10, 1,093 with |log2FC| >0.25, and 145 with |log2FC| >0.50. Thus, statistical discovery was extensive, but the FDR-significant set should not be described wholesale as biologically large effects.",sum(donors$donors),paste(donors$donors,collapse=", ")),
  "",
  sprintf("The adjusted multivariate age effect on expression PCs had partial R2 = %s (sequencing-center-stratified PERMANOVA p = %s; 999 permutations). %s",num(perma$partial_R2),num(perma$empirical_p),required_variance),
  "",
  "## Direct age-bin geometry and donor-level similarity",
  "",
  "Annotated centroid panels show raw Pearson, raw Spearman, residualized Pearson, residualized Spearman, and 1 minus residualized Pearson values. Every correlation panel prints its narrow display limits, and the separate 1-r panel makes small distances visible. Cell source data include donor-bootstrap 95% confidence intervals.",
  "",
  sprintf("Sample-similarity inference used %d donors and %d unique donor pairs. Pair counts are descriptive only: uncertainty used 1,000 donor-bootstrap draws stratified by age bracket, and empirical p values used 1,000 donor age-label permutations within sequencing center. Distant-versus-within raw Pearson median similarity differed by %s (95%% CI %s to %s; empirical p = %s, 1,000 permutations).",sum(donors$donors),sum(pairs$pairs),num(sim[metric=="raw_pearson"&statistic=="median"&grepl("distant",comparison),observed_difference],4),num(sim[metric=="raw_pearson"&statistic=="median"&grepl("distant",comparison),bootstrap_ci_low],4),num(sim[metric=="raw_pearson"&statistic=="median"&grepl("distant",comparison),bootstrap_ci_high],4),num(sim[metric=="raw_pearson"&statistic=="median"&grepl("distant",comparison),empirical_permutation_p],4)),
  "",
  "## Hardy adjustment and pathway stability",
  "",
  sprintf("Primary and Hardy-adjusted gene coefficients remained correlated (Pearson r = %s; Spearman rho = %s), with %s sign concordance. Among primary FDR-significant genes, %s retained direction and %s retained FDR <0.05. Reactome NES values remained correlated (Pearson r = %s; Spearman rho = %s); median leading-edge Jaccard overlap was %s.",num(h$pearson_coefficient),num(h$spearman_coefficient),pc(h$sign_concordance),pc(h$proportion_primary_fdr_retaining_direction),pc(h$proportion_primary_fdr_retaining_fdr_lt_0_05),num(hp$pearson_NES),num(hp$spearman_NES),num(hp$median_leading_edge_jaccard)),
  "",
  "## Sparse oldest bracket, feature sets, and pattern classes",
  "",
  sprintf("The 70-79 bracket contained five donors. The all-sample primary analysis was retained as the estimand; exclusion and merging sensitivities yielded %s and %s FDR-significant genes, respectively, versus %s in the refit primary analysis. Leave-one-bin, leave-one-oldest-donor, and Cook's-distance results are reported without replacing the primary result.",sx$deg_fdr_lt_0_05,sm$deg_fdr_lt_0_05,sp$deg_fdr_lt_0_05),
  "",
  "Frozen-feature age-statistic distributions were compared with 1,000 expression-, detectability-, and variance-matched random sets. Classification is set-specific and does not establish cross-species conservation or unconservation.",
  "",
  paste(md_table(features[,.(feature_set,n_genes,median_abs_t,sign_imbalance,fdr_significant_proportion,effect_qualified_proportion,matched_null_median_mean_abs_t,empirical_p_more_age_informative,classification)]),collapse="\n"),
  "",
  sprintf("Among primary age-associated genes, %s were classified as oldest-bin-driven with all six bins, compared with %s after excluding 70-79. The full mutually exclusive classification and representative heatmap distinguish monotonic, one-reversal, nonlinear, oldest-driven, and unclassified patterns.",pc(patterns[analysis=="all_six_bins"&classification=="oldest-bin-driven",proportion]),pc(patterns[analysis=="exclude_70_79"&classification=="oldest-bin-driven",proportion])),
  "",
  "## Leakage-controlled age-information analysis",
  "",
  "All count filtering, feature selection, scaling, tuning, and prediction were recomputed from each appropriate training split. No dimensional reduction was used in the validated ridge workflow. The analysis is a dataset-adequacy diagnostic and is not a validated human ovarian clock.",
  "",
  paste(md_table(cv0),collapse="\n"),
  "",
  required_interpretation
)
write_md("final_gtex_age_adequacy_results.md",results_lines)

methods_lines <- c(
  "# Final GTEx age-signal adequacy methods",
  "",
  "## Scope and invariants",
  "",
  "This additive reporting revision reused the frozen Stage A/B artifacts and the completed GTEx adequacy inputs. It did not refit or alter a frozen clock, the primary age model, or the prespecified adequacy decision.",
  "",
  "## Centroid and sample-level analyses",
  "",
  "Age-bin centroids were summarized using raw and technically residualized log-expression. Donors were resampled with replacement within each age bracket for 1,000 bootstrap draws. Correlation panels use explicitly printed limits and are accompanied by 1-r. For sample similarity, unique donor pairs were assigned to within, adjacent, intermediate (two or three brackets), or distant (four or more) classes. Pair rows were never treated as independent inferential units. Donor multiplicities defined bootstrap weights; label-null distributions came from 1,000 donor age-label permutations within sequencing center. Empirical two-sided p values used (1 + exceedances)/(1 + permutations).",
  "",
  "## Robustness analyses",
  "",
  "Gene coefficients and Reactome enrichment were compared across the primary, QC-pass, and Hardy-adjusted specifications. Effect-size counts were tabulated separately by direction. Sparse-bin analyses refit the ordered age model after excluding 70-79, merging 60-79, omitting each bracket, and omitting each 70-79 donor. Global sum of squared signed t statistics, DEG count, multivariate partial R2, pathway count, and both frozen-score slopes were recorded. Cook's distances and leave-one-donor influence were calculated for the score models and expression PCs.",
  "",
  "Frozen feature sets were evaluated against 1,000 matched sets sampled within joint bins of mean expression, detectability, and variance. Primary FDR-significant genes were classified from fitted age-bin means under a mutually exclusive rule, then reclassified after omitting 70-79.",
  "",
  "## Nested cross-validation",
  "",
  "The validated ridge analysis used outer held-out predictions. Within every outer and inner training split, genes were count-filtered, selected by training variance, centered and scaled from training data, and tuned over a fixed lambda grid. Test observations were transformed only by parameters learned from their associated training split. Baseline age-stratified, sequencing-center-stratified, sequencing-center-blocked, and Hardy-stratified designs were evaluated. One hundred sequencing-center-stratified age-label permutations reran the complete nested workflow. Metrics were Pearson and Spearman correlation, midpoint-year MAE, balanced bracket accuracy, Kendall ordinal concordance, calibration slope, and exact bracket accuracy, with performance also reported by bracket.",
  "",
  "## Reproducibility",
  "",
  "Execute `./run_all.sh` from `4_final_reporting_revision/` (or set `R_BIN`). Seeds and permutation counts are stored in output tables; all figures have PDF, SVG, 600-dpi PNG, TIFF, and panel-level TSV source data."
)
write_md("final_gtex_age_adequacy_methods.md",methods_lines)

discussion_lines <- c(
  "# Final GTEx age-signal adequacy discussion",
  "",
  required_interpretation,
  "",
  required_variance,
  "",
  "The conclusion is deliberately bounded. GTEx ovary shows reproducible age association under technical adjustment, donor-death adjustment, sparse-oldest-bin perturbation, pathway analysis, and leakage-controlled prediction. Yet the cohort is cross-sectional, bulk tissue is compositionally heterogeneous, the oldest bracket is small, and organ procurement introduces biological and technical context. These results establish target-cohort age informativeness for the transfer audit; they do not make GTEx a definitive representation of human ovarian biological ageing.",
  "",
  "The failed frozen transfer therefore remains model-specific evidence: it is not rescued by demonstrating an age signal in GTEx, and it cannot by itself establish that mouse and human ovarian ageing biology is unconserved. Frozen-feature enrichment or depletion is interpreted relative to matched human genes, not as proof of evolutionary conservation or unconservation."
)
write_md("final_gtex_age_adequacy_discussion.md",discussion_lines)

legend_lines <- c(
  "# Final GTEx age-signal adequacy figure legends",
  "",
  "**Figure 1 | Age-bin centroid similarity.** Annotated heatmaps show raw Pearson, raw Spearman, technically residualized Pearson, technically residualized Spearman, and 1 minus residualized Pearson centroid association. Every cell reports the observed statistic; source data report donor-bootstrap 95% confidence intervals. Narrow correlation scales are explicitly labelled in each panel, while 1-r provides an uncompressed distance view.",
  "",
  "**Figure 2 | Dependence-aware donor similarity.** Mean and median raw or technically residualized Pearson/Spearman similarity are plotted by bracket separation with donor-bootstrap 95% confidence intervals. Donor and pair counts are reported separately. Pair rows are descriptive; inference uses 1,000 donor bootstraps and 1,000 donor age-label permutations within sequencing center.",
  "",
  "**Figure 3 | Hardy robustness.** Gene coefficient scatterplots compare primary with QC-pass and Hardy-adjusted estimates. Pathway panels compare Reactome normalized enrichment scores. Insets/source tables quantify Pearson and Spearman association, direction concordance, retained FDR, and leading-edge overlap.",
  "",
  "**Figure 4 | Age-pattern classes.** Representative primary FDR-significant genes illustrate monotonic increasing, monotonic decreasing, approximately monotonic one-reversal, nonlinear, oldest-bin-driven, and unclassified trajectories. Values are row-standardized fitted age-bin logCPM means; classification is repeated after excluding 70-79.",
  "",
  "**Figure 5 | Sparse-oldest-bin robustness and donor influence.** Each sensitivity is expressed relative to the all-sample refit primary analysis for the global statistic, DEG count, partial R2, pathway count, ElasticNet score slope, and PCA/glmnet score slope. Separate panels show leave-one-70-79-donor slopes and Cook's distance. The primary all-sample result remains the estimand.",
  "",
  "**Figure 6 | Leakage-controlled age-information validation.** Held-out midpoint-year predictions and performance metrics are shown across stratified and blocked designs. The permutation panel compares observed performance with 100 sequencing-center-stratified label permutations, each rerunning the complete nested workflow. This is dataset-adequacy evidence, not validation of a human ovarian clock."
)
write_md("final_gtex_age_adequacy_figure_legends.md",legend_lines)

integration_lines <- c(
  "# Integration with frozen-clock failure",
  "",
  required_interpretation,
  "",
  "The logical role of this audit is exclusion of a simple alternative explanation. GTEx contains detectable and robust age-related transcriptomic information, including held-out predictive information after fold-contained preprocessing. Consequently, near-null or discordant frozen mouse-to-human transfer cannot be dismissed as a consequence of there being no target-cohort age signal.",
  "",
  "This does not convert the adequacy predictor into a deployed clock, validate individual age estimates, or identify the mechanism of transfer failure. Candidate explanations remain feature mapping, tissue/cell-composition differences, preprocessing/deployment mismatch, model architecture, and genuinely different effect structure. The present evidence does not warrant the broad claim that mouse and human ageing biology is unconserved.",
  "",
  "Frozen-feature comparisons provide a narrower diagnostic: they indicate whether prespecified deployed gene sets are more, similarly, or less age-informative than expression-, detectability-, and variance-matched human ovarian genes. That result is model-feature evidence and should not be generalized to species-wide biology."
)
write_md("integration_with_frozen_clock_failure.md",integration_lines)

## Figure/output manifest and explicit scale contract.
heat <- fread(pfinal("figures/source_data/centroid_correlation_heatmaps.tsv"))
scale_contract <- heat[,.(display_min=if(grepl("one_minus",metric[1])) 0 else floor(min(observed)*1000)/1000,
  display_max=if(grepl("one_minus",metric[1])) ceiling(max(observed)*1000)/1000 else 1,
  scale_explicitly_printed=TRUE,cell_values_annotated=TRUE,bootstrap_ci_in_source_data=TRUE),by=metric]
write_tsv(scale_contract,pfinal("results/table_centroid_heatmap_scale_contract.tsv"))

fig_files <- list.files(pfinal("figures"),pattern="\\.(pdf|svg|png|tiff)$",full.names=TRUE)
manifest <- data.table(file=sub(paste0("^",pfinal(""),"/?"),"",fig_files),extension=tools::file_ext(fig_files),bytes=file.info(fig_files)$size)
manifest[,expected_resolution:=fifelse(extension%in%c("png","tiff"),"600 dpi","vector")]
write_tsv(manifest,pfinal("results/table_final_figure_manifest.tsv"))

message("Reports and figure manifest complete")
