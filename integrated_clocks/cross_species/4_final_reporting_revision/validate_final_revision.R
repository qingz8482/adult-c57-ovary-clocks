#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "scripts/common_final.R"))

checks <- list()
add <- function(id, pass, detail) checks[[length(checks)+1L]] <<- data.table(check_id=id,status=ifelse(isTRUE(pass),"PASS","FAIL"),detail=detail)
md_table <- function(d) {
  d <- as.data.table(d); d[] <- lapply(d,function(z)gsub("\\|","\\\\|",as.character(z)))
  c(paste0("| ",paste(names(d),collapse=" | ")," |"),paste0("| ",paste(rep("---",ncol(d)),collapse=" | ")," |"),
    apply(d,1,function(z)paste0("| ",paste(z,collapse=" | ")," |")))
}
all_reports <- c("final_gtex_age_adequacy_results.md","final_gtex_age_adequacy_methods.md","final_gtex_age_adequacy_discussion.md",
  "final_gtex_age_adequacy_figure_legends.md","integration_with_frozen_clock_failure.md")
report_paths <- pfinal("reports",all_reports)
txt <- if(all(file.exists(report_paths))) unlist(lapply(report_paths,readLines,warn=FALSE)) else character()
add("relocated_path_layout",basename(final_dir)=="4_final_reporting_revision" && basename(base_dir)=="3_gtex_age_signal_adequacy" && identical(module_dir,dirname(final_dir)),
  paste("module=",module_dir,"; adequacy=",base_dir,"; final=",final_dir))

## Required outputs and figure variants.
add("reports_present",all(file.exists(report_paths)),paste(all_reports[file.exists(report_paths)],collapse=", "))
bases <- c("01_centroid_correlation_heatmaps_complete","02_centroid_heatmap_raw_pearson","03_centroid_heatmap_raw_spearman",
  "04_centroid_heatmap_residualized_pearson","05_centroid_heatmap_residualized_spearman","06_centroid_heatmap_residualized_one_minus_pearson",
  "07_sample_similarity_dependence_aware","08_hardy_gene_coefficient_stability","09_hardy_pathway_stability","13_representative_age_pattern_heatmap",
  "12_sparse_oldest_sensitivity_metrics","11_leave_one_oldest_donor_clock_slopes","10_donor_cooks_distance",
  "14_nested_cv_validation_predictions","15_nested_cv_validation_metrics","16_nested_cv_stratified_permutation")
formats <- c("pdf","svg","png","tiff")
fp <- unlist(lapply(bases,function(b)pfinal("figures",paste0(b,".",formats))))
add("figure_formats_complete",all(file.exists(fp)),sprintf("%d/%d required format files present",sum(file.exists(fp)),length(fp)))
src_need <- c("centroid_correlation_heatmaps.tsv","sample_similarity_dependence_aware.tsv","hardy_gene_coefficient_stability.tsv",
  "hardy_pathway_stability.tsv","representative_age_pattern_heatmap.tsv","sparse_oldest_sensitivity_metrics.tsv",
  "nested_cv_validation_predictions.tsv","nested_cv_validation_metrics.tsv")
add("panel_source_data_present",all(file.exists(pfinal("figures/source_data",src_need))),paste(src_need,collapse=", "))

## Guardrail: pair rows are not inferential units.
sim <- fread(pfinal("results/table_sample_similarity_dependence_aware_contrasts.tsv"))
pair_ok <- all(sim$n_permutations==1000L) && all(sim$bootstrap_replicates==1000L) &&
  all(grepl("donor",sim$inference_unit,ignore.case=TRUE)) && !any(grepl("independent pair|pairs? as independent",txt,ignore.case=TRUE))
add("no_independent_pair_inference",pair_ok,"Donor bootstrap and within-center donor-label permutations; pair counts descriptive only")
add("permutation_p_has_count",all(is.finite(sim$empirical_permutation_p)) && all(sim$n_permutations>0L) &&
  any(grepl("empirical p.*1,000 permutations|1,000.*empirical p",txt,ignore.case=TRUE)),"Similarity p values carry n_permutations=1000 and reports state the count")

## Guardrail: effect size language and tables.
eff <- fread(pfinal("results/table_effect_size_thresholds_by_specification.tsv"))
need_crit <- c("fdr_lt_0_05","fdr_and_abs_log2fc_gt_0_10","fdr_and_abs_log2fc_gt_0_25","fdr_and_abs_log2fc_gt_0_50",
  "fdr_and_abs_standardized_beta_gt_0_10","fdr_and_abs_standardized_beta_gt_0_25")
bad_large <- grepl("all FDR[- ]significant genes (are|were) (biologically )?large|FDR[- ]significant genes all had large",txt,ignore.case=TRUE)
add("effect_size_thresholds_and_language",all(need_crit%in%eff$criterion)&&all(c("positive","negative")%in%names(eff))&&!any(bad_large),
  "All requested thresholds and directions present; no wholesale large-effect claim")

## Guardrail: sparse-bin sensitivity.
sparse <- fread(pfinal("results/table_sparse_oldest_sensitivity_summary.tsv"))
required_sparse <- c("primary_all_bins","exclude_70_79","merge_60_79",paste0("leave_bin_",gsub("-","_",age_levels)))
old_n <- sum(sparse$sensitivity_type=="leave_one_70_79_donor")
add("sparse_70_79_sensitivity",all(required_sparse%in%sparse$sensitivity)&&old_n==5L&&file.exists(pfinal("results/table_donor_cooks_distance.tsv")),
  sprintf("%d specifications; %d leave-one-oldest-donor fits; Cook's-distance table present",nrow(sparse),old_n))

## Guardrail: all fold-dependent preprocessing lives inside training folds.
leak <- fread(pfinal("results/table_nested_cv_leakage_audit.tsv"))
leak_ok <- all(leak$leakage_status%in%c("PASS","NOT_APPLICABLE","fixed upstream scope")) && all(vapply(c("count filtering","feature selection","scaling","tuning"),
  function(z)any(grepl(z,leak$step,ignore.case=TRUE)),logical(1)))
add("nested_cv_no_leakage",leak_ok,paste(leak$step,leak$leakage_status,sep="=",collapse="; "))
cv <- fread(pfinal("results/table_nested_cv_performance_validated.tsv"))
metrics <- c("pearson","spearman","mae_midpoint_years","balanced_age_bracket_accuracy","ordinal_concordance_kendall","calibration_slope")
add("nested_cv_metrics_sensitivities",all(metrics%in%names(cv))&&all(c("baseline_age_stratified_5fold","sequencing_center_blocked","batch_stratified_5fold","hardy_age_stratified_5fold")%in%cv$specification),
  paste(cv$specification,collapse=", "))
add("nested_cv_permutation",file.exists(pfinal("results/table_nested_cv_stratified_permutation_null.tsv"))&&all(fread(pfinal("results/table_nested_cv_permutation_summary.tsv"))$n_permutations==100L),
  "100 sequencing-center-stratified full-workflow label permutations")

## Guardrail: adequacy predictor is not represented as a validated human clock.
clock_lines <- txt[grepl("human ovarian clock",txt,ignore.case=TRUE)]
clock_ok <- length(clock_lines)>0L && all(grepl("not|no ",clock_lines,ignore.case=TRUE))
add("adequacy_not_human_clock",clock_ok,paste(clock_lines,collapse=" | "))

## Guardrail: no hidden heatmap truncation.
sc <- fread(pfinal("results/table_centroid_heatmap_scale_contract.tsv"))
heat_ok <- nrow(sc)==5L && all(sc$scale_explicitly_printed) && all(sc$cell_values_annotated) &&
  all(is.finite(sc$display_min)) && all(is.finite(sc$display_max)) && any(grepl("one_minus",sc$metric))
add("heatmap_scale_explicit",heat_ok,"Five scales tabulated; every cell annotated; separate 1-r panel")

## Guardrail: Hardy stability quantified.
hardy <- fread(pfinal("results/table_hardy_gene_coefficient_stability.tsv")); hp <- fread(pfinal("results/table_hardy_pathway_stability.tsv"))
hardy_ok <- "primary vs hardy_adjusted"%in%hardy$comparison && "primary vs hardy_adjusted"%in%hp$comparison &&
  all(c("pearson_coefficient","spearman_coefficient","sign_concordance","proportion_primary_fdr_retaining_fdr_lt_0_05")%in%names(hardy)) &&
  all(c("pearson_NES","spearman_NES","median_leading_edge_jaccard")%in%names(hp))
add("hardy_stability_quantified",hardy_ok,"Gene coefficients, directions, retained FDR, Reactome NES, and leading-edge overlap quantified")

## Guardrail: bounded biological interpretation and exact required statements.
unconserved_lines <- txt[grepl("unconserved",txt,ignore.case=TRUE)]
unconserved_ok <- length(unconserved_lines)>0L && all(grepl("cannot|does not|not warrant|should not",unconserved_lines,ignore.case=TRUE))
add("no_unconserved_biology_claim",unconserved_ok,paste(unconserved_lines,collapse=" | "))
req1 <- "GTEx V11 ovarian tissue contains substantial age-associated transcriptomic structure, so the absence of convergent and model-specific frozen-clock transfer cannot be attributed simply to an age-uninformative target cohort."
req2 <- "Age explained a measurable but not dominant fraction of whole-transcriptome variation."
add("required_interpretation_present",any(grepl(req1,txt,fixed=TRUE))&&any(grepl(req2,txt,fixed=TRUE)),"Both required statements present verbatim")
add("no_definitive_representation_claim",!any(grepl("GTEx provides a definitive representation of human ovarian biological ageing",txt,fixed=TRUE)),"Prohibited affirmative claim absent")

out <- rbindlist(checks)
write_tsv(out,pfinal("results/validation_final_reporting_revision.tsv"))

commands <- c(
  "Rscript scripts/01_centroid_similarity.R",
  "Rscript scripts/02_hardy_features_patterns.R",
  "Rscript scripts/03_sparse_oldest_sensitivity.R",
  "Rscript scripts/04_nested_cv_validation.R",
  "Rscript scripts/05_reports_figures.R",
  "Rscript scripts/06_visual_qa_repairs.R",
  "Rscript validate_final_revision.R")
md <- c("# Final reporting revision validation","",sprintf("Overall: **%s**",ifelse(all(out$status=="PASS"),"PASS","FAIL")),"",
  "## Exact execution commands","","```bash",commands,"```","","## Output paths","",
  "- Figures: `4_final_reporting_revision/figures/`","- Panel source data: `4_final_reporting_revision/figures/source_data/`",
  "- Numerical results: `4_final_reporting_revision/results/`","- Manuscript reports: `4_final_reporting_revision/reports/`","",
  "## Checks","",md_table(out))
writeLines(md,pfinal("reports/validation_final_reporting_revision.md"),useBytes=TRUE)
message(sprintf("Validation: %d PASS, %d FAIL",sum(out$status=="PASS"),sum(out$status=="FAIL")))
if(any(out$status=="FAIL")) stop("Final reporting revision validation failed")
