#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=", "", args[grepl("^--file=", args)][1]), mustWork = TRUE)
root <- dirname(script)
p <- function(...) file.path(root, ...)

checks <- list()
add <- function(name, pass, detail) checks[[length(checks) + 1L]] <<- data.table(check = name, status = ifelse(isTRUE(pass), "PASS", "FAIL"), detail = detail)

required <- c(
  "results/table_gtex_ovary_donor_level_audit.tsv",
  "results/table_gtex_ovary_age_bin_counts.tsv",
  "results/table_gtex_ovary_qc_by_age_bin.tsv",
  "results/table_gtex_ovary_metadata_balance.tsv",
  "results/table_ordered_age_trend_full.tsv",
  "results/table_ordered_age_trend_summary.tsv",
  "results/table_age_omnibus_nonlinear_full.tsv",
  "results/table_pairwise_DEG_counts.tsv",
  "results/table_pairwise_DEG_full.tsv",
  "results/table_age_label_permutation_adequacy_summary.tsv",
  "results/table_age_bin_centroid_correlations_distances.tsv",
  "results/table_permanova_age_structure.tsv",
  "results/table_multivariate_dispersion_test.tsv",
  "results/table_reactome_preranked_pathway_enrichment.tsv",
  "results/table_clock_feature_set_age_information.tsv",
  "reports/gtex_age_signal_adequacy_verdict.md",
  "reports/gtex_age_signal_manuscript_revision.md",
  "CODEX_GTEX_AGE_SIGNAL_QC_SUMMARY.md"
)
for (f in required) add(paste0("exists:", f), file.exists(p(f)) && file.info(p(f))$size > 0, as.character(file.info(p(f))$size))

audit <- fread(p("results/table_gtex_ovary_donor_level_audit.tsv"))
counts <- fread(p("results/table_gtex_ovary_age_bin_counts.tsv"))
ordered <- fread(p("results/table_ordered_age_trend_full.tsv"))
os <- fread(p("results/table_ordered_age_trend_summary.tsv"))
perm <- fread(p("results/table_age_label_permutation_adequacy_summary.tsv"))
draw <- fread(p("results/table_age_label_permutation_draws.tsv"))
pair <- fread(p("results/table_pairwise_DEG_full.tsv"))
perma <- fread(p("results/table_permanova_age_structure.tsv"))
vp <- fread(p("results/table_multivariate_variance_partition.tsv"))
pred <- fread(p("results/table_nested_ridge_age_prediction_summary.tsv"))

add("193 unique donors", nrow(audit) == 193L && uniqueN(audit$donor_id) == 193L, sprintf("rows=%d donors=%d", nrow(audit), uniqueN(audit$donor_id)))
add("six exact age bins sum to 193", nrow(counts) == 6L && sum(counts$n_samples) == 193L, paste(counts$n_samples, collapse = ","))
add("ordered table has two complete 25,927-gene analyses", nrow(ordered) == 2L * 25927L && all(ordered[, .N, by = analysis_set]$N == 25927L), as.character(nrow(ordered)))
add("primary count-aware FDR signal present", os[analysis_set == "primary_all_193", fdr_lt_0_05] > 0, as.character(os[analysis_set == "primary_all_193", fdr_lt_0_05]))
add("1,000 permutation draws", nrow(draw) == 1000L, as.character(nrow(draw)))
add("permutation summary finite", all(is.finite(perm$observed)) && all(is.finite(perm$empirical_p)), paste(range(perm$empirical_p), collapse = ".."))
add("pairwise contrasts include both extreme comparisons", all(c("70-79_vs_20-29", "70-79_vs_20-39") %in% pair$contrast), paste(unique(pair$contrast), collapse = ";"))
add("PERMANOVA valid", perma$partial_R2 >= 0 && perma$partial_R2 <= 1 && perma$empirical_p > 0 && perma$empirical_p <= 1, sprintf("R2=%.4f p=%.4g", perma$partial_R2, perma$empirical_p))
add("variance fractions valid", all(vp[!grepl("marginal", component), unique_partial_variance_fraction] >= 0) && sum(vp[!grepl("marginal", component), unique_partial_variance_fraction]) < 1.0001, sprintf("sum=%.6f", sum(vp[!grepl("marginal", component), unique_partial_variance_fraction])))
add("nested prediction metrics finite", all(is.finite(pred[metric %in% c("correlation", "mae_years"), observed])) && all(is.finite(pred[metric %in% c("correlation", "mae_years"), empirical_p])), paste(pred$metric, collapse = ";"))

figure_stems <- c(
  "age_bin_DEG_count_matrix", "age_bin_DEG_upset", "age_bin_logFC_heatmap",
  "expression_pca_qc_panel", "ordered_age_primary_panel", "age_label_permutation_adequacy",
  "age_bin_centroid_distance", "variance_explained_panel", "reactome_preranked_age_pathways",
  "nested_ridge_age_prediction", "gtex_age_signal_adequacy_integrated"
)
for (stem in figure_stems) for (ext in c("pdf", "svg", "png", "tiff")) {
  f <- p("figures", paste0(stem, ".", ext))
  add(paste0("figure:", stem, ".", ext), file.exists(f) && file.info(f)$size > 1000, as.character(file.info(f)$size))
}

res <- rbindlist(checks)
fwrite(res, p("results/validation_gtex_age_signal_adequacy.tsv"), sep = "\t")
if (any(res$status == "FAIL")) {
  print(res[status == "FAIL"])
  stop("Validation failed: ", sum(res$status == "FAIL"), " checks")
}
message("Validation PASS: ", nrow(res), " checks")
