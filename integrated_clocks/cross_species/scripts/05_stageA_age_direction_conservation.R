#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))

mouse <- fread(path_module("derived_data/mouse_models", "mouse_gene_direction_and_weight_master.tsv"))
orth <- fread(path_module("derived_data/orthology", "mouse_human_orthology_primary_one2one.tsv"))
gtex <- readRDS(path_module("derived_data/gtex", "gtex_ovary_prepared.rds"))
hpa <- readRDS(path_module("derived_data/hpa", "hpa_ovary_prepared.rds"))

orth <- orth[!duplicated(mouse_gene_symbol) & !duplicated(human_gene_id_stable)]
master <- merge(mouse, orth, by = "mouse_gene_symbol", all.x = TRUE)
master <- merge(master, gtex$age_effects, by = "human_gene_id_stable", all.x = TRUE)
master <- merge(master, hpa$age_effects, by = "human_gene_id_stable", all.x = TRUE)
master[, human_testable_background := integrated_gene_universe %in% TRUE & !is.na(human_age_t)]
write_tsv(master, path_module("derived_data/stageA", "mouse_human_gene_effect_master.tsv"))

coverage_one <- function(flag, label, weight_col = NULL) {
  d <- master[get(flag) %in% TRUE]
  mapped <- !is.na(d$human_gene_id_stable)
  detected <- mapped & !is.na(d$human_age_t)
  w <- if (is.null(weight_col)) rep(1, nrow(d)) else abs(d[[weight_col]])
  w[!is.finite(w)] <- 0
  data.table(gene_set = label, n_total = nrow(d), n_one2one = sum(mapped), n_human_detected = sum(detected),
             raw_gene_coverage = mean(mapped), human_detected_coverage = mean(detected),
             absolute_model_weight_coverage = if (sum(w) > 0) sum(w[detected]) / sum(w) else NA_real_,
             positive_weight_coverage = if (!is.null(weight_col) && sum(pmax(d[[weight_col]], 0), na.rm=TRUE)>0) sum(pmax(d[[weight_col]][detected],0),na.rm=TRUE)/sum(pmax(d[[weight_col]],0),na.rm=TRUE) else NA_real_,
             negative_weight_coverage = if (!is.null(weight_col) && sum(abs(pmin(d[[weight_col]], 0)),na.rm=TRUE)>0) sum(abs(pmin(d[[weight_col]][detected],0)),na.rm=TRUE)/sum(abs(pmin(d[[weight_col]],0)),na.rm=TRUE) else NA_real_)
}
coverage <- rbindlist(list(
  coverage_one("integrated_gene_universe", "integrated_universe"),
  coverage_one("elasticnet_hvg_2000", "elasticnet_hvg_2000", "elasticnet_expression_space_weight"),
  coverage_one("pca_versionB_2000", "pca_versionB_2000", "pca_expression_space_weight"),
  coverage_one("deployment_shared_artifact", paste0("deployment_shared_artifact_", unique(master$deployment_shared_artifact_n))),
  coverage_one("elasticnet_nonzero", "elasticnet_nonzero", "elasticnet_expression_space_weight")
), fill = TRUE)
write_tsv(coverage, path_module("results/tables", "table_orthology_coverage.tsv"))

direction_one <- function(d, label, x, y, w = NULL, comparison) {
  ok <- is.finite(d[[x]]) & is.finite(d[[y]])
  dd <- d[ok]
  ww <- if (is.null(w)) rep(1, nrow(dd)) else abs(dd[[w]])
  data.table(analysis = label, comparison = comparison, n_genes = nrow(dd),
             spearman_rho = cor(dd[[x]], dd[[y]], method = "spearman"),
             weighted_spearman_rho = weighted_cor(dd[[x]], dd[[y]], ww),
             sign_concordance = mean(sign(dd[[x]]) == sign(dd[[y]])),
             weighted_sign_concordance = weighted.mean(sign(dd[[x]]) == sign(dd[[y]]), ww),
             older_n = sum(dd[[x]] > 0), younger_n = sum(dd[[x]] < 0))
}
dirs <- rbindlist(list(
  direction_one(master[human_testable_background == TRUE], "integrated_testable_background", "mouse_age_beta", "human_age_beta_per_decade", comparison = "observed_expression_age"),
  direction_one(master[elasticnet_hvg_2000 == TRUE], "elasticnet_hvg_2000", "mouse_age_beta", "human_age_beta_per_decade", comparison = "observed_expression_age"),
  direction_one(master[pca_versionB_2000 == TRUE], "pca_versionB_2000", "mouse_age_beta", "human_age_beta_per_decade", comparison = "observed_expression_age"),
  direction_one(master[elasticnet_nonzero == TRUE], "elasticnet_nonzero", "elasticnet_expression_space_weight", "human_age_beta_per_decade", "elasticnet_expression_space_weight", "model_weight_vs_human_age"),
  direction_one(master[pca_versionB_2000 == TRUE], "pca_versionB_2000", "pca_expression_space_weight", "human_age_beta_per_decade", "pca_expression_space_weight", "model_weight_vs_human_age")
), fill = TRUE)
write_tsv(dirs, path_module("results/tables", "table_stageA_direction_conservation.tsv"))

enrich_sets <- list(
  elasticnet_older = master[elasticnet_nonzero == TRUE & elasticnet_expression_space_weight > 0],
  elasticnet_younger = master[elasticnet_nonzero == TRUE & elasticnet_expression_space_weight < 0],
  pca_older = master[pca_versionB_2000 == TRUE & pca_expression_space_weight > 0],
  pca_younger = master[pca_versionB_2000 == TRUE & pca_expression_space_weight < 0]
)
enrich <- rbindlist(lapply(names(enrich_sets), function(nm) {
  d <- enrich_sets[[nm]][is.finite(human_age_t)]
  wc <- if (grepl("elasticnet", nm)) "elasticnet_expression_space_weight" else "pca_expression_space_weight"
  expected <- if (grepl("older", nm)) 1 else -1
  data.table(signature = nm, n_genes = nrow(d), expected_direction = expected,
             observed_signed_mean_t = mean(expected * d$human_age_t),
             observed_weighted_signed_t = weighted.mean(expected * d$human_age_t, abs(d[[wc]])),
             hpa_exploratory_signed_rho = mean(expected * d$hpa_age_spearman, na.rm = TRUE), dataset = "GTEx_v11_primary;HPA_v25.1_exploratory")
}))
write_tsv(enrich, path_module("results/tables", "table_stageA_enrichment.tsv"))
message("Stage A observed conservation statistics complete.")
