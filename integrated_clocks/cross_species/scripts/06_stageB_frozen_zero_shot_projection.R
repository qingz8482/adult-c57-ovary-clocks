#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))
suppressPackageStartupMessages({ library(glmnet) })

arts <- readRDS(path_module("derived_data/mouse_models", "authoritative_artifacts.rds"))
mouse_ref <- readRDS(path_module("derived_data/mouse_models", "mouse_deployment_reference_expression.rds"))
orth <- fread(path_module("derived_data/orthology", "mouse_human_orthology_primary_one2one.tsv"))
orth <- orth[!duplicated(mouse_gene_symbol) & !duplicated(human_gene_id_stable)]
gtex <- readRDS(path_module("derived_data/gtex", "gtex_ovary_prepared.rds"))
meta <- copy(gtex$metadata)

all_proj <- list(); projection_objects <- list(); pc_out <- list()
for (family in names(arts)) {
  a <- arts[[family]]; genes <- a$selected_genes
  map <- orth[match(genes, mouse_gene_symbol)]
  detected <- !is.na(map$human_gene_id_stable) & map$human_gene_id_stable %in% rownames(gtex$logcpm)
  raw <- matrix(rep(a$gene_center[genes], each = nrow(meta)), nrow = nrow(meta), dimnames = list(meta$sample_id, genes))
  raw[, detected] <- t(gtex$logcpm[map$human_gene_id_stable[detected], meta$sample_id, drop = FALSE])
  sc <- a$gene_scale[genes]; sc[!is.finite(sc) | sc == 0] <- 1
  x <- sweep(sweep(raw, 2, a$gene_center[genes], "-"), 2, sc, "/")
  x[!is.finite(x)] <- 0
  pred_raw <- predict_artifact(a, x, calibrate = FALSE)
  pred_week <- a$calibration$intercept + a$calibration$slope * pred_raw
  score <- as.numeric(scale(pred_raw))
  w <- model_gene_weights(a); aw <- abs(w)
  weight_cov <- sum(aw[detected]) / sum(aw)

  ref_x <- t(mouse_ref$expression[genes, a$training_sample_keys, drop = FALSE])
  ref_x <- sweep(sweep(ref_x, 2, a$gene_center[genes], "-"), 2, sc, "/")
  frac3 <- rowMeans(abs(x[, detected, drop = FALSE]) > 3)
  frac5 <- rowMeans(abs(x[, detected, drop = FALSE]) > 5)
  outside <- rowMeans(sweep(x[, detected, drop=FALSE], 2, apply(ref_x[, detected, drop=FALSE],2,min), "<") |
                        sweep(x[, detected, drop=FALSE], 2, apply(ref_x[, detected, drop=FALSE],2,max), ">"))
  if (family == "pca_glmnet") {
    ref_low <- predict(a$model$pca_fit, newdata = ref_x)[, a$model$selected_pc_names, drop=FALSE]
    hum_low <- predict(a$model$pca_fit, newdata = x)[, a$model$selected_pc_names, drop=FALSE]
    pcb <- as.numeric(coef(a$model$glmnet_fit))[-1]; names(pcb) <- rownames(coef(a$model$glmnet_fit))[-1]
    pc_contrib <- sweep(hum_low, 2, pcb[a$model$selected_pc_names], "*")
    pc_out[[family]] <- cbind(data.table(sample_id=meta$sample_id), as.data.table(hum_low), as.data.table(pc_contrib, keep.rownames=FALSE))
  } else {
    pfit <- prcomp(ref_x, center = TRUE, scale. = FALSE, rank. = min(10, nrow(ref_x)-1))
    ref_low <- predict(pfit, ref_x); hum_low <- predict(pfit, x)
  }
  cen <- colMeans(ref_low); cv <- cov(ref_low); diag(cv) <- diag(cv) + 1e-6
  ood <- mahalanobis(hum_low, cen, cv)
  nearest <- vapply(seq_len(nrow(hum_low)), function(i) min(sqrt(rowSums((ref_low - matrix(hum_low[i,], nrow(ref_low), ncol(ref_low), byrow=TRUE))^2))), numeric(1))
  qc <- if (weight_cov >= config$analysis$high_weighted_coverage) "PASS_HIGH_COVERAGE_OOD_QUANTIFIED" else if (weight_cov >= config$analysis$moderate_weighted_coverage) "REVIEW_MODERATE_WEIGHTED_COVERAGE" else "FAIL_LOW_WEIGHTED_COVERAGE"
  all_proj[[family]] <- data.table(donor_id=meta$donor_id, sample_id=meta$sample_id, age_bracket=meta$age_bracket,
    age_midpoint_sensitivity=meta$age_midpoint_sensitivity, model_family=family, mapped_gene_count=sum(detected),
    raw_gene_coverage=mean(detected), absolute_model_weight_coverage=weight_cov,
    positive_weight_coverage=sum(pmax(w[detected],0))/sum(pmax(w,0)), negative_weight_coverage=sum(abs(pmin(w[detected],0)))/sum(abs(pmin(w,0))),
    raw_linear_predictor=pred_raw, mouse_week_equivalent_output_audit_only=pred_week,
    cross_species_mouse_ovarian_clock_score=score, cross_species_score_percentile=frank(score)/length(score),
    projection_ood_distance=ood, fraction_features_abs_z_gt3=frac3, fraction_features_abs_z_gt5=frac5,
    fraction_mapped_features_outside_mouse_range=outside, nearest_mouse_reference_distance=nearest, projection_qc_class=qc)
  projection_objects[[family]] <- list(x=x, detected=detected, map=map, weights=w, reference_x=ref_x, human_low=hum_low, reference_low=ref_low)
}
proj <- rbindlist(all_proj)
write_tsv(proj, path_module("derived_data/stageB", "gtex_frozen_zero_shot_scores.tsv"))
saveRDS(projection_objects, path_module("derived_data/stageB", "projection_objects.rds"))
if (length(pc_out)) write_tsv(rbindlist(pc_out, idcol="model_family", fill=TRUE), path_module("derived_data/stageB", "gtex_frozen_pca_scores_and_contributions.tsv"))

assoc <- list(); ood_assoc <- list()
for (family in unique(proj$model_family)) {
  d <- merge(proj[model_family == family], meta, by=c("donor_id","sample_id","age_bracket","age_midpoint_sensitivity"))
  d[, age_decade := age_midpoint_sensitivity/10]
  d[, rin_z := as.numeric(scale(rin))]
  f <- cross_species_mouse_ovarian_clock_score ~ age_decade + rin_z
  a <- hc3_lm(f, d)[term == "age_decade"]
  a[, `:=`(model_family=family, analysis="primary_technical_adjusted", spearman=cor(d$cross_species_mouse_ovarian_clock_score,d$age_midpoint_sensitivity,method="spearman"))]
  b <- hc3_lm(cross_species_mouse_ovarian_clock_score ~ age_decade, d)[term == "age_decade"]
  b[, `:=`(model_family=family, analysis="sensitivity_age_only", spearman=cor(d$cross_species_mouse_ovarian_clock_score,d$age_midpoint_sensitivity,method="spearman"))]
  assoc[[family]] <- rbind(a,b,fill=TRUE)
  oa <- hc3_lm(projection_ood_distance ~ age_decade + rin_z, d)[term == "age_decade"]
  oa[, model_family:=family]; ood_assoc[[family]] <- oa
}
assoc <- rbindlist(assoc); assoc[analysis=="primary_technical_adjusted", fdr:=p.adjust(p,"BH")]
write_tsv(assoc, path_module("results/tables", "table_stageB_score_age_association.tsv"))
write_tsv(rbindlist(ood_assoc), path_module("results/tables", "table_stageB_ood_qc.tsv"))
write_tsv(proj[, .(n=.N, score_mean=mean(cross_species_mouse_ovarian_clock_score), score_se=sd(cross_species_mouse_ovarian_clock_score)/sqrt(.N)), by=.(model_family,age_bracket)], path_module("results/stageB", "gtex_score_by_age_bracket.tsv"))

wide <- dcast(proj, donor_id + age_bracket ~ model_family, value.var="cross_species_mouse_ovarian_clock_score")
cross <- data.table(n=nrow(wide), pearson=cor(wide$elasticnet,wide$pca_glmnet), spearman=cor(wide$elasticnet,wide$pca_glmnet,method="spearman"))
write_tsv(cross, path_module("results/stageB", "cross_model_concordance.tsv"))
message("Stage B frozen projection and OOD analyses complete.")
