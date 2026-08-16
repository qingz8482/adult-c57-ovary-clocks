#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))

master <- fread(path_module("derived_data/stageA", "mouse_human_gene_effect_master.tsv"))
enrich <- fread(path_module("results/tables", "table_stageA_enrichment.tsv"))
bg <- master[human_testable_background == TRUE]
bg[, expr_bin := pmin(10L, pmax(1L, ceiling(frank(human_mean_expression, ties.method="average")/.N*10))), by=NULL]
bg[, det_bin := pmin(5L, pmax(1L, ceiling(frank(human_detection_fraction, ties.method="average")/.N*5))), by=NULL]
bg[, var_bin := pmin(10L, pmax(1L, ceiling(frank(human_expression_variance, ties.method="average")/.N*10))), by=NULL]
bg[, stratum := paste(expr_bin,det_bin,var_bin,sep="_")]
pools <- split(seq_len(nrow(bg)), bg$stratum)
B <- as.integer(config$analysis$random_sets)
null_rows <- list()
for (sig in enrich$signature) {
  is_el <- grepl("elasticnet",sig); older <- grepl("older",sig)
  d <- if (is_el) master[elasticnet_nonzero==TRUE] else master[pca_versionB_2000==TRUE]
  wc <- if (is_el) "elasticnet_expression_space_weight" else "pca_expression_space_weight"
  d <- d[is.finite(human_age_t) & ((get(wc)>0)==older)]
  idx_bg <- match(d$human_gene_id_stable,bg$human_gene_id_stable)
  strata <- bg$stratum[idx_bg]; weights <- abs(d[[wc]]); weights <- weights/sum(weights)
  null <- numeric(B)
  for (j in seq_along(strata)) {
    pool <- pools[[strata[j]]]; if (is.null(pool)||!length(pool)) pool <- seq_len(nrow(bg))
    null <- null + weights[j] * sample(bg$human_age_t[pool], B, replace=TRUE) * ifelse(older,1,-1)
  }
  obs <- enrich[signature==sig, observed_weighted_signed_t]
  null_rows[[sig]] <- data.table(signature=sig, observed=obs, predicted_direction_pass=obs>0,
                                 null_mean=mean(null), null_sd=sd(null),
                                 empirical_p=if (obs>0) (1+sum(null>=obs))/(B+1) else 1,
                                 null_q95=quantile(null,.95), n_permutations=B)
  write_tsv(data.table(iteration=seq_len(B),null_statistic=null), path_module("results/null_tests",paste0("stageA_",sig,"_matched_null.tsv")))
}
stagea_null <- rbindlist(null_rows); stagea_null[, fdr:=p.adjust(empirical_p,"BH")]
write_tsv(stagea_null, path_module("results/tables","table_random_signature_benchmarks.tsv"))

# Stage B null: preserve the observed mapped genes and absolute-weight distribution,
# but permute weights/signs across expression-matched mapped genes. This is a frozen
# signature null and avoids any age-based fitting.
gtex <- readRDS(path_module("derived_data/gtex","gtex_ovary_prepared.rds"))
objs <- readRDS(path_module("derived_data/stageB","projection_objects.rds"))
assoc <- fread(path_module("results/tables","table_stageB_score_age_association.tsv"))[analysis=="primary_technical_adjusted"]
stageb <- list()
for (family in names(objs)) {
  o <- objs[[family]]; X <- o$x[,o$detected,drop=FALSE]; w <- o$weights[o$detected]
  age <- gtex$metadata$age_midpoint_sensitivity/10
  observed <- assoc[model_family==family,estimate]
  null <- numeric(B); chunk <- 250L
  for (start in seq(1L,B,by=chunk)) {
    ids <- start:min(B,start+chunk-1L); S <- matrix(NA_real_,nrow=nrow(X),ncol=length(ids))
    for (k in seq_along(ids)) {
      wp <- sample(w,length(w),replace=FALSE)
      s <- as.numeric(X %*% wp); S[,k] <- as.numeric(scale(s))
    }
    null[ids] <- apply(S,2,function(s) coef(lm(s~age))[2])
  }
  stageb[[family]] <- data.table(model_family=family,observed=observed,predicted_direction_pass=observed>0,
                                  null_mean=mean(null),null_sd=sd(null),
                                  empirical_p=if (observed>0) (1+sum(null>=observed))/(B+1) else 1,
                                  null_q95=quantile(null,.95),n_permutations=B)
  write_tsv(data.table(iteration=seq_len(B),null_age_trend=null),path_module("results/null_tests",paste0("stageB_",family,"_matched_frozen_signature_null.tsv")))
}
stageb <- rbindlist(stageb); stageb[,fdr:=p.adjust(empirical_p,"BH")]
write_tsv(stageb,path_module("results/null_tests","stageB_random_signature_benchmarks.tsv"))
message("10,000-draw Stage A and Stage B matched-null benchmarks complete.")
