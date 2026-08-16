#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))
suppressPackageStartupMessages({ library(glmnet); library(limma) })

artifact_paths <- list.files(path_repo(config$mouse$artifact_dir), pattern = "integrated_bulk_pseudobulk__.*deployment_all_control_like[.]rds$", full.names = TRUE)
if (length(artifact_paths) != 2L) stop("Artifact cardinality gate failed.")
arts <- lapply(artifact_paths, readRDS)
names(arts) <- vapply(arts, `[[`, character(1), "model_family")
if (!identical(sort(names(arts)), c("elasticnet", "pca_glmnet"))) stop("Expected ElasticNet and PCA/glmnet artifacts.")

manifest <- rbindlist(lapply(seq_along(arts), function(i) {
  a <- arts[[i]]
  if (a$schema_version != "integrated_ovary_clock_artifact_v1" ||
      a$artifact_scope != "deployment_all_control_like" ||
      a$missing_gene_policy != "impute_training_center_then_scale_to_zero") stop("Artifact schema/policy gate failed.")
  if (length(a$selected_genes) != 2000L || anyDuplicated(a$selected_genes)) stop("2,000 unique-gene gate failed.")
  if (!identical(names(a$gene_center), a$selected_genes) || !identical(names(a$gene_scale), a$selected_genes)) stop("Center/scale order gate failed.")
  if (a$model_family == "pca_glmnet" && !identical(a$model$selected_pc_names, paste0("PC", 1:7))) stop("PC1-PC7 gate failed.")
  data.table(model_family = a$model_family, selector = a$selector, artifact_path = normalizePath(artifact_paths[i]),
             schema_version = a$schema_version, n_selected_genes = length(a$selected_genes),
             n_reference_samples = length(a$training_sample_keys), missing_gene_policy = a$missing_gene_policy,
             calibration_intercept = a$calibration$intercept, calibration_slope = a$calibration$slope)
}))
write_tsv(manifest, path_module("data_manifest", "model_artifact_manifest.tsv"))

ref <- load_mouse_reference(arts$elasticnet)
if (!identical(ref$metadata$sample_key, arts$pca_glmnet$training_sample_keys)) stop("Mouse model reference sample mismatch.")
saveRDS(ref, path_module("derived_data/mouse_models", "mouse_deployment_reference_expression.rds"))

repro <- list()
weights <- list()
for (family in names(arts)) {
  a <- arts[[family]]
  if (!all(a$selected_genes %in% rownames(ref$expression))) stop("Reference matrix lacks model genes for ", family)
  x <- scale_with_artifact(ref$expression, a)
  raw_predict <- predict_artifact(a, x, calibrate = FALSE)
  if (family == "elasticnet") {
    co <- as.data.table(a$model$coefficient_table)
    b <- setNames(co$coefficient, co$gene_id)[a$selected_genes]
    raw_formula <- a$model$glmnet_intercept + as.numeric(x %*% b)
    pc_gene_error <- NA_real_
  } else {
    pc <- predict(a$model$pca_fit, newdata = x)[, a$model$selected_pc_names, drop = FALSE]
    pcb <- as.numeric(coef(a$model$glmnet_fit))[-1]
    names(pcb) <- rownames(coef(a$model$glmnet_fit))[-1]
    raw_formula <- as.numeric(coef(a$model$glmnet_fit)[1]) + as.numeric(pc %*% pcb[a$model$selected_pc_names])
    rot <- a$model$pca_fit$rotation[a$selected_genes, a$model$selected_pc_names, drop = FALSE]
    wscaled <- as.numeric(rot %*% pcb[a$model$selected_pc_names]); names(wscaled) <- a$selected_genes
    pca_center_offset <- sum(a$model$pca_fit$center[a$selected_genes] * wscaled)
    raw_gene <- as.numeric(coef(a$model$glmnet_fit)[1]) - pca_center_offset + as.numeric(x %*% wscaled)
    pc_gene_error <- max(abs(raw_formula - raw_gene))
  }
  err <- max(abs(raw_predict - raw_formula))
  tol <- config$analysis$hard_reproduction_tolerance
  if (!is.finite(err) || err > tol || (!is.na(pc_gene_error) && pc_gene_error > tol)) stop("Frozen reproduction gate failed for ", family)
  repro[[family]] <- data.table(model_family = family, n_checked = length(raw_predict), max_abs_predict_vs_formula = err,
                                max_abs_pc_vs_gene_space = pc_gene_error, tolerance = tol, status = "PASS")
  w <- model_gene_weights(a)
  weights[[family]] <- data.table(mouse_gene_symbol = names(w), model_family = family,
                                 expression_space_weight = as.numeric(w), abs_expression_space_weight = abs(as.numeric(w)),
                                 selected_gene = TRUE, nonzero = as.numeric(w) != 0)
}
write_tsv(rbindlist(repro, fill = TRUE), path_module("results/qc", "mouse_artifact_reproduction.tsv"))

shared <- intersect(arts$elasticnet$selected_genes, arts$pca_glmnet$selected_genes)
write_tsv(data.table(prompt_expected_shared_gene_count = 277L,
                     artifact_observed_shared_gene_count = length(shared),
                     authoritative_source = "current frozen top-level deployment artifacts",
                     resolution = "artifact-derived count used; prompt count recorded as stale"),
          path_module("results/qc", "prompt_vs_artifact_shared_gene_discrepancy.tsv"))

design <- model.matrix(~ age_week + source_batch, data = as.data.frame(ref$metadata))
fit <- eBayes(lmFit(ref$expression, design))
age_col <- which(colnames(design) == "age_week")
age_tab <- data.table(mouse_gene_symbol = rownames(ref$expression), mouse_age_beta = fit$coefficients[, age_col],
                      mouse_age_t = fit$t[, age_col], mouse_age_p = fit$p.value[, age_col],
                      mouse_age_fdr = p.adjust(fit$p.value[, age_col], "BH"), n_mouse_samples = ncol(ref$expression),
                      source_table = "recomputed_frozen_117_reference_limma")

master <- Reduce(function(x, y) merge(x, y, by = "mouse_gene_symbol", all = TRUE),
                 list(age_tab,
                      weights$elasticnet[, .(mouse_gene_symbol, elasticnet_expression_space_weight = expression_space_weight, elasticnet_nonzero = nonzero)],
                      weights$pca_glmnet[, .(mouse_gene_symbol, pca_expression_space_weight = expression_space_weight)]))
master[, `:=`(integrated_gene_universe = TRUE,
              elasticnet_hvg_2000 = mouse_gene_symbol %in% arts$elasticnet$selected_genes,
              pca_versionB_2000 = mouse_gene_symbol %in% arts$pca_glmnet$selected_genes,
              deployment_shared_artifact = mouse_gene_symbol %in% shared,
              deployment_shared_artifact_n = length(shared))]
setcolorder(master, c("mouse_gene_symbol", "integrated_gene_universe", "elasticnet_hvg_2000", "pca_versionB_2000", "deployment_shared_artifact", "deployment_shared_artifact_n"))
write_tsv(master, path_module("derived_data/mouse_models", "mouse_gene_direction_and_weight_master.tsv"))
saveRDS(arts, path_module("derived_data/mouse_models", "authoritative_artifacts.rds"))
message("Mouse artifacts reproduced and gene weights extracted.")
