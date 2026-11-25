#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(glmnet)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260707)

root_candidates <- c(
  Sys.getenv("OVARY_CLOCK_ROOT"),
  Sys.getenv("OVARY_CLOCK_ROOT")
)
root_dir <- normalizePath(root_candidates[file.exists(root_candidates)][[1]])
analysis_dir <- file.path(root_dir, "adult_c57_ovary_clocks")
scrna_results_dir <- file.path(analysis_dir, "7 scRNAseq_clock_results")
level_dir <- file.path(scrna_results_dir, "clocks", "metacell_level_scRNA")
# v2: read v2 metacells, write v2 models
V2 <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
previous_root <- file.path(V2, "06_clock_run", Sys.getenv("V2_BUILD_DIR", "build"))  # v2 sensitivity
run_root <- file.path(V2, "06_clock_run", Sys.getenv("V2_CV_DIR", "cv"))  # v2 sensitivity
fold_mode <- Sys.getenv("V2_FOLD_MODE", "metacell")  # v2 sensitivity: "dataset" = leave-one-dataset-out

step01_dir <- file.path(run_root, "step_01_metadata")
step02_dir <- file.path(run_root, "step_02_preprocessing_filtering")
step03_dir <- file.path(run_root, "step_03_metacell_cv_split")
model_dir <- file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet")
validation_dir <- file.path(run_root, "step_05_internal_cv_validation")
projection_dir <- file.path(run_root, "step_06_projection_stress_test")
check_dir <- file.path(run_root, "step_07_hard_checks")
release_dir <- file.path(run_root, "step_08_figures_tables_release")
figure_dir <- file.path(release_dir, "figures")
table_dir <- file.path(release_dir, "tables")
plot_data_dir <- file.path(release_dir, "plot_data")

for (dd in c(step01_dir, step02_dir, step03_dir, model_dir, validation_dir,
             projection_dir, check_dir, release_dir, figure_dir, table_dir, plot_data_dir)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

old_preprocess_dir <- file.path(previous_root, "step_02_preprocessing_filtering")
old_counts_path <- file.path(old_preprocess_dir, "filtered_metacell_counts.rds")
old_meta_path <- file.path(old_preprocess_dir, "filtered_metacell_metadata.tsv")
if (!file.exists(old_counts_path) || !file.exists(old_meta_path)) {
  stop("Missing reusable preprocessing outputs in previous run: ", old_preprocess_dir, call. = FALSE)
}

writeLines(
  c(
    "# Stale/superseded run",
    "",
    "This run used biological-sample-level train/test splitting and is not canonical.",
    "It is superseded by:",
    file.path(level_dir, "scRNA_metacell_elasticnet_metacell_cv_20260707"),
    "",
    "Canonical evaluation is now metacell-level internal out-of-fold CV."
  ),
  file.path(previous_root, "STALE_SUPERSEDED_BY_scRNA_metacell_elasticnet_metacell_cv_20260707.md")
)

age_bin_levels <- c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
age_bin_labels <- c(
  young_adult = "Young adult",
  mature_adult = "Mature adult",
  advanced_maternal_age = "Advanced maternal age",
  post_reproductive_ageing = "Post-reproductive"
)
age_bin_shapes <- c(
  young_adult = 16,
  mature_adult = 17,
  advanced_maternal_age = 15,
  post_reproductive_ageing = 18
)
dataset_palette <- c(
  "2022_eLife_Morris_cycling_murine_ovary_atlas" = "#0072B2",
  "2024_Cell_Winkler_FRT_scRNA_spatial" = "#D55E00",
  "2024_LifeMedicine_Jiang_whole_life_mouse_ovary" = "#009E73",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "#CC79A7",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary" = "#7F3C8D",
  "2025_CommBio_Zhao_developing_ovary_scRNA_spatial" = "#11A579",
  "2026_CommBio_Galligos_ovarian_inflammaging" = "#E69F00",
  "2025_Kim_PRJNA863443_mouse_menopause_models" = "#56B4E9"
)

min_cv_folds <- 3L
preferred_cv_folds <- 5L
min_metacells_per_fold <- 10L
min_age_points <- 3L
min_model_metacells <- 30L
somatic_top_candidate_genes <- 1000L
somatic_feature_count_grid <- c(20L, 40L, 80L, 160L, 320L)
somatic_min_features <- 8L
somatic_cor_threshold <- 0.25
alpha_grid <- c(0, 0.5, 1)

theme_clock_plot <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold")
    )
}

logcpm_from_counts <- function(count_mat) {
  lib <- Matrix::colSums(count_mat)
  lib[lib == 0] <- NA_real_
  cpm <- t(t(count_mat) / lib) * 1e6
  cpm[is.na(cpm)] <- 0
  as.matrix(log1p(cpm))
}

read_prjna863443_gene_id_map <- function() {
  prjna_counts_dir <- file.path(
    analysis_dir,
    "data_input",
    "scRNA_dataset",
    "PRJNA863443",
    "cellranger_counts"
  )
  if (!dir.exists(prjna_counts_dir)) return(tibble(gene_id = character(), gene_symbol = character()))
  sample_dirs <- list.dirs(prjna_counts_dir, full.names = TRUE, recursive = FALSE)
  feature_paths <- file.path(sample_dirs, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")
  feature_paths <- feature_paths[file.exists(feature_paths)]
  if (length(feature_paths) == 0L) return(tibble(gene_id = character(), gene_symbol = character()))
  bind_rows(lapply(feature_paths, function(path) {
    fread(
      cmd = paste("gzip -cd", shQuote(path)),
      header = FALSE,
      col.names = c("gene_id", "gene_symbol", "feature_type")
    ) %>%
      filter(feature_type == "Gene Expression")
  })) %>%
    filter(
      grepl("^ENSMUSG", gene_id),
      !is.na(gene_symbol),
      gene_symbol != "",
      !grepl("^ENSMUSG", gene_symbol)
    ) %>%
    count(gene_id, gene_symbol, name = "n_feature_files") %>%
    arrange(gene_id, desc(n_feature_files), gene_symbol) %>%
    distinct(gene_id, .keep_all = TRUE) %>%
    select(gene_id, gene_symbol)
}

collapse_count_rows <- function(count_mat, target_gene_names) {
  if (length(target_gene_names) != nrow(count_mat)) {
    stop("target_gene_names/count_mat mismatch", call. = FALSE)
  }
  row_group <- factor(target_gene_names, levels = unique(target_gene_names))
  row_map <- sparseMatrix(
    i = seq_along(row_group),
    j = as.integer(row_group),
    x = 1,
    dims = c(length(row_group), nlevels(row_group))
  )
  collapsed <- t(row_map) %*% count_mat
  rownames(collapsed) <- levels(row_group)
  colnames(collapsed) <- colnames(count_mat)
  collapsed
}

harmonize_count_gene_namespace <- function(count_mat) {
  feature_map <- read_prjna863443_gene_id_map()
  if (nrow(feature_map) == 0L) {
    return(list(
      counts = count_mat,
      audit = tibble(
        metric = c("gene_namespace_harmonization_applied", "prjna863443_feature_map_rows"),
        value = c("FALSE", "0")
      )
    ))
  }
  map_vec <- setNames(feature_map$gene_symbol, feature_map$gene_id)
  old_names <- rownames(count_mat)
  new_names <- unname(map_vec[old_names])
  new_names <- ifelse(is.na(new_names) | new_names == "", old_names, new_names)
  n_ensembl_rows_before <- sum(grepl("^ENSMUSG", old_names))
  n_ensembl_rows_mapped <- sum(grepl("^ENSMUSG", old_names) & new_names != old_names)
  n_rows_before <- nrow(count_mat)
  n_duplicate_targets <- sum(duplicated(new_names))
  collapsed <- collapse_count_rows(count_mat, new_names)
  audit <- tibble(
    metric = c(
      "gene_namespace_harmonization_applied",
      "prjna863443_feature_map_rows",
      "count_rows_before_harmonization",
      "count_rows_after_harmonization",
      "ensembl_like_rows_before_harmonization",
      "ensembl_like_rows_mapped_to_symbol",
      "duplicate_target_symbol_rows_collapsed"
    ),
    value = as.character(c(
      TRUE,
      nrow(feature_map),
      n_rows_before,
      nrow(collapsed),
      n_ensembl_rows_before,
      n_ensembl_rows_mapped,
      n_duplicate_targets
    ))
  )
  list(counts = collapsed, audit = audit)
}

metric_table <- function(df, truth = "age_weeks", pred = "predicted_age_weeks", n_col = "metacell_id") {
  if (nrow(df) == 0) {
    return(tibble(n_metacells = 0L, mae_weeks = NA_real_, rmse_weeks = NA_real_,
                  mean_delta_age_weeks = NA_real_, slope = NA_real_, r2 = NA_real_,
                  pearson_r = NA_real_, spearman_r = NA_real_))
  }
  y <- df[[truth]]
  p <- df[[pred]]
  ok <- is.finite(y) & is.finite(p)
  fit <- if (sum(ok) >= 2 && dplyr::n_distinct(y[ok]) >= 2) lm(p[ok] ~ y[ok]) else NULL
  tibble(
    n_metacells = dplyr::n_distinct(df[[n_col]]),
    mae_weeks = mean(abs(p - y), na.rm = TRUE),
    rmse_weeks = sqrt(mean((p - y)^2, na.rm = TRUE)),
    mean_delta_age_weeks = mean(p - y, na.rm = TRUE),
    slope = if (is.null(fit)) NA_real_ else unname(coef(fit)[[2]]),
    r2 = if (is.null(fit)) NA_real_ else summary(fit)$r.squared,
    pearson_r = suppressWarnings(cor(y, p, method = "pearson", use = "complete.obs")),
    spearman_r = suppressWarnings(cor(y, p, method = "spearman", use = "complete.obs"))
  )
}

choose_outer_k <- function(df) {
  n <- nrow(df)
  n_age_points <- dplyr::n_distinct(df$age_weeks)
  n_age_bins <- dplyr::n_distinct(df$reproductive_age_bin)
  if (n < min_model_metacells || n_age_points < min_age_points || n_age_bins < 2) return(NA_integer_)
  age_counts <- df %>% count(reproductive_age_bin, name = "n")
  if (n >= preferred_cv_folds * min_metacells_per_fold && min(age_counts$n) >= preferred_cv_folds) {
    preferred_cv_folds
  } else if (n >= min_cv_folds * min_metacells_per_fold && min(age_counts$n) >= min_cv_folds) {
    min_cv_folds
  } else {
    NA_integer_
  }
}

assign_balanced_folds <- function(df, k) {
  fold_sizes <- setNames(rep(0L, k), paste0("fold", seq_len(k)))
  df <- df %>%
    mutate(stratum = paste(reproductive_age_bin, dataset_id, sep = "__")) %>%
    arrange(stratum, sample_key, metacell_id)
  pieces <- list()
  strata <- df %>% count(stratum, sort = TRUE)
  for (ss in strata$stratum) {
    cur <- df %>% filter(stratum == ss) %>% arrange(sample_key, metacell_id)
    if (nrow(cur) > 1) cur <- cur[sample(seq_len(nrow(cur))), , drop = FALSE]
    cur_fold <- character(nrow(cur))
    for (ii in seq_len(nrow(cur))) {
      target <- names(fold_sizes)[which.min(fold_sizes)]
      cur_fold[[ii]] <- target
      fold_sizes[[target]] <- fold_sizes[[target]] + 1L
    }
    cur$fold_id <- cur_fold
    pieces[[ss]] <- cur
  }
  bind_rows(pieces) %>%
    select(-stratum) %>%
    arrange(harmonized_cell_class, fold_id, metacell_id)
}

make_inner_foldid <- function(train_meta) {
  k <- if (nrow(train_meta) >= preferred_cv_folds * min_metacells_per_fold) preferred_cv_folds else min_cv_folds
  if (nrow(train_meta) < k * 2) return(NULL)
  inner <- assign_balanced_folds(train_meta, k)
  as.integer(factor(inner$fold_id[match(train_meta$metacell_id, inner$metacell_id)], levels = paste0("fold", seq_len(k))))
}

fit_alpha_grid <- function(x, y, foldid = NULL) {
  best <- NULL
  records <- list()
  nfolds <- if (is.null(foldid)) min(preferred_cv_folds, nrow(x)) else length(unique(foldid))
  if (nrow(x) < 3 || nfolds < 2) return(NULL)
  for (a in alpha_grid) {
    fit_cv <- tryCatch(
      cv.glmnet(
        x = x,
        y = y,
        alpha = a,
        family = "gaussian",
        foldid = foldid,
        nfolds = nfolds,
        type.measure = "mae",
        standardize = TRUE
      ),
      error = function(e) NULL
    )
    if (is.null(fit_cv)) next
    mae <- min(fit_cv$cvm, na.rm = TRUE)
    records[[length(records) + 1L]] <- tibble(alpha = a, lambda_min = fit_cv$lambda.min, cv_mae_weeks = mae)
    if (is.null(best) || mae < best$cv_mae_weeks) {
      best <- list(alpha = a, lambda = fit_cv$lambda.min, cv_mae_weeks = mae)
    }
  }
  if (is.null(best)) return(NULL)
  list(best = best, cv_table = bind_rows(records))
}

fit_metacell_elasticnet <- function(train_mat, train_meta, predict_mat = NULL) {
  if (ncol(train_mat) != nrow(train_meta)) stop("train_mat/train_meta mismatch", call. = FALSE)
  if (nrow(train_meta) < min_model_metacells || dplyr::n_distinct(train_meta$age_weeks) < min_age_points) return(NULL)

  feature_values <- as.matrix(train_mat)
  gene_stats <- tibble(
    gene = rownames(feature_values),
    mean_expression = rowMeans(feature_values),
    expression_fraction = rowMeans(feature_values > 0),
    training_sd = apply(t(feature_values), 2, sd)
  ) %>%
    mutate(
      passes_variance_filter = is.finite(training_sd) & training_sd > 0,
      passes_low_expression_filter = is.finite(mean_expression) & expression_fraction > 0,
      passes_prefilter = passes_variance_filter & passes_low_expression_filter
    )
  feature_diagnostics <- gene_stats %>%
    filter(passes_prefilter) %>%
    arrange(desc(training_sd), gene) %>%
    mutate(variability_rank = row_number())
  if (nrow(feature_diagnostics) < somatic_min_features) return(NULL)

  candidates <- feature_diagnostics$gene[seq_len(min(somatic_top_candidate_genes, nrow(feature_diagnostics)))]
  rho <- suppressWarnings(cor(t(feature_values[candidates, , drop = FALSE]), train_meta$age_weeks, method = "spearman"))
  feature_tbl <- feature_diagnostics %>%
    filter(gene %in% candidates) %>%
    left_join(tibble(gene = candidates, rho = as.numeric(rho)), by = "gene") %>%
    filter(is.finite(rho)) %>%
    mutate(abs_rho = abs(rho)) %>%
    arrange(desc(abs_rho), gene)
  if (nrow(feature_tbl) < somatic_min_features) return(NULL)

  threshold_tbl <- feature_tbl %>% filter(abs_rho >= somatic_cor_threshold)
  feature_cap_grid <- sort(unique(pmin(somatic_feature_count_grid[somatic_feature_count_grid >= somatic_min_features], nrow(feature_tbl))))
  tuning_records <- list()
  best_tuning <- NULL
  foldid <- make_inner_foldid(train_meta)
  for (cap in feature_cap_grid) {
    if (nrow(threshold_tbl) >= somatic_min_features) {
      selected_cap <- threshold_tbl %>% slice_head(n = min(cap, nrow(threshold_tbl)))
      selection_mode <- "rho_threshold_then_top_abs_rho"
    } else {
      selected_cap <- feature_tbl %>% slice_head(n = min(cap, nrow(feature_tbl)))
      selection_mode <- "fallback_top_abs_rho_insufficient_threshold_features"
    }
    x_train_cap <- t(as.matrix(train_mat[selected_cap$gene, , drop = FALSE]))
    keep <- apply(x_train_cap, 2, sd) > 0
    x_train_cap <- x_train_cap[, keep, drop = FALSE]
    features_cap <- colnames(x_train_cap)
    if (length(features_cap) < somatic_min_features) next
    alpha_fit <- fit_alpha_grid(x_train_cap, train_meta$age_weeks, foldid = foldid)
    if (is.null(alpha_fit)) next
    tuning_records[[length(tuning_records) + 1L]] <- alpha_fit$cv_table %>%
      mutate(feature_cap = cap, n_features_after_variance_filter = length(features_cap),
             selection_mode = selection_mode, .before = alpha)
    local_best <- alpha_fit$best
    local_best$feature_cap <- cap
    local_best$features <- features_cap
    local_best$selected <- selected_cap %>% filter(gene %in% features_cap)
    local_best$selection_mode <- selection_mode
    if (is.null(best_tuning) || local_best$cv_mae_weeks < best_tuning$cv_mae_weeks) best_tuning <- local_best
  }
  if (is.null(best_tuning)) return(NULL)

  features <- best_tuning$features
  x_train <- t(as.matrix(train_mat[features, , drop = FALSE]))
  final_fit <- glmnet(
    x_train,
    train_meta$age_weeks,
    alpha = best_tuning$alpha,
    lambda = best_tuning$lambda,
    family = "gaussian",
    standardize = TRUE
  )
  coef_vec <- as.numeric(coef(final_fit))[-1]
  names(coef_vec) <- features
  train_pred <- as.numeric(predict(final_fit, newx = x_train))
  predict_pred <- numeric(0)
  if (!is.null(predict_mat) && ncol(predict_mat) > 0) {
    x_predict <- t(as.matrix(predict_mat[features, , drop = FALSE]))
    predict_pred <- as.numeric(predict(final_fit, newx = x_predict))
  }

  list(
    fit = final_fit,
    features = features,
    train_prediction = train_pred,
    predict_prediction = predict_pred,
    best_alpha = best_tuning$alpha,
    best_lambda = best_tuning$lambda,
    cv_mae_weeks = best_tuning$cv_mae_weeks,
    selected_feature_cap = best_tuning$feature_cap,
    cv_grid = bind_rows(tuning_records),
    feature_table = best_tuning$selected %>%
      mutate(
        selected_feature_cap = best_tuning$feature_cap,
        selection_mode = best_tuning$selection_mode,
        feature_source = "metacell_training_only"
      ),
    coefficient_table = tibble(gene = features, coefficient = coef_vec) %>%
      mutate(
        abs_coefficient = abs(coefficient),
        nonzero_coefficient = coefficient != 0,
        coefficient_direction = case_when(coefficient > 0 ~ "positive", coefficient < 0 ~ "negative", TRUE ~ "zero")
      ) %>%
      arrange(desc(abs_coefficient), gene) %>%
      mutate(abs_coefficient_rank = row_number()) %>%
      left_join(feature_tbl %>% select(gene, mean_expression, expression_fraction, training_sd, rho, abs_rho), by = "gene") %>%
      mutate(
        best_alpha = best_tuning$alpha,
        best_lambda = best_tuning$lambda,
        fold_cv_mae_weeks = best_tuning$cv_mae_weeks,
        selected_feature_cap = best_tuning$feature_cap,
        feature_source = "metacell_training_only"
      )
  )
}

message("Loading reusable filtered metacells")
metacell_counts <- readRDS(old_counts_path)
metacell_meta <- fread(old_meta_path) %>%
  as_tibble() %>%
  mutate(
    legacy_sample_split_set_not_used = if ("split_set" %in% colnames(.)) split_set else NA_character_,
    cv_role = case_when(
      model_role == "clock_candidate" &
        intervention_class == "natural_aging" &
        !is.na(age_weeks) &
        reproductive_age_bin != "out_of_scope" ~ "metacell_cv_candidate",
      model_role == "projection_only" ~ "projection_only",
      TRUE ~ "diagnostic_only"
    ),
    reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels)
  ) %>%
  select(-any_of("split_set"))
stopifnot(identical(metacell_meta$metacell_id, colnames(metacell_counts)))

message("Harmonizing count-matrix gene namespace")
gene_namespace_harmonization <- harmonize_count_gene_namespace(metacell_counts)
metacell_counts <- gene_namespace_harmonization$counts

fwrite(
  metacell_meta %>%
    distinct(dataset_folder, dataset_id, sample_id, sample_key, age_weeks, reproductive_age_bin,
             intervention_class, model_role, include_clock_candidate, projection_role, cv_role),
  file.path(step01_dir, "scRNA_sample_role_metadata.tsv"),
  sep = "\t"
)
fwrite(
  metacell_meta %>%
    filter(cv_role == "metacell_cv_candidate") %>%
    count(reproductive_age_bin, name = "n_metacells"),
  file.path(step01_dir, "model_candidate_metacell_age_bin_counts.tsv"),
  sep = "\t"
)

message("Writing preprocessing carry-forward outputs")
copy_candidates <- c(
  "confident_celltype_filter_summary.tsv",
  "confident_celltype_filter_audit.tsv"
)
for (ff in copy_candidates) {
  src <- file.path(old_preprocess_dir, ff)
  if (file.exists(src)) file.copy(src, file.path(step02_dir, ff), overwrite = TRUE)
}
saveRDS(metacell_counts, file.path(step02_dir, "filtered_metacell_counts.rds"))
fwrite(
  gene_namespace_harmonization$audit,
  file.path(step02_dir, "gene_namespace_harmonization_audit.tsv"),
  sep = "\t"
)
fwrite(metacell_meta, file.path(step02_dir, "filtered_metacell_metadata.tsv"), sep = "\t")
writeLines(
  c(
    "# scRNA metacell ElasticNet preprocessing",
    "",
    "This step carries forward the confident-cell-type metacells from the superseded run.",
    "Only the split/modeling layer was superseded; the confident-cell-type filtering and metacell construction are reused.",
    "",
    "Filters applied before model building:",
    "- curation_status in confirmed, probable",
    "- label_assignment_confidence >= 0.70",
    "- exclude Unknown and Erythroid harmonized_cell_class",
    "- metacell_size >= 25 cells",
    "- retain harmonized_cell_class labels only",
    "- model-building metacells require nonmissing age_weeks and reproductive_age_bin != out_of_scope",
    "- PRJNA863443 VCD/Foxl2/control samples are retained only as projection_only rows",
    "",
    "The legacy biological-sample split column from the superseded run is not used for modeling.",
    "Canonical evaluation starts in step_03_metacell_cv_split and uses metacell_id-level folds."
  ),
  file.path(step02_dir, "README.md")
)

message("Assigning metacell-level CV folds")
cv_candidates <- metacell_meta %>%
  filter(cv_role == "metacell_cv_candidate") %>%
  mutate(reproductive_age_bin = as.character(reproductive_age_bin))

celltype_eligibility <- cv_candidates %>%
  group_by(harmonized_cell_class) %>%
  group_modify(~{
    kk <- choose_outer_k(.x)
    tibble(
      n_metacells = nrow(.x),
      n_samples = n_distinct(.x$sample_key),
      n_datasets = n_distinct(.x$dataset_id),
      n_age_bins = n_distinct(.x$reproductive_age_bin),
      n_age_points = n_distinct(.x$age_weeks),
      n_cells = sum(.x$n_cells),
      outer_cv_folds = kk,
      eligible_for_metacell_cv = !is.na(kk),
      ineligibility_reason = case_when(
        !is.na(kk) ~ NA_character_,
        nrow(.x) < min_model_metacells ~ "insufficient_model_metacells",
        n_distinct(.x$age_weeks) < min_age_points ~ "insufficient_age_points",
        n_distinct(.x$reproductive_age_bin) < 2 ~ "insufficient_age_bin_coverage",
        TRUE ~ "cannot_support_at_least_3_fold_metacell_cv"
      )
    )
  }) %>%
  ungroup()
fwrite(celltype_eligibility, file.path(step03_dir, "celltype_metacell_cv_eligibility.tsv"), sep = "\t")

fold_assignments <- bind_rows(lapply(celltype_eligibility$harmonized_cell_class[celltype_eligibility$eligible_for_metacell_cv], function(ct) {
  ct_meta <- cv_candidates %>% filter(harmonized_cell_class == ct)
  kk <- celltype_eligibility$outer_cv_folds[celltype_eligibility$harmonized_cell_class == ct][[1]]
  if (fold_mode == "dataset") {  # v2 sensitivity: one fold per dataset
    ds_levels <- sort(unique(ct_meta$dataset_id))
    kk <- length(ds_levels)
    return(ct_meta %>% mutate(fold_id = factor(paste0("fold", match(dataset_id, ds_levels)), levels = paste0("fold", seq_len(kk))),
                              heldout_dataset = dataset_id, split_id = "leave_one_dataset_out", split_unit = "dataset_id"))
  }
  assign_balanced_folds(ct_meta, kk) %>%
    mutate(
      split_id = "metacell_internal_cv",
      split_unit = "metacell_id",
      fold_id = factor(fold_id, levels = paste0("fold", seq_len(kk)))
    )
}))
fwrite(fold_assignments, file.path(step03_dir, "metacell_fold_assignments.tsv"), sep = "\t")

fold_summary <- fold_assignments %>%
  count(harmonized_cell_class, fold_id, name = "n_metacells") %>%
  group_by(harmonized_cell_class) %>%
  mutate(
    min_fold_metacells = min(n_metacells),
    max_fold_metacells = max(n_metacells),
    fold_size_range = max_fold_metacells - min_fold_metacells
  ) %>%
  ungroup()
fold_summary_age <- fold_assignments %>%
  count(harmonized_cell_class, fold_id, reproductive_age_bin, name = "n_metacells")
fold_summary_dataset <- fold_assignments %>%
  count(harmonized_cell_class, fold_id, dataset_id, name = "n_metacells")
sample_overlap_expected <- fold_assignments %>%
  distinct(harmonized_cell_class, sample_key, sample_id, dataset_id, fold_id) %>%
  group_by(harmonized_cell_class, sample_key, sample_id, dataset_id) %>%
  summarise(
    n_cv_folds_with_sample_metacells = n_distinct(fold_id),
    cv_folds = paste(sort(unique(as.character(fold_id))), collapse = ";"),
    expected_overlap_across_metacell_folds = n_cv_folds_with_sample_metacells > 1,
    .groups = "drop"
  )
fwrite(fold_summary, file.path(step03_dir, "fold_summary_by_celltype.tsv"), sep = "\t")
fwrite(fold_summary_age, file.path(step03_dir, "fold_summary_by_celltype_age_bin.tsv"), sep = "\t")
fwrite(fold_summary_dataset, file.path(step03_dir, "fold_summary_by_celltype_dataset.tsv"), sep = "\t")
fwrite(sample_overlap_expected, file.path(step03_dir, "sample_overlap_expected_audit.tsv"), sep = "\t")

message("Fitting metacell-level cell-type ElasticNet CV clocks")
oof_predictions <- list()
selected_features <- list()
coefficients <- list()
cv_grid <- list()
training_manifest <- list()
final_models <- list()
projection_predictions <- list()

for (ct in celltype_eligibility$harmonized_cell_class[celltype_eligibility$eligible_for_metacell_cv]) {
  message("Cell type: ", ct)
  ct_meta_all <- metacell_meta %>% filter(harmonized_cell_class == ct)
  ct_counts <- metacell_counts[, ct_meta_all$metacell_id, drop = FALSE]
  ct_log <- logcpm_from_counts(ct_counts)
  rownames(ct_log) <- rownames(ct_counts)
  colnames(ct_log) <- colnames(ct_counts)

  ct_folds <- fold_assignments %>% filter(harmonized_cell_class == ct)
  for (fold in sort(unique(as.character(ct_folds$fold_id)))) {
    test_meta <- ct_folds %>% filter(as.character(fold_id) == fold)
    train_meta <- ct_folds %>% filter(as.character(fold_id) != fold)
    train_mat <- ct_log[, train_meta$metacell_id, drop = FALSE]
    test_mat <- ct_log[, test_meta$metacell_id, drop = FALSE]
    fit <- fit_metacell_elasticnet(train_mat, train_meta, test_mat)
    if (is.null(fit)) next

    oof_predictions[[paste(ct, fold, sep = "__")]] <- test_meta %>%
      transmute(
        split_id = "metacell_internal_cv",
        fold_id = fold,
        clock_type = "celltype_metacell_elasticnet_clock",
        model_evaluation_role = "metacell_oof_cv",
        metacell_id, sample_key, sample_id, dataset_folder, dataset_id,
        intervention_class, projection_role, harmonized_cell_class,
        age_weeks, reproductive_age_bin, n_cells,
        predicted_age_weeks_raw = fit$predict_prediction,
        predicted_age_weeks = fit$predict_prediction,
        delta_age_weeks = predicted_age_weeks - age_weeks,
        n_model_features = length(fit$features),
        best_alpha = fit$best_alpha,
        best_lambda = fit$best_lambda,
        inner_cv_mae_weeks = fit$cv_mae_weeks
      )
    selected_features[[paste(ct, fold, sep = "__")]] <- fit$feature_table %>%
      mutate(split_id = "metacell_internal_cv", fold_id = fold,
             clock_type = "celltype_metacell_elasticnet_clock",
             harmonized_cell_class = ct, .before = gene)
    coefficients[[paste(ct, fold, sep = "__")]] <- fit$coefficient_table %>%
      mutate(split_id = "metacell_internal_cv", fold_id = fold,
             clock_type = "celltype_metacell_elasticnet_clock",
             harmonized_cell_class = ct, .before = gene)
    cv_grid[[paste(ct, fold, sep = "__")]] <- fit$cv_grid %>%
      mutate(split_id = "metacell_internal_cv", fold_id = fold,
             clock_type = "celltype_metacell_elasticnet_clock",
             harmonized_cell_class = ct, .before = alpha)
    training_manifest[[paste(ct, fold, sep = "__")]] <- tibble(
      model_fit_role = "outer_cv_fold",
      split_id = "metacell_internal_cv",
      fold_id = fold,
      clock_type = "celltype_metacell_elasticnet_clock",
      harmonized_cell_class = ct,
      split_unit = "metacell_id",
      n_train_metacells = nrow(train_meta),
      n_test_metacells = nrow(test_meta),
      n_train_samples = n_distinct(train_meta$sample_key),
      n_test_samples = n_distinct(test_meta$sample_key),
      n_train_age_points = n_distinct(train_meta$age_weeks),
      n_model_features = length(fit$features),
      best_alpha = fit$best_alpha,
      best_lambda = fit$best_lambda,
      inner_cv_mae_weeks = fit$cv_mae_weeks
    )
  }

  final_train_meta <- ct_folds
  final_train_mat <- ct_log[, final_train_meta$metacell_id, drop = FALSE]
  proj_meta <- ct_meta_all %>% filter(model_role == "projection_only")
  proj_mat <- if (nrow(proj_meta) > 0) ct_log[, proj_meta$metacell_id, drop = FALSE] else NULL
  final_fit <- fit_metacell_elasticnet(final_train_mat, final_train_meta, proj_mat)
  if (!is.null(final_fit)) {
    final_models[[ct]] <- final_fit
    training_manifest[[paste(ct, "final", sep = "__")]] <- tibble(
      model_fit_role = "final_all_cv_metacells_for_projection",
      split_id = "metacell_internal_cv",
      fold_id = "final_all_metacells",
      clock_type = "celltype_metacell_elasticnet_clock",
      harmonized_cell_class = ct,
      split_unit = "metacell_id",
      n_train_metacells = nrow(final_train_meta),
      n_test_metacells = nrow(proj_meta),
      n_train_samples = n_distinct(final_train_meta$sample_key),
      n_test_samples = n_distinct(proj_meta$sample_key),
      n_train_age_points = n_distinct(final_train_meta$age_weeks),
      n_model_features = length(final_fit$features),
      best_alpha = final_fit$best_alpha,
      best_lambda = final_fit$best_lambda,
      inner_cv_mae_weeks = final_fit$cv_mae_weeks
    )
    if (nrow(proj_meta) > 0 && length(final_fit$predict_prediction) == nrow(proj_meta)) {
      projection_predictions[[ct]] <- proj_meta %>%
        transmute(
          split_id = "projection_only",
          fold_id = "final_all_metacells",
          clock_type = "celltype_metacell_elasticnet_clock",
          model_evaluation_role = "projection_stress_test",
          metacell_id, sample_key, sample_id, dataset_folder, dataset_id,
          intervention_class, projection_role, harmonized_cell_class,
          age_weeks, reproductive_age_bin, n_cells,
          predicted_age_weeks_raw = final_fit$predict_prediction,
          predicted_age_weeks = final_fit$predict_prediction,
          delta_age_weeks = predicted_age_weeks - age_weeks,
          n_model_features = length(final_fit$features),
          best_alpha = final_fit$best_alpha,
          best_lambda = final_fit$best_lambda,
          inner_cv_mae_weeks = final_fit$cv_mae_weeks
        )
    }
  }
  rm(ct_log)
  gc()
}

oof_tbl <- bind_rows(oof_predictions)
selected_features_tbl <- bind_rows(selected_features)
cv_grid_tbl <- bind_rows(cv_grid)
training_manifest_tbl <- bind_rows(training_manifest)
coefficients_tbl <- bind_rows(coefficients)
projection_tbl <- bind_rows(projection_predictions)

cv_eval <- oof_tbl %>%
  group_by(split_id, clock_type, harmonized_cell_class) %>%
  group_modify(~metric_table(.x)) %>%
  ungroup()

fwrite(oof_tbl, file.path(model_dir, "celltype_metacell_elasticnet_oof_predictions.tsv"), sep = "\t")
fwrite(cv_eval, file.path(model_dir, "celltype_metacell_elasticnet_cv_evaluation_summary.tsv"), sep = "\t")
fwrite(selected_features_tbl, file.path(model_dir, "celltype_metacell_elasticnet_selected_features_by_fold.tsv"), sep = "\t")
fwrite(cv_grid_tbl, file.path(model_dir, "celltype_metacell_elasticnet_cv_grid.tsv"), sep = "\t")
fwrite(training_manifest_tbl, file.path(model_dir, "model_training_manifest.tsv"), sep = "\t")
fwrite(coefficients_tbl, file.path(model_dir, "celltype_metacell_elasticnet_coefficients_by_fold.tsv"), sep = "\t")
saveRDS(final_models, file.path(model_dir, "celltype_metacell_elasticnet_final_models.rds"))

message("Writing validation outputs")
metacell_prediction_summary <- oof_tbl %>%
  select(split_id, fold_id, clock_type, model_evaluation_role, metacell_id, sample_key, sample_id,
         dataset_folder, dataset_id, intervention_class, projection_role, harmonized_cell_class,
         age_weeks, reproductive_age_bin, n_cells, predicted_age_weeks, delta_age_weeks,
         n_model_features, best_alpha, best_lambda, inner_cv_mae_weeks)
residual_bias_by_dataset <- oof_tbl %>%
  group_by(dataset_id) %>%
  group_modify(~metric_table(.x)) %>%
  ungroup()
residual_bias_by_age_bin <- oof_tbl %>%
  group_by(reproductive_age_bin) %>%
  group_modify(~metric_table(.x)) %>%
  ungroup()
residual_bias_by_celltype <- cv_eval
optional_sample_level_sanity <- oof_tbl %>%
  group_by(sample_key, sample_id, dataset_id, harmonized_cell_class, age_weeks, reproductive_age_bin) %>%
  summarise(
    predicted_age_weeks = weighted.mean(predicted_age_weeks, n_cells, na.rm = TRUE),
    delta_age_weeks = predicted_age_weeks - first(age_weeks),
    n_predicted_metacells = n(),
    n_predicted_cells = sum(n_cells),
    .groups = "drop"
  )
fwrite(metacell_prediction_summary, file.path(validation_dir, "metacell_level_prediction_summary.tsv"), sep = "\t")
fwrite(cv_eval, file.path(validation_dir, "celltype_clock_cv_evaluation_summary.tsv"), sep = "\t")
fwrite(residual_bias_by_dataset, file.path(validation_dir, "residual_bias_by_dataset.tsv"), sep = "\t")
fwrite(residual_bias_by_age_bin, file.path(validation_dir, "residual_bias_by_age_bin.tsv"), sep = "\t")
fwrite(residual_bias_by_celltype, file.path(validation_dir, "residual_bias_by_celltype.tsv"), sep = "\t")
fwrite(optional_sample_level_sanity, file.path(validation_dir, "optional_sample_level_sanity_summary.tsv"), sep = "\t")

message("Writing projection outputs")
if (nrow(projection_tbl) > 0) {
  projection_sample_celltype <- projection_tbl %>%
    group_by(split_id, fold_id, clock_type, model_evaluation_role, sample_key, sample_id, dataset_folder, dataset_id,
             intervention_class, projection_role, reproductive_age_bin, harmonized_cell_class, age_weeks) %>%
    summarise(
      predicted_age_weeks = weighted.mean(predicted_age_weeks, n_cells, na.rm = TRUE),
      delta_age_weeks = predicted_age_weeks - first(age_weeks),
      n_projected_metacells = n(),
      n_projected_cells = sum(n_cells),
      .groups = "drop"
    )
  projection_delta_by_celltype <- projection_tbl %>%
    group_by(intervention_class, projection_role, harmonized_cell_class) %>%
    summarise(
      n_projected_metacells = n(),
      n_projected_samples = n_distinct(sample_key),
      mean_delta_age_weeks = mean(delta_age_weeks, na.rm = TRUE),
      median_delta_age_weeks = median(delta_age_weeks, na.rm = TRUE),
      .groups = "drop"
    )
  vcd_summary <- projection_sample_celltype %>%
    filter(grepl("vcd", projection_role, ignore.case = TRUE)) %>%
    group_by(projection_role, intervention_class, harmonized_cell_class) %>%
    summarise(
      n_projected_samples = n_distinct(sample_key),
      n_projected_metacells = sum(n_projected_metacells),
      mean_delta_age_weeks = mean(delta_age_weeks, na.rm = TRUE),
      median_delta_age_weeks = median(delta_age_weeks, na.rm = TRUE),
      .groups = "drop"
    )
  foxl2_summary <- projection_sample_celltype %>%
    filter(grepl("foxl2", projection_role, ignore.case = TRUE)) %>%
    group_by(projection_role, intervention_class, harmonized_cell_class) %>%
    summarise(
      n_projected_samples = n_distinct(sample_key),
      n_projected_metacells = sum(n_projected_metacells),
      mean_delta_age_weeks = mean(delta_age_weeks, na.rm = TRUE),
      median_delta_age_weeks = median(delta_age_weeks, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  projection_sample_celltype <- tibble()
  projection_delta_by_celltype <- tibble()
  vcd_summary <- tibble()
  foxl2_summary <- tibble()
}
fwrite(projection_tbl, file.path(projection_dir, "prjna863443_projection_metacell_predictions.tsv"), sep = "\t")
fwrite(projection_sample_celltype, file.path(projection_dir, "prjna863443_projection_sample_celltype_summary.tsv"), sep = "\t")
fwrite(projection_delta_by_celltype, file.path(projection_dir, "prjna863443_projection_delta_by_celltype.tsv"), sep = "\t")
fwrite(vcd_summary, file.path(projection_dir, "vcd_projection_summary_by_celltype.tsv"), sep = "\t")
fwrite(foxl2_summary, file.path(projection_dir, "foxl2_projection_summary_by_celltype.tsv"), sep = "\t")

message("Building release figure")
normalize_dataset_id <- function(dataset_id, dataset_folder) {
  dataset_id <- as.character(dataset_id)
  dataset_folder <- as.character(dataset_folder)
  if_else(is.na(dataset_id) | dataset_id == "", dataset_folder, dataset_id)
}

facet_metrics <- cv_eval %>%
  mutate(
    facet_label = harmonized_cell_class,
    metric_label = sprintf("n=%s\nMAE=%.1f wk\nslope=%.2f\nR2=%.2f", n_metacells, mae_weeks, slope, r2)
  ) %>%
  select(harmonized_cell_class, facet_label, metric_label)
plot_data <- oof_tbl %>%
  left_join(facet_metrics, by = "harmonized_cell_class") %>%
  mutate(
    dataset_id = normalize_dataset_id(dataset_id, dataset_folder),
    reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels),
    dataset_id = factor(dataset_id, levels = names(dataset_palette))
  )
fwrite(plot_data, file.path(plot_data_dir, "filtered_elasticnet_metacell_internal_grouped_cv_plot_data.tsv"), sep = "\t")

metric_label_data <- facet_metrics %>%
  mutate(
    age_weeks = -Inf,
    predicted_age_weeks = Inf
  )
facet_labels <- setNames(facet_metrics$facet_label, facet_metrics$harmonized_cell_class)
p_cv <- ggplot(plot_data, aes(age_weeks, predicted_age_weeks, color = dataset_id, shape = reproductive_age_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
  geom_point(alpha = 0.55, size = 1.05) +
  geom_text(
    data = metric_label_data,
    aes(age_weeks, predicted_age_weeks, label = metric_label),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.04,
    size = 2.35,
    lineheight = 0.92,
    color = "grey20"
  ) +
  facet_wrap(~ harmonized_cell_class, scales = "free", labeller = labeller(harmonized_cell_class = facet_labels)) +
  scale_color_manual(values = dataset_palette, drop = FALSE) +
  scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
  theme_clock_plot(base_size = 8.5) +
  labs(
    title = "Filtered ElasticNet metacell clocks: internal grouped CV",
    x = "Observed age (weeks)",
    y = "Predicted age (weeks)",
    color = "Dataset",
    shape = "Age bin"
  )
ggsave(file.path(figure_dir, "filtered_elasticnet_metacell_internal_grouped_cv.pdf"), p_cv, width = 11, height = 8)
ggsave(file.path(figure_dir, "filtered_elasticnet_metacell_internal_grouped_cv.png"), p_cv, width = 11, height = 8, dpi = 240)

release_tables <- c(
  file.path(step01_dir, "scRNA_sample_role_metadata.tsv"),
  file.path(step01_dir, "model_candidate_metacell_age_bin_counts.tsv"),
  file.path(step02_dir, "confident_celltype_filter_summary.tsv"),
  file.path(step02_dir, "filtered_metacell_metadata.tsv"),
  file.path(step03_dir, "celltype_metacell_cv_eligibility.tsv"),
  file.path(step03_dir, "metacell_fold_assignments.tsv"),
  file.path(step03_dir, "fold_summary_by_celltype.tsv"),
  file.path(step03_dir, "fold_summary_by_celltype_age_bin.tsv"),
  file.path(step03_dir, "fold_summary_by_celltype_dataset.tsv"),
  file.path(step03_dir, "sample_overlap_expected_audit.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_oof_predictions.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_cv_evaluation_summary.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_selected_features_by_fold.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_coefficients_by_fold.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_cv_grid.tsv"),
  file.path(model_dir, "model_training_manifest.tsv"),
  file.path(validation_dir, "metacell_level_prediction_summary.tsv"),
  file.path(validation_dir, "celltype_clock_cv_evaluation_summary.tsv"),
  file.path(validation_dir, "residual_bias_by_dataset.tsv"),
  file.path(validation_dir, "residual_bias_by_age_bin.tsv"),
  file.path(validation_dir, "residual_bias_by_celltype.tsv"),
  file.path(validation_dir, "optional_sample_level_sanity_summary.tsv"),
  file.path(projection_dir, "prjna863443_projection_metacell_predictions.tsv"),
  file.path(projection_dir, "prjna863443_projection_sample_celltype_summary.tsv"),
  file.path(projection_dir, "prjna863443_projection_delta_by_celltype.tsv"),
  file.path(projection_dir, "vcd_projection_summary_by_celltype.tsv"),
  file.path(projection_dir, "foxl2_projection_summary_by_celltype.tsv")
)
file.copy(release_tables[file.exists(release_tables)], table_dir, overwrite = TRUE)

message("Running hard checks")
all_paths <- list.files(run_root, recursive = TRUE, full.names = FALSE, all.files = FALSE)
dir_paths <- list.dirs(run_root, recursive = TRUE, full.names = FALSE)
check_tbl <- tibble(check_id = character(), check_name = character(), passed = logical(), detail = character())
add_check <- function(id, name, passed, detail = "") {
  check_tbl <<- bind_rows(check_tbl, tibble(check_id = id, check_name = name, passed = isTRUE(passed), detail = detail))
}

old_stale_file <- file.path(previous_root, "STALE_SUPERSEDED_BY_scRNA_metacell_elasticnet_metacell_cv_20260707.md")
add_check("01", "previous_sample_split_rerun_marked_stale", file.exists(old_stale_file), old_stale_file)
add_check("02", "new_canonical_run_uses_metacell_id_split_unit",
          nrow(fold_assignments) > 0 && all(fold_assignments$split_unit == "metacell_id"),
          paste(unique(fold_assignments$split_unit), collapse = ";"))
fold_balance_check <- fold_summary %>%
  group_by(harmonized_cell_class) %>%
  summarise(max_range = max(fold_size_range), .groups = "drop")
add_check("03", "fold_sizes_balanced_by_metacell_counts_within_celltype",
          nrow(fold_balance_check) > 0 && all(fold_balance_check$max_range <= 1),
          paste(fold_balance_check$harmonized_cell_class[fold_balance_check$max_range > 1], collapse = ";"))
bad_prjna_cv <- fold_assignments %>%
  filter(dataset_folder == "PRJNA863443", grepl("^VCD|^Foxl2", sample_id))
add_check("04", "no_prjna863443_vcd_or_foxl2_metacell_in_cv_folds",
          nrow(bad_prjna_cv) == 0,
          paste(unique(bad_prjna_cv$sample_id), collapse = ";"))
vcd_foxl2_proj <- projection_tbl %>%
  filter(dataset_folder == "PRJNA863443", grepl("^VCD|^Foxl2", sample_id))
vcd_foxl2_nonproj <- metacell_meta %>%
  filter(dataset_folder == "PRJNA863443", grepl("^VCD|^Foxl2", sample_id), cv_role != "projection_only")
add_check("05", "vcd_foxl2_metacells_only_projection_outputs",
          nrow(vcd_foxl2_proj) > 0 && nrow(vcd_foxl2_nonproj) == 0,
          paste(unique(vcd_foxl2_nonproj$sample_id), collapse = ";"))
add_check("06", "no_sample_level_train_test_split_used_as_primary_split",
          !file.exists(file.path(step01_dir, "global_train_test_split.tsv")) &&
            all(!grepl("global_age_bin_train_test|heldout_test|train_fit", all_paths)),
          paste(head(all_paths[grepl("global_age_bin_train_test|heldout_test|train_fit", all_paths)], 20), collapse = ";"))
add_check("07", "sample_leakage_check_not_failure_condition",
          file.exists(file.path(step03_dir, "sample_overlap_expected_audit.tsv")) &&
            any(sample_overlap_expected$expected_overlap_across_metacell_folds),
          "sample overlap across metacell folds is expected and audited")
add_check("08", "primary_evaluation_metacell_oof_cv",
          nrow(oof_tbl) > 0 && all(oof_tbl$model_evaluation_role == "metacell_oof_cv"),
          paste(unique(oof_tbl$model_evaluation_role), collapse = ";"))
add_check("09", "sample_level_summaries_optional_sanity_only",
          file.exists(file.path(validation_dir, "optional_sample_level_sanity_summary.tsv")) &&
            !any(grepl("sample_prediction_summary.tsv$", all_paths)),
          paste(head(all_paths[grepl("sample_prediction_summary.tsv$", all_paths)], 20), collapse = ";"))
add_check("10", "no_forbidden_model_outputs_in_new_canonical_root",
          !any(grepl("celltype_metacell_elasticnet_filtered|celltype_metacell_pca_glmnet|celltype_weighted_elasticnet|celltype_weighted_pca_glmnet|sample_pseudobulk|weighted", all_paths, ignore.case = TRUE)),
          paste(head(all_paths[grepl("celltype_metacell_elasticnet_filtered|celltype_metacell_pca_glmnet|celltype_weighted_elasticnet|celltype_weighted_pca_glmnet|sample_pseudobulk|weighted", all_paths, ignore.case = TRUE)], 20), collapse = ";"))
release_files <- list.files(release_dir, recursive = TRUE, full.names = FALSE)
add_check("11", "release_tables_plot_data_and_figures_copied_to_step_08",
          any(endsWith(release_files, ".tsv")) &&
            any(grepl("^plot_data/", release_files)) &&
            any(endsWith(release_files, ".pdf")) &&
            any(endsWith(release_files, ".png")),
          paste(release_files, collapse = ";"))

fwrite(check_tbl, file.path(check_dir, "hard_check_summary.tsv"), sep = "\t")
writeLines(
  c(
    "# scRNA metacell ElasticNet metacell-CV run",
    "",
    paste0("- Output root: ", run_root),
    "- Canonical object: cell-type-specific metacell-level ElasticNet aging clocks.",
    "- Primary split/evaluation: internal metacell_id-level out-of-fold CV.",
    "- Each validation point is one held-out metacell.",
    "- Performance n is the number of held-out metacells, not biological samples.",
    "- Biological sample IDs are allowed to appear in multiple folds through different metacells; this is expected and audited.",
    "- The previous sample-level train/test rerun is stale and noncanonical.",
    "- Final all-metacell models are used only for PRJNA863443 projection stress tests.",
    paste0("- Eligible cell classes: ", sum(celltype_eligibility$eligible_for_metacell_cv)),
    paste0("- OOF metacell predictions: ", nrow(oof_tbl)),
    paste0("- Hard checks passed: ", sum(check_tbl$passed), "/", nrow(check_tbl))
  ),
  file.path(run_root, "README.md")
)

if (!all(check_tbl$passed)) {
  warning("Some hard checks failed: ", paste(check_tbl$check_id[!check_tbl$passed], collapse = ", "))
}
message("Metacell CV run complete: ", run_root)
