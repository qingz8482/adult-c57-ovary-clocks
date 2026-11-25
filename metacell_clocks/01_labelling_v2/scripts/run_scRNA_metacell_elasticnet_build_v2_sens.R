#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(glmnet)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260707)

root_dir <- normalizePath(Sys.getenv("OVARY_CLOCK_ROOT"))
analysis_dir <- file.path(root_dir, "adult_c57_ovary_clocks")
scrna_results_dir <- file.path(analysis_dir, "7 scRNAseq_clock_results")
level_dir <- file.path(scrna_results_dir, "clocks", "metacell_level_scRNA")
# v2: outputs and label inputs redirected to the v2 rebuild (only changes vs 20260707 are marked "v2")
V2 <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
rerun_root <- file.path(V2, "06_clock_run", Sys.getenv("V2_BUILD_DIR", "build"))  # v2 sensitivity
metadata_out_dir <- file.path(rerun_root, "step_01_metadata")
preprocess_dir <- file.path(rerun_root, "step_02_preprocessing_filtering")
model_dir <- file.path(rerun_root, "step_03_model_building", "celltype_metacell_elasticnet")
projection_dir <- file.path(rerun_root, "step_04_projection_stress_test")
check_dir <- file.path(rerun_root, "step_07_hard_checks")
release_dir <- file.path(rerun_root, "step_08_figures_tables_release")
plot_dir <- file.path(release_dir, "figures")
table_dir <- file.path(release_dir, "tables")
for (dd in c(metadata_out_dir, preprocess_dir, model_dir, projection_dir, check_dir, release_dir, plot_dir, table_dir)) {
  dir.create(dd, recursive = TRUE, showWarnings = FALSE)
}

metadata_path <- file.path(analysis_dir, "data_input", "metadata_scRNA.tsv")
manifest_path <- file.path(analysis_dir, "data_input", "standardized_seurat", "00_metadata", "standardized_processed_manifest.tsv")
consensus_path <- file.path(V2, "04_consensus", "harmonized_annotation_consensus_v2.tsv.gz")  # v2
age_bin_reference_path <- file.path(analysis_dir, "bulk-pseudobulk_clock", "methods_tables", "15_original_age_bin_summary.tsv")

if (!file.exists(metadata_path)) stop("Missing metadata_scRNA.tsv: ", metadata_path, call. = FALSE)
if (!file.exists(manifest_path)) stop("Missing standardized manifest: ", manifest_path, call. = FALSE)
if (!file.exists(consensus_path)) stop("Missing harmonized annotation consensus: ", consensus_path, call. = FALSE)
if (!file.exists(age_bin_reference_path)) stop("Missing age-bin reference: ", age_bin_reference_path, call. = FALSE)

allowed_prjna_clock_samples <- c("AC_13weeks", "YF_1", "YF_2", "AC_76weeks", "AC_86weeks", "OF_1", "OF_2")
target_test_by_bin <- c(
  young_adult = 6L,
  mature_adult = 4L,
  advanced_maternal_age = 2L,
  post_reproductive_ageing = 1L
)
target_n_by_bin <- c(
  young_adult = 32L,
  mature_adult = 18L,
  advanced_maternal_age = 11L,
  post_reproductive_ageing = 6L
)

metacell_target_cells <- 50L
metacell_min_cells <- 25L
min_label_confidence <- 0.70
clock_label_statuses <- strsplit(Sys.getenv("V2_STATUSES", "confirmed,probable"), ",")[[1]]  # v2 sensitivity
exclude_granulosa_luteal <- Sys.getenv("V2_EXCLUDE_GRANULOSA_LUTEAL", "0") == "1"  # v2 sensitivity
clock_excluded_labels <- c("Unknown", "Erythroid")
min_train_samples <- 8L
min_train_age_points <- 3L
min_train_metacells <- 20L
somatic_top_candidate_genes <- 1000L
somatic_feature_count_grid <- c(20L, 40L, 80L, 160L, 320L)
somatic_min_features <- 8L
somatic_cor_threshold <- 0.25
alpha_grid <- c(0, 0.5, 1)

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == "" | x == "NA"] <- NA_character_
  x
}

parse_age_weeks <- function(age, stage = "") {
  txt <- tolower(paste(na_chr(age), na_chr(stage)))
  out <- rep(NA_real_, length(txt))
  hit <- regmatches(txt, regexec("(?:postnatal\\s*day|pd|p)\\s*(\\d+(?:\\.\\d+)?)", txt, perl = TRUE))
  has <- lengths(hit) >= 2
  out[has] <- vapply(hit[has], function(z) as.numeric(z[[2]]) / 7, numeric(1))
  hit <- regmatches(txt, regexec("(\\d+(?:\\.\\d+)?)\\s*-?\\s*day\\s*-?\\s*old", txt, perl = TRUE))
  has <- is.na(out) & lengths(hit) >= 2
  out[has] <- vapply(hit[has], function(z) as.numeric(z[[2]]) / 7, numeric(1))
  hit <- regmatches(txt, regexec("(\\d+(?:\\.\\d+)?)\\s*-\\s*(\\d+(?:\\.\\d+)?)\\s*weeks?", txt, perl = TRUE))
  has <- is.na(out) & lengths(hit) >= 3
  out[has] <- vapply(hit[has], function(z) mean(as.numeric(z[2:3])), numeric(1))
  hit <- regmatches(txt, regexec("(\\d+(?:\\.\\d+)?)\\s*-?\\s*weeks?", txt, perl = TRUE))
  has <- is.na(out) & lengths(hit) >= 2
  out[has] <- vapply(hit[has], function(z) as.numeric(z[[2]]), numeric(1))
  hit <- regmatches(txt, regexec("(\\d+(?:\\.\\d+)?)\\s*-?\\s*(?:m|months?)", txt, perl = TRUE))
  has <- is.na(out) & lengths(hit) >= 2
  out[has] <- vapply(hit[has], function(z) as.numeric(z[[2]]) * 30.4375 / 7, numeric(1))
  out
}

assign_reproductive_age_bin <- function(age_weeks) {
  dplyr::case_when(
    is.na(age_weeks) | age_weeks < 7 ~ "out_of_scope",
    age_weeks >= 7 & age_weeks < 16 ~ "young_adult",
    age_weeks >= 16 & age_weeks < 40 ~ "mature_adult",
    age_weeks >= 40 & age_weeks < 65 ~ "advanced_maternal_age",
    age_weeks >= 65 ~ "post_reproductive_ageing",
    TRUE ~ "out_of_scope"
  )
}

prefer_joined_column <- function(df, name) {
  y_name <- paste0(name, ".y")
  x_name <- paste0(name, ".x")
  fill_missing <- function(primary, fallback) {
    out <- primary
    miss <- is.na(out)
    if (!any(miss)) return(out)
    fill <- fallback
    if (is.numeric(out)) {
      fill <- suppressWarnings(as.numeric(fill))
    } else if (is.logical(out)) {
      fill <- suppressWarnings(as.logical(fill))
    } else {
      fill <- as.character(fill)
    }
    out[miss] <- fill[miss]
    out
  }
  if (y_name %in% colnames(df)) {
    out <- df[[y_name]]
    if (x_name %in% colnames(df)) out <- fill_missing(out, df[[x_name]])
    if (name %in% colnames(df)) out <- fill_missing(out, df[[name]])
    df[[name]] <- out
  } else if (name %in% colnames(df) && x_name %in% colnames(df)) {
    df[[name]] <- fill_missing(df[[name]], df[[x_name]])
  } else if (!name %in% colnames(df) && x_name %in% colnames(df)) {
    df[[name]] <- df[[x_name]]
  }
  df
}

standardized_dataset_id <- function(dataset_folder) {
  dplyr::if_else(
    dataset_folder == "PRJNA863443",
    "2025_Kim_PRJNA863443_mouse_menopause_models",
    dataset_folder
  )
}

metric_table <- function(df, truth = "age_weeks", pred = "predicted_age_weeks") {
  if (nrow(df) == 0) {
    return(tibble(n_samples = 0, mae_weeks = NA_real_, rmse_weeks = NA_real_,
                  mean_delta_age_weeks = NA_real_, pearson_r = NA_real_, spearman_r = NA_real_))
  }
  y <- df[[truth]]
  p <- df[[pred]]
  tibble(
    n_samples = dplyr::n_distinct(df$sample_key),
    mae_weeks = mean(abs(p - y), na.rm = TRUE),
    rmse_weeks = sqrt(mean((p - y)^2, na.rm = TRUE)),
    mean_delta_age_weeks = mean(p - y, na.rm = TRUE),
    pearson_r = suppressWarnings(cor(y, p, method = "pearson", use = "complete.obs")),
    spearman_r = suppressWarnings(cor(y, p, method = "spearman", use = "complete.obs"))
  )
}

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
  "2026_CommBio_Galligos_ovarian_inflammaging" = "#E69F00",
  "2025_Kim_PRJNA863443_mouse_menopause_models" = "#56B4E9"
)

format_for_age_bin_plots <- function(df) {
  df %>%
    mutate(
      reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels),
      dataset_id = factor(dataset_id, levels = names(dataset_palette))
    )
}

performance_table <- function(df) {
  df %>%
    group_by(harmonized_cell_class) %>%
    summarise(
      n_samples = n_distinct(sample_key),
      n_prediction_rows = n(),
      mae_weeks = mean(abs(predicted_age_weeks - age_weeks), na.rm = TRUE),
      slope = {
        ok <- is.finite(age_weeks) & is.finite(predicted_age_weeks)
        if (sum(ok) >= 2 && dplyr::n_distinct(age_weeks[ok]) >= 2) {
          unname(coef(lm(predicted_age_weeks[ok] ~ age_weeks[ok]))[[2]])
        } else {
          NA_real_
        }
      },
      r2 = {
        ok <- is.finite(age_weeks) & is.finite(predicted_age_weeks)
        if (sum(ok) >= 2 && dplyr::n_distinct(age_weeks[ok]) >= 2) {
          summary(lm(predicted_age_weeks[ok] ~ age_weeks[ok]))$r.squared
        } else {
          NA_real_
        }
      },
      .groups = "drop"
    ) %>%
    mutate(
      performance_label = sprintf(
        "n=%s; MAE=%.1f wk; slope=%.2f; R2=%.2f",
        n_samples, mae_weeks, slope, r2
      )
    )
}

theme_clock_plot <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      legend.key = element_blank()
    )
}

get_assay_matrix <- function(obj) {
  assay <- DefaultAssay(obj)
  if (exists("Layers", mode = "function") && exists("JoinLayers", mode = "function")) {
    layer_names <- tryCatch(Layers(obj[[assay]]), error = function(e) character())
    if (sum(grepl("^counts", layer_names)) > 1 || sum(grepl("^data", layer_names)) > 1) {
      obj <- JoinLayers(obj, assay = assay)
    }
  }
  mat <- tryCatch(GetAssayData(obj, assay = assay, layer = "counts"), error = function(e) NULL)
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    mat <- tryCatch(GetAssayData(obj, assay = assay, slot = "counts"), error = function(e) NULL)
  }
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) stop("Could not extract count matrix", call. = FALSE)
  as(mat, "dgCMatrix")
}

aggregate_columns <- function(mat, groups) {
  groups <- as.character(groups)
  levels <- unique(groups)
  design <- sparseMatrix(
    i = seq_len(ncol(mat)),
    j = match(groups, levels),
    x = 1,
    dims = c(ncol(mat), length(levels))
  )
  out <- mat %*% design
  colnames(out) <- levels
  out
}

align_sparse_rows <- function(mat, genes) {
  out <- Matrix(0, nrow = length(genes), ncol = ncol(mat), sparse = TRUE)
  rownames(out) <- genes
  colnames(out) <- colnames(mat)
  idx <- match(rownames(mat), genes)
  out[idx, ] <- mat
  out
}

combine_sparse_by_rows <- function(mats) {
  mats <- mats[!vapply(mats, is.null, logical(1))]
  genes <- sort(unique(unlist(lapply(mats, rownames))))
  n_cols <- vapply(mats, ncol, integer(1))
  i_parts <- vector("list", length(mats))
  j_parts <- vector("list", length(mats))
  x_parts <- vector("list", length(mats))
  col_start <- 1L
  for (idx in seq_along(mats)) {
    mat <- mats[[idx]]
    trip <- summary(mat)
    row_lookup <- match(rownames(mat), genes)
    i_parts[[idx]] <- row_lookup[trip$i]
    j_parts[[idx]] <- trip$j + col_start - 1L
    x_parts[[idx]] <- trip$x
    col_start <- col_start + ncol(mat)
  }
  out <- sparseMatrix(
    i = unlist(i_parts, use.names = FALSE),
    j = unlist(j_parts, use.names = FALSE),
    x = unlist(x_parts, use.names = FALSE),
    dims = c(length(genes), sum(n_cols))
  )
  rownames(out) <- genes
  colnames(out) <- unlist(lapply(mats, colnames), use.names = FALSE)
  out
}

logcpm_from_counts <- function(count_mat) {
  lib <- Matrix::colSums(count_mat)
  lib[lib == 0] <- NA_real_
  cpm <- t(t(count_mat) / lib) * 1e6
  cpm[is.na(cpm)] <- 0
  as.matrix(log1p(cpm))
}

make_group_foldid <- function(sample_ids, max_folds = 5L) {
  sample_ids <- as.character(sample_ids)
  sample_levels <- sort(unique(sample_ids))
  n_folds <- min(max_folds, length(sample_levels))
  if (n_folds < 2) return(NULL)
  fold_lookup <- setNames(rep(seq_len(n_folds), length.out = length(sample_levels)), sample_levels)
  unname(fold_lookup[sample_ids])
}

fit_alpha_grid <- function(x, y, foldid = NULL) {
  best <- NULL
  records <- list()
  nfolds <- if (is.null(foldid)) min(5L, nrow(x)) else length(unique(foldid))
  if (nrow(x) < 3 || nfolds < 2) return(NULL)
  for (a in alpha_grid) {
    fit_cv <- tryCatch(
      cv.glmnet(
        x = x, y = y, alpha = a, family = "gaussian", foldid = foldid,
        nfolds = nfolds, type.measure = "mae", standardize = TRUE
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

calibrate_predictions <- function(y_train, pred_train_raw, pred_new_raw) {
  if (length(unique(pred_train_raw)) < 2 || length(unique(y_train)) < 2) {
    return(list(train_prediction = pred_train_raw, new_prediction = pred_new_raw,
                calibration_method = "identity_insufficient_training_variation", calibration_fit = NULL))
  }
  fit <- lm(age_weeks ~ raw_prediction, data = data.frame(age_weeks = y_train, raw_prediction = pred_train_raw))
  list(
    train_prediction = as.numeric(predict(fit, newdata = data.frame(raw_prediction = pred_train_raw))),
    new_prediction = as.numeric(predict(fit, newdata = data.frame(raw_prediction = pred_new_raw))),
    calibration_method = "training_only_lm_age_weeks_on_raw_prediction",
    calibration_fit = fit
  )
}

fit_direct_elasticnet <- function(train_mat, train_meta, test_mat, test_meta, projection_mat = NULL, projection_meta = NULL) {
  train_units <- train_meta$sample_key
  if (n_distinct(train_units) < min_train_samples || n_distinct(train_meta$age_weeks) < min_train_age_points) return(NULL)
  if (ncol(train_mat) < min_train_metacells || nrow(test_meta) == 0) return(NULL)

  sample_means <- t(rowsum(t(as.matrix(train_mat)), group = train_units, reorder = FALSE))
  sample_counts <- as.numeric(table(factor(train_units, levels = colnames(sample_means))))
  sample_means <- t(t(sample_means) / sample_counts)
  y_feature <- train_meta %>%
    mutate(sample_unit = train_units) %>%
    distinct(sample_unit, age_weeks) %>%
    arrange(match(sample_unit, colnames(sample_means))) %>%
    pull(age_weeks)

  feature_values <- as.matrix(sample_means)
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
  rho <- suppressWarnings(cor(t(feature_values[candidates, , drop = FALSE]), y_feature, method = "spearman"))
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
  foldid <- make_group_foldid(train_units)
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
  final_fit <- glmnet(x_train, train_meta$age_weeks, alpha = best_tuning$alpha,
                      lambda = best_tuning$lambda, family = "gaussian", standardize = TRUE)
  coef_vec <- as.numeric(coef(final_fit))[-1]
  names(coef_vec) <- features
  pred_train_raw <- as.numeric(predict(final_fit, newx = x_train))
  x_test <- t(as.matrix(test_mat[features, , drop = FALSE]))
  pred_test_raw <- as.numeric(predict(final_fit, newx = x_test))
  calibration <- calibrate_predictions(train_meta$age_weeks, pred_train_raw, pred_test_raw)

  projection_raw <- numeric(0)
  projection_pred <- numeric(0)
  if (!is.null(projection_mat) && ncol(projection_mat) > 0) {
    x_proj <- t(as.matrix(projection_mat[features, , drop = FALSE]))
    projection_raw <- as.numeric(predict(final_fit, newx = x_proj))
    if (is.null(calibration$calibration_fit)) {
      projection_pred <- projection_raw
    } else {
      projection_pred <- as.numeric(predict(calibration$calibration_fit, newdata = data.frame(raw_prediction = projection_raw)))
    }
  }

  list(
    fit = final_fit,
    calibration_fit = calibration$calibration_fit,
    features = features,
    train_raw_prediction = pred_train_raw,
    train_prediction = calibration$train_prediction,
    test_raw_prediction = pred_test_raw,
    test_prediction = calibration$new_prediction,
    projection_raw_prediction = projection_raw,
    projection_prediction = projection_pred,
    calibration_method = calibration$calibration_method,
    best_alpha = best_tuning$alpha,
    best_lambda = best_tuning$lambda,
    cv_mae_weeks = best_tuning$cv_mae_weeks,
    selected_feature_cap = best_tuning$feature_cap,
    cv_grid = bind_rows(tuning_records),
    feature_table = best_tuning$selected %>% mutate(feature_source = "celltype_metacell_training_only"),
    feature_selection_diagnostics = feature_diagnostics %>%
      left_join(feature_tbl %>% select(gene, rho, abs_rho), by = "gene") %>%
      mutate(
        in_top_candidate_set = gene %in% candidates,
        spearman_rank = if_else(is.finite(abs_rho), min_rank(desc(abs_rho)), NA_integer_),
        passes_cor_threshold = is.finite(abs_rho) & abs_rho >= somatic_cor_threshold,
        selected_for_model = gene %in% features,
        selection_mode = best_tuning$selection_mode,
        somatic_top_candidate_genes = somatic_top_candidate_genes,
        somatic_cor_threshold = somatic_cor_threshold,
        selected_feature_cap = best_tuning$feature_cap,
        feature_count_grid = paste(feature_cap_grid, collapse = ","),
        n_total_genes_before_filter = nrow(gene_stats),
        n_prefilter_pass = nrow(feature_diagnostics),
        n_training_metacells = nrow(train_meta),
        n_feature_selection_samples = ncol(sample_means),
        feature_source = "celltype_metacell_training_only"
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
      mutate(best_alpha = best_tuning$alpha, best_lambda = best_tuning$lambda,
             fold_cv_mae_weeks = best_tuning$cv_mae_weeks,
             selected_feature_cap = best_tuning$feature_cap,
             feature_source = "celltype_metacell_training_only")
  )
}

choose_split <- function(candidate_samples) {
  pieces <- list()
  for (bin in names(target_test_by_bin)) {
    bin_df <- candidate_samples %>% filter(reproductive_age_bin == bin) %>% arrange(dataset_id, age_weeks, sample_id)
    test_n <- min(target_test_by_bin[[bin]], nrow(bin_df))
    ds_counts <- bin_df %>% count(dataset_id, name = "n_dataset")
    ds_counts <- ds_counts %>%
      mutate(raw_quota = n_dataset / sum(n_dataset) * test_n,
             quota = floor(raw_quota),
             remainder = raw_quota - quota) %>%
      arrange(desc(remainder), desc(n_dataset), dataset_id)
    while (sum(ds_counts$quota) < test_n) {
      ds_counts$quota[[which.min(ds_counts$quota / pmax(ds_counts$n_dataset, 1))]] <-
        ds_counts$quota[[which.min(ds_counts$quota / pmax(ds_counts$n_dataset, 1))]] + 1L
    }
    selected <- bind_rows(lapply(seq_len(nrow(ds_counts)), function(i) {
      qq <- ds_counts$quota[[i]]
      if (qq <= 0) return(tibble())
      cur <- bin_df %>% filter(dataset_id == ds_counts$dataset_id[[i]]) %>% arrange(age_weeks, sample_id)
      idx <- unique(pmax(1L, pmin(nrow(cur), round(seq(1, nrow(cur), length.out = qq)))))
      cur[idx, , drop = FALSE]
    })) %>%
      distinct(sample_key, .keep_all = TRUE)
    if (nrow(selected) > test_n) selected <- selected %>% arrange(dataset_id, age_weeks, sample_id) %>% slice_head(n = test_n)
    pieces[[bin]] <- bin_df %>%
      mutate(split_set = if_else(sample_key %in% selected$sample_key, "test", "train"))
  }
  bind_rows(pieces) %>% arrange(reproductive_age_bin, split_set, dataset_id, age_weeks, sample_id)
}

metadata <- fread(metadata_path) %>%
  as_tibble() %>%
  mutate(
    dataset_id = standardized_dataset_id(dataset_folder),
    sample_id = as.character(sample_id),
    age_weeks = parse_age_weeks(age, stage),
    reproductive_age_bin = assign_reproductive_age_bin(age_weeks),
    sample_key = paste(dataset_id, sample_id, sep = "__"),
    intervention_class = "natural_aging",
    model_role = "clock_candidate",
    include_clock_candidate = TRUE,
    projection_role = "not_projection"
  )

metadata <- metadata %>%
  mutate(
    intervention_class = case_when(
      dataset_folder == "PRJNA863443" & sample_id %in% allowed_prjna_clock_samples ~ "natural_aging",
      dataset_folder == "PRJNA863443" & grepl("^VCD", sample_id) ~ "vcd_model",
      dataset_folder == "PRJNA863443" & grepl("^CTL", sample_id) ~ "vcd_control",
      dataset_folder == "PRJNA863443" & grepl("^Foxl2_het", sample_id) ~ "foxl2_model",
      dataset_folder == "PRJNA863443" & grepl("^Foxl2_wt", sample_id) ~ "foxl2_control",
      TRUE ~ intervention_class
    ),
    model_role = case_when(
      dataset_folder == "PRJNA863443" & sample_id %in% allowed_prjna_clock_samples ~ "clock_candidate",
      dataset_folder == "PRJNA863443" ~ "projection_only",
      TRUE ~ model_role
    ),
    include_clock_candidate = model_role == "clock_candidate",
    projection_role = case_when(
      dataset_folder == "PRJNA863443" & sample_id %in% allowed_prjna_clock_samples ~ "natural_control_reference",
      dataset_folder == "PRJNA863443" & grepl("^VCD", sample_id) ~ "stress_test_vcd",
      dataset_folder == "PRJNA863443" & grepl("^CTL", sample_id) ~ "stress_test_vcd_control",
      dataset_folder == "PRJNA863443" & grepl("^Foxl2_het", sample_id) ~ "stress_test_foxl2",
      dataset_folder == "PRJNA863443" & grepl("^Foxl2_wt", sample_id) ~ "stress_test_foxl2_control",
      TRUE ~ projection_role
    ),
    reproductive_age_bin = case_when(
      dataset_folder == "PRJNA863443" & sample_id == "AC_13weeks" ~ "young_adult",
      dataset_folder == "PRJNA863443" & sample_id %in% c("YF_1", "YF_2") ~ "mature_adult",
      dataset_folder == "PRJNA863443" & sample_id %in% c("AC_76weeks", "AC_86weeks", "OF_1", "OF_2") ~ "post_reproductive_ageing",
      TRUE ~ reproductive_age_bin
    )
  )

bad_prjna <- metadata %>%
  filter(dataset_folder == "PRJNA863443") %>%
  filter(is.na(intervention_class) | is.na(model_role) | is.na(projection_role))
if (nrow(bad_prjna) > 0) {
  stop("Unclassified PRJNA863443 samples: ", paste(bad_prjna$sample_id, collapse = ", "), call. = FALSE)
}

candidate_samples <- metadata %>%
  filter(model_role == "clock_candidate", reproductive_age_bin != "out_of_scope", !is.na(age_weeks)) %>%
  distinct(dataset_folder, dataset_id, sample_id, sample_key, age_weeks, reproductive_age_bin,
           intervention_class, model_role, include_clock_candidate, projection_role)
split_plan <- choose_split(candidate_samples)
projection_samples <- metadata %>%
  filter(model_role == "projection_only", !is.na(age_weeks)) %>%
  distinct(dataset_folder, dataset_id, sample_id, sample_key, age_weeks, reproductive_age_bin,
           intervention_class, model_role, include_clock_candidate, projection_role)

fwrite(metadata, file.path(metadata_out_dir, "scRNA_sample_role_metadata.tsv"), sep = "\t")
fwrite(split_plan, file.path(metadata_out_dir, "global_train_test_split.tsv"), sep = "\t")
fwrite(
  candidate_samples %>%
    count(reproductive_age_bin, name = "observed_candidate_n") %>%
    mutate(target_n_from_15_original_age_bin_summary = target_n_by_bin[reproductive_age_bin]),
  file.path(metadata_out_dir, "candidate_age_bin_count_diagnostic.tsv"),
  sep = "\t"
)
fwrite(
  split_plan %>% count(reproductive_age_bin, split_set, name = "n_samples"),
  file.path(metadata_out_dir, "split_balance_by_age_bin.tsv"),
  sep = "\t"
)
fwrite(
  split_plan %>% count(reproductive_age_bin, dataset_id, split_set, name = "n_samples"),
  file.path(metadata_out_dir, "split_balance_by_age_bin_dataset.tsv"),
  sep = "\t"
)

manifest <- fread(manifest_path) %>% as_tibble() %>%
  mutate(processed_seurat_path = file.path(V2, "00_input_seurat", basename(processed_seurat_path)))  # v2: copied objects
annotation_consensus <- fread(cmd = paste("zcat", shQuote(consensus_path))) %>%  # v2
  as_tibble() %>%
  transmute(
    cell_id = as.character(cell_id),
    dataset_id = as.character(dataset_id),
    sample_id = as.character(sample_id),
    harmonized_cell_class = as.character(annotation_label_broad),
    annotation_label_fine = as.character(annotation_label_fine),
    curation_status = as.character(curation_status),
    label_assignment_source = as.character(annotation_source),
    label_assignment_confidence = as.numeric(annotation_confidence),
    evidence_notes = as.character(evidence_notes)
  )

sample_roles <- bind_rows(
  split_plan %>% mutate(split_set = split_set),
  projection_samples %>% mutate(split_set = "projection_only")
) %>%
  distinct(sample_key, .keep_all = TRUE)
needed_datasets <- unique(sample_roles$dataset_id)
manifest_needed <- manifest %>% filter(dataset_id %in% needed_datasets)
if (nrow(manifest_needed) != length(needed_datasets)) {
  missing_ds <- setdiff(needed_datasets, manifest_needed$dataset_id)
  stop("Missing standardized Seurat manifest rows for: ", paste(missing_ds, collapse = ", "), call. = FALSE)
}

metacell_count_mats <- list()
metacell_meta_list <- list()
filter_audit_list <- list()
cell_filter_summary_list <- list()

message("Building confident filtered metacells from standardized Seurat objects")
for (i in seq_len(nrow(manifest_needed))) {
  ds <- manifest_needed$dataset_id[[i]]
  rds_path <- manifest_needed$processed_seurat_path[[i]]
  ds_roles <- sample_roles %>% filter(dataset_id == ds)
  if (nrow(ds_roles) == 0) next
  message("Reading ", ds)
  obj <- readRDS(rds_path)
  md <- obj@meta.data %>% rownames_to_column("cell_id") %>% as_tibble()
  if (!"sample_id" %in% colnames(md)) stop("Missing sample_id in Seurat metadata for ", ds, call. = FALSE)
  md <- md %>%
    mutate(sample_id = as.character(sample_id), dataset_id = ds, sample_key = paste(dataset_id, sample_id, sep = "__")) %>%
    filter(sample_key %in% ds_roles$sample_key) %>%
    left_join(ds_roles, by = c("dataset_id", "sample_id", "sample_key"))
  for (joined_name in c(
    "dataset_folder", "age_weeks", "reproductive_age_bin", "intervention_class",
    "model_role", "include_clock_candidate", "projection_role", "split_set"
  )) {
    md <- prefer_joined_column(md, joined_name)
  }
  if (nrow(md) == 0) {
    rm(obj)
    gc()
    next
  }
  counts <- get_assay_matrix(obj)[, md$cell_id, drop = FALSE]
  rm(obj)
  gc()

  md <- md %>%
    left_join(annotation_consensus, by = c("cell_id", "dataset_id", "sample_id"))
  for (joined_name in c(
    "harmonized_cell_class", "annotation_label_fine", "curation_status",
    "label_assignment_source", "label_assignment_confidence", "evidence_notes"
  )) {
    md <- prefer_joined_column(md, joined_name)
  }
  # v2: PRJNA863443 broad_cell_type fallback removed; all datasets are labelled by the v2 consensus
  md <- md %>%
    mutate(
      harmonized_cell_class = if_else(is.na(harmonized_cell_class), "Unknown", harmonized_cell_class),
      curation_status = if_else(is.na(curation_status), "ambiguous", curation_status),
      label_assignment_confidence = if_else(is.na(label_assignment_confidence), 0, label_assignment_confidence),
      confident_celltype_filter_pass = curation_status %in% clock_label_statuses &
        label_assignment_confidence >= min_label_confidence &
        !(harmonized_cell_class %in% clock_excluded_labels),
      model_candidate_filter_pass = model_role == "clock_candidate" &
        !is.na(age_weeks) &
        reproductive_age_bin != "out_of_scope",
      retained_for_metacell = confident_celltype_filter_pass &
        !(exclude_granulosa_luteal & harmonized_cell_class == "Granulosa" & annotation_label_fine %in% "luteal") &  # v2 sensitivity
        (model_candidate_filter_pass | model_role == "projection_only")
    )

  cell_filter_summary_list[[ds]] <- md %>%
    summarise(
      dataset_id = ds,
      n_cells_seen = n(),
      n_confident_celltype_pass = sum(confident_celltype_filter_pass, na.rm = TRUE),
      n_model_candidate_cell_pass = sum(confident_celltype_filter_pass & model_candidate_filter_pass, na.rm = TRUE),
      n_projection_cell_pass = sum(confident_celltype_filter_pass & model_role == "projection_only", na.rm = TRUE),
      n_retained_for_metacell = sum(retained_for_metacell, na.rm = TRUE)
    )

  filter_audit_list[[ds]] <- md %>%
    count(dataset_id, sample_id, sample_key, split_set, model_role, intervention_class, projection_role,
          reproductive_age_bin, harmonized_cell_class, curation_status, confident_celltype_filter_pass,
          model_candidate_filter_pass, retained_for_metacell, name = "n_cells")

  md_keep <- md %>% filter(retained_for_metacell)
  if (nrow(md_keep) == 0) {
    rm(counts)
    gc()
    next
  }
  md_keep$original_col_index <- match(md_keep$cell_id, md$cell_id)
  group <- paste(md_keep$sample_key, md_keep$harmonized_cell_class, sep = "__")
  metacell_cell_index_list <- list()
  metacell_cell_group_list <- list()
  metacell_meta_ds <- list()
  chunk_counter <- 0L
  for (grp in unique(group)) {
    idx_keep <- which(group == grp)
    if (length(idx_keep) < metacell_min_cells) next
    idx_keep <- sample(idx_keep)
    chunk <- ceiling(seq_along(idx_keep) / metacell_target_cells)
    valid_chunks <- as.integer(names(table(chunk))[table(chunk) >= metacell_min_cells])
    if (length(valid_chunks) == 0 && length(idx_keep) >= metacell_min_cells) valid_chunks <- 1L
    for (cc in valid_chunks) {
      local_idx <- idx_keep[chunk == cc]
      row_md <- md_keep[local_idx[[1]], ]
      metacell_id <- paste(row_md$sample_key, row_md$harmonized_cell_class, sprintf("mc%04d", cc), sep = "__")
      original_idx <- md_keep$original_col_index[local_idx]
      chunk_counter <- chunk_counter + 1L
      metacell_cell_index_list[[chunk_counter]] <- original_idx
      metacell_cell_group_list[[chunk_counter]] <- rep(metacell_id, length(original_idx))
      metacell_meta_ds[[metacell_id]] <- tibble(
        metacell_id = metacell_id,
        sample_key = row_md$sample_key,
        dataset_folder = row_md$dataset_folder,
        dataset_id = row_md$dataset_id,
        sample_id = row_md$sample_id,
        split_set = row_md$split_set,
        model_role = row_md$model_role,
        include_clock_candidate = row_md$include_clock_candidate,
        projection_role = row_md$projection_role,
        intervention_class = row_md$intervention_class,
        age_weeks = row_md$age_weeks,
        reproductive_age_bin = row_md$reproductive_age_bin,
        harmonized_cell_class = row_md$harmonized_cell_class,
        label_assignment_source = paste(sort(unique(md_keep$label_assignment_source[local_idx])), collapse = ";"),
        label_assignment_confidence = mean(md_keep$label_assignment_confidence[local_idx], na.rm = TRUE),
        min_label_assignment_confidence = min(md_keep$label_assignment_confidence[local_idx], na.rm = TRUE),
        curation_status_summary = paste(paste0(names(table(md_keep$curation_status[local_idx])),
                                               "=", as.integer(table(md_keep$curation_status[local_idx]))), collapse = ";"),
        n_confirmed_cells = sum(md_keep$curation_status[local_idx] == "confirmed"),
        n_probable_cells = sum(md_keep$curation_status[local_idx] == "probable"),
        n_cells = length(local_idx),
        n_luteal_cells = sum(md_keep$annotation_label_fine[local_idx] == "luteal", na.rm = TRUE),  # v2
        frac_luteal = mean(md_keep$annotation_label_fine[local_idx] == "luteal", na.rm = TRUE),  # v2
        frac_cycling_state = mean(md_keep$annotation_label_fine[local_idx] == "cycling", na.rm = TRUE),  # v2
        frac_confirmed = mean(md_keep$curation_status[local_idx] == "confirmed", na.rm = TRUE)  # v2
      )
    }
  }
  if (length(metacell_meta_ds) > 0) {
    levels <- names(metacell_meta_ds)
    metacell_cell_index <- unlist(metacell_cell_index_list, use.names = FALSE)
    metacell_cell_group <- unlist(metacell_cell_group_list, use.names = FALSE)
    design <- sparseMatrix(
      i = seq_along(metacell_cell_index),
      j = match(metacell_cell_group, levels),
      x = 1,
      dims = c(length(metacell_cell_index), length(levels))
    )
    mat <- counts[, metacell_cell_index, drop = FALSE] %*% design
    colnames(mat) <- levels
    metacell_count_mats[[ds]] <- as(mat, "dgCMatrix")
    metacell_meta_list[[ds]] <- bind_rows(metacell_meta_ds)
  }
  rm(counts)
  gc()
}

message("Combining dataset metacell matrices")
metacell_counts <- combine_sparse_by_rows(metacell_count_mats)
metacell_meta <- bind_rows(metacell_meta_list) %>% arrange(match(metacell_id, colnames(metacell_counts)))
stopifnot(identical(metacell_meta$metacell_id, colnames(metacell_counts)))

message("Writing preprocessing outputs")
fwrite(bind_rows(cell_filter_summary_list), file.path(preprocess_dir, "confident_celltype_filter_summary.tsv"), sep = "\t")
fwrite(bind_rows(filter_audit_list), file.path(preprocess_dir, "confident_celltype_filter_audit.tsv"), sep = "\t")
fwrite(metacell_meta, file.path(preprocess_dir, "filtered_metacell_metadata.tsv"), sep = "\t")
saveRDS(metacell_counts, file.path(preprocess_dir, "filtered_metacell_counts.rds"))

preprocess_readme <- c(
  "# scRNA metacell ElasticNet preprocessing",
  "",
  "Filters applied before model building:",
  "- curation_status in confirmed, probable",
  "- label_assignment_confidence >= 0.70",
  "- exclude Unknown and Erythroid harmonized_cell_class",
  "- metacell_size >= 25 cells",
  "- retain harmonized_cell_class labels only",
  "- model-building candidates require nonmissing age_weeks and reproductive_age_bin != out_of_scope",
  "- PRJNA863443 VCD/Foxl2/control samples are retained only as projection_only rows",
  "",
  "Train-dependent gene filters are fit inside each cell-type model using training samples only."
)
writeLines(preprocess_readme, file.path(preprocess_dir, "README.md"))

# v2: stop after metacell construction; models are fit by run_scRNA_metacell_elasticnet_metacell_cv_v2.R
message("v2 metacell build complete")
quit(save = "no", status = 0)

message("Fitting cell-type-specific metacell ElasticNet clocks")
celltype_eligibility <- metacell_meta %>%
  filter(model_role == "clock_candidate") %>%
  group_by(harmonized_cell_class) %>%
  summarise(
    n_metacells = n(),
    n_samples = n_distinct(sample_key),
    n_train_samples = n_distinct(sample_key[split_set == "train"]),
    n_test_samples = n_distinct(sample_key[split_set == "test"]),
    n_train_age_points = n_distinct(age_weeks[split_set == "train"]),
    n_test_age_points = n_distinct(age_weeks[split_set == "test"]),
    n_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  mutate(
    eligible_for_celltype_metacell_elasticnet = n_train_samples >= min_train_samples &
      n_train_age_points >= min_train_age_points &
      n_metacells >= min_train_metacells &
      n_test_samples > 0,
    ineligibility_reason = case_when(
      eligible_for_celltype_metacell_elasticnet ~ NA_character_,
      n_train_samples < min_train_samples ~ "insufficient_training_samples_after_global_split",
      n_train_age_points < min_train_age_points ~ "insufficient_training_age_points_after_global_split",
      n_metacells < min_train_metacells ~ "insufficient_metacells",
      n_test_samples == 0 ~ "no_test_sample_coverage_after_global_split",
      TRUE ~ "unknown"
    )
  )
fwrite(celltype_eligibility, file.path(preprocess_dir, "celltype_metacell_elasticnet_eligibility.tsv"), sep = "\t")

predictions <- list()
train_predictions <- list()
projection_predictions <- list()
features <- list()
feature_diagnostics <- list()
coefficients <- list()
cv_grids <- list()
model_objects <- list()
gene_filter_rows <- list()

for (ct in celltype_eligibility$harmonized_cell_class[celltype_eligibility$eligible_for_celltype_metacell_elasticnet]) {
  ct_meta <- metacell_meta %>% filter(harmonized_cell_class == ct)
  ct_counts <- metacell_counts[, ct_meta$metacell_id, drop = FALSE]
  ct_log <- logcpm_from_counts(ct_counts)
  rownames(ct_log) <- rownames(ct_counts)
  colnames(ct_log) <- colnames(ct_counts)

  train_meta <- ct_meta %>% filter(model_role == "clock_candidate", split_set == "train")
  test_meta <- ct_meta %>% filter(model_role == "clock_candidate", split_set == "test")
  proj_meta <- ct_meta %>% filter(model_role == "projection_only")
  train_mat <- ct_log[, train_meta$metacell_id, drop = FALSE]
  test_mat <- ct_log[, test_meta$metacell_id, drop = FALSE]
  proj_mat <- if (nrow(proj_meta) > 0) ct_log[, proj_meta$metacell_id, drop = FALSE] else NULL

  gene_filter_rows[[ct]] <- tibble(
    harmonized_cell_class = ct,
    gene = rownames(train_mat),
    train_mean_logcpm = rowMeans(train_mat),
    train_expression_fraction = rowMeans(train_mat > 0),
    train_expression_sd = apply(train_mat, 1, sd),
    n_train_metacells = ncol(train_mat),
    n_train_samples = n_distinct(train_meta$sample_key)
  ) %>%
    mutate(pass_pre_elasticnet_expression_filter = is.finite(train_expression_sd) & train_expression_sd > 0 & train_expression_fraction > 0)

  fit <- fit_direct_elasticnet(train_mat, train_meta, test_mat, test_meta, proj_mat, proj_meta)
  if (is.null(fit)) next

  train_predictions[[ct]] <- train_meta %>%
    transmute(
      split_id = "global_age_bin_train_test",
      fold_id = "train_fit",
      clock_type = "celltype_metacell_elasticnet_clock",
      model_evaluation_role = "training_fit",
      metacell_id, sample_key, sample_id, dataset_folder, dataset_id,
      intervention_class, projection_role, harmonized_cell_class, age_weeks, reproductive_age_bin, n_cells,
      predicted_age_weeks_raw = fit$train_raw_prediction,
      predicted_age_weeks = fit$train_prediction,
      delta_age_weeks = predicted_age_weeks - age_weeks,
      calibration_method = fit$calibration_method,
      n_model_features = length(fit$features),
      best_alpha = fit$best_alpha,
      best_lambda = fit$best_lambda,
      fold_cv_mae_weeks = fit$cv_mae_weeks
    )
  predictions[[ct]] <- test_meta %>%
    transmute(
      split_id = "global_age_bin_train_test",
      fold_id = "test",
      clock_type = "celltype_metacell_elasticnet_clock",
      model_evaluation_role = "heldout_test",
      metacell_id, sample_key, sample_id, dataset_folder, dataset_id,
      intervention_class, projection_role, harmonized_cell_class, age_weeks, reproductive_age_bin, n_cells,
      predicted_age_weeks_raw = fit$test_raw_prediction,
      predicted_age_weeks = fit$test_prediction,
      delta_age_weeks = predicted_age_weeks - age_weeks,
      calibration_method = fit$calibration_method,
      n_model_features = length(fit$features),
      best_alpha = fit$best_alpha,
      best_lambda = fit$best_lambda,
      fold_cv_mae_weeks = fit$cv_mae_weeks
    )
  if (nrow(proj_meta) > 0 && length(fit$projection_prediction) == nrow(proj_meta)) {
    projection_predictions[[ct]] <- proj_meta %>%
      transmute(
        split_id = "projection_only",
        fold_id = "frozen_train_clock",
        clock_type = "celltype_metacell_elasticnet_clock",
        model_evaluation_role = "projection_stress_test",
        metacell_id, sample_key, sample_id, dataset_folder, dataset_id,
        intervention_class, projection_role, harmonized_cell_class, age_weeks, reproductive_age_bin, n_cells,
        predicted_age_weeks_raw = fit$projection_raw_prediction,
        predicted_age_weeks = fit$projection_prediction,
        delta_age_weeks = predicted_age_weeks - age_weeks,
        calibration_method = fit$calibration_method,
        n_model_features = length(fit$features),
        best_alpha = fit$best_alpha,
        best_lambda = fit$best_lambda,
        fold_cv_mae_weeks = fit$cv_mae_weeks
      )
  }
  features[[ct]] <- fit$feature_table %>%
    mutate(split_id = "global_age_bin_train_test", fold_id = "test",
           clock_type = "celltype_metacell_elasticnet_clock", harmonized_cell_class = ct, .before = gene)
  feature_diagnostics[[ct]] <- fit$feature_selection_diagnostics %>%
    mutate(split_id = "global_age_bin_train_test", fold_id = "test",
           clock_type = "celltype_metacell_elasticnet_clock", harmonized_cell_class = ct, .before = gene)
  coefficients[[ct]] <- fit$coefficient_table %>%
    mutate(split_id = "global_age_bin_train_test", fold_id = "test",
           clock_type = "celltype_metacell_elasticnet_clock", harmonized_cell_class = ct, .before = gene)
  cv_grids[[ct]] <- fit$cv_grid %>%
    mutate(split_id = "global_age_bin_train_test", fold_id = "test",
           clock_type = "celltype_metacell_elasticnet_clock", harmonized_cell_class = ct, .before = alpha)
  model_objects[[ct]] <- fit
  rm(ct_log)
  gc()
}

prediction_tbl <- bind_rows(predictions)
train_prediction_tbl <- bind_rows(train_predictions)
projection_tbl <- bind_rows(projection_predictions)
sample_prediction_summary <- prediction_tbl %>%
  group_by(split_id, fold_id, clock_type, model_evaluation_role, sample_key, sample_id, dataset_folder,
           dataset_id, intervention_class, reproductive_age_bin, harmonized_cell_class, age_weeks) %>%
  summarise(
    predicted_age_weeks_raw = weighted.mean(predicted_age_weeks_raw, n_cells, na.rm = TRUE),
    predicted_age_weeks = weighted.mean(predicted_age_weeks, n_cells, na.rm = TRUE),
    delta_age_weeks = predicted_age_weeks - first(age_weeks),
    n_predicted_metacells = n(),
    n_predicted_cells = sum(n_cells),
    .groups = "drop"
  )
train_sample_prediction_summary <- train_prediction_tbl %>%
  group_by(split_id, fold_id, clock_type, model_evaluation_role, sample_key, sample_id, dataset_folder,
           dataset_id, intervention_class, reproductive_age_bin, harmonized_cell_class, age_weeks) %>%
  summarise(
    predicted_age_weeks_raw = weighted.mean(predicted_age_weeks_raw, n_cells, na.rm = TRUE),
    predicted_age_weeks = weighted.mean(predicted_age_weeks, n_cells, na.rm = TRUE),
    delta_age_weeks = predicted_age_weeks - first(age_weeks),
    n_predicted_metacells = n(),
    n_predicted_cells = sum(n_cells),
    .groups = "drop"
  )
training_performance_summary <- performance_table(train_sample_prediction_summary)
evaluation_summary <- sample_prediction_summary %>%
  group_by(split_id, clock_type, harmonized_cell_class) %>%
  group_modify(~metric_table(.x)) %>%
  ungroup()

if (nrow(projection_tbl) > 0) {
  projection_sample_summary <- projection_tbl %>%
    group_by(split_id, fold_id, clock_type, model_evaluation_role, sample_key, sample_id, dataset_folder,
             dataset_id, intervention_class, projection_role, reproductive_age_bin, harmonized_cell_class, age_weeks) %>%
    summarise(
      predicted_age_weeks_raw = weighted.mean(predicted_age_weeks_raw, n_cells, na.rm = TRUE),
      predicted_age_weeks = weighted.mean(predicted_age_weeks, n_cells, na.rm = TRUE),
      delta_age_weeks = predicted_age_weeks - first(age_weeks),
      n_predicted_metacells = n(),
      n_predicted_cells = sum(n_cells),
      .groups = "drop"
    )
} else {
  projection_sample_summary <- tibble()
}

fwrite(prediction_tbl, file.path(model_dir, "celltype_metacell_elasticnet_predictions.tsv"), sep = "\t")
fwrite(train_prediction_tbl, file.path(model_dir, "celltype_metacell_elasticnet_training_fit_predictions.tsv"), sep = "\t")
fwrite(sample_prediction_summary, file.path(model_dir, "celltype_metacell_elasticnet_sample_prediction_summary.tsv"), sep = "\t")
fwrite(train_sample_prediction_summary, file.path(model_dir, "celltype_metacell_elasticnet_training_sample_prediction_summary.tsv"), sep = "\t")
fwrite(training_performance_summary, file.path(model_dir, "celltype_metacell_elasticnet_training_performance_summary.tsv"), sep = "\t")
fwrite(evaluation_summary, file.path(model_dir, "celltype_metacell_elasticnet_evaluation_summary.tsv"), sep = "\t")
fwrite(bind_rows(features), file.path(model_dir, "celltype_metacell_elasticnet_selected_features.tsv"), sep = "\t")
fwrite(bind_rows(feature_diagnostics), file.path(model_dir, "celltype_metacell_elasticnet_feature_selection_diagnostics.tsv"), sep = "\t")
fwrite(bind_rows(coefficients), file.path(model_dir, "celltype_metacell_elasticnet_coefficient_diagnostics.tsv"), sep = "\t")
fwrite(bind_rows(cv_grids), file.path(model_dir, "celltype_metacell_elasticnet_cv_grid.tsv"), sep = "\t")
saveRDS(model_objects, file.path(model_dir, "celltype_metacell_elasticnet_split_models.rds"))
fwrite(bind_rows(gene_filter_rows), file.path(preprocess_dir, "celltype_metacell_elasticnet_training_gene_filter.tsv"), sep = "\t")

fwrite(projection_tbl, file.path(projection_dir, "celltype_metacell_elasticnet_prjna863443_projection_predictions.tsv"), sep = "\t")
fwrite(projection_sample_summary, file.path(projection_dir, "celltype_metacell_elasticnet_prjna863443_projection_sample_summary.tsv"), sep = "\t")

if (nrow(sample_prediction_summary) > 0) {
  p1_data <- format_for_age_bin_plots(sample_prediction_summary)
  p1 <- ggplot(p1_data, aes(age_weeks, predicted_age_weeks, color = dataset_id, shape = reproductive_age_bin)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 1.9, alpha = 0.85) +
    facet_wrap(~ harmonized_cell_class, scales = "free") +
    scale_color_manual(values = dataset_palette, drop = FALSE) +
    scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
    theme_clock_plot(base_size = 9) +
    labs(
      x = "Observed age (weeks)",
      y = "Predicted age (weeks)",
      color = "Dataset",
      shape = "Age bin"
    )
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_predicted_vs_observed.pdf"), p1, width = 11, height = 8)
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_predicted_vs_observed.png"), p1, width = 11, height = 8, dpi = 220)

  p_train_data <- format_for_age_bin_plots(train_sample_prediction_summary)
  train_label_pos <- p_train_data %>%
    group_by(harmonized_cell_class) %>%
    summarise(
      x = min(age_weeks, na.rm = TRUE),
      y = max(predicted_age_weeks, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(training_performance_summary, by = "harmonized_cell_class") %>%
    mutate(performance_label_plot = gsub("; ", "\n", performance_label, fixed = TRUE))
  p_train <- ggplot(p_train_data, aes(age_weeks, predicted_age_weeks, color = dataset_id, shape = reproductive_age_bin)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 1.35, alpha = 0.45) +
    geom_text(
      data = train_label_pos,
      aes(x = x, y = y, label = performance_label_plot),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.2,
      lineheight = 0.9
    ) +
    facet_wrap(~ harmonized_cell_class, scales = "free") +
    scale_color_manual(values = dataset_palette, drop = FALSE) +
    scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
    theme_clock_plot(base_size = 9) +
    labs(
      x = "Observed age (weeks)",
      y = "Training-fit predicted age (weeks)",
      color = "Dataset",
      shape = "Age bin"
    )
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_training_performance_all_samples.pdf"), p_train, width = 11, height = 8)
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_training_performance_all_samples.png"), p_train, width = 11, height = 8, dpi = 220)

  p2 <- ggplot(evaluation_summary, aes(reorder(harmonized_cell_class, mae_weeks), mae_weeks)) +
    geom_col(fill = "#2f6f73") +
    coord_flip() +
    theme_clock_plot(base_size = 9) +
    labs(x = "Cell class", y = "Held-out sample-level MAE (weeks)")
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_mae_by_cell_class.pdf"), p2, width = 7, height = 5)
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_mae_by_cell_class.png"), p2, width = 7, height = 5, dpi = 220)
}

split_plot_data <- split_plan %>%
  mutate(reproductive_age_bin = factor(reproductive_age_bin, levels = age_bin_levels)) %>%
  count(reproductive_age_bin, split_set, name = "n_samples")
fwrite(split_plot_data, file.path(table_dir, "celltype_metacell_elasticnet_split_plot_data.tsv"), sep = "\t")
p3 <- ggplot(split_plot_data, aes(reproductive_age_bin, n_samples, fill = split_set)) +
  geom_col(position = "dodge") +
  scale_x_discrete(labels = age_bin_labels, drop = FALSE) +
  theme_clock_plot(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(x = "Reproductive age bin", y = "Samples", fill = "Split")
ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_split_balance.pdf"), p3, width = 7, height = 4)
ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_split_balance.png"), p3, width = 7, height = 4, dpi = 220)

if (nrow(projection_sample_summary) > 0) {
  p4_data <- format_for_age_bin_plots(projection_sample_summary)
  p4 <- ggplot(p4_data, aes(intervention_class, delta_age_weeks, color = projection_role, shape = reproductive_age_bin)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
    geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.75, size = 1.6) +
    facet_wrap(~ harmonized_cell_class, scales = "free_y") +
    scale_shape_manual(values = age_bin_shapes, labels = age_bin_labels, drop = FALSE) +
    theme_clock_plot(base_size = 9) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
    labs(
      x = "PRJNA863443 projection class",
      y = "Predicted - observed age (weeks)",
      color = "Projection role",
      shape = "Age bin"
    )
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_prjna863443_projection_delta.pdf"), p4, width = 11, height = 8)
  ggsave(file.path(plot_dir, "celltype_metacell_elasticnet_prjna863443_projection_delta.png"), p4, width = 11, height = 8, dpi = 220)
}

release_tables <- c(
  file.path(metadata_out_dir, "scRNA_sample_role_metadata.tsv"),
  file.path(metadata_out_dir, "global_train_test_split.tsv"),
  file.path(metadata_out_dir, "candidate_age_bin_count_diagnostic.tsv"),
  file.path(metadata_out_dir, "split_balance_by_age_bin.tsv"),
  file.path(preprocess_dir, "confident_celltype_filter_summary.tsv"),
  file.path(preprocess_dir, "filtered_metacell_metadata.tsv"),
  file.path(preprocess_dir, "celltype_metacell_elasticnet_eligibility.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_sample_prediction_summary.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_training_sample_prediction_summary.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_training_performance_summary.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_evaluation_summary.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_selected_features.tsv"),
  file.path(model_dir, "celltype_metacell_elasticnet_coefficient_diagnostics.tsv"),
  file.path(projection_dir, "celltype_metacell_elasticnet_prjna863443_projection_sample_summary.tsv")
)
file.copy(release_tables[file.exists(release_tables)], table_dir, overwrite = TRUE)

all_paths <- list.files(rerun_root, recursive = TRUE, full.names = FALSE, all.files = FALSE)
dir_paths <- list.dirs(rerun_root, recursive = TRUE, full.names = FALSE)
check_tbl <- tibble(
  check_id = character(),
  check_name = character(),
  passed = logical(),
  detail = character()
)
add_check <- function(id, name, passed, detail = "") {
  check_tbl <<- bind_rows(check_tbl, tibble(check_id = id, check_name = name, passed = isTRUE(passed), detail = detail))
}

add_check("01", "no_celltype_metacell_elasticnet_filtered_folder",
          !any(grepl("celltype_metacell_elasticnet_filtered", dir_paths, fixed = TRUE)),
          paste(dir_paths[grepl("celltype_metacell_elasticnet_filtered", dir_paths, fixed = TRUE)], collapse = ";"))
add_check("02", "no_pca_glmnet_outputs",
          !any(grepl("pca|glmnet", all_paths, ignore.case = TRUE)),
          paste(head(all_paths[grepl("pca|glmnet", all_paths, ignore.case = TRUE)], 20), collapse = ";"))
add_check("03", "no_weighted_clock_outputs",
          !any(grepl("weighted", all_paths, ignore.case = TRUE)),
          paste(head(all_paths[grepl("weighted", all_paths, ignore.case = TRUE)], 20), collapse = ";"))
add_check("04", "no_sample_pseudobulk_clock_outputs",
          !any(grepl("sample_pseudobulk", all_paths, ignore.case = TRUE)),
          paste(head(all_paths[grepl("sample_pseudobulk", all_paths, ignore.case = TRUE)], 20), collapse = ";"))
train_test_prjna_bad <- split_plan %>%
  filter(dataset_folder == "PRJNA863443", grepl("^VCD|^Foxl2", sample_id))
add_check("05", "no_vcd_or_foxl2_sample_in_train_test",
          nrow(train_test_prjna_bad) == 0,
          paste(train_test_prjna_bad$sample_id, collapse = ";"))
leak <- split_plan %>%
  count(sample_key, name = "n_rows") %>%
  filter(n_rows > 1)
add_check("06", "no_sample_leakage_across_train_test",
          nrow(leak) == 0 && length(intersect(split_plan$sample_key[split_plan$split_set == "train"],
                                              split_plan$sample_key[split_plan$split_set == "test"])) == 0,
          paste(leak$sample_key, collapse = ";"))
age_ref <- fread(age_bin_reference_path) %>%
  as_tibble() %>%
  filter(modality == "scrna_pseudobulk") %>%
  select(reproductive_age_bin, age_range_week, N)
expected_ranges <- tibble(
  reproductive_age_bin = names(target_n_by_bin),
  age_range_week = c("7 to <16", "16 to <40", "40 to <65", ">=65"),
  expected_original_N = as.integer(target_n_by_bin)
)
range_check <- expected_ranges %>%
  left_join(age_ref, by = "reproductive_age_bin") %>%
  mutate(range_matches = age_range_week.x == age_range_week.y)
fwrite(range_check, file.path(check_dir, "age_bin_mapping_reference_check.tsv"), sep = "\t")
add_check("07", "age_bin_mapping_ranges_match_15_original_age_bin_summary",
          all(range_check$range_matches),
          paste(range_check$reproductive_age_bin[!range_check$range_matches], collapse = ";"))
add_check("08", "filtering_documented_in_step_02_preprocessing_filtering",
          file.exists(file.path(preprocess_dir, "README.md")) &&
            file.exists(file.path(preprocess_dir, "confident_celltype_filter_summary.tsv")) &&
            file.exists(file.path(preprocess_dir, "celltype_metacell_elasticnet_training_gene_filter.tsv")),
          preprocess_dir)
release_files <- list.files(release_dir, recursive = TRUE, full.names = FALSE)
add_check("09", "final_tables_plot_data_and_figures_copied_to_step_08",
          length(release_files) > 0 &&
            any(grepl("\\.(pdf|png)$", release_files)) &&
            any(grepl("\\.tsv$", release_files)),
          paste(release_files, collapse = ";"))

fwrite(check_tbl, file.path(check_dir, "hard_check_summary.tsv"), sep = "\t")
writeLines(
  c(
    "# scRNA metacell ElasticNet rerun",
    "",
    paste0("- Output root: ", rerun_root),
    paste0("- Candidate samples: ", nrow(candidate_samples)),
    paste0("- Train samples: ", sum(split_plan$split_set == "train")),
    paste0("- Test samples: ", sum(split_plan$split_set == "test")),
    paste0("- Projection-only samples: ", nrow(projection_samples)),
    paste0("- Eligible cell classes: ", sum(celltype_eligibility$eligible_for_celltype_metacell_elasticnet)),
    paste0("- Hard checks passed: ", sum(check_tbl$passed), "/", nrow(check_tbl)),
    "",
    "Note: `candidate_age_bin_count_diagnostic.tsv` records the candidate-count difference introduced by applying the explicit PRJNA863443 clock-candidate rule to the original 67-sample scRNA age-bin reference."
  ),
  file.path(rerun_root, "README.md")
)

message("Rerun complete: ", rerun_root)
