#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

age_dataset_penalized_defaults <- list(
  n_genes = 2000L,
  dataset_penalty = 1,
  max_dataset_partial_r2 = 0.90,
  balance_direction = TRUE,
  epsilon = 1e-8,
  min_residual_df = 3L,
  min_nonzero_variance_datasets = 2L
)

age_dataset_penalized_tuning_grid <- CJ(
  n_genes = c(1000L, 2000L, 3000L),
  dataset_penalty = c(0.25, 0.5, 1, 2),
  max_dataset_partial_r2 = c(0.75, 0.90, 0.98)
)

safe_partial_r2 <- function(rss_reduced, rss_full) {
  if (!is.finite(rss_reduced) || !is.finite(rss_full) || rss_reduced <= 0) {
    return(NA_real_)
  }
  max(0, min(1, (rss_reduced - rss_full) / rss_reduced))
}

safe_lm_gene_score <- function(expr, meta, age_col, dataset_col, covariates = character(0),
                               min_residual_df = 3L) {
  df <- data.frame(expr_gene = as.numeric(expr), meta, check.names = FALSE)
  keep <- is.finite(df$expr_gene) & is.finite(df[[age_col]]) & !is.na(df[[dataset_col]])
  for (covar in covariates) {
    keep <- keep & !is.na(df[[covar]])
    if (is.numeric(df[[covar]])) {
      keep <- keep & is.finite(df[[covar]])
    }
  }
  df <- df[keep, , drop = FALSE]

  if (nrow(df) <= min_residual_df || stats::var(df$expr_gene) == 0) {
    return(NULL)
  }

  df[[dataset_col]] <- factor(df[[dataset_col]])
  rhs_covars <- c(age_col, dataset_col, covariates)
  full_formula <- stats::as.formula(paste("expr_gene ~", paste(rhs_covars, collapse = " + ")))
  no_age_formula <- stats::as.formula(paste(
    "expr_gene ~",
    paste(c(dataset_col, covariates), collapse = " + ")
  ))
  no_dataset_formula <- stats::as.formula(paste(
    "expr_gene ~",
    paste(c(age_col, covariates), collapse = " + ")
  ))

  full_fit <- tryCatch(stats::lm(full_formula, data = df), error = function(e) NULL)
  no_age_fit <- tryCatch(stats::lm(no_age_formula, data = df), error = function(e) NULL)
  no_dataset_fit <- tryCatch(stats::lm(no_dataset_formula, data = df), error = function(e) NULL)
  if (is.null(full_fit) || is.null(no_age_fit) || is.null(no_dataset_fit)) {
    return(NULL)
  }
  if (stats::df.residual(full_fit) < min_residual_df || anyNA(stats::coef(full_fit))) {
    return(NULL)
  }

  coef_table <- tryCatch(summary(full_fit)$coefficients, error = function(e) NULL)
  if (is.null(coef_table) || !(age_col %in% rownames(coef_table))) {
    return(NULL)
  }
  age_t <- coef_table[age_col, "t value"]
  age_beta <- coef_table[age_col, "Estimate"]
  if (!is.finite(age_t) || !is.finite(age_beta)) {
    return(NULL)
  }

  rss_full <- sum(stats::residuals(full_fit)^2)
  rss_no_age <- sum(stats::residuals(no_age_fit)^2)
  rss_no_dataset <- sum(stats::residuals(no_dataset_fit)^2)

  data.table(
    age_beta = age_beta,
    age_t = age_t,
    partial_r2_age = safe_partial_r2(rss_no_age, rss_full),
    partial_r2_dataset = safe_partial_r2(rss_no_dataset, rss_full),
    n_samples_used = nrow(df),
    residual_df = stats::df.residual(full_fit)
  )
}

count_nonzero_variance_datasets <- function(expr, dataset) {
  by_dataset <- split(as.numeric(expr), dataset)
  sum(vapply(by_dataset, function(x) {
    x <- x[is.finite(x)]
    length(x) > 1 && is.finite(stats::var(x)) && stats::var(x) > 0
  }, logical(1)))
}

select_age_dataset_penalized_genes <- function(logcpm_train, meta_train,
                                               age_col = "age_week",
                                               dataset_col = "dataset_folder",
                                               covariates = character(0),
                                               n_genes = age_dataset_penalized_defaults$n_genes,
                                               dataset_penalty = age_dataset_penalized_defaults$dataset_penalty,
                                               max_dataset_partial_r2 = age_dataset_penalized_defaults$max_dataset_partial_r2,
                                               balance_direction = age_dataset_penalized_defaults$balance_direction,
                                               epsilon = age_dataset_penalized_defaults$epsilon,
                                               min_residual_df = age_dataset_penalized_defaults$min_residual_df,
                                               min_nonzero_variance_datasets = age_dataset_penalized_defaults$min_nonzero_variance_datasets) {
  stopifnot(is.matrix(logcpm_train) || is.data.frame(logcpm_train))
  logcpm_train <- as.matrix(logcpm_train)
  meta_train <- as.data.frame(meta_train, check.names = FALSE)
  stopifnot(age_col %in% colnames(meta_train))
  stopifnot(dataset_col %in% colnames(meta_train))
  stopifnot(all(covariates %in% colnames(meta_train)))

  if (!is.null(colnames(logcpm_train)) && !is.null(rownames(meta_train))) {
    common_samples <- intersect(colnames(logcpm_train), rownames(meta_train))
    logcpm_train <- logcpm_train[, common_samples, drop = FALSE]
    meta_train <- meta_train[common_samples, , drop = FALSE]
  }
  stopifnot(ncol(logcpm_train) == nrow(meta_train))

  scores <- rbindlist(lapply(rownames(logcpm_train), function(gene_id) {
    expr <- logcpm_train[gene_id, ]
    n_var_datasets <- count_nonzero_variance_datasets(expr, meta_train[[dataset_col]])
    if (n_var_datasets < min_nonzero_variance_datasets) {
      return(NULL)
    }
    stat <- safe_lm_gene_score(
      expr = expr,
      meta = meta_train,
      age_col = age_col,
      dataset_col = dataset_col,
      covariates = covariates,
      min_residual_df = min_residual_df
    )
    if (is.null(stat)) {
      return(NULL)
    }
    stat[, gene_id := gene_id]
    stat[, n_nonzero_variance_datasets := n_var_datasets]
    stat
  }), fill = TRUE)

  if (nrow(scores) == 0) {
    stop("No genes passed age/dataset-penalized feature-selection guards.")
  }

  scores[, coverage_weight := n_nonzero_variance_datasets / uniqueN(meta_train[[dataset_col]])]
  scores[, feature_score := (
    abs(age_t) * sqrt(pmax(partial_r2_age, 0)) * coverage_weight
  ) / (
    1 + dataset_penalty * partial_r2_dataset / (partial_r2_age + epsilon)
  )]
  scores[!is.finite(feature_score), feature_score := NA_real_]
  scores[partial_r2_dataset > max_dataset_partial_r2, feature_score := NA_real_]

  eligible <- scores[is.finite(feature_score)]
  if (balance_direction) {
    n_up <- ceiling(n_genes / 2)
    n_down <- floor(n_genes / 2)
    selected <- rbindlist(list(
      head(eligible[age_beta >= 0][order(-feature_score)], n_up),
      head(eligible[age_beta < 0][order(-feature_score)], n_down)
    ), use.names = TRUE, fill = TRUE)
    if (nrow(selected) < n_genes) {
      selected <- unique(rbindlist(list(
        selected,
        head(eligible[order(-feature_score)], n_genes)
      ), use.names = TRUE, fill = TRUE), by = "gene_id")
    }
  } else {
    selected <- head(eligible[order(-feature_score)], n_genes)
  }

  setcolorder(scores, c(
    "gene_id",
    "feature_score",
    "age_beta",
    "age_t",
    "partial_r2_age",
    "partial_r2_dataset",
    "coverage_weight",
    "n_nonzero_variance_datasets",
    "n_samples_used",
    "residual_df"
  ))

  list(
    selected_genes = selected$gene_id,
    selected_gene_table = selected,
    feature_score_table = scores,
    settings = list(
      n_genes = n_genes,
      dataset_penalty = dataset_penalty,
      max_dataset_partial_r2 = max_dataset_partial_r2,
      balance_direction = balance_direction,
      covariates = covariates,
      min_residual_df = min_residual_df,
      min_nonzero_variance_datasets = min_nonzero_variance_datasets
    )
  )
}
