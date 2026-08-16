suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
module_dir <- if (basename(script_dir) == "scripts") normalizePath(file.path(script_dir, ".."), mustWork = TRUE) else normalizePath(script_dir, mustWork = TRUE)
repo_dir <- normalizePath(file.path(module_dir, "..", ".."), mustWork = TRUE)
config <- yaml::read_yaml(file.path(module_dir, "config", "analysis_config.yaml"))
set.seed(config$seed)

path_module <- function(...) file.path(module_dir, ...)
path_repo <- function(...) file.path(repo_dir, ...)
ensure_dir <- function(path) dir.create(path, recursive = TRUE, showWarnings = FALSE)
write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  data.table::fwrite(as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE)
  invisible(path)
}
read_tsv <- function(path, ...) data.table::fread(path, sep = "\t", ...)
strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))
safe_num <- function(x) suppressWarnings(as.numeric(x))
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

sha256_one <- function(path) {
  unname(tools::md5sum(path)) # supplemented by sha256sum in inventory script
}

read_expression_matrix <- function(path) {
  dt <- fread(path)
  gene_col <- names(dt)[1]
  genes <- trimws(as.character(dt[[gene_col]]))
  dt[, (gene_col) := NULL]
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- genes
  keep <- !is.na(genes) & genes != ""
  mat[keep, , drop = FALSE]
}

collapse_rows_mean <- function(mat) {
  ids <- rownames(mat)
  if (!anyDuplicated(ids)) return(mat)
  groups <- split(seq_along(ids), ids)
  out <- t(vapply(groups, function(i) colMeans(mat[i, , drop = FALSE], na.rm = TRUE), numeric(ncol(mat))))
  colnames(out) <- colnames(mat)
  out
}

load_mouse_reference <- function(artifact) {
  meta <- as.data.table(artifact$reference_metadata)
  align <- fread(path_repo(config$mouse$alignment_file))
  align <- align[sample_key %in% artifact$training_sample_keys]
  if (nrow(align) != length(artifact$training_sample_keys)) {
    stop("Frozen alignment does not contain exactly one row per mouse reference sample.")
  }
  pieces <- list()
  for (p in unique(meta$expression_logcpm_path)) {
    rows <- meta[expression_logcpm_path == p]
    a <- align[sample_key %in% rows$sample_key]
    mat <- collapse_rows_mean(read_expression_matrix(p))
    idx <- match(rows$sample_key, a$sample_key)
    cols <- a$expression_column[idx]
    if (anyNA(cols) || any(!cols %in% colnames(mat))) {
      stop("Cannot reconstruct mouse reference expression for ", basename(p))
    }
    z <- mat[, cols, drop = FALSE]
    colnames(z) <- rows$sample_key
    pieces[[p]] <- z
  }
  common <- Reduce(intersect, lapply(pieces, rownames))
  expr <- do.call(cbind, lapply(pieces, function(x) x[common, , drop = FALSE]))
  expr <- expr[, artifact$training_sample_keys, drop = FALSE]
  meta <- meta[match(artifact$training_sample_keys, sample_key)]
  stopifnot(identical(colnames(expr), meta$sample_key))
  list(expression = expr, metadata = meta)
}

scale_with_artifact <- function(expr, artifact, genes = artifact$selected_genes) {
  x <- t(expr[genes, artifact$training_sample_keys, drop = FALSE])
  sc <- artifact$gene_scale[genes]
  sc[!is.finite(sc) | sc == 0] <- 1
  x <- sweep(x, 2, artifact$gene_center[genes], "-")
  x <- sweep(x, 2, sc, "/")
  x[!is.finite(x)] <- 0
  x
}

predict_artifact <- function(artifact, x_scaled, calibrate = TRUE) {
  if (artifact$model_family == "elasticnet") {
    raw <- as.numeric(predict(artifact$model$glmnet_fit, newx = x_scaled))
  } else {
    pcs <- predict(artifact$model$pca_fit, newdata = x_scaled)
    pcs <- pcs[, artifact$model$selected_pc_names, drop = FALSE]
    raw <- as.numeric(predict(artifact$model$glmnet_fit, newx = pcs))
  }
  if (!calibrate) return(raw)
  artifact$calibration$intercept + artifact$calibration$slope * raw
}

model_gene_weights <- function(artifact) {
  genes <- artifact$selected_genes
  if (artifact$model_family == "elasticnet") {
    tab <- as.data.table(artifact$model$coefficient_table)
    b <- setNames(tab$coefficient, tab$gene_id)[genes]
  } else {
    pc_b <- as.numeric(coef(artifact$model$glmnet_fit))[-1]
    names(pc_b) <- rownames(coef(artifact$model$glmnet_fit))[-1]
    pc_b <- pc_b[artifact$model$selected_pc_names]
    rot <- artifact$model$pca_fit$rotation[genes, artifact$model$selected_pc_names, drop = FALSE]
    b <- as.numeric(rot %*% pc_b)
    names(b) <- genes
  }
  sc <- artifact$gene_scale[genes]
  sc[!is.finite(sc) | sc == 0] <- 1
  artifact$calibration$slope * b / sc
}

weighted_cor <- function(x, y, w) {
  ok <- is.finite(x) & is.finite(y) & is.finite(w) & w > 0
  x <- rank(x[ok], ties.method = "average")
  y <- rank(y[ok], ties.method = "average")
  w <- w[ok] / sum(w[ok])
  mx <- sum(w * x); my <- sum(w * y)
  sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
}

hc3_lm <- function(formula, data) {
  fit <- lm(formula, data = data)
  X <- model.matrix(fit)
  e <- residuals(fit)
  h <- hatvalues(fit)
  meat <- crossprod(X, X * as.numeric((e / pmax(1 - h, 1e-8))^2))
  bread <- solve(crossprod(X))
  vc <- bread %*% meat %*% bread
  se <- sqrt(diag(vc))
  b <- coef(fit)
  t <- b / se
  p <- 2 * pt(abs(t), df = df.residual(fit), lower.tail = FALSE)
  data.table(term = names(b), estimate = b, std_error_hc3 = se,
             conf_low = b - qt(.975, df.residual(fit)) * se,
             conf_high = b + qt(.975, df.residual(fit)) * se,
             t = t, p = p, n = nobs(fit))
}

age_midpoint <- function(x) {
  z <- gsub("[^0-9-]", "", as.character(x))
  vapply(strsplit(z, "-", fixed = TRUE), function(a) mean(as.numeric(a), na.rm = TRUE), numeric(1))
}

make_data_dictionary <- function(path, descriptions = NULL) {
  x <- fread(path, nrows = 1)
  d <- data.table(column = names(x), description = names(x))
  if (!is.null(descriptions)) d[match(names(descriptions), column), description := unname(descriptions)]
  write_tsv(d, paste0(path, ".data_dictionary.tsv"))
}

init_dirs <- function() {
  dirs <- c("data_manifest", "derived_data/mouse_models", "derived_data/orthology",
            "derived_data/gtex/raw", "derived_data/hpa/raw", "derived_data/stageA",
            "derived_data/stageB", "results/tables", "results/qc", "results/stageA",
            "results/stageB", "results/null_tests", "1_publication/panels",
            "1_publication/assembled/human_cross_species_main", "1_publication/source_data",
            "1_publication/captions", "reports", "logs")
  invisible(lapply(path_module(dirs), ensure_dir))
}

init_dirs()
