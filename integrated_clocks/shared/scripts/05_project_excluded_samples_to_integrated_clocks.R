#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(glmnet)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_input_dir <- normalizePath(file.path(project_dir, "..", "data_input"), mustWork = TRUE)
results_dir <- Sys.getenv("RESULTS_ROOT", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = TRUE)
if (grepl("results_groupmapped_v2", results_dir, fixed = TRUE) || grepl("_v2", basename(results_dir), fixed = TRUE)) {
  stop("Refusing to run v1 external projection against a v2 results root: ", results_dir, call. = FALSE)
}
metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
out_dir <- file.path(results_dir, "external_projection")
plot_dir <- file.path(out_dir, "plots")

archive_existing_external_projection <- function() {
  overwrite_requested <- identical(Sys.getenv("OVERWRITE_EXTERNAL_PROJECTION", unset = "FALSE"), "TRUE")
  overwrite_allowed <- identical(Sys.getenv("ALLOW_EXTERNAL_PROJECTION_OVERWRITE", unset = "FALSE"), "TRUE")
  if (!dir.exists(out_dir) || length(list.files(out_dir, all.files = TRUE, no.. = TRUE)) == 0) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    return(invisible(NA_character_))
  }
  if (!overwrite_requested || !overwrite_allowed) {
    stop(
      "External projection output already exists. Set OVERWRITE_EXTERNAL_PROJECTION=TRUE ",
      "and ALLOW_EXTERNAL_PROJECTION_OVERWRITE=TRUE to archive and overwrite.",
      call. = FALSE
    )
  }
  archive_root <- file.path(results_dir, "external_projection_archive")
  dir.create(archive_root, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  archive_dir <- file.path(archive_root, paste0("external_projection_before_reference_control_patch_", stamp))
  if (!file.rename(out_dir, archive_dir)) {
    stop("Could not archive existing external projection directory to: ", archive_dir, call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old_alignment <- file.path(archive_dir, "alignment_audit")
  if (dir.exists(old_alignment)) {
    file.copy(old_alignment, out_dir, recursive = TRUE)
  }
  message("Archived previous external projection to: ", archive_dir)
  invisible(archive_dir)
}

projection_archive_dir <- archive_existing_external_projection()
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

missing_artifact_path <- file.path(out_dir, "projection_missing_artifacts.txt")
manual_projection_mapping_path <- file.path(
  out_dir,
  "alignment_audit",
  "manual_projection_sample_mapping_APPROVED_AUTO.tsv"
)

retained_clocks <- data.table(
  model_family = c("elasticnet", "pca_glmnet"),
  selector = c("hvg_baseline", "version_b_age_dataset_penalized"),
  artifact_dir = file.path(
    results_dir,
    "model_artifacts",
    "deployment"
  ),
  model_path = file.path(
    results_dir,
    "model_artifacts",
    "deployment",
    c(
      "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds",
      "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
    )
  ),
  selected_genes_path = file.path(
    results_dir,
    "model_artifacts",
    "deployment",
    c(
      "integrated_bulk_pseudobulk__elasticnet__hvg_baseline__deployment_all_control_like.rds",
      "integrated_bulk_pseudobulk__pca_glmnet__version_b_age_dataset_penalized__deployment_all_control_like.rds"
    )
  )
)

required_external_projection_datasets <- c(
  "CRA003645", "CRA017492", "GSE194196", "GSE252791", "GSE273724",
  "GSE279446", "GSE279573", "GSE329024", "GSE228975", "GSE240547",
  "PRJNA1076509", "GSE235021", "GSE154890", "GSE249151", "GSE106339",
  "GSE128240", "GSE98497", "PRJNA787028", "PRJNA1002222"
)

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

gene_id_symbol_map_cache <- NULL

load_gene_id_symbol_map <- function() {
  if (!is.null(gene_id_symbol_map_cache)) return(gene_id_symbol_map_cache)

  feature_paths <- unique(c(
    Sys.glob(file.path(
      data_input_dir, "standardized_seurat", "01_count_matrices",
      "*", "*", "filtered_feature_bc_matrix", "features.tsv.gz"
    )),
    Sys.glob(file.path(
      data_input_dir, "scRNA_dataset", "*", "cellranger_counts",
      "*", "outs", "filtered_feature_bc_matrix", "features.tsv.gz"
    ))
  ))
  feature_paths <- feature_paths[file.exists(feature_paths)]

  maps <- list()
  for (feature_path in feature_paths) {
    dt <- tryCatch(
      fread(cmd = paste("gzip -cd", shQuote(feature_path)), header = FALSE),
      error = function(e) NULL
    )
    if (is.null(dt) || ncol(dt) < 2) next
    maps[[length(maps) + 1L]] <- dt[, .(
      gene_id = na_chr(V1),
      gene_symbol = na_chr(V2)
    )]
  }

  if (length(maps) == 0) {
    gene_id_symbol_map_cache <<- data.table(gene_id = character(), gene_symbol = character())
    return(gene_id_symbol_map_cache)
  }

  mapping <- rbindlist(maps, use.names = TRUE, fill = TRUE)
  mapping <- mapping[
    grepl("^ENSMUSG", gene_id) &
      !is.na(gene_symbol) &
      gene_symbol != "" &
      !grepl("^ENSMUSG", gene_symbol)
  ]
  mapping[, gene_id_base := sub("\\.[0-9]+$", "", gene_id)]
  mapping <- mapping[, .N, by = .(gene_id_base, gene_symbol)]
  setorder(mapping, gene_id_base, -N, gene_symbol)
  mapping <- mapping[, .SD[1], by = gene_id_base][, .(gene_id_base, gene_symbol)]
  gene_id_symbol_map_cache <<- mapping
  mapping
}

map_ensembl_rows_to_symbols <- function(mat, path) {
  gene_ids <- rownames(mat)
  if (!any(grepl("^ENSMUSG", gene_ids))) return(mat)

  mapping <- load_gene_id_symbol_map()
  if (nrow(mapping) == 0) return(mat)

  symbol_lookup <- mapping$gene_symbol
  names(symbol_lookup) <- mapping$gene_id_base
  gene_id_base <- sub("\\.[0-9]+$", "", gene_ids)
  mapped_symbols <- unname(symbol_lookup[gene_id_base])
  mapped <- !is.na(mapped_symbols) & mapped_symbols != ""
  if (!any(mapped)) return(mat)

  new_gene_ids <- gene_ids
  new_gene_ids[mapped] <- mapped_symbols[mapped]
  if (!anyDuplicated(new_gene_ids)) {
    rownames(mat) <- new_gene_ids
    message(
      "Mapped ", sum(mapped), " Ensembl gene rows to symbols for ",
      basename(path), "."
    )
    return(mat)
  }

  groups <- split(seq_along(new_gene_ids), new_gene_ids)
  collapsed <- t(vapply(
    groups,
    function(idx) colMeans(mat[idx, , drop = FALSE], na.rm = TRUE),
    numeric(ncol(mat))
  ))
  colnames(collapsed) <- colnames(mat)
  message(
    "Mapped ", sum(mapped), " Ensembl gene rows to symbols and collapsed ",
    length(new_gene_ids) - length(groups), " duplicate rows for ", basename(path), "."
  )
  collapsed
}

detect_col <- function(dt, candidates, context, required = TRUE) {
  nms <- names(dt)
  direct <- candidates[candidates %in% nms]
  if (length(direct) > 0) return(direct[[1]])
  lower_hit <- nms[tolower(nms) %in% tolower(candidates)]
  if (length(lower_hit) > 0) return(lower_hit[[1]])
  regex <- paste(candidates, collapse = "|")
  fuzzy <- nms[grepl(regex, nms, ignore.case = TRUE)]
  if (length(fuzzy) > 0) return(fuzzy[[1]])
  if (required) stop("Could not detect ", context, " column. Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  NA_character_
}

read_expression_matrix <- function(path) {
  if (!file.exists(path)) stop("Missing expression matrix: ", path, call. = FALSE)
  dt <- fread(path)
  gene_col <- names(dt)[1]
  genes <- na_chr(dt[[gene_col]])
  dt[, (gene_col) := NULL]
  mat <- as.matrix(dt)
  mode(mat) <- "numeric"
  rownames(mat) <- genes
  keep <- !is.na(rownames(mat)) & rownames(mat) != "" & !duplicated(rownames(mat))
  mat <- mat[keep, , drop = FALSE]
  map_ensembl_rows_to_symbols(mat, path)
}

load_manual_projection_mapping <- function() {
  if (!file.exists(manual_projection_mapping_path)) return(data.table())
  mapping <- fread(manual_projection_mapping_path, sep = "\t")
  required <- c(
    "metadata_sample_key", "candidate_expression_column",
    "use_for_projection_TRUE_FALSE", "match_status"
  )
  missing <- setdiff(required, names(mapping))
  if (length(missing) > 0) {
    stop(
      "Manual projection mapping is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  mapping[, use_for_projection_TRUE_FALSE := as.logical(use_for_projection_TRUE_FALSE)]
  mapping <- mapping[use_for_projection_TRUE_FALSE == TRUE]
  mapping[grepl("^safe_", match_status)]
}

manual_projection_mapping <- load_manual_projection_mapping()

align_dataset_matrix <- function(mat, meta_ds, modality) {
  expr_cols <- colnames(mat)
  candidate_values <- list()
  for (candidate_name in c("sample_id", "accession", "sample_key")) {
    if (candidate_name %in% names(meta_ds)) candidate_values[[candidate_name]] <- meta_ds[[candidate_name]]
  }
  if (modality != "bulk_rnaseq") {
    candidate_values$accession <- NULL
  }

  matched_col <- rep(NA_character_, nrow(meta_ds))
  matched_by <- rep(NA_character_, nrow(meta_ds))
  replicate_cols <- vector("list", nrow(meta_ds))

  if (nrow(manual_projection_mapping) > 0) {
    map_ds <- manual_projection_mapping[metadata_sample_key %chin% meta_ds$sample_key]
    if (nrow(map_ds) > 0) {
      for (i in seq_len(nrow(meta_ds))) {
        map_row <- map_ds[metadata_sample_key == meta_ds$sample_key[i]]
        if (nrow(map_row) == 0) next
        manual_cols <- unique(unlist(strsplit(as.character(map_row$candidate_expression_column[[1]]), ";", fixed = TRUE)))
        manual_cols <- trimws(manual_cols)
        manual_cols <- manual_cols[manual_cols %chin% expr_cols]
        expected_cols <- unique(unlist(strsplit(as.character(map_row$candidate_expression_column[[1]]), ";", fixed = TRUE)))
        expected_cols <- trimws(expected_cols)
        missing_cols <- setdiff(expected_cols, manual_cols)
        if (length(missing_cols) > 0) {
          stop(
            "Manual projection mapping for ", meta_ds$sample_key[i],
            " references expression columns not present in matrix: ",
            paste(missing_cols, collapse = ", "),
            call. = FALSE
          )
        }
        if (length(manual_cols) == 0) next
        replicate_cols[[i]] <- manual_cols
        matched_col[i] <- paste(manual_cols, collapse = ";")
        matched_by[i] <- paste0("manual_", map_row$match_tier[[1]])
      }
    }
  }

  for (candidate_name in names(candidate_values)) {
    idx <- match(as.character(candidate_values[[candidate_name]]), expr_cols)
    use <- is.na(matched_col) & !is.na(idx)
    matched_col[use] <- expr_cols[idx[use]]
    matched_by[use] <- candidate_name
  }

  norm_token <- function(x) {
    x <- gsub("[^a-z0-9]+", "", tolower(na_chr(x)))
    x[!is.na(x) & nchar(x) >= 3]
  }
  norm_keep <- function(x) gsub("[^a-z0-9]+", "", tolower(na_chr(x)))
  replicate_hits_for_row <- function(row) {
    txt <- c(
      row$sample_id %||% NA_character_,
      row$accession %||% NA_character_,
      row$treatment %||% NA_character_,
      row$sample_title %||% NA_character_
    )
    tokens <- unique(unlist(strsplit(tolower(paste(txt, collapse = " ")), "[^a-z0-9]+")))
    tokens <- tokens[!is.na(tokens) & nchar(tokens) >= 3]
    if (any(tokens %chin% c("control", "ctrl"))) tokens <- unique(c(tokens, "control", "ctrl"))
    tokens <- tokens[!tokens %chin% c("diet", "sample", "rna", "seq", "gse")]
    expr_norm <- norm_keep(expr_cols)
    names(expr_norm) <- expr_cols
    expr_norm <- expr_norm[!is.na(expr_norm) & nchar(expr_norm) >= 3]
    hit_cols <- character()
    for (token in tokens) {
      token_norm <- norm_token(token)
      hit_cols <- c(hit_cols, names(expr_norm)[grepl(token_norm, expr_norm, fixed = TRUE)])
    }
    unique(hit_cols)
  }

  if (any(is.na(matched_col))) {
    for (i in which(is.na(matched_col))) {
      hits <- replicate_hits_for_row(meta_ds[i])
      if (length(hits) > 0) {
        replicate_cols[[i]] <- hits
        matched_col[i] <- paste(hits, collapse = ";")
        matched_by[i] <- "replicate_pattern_mean"
      }
    }
  }

  if (any(is.na(matched_col))) {
    if (ncol(mat) == nrow(meta_ds)) {
      missing_idx <- which(is.na(matched_col))
      matched_col[missing_idx] <- expr_cols[missing_idx]
      matched_by[missing_idx] <- "dataset_position"
    } else {
      missing_samples <- meta_ds[is.na(matched_col), sample_key]
      stop(
        "Could not align all samples for ", unique(meta_ds$dataset_folder),
        ". Missing: ", paste(missing_samples, collapse = ", "),
        call. = FALSE
      )
    }
  }

  out <- matrix(NA_real_, nrow = nrow(mat), ncol = nrow(meta_ds))
  rownames(out) <- rownames(mat)
  for (i in seq_len(nrow(meta_ds))) {
    if (length(replicate_cols[[i]]) > 0) {
      out[, i] <- rowMeans(mat[, replicate_cols[[i]], drop = FALSE], na.rm = TRUE)
    } else {
      out[, i] <- mat[, matched_col[[i]], drop = TRUE]
    }
  }
  colnames(out) <- meta_ds$sample_key
  list(
    matrix = out,
    alignment = data.table(
      sample_key = meta_ds$sample_key,
      dataset_folder = meta_ds$dataset_folder,
      expression_column = matched_col,
      matched_by = matched_by
    )
  )
}

load_expression_for_keys <- function(meta, keys) {
  meta_use <- copy(meta[sample_key %chin% keys])
  if (nrow(meta_use) == 0) stop("No metadata rows selected for expression loading.", call. = FALSE)
  pieces <- list()
  alignments <- list()
  alignment_issues <- list()
  for (mod in unique(meta_use$modality)) {
    meta_mod <- meta_use[modality == mod]
    for (ds in unique(meta_mod$dataset_folder)) {
      meta_ds <- meta_mod[dataset_folder == ds]
      paths <- unique(na_chr(meta_ds$expression_logcpm_path))
      paths <- paths[!is.na(paths)]
      if (length(paths) != 1) {
        alignment_issues[[length(alignment_issues) + 1L]] <- data.table(
          modality = mod,
          dataset_folder = ds,
          sample_key = meta_ds$sample_key,
          issue = paste0("Expected one expression_logcpm_path; found ", length(paths))
        )
        next
      }
      aligned <- tryCatch(
        {
          mat <- read_expression_matrix(paths[[1]])
          align_dataset_matrix(mat, meta_ds, mod)
        },
        error = function(e) {
          alignment_issues[[length(alignment_issues) + 1L]] <<- data.table(
            modality = mod,
            dataset_folder = ds,
            sample_key = meta_ds$sample_key,
            issue = conditionMessage(e)
          )
          NULL
        }
      )
      if (is.null(aligned)) next
      pieces[[paste(mod, ds, sep = "\r")]] <- aligned$matrix
      alignments[[paste(mod, ds, sep = "\r")]] <- aligned$alignment
    }
  }
  if (length(pieces) == 0) stop("No expression matrices could be aligned for projection.", call. = FALSE)
  list(
    expression_by_piece = pieces,
    alignment = rbindlist(alignments, fill = TRUE),
    alignment_issues = if (length(alignment_issues) > 0) rbindlist(alignment_issues, fill = TRUE) else data.table(),
    aligned_sample_keys = unique(unlist(lapply(pieces, colnames), use.names = FALSE))
  )
}

extract_calibration <- function(obj) {
  candidates <- c("calibration", "calibration_model", "calibration_fit", "calibration_params", "calibration_parameters")
  hit <- candidates[candidates %in% names(obj)]
  if (length(hit) > 0) {
    cal <- obj[[hit[[1]]]]
    if (inherits(cal, "lm")) return(list(type = "lm", fit = cal))
    if (is.list(cal)) {
      intercept <- cal$intercept %||% cal$calibration_intercept %||% cal$a
      slope <- cal$slope %||% cal$calibration_slope %||% cal$b
      if (length(intercept) == 1 && length(slope) == 1 && is.finite(intercept) && is.finite(slope)) {
        return(list(type = "linear", intercept = as.numeric(intercept), slope = as.numeric(slope)))
      }
    }
    if (is.numeric(cal) && all(c("intercept", "slope") %in% names(cal))) {
      return(list(type = "linear", intercept = unname(cal[["intercept"]]), slope = unname(cal[["slope"]])))
    }
  }
  intercept <- obj$calibration_intercept %||% obj$intercept
  slope <- obj$calibration_slope %||% obj$slope
  if (length(intercept) == 1 && length(slope) == 1 && is.finite(intercept) && is.finite(slope)) {
    return(list(type = "linear", intercept = as.numeric(intercept), slope = as.numeric(slope)))
  }
  NULL
}

apply_calibration <- function(raw_pred, cal) {
  if (is.null(cal)) stop("Cannot apply missing calibration.", call. = FALSE)
  if (identical(cal$type, "lm")) {
    return(as.numeric(predict(cal$fit, newdata = data.frame(raw_train = raw_pred))))
  }
  if (identical(cal$type, "linear")) {
    return(cal$intercept + cal$slope * raw_pred)
  }
  stop("Unsupported calibration object type: ", cal$type, call. = FALSE)
}

normalize_artifact_object <- function(obj, model_family) {
  if (!is.null(obj$schema_version) && obj$schema_version == "integrated_ovary_clock_artifact_v1") {
    out <- list(
      model = obj$model$glmnet_fit,
      genes = obj$selected_genes,
      center = obj$gene_center,
      scale = obj$gene_scale,
      training_scope = obj$training_scope,
      calibration = obj$calibration,
      source_artifact = obj
    )
    if (model_family == "pca_glmnet") {
      out$pca <- obj$model$pca_fit
      out$pc_count <- obj$model$selected_pc_count
      out$pc_coefficient_table <- obj$model$pc_coefficient_table
      out$calibrated_intercept <- obj$model$calibrated_intercept
    }
    return(out)
  }
  obj
}

audit_artifacts <- function(clock_dt) {
  rows <- list()
  loaded <- list()
  for (i in seq_len(nrow(clock_dt))) {
    family <- clock_dt$model_family[i]
    selector <- clock_dt$selector[i]
    model_path <- clock_dt$model_path[i]
    selected_genes_path <- clock_dt$selected_genes_path[i]

    add <- function(requirement, present, detail) {
      rows[[length(rows) + 1L]] <<- data.table(
        model_family = family,
        selector = selector,
        requirement = requirement,
        present = present,
        detail = detail
      )
    }

    add("locked/deployment model RDS", file.exists(model_path), model_path)
    if (!file.exists(model_path)) next

    raw_obj <- readRDS(model_path)
    obj <- normalize_artifact_object(raw_obj, family)
    cal <- extract_calibration(obj)
    loaded[[paste(family, selector, sep = "\r")]] <- obj
    add("selected genes", !is.null(obj$genes) && length(obj$genes) > 0, "selected genes in deployment model RDS")
    add("training/reference metadata", !is.null(obj$training_scope) && obj$training_scope == "integrated_bulk_pseudobulk", "training_scope in model RDS")
    add("artifact scope is deployment", !is.null(raw_obj$artifact_scope) && grepl("deployment", raw_obj$artifact_scope), "artifact_scope in model RDS")
    add("gene center/scale or preprocessing object", all(c("center", "scale") %in% names(obj)) && length(obj$center) > 0 && length(obj$scale) > 0, "center and scale vectors in model RDS")
    add("glmnet model or coefficient table", "model" %in% names(obj) && !is.null(obj$model), "model object in deployment model RDS")
    add("calibration intercept", !is.null(cal) && !is.null(cal$intercept) && is.finite(cal$intercept), "calibration$intercept in model RDS")
    add("calibration slope", !is.null(cal) && !is.null(cal$slope) && is.finite(cal$slope), "calibration$slope in model RDS")

    if (family == "pca_glmnet") {
      add("PCA rotation/loadings", "pca" %in% names(obj) && "rotation" %in% names(obj$pca) && nrow(obj$pca$rotation) > 0, "pca$rotation in model RDS")
      add("selected PC list", "pc_count" %in% names(obj) && is.finite(obj$pc_count) && obj$pc_count > 0, "pc_count in model RDS")
    }
  }
  list(audit = rbindlist(rows, fill = TRUE), objects = loaded)
}

write_missing_and_stop <- function(audit_dt) {
  missing <- audit_dt[present != TRUE]
  if (nrow(missing) == 0) return(invisible(FALSE))
  lines <- c(
    "External projection stopped because required locked/deployment artifacts are missing.",
    "",
    "Missing artifacts/objects:",
    paste0(
      "- ", missing$model_family, " / ", missing$selector, ": ",
      missing$requirement, " (", missing$detail, ")"
    ),
    "",
    "No approximate projection was attempted."
  )
  writeLines(lines, missing_artifact_path)
  stop(paste(lines, collapse = "\n"), call. = FALSE)
}

projection_class_one <- function(row) {
  txt <- tolower(paste(
    row$exclusion_reason %||% "",
    row$sample_context %||% "",
    row$treatment %||% "",
    row$genotype %||% "",
    row$sample_title %||% "",
    row$dev_stage %||% "",
    row$bio_age_bin %||% "",
    collapse = " "
  ))
  age <- safe_num(row$age_week)
  if (grepl("postpartum|lactat", txt)) return("postpartum_lactation")
  if (grepl("microbi|fmt|fecal", txt)) return("microbiome")
  if (grepl("diet|high-fat|hfhs|sucrose|hfd", txt)) return("diet")
  if (grepl("stress|cums", txt)) return("stress")
  if (grepl("ko|knock|mutant|cko|tr/tr|heterozygous|zeb1|abcg5|bnc1|srsf1|rnf111", txt)) return("genetic_model")
  if (grepl("developmental|prepubertal|neonatal|embryo|juvenile", txt) || (!is.na(age) && age < 7)) return("developmental_or_prepubertal")
  if (grepl("\\b(nmn|hep14|rescue|vehicle|treat(ed|ment)?|intervention|drug)\\b", txt)) return("intervention")
  if (!is.na(age) && (age < 7 || age > 110)) return("age_range_OOD")
  "other_excluded"
}

row_text <- function(dt) {
  cols <- intersect(
    c("sample_key", "sample_id", "accession", "treatment", "genotype", "sample_context", "sample_title"),
    names(dt)
  )
  if (length(cols) == 0) return(rep("", nrow(dt)))
  do.call(paste, c(dt[, ..cols], sep = " "))
}

gse279446_cohort_one <- function(txt) {
  txt <- tolower(txt)
  if (grepl("youngcohort|young cohort", txt)) return("young_cohort")
  "adult_cohort"
}

gse279446_platform_one <- function(cohort) {
  if (identical(cohort, "young_cohort")) "GPL24247_NovaSeq6000" else "GPL13112_HiSeq2000"
}

reference_control_group_label_one <- function(group_id) {
  if (is.na(group_id) || !nzchar(group_id)) return(NA_character_)
	labels <- c(
	  GSE194196_WT_1_4_7 = "GSE194196 WT-1/WT-4/WT-7 controls",
	  GSE252791_WT_1_4_7 = "GSE252791 WT1/WT2/WT3 controls",
	  GSE279446_adult_WT_ovary_controls = "GSE279446 adult WT ovary controls",
	  GSE279446_young_WT_controls = "GSE279446 young WT ovary controls",
	  GSE279573_WT_chow_fertile_controls = "GSE279573 WT chow fertile controls",
	  GSE329024_2M_8W_control_meta = "GSE329024 2M-8W same-study controls",
	  GSE329024_2M_12W_control_meta = "GSE329024 2M-12W same-study controls",
	  GSE329024_6M_8W_control_meta = "GSE329024 6M-8W same-study controls",
	  GSE329024_6M_12W_control_meta = "GSE329024 6M-12W same-study controls",
	  GSE240547_control_diet = "GSE240547 control diet controls",
	  PRJNA1076509_young_female_donor_microbiota = "PRJNA1076509 young female donor microbiota controls",
	  GSE235021_polished_rice_6m = "GSE235021 6-month polished-rice controls",
	  GSE235021_polished_rice_9m = "GSE235021 9-month polished-rice controls",
	  GSE249151_sesame_oil_controls = "GSE249151 sesame-oil vehicle controls",
	  GSE128240_same_study_controls = "GSE128240 same-study controls",
	  GSE98497_WT_controls = "GSE98497 WT controls",
	  PRJNA787028_two_per_cage_controls = "PRJNA787028 two-mice-per-cage controls",
	  no_matched_control = "No matched control"
	)
  hit <- unname(labels[group_id])
	if (is.na(hit) || !nzchar(hit)) paste0(group_id, " controls") else hit
}

gse329024_control_group_id_one <- function(txt) {
	txt <- tolower(txt)
	start_age <- if (grepl("(^|[^a-z0-9])2m([^a-z0-9]|$)", txt)) {
		"2M"
	} else if (grepl("(^|[^a-z0-9])6m([^a-z0-9]|$)", txt)) {
		"6M"
	} else {
		NA_character_
	}
	duration <- if (grepl("(^|[^a-z0-9])12w([^a-z0-9]|$)", txt)) {
		"12W"
	} else if (grepl("(^|[^a-z0-9])8w([^a-z0-9]|$)", txt)) {
		"8W"
	} else {
		NA_character_
	}
	if (is.na(start_age) || is.na(duration)) return("GSE329024_same_dataset_control_meta")
	paste0("GSE329024_", start_age, "_", duration, "_control_meta")
}

default_control_group_id <- function(dt) {
	vapply(seq_len(nrow(dt)), function(i) {
	  row <- dt[i]
	  ds <- as.character(row$dataset_folder)
	  txt <- tolower(paste(row$sample_key, row$sample_id, row$accession, row$treatment, row$genotype, row$sample_context, row$sample_title, sep = " "))
    if (ds == "GSE194196") return("GSE194196_WT_1_4_7")
	  if (ds == "GSE252791") return("GSE252791_WT_1_4_7")
	  if (ds == "GSE279573") return("GSE279573_WT_chow_fertile_controls")
	  if (ds == "GSE329024") return(gse329024_control_group_id_one(txt))
    if (ds == "GSE240547") return("GSE240547_control_diet")
    if (ds == "PRJNA1076509") return("PRJNA1076509_young_female_donor_microbiota")
    if (ds == "GSE235021") {
      if (is.finite(safe_num(row$age_week)) && safe_num(row$age_week) < 30) {
        return("GSE235021_polished_rice_6m")
      }
      return("GSE235021_polished_rice_9m")
    }
    if (ds == "GSE249151") return("GSE249151_sesame_oil_controls")
    if (ds == "GSE128240") return("GSE128240_same_study_controls")
    if (ds == "GSE98497") return("GSE98497_WT_controls")
    if (ds == "PRJNA787028") return("PRJNA787028_two_per_cage_controls")
    if (ds == "GSE106339") return("GSE106339_prepubertal_WT_controls")
    if (ds == "GSE154890") return("no_matched_control")
    if (ds == "GSE279446") {
      if (gse279446_cohort_one(txt) == "young_cohort") return("GSE279446_young_WT_controls")
      return("GSE279446_adult_WT_ovary_controls")
    }
    if (ds == "GSE228975") return("no_matched_control")
    paste0(ds, "_same_dataset_control_meta")
  }, character(1))
}

assign_projection_roles <- function(dt, default_role = "external_projected_treatment") {
  dt <- copy(dt)
  if (!"projection_role" %in% names(dt)) dt[, projection_role := default_role]
  if (!"projection_class" %in% names(dt)) dt[, projection_class := NA_character_]
  if (!"projection_reason" %in% names(dt)) dt[, projection_reason := NA_character_]
  if (!"projection_dataset_note" %in% names(dt)) dt[, projection_dataset_note := NA_character_]
  if (!"control_pool_source" %in% names(dt)) dt[, control_pool_source := NA_character_]
  if (!"control_like_clock_candidate" %in% names(dt)) dt[, control_like_clock_candidate := FALSE]

  txt <- row_text(dt)
  dt[, `:=`(
    matched_reference_control_flag = FALSE,
    matched_reference_control_source = NA_character_,
    primary_external_projection_include = TRUE,
    primary_effect_interpretation_include = TRUE,
    matched_control_allowed = TRUE,
    gse279446_cohort = NA_character_,
    gse279446_platform = NA_character_,
    reference_control_group_id = NA_character_,
    reference_control_group_label = NA_character_
  )]

  dt[dataset_folder == "CRA003645" & age_week < 7, `:=`(
    matched_control_allowed = FALSE,
    matched_reference_control_source = "none",
    reference_control_group_id = "no_matched_control",
    reference_control_group_label = "No matched control",
    projection_dataset_note = paste(
      "Developmental/OOD placement only; the 16 adult CRA003645 samples are",
      "reserved exclusively for the frozen held-out test and are not projected."
    )
  )]

  dt[dataset_folder == "GSE228975", `:=`(
    projection_role = "OOD_supplementary_projection",
    projection_class = "age_bin_OOD_not_primary",
    projection_reason = "outside_current_age_bin_frame",
    primary_external_projection_include = FALSE,
    primary_effect_interpretation_include = FALSE,
    matched_control_allowed = FALSE,
    matched_reference_control_source = "none",
    reference_control_group_id = "no_matched_control",
    reference_control_group_label = "No matched control",
    projection_dataset_note = "Outside current age-bin frame; excluded from primary interpretation"
  )]

  gse240_ctrl <- dt$dataset_folder == "GSE240547" &
    grepl("control diet|control_diet|\\bctrl\\b", txt, ignore.case = TRUE)
  dt[gse240_ctrl, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "control_diet",
    projection_reason = "GSE240547_control_diet_projection_internal_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE,
    reference_control_group_id = "GSE240547_control_diet",
    reference_control_group_label = "GSE240547 control diet controls",
    projection_dataset_note = "HFHS diet compared against same-study control_diet projection-internal controls"
  )]
  gse240_hfhs <- dt$dataset_folder == "GSE240547" &
    grepl("hfhs|high-fat|high-sucrose", txt, ignore.case = TRUE)
  dt[gse240_hfhs, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "diet_intervention",
    reference_control_group_id = "GSE240547_control_diet",
    reference_control_group_label = "GSE240547 control diet controls",
    projection_dataset_note = "HFHS diet compared against same-study control_diet projection-internal controls"
  )]

  prjna_young <- dt$dataset_folder == "PRJNA1076509" &
    grepl("young female donor microbiota|fmt-yf", txt, ignore.case = TRUE)
  prjna_est <- dt$dataset_folder == "PRJNA1076509" &
    grepl("estropausal female donor microbiota|fmt-of", txt, ignore.case = TRUE)
  dt[prjna_young, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "young_female_donor_microbiota_control",
    projection_reason = "PRJNA1076509_young_female_donor_microbiota_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE,
    reference_control_group_id = "PRJNA1076509_young_female_donor_microbiota",
    reference_control_group_label = "PRJNA1076509 young female donor microbiota controls",
    projection_dataset_note = "Estropausal female donor microbiota compared against young female donor microbiota controls"
  )]
  dt[prjna_est, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "FMT_estropausal_female_donor_microbiota",
    projection_reason = "PRJNA1076509_estropausal_female_donor_microbiota",
    reference_control_group_id = "PRJNA1076509_young_female_donor_microbiota",
    reference_control_group_label = "PRJNA1076509 young female donor microbiota controls",
    projection_dataset_note = "Estropausal female donor microbiota compared against young female donor microbiota controls"
  )]

  gse279446 <- dt$dataset_folder == "GSE279446"
  if (any(gse279446)) {
    dt[gse279446, gse279446_cohort := vapply(txt[gse279446], gse279446_cohort_one, character(1))]
    dt[gse279446, gse279446_platform := vapply(gse279446_cohort, gse279446_platform_one, character(1))]
    gse279446_wt <- gse279446 & grepl("(^|[^a-z0-9])wt[_ ]ovary|zeb1[_ ]wt[_ ]ovary|\\bwt ovary", txt, ignore.case = TRUE)
    gse279446_h <- gse279446 & grepl("zeb1[_ ]200h", txt, ignore.case = TRUE)
    gse279446_m <- gse279446 & grepl("zeb1[_ ]200m", txt, ignore.case = TRUE)
    dt[gse279446_wt, `:=`(
      projection_role = "projection_internal_control",
      projection_class = "same_study_WT_control",
      projection_reason = "GSE279446_WT_same_cohort_control",
      matched_reference_control_flag = TRUE,
      matched_reference_control_source = "control_meta",
      primary_effect_interpretation_include = FALSE
    )]
    dt[gse279446_h, `:=`(
      projection_role = "external_projected_treatment",
      projection_class = fifelse(gse279446_cohort == "young_cohort", "genetic_model_heterozygous_young", "genetic_model_heterozygous"),
      projection_reason = "GSE279446_Zeb1_200H_genetic_model"
    )]
    dt[gse279446_m, `:=`(
      projection_role = "external_projected_treatment",
      projection_class = fifelse(gse279446_cohort == "young_cohort", "genetic_model_homozygous_young", "genetic_model_homozygous"),
      projection_reason = "GSE279446_Zeb1_200M_genetic_model"
    )]
    dt[gse279446 & gse279446_cohort == "young_cohort", `:=`(
      reference_control_group_id = "GSE279446_young_WT_controls",
      reference_control_group_label = "GSE279446 young WT ovary controls"
    )]
    dt[gse279446 & gse279446_cohort == "adult_cohort", `:=`(
      reference_control_group_id = "GSE279446_adult_WT_ovary_controls",
      reference_control_group_label = "GSE279446 adult WT ovary controls"
    )]
  }

  gse279573_wt <- dt$dataset_folder == "GSE279573" & grepl("wt chow fertile", txt, ignore.case = TRUE)
  dt[gse279573_wt, `:=`(
    projection_role = "control_meta_reference_control",
    projection_class = "same_study_WT_chow_fertile_control",
    projection_reason = "GSE279573_WT_chow_fertile_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "control_meta",
    primary_effect_interpretation_include = FALSE,
    reference_control_group_id = "GSE279573_WT_chow_fertile_controls",
    reference_control_group_label = "GSE279573 WT chow fertile controls",
    projection_dataset_note = "KO groups compared against same-study WT chow fertile controls; n_controls = 2, descriptive interpretation"
  )]
  dt[dataset_folder == "GSE279573" & grepl("ko chow infertile", txt, ignore.case = TRUE), `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "genetic_model",
    projection_reason = "GSE279573_Abcg5g8_KO_chow_infertile",
    reference_control_group_id = "GSE279573_WT_chow_fertile_controls",
    reference_control_group_label = "GSE279573 WT chow fertile controls"
  )]
  dt[dataset_folder == "GSE279573" & grepl("vegetable-oil-free|vof", txt, ignore.case = TRUE), `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "diet_rescue_or_diet_intervention",
    projection_reason = "GSE279573_Abcg5g8_KO_vegetable_oil_free_rescue",
    reference_control_group_id = "GSE279573_WT_chow_fertile_controls",
    reference_control_group_label = "GSE279573 WT chow fertile controls"
  )]

  bnc_wt <- dt$dataset_folder %chin% c("GSE194196", "GSE252791") &
    grepl("\\bwt[-_ ]?[147123]?\\b|\\bwt ovary\\b|estrus wt", txt, ignore.case = TRUE) &
    !grepl("bnc1|het|ko|heterozygous", txt, ignore.case = TRUE)
  dt[bnc_wt, `:=`(
    projection_role = "control_meta_reference_control",
    projection_class = "same_study_WT_control",
    projection_reason = paste0(dataset_folder, "_WT_reference_control"),
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "control_meta",
    primary_effect_interpretation_include = FALSE
  )]
  dt[dataset_folder == "GSE194196", `:=`(
    reference_control_group_id = "GSE194196_WT_1_4_7",
    reference_control_group_label = "GSE194196 WT-1/WT-4/WT-7 controls"
  )]
	  dt[dataset_folder == "GSE252791", `:=`(
	    reference_control_group_id = "GSE252791_WT_1_4_7",
	    reference_control_group_label = "GSE252791 WT1/WT2/WT3 controls"
	  )]

	  gse329024 <- dt$dataset_folder == "GSE329024"
	  if (any(gse329024)) {
	    gse329024_group_id <- vapply(txt[gse329024], gse329024_control_group_id_one, character(1))
	    dt[gse329024, `:=`(
	      reference_control_group_id = gse329024_group_id,
	      reference_control_group_label = vapply(gse329024_group_id, reference_control_group_label_one, character(1)),
	      projection_dataset_note = "GSE329024 CUMS samples are compared against controls matched by starting age and exposure duration"
	    )]
	    dt[gse329024 & control_like_clock_candidate == TRUE, `:=`(
	      projection_role = "control_meta_reference_control",
	      projection_class = "same_study_reference_control",
	      projection_reason = "GSE329024_design_matched_reference_control",
	      matched_reference_control_flag = TRUE,
	      matched_reference_control_source = "control_meta",
	      primary_effect_interpretation_include = FALSE
	    )]
	    dt[gse329024 & control_like_clock_candidate != TRUE, `:=`(
	      projection_role = "external_projected_treatment",
	      projection_class = "stress",
	      projection_reason = "GSE329024_CUMS_stress_model"
	    )]
	  }

	  std_control <- dt$dataset_folder %chin% c("CRA003645", "CRA017492", "GSE273724") &
	    dt$control_like_clock_candidate == TRUE
	  dt[std_control, `:=`(
    projection_role = "control_meta_reference_control",
    projection_class = "same_study_reference_control",
    projection_reason = paste0(dataset_folder, "_control_meta_reference_control"),
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "control_meta",
    primary_effect_interpretation_include = FALSE,
    reference_control_group_id = paste0(dataset_folder, "_same_dataset_control_meta"),
    reference_control_group_label = paste0(dataset_folder, " same-dataset control_meta controls")
  )]

  gse235021 <- dt$dataset_folder == "GSE235021"
  gse235021_control <- gse235021 & grepl("polished rice", txt, ignore.case = TRUE)
  dt[gse235021, `:=`(
    reference_control_group_id = fifelse(
      age_week < 30,
      "GSE235021_polished_rice_6m",
      "GSE235021_polished_rice_9m"
    ),
    reference_control_group_label = fifelse(
      age_week < 30,
      "GSE235021 6-month polished-rice controls",
      "GSE235021 9-month polished-rice controls"
    ),
    projection_dataset_note = "Black-rice and white-brown-rice groups are compared with pooled bulk and scRNA-pseudobulk polished-rice controls within age."
  )]
  dt[gse235021_control, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "control_diet",
    projection_reason = "GSE235021_polished_rice_projection_internal_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE
  )]
  dt[gse235021 & !gse235021_control, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "diet_intervention",
    projection_reason = "GSE235021_rice_diet_intervention"
  )]

  dt[dataset_folder == "GSE154890", `:=`(
    projection_role = "external_reference_placement",
    projection_class = "natural_age_reference_control",
    projection_reason = "GSE154890_natural_age_external_validation",
    matched_reference_control_flag = FALSE,
    matched_reference_control_source = "none",
    primary_effect_interpretation_include = FALSE,
    matched_control_allowed = FALSE,
    reference_control_group_id = "no_matched_control",
    reference_control_group_label = "No matched control",
    projection_dataset_note = "Natural-age external placement cohort; absolute clock placement only."
  )]

  gse249151 <- dt$dataset_folder == "GSE249151"
  gse249151_control <- gse249151 & grepl("sesame oil control|vehicle control", txt, ignore.case = TRUE)
  dt[gse249151, `:=`(
    reference_control_group_id = "GSE249151_sesame_oil_controls",
    reference_control_group_label = "GSE249151 sesame-oil vehicle controls",
    projection_dataset_note = "VCD day-15 is age matched; VCD day-30 uses the available younger vehicle group and is age-confounded."
  )]
  dt[gse249151_control, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "same_study_reference_control",
    projection_reason = "GSE249151_sesame_oil_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE
  )]
  dt[gse249151 & !gse249151_control, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "intervention",
    projection_reason = "GSE249151_VCD_ovarian_insufficiency_intervention"
  )]

  gse106339 <- dt$dataset_folder == "GSE106339"
  gse106339_wt <- gse106339 & toupper(trimws(dt$genotype)) == "WT"
  dt[gse106339, `:=`(
    projection_reason = "GSE106339_GOAT_genetic_model_with_prepubertal_OOD_flag",
    primary_external_projection_include = FALSE,
    primary_effect_interpretation_include = FALSE,
    matched_control_allowed = TRUE,
    reference_control_group_id = "GSE106339_prepubertal_WT_controls",
    reference_control_group_label = "GSE106339 prepubertal WT controls",
    projection_dataset_note = "GOAT-KO samples are genetic models and WT samples are reference controls; all remain prepubertal/OOD and descriptive."
  )]
  dt[gse106339_wt, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "same_study_reference_control",
    projection_reason = "GSE106339_prepubertal_WT_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control"
  )]
  dt[gse106339 & !gse106339_wt, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "genetic_model",
    matched_reference_control_flag = FALSE,
    matched_reference_control_source = "projection_internal_control"
  )]

  gse128240 <- dt$dataset_folder == "GSE128240"
  gse128240_control <- gse128240 & grepl("(^|[^a-z])control([^a-z]|$)", txt, ignore.case = TRUE) &
    !grepl("cyclophosphamide|ctx", txt, ignore.case = TRUE)
  dt[gse128240, `:=`(
    reference_control_group_id = "GSE128240_same_study_controls",
    reference_control_group_label = "GSE128240 same-study controls",
    projection_dataset_note = "Cyclophosphamide ovarian-injury samples are compared with same-study controls."
  )]
  dt[gse128240_control, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "same_study_reference_control",
    projection_reason = "GSE128240_same_study_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE
  )]
  dt[gse128240 & !gse128240_control, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "intervention",
    projection_reason = "GSE128240_cyclophosphamide_intervention"
  )]

  gse98497 <- dt$dataset_folder == "GSE98497"
  gse98497_wt <- gse98497 & grepl("(^|[^a-z0-9])wt([^a-z0-9]|$)|health state control", txt, ignore.case = TRUE) &
    !grepl("cko|knockout", txt, ignore.case = TRUE)
  dt[gse98497, `:=`(
    reference_control_group_id = "GSE98497_WT_controls",
    reference_control_group_label = "GSE98497 WT controls",
    projection_dataset_note = "Mtor oocyte-conditional knockout samples are compared with same-study WT controls; PMSG is common to both groups."
  )]
  dt[gse98497_wt, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "same_study_reference_control",
    projection_reason = "GSE98497_WT_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE
  )]
  dt[gse98497 & !gse98497_wt, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "genetic_model",
    projection_reason = "GSE98497_Mtor_oocyte_conditional_KO"
  )]

  prjna787028 <- dt$dataset_folder == "PRJNA787028"
  prjna787028_control <- prjna787028 & grepl("2 mice per cage|2pc", txt, ignore.case = TRUE)
  dt[prjna787028, `:=`(
    reference_control_group_id = "PRJNA787028_two_per_cage_controls",
    reference_control_group_label = "PRJNA787028 two-mice-per-cage controls",
    projection_dataset_note = "Eight-mice-per-cage housing-density stress is compared with two-mice-per-cage controls."
  )]
  dt[prjna787028_control, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "same_study_reference_control",
    projection_reason = "PRJNA787028_two_per_cage_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE
  )]
  dt[prjna787028 & !prjna787028_control, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "stress",
    projection_reason = "PRJNA787028_high_density_housing_stress"
  )]

  prjna1002222 <- dt$dataset_folder == "PRJNA1002222"
  prjna1002222_wt <- prjna1002222 & grepl("(^|[^a-z0-9])wt([^a-z0-9]|$)", txt, ignore.case = TRUE) &
    !grepl("cd38|knockout", txt, ignore.case = TRUE)
  prjna1002222_2m <- prjna1002222 & abs(dt$age_week - 2 * 52.1775 / 12) < 1
  prjna1002222_8m <- prjna1002222 & abs(dt$age_week - 8 * 52.1775 / 12) < 1
  prjna1002222_12m <- prjna1002222 & abs(dt$age_week - 12 * 52.1775 / 12) < 1
  for (age_info in list(
    list(index = prjna1002222_8m, label = "8M"),
    list(index = prjna1002222_12m, label = "12M")
  )) {
    idx <- age_info$index
    label <- age_info$label
    dt[idx, `:=`(
      reference_control_group_id = paste0("PRJNA1002222_Cd38KO_vs_WT_", label),
      reference_control_group_label = paste0("PRJNA1002222 WT ", label, " controls"),
      projection_dataset_note = "Cd38-knockout whole ovaries are compared with same-age WT controls."
    )]
    dt[idx & prjna1002222_wt, `:=`(
      projection_role = "projection_internal_control",
      projection_class = "same_study_reference_control",
      projection_reason = "PRJNA1002222_same_age_WT_reference_control",
      matched_reference_control_flag = TRUE,
      matched_reference_control_source = "projection_internal_control",
      primary_effect_interpretation_include = FALSE
    )]
    dt[idx & !prjna1002222_wt, `:=`(
      projection_role = "external_projected_treatment",
      projection_class = "genetic_model",
      projection_reason = "PRJNA1002222_Cd38_knockout_genetic_model"
    )]
  }
  dt[prjna1002222_2m, `:=`(
    projection_role = "control_meta_reference_control",
    projection_class = "young_natural_age_reference_control",
    projection_reason = "PRJNA1002222_2M_WT_natural_age_reference",
    matched_reference_control_flag = FALSE,
    matched_reference_control_source = "none",
    primary_effect_interpretation_include = FALSE,
    matched_control_allowed = FALSE,
    reference_control_group_id = "no_matched_control",
    reference_control_group_label = "No matched control",
    projection_dataset_note = "Two-month WT ovaries are placement-only young references; no 2M Cd38-knockout group exists."
  )]

  cra017492 <- dt$dataset_folder == "CRA017492"
  cra017492_young <- cra017492 & grepl("young 2 months control", txt, ignore.case = TRUE)
  cra017492_aged_water <- cra017492 & grepl("aged 12 months water control|12 months water control", txt, ignore.case = TRUE)
  cra017492_aged_nmn <- cra017492 & grepl("aged 12 months plus NMN|12 months.*NMN", txt, ignore.case = TRUE)
  dt[cra017492_young, `:=`(
    projection_role = "control_meta_reference_control",
    projection_class = "young_natural_age_reference_control",
    projection_reason = "CRA017492_young_natural_age_reference_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "control_meta",
    primary_effect_interpretation_include = FALSE,
    reference_control_group_id = "CRA017492_young_control_reference",
    reference_control_group_label = "CRA017492 young 2-month natural-age reference controls",
    projection_dataset_note = "Young controls are retained as a natural-age reference, not as the NMN treatment comparator"
  )]
  dt[cra017492_aged_water, `:=`(
    projection_role = "projection_internal_control",
    projection_class = "aged_water_control",
    projection_reason = "CRA017492_aged_water_age_matched_control",
    matched_reference_control_flag = TRUE,
    matched_reference_control_source = "projection_internal_control",
    primary_effect_interpretation_include = FALSE,
    reference_control_group_id = "CRA017492_aged_water_controls",
    reference_control_group_label = "CRA017492 aged 12-month water controls",
    projection_dataset_note = "Aged water controls are the age-matched comparator for aged NMN"
  )]
  dt[cra017492_aged_nmn, `:=`(
    projection_role = "external_projected_treatment",
    projection_class = "NMN_intervention",
    projection_reason = "CRA017492_aged_NMN_intervention",
    matched_reference_control_flag = FALSE,
    matched_reference_control_source = NA_character_,
    reference_control_group_id = "CRA017492_aged_water_controls",
    reference_control_group_label = "CRA017492 aged 12-month water controls",
    projection_dataset_note = "Aged NMN is compared against age-matched aged water controls; young controls are shown only as a natural-age reference"
  )]

  dt[is.na(reference_control_group_id), reference_control_group_id := default_control_group_id(.SD)]
  dt[is.na(reference_control_group_label), reference_control_group_label := vapply(reference_control_group_id, reference_control_group_label_one, character(1))]
  dt[projection_role %chin% c("control_meta_reference_control", "projection_internal_control") &
       is.na(matched_reference_control_source),
     matched_reference_control_source := fifelse(control_pool_source == "projection_internal_control", "projection_internal_control", "control_meta")]
  dt
}

make_projection_metadata <- function(meta) {
  meta <- copy(meta)
  control_col <- detect_col(meta, c("control_like_clock_candidate", "control_like", "clock_candidate"), "control-like flag")
  exclusion_col <- detect_col(meta, c("exclusion_reason", "exclude_reason", "exclusion"), "exclusion reason", required = FALSE)
  if (control_col != "control_like_clock_candidate") setnames(meta, control_col, "control_like_clock_candidate")
  if (!is.na(exclusion_col) && exclusion_col != "exclusion_reason") setnames(meta, exclusion_col, "exclusion_reason")
  if (!"exclusion_reason" %in% names(meta)) meta[, exclusion_reason := NA_character_]
  meta[, control_like_clock_candidate := as.logical(control_like_clock_candidate)]
  meta[, exclusion_reason := na_chr(exclusion_reason)]

  projection_flag <- meta$control_like_clock_candidate != TRUE | !is.na(meta$exclusion_reason)
  if ("projection_candidate" %in% names(meta)) {
    projection_flag <- projection_flag | as.logical(meta$projection_candidate) == TRUE
  }
  projection_flag[is.na(projection_flag)] <- FALSE
  proj <- copy(meta[projection_flag])
  if (nrow(proj) == 0) stop("No external projection samples detected from metadata flags.", call. = FALSE)
  if (any(proj$control_like_clock_candidate == TRUE & is.na(proj$exclusion_reason))) {
    stop("Control-like samples were mixed into the external projection set without an exclusion reason.", call. = FALSE)
  }
  proj[, projection_role := "external_projected_treatment"]
  proj[, projection_class := vapply(seq_len(.N), function(i) projection_class_one(.SD[i]), character(1))]
  proj[, projection_reason := fifelse(!is.na(exclusion_reason), exclusion_reason, sample_context %||% "metadata_flagged_projection")]
  assign_projection_roles(proj)
}

model_matrix_for_genes <- function(expr_pieces, meta_dt, model_genes, center, scale) {
  out <- matrix(NA_real_, nrow = nrow(meta_dt), ncol = length(model_genes))
  rownames(out) <- meta_dt$sample_key
  colnames(out) <- model_genes
  observed <- vector("list", nrow(meta_dt))
  names(observed) <- meta_dt$sample_key

  for (piece_name in names(expr_pieces)) {
    mat <- expr_pieces[[piece_name]]
    keys <- intersect(colnames(mat), meta_dt$sample_key)
    if (length(keys) == 0) next
    genes_present <- intersect(model_genes, rownames(mat))
    if (length(genes_present) > 0) {
      out[keys, genes_present] <- t(mat[genes_present, keys, drop = FALSE])
    }
    for (key in keys) observed[[key]] <- genes_present
  }

  missing <- is.na(out)
  raw <- out
  for (gene in model_genes) {
    raw[missing[, gene], gene] <- center[[gene]]
  }
  scaled <- sweep(raw, 2, center[model_genes], "-")
  scaled <- sweep(scaled, 2, scale[model_genes], "/")
  scaled[!is.finite(scaled)] <- 0

  coverage <- data.table(
    sample_key = rownames(out),
    n_model_genes = length(model_genes),
    n_observed_model_genes = rowSums(!missing),
    gene_coverage_fraction = rowSums(!missing) / length(model_genes),
    missing_model_genes = apply(missing, 1, function(x) paste(model_genes[x], collapse = ";"))
  )
  coverage[, projection_gene_coverage_flag := fifelse(gene_coverage_fraction < 0.80, "low_coverage", "ok")]
  list(x = scaled, coverage = coverage)
}

robust_mahalanobis <- function(x, ref_x) {
  if (nrow(ref_x) < 3 || ncol(ref_x) < 1) return(rep(NA_real_, nrow(x)))
  center <- colMeans(ref_x, na.rm = TRUE)
  cov_mat <- stats::cov(ref_x, use = "pairwise.complete.obs")
  cov_mat <- as.matrix(cov_mat)
  cov_mat[!is.finite(cov_mat)] <- 0
  diag(cov_mat) <- diag(cov_mat) + 1e-6
  as.numeric(stats::mahalanobis(x, center = center, cov = cov_mat))
}

distance_percentile <- function(distance, ref_distance) {
  vapply(distance, function(d) {
    if (!is.finite(d) || all(!is.finite(ref_distance))) return(NA_real_)
    mean(ref_distance <= d, na.rm = TRUE)
  }, numeric(1))
}

clock_predict <- function(clock, obj, cal, expr_pieces, meta_apply, meta_reference) {
  genes <- obj$genes
  if (is.null(genes) || length(genes) == 0) {
    genes <- fread(clock$selected_genes_path)$gene_id
  }
  genes <- as.character(genes)
  center <- obj$center[genes]
  scale <- obj$scale[genes]
  scale[!is.finite(scale) | scale == 0] <- 1

  apply_x <- model_matrix_for_genes(expr_pieces, meta_apply, genes, center, scale)
  ref_x <- model_matrix_for_genes(expr_pieces, meta_reference, genes, center, scale)

  if (clock$model_family == "elasticnet") {
    raw_apply <- as.numeric(predict(obj$model, newx = apply_x$x))
    raw_ref <- as.numeric(predict(obj$model, newx = ref_x$x))
    pred <- apply_calibration(raw_apply, cal)
    pred_ref <- apply_calibration(raw_ref, cal)
    centroid <- colMeans(ref_x$x, na.rm = TRUE)
    dist_apply <- sqrt(rowSums(sweep(apply_x$x, 2, centroid, "-")^2))
    dist_ref <- sqrt(rowSums(sweep(ref_x$x, 2, centroid, "-")^2))
    pctl <- distance_percentile(dist_apply, dist_ref)
    pred_dt <- data.table(
      sample_key = meta_apply$sample_key,
      model_family = clock$model_family,
      selector = clock$selector,
      raw_predicted_age_week = raw_apply,
      predicted_age_week = pred,
      projection_distance = dist_apply,
      projection_distance_percentile = pctl,
      OOD_flag = fifelse(pctl > 0.99, "OOD_gt_99th", fifelse(pctl > 0.95, "OOD_gt_95th", "in_reference_range"))
    )
    ref_dt <- data.table(
      sample_key = meta_reference$sample_key,
      model_family = clock$model_family,
      selector = clock$selector,
      predicted_age_week = pred_ref,
      delta_age_week = pred_ref - safe_num(meta_reference$age_week),
      projection_distance = dist_ref
    )
    return(list(predictions = merge(pred_dt, apply_x$coverage, by = "sample_key"), reference_predictions = ref_dt, pca_scores = NULL))
  }

  pc_count <- as.integer(obj$pc_count)
  pc_apply_all <- predict(obj$pca, newdata = apply_x$x)
  pc_ref_all <- predict(obj$pca, newdata = ref_x$x)
  pc_cols <- paste0("PC", seq_len(pc_count))
  pc_apply <- pc_apply_all[, seq_len(pc_count), drop = FALSE]
  pc_ref <- pc_ref_all[, seq_len(pc_count), drop = FALSE]
  colnames(pc_apply) <- pc_cols
  colnames(pc_ref) <- pc_cols
  raw_apply <- as.numeric(predict(obj$model, newx = pc_apply))
  raw_ref <- as.numeric(predict(obj$model, newx = pc_ref))
  pred <- apply_calibration(raw_apply, cal)
  pred_ref <- apply_calibration(raw_ref, cal)
  dist_apply <- robust_mahalanobis(pc_apply, pc_ref)
  dist_ref <- robust_mahalanobis(pc_ref, pc_ref)
  pctl <- distance_percentile(dist_apply, dist_ref)
  pred_dt <- data.table(
    sample_key = meta_apply$sample_key,
    model_family = clock$model_family,
    selector = clock$selector,
    raw_predicted_age_week = raw_apply,
    predicted_age_week = pred,
    projection_distance = dist_apply,
    projection_distance_percentile = pctl,
    OOD_flag = fifelse(pctl > 0.99, "OOD_gt_99th", fifelse(pctl > 0.95, "OOD_gt_95th", "in_reference_range"))
  )
  beta <- as.numeric(coef(obj$model))[-1]
  beta <- beta[seq_len(pc_count)]
  score_dt <- as.data.table(pc_apply)
  score_dt[, sample_key := meta_apply$sample_key]
  setcolorder(score_dt, "sample_key")
  for (j in seq_len(pc_count)) {
    score_dt[, paste0(pc_cols[[j]], "_contribution_week") := pc_apply[, j] * beta[[j]]]
  }
  score_dt[, `:=`(
    model_family = clock$model_family,
    selector = clock$selector,
    total_calibrated_predicted_age_week = pred,
    delta_age_week = pred - safe_num(meta_apply$age_week)
  )]
  ref_dt <- data.table(
    sample_key = meta_reference$sample_key,
    model_family = clock$model_family,
    selector = clock$selector,
    predicted_age_week = pred_ref,
    delta_age_week = pred_ref - safe_num(meta_reference$age_week),
    projection_distance = dist_ref
  )
  ref_scores <- as.data.table(pc_ref[, seq_len(min(2L, pc_count)), drop = FALSE])
  ref_scores[, sample_key := meta_reference$sample_key]
  list(
    predictions = merge(pred_dt, apply_x$coverage, by = "sample_key"),
    reference_predictions = ref_dt,
    pca_scores = score_dt,
    reference_pca_scores = ref_scores
  )
}

make_matched_control_metadata <- function(meta, control_meta) {
  matched_control_meta <- copy(control_meta)
  matched_control_meta[, `:=`(
    control_pool_source = "control_meta",
    control_candidate_status = "accepted_control_meta",
    control_candidate_note = "control-like deployment/reference metadata"
  )]
  gse240_control_diet <- meta[
    dataset_folder == "GSE240547" &
      grepl("control diet|control_diet|\\bctrl\\b", paste(sample_id, accession, treatment, sample_title, sep = " "), ignore.case = TRUE)
  ]
  if (nrow(gse240_control_diet) > 0) {
    gse240_control_diet[, `:=`(
      control_pool_source = "projection_internal_control",
      control_candidate_status = "accepted_projection_internal_control",
      control_candidate_note = "GSE240547 control_diet projection-internal same-study controls"
    )]
    matched_control_meta <- rbindlist(
      list(matched_control_meta, gse240_control_diet),
      use.names = TRUE,
      fill = TRUE
    )
    matched_control_meta <- unique(matched_control_meta, by = "sample_key")
  }
  prjna_young_donor <- meta[
    dataset_folder == "PRJNA1076509" &
      grepl(
        "young_donor_microbiota_reference|young female donor microbiota",
        paste(sample_context, treatment, sample_title, sep = " "),
        ignore.case = TRUE
      )
  ]
  if (nrow(prjna_young_donor) > 0) {
    prjna_young_donor[, `:=`(
      control_pool_source = "projection_internal_control",
      control_candidate_status = "accepted_projection_internal_control",
      control_candidate_note = "PRJNA1076509 young female donor microbiota controls"
    )]
    matched_control_meta <- rbindlist(
      list(matched_control_meta, prjna_young_donor),
      use.names = TRUE,
      fill = TRUE
    )
    matched_control_meta <- unique(matched_control_meta, by = "sample_key")
  }
  cra017492_aged_water <- meta[
    dataset_folder == "CRA017492" &
      grepl(
        "aged 12 months water control|12 months water control",
        paste(sample_context, treatment, sample_title, sep = " "),
        ignore.case = TRUE
      )
  ]
  if (nrow(cra017492_aged_water) > 0) {
    cra017492_aged_water[, `:=`(
      control_pool_source = "projection_internal_control",
      control_candidate_status = "accepted_projection_internal_control",
      control_candidate_note = "CRA017492 age-matched aged water controls for aged NMN"
    )]
    matched_control_meta <- rbindlist(
      list(matched_control_meta, cra017492_aged_water),
      use.names = TRUE,
      fill = TRUE
    )
    matched_control_meta <- unique(matched_control_meta, by = "sample_key")
  }
  new_projection_internal_controls <- meta[
    (dataset_folder == "GSE235021" &
       grepl("polished rice", paste(treatment, sample_title, sample_id, sep = " "), ignore.case = TRUE)) |
    (dataset_folder == "GSE249151" &
       grepl("sesame oil control|vehicle control", paste(treatment, sample_title, sample_id, sep = " "), ignore.case = TRUE)) |
    (dataset_folder == "GSE128240" &
       grepl("(^|[^a-z])control([^a-z]|$)", paste(treatment, sample_title, sample_id, sep = " "), ignore.case = TRUE) &
       !grepl("cyclophosphamide|ctx", paste(treatment, sample_title, sample_id, sep = " "), ignore.case = TRUE)) |
    (dataset_folder == "GSE98497" &
       grepl("(^|[^a-z0-9])wt([^a-z0-9]|$)", paste(genotype, sample_title, sample_id, sep = " "), ignore.case = TRUE) &
       !grepl("cko|knockout", paste(genotype, sample_title, sample_id, sep = " "), ignore.case = TRUE)) |
    (dataset_folder == "PRJNA787028" &
       grepl("2 mice per cage|2pc", paste(treatment, sample_title, sample_id, sep = " "), ignore.case = TRUE))
  ]
  if (nrow(new_projection_internal_controls) > 0) {
    new_projection_internal_controls[, `:=`(
      control_pool_source = "projection_internal_control",
      control_candidate_status = "accepted_projection_internal_control",
      control_candidate_note = "Prespecified same-study control for added external projection cohort"
    )]
    matched_control_meta <- rbindlist(
      list(matched_control_meta, new_projection_internal_controls),
      use.names = TRUE,
      fill = TRUE
    )
    matched_control_meta <- unique(matched_control_meta, by = "sample_key")
  }
  matched_control_meta
}

compute_relative_to_controls <- function(pred, control_pred, meta_controls) {
  controls <- merge(
    control_pred,
    meta_controls[, .(
      sample_key, dataset_folder, modality, age_week, bio_age_bin, sample_context,
      treatment, sample_id, accession, sample_title, control_pool_source,
      projection_role, matched_reference_control_source, reference_control_group_id,
      reference_control_group_label, gse279446_cohort
    )],
    by = "sample_key",
    all.x = TRUE
  )
  pred[, row_id := seq_len(.N)]
  out <- pred[, {
    use <- data.table()
    match_type <- "no_matched_control"
    match_source <- "none"
    warning <- NA_character_
    ctrl <- controls[
      model_family == .BY$model_family &
        selector == .BY$selector &
        dataset_folder == .BY$dataset_folder
    ]
    if (.BY$dataset_folder != "GSE235021") {
      ctrl <- ctrl[modality == .BY$modality]
    }
    if (.BY$matched_control_allowed == FALSE || .BY$primary_external_projection_include == FALSE) {
      ctrl <- data.table()
      warning <- "matched controls disabled for non-primary/OOD or projection-internal control row"
    }
    if (.BY$matched_control_allowed == TRUE && .BY$primary_external_projection_include == TRUE) {
      group_id <- .BY$reference_control_group_id %||% NA_character_
      if (!is.na(group_id) && group_id != "" && group_id != "no_matched_control") {
        use <- ctrl[reference_control_group_id == group_id]
        if (nrow(use) > 0) {
          match_type <- group_id
          match_source <- unique(na_chr(use$matched_reference_control_source %||% use$control_pool_source))[1] %||% "control_meta"
        }
      }
	      if (nrow(use) == 0 && .BY$dataset_folder %chin% c(
	        "GSE194196", "GSE252791", "GSE279573", "GSE279446", "GSE329024",
	        "GSE240547", "PRJNA1076509", "GSE235021", "GSE249151",
	        "GSE128240", "GSE98497", "PRJNA787028", "PRJNA1002222"
	      )) {
	        ctrl <- data.table()
	        warning <- paste0(.BY$dataset_folder, " required matched control group not detected: ", group_id)
	      }
    }
    if (nrow(use) == 0 && nrow(ctrl) > 0) {
      exact <- ctrl[is.finite(age_week) & age_week == .BY$age_week]
      if (nrow(exact) > 0) {
        use <- exact
        match_type <- "same_dataset_same_age_control_meta"
      } else {
        nearest_age <- min(abs(ctrl$age_week - .BY$age_week), na.rm = TRUE)
        if (is.finite(nearest_age)) {
          use <- ctrl[abs(age_week - .BY$age_week) == nearest_age]
          match_type <- "same_dataset_nearest_age_control_meta"
        } else {
          use <- ctrl[bio_age_bin == .BY$bio_age_bin]
          if (nrow(use) > 0) match_type <- "same_dataset_same_biological_age_bin_control_meta"
        }
      }
      if (nrow(use) > 0) {
        match_source <- unique(na_chr(use$control_pool_source))[1] %||% "control_meta"
      }
    }
    if (nrow(use) == 0) {
      data.table(
        control_expected_predicted_age_week = NA_real_,
        control_expected_delta_age_week = NA_real_,
        relative_delta_age_week = NA_real_,
        matched_control_method = match_type,
        matched_control_source = "none",
        n_matched_controls = 0L,
        matched_control_sample_keys = NA_character_,
        relative_control_reliability = "not_available",
        relative_control_warning = warning %||% "no same-dataset matched control pool available"
      )
    } else {
      control_expected_delta <- mean(use$delta_age_week, na.rm = TRUE)
      data.table(
        control_expected_predicted_age_week = mean(use$predicted_age_week, na.rm = TRUE),
        control_expected_delta_age_week = control_expected_delta,
        relative_delta_age_week = .BY$delta_age_week - control_expected_delta,
        matched_control_method = match_type,
        matched_control_source = match_source,
        n_matched_controls = nrow(use),
        matched_control_sample_keys = paste(sort(unique(use$sample_key)), collapse = ";"),
        relative_control_reliability = fifelse(nrow(use) >= 3, "valid_same_study_control_pool", "limited_control_pool"),
        relative_control_warning = NA_character_
      )
    }
  }, by = .(
    row_id, sample_key, model_family, selector, dataset_folder, modality, age_week,
    bio_age_bin, delta_age_week, primary_external_projection_include, matched_control_allowed,
    projection_role, reference_control_group_id
  )]
  pred[, row_id := NULL]
  out[, c("row_id", "delta_age_week") := NULL]
  out[, `:=`(
    control_expected_pred_age = control_expected_predicted_age_week,
    control_expected_delta_age = control_expected_delta_age_week,
    relative_delta_age = relative_delta_age_week
  )]
  out
}

save_plot <- function(p, filename, width = 10, height = 7) {
  png_path <- file.path(plot_dir, paste0(filename, ".png"))
  pdf_path <- file.path(plot_dir, paste0(filename, ".pdf"))
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(pdf_path, p, width = width, height = height, bg = "white", limitsize = FALSE)
}

theme_projection <- function() {
  theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

label_projection_class <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & x == "stress"] <- "stress_model"
  x
}

add_projection_class_label <- function(dt) {
  if (!is.null(dt) && "projection_class" %in% names(dt)) {
    dt[, projection_class_label := label_projection_class(projection_class)]
  }
  dt
}

write_dataset_inclusion_notes <- function() {
  notes_path <- file.path(project_dir, "metadata", "external_projection_dataset_inclusion_notes.md")
  lines <- c(
    "# External projection dataset inclusion notes",
    "",
    "## Primary external-projection frame",
    "The chronological clocks were trained on control-like reference samples only. External projection positions excluded intervention, genotype, diet, microbiome/FMT, stress, perturbation, postpartum/lactation, developmental, or OOD samples into the frozen ElasticNet and PCA/glmnet clock spaces. Projected samples are never used for fitting, scaling, PCA, glmnet, or calibration.",
    "",
    "## Age-bin frame",
    "Current biological age bins:",
    "- young_adult: 7 to <16 weeks",
    "- mature_adult: 16 to <40 weeks",
    "- advanced_maternal_age: 40 to <65 weeks",
    "- post_reproductive_ageing: >=65 weeks",
    "",
    "Rows outside this frame are OOD/supplementary unless explicitly approved.",
    "",
    "## Matched-control interpretation",
    "Raw delta age is computed for projected samples with sufficient gene coverage.",
    "Matched-control relative delta is computed only when a valid same-study control pool is available.",
    "",
    "## Dataset-specific control logic",
    "",
    "### CRA003645",
    "Use same-dataset control_meta controls when available.",
    "",
    "### CRA017492",
    "Use same-dataset control_meta controls when available.",
    "",
    "### GSE194196",
    "Use same-dataset control_meta controls when available.",
    "",
    "### GSE252791",
    "Use same-dataset control_meta controls when available.",
    "",
    "### GSE273724",
    "Use same-dataset control_meta controls when available.",
    "",
    "### GSE279446",
    "Use same-dataset control_meta controls when available.",
    "",
    "### GSE279573",
    "Use same-dataset control_meta controls when available.",
    "",
	    "### GSE329024",
	    "CUMS samples are compared against same-study controls matched by starting age and exposure duration.",
	    "Expected groups: 2M-8W, 2M-12W, 6M-8W, and 6M-12W.",
	    "",
    "### GSE240547",
    "Metadata was expanded to sample level.",
    "control_diet: Ctrl_1, Ctrl_2, Ctrl_3, Ctrl_4.",
    "hfhs_diet: HFHS_1, HFHS_2, HFHS_3, HFHS_4.",
    "control_diet samples are projection-internal same-study controls.",
    "hfhs_diet samples are diet-intervention projection samples.",
    "HFHS relative delta is computed against the control_diet group.",
    "",
    "### PRJNA1076509",
    "FMT from estropausal female donor microbiota is compared against FMT from young female donor microbiota.",
    "The young female donor microbiota group is the explicit matched-control pool.",
    "Expected n_matched_controls = 7 for every PRJNA1076509 row.",
    "",
    "### GSE228975",
    "Excluded from primary external-projection interpretation because it is outside the current age-bin frame.",
    "Predictions may be retained as all-projection/OOD output only.",
    "No primary matched-control relative delta is computed.",
    "",
    "## Caution",
    "External projections position samples along frozen transcriptomic ageing references. They do not prove rejuvenation, accelerated ageing, or causality without independent functional evidence."
  )
  dir.create(dirname(notes_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, notes_path)
  file.copy(notes_path, file.path(out_dir, basename(notes_path)), overwrite = TRUE)
  invisible(notes_path)
}

invalid_relative_rows <- function(rel) {
  rel[
    (matched_control_method == "no_matched_control" & is.finite(relative_delta_age_week)) |
      (matched_control_method == "no_matched_control" & n_matched_controls != 0) |
      (n_matched_controls == 0 & is.finite(relative_delta_age_week)) |
      (primary_external_projection_include == FALSE & is.finite(relative_delta_age_week)) |
      (matched_control_allowed == FALSE & is.finite(relative_delta_age_week)) |
      (dataset_folder == "GSE228975" & is.finite(relative_delta_age_week)) |
      (dataset_folder == "PRJNA1076509" & matched_control_method != "PRJNA1076509_young_female_donor_microbiota") |
      (dataset_folder == "PRJNA1076509" & n_matched_controls != 7) |
      (dataset_folder == "GSE194196" & primary_external_projection_include == TRUE & matched_control_method != "GSE194196_WT_1_4_7") |
      (dataset_folder == "GSE194196" & primary_external_projection_include == TRUE & n_matched_controls != 3) |
      (dataset_folder == "GSE252791" & primary_external_projection_include == TRUE & matched_control_method != "GSE252791_WT_1_4_7") |
      (dataset_folder == "GSE252791" & primary_external_projection_include == TRUE & n_matched_controls != 3) |
      (dataset_folder == "GSE279573" & projection_role == "external_projected_treatment" & n_matched_controls != 2) |
      (dataset_folder == "GSE279446" & grepl("Zeb1_200H|Zeb1_200M|genetic_model", paste(projection_class, treatment, genotype, sample_title, sep = " "), ignore.case = TRUE) &
        gse279446_cohort == "adult_cohort" & matched_control_method != "GSE279446_adult_WT_ovary_controls") |
      (dataset_folder == "GSE279446" & grepl("Zeb1_200H|Zeb1_200M|genetic_model", paste(projection_class, treatment, genotype, sample_title, sep = " "), ignore.case = TRUE) &
        gse279446_cohort == "young_cohort" & matched_control_method != "GSE279446_young_WT_controls") |
	      (dataset_folder == "GSE279446" & grepl("Zeb1_WT|WT ovary", paste(sample_key, treatment, genotype, sample_title, sep = " "), ignore.case = TRUE) &
	        grepl("heterozygous", projection_class, ignore.case = TRUE)) |
	      (dataset_folder == "GSE329024" & primary_external_projection_include == TRUE &
	        matched_control_method != reference_control_group_id) |
	      (dataset_folder == "GSE240547" &
	        grepl("hfhs|high-fat|high-sucrose", paste(sample_key, treatment, sample_title, sep = " "), ignore.case = TRUE) &
	        matched_control_method != "GSE240547_control_diet") |
      (dataset_folder == "GSE240547" &
        grepl("hfhs|high-fat|high-sucrose", paste(sample_key, treatment, sample_title, sep = " "), ignore.case = TRUE) &
        n_matched_controls != 4) |
      (dataset_folder == "GSE235021" &
        !grepl("polished rice", paste(treatment, sample_title, sep = " "), ignore.case = TRUE) &
        (matched_control_method != reference_control_group_id |
          n_matched_controls != fifelse(age_week < 30, 6L, 5L))) |
      (dataset_folder == "GSE249151" &
        grepl("VCD", paste(treatment, sample_title, sep = " "), ignore.case = TRUE) &
        (matched_control_method != "GSE249151_sesame_oil_controls" | n_matched_controls != 5)) |
      (dataset_folder == "GSE128240" &
        grepl("cyclophosphamide|CTX", paste(treatment, sample_title, sep = " "), ignore.case = TRUE) &
        (matched_control_method != "GSE128240_same_study_controls" | n_matched_controls != 3)) |
      (dataset_folder == "GSE98497" &
        grepl("cKO|knockout", paste(genotype, sample_title, sep = " "), ignore.case = TRUE) &
        (matched_control_method != "GSE98497_WT_controls" | n_matched_controls != 7)) |
      (dataset_folder == "PRJNA787028" &
        grepl("8 mice per cage|8pc", paste(treatment, sample_title, sep = " "), ignore.case = TRUE) &
        (matched_control_method != "PRJNA787028_two_per_cage_controls" | n_matched_controls != 6))
  ]
}

write_summary <- function(pred_all, pred_primary, rel_all, pca_scores, invalid_rel) {
  pred <- pred_primary
  class_summary <- pred[, .(
    n_samples = uniqueN(sample_key),
    median_delta_age_week = median(delta_age_week, na.rm = TRUE),
    mean_delta_age_week = mean(delta_age_week, na.rm = TRUE)
  ), by = .(model_family, projection_class)][order(model_family, mean_delta_age_week)]
  agreement <- dcast(pred, sample_key + projection_class ~ model_family, value.var = "delta_age_week")
  agreement_text <- if (all(c("elasticnet", "pca_glmnet") %in% names(agreement))) {
    ok <- complete.cases(agreement[, .(elasticnet, pca_glmnet)])
    if (sum(ok) >= 3) {
      paste0("ElasticNet/PCA delta-age correlation: ", round(cor(agreement$elasticnet[ok], agreement$pca_glmnet[ok]), 3), ".")
    } else {
      "Too few paired projections to estimate ElasticNet/PCA agreement."
    }
  } else {
    "Both retained clocks were not present in the final projection table."
  }

  pc_contrib_cols <- grep("_contribution_week$", names(pca_scores), value = TRUE)
  top_pc_text <- "PC contribution coefficients were not available."
  if (length(pc_contrib_cols) > 0 && nrow(pca_scores) > 0) {
    pc_abs <- pca_scores[, lapply(.SD, function(x) mean(abs(x), na.rm = TRUE)), .SDcols = pc_contrib_cols]
    ord <- names(sort(unlist(pc_abs), decreasing = TRUE))
    top_pc_text <- paste("Largest average PCA contributions:", paste(head(gsub("_contribution_week$", "", ord), 5), collapse = ", "), ".")
  }

  lines <- c(
    "# External projection summary",
    "",
    "## Run and model scope",
    "",
    "- v1 deployment artifacts were used from `results/model_artifacts/deployment/`.",
    "- No model was refit: projected/reference rows did not modify fitting, scaling, PCA, glmnet, or calibration.",
    paste0("- All projected/reference samples: ", uniqueN(pred_all$sample_key), "."),
    paste0("- Primary projected/reference samples: ", uniqueN(pred_primary$sample_key), "."),
    paste0("- Prediction rows: ", nrow(pred_all), "."),
    paste0("- Relative-control rows: ", nrow(rel_all), "."),
    paste0("- Invalid relative-control rows: ", nrow(invalid_rel), "."),
    paste0("- Datasets: ", uniqueN(pred_all$dataset_folder), "."),
    "",
    "## Delta-age definitions",
    "",
    "- Raw delta age = predicted age - chronological age.",
    "- Matched-control relative delta age = sample raw delta age - mean raw delta age of the matched same-study control group.",
    "- Raw delta age is absolute clock placement; relative delta age is the genotype/intervention/diet/FMT effect relative to matched controls.",
    "",
    "## Dataset-specific matched controls",
    "",
    "- GSE194196: Bnc1tr/tr samples are compared against WT-1/WT-4/WT-7; raw plots include these WT controls as reference rows.",
    "- GSE252791: genetic-model treatment samples are compared against WT1/WT2/WT3 reference controls; raw plots include WT controls.",
	    "- GSE279446: WT, Zeb1_200H, and Zeb1_200M samples from adult and young cohorts are projected; adult and young cohorts are kept separate; Zeb1_WT is treated as WT, not heterozygous.",
	    "- GSE279573: WT chow fertile controls are included as reference-control projection rows; KO chow infertile and KO vegetable-oil-free diet/rescue samples are compared to WT chow fertile controls; n_controls = 2, so interpretation is descriptive.",
	    "- GSE329024: CUMS samples are compared to same-study controls matched by starting age and exposure duration.",
	    "- GSE240547: HFHS diet samples are compared to same-study control_diet projection-internal controls.",
    "- PRJNA1076509: estropausal female donor FMT samples are compared to young female donor FMT controls; n_controls = 7.",
    "- GSE228975: outside the current age-bin frame; retained only as OOD/supplementary all-projection output with no relative delta.",
    "",
    "## Reference-control rows",
    "",
    paste0(
      "- Reference-control rows: ",
      nrow(pred_all[matched_reference_control_flag == TRUE]),
      " sample x model rows across ",
      uniqueN(pred_all[matched_reference_control_flag == TRUE]$sample_key),
      " reference-control samples."
    ),
    "",
    "## Class-level shifts",
    "",
    paste(capture.output(print(class_summary)), collapse = "\n"),
    "",
    "## Clock agreement",
    "",
    agreement_text,
    "",
    "## PCA axes",
    "",
    top_pc_text,
    "",
    "## Caution",
    "",
    "These are model-level transcriptomic projections, not proof of rejuvenation or causality."
  )
  writeLines(lines, file.path(out_dir, "external_projection_summary.md"))
}

message("Auditing locked/deployment artifacts.")
artifact_audit <- audit_artifacts(retained_clocks)
fwrite(artifact_audit$audit, file.path(out_dir, "projection_artifact_audit.tsv"), sep = "\t")
write_missing_and_stop(artifact_audit$audit)

message("Loading metadata and external projection samples.")
meta <- fread(metadata_path, na.strings = c("", "NA", "NaN"))
required_cols <- c("sample_key", "dataset_folder", "modality", "age_week", "expression_logcpm_path")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
meta[, age_week := safe_num(age_week)]
setkey(meta, sample_key)
cra003645_heldout_keys <- meta[
  dataset_folder == "CRA003645" & split_role == "test",
  unique(sample_key)
]
if (length(cra003645_heldout_keys) != 16L) {
  stop(
    "Expected 16 frozen held-out CRA003645 samples, observed ",
    length(cra003645_heldout_keys), ".",
    call. = FALSE
  )
}
proj_meta <- make_projection_metadata(meta)[
  dataset_folder %chin% required_external_projection_datasets &
    !sample_key %chin% cra003645_heldout_keys
]
control_meta <- copy(meta[control_like_clock_candidate == TRUE])
control_meta[, control_like_clock_candidate := as.logical(control_like_clock_candidate)]
matched_control_meta <- make_matched_control_metadata(meta, control_meta)
matched_control_meta <- assign_projection_roles(matched_control_meta, default_role = "control_meta_reference_control")
matched_control_meta <- matched_control_meta[
  !sample_key %chin% cra003645_heldout_keys
]
matched_control_meta[, reference_control_group_label := vapply(reference_control_group_id, reference_control_group_label_one, character(1))]
reference_meta <- control_meta[split_role %chin% c("train", "validation")]
if (nrow(reference_meta) == 0) stop("No train/validation control-like reference samples detected.", call. = FALSE)

reference_control_meta <- matched_control_meta[
  dataset_folder %chin% unique(proj_meta$dataset_folder) &
    dataset_folder != "GSE228975" &
    control_pool_source == "control_meta"
]
projection_apply_meta <- rbindlist(
  list(proj_meta, reference_control_meta),
  use.names = TRUE,
  fill = TRUE
)
setorder(projection_apply_meta, dataset_folder, projection_role, sample_key)
projection_apply_meta <- unique(projection_apply_meta, by = "sample_key")

all_keys <- unique(c(projection_apply_meta$sample_key, matched_control_meta$sample_key))
loaded_expr <- load_expression_for_keys(meta, all_keys)
alignment_issues_out <- loaded_expr$alignment_issues
if (nrow(alignment_issues_out) == 0) {
  alignment_issues_out <- data.table(
    modality = character(),
    dataset_folder = character(),
    sample_key = character(),
    issue = character()
  )
}
fwrite(alignment_issues_out, file.path(out_dir, "external_projection_alignment_issues.tsv"), sep = "\t")
aligned_keys <- loaded_expr$aligned_sample_keys
skipped_projection <- proj_meta[!sample_key %chin% aligned_keys]
skipped_projection_out <- skipped_projection[, .(
  sample_key, dataset_folder, modality, age_week, projection_class, projection_reason
)]
if (nrow(skipped_projection_out) == 0) {
  skipped_projection_out <- data.table(
    sample_key = character(),
    dataset_folder = character(),
    modality = character(),
    age_week = numeric(),
    projection_class = character(),
    projection_reason = character()
  )
}
fwrite(
  skipped_projection_out,
  file.path(out_dir, "external_projection_skipped_unaligned_samples.tsv"),
  sep = "\t"
)
proj_meta <- proj_meta[sample_key %chin% aligned_keys]
projection_apply_meta <- projection_apply_meta[sample_key %chin% aligned_keys]
control_meta <- control_meta[sample_key %chin% aligned_keys]
matched_control_meta <- matched_control_meta[sample_key %chin% aligned_keys]
reference_meta <- reference_meta[sample_key %chin% aligned_keys]
if (nrow(projection_apply_meta) == 0) stop("No projection/reference samples remained after expression alignment.", call. = FALSE)
if (nrow(reference_meta) < 3) stop("Too few aligned reference controls remained for OOD distance calculations.", call. = FALSE)

deployment_training_keys <- unique(unlist(lapply(artifact_audit$objects, function(obj) {
  src <- obj$source_artifact %||% obj
  as.character(src$training_sample_keys %||% src$train_sample_keys %||% character())
}), use.names = FALSE))
projection_apply_meta[, used_in_deployment_training_flag := sample_key %chin% deployment_training_keys]
matched_control_meta[, used_in_deployment_training_flag := sample_key %chin% deployment_training_keys]

prediction_list <- list()
control_prediction_list <- list()
pca_score_list <- list()
reference_pca_score_list <- list()

for (i in seq_len(nrow(retained_clocks))) {
  clock <- retained_clocks[i]
  key <- paste(clock$model_family, clock$selector, sep = "\r")
  obj <- artifact_audit$objects[[key]]
  cal <- extract_calibration(obj)
  message("Projecting ", clock$model_family, " / ", clock$selector, ".")
  res <- clock_predict(clock, obj, cal, loaded_expr$expression_by_piece, projection_apply_meta, reference_meta)
  ctl <- clock_predict(clock, obj, cal, loaded_expr$expression_by_piece, matched_control_meta, reference_meta)
  prediction_list[[key]] <- res$predictions
  control_prediction_list[[key]] <- ctl$predictions[, .(
    sample_key, model_family, selector, predicted_age_week,
    delta_age_week = predicted_age_week - matched_control_meta$age_week[match(sample_key, matched_control_meta$sample_key)],
    projection_distance
  )]
  if (!is.null(res$pca_scores)) pca_score_list[[key]] <- res$pca_scores
  if (!is.null(res$reference_pca_scores)) reference_pca_score_list[[key]] <- res$reference_pca_scores
}

pred <- rbindlist(prediction_list, fill = TRUE)
pred <- merge(pred, projection_apply_meta, by = "sample_key", all.x = TRUE)
pred[, delta_age_week := predicted_age_week - age_week]
pred[, primary_external_projection_include := as.logical(primary_external_projection_include)]
pred[is.na(primary_external_projection_include), primary_external_projection_include := TRUE]
pred[, matched_control_allowed := as.logical(matched_control_allowed)]
pred[is.na(matched_control_allowed), matched_control_allowed := TRUE]
pred[, same_dataset_control_available := dataset_folder %chin% unique(matched_control_meta$dataset_folder)]
pred[, matched_control_available := FALSE]

treatment_cols <- intersect(c("treatment", "intervention", "genotype", "diet", "stress", "microbiome", "sample_context", "sample_title"), names(pred))
base_cols <- c(
  "sample_key", "dataset_folder", "modality", "age_week", "bio_age_bin",
  "projection_role", "projection_class", "projection_reason",
  "matched_reference_control_flag", "matched_reference_control_source",
  "used_in_deployment_training_flag", "primary_external_projection_include",
  "primary_effect_interpretation_include", "matched_control_allowed",
  "reference_control_group_id", "reference_control_group_label",
  "gse279446_cohort", "gse279446_platform", "projection_dataset_note",
  treatment_cols, "model_family", "selector", "raw_predicted_age_week", "predicted_age_week", "delta_age_week",
  "gene_coverage_fraction", "n_model_genes", "n_observed_model_genes", "missing_model_genes",
  "projection_gene_coverage_flag", "projection_distance", "projection_distance_percentile", "OOD_flag"
)
base_cols <- unique(base_cols[base_cols %in% names(pred)])

control_pred <- rbindlist(control_prediction_list, fill = TRUE)
relative_input_cols <- c(
  "sample_key", "model_family", "selector", "dataset_folder", "modality", "age_week",
  "bio_age_bin", "delta_age_week", "projection_role", "reference_control_group_id",
  "primary_external_projection_include", "matched_control_allowed"
)
relative <- merge(
  pred[, ..relative_input_cols],
  compute_relative_to_controls(pred[, ..relative_input_cols], control_pred, matched_control_meta),
  by = c("sample_key", "model_family", "selector", "dataset_folder", "modality", "age_week", "bio_age_bin", "projection_role", "reference_control_group_id"),
  all.x = TRUE
)
relative <- merge(
  relative,
  unique(pred[, c(
    "sample_key", "projection_class", "projection_reason", "matched_reference_control_flag",
    "matched_reference_control_source", "used_in_deployment_training_flag",
    "primary_external_projection_include", "primary_effect_interpretation_include",
    "matched_control_allowed", "reference_control_group_label", "gse279446_cohort",
    "gse279446_platform", treatment_cols, "sample_id", "accession", "projection_dataset_note"
  ), with = FALSE]),
  by = "sample_key",
  all.x = TRUE,
  suffixes = c("", ".pred")
)
for (cc in c("primary_external_projection_include", "matched_control_allowed")) {
  alt <- paste0(cc, ".pred")
  if (alt %in% names(relative)) {
    relative[is.na(get(cc)), (cc) := get(alt)]
    relative[, (alt) := NULL]
  }
}
relative[, primary_external_projection_include := as.logical(primary_external_projection_include)]
relative[, matched_control_allowed := as.logical(matched_control_allowed)]
relative[matched_control_method == "no_matched_control", `:=`(
  relative_delta_age_week = NA_real_,
  relative_delta_age = NA_real_,
  control_expected_delta_age_week = NA_real_,
  control_expected_delta_age = NA_real_,
  control_expected_predicted_age_week = NA_real_,
  control_expected_pred_age = NA_real_,
  n_matched_controls = 0L,
  matched_control_sample_keys = NA_character_,
  matched_control_source = "none"
)]
relative[matched_control_allowed == FALSE | primary_external_projection_include == FALSE, `:=`(
  relative_delta_age_week = NA_real_,
  relative_delta_age = NA_real_,
  control_expected_delta_age_week = NA_real_,
  control_expected_delta_age = NA_real_,
  control_expected_predicted_age_week = NA_real_,
  control_expected_pred_age = NA_real_,
  n_matched_controls = 0L,
  matched_control_sample_keys = NA_character_,
  matched_control_source = "none",
  matched_control_method = "no_matched_control"
)]
relative_all_cols <- c(
  "sample_key", "dataset_folder", "modality", "age_week", "bio_age_bin",
  "projection_role", "projection_class", "projection_reason",
  "matched_reference_control_flag", "matched_reference_control_source",
  "used_in_deployment_training_flag", "primary_external_projection_include",
  "primary_effect_interpretation_include", "matched_control_allowed",
  "reference_control_group_id", "reference_control_group_label",
  "gse279446_cohort", "gse279446_platform",
  treatment_cols, "sample_id", "accession", "model_family", "selector",
  "control_expected_predicted_age_week", "control_expected_delta_age_week",
  "relative_delta_age_week", "control_expected_pred_age", "control_expected_delta_age",
  "relative_delta_age", "matched_control_method", "matched_control_source",
  "n_matched_controls", "matched_control_sample_keys", "relative_control_reliability",
  "relative_control_warning", "projection_dataset_note"
)
relative_all_cols <- unique(relative_all_cols[relative_all_cols %in% names(relative)])
relative_all_out <- relative[, ..relative_all_cols]
relative_primary_out <- relative[primary_external_projection_include == TRUE, ..relative_all_cols]
invalid_rel <- invalid_relative_rows(relative)
fwrite(invalid_rel, file.path(out_dir, "external_projection_unreliable_relative_controls_audit.tsv"), sep = "\t")
if (nrow(invalid_rel) > 0) {
  stop(
    "Invalid matched-control rows detected; see external_projection_unreliable_relative_controls_audit.tsv",
    call. = FALSE
  )
}
fwrite(relative_all_out, file.path(out_dir, "external_projection_relative_to_controls_all.tsv"), sep = "\t")
fwrite(relative_primary_out, file.path(out_dir, "external_projection_relative_to_controls.tsv"), sep = "\t")

pred_match_cols <- c(
  "sample_key", "model_family", "selector", "matched_control_method",
  "matched_control_source", "n_matched_controls", "relative_control_reliability"
)
pred <- merge(pred, unique(relative[, ..pred_match_cols]), by = c("sample_key", "model_family", "selector"), all.x = TRUE)
base_cols <- unique(c(
  base_cols,
  "matched_control_method", "matched_control_source", "n_matched_controls",
  "relative_control_reliability"
))
base_cols <- unique(base_cols[base_cols %in% names(pred)])
pred_all_out <- pred[, ..base_cols]
pred_primary_out <- pred[primary_external_projection_include == TRUE, ..base_cols]
fwrite(pred_all_out, file.path(out_dir, "external_projection_predictions_all.tsv"), sep = "\t")
fwrite(pred_primary_out, file.path(out_dir, "external_projection_predictions.tsv"), sep = "\t")

reference_controls_out <- pred[
  matched_reference_control_flag == TRUE,
  .(
    dataset_folder, sample_key, age_week, model_family, selector,
    projection_role, matched_reference_control_source,
    used_in_deployment_training_flag, predicted_age_week, delta_age_week,
    reference_control_group_id, reference_control_group_label
  )
]
fwrite(reference_controls_out, file.path(out_dir, "external_projection_reference_controls.tsv"), sep = "\t")

control_pool_candidates <- rbindlist(list(
  matched_control_meta[, .(
    dataset_folder, sample_key, sample_id, accession, modality, age_week, bio_age_bin,
    treatment, sample_context, sample_title, projection_role, projection_class,
    matched_reference_control_flag, matched_reference_control_source,
    reference_control_group_id, reference_control_group_label, control_pool_source,
    control_candidate_status, control_candidate_note
  )],
  projection_apply_meta[
    !(sample_key %chin% matched_control_meta$sample_key),
    .(
      dataset_folder, sample_key, sample_id, accession, modality, age_week, bio_age_bin,
      treatment, sample_context, sample_title, projection_role, projection_class,
      matched_reference_control_flag, matched_reference_control_source,
      reference_control_group_id, reference_control_group_label,
      control_pool_source = "none",
      control_candidate_status = "rejected_or_projection_sample",
      control_candidate_note = projection_reason
    )
  ]
), use.names = TRUE, fill = TRUE)
fwrite(control_pool_candidates, file.path(out_dir, "external_projection_control_pool_candidates.tsv"), sep = "\t")

control_matching_audit <- relative[, .(
  dataset_folder, sample_key, model_family, matched_control_method, matched_control_source,
  n_matched_controls, matched_control_sample_keys, relative_control_reliability,
  relative_control_warning
)]
fwrite(control_matching_audit, file.path(out_dir, "external_projection_control_matching_audit.tsv"), sep = "\t")

dataset_notes <- data.table(dataset_folder = required_external_projection_datasets)
dataset_notes[, `:=`(
  expected_relative_method = "same_dataset_control_meta_when_available",
  dataset_note = "Use same-dataset control metadata when available."
)]
dataset_notes[dataset_folder == "CRA017492", `:=`(
  expected_relative_method = "CRA017492_young_control_reference;CRA017492_aged_water_controls",
  dataset_note = "Young controls are a natural-age reference; aged NMN is matched to aged water controls."
)]
dataset_notes[dataset_folder == "GSE194196", `:=`(
  expected_relative_method = "GSE194196_WT_1_4_7",
  dataset_note = "Bnc1tr/tr samples are compared with WT-1/WT-4/WT-7 controls."
)]
dataset_notes[dataset_folder == "GSE252791", `:=`(
  expected_relative_method = "GSE252791_WT_1_4_7",
  dataset_note = "Rnf111 genetic-model samples are compared with WT1/WT2/WT3 controls."
)]
dataset_notes[dataset_folder == "GSE279446", `:=`(
  expected_relative_method = "GSE279446_adult_WT_ovary_controls;GSE279446_young_WT_controls",
  dataset_note = "Adult and young Zeb1 cohorts remain separate and use cohort-matched WT controls."
)]
dataset_notes[dataset_folder == "GSE279573", `:=`(
  expected_relative_method = "GSE279573_WT_chow_fertile_controls",
  dataset_note = "KO groups are compared with two WT chow fertile controls."
)]
dataset_notes[dataset_folder == "GSE329024", `:=`(
  expected_relative_method = "GSE329024_design_matched_controls",
  dataset_note = "CUMS samples use controls matched by starting age and exposure duration."
)]
dataset_notes[dataset_folder == "GSE228975", `:=`(
  expected_relative_method = "no_matched_control",
  dataset_note = "Outside the current age-bin frame and excluded from primary interpretation."
)]
dataset_notes[dataset_folder == "GSE240547", `:=`(
  expected_relative_method = "GSE240547_control_diet",
  dataset_note = "HFHS diet is compared with same-study control-diet samples."
)]
dataset_notes[dataset_folder == "PRJNA1076509", `:=`(
  expected_relative_method = "PRJNA1076509_young_female_donor_microbiota",
  dataset_note = "Estropausal-donor FMT is compared with young-donor FMT."
)]
dataset_notes[dataset_folder == "GSE235021", `:=`(
  expected_relative_method = "GSE235021_polished_rice_6m;GSE235021_polished_rice_9m",
  dataset_note = "Black-rice and white-brown-rice samples use pooled bulk and scRNA-pseudobulk polished-rice controls within age."
)]
dataset_notes[dataset_folder == "GSE154890", `:=`(
  expected_relative_method = "no_matched_control",
  dataset_note = "Natural-age external placement cohort; no treatment-relative effect is calculated."
)]
dataset_notes[dataset_folder == "GSE249151", `:=`(
  expected_relative_method = "GSE249151_sesame_oil_controls",
  dataset_note = "VCD day-15 is age matched; VCD day-30 is retained with an age-confounding warning."
)]
dataset_notes[dataset_folder == "GSE106339", `:=`(
  expected_relative_method = "GSE106339_prepubertal_WT_controls",
  dataset_note = "Three-week GOAT-KO samples are genetic models and WT samples are reference controls; all retain a cross-cutting prepubertal/OOD eligibility flag."
)]
dataset_notes[dataset_folder == "GSE128240", `:=`(
  expected_relative_method = "GSE128240_same_study_controls",
  dataset_note = "Cyclophosphamide samples use three same-study controls."
)]
dataset_notes[dataset_folder == "GSE98497", `:=`(
  expected_relative_method = "GSE98497_WT_controls",
  dataset_note = "Mtor oocyte-conditional KO samples use seven WT controls."
)]
dataset_notes[dataset_folder == "PRJNA787028", `:=`(
  expected_relative_method = "PRJNA787028_two_per_cage_controls",
  dataset_note = "Eight-mice-per-cage stress samples use six two-mice-per-cage controls."
)]
dataset_notes[dataset_folder == "PRJNA1002222", `:=`(
  expected_relative_method = "PRJNA1002222_Cd38KO_vs_WT_8M;PRJNA1002222_Cd38KO_vs_WT_12M;no_matched_control",
  dataset_note = "Cd38-knockout samples use same-age WT controls at 8M and 12M; 2M WT is placement-only."
)]
dataset_availability <- data.table(dataset_folder = required_external_projection_datasets)
dataset_availability <- merge(
  dataset_availability,
  unique(pred[, .(sample_key, dataset_folder, primary_external_projection_include)])[, .(
    n_projected_samples_all = .N,
    n_projected_samples_primary = sum(primary_external_projection_include == TRUE, na.rm = TRUE),
    primary_external_projection_include = any(primary_external_projection_include == TRUE, na.rm = TRUE)
  ), by = dataset_folder],
  by = "dataset_folder",
  all.x = TRUE
)
candidate_counts <- control_pool_candidates[, .(
  n_control_meta_candidates = sum(control_pool_source == "control_meta", na.rm = TRUE),
  n_projection_internal_control_candidates = sum(control_pool_source == "projection_internal_control", na.rm = TRUE),
  has_same_dataset_control_meta = any(control_pool_source == "control_meta", na.rm = TRUE),
  has_projection_internal_control = any(control_pool_source == "projection_internal_control", na.rm = TRUE)
), by = dataset_folder]
dataset_availability <- merge(dataset_availability, candidate_counts, by = "dataset_folder", all.x = TRUE)
observed_methods <- relative[, .(
  observed_relative_method = paste(sort(unique(matched_control_method)), collapse = ";")
), by = dataset_folder]
dataset_availability <- merge(dataset_availability, observed_methods, by = "dataset_folder", all.x = TRUE)
dataset_availability <- merge(dataset_availability, dataset_notes, by = "dataset_folder", all.x = TRUE)
dataset_availability[is.na(n_projected_samples_all), n_projected_samples_all := 0L]
dataset_availability[is.na(n_projected_samples_primary), n_projected_samples_primary := 0L]
for (cc in c("n_control_meta_candidates", "n_projection_internal_control_candidates")) {
  dataset_availability[is.na(get(cc)), (cc) := 0L]
}
for (cc in c("has_same_dataset_control_meta", "has_projection_internal_control", "primary_external_projection_include")) {
  dataset_availability[is.na(get(cc)), (cc) := FALSE]
}
dataset_availability[, recommended_relative_method := fifelse(
  dataset_folder == "GSE228975", "no_matched_control",
	  fifelse(dataset_folder == "CRA017492", "CRA017492_young_control_reference;CRA017492_aged_water_controls",
	  fifelse(dataset_folder == "GSE240547", "GSE240547_control_diet",
	    fifelse(dataset_folder == "PRJNA1076509", "PRJNA1076509_young_female_donor_microbiota",
	      fifelse(dataset_folder == "GSE329024", "GSE329024_design_matched_controls",
	        fifelse(dataset_folder == "GSE194196", "GSE194196_WT_1_4_7",
	          fifelse(dataset_folder == "GSE252791", "GSE252791_WT_1_4_7",
	            fifelse(dataset_folder == "GSE279446", "GSE279446_adult_WT_ovary_controls;GSE279446_young_WT_controls",
	              fifelse(dataset_folder == "GSE279573", "GSE279573_WT_chow_fertile_controls",
	                "same_dataset_control_meta_when_available"
	              )
	            )
	          )
	        )
      )
    )
  ))
)]
dataset_availability[, pass_fail_status := fifelse(
  dataset_folder == "GSE228975" & grepl("no_matched_control", observed_relative_method %||% ""), "PASS",
  fifelse(dataset_folder == "CRA017492" &
            grepl("CRA017492_young_control_reference", observed_relative_method %||% "") &
            grepl("CRA017492_aged_water_controls", observed_relative_method %||% ""), "PASS",
	  fifelse(dataset_folder == "GSE240547" & grepl("GSE240547_control_diet", observed_relative_method %||% ""), "PASS",
	    fifelse(dataset_folder == "PRJNA1076509" & observed_relative_method == "PRJNA1076509_young_female_donor_microbiota", "PASS",
	      fifelse(dataset_folder == "GSE329024" &
	                grepl("GSE329024_2M_8W_control_meta", observed_relative_method %||% "") &
	                grepl("GSE329024_2M_12W_control_meta", observed_relative_method %||% "") &
	                grepl("GSE329024_6M_8W_control_meta", observed_relative_method %||% "") &
	                grepl("GSE329024_6M_12W_control_meta", observed_relative_method %||% ""), "PASS",
	        fifelse(dataset_folder == "GSE194196" & grepl("GSE194196_WT_1_4_7", observed_relative_method %||% ""), "PASS",
	          fifelse(dataset_folder == "GSE252791" & grepl("GSE252791_WT_1_4_7", observed_relative_method %||% ""), "PASS",
	            fifelse(dataset_folder == "GSE279446" &
	                      grepl("GSE279446_adult_WT_ovary_controls", observed_relative_method %||% "") &
	                      grepl("GSE279446_young_WT_controls", observed_relative_method %||% ""), "PASS",
	              fifelse(dataset_folder == "GSE279573" & grepl("GSE279573_WT_chow_fertile_controls", observed_relative_method %||% ""), "PASS",
	                fifelse(!dataset_folder %chin% c("GSE228975", "GSE240547", "PRJNA1076509", "GSE329024", "GSE194196", "GSE252791", "GSE279446", "GSE279573") & n_projected_samples_all > 0, "PASS", "FAIL")
	              )
	            )
	          )
	        )
      )
    ))
  )
)]
dataset_availability[dataset_folder %chin% c(
  "GSE235021", "GSE154890", "GSE249151", "GSE106339",
  "GSE128240", "GSE98497", "PRJNA787028", "PRJNA1002222"
), recommended_relative_method := expected_relative_method]
dataset_availability[dataset_folder == "GSE235021", pass_fail_status := fifelse(
  grepl("GSE235021_polished_rice_6m", observed_relative_method %||% "") &
    grepl("GSE235021_polished_rice_9m", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "GSE154890", pass_fail_status := fifelse(
  n_projected_samples_all == 19L & observed_relative_method == "no_matched_control",
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "PRJNA1002222", pass_fail_status := fifelse(
  n_projected_samples_all == 21L &
    grepl("PRJNA1002222_Cd38KO_vs_WT_8M", observed_relative_method %||% "") &
    grepl("PRJNA1002222_Cd38KO_vs_WT_12M", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "GSE249151", pass_fail_status := fifelse(
  grepl("GSE249151_sesame_oil_controls", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "GSE106339", pass_fail_status := fifelse(
  n_projected_samples_all == 8L & n_projected_samples_primary == 0L &
    grepl("GSE106339_prepubertal_WT_controls", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "GSE128240", pass_fail_status := fifelse(
  grepl("GSE128240_same_study_controls", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "GSE98497", pass_fail_status := fifelse(
  grepl("GSE98497_WT_controls", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
dataset_availability[dataset_folder == "PRJNA787028", pass_fail_status := fifelse(
  grepl("PRJNA787028_two_per_cage_controls", observed_relative_method %||% ""),
  "PASS", "FAIL"
)]
setcolorder(dataset_availability, c(
  "dataset_folder", "n_projected_samples_all", "n_projected_samples_primary",
  "n_control_meta_candidates", "n_projection_internal_control_candidates",
  "has_same_dataset_control_meta", "has_projection_internal_control",
  "recommended_relative_method", "expected_relative_method", "observed_relative_method",
  "primary_external_projection_include", "dataset_note", "pass_fail_status"
))
fwrite(dataset_availability, file.path(out_dir, "external_projection_dataset_control_availability.tsv"), sep = "\t")

active_group_col <- if ("treatment" %in% names(pred) && any(!is.na(pred$treatment) & pred$treatment != "")) "treatment" else "projection_reason"
group_summary <- pred[, .(
  n_samples = uniqueN(sample_key),
  mean_age_week = mean(age_week, na.rm = TRUE),
  mean_predicted_age_week = mean(predicted_age_week, na.rm = TRUE),
  mean_delta_age_week = mean(delta_age_week, na.rm = TRUE),
  median_delta_age_week = median(delta_age_week, na.rm = TRUE)
), by = .(dataset_folder, projection_class, group = get(active_group_col), model_family, selector)]
rel_summary <- relative[, .(
  mean_relative_delta_age_week = mean(relative_delta_age_week, na.rm = TRUE),
  median_relative_delta_age_week = median(relative_delta_age_week, na.rm = TRUE),
  mean_n_matched_controls = mean(n_matched_controls, na.rm = TRUE)
), by = .(dataset_folder, projection_class, group = get(active_group_col), model_family, selector)]
group_summary <- merge(group_summary, rel_summary, by = c("dataset_folder", "projection_class", "group", "model_family", "selector"), all.x = TRUE)
fwrite(group_summary, file.path(out_dir, "external_projection_group_summary.tsv"), sep = "\t")

model_concordance <- dcast(
  pred[, .(sample_key, dataset_folder, projection_class, model_family, delta_age_week)],
  sample_key + dataset_folder + projection_class ~ model_family,
  value.var = "delta_age_week"
)
if (all(c("elasticnet", "pca_glmnet") %in% names(model_concordance))) {
  model_concordance[, `:=`(
    delta_age_difference_pca_minus_elasticnet = pca_glmnet - elasticnet,
    absolute_delta_age_disagreement = abs(pca_glmnet - elasticnet)
  )]
}
fwrite(model_concordance, file.path(out_dir, "external_projection_model_concordance.tsv"), sep = "\t")

gene_coverage <- pred[, .(
  sample_key, dataset_folder, model_family, selector, n_model_genes,
  n_observed_model_genes, gene_coverage_fraction, projection_gene_coverage_flag,
  missing_model_genes
)]
fwrite(gene_coverage, file.path(out_dir, "external_projection_gene_coverage.tsv"), sep = "\t")
qc_flags <- pred[, .(
  sample_key, dataset_folder, model_family, selector, projection_gene_coverage_flag,
  OOD_flag, projection_distance, projection_distance_percentile,
  primary_external_projection_include, matched_control_allowed
)]
fwrite(qc_flags, file.path(out_dir, "external_projection_qc_flags.tsv"), sep = "\t")

pca_scores <- rbindlist(pca_score_list, fill = TRUE)
if (nrow(pca_scores) > 0) {
  fwrite(pca_scores, file.path(out_dir, "pca_projection_scores.tsv"), sep = "\t")
  contrib_cols <- grep("_contribution_week$", names(pca_scores), value = TRUE)
  if (length(contrib_cols) > 0) {
    contrib <- melt(
      pca_scores,
      id.vars = intersect(c("sample_key", "model_family", "selector"), names(pca_scores)),
      measure.vars = contrib_cols,
      variable.name = "pc",
      value.name = "pc_contribution_week"
    )
    contrib[, pc := gsub("_contribution_week$", "", pc)]
    fwrite(contrib, file.path(out_dir, "pca_projection_pc_contributions.tsv"), sep = "\t")
  }
}
reference_pca_scores <- rbindlist(reference_pca_score_list, fill = TRUE)
if (nrow(reference_pca_scores) > 0) {
  reference_pca_scores <- merge(
    reference_pca_scores,
    reference_meta[, .(sample_key, dataset_folder, modality, age_week)],
    by = "sample_key",
    all.x = TRUE
  )
  fwrite(reference_pca_scores, file.path(out_dir, "reference_pca_projection_scores.tsv"), sep = "\t")
}

write_dataset_inclusion_notes()

run_manifest <- data.table(
  field = c(
    "run_label", "results_root", "output_dir", "archive_dir", "n_projected_samples_all",
    "n_projected_samples_primary", "n_prediction_rows_all", "n_relative_control_rows_all",
    "n_invalid_relative_control_rows", "created_at"
  ),
  value = c(
    Sys.getenv("RUN_LABEL", unset = "v1"), results_dir, out_dir, projection_archive_dir %||% NA_character_,
    uniqueN(pred$sample_key), uniqueN(pred[primary_external_projection_include == TRUE]$sample_key),
    nrow(pred), nrow(relative), nrow(invalid_rel), as.character(Sys.time())
  )
)
fwrite(run_manifest, file.path(out_dir, "external_projection_run_manifest.tsv"), sep = "\t")
writeLines(
  c(
    "# External projection run manifest",
    "",
    paste0("- Run label: ", Sys.getenv("RUN_LABEL", unset = "v1")),
    paste0("- Results root: ", results_dir),
    paste0("- Output directory: ", out_dir),
    paste0("- Archive directory: ", projection_archive_dir %||% "none"),
    paste0("- All projected samples: ", uniqueN(pred$sample_key)),
    paste0("- Primary projected samples: ", uniqueN(pred[primary_external_projection_include == TRUE]$sample_key)),
    paste0("- Prediction rows: ", nrow(pred)),
    paste0("- Relative-control rows: ", nrow(relative)),
    paste0("- Invalid relative-control rows: ", nrow(invalid_rel))
  ),
  file.path(out_dir, "external_projection_run_manifest.md")
)
warnings_lines <- unique(na.omit(relative$relative_control_warning))
if (length(warnings_lines) == 0) warnings_lines <- "No external projection warnings."
writeLines(warnings_lines, file.path(out_dir, "external_projection_warnings.txt"))

plot_pred <- pred[primary_external_projection_include == TRUE & projection_gene_coverage_flag != "low_coverage"]
plot_pred <- add_projection_class_label(plot_pred)
if (nrow(plot_pred) > 0) {
  p <- ggplot(plot_pred, aes(age_week, predicted_age_week, color = projection_class_label, shape = modality)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey45") +
    geom_point(size = 2.4, alpha = 0.85) +
    facet_wrap(~model_family) +
    labs(x = "Chronological age (weeks)", y = "Projected transcriptomic age (weeks)", color = "Projection class", shape = "Modality") +
    theme_projection()
  save_plot(p, "01_external_projection_predicted_vs_chronological_age")

  p <- ggplot(plot_pred, aes(projection_class_label, delta_age_week, color = projection_class_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.15) +
    geom_jitter(width = 0.18, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~model_family, scales = "free_x") +
    labs(x = NULL, y = "Delta age (projected - chronological weeks)", color = "Projection class") +
    theme_projection()
  save_plot(p, "02_delta_age_by_projection_class", width = 11, height = 7)
}

plot_rel <- merge(plot_pred, relative[, .(sample_key, model_family, selector, relative_delta_age, matched_control_method)], by = c("sample_key", "model_family", "selector"), all.x = TRUE)
plot_rel <- plot_rel[is.finite(relative_delta_age)]
if (nrow(plot_rel) > 0) {
  p <- ggplot(plot_rel, aes(projection_class_label, relative_delta_age, color = projection_class_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.15) +
    geom_jitter(width = 0.18, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~model_family, scales = "free_x") +
    labs(x = NULL, y = "Matched-control relative delta age (weeks)", color = "Projection class") +
    theme_projection()
  save_plot(p, "07_matched_control_relative_delta", width = 11, height = 7)
}

if (length(treatment_cols) > 0 && nrow(plot_rel) > 0) {
  treatment_col <- treatment_cols[[1]]
  dt <- plot_rel[!is.na(get(treatment_col)) & same_dataset_control_available == TRUE]
  if (nrow(dt) > 0) {
    p <- ggplot(dt, aes(get(treatment_col), delta_age_week, color = projection_class_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
      geom_point(position = position_jitter(width = 0.15, height = 0), alpha = 0.8, size = 1.8) +
      facet_grid(model_family ~ dataset_folder, scales = "free_x", space = "free_x") +
      labs(x = treatment_col, y = "Delta age (weeks)", color = "Projection class") +
      theme_projection()
    save_plot(p, "03_delta_age_by_treatment_within_dataset", width = 14, height = 8)
  }
}

if (nrow(pca_scores) > 0 && all(c("PC1", "PC2") %in% names(pca_scores))) {
  ref_scores <- rbindlist(reference_pca_score_list, fill = TRUE)
  if (nrow(ref_scores) > 0 && all(c("PC1", "PC2") %in% names(ref_scores))) {
    ref_scores <- merge(ref_scores, reference_meta[, .(sample_key, modality)], by = "sample_key", all.x = TRUE)
    proj_scores <- merge(pca_scores, projection_apply_meta[, .(sample_key, projection_role, projection_class, modality, sample_context, treatment)], by = "sample_key", all.x = TRUE)
    proj_scores <- add_projection_class_label(proj_scores)
    p <- ggplot() +
      geom_point(data = ref_scores, aes(PC1, PC2, shape = modality), color = "grey70", alpha = 0.55, size = 1.8) +
      geom_point(data = proj_scores, aes(PC1, PC2, color = projection_class_label, shape = modality), size = 2.4, alpha = 0.9) +
      labs(x = "Frozen PC1 score", y = "Frozen PC2 score", color = "Projection class", shape = "Modality") +
      theme_projection()
    save_plot(p, "04_pca_pc1_pc2_projection_map", width = 9, height = 7)
  }

  contrib_cols <- grep("_contribution_week$", names(pca_scores), value = TRUE)
  if (length(contrib_cols) > 0) {
    heat_input <- merge(pca_scores, projection_apply_meta[, .(sample_key, projection_role, projection_class)], by = "sample_key", all.x = TRUE)
    heat_input <- add_projection_class_label(heat_input)
    heat <- melt(
      heat_input,
      id.vars = c("sample_key", "projection_role", "projection_class", "projection_class_label"),
      measure.vars = contrib_cols,
      variable.name = "PC",
      value.name = "contribution_week"
    )
    heat[, PC := gsub("_contribution_week$", "", PC)]
    p <- ggplot(heat, aes(PC, sample_key, fill = contribution_week)) +
      geom_tile() +
      scale_fill_gradient2(low = "#3B6EA8", mid = "white", high = "#C44E52", midpoint = 0) +
      facet_grid(projection_class_label ~ ., scales = "free_y", space = "free_y") +
      labs(x = "Frozen selected PC", y = NULL, fill = "Weeks") +
      theme_projection() +
      theme(axis.text.y = element_text(size = 5))
    save_plot(p, "05_pca_selected_pc_contribution_heatmap", width = 9, height = 10)
  }
}

if (nrow(pred) > 0) {
  p <- ggplot(pred, aes(gene_coverage_fraction, projection_distance, color = OOD_flag, shape = projection_gene_coverage_flag)) +
    geom_vline(xintercept = 0.80, linetype = "dashed", color = "grey45") +
    geom_point(alpha = 0.85, size = 2.2) +
    facet_wrap(~model_family, scales = "free_y") +
    labs(x = "Model-gene coverage fraction", y = "Projection distance", color = "OOD flag", shape = "Coverage flag") +
    theme_projection()
  save_plot(p, "06_ood_gene_coverage")
}

write_summary(pred, pred[primary_external_projection_include == TRUE], relative, pca_scores, invalid_rel)
cat("All projected sample count:", uniqueN(pred$sample_key), "\n")
cat("Primary projected sample count:", uniqueN(pred[primary_external_projection_include == TRUE]$sample_key), "\n")
cat("Total prediction rows:", nrow(pred), "\n")
cat("Total relative-control rows:", nrow(relative), "\n")
cat("Invalid relative-control rows:", nrow(invalid_rel), "\n")
cat(
  "GSE240547 matched-control status:",
  dataset_availability[dataset_folder == "GSE240547", pass_fail_status],
  dataset_availability[dataset_folder == "GSE240547", observed_relative_method],
  "\n"
)
cat(
  "PRJNA1076509 matched-control status:",
  dataset_availability[dataset_folder == "PRJNA1076509", pass_fail_status],
  dataset_availability[dataset_folder == "PRJNA1076509", observed_relative_method],
  "\n"
)
cat(
  "GSE228975 OOD/exclusion status:",
  dataset_availability[dataset_folder == "GSE228975", pass_fail_status],
  dataset_availability[dataset_folder == "GSE228975", dataset_note],
  "\n"
)
cat("Summary markdown:", file.path(out_dir, "external_projection_summary.md"), "\n")
message("External projection complete: ", out_dir)
