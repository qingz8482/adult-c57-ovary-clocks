#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

bool_env <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = if (default) "TRUE" else "FALSE")
  toupper(trimws(value)) %in% c("TRUE", "T", "1", "YES", "Y")
}

int_env <- function(name, default) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(value)) default else value
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
results_dir <- Sys.getenv("CLOCK_RESULTS_DIR", unset = file.path(project_dir, "results"))
results_dir <- normalizePath(results_dir, mustWork = FALSE)

out_dir <- file.path(results_dir, "integrated_diagnostics", "12_PCA deployment model GSEA")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

diag_dir <- file.path(results_dir, "integrated_diagnostics", "10_pca_diagnostics")
pc_summary_path <- file.path(diag_dir, "final_pca_clock_selected_pc_summary.tsv")
loadings_path <- file.path(diag_dir, "final_pca_clock_gene_pc_loadings_with_feature_scores.tsv")
locked_path <- file.path(results_dir, "locked_model_selection.tsv")
if (!file.exists(locked_path)) stop("Missing locked model selection table: ", locked_path)
locked <- fread(locked_path, sep = "\t")
pca_selector <- locked[model_family == "pca_glmnet", selector][1]
if (is.na(pca_selector) || !nzchar(pca_selector)) stop("Could not determine locked PCA/glmnet selector.")
deployment_artifact_path <- file.path(
  results_dir,
  "model_artifacts",
  "deployment",
  paste0("integrated_bulk_pseudobulk__pca_glmnet__", pca_selector, "__deployment_all_control_like.rds")
)
coefficient_fallback_path <- file.path(results_dir, "model_artifacts", "pca_calibrated_pc_coefficients.tsv")

required_files <- c(pc_summary_path, loadings_path, deployment_artifact_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"))
}

deployment_scope <- "deployment_all_control_like_external_projection"
artifact_scope_short <- "deployment_all_control_like"
seed <- int_env("GSEA_SEED", 20260630)
n_perm <- int_env("N_PERM", 500)
n_perm_refined <- int_env("N_PERM_REFINED", 10000)
run_refined <- bool_env("RUN_REFINED_PERMUTATIONS", FALSE)
run_all_collections <- bool_env("RUN_ALL_COLLECTIONS", FALSE)
min_size <- int_env("MIN_SIZE_IN_UNIVERSE", 10)
max_size <- int_env("MAX_SIZE_IN_UNIVERSE", 300)
max_size_sensitivity <- int_env("MAX_SIZE_SENSITIVITY", 500)
top_n_plot <- int_env("TOP_N_PLOT", 8)
score_type <- 1

if (n_perm < 10) stop("N_PERM must be at least 10.")
if (min_size < 1 || max_size <= min_size) stop("Invalid gene-set size filters.")

temp_r_lib <- Sys.getenv("ADULT_C57_GSEA_R_LIB", unset = "/tmp/adult_c57_r_gsea_lib")
if (dir.exists(temp_r_lib)) {
  .libPaths(c(temp_r_lib, .libPaths()))
}

clean_gene <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "NaN")] <- NA_character_
  x
}

artifact <- readRDS(deployment_artifact_path)
identity_checks <- list(
  training_scope = artifact$training_scope %||% NA_character_,
  artifact_scope = artifact$artifact_scope %||% NA_character_,
  model_family = artifact$model_family %||% NA_character_,
  selector = artifact$selector %||% NA_character_
)
if (!identical(identity_checks$training_scope, "integrated_bulk_pseudobulk") ||
  !identical(identity_checks$artifact_scope, artifact_scope_short) ||
  !identical(identity_checks$model_family, "pca_glmnet") ||
  !identical(identity_checks$selector, pca_selector)) {
  stop(
    "Deployment artifact identity check failed. Found: ",
    paste(names(identity_checks), unlist(identity_checks), sep = "=", collapse = ", ")
  )
}

selected_genes <- clean_gene(artifact$selected_genes)
selected_genes <- selected_genes[!is.na(selected_genes)]
if (length(selected_genes) < 1900) {
  stop("Fewer than 1900 selected PCA genes found in deployment artifact: ", length(selected_genes))
}
if (anyDuplicated(selected_genes)) {
  stop("Deployment artifact selected_genes contains duplicated gene names after cleaning.")
}
gene_universe <- selected_genes

artifact_pcs <- artifact$model$selected_pc_names %||% character()
if (length(artifact_pcs) == 0) {
  stop("Deployment artifact selected PCs are missing.")
}
expected_pcs <- artifact_pcs

pc_coef <- as.data.table(artifact$model$pc_coefficient_table %||% data.frame())
if (nrow(pc_coef) == 0 || !"calibrated_pc_coefficient" %in% names(pc_coef)) {
  if (!file.exists(coefficient_fallback_path)) {
    stop("Calibrated PC coefficients are missing from artifact and fallback table is absent.")
  }
  pc_coef <- fread(coefficient_fallback_path)
  pc_coef <- pc_coef[artifact_scope == artifact_scope_short]
}
if (!all(c("pc", "calibrated_pc_coefficient") %in% names(pc_coef))) {
  stop("Calibrated PC coefficient table lacks required columns pc and calibrated_pc_coefficient.")
}
pc_coef <- pc_coef[pc %in% expected_pcs, .(
  pc = as.character(pc),
  calibrated_pc_coefficient = as.numeric(calibrated_pc_coefficient)
)]
if (!identical(pc_coef$pc[order(match(pc_coef$pc, expected_pcs))], expected_pcs) ||
  any(!is.finite(pc_coef$calibrated_pc_coefficient))) {
  stop("Calibrated PC coefficients are missing for one or more selected deployment PCs: ", paste(expected_pcs, collapse = ", "))
}

pc_summary <- fread(pc_summary_path)
if (!"artifact_scope" %in% names(pc_summary)) stop("PC summary lacks artifact_scope.")
pc_summary <- pc_summary[artifact_scope == deployment_scope]
if (nrow(pc_summary) == 0) stop("No deployment rows found in PC summary.")
if (!setequal(pc_summary$pc, expected_pcs) || uniqueN(pc_summary$pc) != length(expected_pcs)) {
  stop("Deployment PC summary must contain exactly the selected deployment PCs: ", paste(expected_pcs, collapse = ", "), ". Found: ", paste(sort(unique(pc_summary$pc)), collapse = ", "))
}
pc_summary <- merge(
  pc_summary,
  pc_coef,
  by = "pc",
  all.x = TRUE,
  suffixes = c("_summary", "")
)
if ("calibrated_pc_coefficient_summary" %in% names(pc_summary)) {
  pc_summary[, calibrated_pc_coefficient_summary := NULL]
}
pc_summary[, pc := factor(as.character(pc), levels = expected_pcs)]
setorder(pc_summary, pc)
pc_summary[, pc := as.character(pc)]

loadings <- fread(loadings_path)
if (!all(c("artifact_scope", "gene_id", "pc", "pc_loading") %in% names(loadings))) {
  stop("PC loading table lacks one or more required columns: artifact_scope, gene_id, pc, pc_loading.")
}
loadings <- loadings[artifact_scope == deployment_scope]
if (nrow(loadings) == 0) stop("No deployment rows found in loading table.")
loadings[, gene_id := clean_gene(gene_id)]
loadings <- loadings[!is.na(gene_id)]
if (!setequal(loadings$pc, expected_pcs) || uniqueN(loadings$pc) != length(expected_pcs)) {
  stop("Deployment loading table must contain exactly the selected deployment PCs: ", paste(expected_pcs, collapse = ", "), ". Found: ", paste(sort(unique(loadings$pc)), collapse = ", "))
}
if (loadings[, anyDuplicated(gene_id) > 0, by = pc][V1 == TRUE, .N] > 0) {
  dup_pcs <- loadings[, anyDuplicated(gene_id) > 0, by = pc][V1 == TRUE, pc]
  stop("Rank vectors contain duplicated gene names after cleaning for: ", paste(dup_pcs, collapse = ", "))
}
loadings <- loadings[gene_id %in% gene_universe]
if (loadings[, uniqueN(gene_id)] < 1900) {
  stop("Fewer than 1900 selected PCA genes found in deployment loading table after universe filtering.")
}

rank_vectors <- merge(
  loadings[, .(gene_id, pc, pc_loading)],
  pc_summary[, .(
    pc,
    pc_variance_explained = variance_explained,
    pc_cumulative_variance = cumulative_variance_through_pc,
    pc_pearson_age_cor = pc_age_pearson,
    pc_spearman_age_cor = pc_age_spearman,
    clock_direction = clock_coefficient_direction,
    calibrated_pc_coefficient
  )],
  by = "pc",
  all.x = TRUE
)
rank_vectors[, rank_score := as.numeric(pc_loading) * as.numeric(calibrated_pc_coefficient)]
rank_vectors[, rank_direction := fifelse(rank_score >= 0, "older_predictive", "younger_predictive")]
rank_vectors[, abs_rank_score := abs(rank_score)]
rank_vectors[, rank_order := frank(-rank_score, ties.method = "first"), by = pc]
rank_vectors[, abs_rank_order := frank(-abs_rank_score, ties.method = "first"), by = pc]
if (rank_vectors[!is.finite(rank_score), .N] > 0) stop("Rank vectors contain non-finite rank scores.")
if (rank_vectors[, uniqueN(gene_id), by = pc][V1 < 1900, .N] > 0) {
  stop("At least one PC has fewer than 1900 ranked selected genes.")
}

integrated_rank_vectors <- rank_vectors[, {
  abs_contrib <- abs(rank_score)
  dominant_i <- which.max(abs_contrib)
  data.table(
    pc = "integrated_PCA_clock",
    rank_score = sum(rank_score),
    abs_rank_score = abs(sum(rank_score)),
    abs_sum_rank_score = sum(abs_contrib),
    n_pcs_used = .N,
    dominant_pc = pc[[dominant_i]],
    dominant_pc_contribution = rank_score[[dominant_i]],
    dominant_pc_abs_contribution = abs_contrib[[dominant_i]],
    dominant_pc_share = if (sum(abs_contrib) > 0) abs_contrib[[dominant_i]] / sum(abs_contrib) else NA_real_
  )
}, by = gene_id]
integrated_rank_vectors[, rank_direction := fifelse(rank_score >= 0, "older_predictive", "younger_predictive")]
integrated_rank_vectors[, rank_order := frank(-rank_score, ties.method = "first")]
integrated_rank_vectors[, abs_rank_order := frank(-abs_rank_score, ties.method = "first")]
if (integrated_rank_vectors[!is.finite(rank_score), .N] > 0) stop("Integrated rank vector contains non-finite rank scores.")
if (integrated_rank_vectors[, uniqueN(gene_id)] < 1900) stop("Integrated rank vector has fewer than 1900 ranked selected genes.")

load_msigdb_mouse <- function() {
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    message("Loading mouse-native MSigDB through msigdbr.")
    msig <- as.data.table(msigdbr::msigdbr(db_species = "MM", species = "Mus musculus"))
    gene_col <- if ("gene_symbol" %in% names(msig)) "gene_symbol" else "db_gene_symbol"
    needed <- c(gene_col, "gs_name", "gs_collection", "gs_subcollection")
    if (!all(needed %in% names(msig))) {
      stop("msigdbr returned an unexpected schema; missing: ", paste(setdiff(needed, names(msig)), collapse = ", "))
    }
    if (!"gs_collection_name" %in% names(msig)) msig[, gs_collection_name := gs_collection]
    if (!"gs_description" %in% names(msig)) msig[, gs_description := ""]
    if (!"db_version" %in% names(msig)) msig[, db_version := NA_character_]
    out <- msig[, .(
      gene_set_name = as.character(gs_name),
      gene_symbol = clean_gene(get(gene_col)),
      collection = fifelse(is.na(gs_collection), "", as.character(gs_collection)),
      subcollection = fifelse(is.na(gs_subcollection), "", as.character(gs_subcollection)),
      collection_name = fifelse(is.na(gs_collection_name), "", as.character(gs_collection_name)),
      gene_set_description = fifelse(is.na(gs_description), "", as.character(gs_description)),
      msigdb_version = fifelse(is.na(db_version), "", as.character(db_version)),
      gene_set_source = "msigdbr_mouse_native"
    )]
    return(unique(out[!is.na(gene_symbol)]))
  }

  gmt_path <- Sys.getenv("MSIGDB_MOUSE_GMT", unset = "")
  if (!nzchar(gmt_path) || !file.exists(gmt_path)) {
    stop(
      "Mouse-native MSigDB gene sets are unavailable. Install msigdbr for this R library ",
      "or set MSIGDB_MOUSE_GMT to a local mouse-native GMT file."
    )
  }
  message("Loading mouse-native MSigDB from GMT: ", gmt_path)
  lines <- readLines(gmt_path, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) stop("MSIGDB_MOUSE_GMT is empty: ", gmt_path)
  base <- basename(gmt_path)
  infer_from_file <- function() {
    b <- tolower(base)
    if (grepl("(^|[._-])mh([._-]|$)", b)) return(c("MH", ""))
    if (grepl("m2.*cp", b)) return(c("M2", "CP"))
    if (grepl("m5.*go.*bp", b)) return(c("M5", "GO:BP"))
    if (grepl("(^|[._-])m7([._-]|$)", b)) return(c("M7", ""))
    c("", "")
  }
  file_collection <- infer_from_file()
  rows <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(parts) < 3) return(NULL)
    name <- parts[[1]]
    description <- parts[[2]]
    genes <- clean_gene(parts[-c(1, 2)])
    genes <- unique(genes[!is.na(genes)])
    collection <- file_collection[[1]]
    subcollection <- file_collection[[2]]
    if (!nzchar(collection)) {
      if (grepl("^HALLMARK_", name)) {
        collection <- "MH"
      } else if (grepl("^GOBP_", name)) {
        collection <- "M5"; subcollection <- "GO:BP"
      } else if (grepl("^GOCC_", name)) {
        collection <- "M5"; subcollection <- "GO:CC"
      } else if (grepl("^GOMF_", name)) {
        collection <- "M5"; subcollection <- "GO:MF"
      } else if (grepl("^(BIOCARTA|KEGG|REACTOME|WP|PID)_", name)) {
        collection <- "M2"; subcollection <- "CP"
      } else {
        collection <- "GMT_UNSPECIFIED"; subcollection <- ""
      }
    }
    data.table(
      gene_set_name = name,
      gene_symbol = genes,
      collection = collection,
      subcollection = subcollection,
      collection_name = collection,
      gene_set_description = description,
      msigdb_version = paste0("GMT:", base),
      gene_set_source = "MSIGDB_MOUSE_GMT"
    )
  })
  unique(rbindlist(rows, fill = TRUE))
}

msig <- load_msigdb_mouse()
if (nrow(msig) == 0) stop("Mouse-native MSigDB table loaded zero gene-set rows.")

msig[, gene_symbol := clean_gene(gene_symbol)]
msig <- unique(msig[!is.na(gene_symbol)])
msig[, collection := fifelse(is.na(collection) | collection == "", "UNSPECIFIED", collection)]
msig[, subcollection := fifelse(is.na(subcollection), "", subcollection)]
msig[, collection_label := fifelse(subcollection == "", collection, paste(collection, subcollection, sep = ":"))]

available_collections <- sort(unique(msig$collection_label))
if (run_all_collections) {
  msig_selected <- copy(msig)
  collection_selection_note <- "RUN_ALL_COLLECTIONS=TRUE; all available mouse-native MSigDB collections were used."
} else {
  msig_selected <- msig[
    collection == "MH" |
      (collection == "M2" & grepl("^CP", subcollection)) |
      (collection == "M5" & subcollection == "GO:BP") |
      collection == "M7"
  ]
  collection_selection_note <- "Primary collection selection: MH/Hallmark, M2:CP canonical pathways, M5:GO:BP, and M7 where available."
}
if (nrow(msig_selected) == 0) {
  stop(
    "No requested mouse-native MSigDB collections were available. Available collection labels: ",
    paste(available_collections, collapse = ", ")
  )
}

gene_set_original <- msig_selected[, .(
  gene_set_size_original = uniqueN(gene_symbol),
  collection_name = collection_name[which.max(nzchar(collection_name))] %||% "",
  gene_set_description = gene_set_description[which.max(nzchar(gene_set_description))] %||% "",
  msigdb_version = msigdb_version[which.max(nzchar(msigdb_version))] %||% "",
  gene_set_source = gene_set_source[1]
), by = .(collection, subcollection, collection_label, gene_set_name)]
gene_set_original[, gene_set_id := paste(collection, subcollection, gene_set_name, sep = "\t")]

gene_sets <- msig_selected[gene_symbol %in% gene_universe, .(
  genes = list(sort(unique(gene_symbol))),
  gene_set_size_in_universe = uniqueN(gene_symbol)
), by = .(collection, subcollection, collection_label, gene_set_name)]
gene_sets <- merge(gene_set_original, gene_sets, by = c("collection", "subcollection", "collection_label", "gene_set_name"), all.x = TRUE)
gene_sets[is.na(gene_set_size_in_universe), gene_set_size_in_universe := 0L]
gene_sets[, primary_eligible := gene_set_size_in_universe >= min_size & gene_set_size_in_universe <= max_size]
gene_sets[, sensitivity_eligible_max500 := gene_set_size_in_universe >= min_size & gene_set_size_in_universe <= max_size_sensitivity]
gene_sets_primary <- gene_sets[primary_eligible == TRUE]
if (nrow(gene_sets_primary) == 0) {
  stop("No gene sets remain after intersecting with the 2000-gene universe and applying size filters.")
}

calc_es_detail <- function(stats, genes) {
  genes <- intersect(genes, names(stats))
  n <- length(stats)
  nh <- length(genes)
  if (nh == 0 || nh >= n) {
    return(list(es = NA_real_, leading_edge = character(), peak_rank = NA_integer_))
  }
  pos <- sort(match(genes, names(stats)))
  weights <- abs(stats[pos])^score_type
  if (!is.finite(sum(weights)) || sum(weights) == 0) weights <- rep(1, length(pos))
  hit_increment <- weights / sum(weights)
  nm <- n - nh
  hit_index <- seq_along(pos)
  cum_hit <- cumsum(hit_increment)
  running_after_hit <- cum_hit - (pos - hit_index) / nm
  running_before_hit <- c(0, head(cum_hit, -1)) - (pos - hit_index) / nm
  max_i <- which.max(running_after_hit)
  min_i <- which.min(running_before_hit)
  max_es <- running_after_hit[[max_i]]
  min_es <- running_before_hit[[min_i]]
  if (abs(max_es) >= abs(min_es)) {
    peak_rank <- pos[[max_i]]
    leading_edge <- names(stats)[pos[pos <= peak_rank]]
    es <- max_es
  } else {
    peak_rank <- pos[[min_i]]
    leading_edge <- names(stats)[pos[pos >= peak_rank]]
    es <- min_es
  }
  list(es = es, leading_edge = leading_edge, peak_rank = peak_rank)
}

calc_es_only_from_positions <- function(stats_abs, pos) {
  n <- length(stats_abs)
  pos <- sort(pos)
  nh <- length(pos)
  if (nh == 0 || nh >= n) return(NA_real_)
  weights <- stats_abs[pos]^score_type
  if (!is.finite(sum(weights)) || sum(weights) == 0) weights <- rep(1, length(pos))
  hit_increment <- weights / sum(weights)
  nm <- n - nh
  hit_index <- seq_along(pos)
  cum_hit <- cumsum(hit_increment)
  running_after_hit <- cum_hit - (pos - hit_index) / nm
  running_before_hit <- c(0, head(cum_hit, -1)) - (pos - hit_index) / nm
  max_es <- max(running_after_hit)
  min_es <- min(running_before_hit)
  if (abs(max_es) >= abs(min_es)) max_es else min_es
}

same_direction_nes <- function(es, null_es) {
  if (es > 0) {
    same <- null_es[null_es > 0]
  } else if (es < 0) {
    same <- null_es[null_es < 0]
  } else {
    same <- numeric()
  }
  denom <- if (length(same) >= 5) mean(abs(same), na.rm = TRUE) else mean(abs(null_es), na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  es / denom
}

empirical_p_value <- function(es, null_es, n_perm_used) {
  if (es > 0) {
    (1 + sum(null_es >= es, na.rm = TRUE)) / (1 + n_perm_used)
  } else if (es < 0) {
    (1 + sum(null_es <= es, na.rm = TRUE)) / (1 + n_perm_used)
  } else {
    1
  }
}

run_gsea_for_pc <- function(pc_rank, terms, n_perm_used, stage, seed_offset) {
  stats <- pc_rank$rank_score
  names(stats) <- pc_rank$gene_id
  stats <- stats[is.finite(stats)]
  stats <- stats[!duplicated(names(stats))]
  stats <- sort(stats, decreasing = TRUE)
  stats_abs <- abs(stats)
  n <- length(stats)
  if (n < 1900) stop("Rank vector has fewer than 1900 finite genes for ", pc_rank$pc[[1]])

  observed <- rbindlist(lapply(seq_len(nrow(terms)), function(i) {
    genes <- intersect(terms$genes[[i]], names(stats))
    set_size <- length(genes)
    if (set_size < min_size || set_size > max_size || set_size >= n) return(NULL)
    detail <- calc_es_detail(stats, genes)
    if (!is.finite(detail$es)) return(NULL)
    data.table(
      pc = pc_rank$pc[[1]],
      gene_set_id = terms$gene_set_id[[i]],
      collection = terms$collection[[i]],
      subcollection = terms$subcollection[[i]],
      collection_label = terms$collection_label[[i]],
      gene_set_name = terms$gene_set_name[[i]],
      gene_set_size_original = terms$gene_set_size_original[[i]],
      gene_set_size_in_universe = set_size,
      ES = detail$es,
      peak_rank = detail$peak_rank,
      n_leading_edge = length(detail$leading_edge),
      leading_edge_genes = paste(detail$leading_edge, collapse = ";")
    )
  }), fill = TRUE)
  if (nrow(observed) == 0) return(list(results = data.table(), null_qc = data.table()))

  set.seed(seed + seed_offset)
  unique_sizes <- sort(unique(observed$gene_set_size_in_universe))
  null_by_size <- vector("list", length(unique_sizes))
  names(null_by_size) <- as.character(unique_sizes)
  null_qc <- rbindlist(lapply(unique_sizes, function(sz) {
    null_es <- replicate(n_perm_used, calc_es_only_from_positions(stats_abs, sample.int(n, sz)))
    null_by_size[[as.character(sz)]] <<- null_es
    data.table(
      pc = pc_rank$pc[[1]],
      permutation_stage = stage,
      gene_set_size_in_universe = sz,
      n_permutations = n_perm_used,
      n_positive_null = sum(null_es > 0, na.rm = TRUE),
      n_negative_null = sum(null_es < 0, na.rm = TRUE),
      mean_abs_positive_null = mean(abs(null_es[null_es > 0]), na.rm = TRUE),
      mean_abs_negative_null = mean(abs(null_es[null_es < 0]), na.rm = TRUE),
      mean_abs_all_null = mean(abs(null_es), na.rm = TRUE),
      min_null_ES = min(null_es, na.rm = TRUE),
      max_null_ES = max(null_es, na.rm = TRUE)
    )
  }), fill = TRUE)

  observed[, `:=`(
    NES = mapply(function(es, sz) same_direction_nes(es, null_by_size[[as.character(sz)]]), ES, gene_set_size_in_universe),
    empirical_p = mapply(function(es, sz) empirical_p_value(es, null_by_size[[as.character(sz)]], n_perm_used), ES, gene_set_size_in_universe),
    n_permutations = n_perm_used,
    permutation_stage = stage
  )]
  list(results = observed, null_qc = null_qc)
}

message("Running first-pass PC-level GSEA with ", n_perm, " permutations.")
first_pass <- lapply(seq_along(expected_pcs), function(i) {
  pc_name <- expected_pcs[[i]]
  run_gsea_for_pc(rank_vectors[pc == pc_name], gene_sets_primary, n_perm, "first_pass", i * 100000)
})
all_results <- rbindlist(lapply(first_pass, `[[`, "results"), fill = TRUE)
null_qc <- rbindlist(lapply(first_pass, `[[`, "null_qc"), fill = TRUE)

if (nrow(all_results) == 0) stop("GSEA returned zero tested gene sets.")
all_results <- merge(
  all_results,
  pc_summary[, .(
    pc,
    pc_variance_explained = variance_explained,
    pc_cumulative_variance = cumulative_variance_through_pc,
    pc_pearson_age_cor = pc_age_pearson,
    pc_spearman_age_cor = pc_age_spearman,
    clock_direction = clock_coefficient_direction,
    calibrated_pc_coefficient
  )],
  by = "pc",
  all.x = TRUE
)
all_results[, direction := fifelse(NES > 0, "older_aligned", "younger_aligned")]

if (run_refined) {
  message("Running refined GSEA for candidate terms with ", n_perm_refined, " permutations.")
  all_results[, abs_NES_rank := frank(-abs(NES), ties.method = "first"), by = pc]
  candidates <- all_results[empirical_p <= 0.05 | abs_NES_rank <= 100, unique(gene_set_id), by = pc]
  refined_list <- lapply(seq_len(nrow(candidates)), function(i) {
    pc_name <- candidates$pc[[i]]
    ids <- candidates$V1[[i]]
    run_gsea_for_pc(
      rank_vectors[pc == pc_name],
      gene_sets_primary[gene_set_id %in% ids],
      n_perm_refined,
      "refined",
      9000000 + match(pc_name, expected_pcs) * 100000
    )
  })
  refined_results <- rbindlist(lapply(refined_list, `[[`, "results"), fill = TRUE)
  refined_null_qc <- rbindlist(lapply(refined_list, `[[`, "null_qc"), fill = TRUE)
  if (nrow(refined_results) > 0) {
    refined_results <- merge(
      refined_results,
      pc_summary[, .(
        pc,
        pc_variance_explained = variance_explained,
        pc_cumulative_variance = cumulative_variance_through_pc,
        pc_pearson_age_cor = pc_age_pearson,
        pc_spearman_age_cor = pc_age_spearman,
        clock_direction = clock_coefficient_direction,
        calibrated_pc_coefficient
      )],
      by = "pc",
      all.x = TRUE
    )
    refined_results[, direction := fifelse(NES > 0, "older_aligned", "younger_aligned")]
    all_results <- all_results[!paste(pc, gene_set_id) %in% paste(refined_results$pc, refined_results$gene_set_id)]
    all_results <- rbindlist(list(all_results, refined_results), fill = TRUE)
    null_qc <- rbindlist(list(null_qc, refined_null_qc), fill = TRUE)
  }
  all_results[, abs_NES_rank := NULL]
}

message("Running integrated PCA-clock GSEA with ", n_perm, " permutations.")
integrated_first_pass <- run_gsea_for_pc(
  integrated_rank_vectors[, .(pc, gene_id, rank_score)],
  gene_sets_primary,
  n_perm,
  "first_pass_integrated",
  8000000
)
integrated_results <- integrated_first_pass$results
integrated_null_qc <- integrated_first_pass$null_qc
if (nrow(integrated_results) == 0) stop("Integrated PCA-clock GSEA returned zero tested gene sets.")
integrated_results[, direction := fifelse(NES > 0, "older_aligned", "younger_aligned")]

if (run_refined) {
  message("Running refined integrated PCA-clock GSEA for candidate terms with ", n_perm_refined, " permutations.")
  integrated_results[, abs_NES_rank := frank(-abs(NES), ties.method = "first")]
  integrated_candidates <- integrated_results[empirical_p <= 0.05 | abs_NES_rank <= 100, unique(gene_set_id)]
  if (length(integrated_candidates) > 0) {
    integrated_refined <- run_gsea_for_pc(
      integrated_rank_vectors[, .(pc, gene_id, rank_score)],
      gene_sets_primary[gene_set_id %in% integrated_candidates],
      n_perm_refined,
      "refined_integrated",
      8800000
    )
    integrated_refined_results <- integrated_refined$results
    integrated_refined_null_qc <- integrated_refined$null_qc
    if (nrow(integrated_refined_results) > 0) {
      integrated_refined_results[, direction := fifelse(NES > 0, "older_aligned", "younger_aligned")]
      integrated_results <- integrated_results[!gene_set_id %in% integrated_refined_results$gene_set_id]
      integrated_results <- rbindlist(list(integrated_results, integrated_refined_results), fill = TRUE)
      integrated_null_qc <- rbindlist(list(integrated_null_qc, integrated_refined_null_qc), fill = TRUE)
    }
  }
  integrated_results[, abs_NES_rank := NULL]
}
integrated_results[, fdr_within_pc := p.adjust(empirical_p, method = "BH")]
integrated_results[, fdr_within_pc_collection := p.adjust(empirical_p, method = "BH"), by = collection_label]
setcolorder(integrated_results, c(
  "pc",
  "collection",
  "subcollection",
  "gene_set_name",
  "gene_set_size_original",
  "gene_set_size_in_universe",
  "ES",
  "NES",
  "direction",
  "empirical_p",
  "fdr_within_pc",
  "fdr_within_pc_collection",
  "n_leading_edge",
  "leading_edge_genes"
))
integrated_results[, abs_NES_for_order := abs(NES)]
setorder(integrated_results, fdr_within_pc, empirical_p, -abs_NES_for_order, gene_set_name)
integrated_results[, abs_NES_for_order := NULL]
integrated_significant_terms <- integrated_results[fdr_within_pc < 0.10 & empirical_p < 0.05]
integrated_older_terms <- integrated_results[direction == "older_aligned"]
integrated_younger_terms <- integrated_results[direction == "younger_aligned"]
integrated_leading_edge_genes <- integrated_results[, .(
  pc, collection, subcollection, gene_set_name, direction, NES, empirical_p,
  fdr_within_pc, n_leading_edge, leading_edge_genes
)]
integrated_leading_edge_long <- integrated_results[nzchar(leading_edge_genes), .(
  gene_id = unlist(strsplit(leading_edge_genes, ";", fixed = TRUE))
), by = .(pc, collection, subcollection, gene_set_name, direction, NES, empirical_p, fdr_within_pc)]

all_results[, fdr_within_pc := p.adjust(empirical_p, method = "BH"), by = pc]
all_results[, fdr_within_pc_collection := p.adjust(empirical_p, method = "BH"), by = .(pc, collection_label)]
setcolorder(all_results, c(
  "pc",
  "pc_variance_explained",
  "pc_cumulative_variance",
  "pc_pearson_age_cor",
  "pc_spearman_age_cor",
  "clock_direction",
  "calibrated_pc_coefficient",
  "collection",
  "subcollection",
  "gene_set_name",
  "gene_set_size_original",
  "gene_set_size_in_universe",
  "ES",
  "NES",
  "direction",
  "empirical_p",
  "fdr_within_pc",
  "fdr_within_pc_collection",
  "n_leading_edge",
  "leading_edge_genes"
))
all_results[, abs_NES_for_order := abs(NES)]
setorder(all_results, pc, fdr_within_pc, empirical_p, -abs_NES_for_order, gene_set_name)
all_results[, abs_NES_for_order := NULL]

significant_terms <- all_results[fdr_within_pc < 0.10 & empirical_p < 0.05]
older_terms <- all_results[direction == "older_aligned"]
younger_terms <- all_results[direction == "younger_aligned"]

leading_edge_genes <- all_results[, .(
  pc, collection, subcollection, gene_set_name, direction, NES, empirical_p,
  fdr_within_pc, n_leading_edge, leading_edge_genes
)]
leading_edge_long <- all_results[nzchar(leading_edge_genes), .(
  gene_id = unlist(strsplit(leading_edge_genes, ";", fixed = TRUE))
), by = .(pc, collection, subcollection, gene_set_name, direction, NES, empirical_p, fdr_within_pc)]

collection_summary <- gene_sets[, .(
  n_gene_sets_available = as.integer(.N),
  n_tested_primary_max300 = as.integer(sum(primary_eligible)),
  n_tested_sensitivity_max500 = as.integer(sum(sensitivity_eligible_max500)),
  n_dropped_below_min = as.integer(sum(gene_set_size_in_universe < min_size)),
  n_dropped_above_max300 = as.integer(sum(gene_set_size_in_universe > max_size)),
  n_dropped_above_max500 = as.integer(sum(gene_set_size_in_universe > max_size_sensitivity)),
  median_size_in_universe = as.numeric(median(gene_set_size_in_universe)),
  max_size_in_universe = as.numeric(max(gene_set_size_in_universe))
), by = .(collection, subcollection, collection_label)]
enrich_summary <- all_results[, .(
  n_results_per_pc = .N,
  n_nominal_p_lt_0_05 = sum(empirical_p < 0.05),
  n_enriched_q_lt_0_10_nominal_p_lt_0_05 = sum(fdr_within_pc < 0.10 & empirical_p < 0.05),
  n_older_q_lt_0_10_nominal_p_lt_0_05 = sum(direction == "older_aligned" & fdr_within_pc < 0.10 & empirical_p < 0.05),
  n_younger_q_lt_0_10_nominal_p_lt_0_05 = sum(direction == "younger_aligned" & fdr_within_pc < 0.10 & empirical_p < 0.05)
), by = .(pc, collection, subcollection, collection_label)]
collection_summary <- merge(collection_summary, enrich_summary, by = c("collection", "subcollection", "collection_label"), all.x = TRUE, allow.cartesian = TRUE)
setorder(collection_summary, pc, collection, subcollection)

config <- data.table(
  parameter = c(
    "analysis",
    "model_identity",
    "benchmark_model_reminder",
    "ranking_formula",
    "integrated_ranking_formula",
    "gene_universe",
    "mouse_msigdb_source",
    "available_collection_labels",
    "collection_selection",
    "run_all_collections",
    "min_size_in_universe",
    "max_size_in_universe_primary",
    "max_size_in_universe_sensitivity",
    "n_perm_first_pass",
    "run_refined_permutations",
    "n_perm_refined",
    "seed",
    "permutation_warning",
    "interpretation_caution"
  ),
  value = c(
    "final_pca_pc_preranked_gsea",
    paste0("Final deployment PCA/glmnet clock; training_scope=integrated_bulk_pseudobulk; selector=", pca_selector, "; deployment_all_control_like; selected PCs=", paste(expected_pcs, collapse = ","), "."),
    "The held-out benchmark model is not interpreted by this analysis.",
    "rank_score_gene_pc = pc_loading * calibrated_pc_coefficient; positive scores push predicted age older through that PC.",
    "integrated_rank_score_gene = sum over selected PCs of pc_loading * calibrated_pc_coefficient; positive scores push the integrated PCA clock prediction older.",
    paste0(length(gene_universe), " final selected PCA genes from the deployment artifact."),
    paste(unique(msig_selected$gene_set_source), collapse = ";"),
    paste(available_collections, collapse = ";"),
    collection_selection_note,
    as.character(run_all_collections),
    as.character(min_size),
    as.character(max_size),
    as.character(max_size_sensitivity),
    as.character(n_perm),
    as.character(run_refined),
    as.character(if (run_refined) n_perm_refined else NA_integer_),
    as.character(seed),
    "With 500 first-pass permutations, the minimum empirical P value is about 0.002; FDR values are conservative/coarse and should not be overstated.",
    "PCA loadings are dense; pathway results describe model-oriented enrichment within the selected 2000-gene clock universe, not direct causal gene effects."
  )
)

fwrite(rank_vectors, file.path(out_dir, "final_pca_pc_gsea_rank_vectors.tsv"), sep = "\t")
setorder(integrated_rank_vectors, rank_order)
fwrite(integrated_rank_vectors, file.path(out_dir, "final_pca_integrated_gsea_rank_vector.tsv"), sep = "\t")
fwrite(all_results, file.path(out_dir, "final_pca_pc_gsea_all_results.tsv"), sep = "\t")
fwrite(significant_terms, file.path(out_dir, "final_pca_pc_gsea_significant_terms.tsv"), sep = "\t")
fwrite(older_terms, file.path(out_dir, "final_pca_pc_gsea_older_terms.tsv"), sep = "\t")
fwrite(younger_terms, file.path(out_dir, "final_pca_pc_gsea_younger_terms.tsv"), sep = "\t")
fwrite(leading_edge_genes, file.path(out_dir, "final_pca_pc_gsea_leading_edge_genes.tsv"), sep = "\t")
fwrite(leading_edge_long, file.path(out_dir, "final_pca_pc_gsea_leading_edge_genes_long.tsv"), sep = "\t")
fwrite(integrated_results, file.path(out_dir, "final_pca_integrated_gsea_all_results.tsv"), sep = "\t")
fwrite(integrated_significant_terms, file.path(out_dir, "final_pca_integrated_gsea_significant_terms.tsv"), sep = "\t")
fwrite(integrated_older_terms, file.path(out_dir, "final_pca_integrated_gsea_older_terms.tsv"), sep = "\t")
fwrite(integrated_younger_terms, file.path(out_dir, "final_pca_integrated_gsea_younger_terms.tsv"), sep = "\t")
fwrite(integrated_leading_edge_genes, file.path(out_dir, "final_pca_integrated_gsea_leading_edge_genes.tsv"), sep = "\t")
fwrite(integrated_leading_edge_long, file.path(out_dir, "final_pca_integrated_gsea_leading_edge_genes_long.tsv"), sep = "\t")
fwrite(collection_summary, file.path(out_dir, "final_pca_pc_gsea_collection_summary.tsv"), sep = "\t")
fwrite(rbindlist(list(null_qc, integrated_null_qc), fill = TRUE), file.path(out_dir, "final_pca_pc_gsea_null_qc.tsv"), sep = "\t")
fwrite(integrated_null_qc, file.path(out_dir, "final_pca_integrated_gsea_null_qc.tsv"), sep = "\t")
fwrite(config, file.path(out_dir, "final_pca_pc_gsea_config.tsv"), sep = "\t")

theme_gsea <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey35"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

pretty_term <- function(x, width = 46) {
  x <- gsub("^(HALLMARK|GOBP|GOCC|GOMF|REACTOME|KEGG|WP|BIOCARTA|PID)_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  fixes <- c(
    "\\bDna\\b" = "DNA", "\\bRna\\b" = "RNA", "\\bMrna\\b" = "mRNA",
    "\\bTnf\\b" = "TNF", "\\bIfn\\b" = "IFN", "\\bIl\\b" = "IL",
    "\\bP53\\b" = "p53", "\\bMtor\\b" = "mTOR", "\\bNfkb\\b" = "NFKB",
    "\\bTca\\b" = "TCA", "\\bAtp\\b" = "ATP"
  )
  for (pattern in names(fixes)) x <- gsub(pattern, fixes[[pattern]], x)
  vapply(strwrap(x, width = width, simplify = FALSE), paste, character(1), collapse = "\n")
}

save_plot <- function(plot, stem, width, height) {
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf)
  ggsave(png_path, plot, width = width, height = height, dpi = 300, bg = "white")
  c(pdf_path, png_path)
}

empty_plot <- function(title, subtitle = "No terms passed the plotting criteria.") {
  ggplot() +
    annotate("text", x = 0, y = 0, label = subtitle, size = 5) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    labs(title = title, x = NULL, y = NULL) +
    theme_void(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
}

plot_files <- character()

pc_plot_dt <- melt(
  pc_summary[, .(
    pc,
    `Variance explained` = 100 * variance_explained,
    `Pearson age correlation` = pc_age_pearson,
    `Calibrated PC coefficient` = calibrated_pc_coefficient
  )],
  id.vars = "pc",
  variable.name = "metric",
  value.name = "value"
)
p <- ggplot(pc_plot_dt, aes(x = factor(pc, levels = expected_pcs), y = value, fill = metric)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
  geom_col(width = 0.68, color = "grey25", show.legend = FALSE) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("#4E79A7", "#59A14F", "#E15759")) +
  labs(
    title = "Final deployment PCA/glmnet selected PC summary",
    subtitle = paste0("Selected deployment PCs only: ", paste(expected_pcs, collapse = ", "), "; coefficients are calibrated deployment coefficients."),
    x = "Selected PC",
    y = NULL
  ) +
  theme_gsea()
plot_files <- c(plot_files, save_plot(p, "01_pc_summary_variance_age_coefficients", 13, 4.5))

count_dt <- all_results[fdr_within_pc < 0.10 & empirical_p < 0.05, .N, by = .(pc, direction)]
if (nrow(count_dt) == 0) {
  count_dt <- CJ(pc = expected_pcs, direction = c("older_aligned", "younger_aligned"))
  count_dt[, N := 0L]
}
p <- ggplot(count_dt, aes(x = factor(pc, levels = expected_pcs), y = N, fill = direction)) +
  geom_col(position = "dodge", color = "grey25") +
  scale_fill_manual(values = c(older_aligned = "#C65A4A", younger_aligned = "#3E7CB1"), labels = c("Older aligned", "Younger aligned")) +
  labs(
    title = "Enriched term counts by PC and direction",
    subtitle = "Counted at within-PC q < 0.10 and nominal empirical P < 0.05.",
    x = "PC",
    y = "Number of enriched terms",
    fill = NULL
  ) +
  theme_gsea()
plot_files <- c(plot_files, save_plot(p, "02_gsea_term_counts_by_pc_direction", 9, 5.5))

select_top_by_direction <- function(dt, direction_value, n = top_n_plot, pcs = expected_pcs) {
  out <- dt[pc %in% pcs & direction == direction_value]
  out[, sort_key := fdr_within_pc + empirical_p / 1000]
  out <- out[order(pc, sort_key, -abs(NES), gene_set_name)]
  out[, head(.SD, n), by = pc]
}

focused_heatmap_gene_sets <- c(
  "REACTOME_METABOLISM_OF_STEROIDS",
  "GOBP_STEROL_METABOLIC_PROCESS",
  "GOBP_STEROID_METABOLIC_PROCESS",
  "GOBP_STEROID_BIOSYNTHETIC_PROCESS",
  "GOBP_SECONDARY_ALCOHOL_METABOLIC_PROCESS",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "GOBP_FATTY_ACID_METABOLIC_PROCESS",
  "REACTOME_FATTY_ACID_METABOLISM",
  "REACTOME_METABOLISM_OF_LIPIDS",
  "REACTOME_METABOLISM",
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
  "WP_ELECTRON_TRANSPORT_CHAIN",
  "WP_TCA_CYCLE",
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE"
)

term_biology_group <- function(x) {
  out <- rep("Other", length(x))
  out[grepl("STEROID|STEROL|SECONDARY_ALCOHOL|CHOLESTEROL|BILE_ACID", x)] <- "Steroid / sterol"
  out[grepl("OXIDATIVE_PHOSPHORYLATION|MITOCHONDR|RESPIRATORY|ELECTRON_TRANSPORT|TCA|CITRIC_ACID|PYRUVATE", x)] <- "Mitochondrial / energy"
  out[grepl("FATTY_ACID|LIPID|ADIPOGENESIS|BETAOXIDATION|PHOSPHOLIPID|SPHINGOLIPID|GLYCEROLIPID", x)] <- "Fatty acid / lipid"
  out[grepl("METABOL", x) & out == "Other"] <- "Other metabolism"
  out[grepl("IMMUNE|LEUKOCYTE|LYMPHOCYTE|T_CELL|B_CELL|NK_CELL|NEUTROPHIL|MACROPHAGE|MONOCYTE|CYTOKINE|INTERFERON|IL[0-9]|TNF|INFLAM", x)] <- "Immune / inflammatory"
  out[grepl("CELL_CYCLE|MITOTIC|CHROMOSOME|CHROMATID|DNA_|DNA$|E2F|G2M|CHECKPOINT|REPLICATION", x)] <- "Cell cycle / DNA"
  out
}

make_top_term_plot <- function(dt, title, subtitle) {
  if (nrow(dt) == 0) return(empty_plot(title))
  dt <- copy(dt)
  dt[, significance_class := fifelse(
    fdr_within_pc < 0.10 & empirical_p < 0.05 & NES > 0,
    "Positive NES; q<0.10 & P<0.05",
    fifelse(
      fdr_within_pc < 0.10 & empirical_p < 0.05 & NES < 0,
      "Negative NES; q<0.10 & P<0.05",
      "Top ranked"
    )
  )]
  dt[, term_label := pretty_term(gene_set_name, width = 40)]
  dt[, term_plot_chr := make.unique(paste(pc, term_label, sep = " | "), sep = " #")]
  dt[, term_plot := factor(term_plot_chr, levels = rev(term_plot_chr))]
  ggplot(dt, aes(x = NES, y = term_plot, fill = significance_class)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    geom_col(color = "grey25", width = 0.72) +
    facet_wrap(~pc, scales = "free_y", ncol = 2) +
    scale_y_discrete(labels = function(x) sub(" #[0-9]+$", "", sub("^.* \\| ", "", x))) +
    scale_fill_manual(values = c(
      "Positive NES; q<0.10 & P<0.05" = "#8B1E1E",
      "Negative NES; q<0.10 & P<0.05" = "#3E7CB1",
      "Top ranked" = "grey70"
    )) +
    labs(title = title, subtitle = subtitle, x = "NES", y = NULL, fill = NULL) +
    theme_gsea() +
    theme(axis.text.y = element_text(size = 9), axis.text.x = element_text(angle = 0))
}

top_older <- select_top_by_direction(all_results, "older_aligned")
p <- make_top_term_plot(top_older, "Top older-aligned GSEA terms by PC", "Positive NES means the older-predictive side of the PCA clock.")
plot_files <- c(plot_files, save_plot(p, "03_top_older_terms_by_pc", 13, 11))

top_younger <- select_top_by_direction(all_results, "younger_aligned")
p <- make_top_term_plot(top_younger, "Top younger-aligned GSEA terms by PC", "Negative NES means the younger-predictive side of the PCA clock.")
plot_files <- c(plot_files, save_plot(p, "04_top_younger_terms_by_pc", 13, 11))

integrated_top_terms <- rbindlist(list(
  select_top_by_direction(integrated_results, "older_aligned", n = 15, pcs = "integrated_PCA_clock"),
  select_top_by_direction(integrated_results, "younger_aligned", n = 15, pcs = "integrated_PCA_clock")
), fill = TRUE)
p <- make_top_term_plot(
  integrated_top_terms,
  "Integrated PCA-clock GSEA terms",
  "Genes are ranked once by summed loading x calibrated coefficient across PC1-PC7.\nNES sign is the integrated PCA-clock direction."
)
plot_files <- c(plot_files, save_plot(p, "12_integrated_pca_clock_gsea_terms", 11, 9))

heat_terms <- rbindlist(list(top_older[, .(gene_set_id)], top_younger[, .(gene_set_id)]), fill = TRUE)
heat_terms <- unique(heat_terms$gene_set_id)
focus_heat_terms <- all_results[
  pc %in% c("PC3", "PC4") & gene_set_name %in% focused_heatmap_gene_sets,
  unique(gene_set_id)
]
heat_terms <- unique(c(heat_terms, focus_heat_terms))
max_heat_terms <- 45L
if (length(heat_terms) > max_heat_terms) {
  term_rank <- all_results[gene_set_id %in% heat_terms, .(
    best_fdr = min(fdr_within_pc),
    best_abs_nes = max(abs(NES)),
    focus_term = gene_set_id[1] %in% focus_heat_terms
  ), by = gene_set_id]
  focus_rank <- term_rank[focus_term == TRUE][order(best_fdr, -best_abs_nes, gene_set_id)]
  fill_rank <- term_rank[focus_term == FALSE][order(best_fdr, -best_abs_nes, gene_set_id)]
  focus_keep <- head(focus_rank$gene_set_id, max_heat_terms)
  fill_keep <- head(fill_rank$gene_set_id, max_heat_terms - length(focus_keep))
  heat_terms <- c(focus_keep, fill_keep)
}
heat_dt <- all_results[gene_set_id %in% heat_terms]
if (nrow(heat_dt) > 0) {
  heat_dt[, term_label := pretty_term(gene_set_name, width = 42)]
  heat_dt[, focus_term := gene_set_id %in% focus_heat_terms]
  biology_group_levels <- c(
    "Steroid / sterol",
    "Mitochondrial / energy",
    "Fatty acid / lipid",
    "Other metabolism",
    "Immune / inflammatory",
    "Cell cycle / DNA",
    "Other"
  )
  heat_dt[, biology_group := factor(term_biology_group(gene_set_name), levels = biology_group_levels)]
  heat_label_map <- heat_dt[, .(
    term_label = term_label[1],
    biology_group = biology_group[1],
    focus_term = any(focus_term),
    best_fdr = min(fdr_within_pc),
    best_abs_nes = max(abs(NES))
  ), by = gene_set_id]
  setorder(heat_label_map, biology_group, -focus_term, best_fdr, -best_abs_nes, term_label)
  heat_label_map[, term_heat_label := make.unique(term_label, sep = " #")]
  heat_dt <- merge(heat_dt, heat_label_map[, .(gene_set_id, term_heat_label)], by = "gene_set_id", all.x = TRUE, sort = FALSE)
  heat_dt[, term_heat_label := factor(term_heat_label, levels = rev(heat_label_map$term_heat_label))]
  p <- ggplot(heat_dt, aes(x = factor(pc, levels = expected_pcs), y = term_heat_label, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(biology_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#3E7CB1", mid = "white", high = "#C65A4A", midpoint = 0) +
    labs(
      title = "NES heatmap for top deployment PCA GSEA terms",
      subtitle = "Rows are top terms per PC plus focused steroid and metabolic terms.\nFacets group related biology together.",
      x = "PC",
      y = NULL,
      fill = "NES"
    ) +
    theme_gsea() +
    theme(
      axis.text.x = element_text(angle = 0),
      axis.text.y = element_text(size = 8),
      panel.spacing.y = unit(0.18, "lines"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9)
    ) +
    scale_y_discrete(labels = function(x) sub(" #[0-9]+$", "", x))
} else {
  p <- empty_plot("NES heatmap for top deployment PCA GSEA terms")
}
plot_files <- c(plot_files, save_plot(p, "05_gsea_nes_heatmap_top_terms", 11, 13))

if (nrow(heat_dt) > 0) {
  integrated_heat_dt <- integrated_results[gene_set_id %in% heat_terms]
  integrated_heat_dt[, `:=`(
    term_label = pretty_term(gene_set_name, width = 42),
    focus_term = gene_set_id %in% focus_heat_terms,
    biology_group = factor(term_biology_group(gene_set_name), levels = biology_group_levels)
  )]
  integrated_heat_dt <- merge(
    integrated_heat_dt,
    heat_label_map[, .(gene_set_id, term_heat_label)],
    by = "gene_set_id",
    all.x = TRUE,
    sort = FALSE
  )
  integrated_heat_dt <- integrated_heat_dt[!is.na(term_heat_label)]
  heat_dt_plus_integrated <- rbindlist(list(
    heat_dt[, .(
      source = "pc_specific",
      pc,
      pc_plot = pc,
      gene_set_id,
      collection,
      subcollection,
      collection_label,
      gene_set_name,
      gene_set_size_in_universe,
      NES,
      direction,
      empirical_p,
      fdr_within_pc,
      term_label,
      term_heat_label,
      biology_group
    )],
    integrated_heat_dt[, .(
      source = "integrated_total",
      pc,
      pc_plot = "Integrated",
      gene_set_id,
      collection,
      subcollection,
      collection_label,
      gene_set_name,
      gene_set_size_in_universe,
      NES,
      direction,
      empirical_p,
      fdr_within_pc,
      term_label,
      term_heat_label,
      biology_group
    )]
  ), fill = TRUE)
  heat_dt_plus_integrated[, pc_plot := factor(pc_plot, levels = c(expected_pcs, "Integrated"))]
  heat_dt_plus_integrated[, term_heat_label := factor(term_heat_label, levels = rev(heat_label_map$term_heat_label))]
  fwrite(
    heat_dt_plus_integrated[order(biology_group, term_heat_label, pc_plot)],
    file.path(out_dir, "final_pca_pc_plus_integrated_gsea_heatmap_source.tsv"),
    sep = "\t"
  )
  p <- ggplot(heat_dt_plus_integrated, aes(x = pc_plot, y = term_heat_label, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_vline(xintercept = length(expected_pcs) + 0.5, color = "grey25", linewidth = 0.45) +
    facet_grid(biology_group ~ ., scales = "free_y", space = "free_y", switch = "y") +
    scale_fill_gradient2(low = "#3E7CB1", mid = "white", high = "#C65A4A", midpoint = 0) +
    labs(
      title = "PC-specific and integrated PCA-clock GSEA NES",
      subtitle = "Rows match the top-term heatmap.\nThe rightmost column ranks each gene once by summed PC1-PC7 contribution.",
      x = "PCA clock component",
      y = NULL,
      fill = "NES"
    ) +
    theme_gsea() +
    theme(
      axis.text.x = element_text(angle = 0),
      axis.text.y = element_text(size = 8),
      panel.spacing.y = unit(0.18, "lines"),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9)
    ) +
    scale_y_discrete(labels = function(x) sub(" #[0-9]+$", "", x))
} else {
  p <- empty_plot("PC-specific and integrated PCA-clock GSEA NES")
}
plot_files <- c(plot_files, save_plot(p, "13_gsea_nes_heatmap_top_terms_with_integrated_total", 13, 13))

running_curve_dt <- function(pc_name, term_rows) {
  if (nrow(term_rows) == 0) return(data.table())
  pc_rank <- rank_vectors[pc == pc_name][order(-rank_score)]
  stats <- pc_rank$rank_score
  names(stats) <- pc_rank$gene_id
  rbindlist(lapply(seq_len(nrow(term_rows)), function(i) {
    genes <- gene_sets_primary[gene_set_id == term_rows$gene_set_id[[i]], genes][[1]]
    genes <- intersect(genes, names(stats))
    pos <- sort(match(genes, names(stats)))
    nh <- length(pos)
    n <- length(stats)
    if (nh == 0 || nh >= n) return(NULL)
    weights <- abs(stats[pos])^score_type
    if (!is.finite(sum(weights)) || sum(weights) == 0) weights <- rep(1, length(pos))
    increments <- rep(-1 / (n - nh), n)
    increments[pos] <- weights / sum(weights)
    data.table(
      rank = seq_len(n),
      running_ES = cumsum(increments),
      hit = seq_len(n) %in% pos,
      gene_set_name = term_rows$gene_set_name[[i]],
      term_label = pretty_term(term_rows$gene_set_name[[i]], width = 36),
      direction = term_rows$direction[[i]]
    )
  }), fill = TRUE)
}

make_running_plot <- function(pc_name, stem_title) {
  terms <- rbindlist(list(
    select_top_by_direction(all_results[pc == pc_name], "older_aligned", n = 3, pcs = pc_name),
    select_top_by_direction(all_results[pc == pc_name], "younger_aligned", n = 3, pcs = pc_name)
  ), fill = TRUE)
  curve_dt <- running_curve_dt(pc_name, terms)
  if (nrow(curve_dt) == 0) return(empty_plot(stem_title))
  ggplot(curve_dt, aes(x = rank, y = running_ES, color = term_label)) +
    geom_hline(yintercept = 0, color = "grey55") +
    geom_line(linewidth = 0.7) +
    facet_wrap(~direction, scales = "free_y", ncol = 1) +
    labs(
      title = stem_title,
      subtitle = "Top older- and younger-aligned terms; genes ranked by loading x calibrated PC coefficient.",
      x = "Rank in preranked 2000-gene universe",
      y = "Running enrichment score",
      color = NULL
    ) +
    theme_gsea() +
    theme(axis.text.x = element_text(angle = 0))
}

p <- make_running_plot("PC1", "Running enrichment curves for top PC1 terms")
plot_files <- c(plot_files, save_plot(p, "06_running_sum_top_pc1_terms", 10, 8))

p <- make_running_plot("PC2", "Running enrichment curves for top PC2 terms")
plot_files <- c(plot_files, save_plot(p, "07_running_sum_top_pc2_terms", 10, 8))

primary <- rbindlist(list(
  select_top_by_direction(all_results, "older_aligned", n = 10, pcs = c("PC1", "PC2", "PC3", "PC4")),
  select_top_by_direction(all_results, "younger_aligned", n = 10, pcs = c("PC1", "PC2", "PC3", "PC4"))
), fill = TRUE)
p <- make_top_term_plot(primary, "Primary PC GSEA terms", "PC1-PC4 are shown as the primary display block: PC1/PC2 are the strongest age axes, while PC3/PC4 add interpretable steroid and metabolic biology.")
plot_files <- c(plot_files, save_plot(p, "11_primary_pc_terms", 13, 14))

secondary <- rbindlist(list(
  select_top_by_direction(all_results, "older_aligned", n = 5, pcs = c("PC5", "PC6", "PC7")),
  select_top_by_direction(all_results, "younger_aligned", n = 5, pcs = c("PC5", "PC6", "PC7"))
), fill = TRUE)
p <- make_top_term_plot(secondary, "Secondary PC GSEA terms", "PC5, PC6, and PC7 are selected deployment modifiers with lower variance and narrower residual biology despite some large coefficients.")
plot_files <- c(plot_files, save_plot(p, "08_secondary_pc_terms", 12, 9))

overlap_base <- significant_terms
if (nrow(overlap_base) == 0) {
  overlap_base <- rbindlist(list(
    select_top_by_direction(all_results, "older_aligned", n = 5),
    select_top_by_direction(all_results, "younger_aligned", n = 5)
  ), fill = TRUE)
}
overlap_sets <- overlap_base[nzchar(leading_edge_genes), .(
  genes = list(unique(unlist(strsplit(leading_edge_genes, ";", fixed = TRUE))))
), by = .(pc, direction)]
if (nrow(overlap_sets) >= 2) {
  overlap_sets[, label := paste(pc, sub("_aligned", "", direction), sep = " ")]
  pairs <- CJ(i = seq_len(nrow(overlap_sets)), j = seq_len(nrow(overlap_sets)))
  overlap_dt <- pairs[, {
    a <- overlap_sets$genes[[i]]
    b <- overlap_sets$genes[[j]]
    data.table(
      label_x = overlap_sets$label[[j]],
      label_y = overlap_sets$label[[i]],
      jaccard = length(intersect(a, b)) / length(union(a, b)),
      overlap_n = length(intersect(a, b))
    )
  }, by = .(i, j)]
  p <- ggplot(overlap_dt, aes(x = label_x, y = label_y, fill = jaccard)) +
    geom_tile(color = "white") +
    geom_text(aes(label = overlap_n), size = 3) +
    scale_fill_gradient(low = "white", high = "#4E79A7") +
    labs(
      title = "Leading-edge gene overlap",
      subtitle = "Numbers show shared leading-edge genes; color shows Jaccard overlap.",
      x = NULL,
      y = NULL,
      fill = "Jaccard"
    ) +
    theme_gsea() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(angle = 0))
} else {
  p <- empty_plot("Leading-edge gene overlap")
}
plot_files <- c(plot_files, save_plot(p, "09_leading_edge_gene_overlap", 9, 8))

coll_plot <- collection_summary[, .(
  n_tested_primary_max300 = max(n_tested_primary_max300, na.rm = TRUE),
  n_enriched = sum(n_enriched_q_lt_0_10_nominal_p_lt_0_05, na.rm = TRUE)
), by = collection_label]
coll_plot <- melt(coll_plot, id.vars = "collection_label", variable.name = "metric", value.name = "count")
p <- ggplot(coll_plot, aes(x = reorder(collection_label, count), y = count, fill = metric)) +
  geom_col(position = "dodge", color = "grey25") +
  coord_flip() +
  scale_fill_manual(values = c(n_tested_primary_max300 = "#4E79A7", n_enriched = "#C65A4A"), labels = c("Tested", "Enriched PC-collection hits")) +
  labs(
    title = "MSigDB collection summary",
    subtitle = "Tested terms use the primary 10-300 genes-in-universe filter.",
    x = NULL,
    y = "Term count",
    fill = NULL
  ) +
  theme_gsea() +
  theme(axis.text.x = element_text(angle = 0))
plot_files <- c(plot_files, save_plot(p, "10_collection_summary", 9, 6))

format_top_terms_md <- function(dt, direction_value, n = 5) {
  rows <- dt[direction == direction_value][order(fdr_within_pc, empirical_p, -abs(NES))]
  rows <- rows[, head(.SD, n), by = pc]
  if (nrow(rows) == 0) return("- No terms met plotting criteria.")
  rows[, sprintf(
    "- %s: %s (NES %.2f, P %.4f, q %.4f)",
    pc, gene_set_name, NES, empirical_p, fdr_within_pc
  )]
}

leading_highlights <- leading_edge_long[, .N, by = gene_id][order(-N, gene_id)][1:min(.N, 20)]
leading_text <- if (nrow(leading_highlights) > 0) {
  paste(sprintf("- %s (%s leading-edge appearances)", leading_highlights$gene_id, leading_highlights$N), collapse = "\n")
} else {
  "- No leading-edge genes were available."
}

all_q_nonsig <- all_results[fdr_within_pc < 0.10 & empirical_p < 0.05, .N] == 0
summary_lines <- c(
  "# Final PCA/glmnet deployment PC-level GSEA",
  "",
  "## Model Identity",
  "",
  paste0(
    "This analysis interprets the final deployment PCA/glmnet ovarian clock: `training_scope = integrated_bulk_pseudobulk`, selector `",
    pca_selector,
    "`, refit on all ",
    length(artifact$training_sample_keys %||% selected_genes),
    " control-like samples, selected PCs `",
    paste(expected_pcs, collapse = ", "),
    "`, and ",
    length(selected_genes),
    " selected PCA genes."
  ),
  "",
  "The held-out benchmark PCA model is not interpreted here.",
  "",
  "## Ranking",
  "",
  "`rank_score = pc_loading x calibrated_pc_coefficient`.",
  "",
  "Positive NES means enrichment on the older-predictive side of the final PCA clock. Negative NES means enrichment on the younger-predictive side. Raw PC loading sign alone was not used as the biological direction.",
  "",
  "An additional integrated PCA-clock GSEA ranks each of the same 2000 PCA input genes once by the summed contribution `sum(pc_loading x calibrated_pc_coefficient)` across all selected PCs. This integrated analysis is better suited for model-level pathway direction than the PC-specific heatmap.",
  "",
  "## Gene Sets",
  "",
  paste0("Gene universe: ", length(gene_universe), " final selected PCA genes."),
  paste0("MSigDB source: ", paste(unique(msig_selected$gene_set_source), collapse = "; "), "."),
  paste0("MSigDB version(s): ", paste(sort(unique(msig_selected$msigdb_version)), collapse = "; "), "."),
  paste0("Collections used: ", paste(sort(unique(msig_selected$collection_label)), collapse = "; "), "."),
  paste0("Primary size filter: ", min_size, "-", max_size, " genes within the 2000-gene universe. Sensitivity counts were also written for max size ", max_size_sensitivity, "."),
  "",
  "## Permutations",
  "",
  paste0("First-pass weighted-KS GSEA used exponent p = 1, seed ", seed, ", and ", n_perm, " random gene-set permutations per PC and set size."),
  "With 500 first-pass permutations, the minimum empirical P value is about 0.002, so FDR values are conservative/coarse and should not be overstated.",
  if (run_refined) paste0("Refined permutations were run for candidate terms with ", n_perm_refined, " permutations.") else "Refined permutations were not run; outputs are marked as first-pass.",
  if (all_q_nonsig) "Warning: no terms passed q < 0.10 and nominal P < 0.05; this may reflect limited first-pass permutation resolution." else "",
  "",
  "## Top Older-Aligned Terms",
  "",
  format_top_terms_md(all_results, "older_aligned", 5),
  "",
  "## Top Younger-Aligned Terms",
  "",
  format_top_terms_md(all_results, "younger_aligned", 5),
  "",
  "## Integrated PCA-Clock GSEA",
  "",
  "These terms use one aggregate rank vector across PC1-PC7 rather than separate PC-specific rank vectors.",
  "",
  "### Older-Aligned Integrated Terms",
  "",
  format_top_terms_md(integrated_results, "older_aligned", 10),
  "",
  "### Younger-Aligned Integrated Terms",
  "",
  format_top_terms_md(integrated_results, "younger_aligned", 10),
  "",
  "## Leading-Edge Highlights",
  "",
  leading_text,
  "",
  "## PC-Specific Interpretation Strength",
  "",
  "- PC1/PC2: strongest direct age-associated axes.",
  "- PC3/PC4: primary biological interpretation block with substantial cumulative variance and interpretable steroid, sterol, mitochondrial, and lipid-metabolism pathway structure.",
  "- PC5/PC6/PC7: selected secondary deployment modifiers; PC5 and PC7 have large calibrated coefficients, but each explains less variance than PC1-PC4 and should be interpreted as narrower residual axes.",
  "",
  "## Caution",
  "",
  "PCA loadings are dense and pathway results describe model-oriented gene-set enrichment within the selected 2000-gene clock universe. They are not direct causal gene effects."
)
writeLines(summary_lines[nzchar(summary_lines) | c(TRUE, head(nzchar(summary_lines), -1))], file.path(out_dir, "final_pca_pc_gsea_summary.md"))

table_count <- length(list.files(out_dir, pattern = "\\.tsv$"))
plot_count <- length(unique(sub("\\.(png|pdf)$", "", basename(plot_files))))
message("Wrote final PCA deployment PC GSEA outputs under:")
message(out_dir)
message("Tables written: ", table_count)
message("Plot stems written: ", plot_count)
message("Tested gene sets per PC: ", paste(all_results[, .N, by = pc][, paste(pc, N, sep = "=")], collapse = ", "))
message("Done.")
