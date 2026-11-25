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

# Figure contract: show signed pathway enrichment for active frozen PCs and the
# integrated standardized gene coefficient. Retain the 13 reference analysis
# views, use current-model labels, and export PDF/PNG plus complete TSV sources.
arg <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
arg <- gsub("~+~", " ", arg, fixed=TRUE)
root <- normalizePath(file.path(dirname(arg), ".."))
development <- normalizePath(file.path(root, "../.."))
setDTthreads(1)
out_dir <- file.path(root, "12_PCA deployment model GSEA")
diag_dir <- file.path(root, "10_pca_diagnostics")
artifact <- readRDS(file.path(development, "locked_models/count_tuning_20261011/VersionB-PCA-glmnet/model.rds"))
stopifnot(artifact$recipe == "VersionB-PCA", artifact$calibration == "none", artifact$target == "age_week")
selected_genes <- artifact$representation$genes
gene_universe <- selected_genes
expected_pcs <- colnames(artifact$representation$rotation)[artifact$pc_coefficients != 0]
pca_selector <- "VersionB_count_tuning_locked_20261011"
pc_summary <- fread(file.path(diag_dir, "final_pca_clock_selected_pc_summary.tsv"))
loadings <- fread(file.path(diag_dir, "final_pca_clock_gene_pc_loadings_with_feature_scores.tsv"))
stopifnot(setequal(pc_summary$pc, expected_pcs), setequal(loadings$pc, expected_pcs),
          setequal(loadings$gene_id, selected_genes), !anyDuplicated(selected_genes))
stopifnot(all(loadings[, uniqueN(gene_id), by=pc]$V1 == length(gene_universe)))
stopifnot(max(abs(pc_summary$pc_coefficient[match(expected_pcs,pc_summary$pc)] - artifact$pc_coefficients[artifact$pc_coefficients!=0])) < 1e-12)
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
stopifnot(n_perm>=10, min_size>=1, max_size>min_size)
clean_gene <- function(x) {x<-trimws(as.character(x));x[x %in% c("","NA","NaN")]<-NA_character_;x}

rank_vectors <- merge(
  loadings[, .(gene_id, pc, pc_loading)],
  pc_summary[, .(
    pc,
    pc_variance_explained = variance_explained,
    pc_cumulative_variance = cumulative_variance_through_pc,
    pc_pearson_age_cor = pc_age_pearson,
    pc_spearman_age_cor = pc_age_spearman,
    clock_direction = clock_coefficient_direction,
    pc_coefficient
  )],
  by = "pc",
  all.x = TRUE
)
rank_vectors[, rank_score := as.numeric(pc_loading) * as.numeric(pc_coefficient)]
rank_vectors[, rank_direction := fifelse(rank_score >= 0, "older_predictive", "younger_predictive")]
rank_vectors[, abs_rank_score := abs(rank_score)]
rank_vectors[, rank_order := frank(-rank_score, ties.method = "first"), by = pc]
rank_vectors[, abs_rank_order := frank(-abs_rank_score, ties.method = "first"), by = pc]
if (rank_vectors[!is.finite(rank_score), .N] > 0) stop("Rank vectors contain non-finite rank scores.")
if (rank_vectors[, uniqueN(gene_id), by = pc][V1 != length(gene_universe), .N] > 0) {
  stop("At least one PC has an incomplete set of ranked selected genes.")
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
if (integrated_rank_vectors[, uniqueN(gene_id)] != length(gene_universe)) stop("Integrated rank vector has an incomplete set of ranked selected genes.")

load_msigdb_mouse <- function() {
  paths <- sort(list.files(file.path(root,"provenance/msigdb"), pattern="^msigdb.*Mm.*rds$",full.names=TRUE))
  stopifnot(length(paths)==7)
  msig <- rbindlist(lapply(paths,readRDS),use.names=TRUE,fill=TRUE)
  stopifnot(all(msig$db_target_species=="MM"),all(msig$db_version=="2026.1.Mm"))
  unique(msig[,.(gene_set_name=gs_name, gene_symbol=clean_gene(db_gene_symbol),
    collection=gs_collection, subcollection=gs_subcollection,
    collection_name=gs_collection_name,gene_set_description=gs_description,
    msigdb_version=db_version,gene_set_source="msigdbr_mouse_native_cached_snapshot")])
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
  stop("No gene sets remain after intersecting with the 3000-gene universe and applying size filters.")
}

fwrite(gene_sets[, .(gene_set_id,collection,subcollection,collection_label,gene_set_name,gene_set_size_original,gene_set_size_in_universe,primary_eligible,sensitivity_eligible_max500,genes=vapply(genes,paste,character(1),collapse=";"))], file.path(out_dir,"gene_set_universe_audit.tsv"),sep="\t")

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
  if (n != length(gene_universe)) stop("Rank vector has an incomplete set of finite genes for ", pc_rank$pc[[1]])

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
    pc_coefficient
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
        pc_coefficient
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
  "pc_coefficient",
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
    paste0("Final deployment PCA/glmnet clock; training_scope=locked_count_tuning_all_136; selector=", pca_selector, "; no calibration; active PCs=", paste(expected_pcs, collapse = ","), "."),
    "The held-out benchmark model is not interpreted by this analysis.",
    "rank_score_gene_pc = pc_loading * pc_coefficient; positive scores push predicted age older through that PC.",
    "integrated_rank_score_gene = sum over selected PCs of pc_loading * pc_coefficient; positive scores push the integrated PCA clock prediction older.",
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
    "With 500 first-pass permutations, the minimum empirical P value is about 0.002; P values and BH q values have coarse Monte Carlo resolution.",
    "PCA loadings are dense; pathway results describe model-oriented enrichment within the selected 3000-gene clock universe, not direct causal gene effects."
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
    `Frozen PC coefficient` = pc_coefficient
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
    subtitle = paste0("Selected deployment PCs only: ", paste(expected_pcs, collapse = ", "), "; coefficients are frozen deployment coefficients (no calibration)."),
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

# Descriptive, name-based display groups; these are not inferred pathway clusters.
# Later, more specific rules take precedence. The complete assignment is exported.
biology_rules <- data.table(
  priority = 1:10,
  biology_group = c("Signalling / development", "Other metabolism", "Protein / RNA processing",
                    "Adhesion / extracellular matrix", "Neural / synaptic", "Fatty acid / lipid",
                    "Mitochondrial / energy", "Steroid / sterol", "Cell cycle / DNA", "Immune / inflammatory"),
  pattern = c(
    "SIGNAL|DEVELOPMENT|DIFFERENTIATION|MORPHOGEN|WNT|NOTCH|HEDGEHOG|ANGIOGEN|VASCU|GTPASE|GPCR|SECOND_MESSENGER|ENDOCRINE|HORMONE",
    "METABOL|AMINO_ACID|GLYCOL|GLUCOSE|GLUCONEOGEN|CARBOHYDRATE|NUCLEOTIDE_BIOSYN|RIBOSE_PHOSPHATE|ACETYL_COA|ACYL_COA|BIOSYNTHETIC_PROCESS",
    "RNA_PROCESS|MRNA|NCRNA|TRNA|RNA_SPLIC|SPLICING|RNA_METABOL|METABOLISM_OF_RNA|TRANSLATION|RIBOSOM|PROTEIN_FOLD|UNFOLDED_PROTEIN|PROTEASOM|PROTEIN_DEGRADATION|AUTOPHAG|LYSOSOM|PROTEIN_STABILITY",
    "EXTRACELLULAR_MATRIX|COLLAGEN|INTEGRIN|CELL_ADHESION|CELL_CELL_ADHESION|CELL_MATRIX_ADHESION|L1CAM|EPITHELIAL_MESENCHYMAL",
    "NEURON|NEURAL|SYNAP|AXON|NEUROTRANSM|NEUROGEN|NERVOUS_SYSTEM|NEUROPEPTIDE",
    "FATTY_ACID|LIPID|ADIPOGENESIS|BETAOXIDATION|ACYL_COA|LIPOPROTEIN",
    "OXIDATIVE_PHOSPHORYLATION|MITOCHONDR|AEROBIC_RESPIRATION|RESPIRATORY_ELECTRON|ELECTRON_TRANSPORT|(^|_)TCA(_|$)|CITRIC_ACID|PYRUVATE|ATP_SYNTHESIS",
    "STEROID|STEROL|SECONDARY_ALCOHOL|CHOLESTEROL|BILE_ACID",
    "CELL_CYCLE|MITOTIC|MITOSIS|CYTOKINESIS|CHROMOSOME|CHROMATID|CHROMATIN|DNA_(REPAIR|REPLICATION|RECOMBINATION|METABOL|SYNTHESIS|DAMAGE|INTEGRITY|TOPOLOG)|DNA$|E2F|G2M|CHECKPOINT|REPLICATION|TELOMERE|HOMOLOGY_DIRECTED_REPAIR|HOMOLOGOUS_RECOMBINATION|RECOMBINATIONAL_REPAIR|STRAND_BREAK_REPAIR|KINETOCHORE|SPINDLE|CELL_DIVISION|NUCLEAR_DIVISION|(^|_)[MS]_PHASE($|_)|G1_S_TRANSITION|G2_M_TRANSITION",
    "^CUI_|IMMUN|LEUKOCYTE|LYMPH|T_CELL|B_CELL|NK_CELL|NATURAL_KILLER|NEUTROPHIL|MACROPHAGE|MONOCYTE|MONONUCLEAR|MYELOID|GRANULOCYTE|CYTOKINES?(_|$)|CHEMOKINE|INTERFERON|INTERLEUKIN|IL[0-9]|TNF|TUMOR_NECROSIS|INFLAM|ANTIGEN|PHAGOCYT|COMPLEMENT|DEFENSE_RESPONSE|RESPIRATORY_BURST|FC_RECEPTOR|FC_GAMMA|MHC_|ALLOGRAFT|HEMOPOIESIS|HEMATOPOIESIS"
  )
)
term_biology_group <- function(x) {
  out <- rep("Other", length(x))
  for (i in seq_len(nrow(biology_rules))) {
    out[grepl(biology_rules$pattern[i], x)] <- biology_rules$biology_group[i]
  }
  out
}
biology_group_levels <- c("Mitochondrial / energy", "Fatty acid / lipid", "Other metabolism",
                          "Protein / RNA processing", "Adhesion / extracellular matrix",
                          "Neural / synaptic", "Signalling / development", "Immune / inflammatory",
                          "Cell cycle / DNA", "Steroid / sterol", "Other")

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
  "Genes ranked by summed loading x frozen coefficient.\nActive PCs: PC1, PC2, PC5 and PC8; NES sign gives integrated clock direction."
)
plot_files <- c(plot_files, save_plot(p, "12_integrated_pca_clock_gsea_terms", 11, 9))

# Figure contract: show the breadth and direction of supported integrated biology
# while guaranteeing visibility of the top five terms in each direction of every PC.
# Keep the numerical GSEA unchanged; export selection reasons, ranks, q and membership.
heat_top_n <- 5L
functional_n <- 4L
ranked_heat_results <- rbindlist(list(copy(all_results), copy(integrated_results)), fill = TRUE)
ranked_heat_results[, pc := as.character(pc)]
ranked_heat_results[, sort_key := fdr_within_pc + empirical_p / 1000]
ranked_heat_results <- ranked_heat_results[order(pc, direction, sort_key, -abs(NES), gene_set_name)]
ranked_heat_results[, directional_rank := seq_len(.N), by = .(pc, direction)]
ranked_heat_results[, biology_group := term_biology_group(gene_set_name)]
# M7 is an immune-signature collection even where its name lacks an immune keyword.
ranked_heat_results[collection == "M7", biology_group := "Immune / inflammatory"]
ranked_heat_results[, `:=`(
  significant_q05 = fdr_within_pc < 0.05 & empirical_p < 0.05,
  significant_q10 = fdr_within_pc < 0.10 & empirical_p < 0.05,
  pc_top_term = pc %in% expected_pcs & directional_rank <= heat_top_n
)]
pc_anchors <- ranked_heat_results[pc_top_term == TRUE]
membership <- setNames(gene_sets_primary$genes, gene_sets_primary$gene_set_id)
jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
# Up to four integrated representatives per named function, with up to two from
# each direction first, then remaining places by the original FDR/P/|NES| ordering.
# Jaccard >=0.5 in the tested-gene universe removes redundant representatives only;
# it never removes mandatory PC top terms.
pick_functional_terms <- function(dt) {
  dt <- dt[order(sort_key, -abs(NES), gene_set_name)]
  chosen <- character()
  accept <- function(id) {
    if (id %in% chosen) return(FALSE)
    all(vapply(chosen, function(prior) jaccard(membership[[id]], membership[[prior]]) < 0.5, logical(1)))
  }
  for (dir in c("older_aligned", "younger_aligned")) {
    added <- 0L
    for (id in dt[direction == dir, gene_set_id]) {
      if (accept(id)) {
        chosen <- c(chosen, id); added <- added + 1L
      }
      if (added == 2L) break
    }
  }
  for (id in dt$gene_set_id) {
    if (length(chosen) >= functional_n) break
    if (accept(id)) chosen <- c(chosen, id)
  }
  dt[match(chosen, gene_set_id)]
}
functional_candidates <- ranked_heat_results[pc == "integrated_PCA_clock" & significant_q05 & biology_group != "Other"]
functional_representatives <- rbindlist(lapply(
  intersect(biology_group_levels, unique(functional_candidates$biology_group)),
  function(group) pick_functional_terms(functional_candidates[biology_group == group])
))
heat_terms <- unique(c(pc_anchors$gene_set_id, functional_representatives$gene_set_id))
ranked_heat_results[, heatmap_selected := gene_set_id %in% heat_terms]
fwrite(biology_rules, file.path(out_dir, "heatmap_function_assignment_rules.tsv"), sep = "\t")
fwrite(ranked_heat_results[, .(pc, gene_set_id, gene_set_name, collection, biology_group,
  direction, directional_rank, NES, empirical_p, fdr_within_pc, significant_q05,
  significant_q10, pc_top_term, heatmap_selected)],
  file.path(out_dir, "heatmap_all_terms_function_audit.tsv"), sep = "\t")
function_counts <- ranked_heat_results[, .(tested = .N, q05 = sum(significant_q05),
  q10 = sum(significant_q10)), by = .(pc, biology_group, direction)]
fwrite(function_counts, file.path(out_dir, "heatmap_function_support_counts.tsv"), sep = "\t")

heat_selection <- ranked_heat_results[pc == "integrated_PCA_clock" & gene_set_id %in% heat_terms,
  .(gene_set_id, gene_set_name, collection_label, biology_group, integrated_NES = NES,
    integrated_p = empirical_p, integrated_q = fdr_within_pc)]
heat_selection[, pc_top_evidence := vapply(gene_set_id, function(id) {
  rows <- pc_anchors[gene_set_id == id]
  paste(paste(rows$pc, rows$direction, paste0("rank", rows$directional_rank), sep = ":"), collapse = ";")
}, character(1))]
heat_selection[, integrated_function_representative := gene_set_id %in% functional_representatives$gene_set_id]
heat_selection[, selection_reason := fifelse(nzchar(pc_top_evidence), "PC top-five", "")]
heat_selection[integrated_function_representative == TRUE,
  selection_reason := fifelse(nzchar(selection_reason), paste0(selection_reason, "; integrated functional representative"),
                              "Integrated functional representative")]
heat_selection[, biology_group := factor(biology_group, levels = biology_group_levels)]
setorder(heat_selection, biology_group, integrated_q, -integrated_NES, gene_set_name)
heat_selection[, term_label := pretty_term(gene_set_name, width = 57)]
heat_selection[, term_heat_label := make.unique(term_label, sep = " #")]
fwrite(heat_selection, file.path(out_dir, "heatmap_term_selection_audit.tsv"), sep = "\t")

heat_dt_plus_integrated <- merge(ranked_heat_results[heatmap_selected == TRUE],
  heat_selection[, .(gene_set_id, term_label, term_heat_label, selection_reason,
                    integrated_function_representative, pc_top_evidence)],
  by = "gene_set_id", all.x = TRUE, sort = FALSE)
heat_dt_plus_integrated[, `:=`(
  source = fifelse(pc == "integrated_PCA_clock", "integrated_total", "pc_specific"),
  pc_plot = factor(fifelse(pc == "integrated_PCA_clock", "Integrated", pc), levels = c(expected_pcs, "Integrated")),
  biology_group = factor(biology_group, levels = biology_group_levels),
  term_heat_label = factor(term_heat_label, levels = rev(heat_selection$term_heat_label)),
  significance_symbol = fifelse(significant_q05, "*", fifelse(significant_q10, "+", ""))
)]
heat_dt_plus_integrated[, cell_label := significance_symbol]
fwrite(heat_dt_plus_integrated[order(biology_group, term_heat_label, pc_plot),
  .(source, pc, pc_plot, gene_set_id, collection, subcollection, collection_label, gene_set_name,
    gene_set_size_in_universe, NES, direction, empirical_p, fdr_within_pc, directional_rank,
    significant_q05, significant_q10, pc_top_term, significance_symbol, cell_label,
    term_label, term_heat_label, biology_group, selection_reason, integrated_function_representative, pc_top_evidence)],
  file.path(out_dir, "final_pca_pc_plus_integrated_gsea_heatmap_source.tsv"), sep = "\t")

# The full top eight in each direction remain available even if absent from the
# combined heatmap. This table includes function, rank and significance for each PC.
fwrite(ranked_heat_results[pc %in% expected_pcs & directional_rank <= top_n_plot,
  .(pc, direction, directional_rank, biology_group, gene_set_id, gene_set_name,
    NES, empirical_p, fdr_within_pc, significant_q05, significant_q10, heatmap_selected)],
  file.path(out_dir, "top_terms_per_pc_by_function.tsv"), sep = "\t")

plot_function_heatmap <- function(dt, integrated = FALSE) {
  title <- if (integrated) "PC-specific and integrated PCA-clock GSEA" else "Functional groups across active PCA-clock components"
  dt <- copy(dt)
  if (!integrated) dt[, pc_plot := factor(pc, levels = expected_pcs)]
  lim <- max(abs(heat_dt_plus_integrated$NES))
  group_labels <- c("Mitochondrial / energy" = "Mitochondrial /\nenergy", "Fatty acid / lipid" = "Fatty acid /\nlipid",
    "Other metabolism" = "Other\nmetabolism", "Protein / RNA processing" = "Protein / RNA\nprocessing",
    "Adhesion / extracellular matrix" = "Adhesion /\nextracellular matrix", "Neural / synaptic" = "Neural /\nsynaptic",
    "Signalling / development" = "Signalling /\ndevelopment", "Immune / inflammatory" = "Immune /\ninflammatory",
    "Cell cycle / DNA" = "Cell cycle /\nDNA", "Steroid / sterol" = "Steroid / sterol", "Other" = "Other")
  p <- ggplot(dt, aes(x = pc_plot, y = term_heat_label, fill = NES)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = cell_label), size = 2.8, color = "#202020") +
    facet_grid(biology_group ~ ., scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(biology_group = group_labels)) +
    scale_fill_gradient2(low = "#3E7CB1", mid = "white", high = "#C65A4A", midpoint = 0,
                         limits = c(-lim, lim), name = "NES") +
    labs(title = title,
      subtitle = paste0("Top 5 terms per PC in each direction, plus up to 4 significant integrated representatives per function.\n",
                        "* q < 0.05;  + 0.05 <= q < 0.10 (both P < 0.05)."),
      caption = paste0("Positive NES: older-predictive side; negative NES: younger-predictive side. Unmarked tiles: q >= 0.10 or P >= 0.05.\n",
        "BH q is calculated across all tested pathways within each rank. Functional groups are descriptive name-based assignments.\n",
        if (integrated) "Integrated column: GSEA of summed PC gene contributions. " else "",
        "Representative sets require q < 0.05; Jaccard < 0.5 within group."),
      x = NULL, y = NULL) +
    theme_gsea() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
      axis.text.y = element_text(size = 8), axis.ticks = element_blank(), panel.grid = element_blank(),
      panel.spacing.y = unit(0.25, "lines"), strip.placement = "outside",
      strip.background = element_rect(fill = "#f3f3f3", color = NA),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 9),
      plot.subtitle = element_text(size = 10, margin = margin(b = 9)),
      plot.caption = element_text(size = 8, hjust = 0), plot.margin = margin(10, 12, 10, 10)) +
    scale_y_discrete(labels = function(x) sub(" #[0-9]+$", "", x)) +
    scale_x_discrete(position = "top")
  if (integrated) p <- p + geom_vline(xintercept = length(expected_pcs) + 0.5, color = "grey25", linewidth = 0.55)
  p
}
heat_height <- max(13, 3.5 + 0.25 * nrow(heat_selection))
heat_dt <- heat_dt_plus_integrated[pc %in% expected_pcs]
p <- plot_function_heatmap(heat_dt)
plot_files <- c(plot_files, save_plot(p, "05_gsea_nes_heatmap_top_terms", 12.5, heat_height))
p <- plot_function_heatmap(heat_dt_plus_integrated, integrated = TRUE)
plot_files <- c(plot_files, save_plot(p, "13_gsea_nes_heatmap_top_terms_with_integrated_total", 13.5, heat_height))
message("Function-balanced heatmaps: ", nrow(heat_selection), " terms; ", uniqueN(heat_selection$biology_group), " groups.")

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
      subtitle = "Top older- and younger-aligned terms; genes ranked by loading x frozen PC coefficient.",
      x = "Rank in preranked 3000-gene universe",
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
  select_top_by_direction(all_results, "older_aligned", n = 10, pcs = intersect(expected_pcs, c("PC1", "PC2"))),
  select_top_by_direction(all_results, "younger_aligned", n = 10, pcs = intersect(expected_pcs, c("PC1", "PC2")))
), fill = TRUE)
p <- make_top_term_plot(primary, "Primary PC GSEA terms", "PC1 and PC2: the first two variance axes, both with nonzero frozen clock coefficients.")
plot_files <- c(plot_files, save_plot(p, "11_primary_pc_terms", 13, 14))

secondary <- rbindlist(list(
  select_top_by_direction(all_results, "older_aligned", n = 5, pcs = setdiff(expected_pcs, c("PC1", "PC2"))),
  select_top_by_direction(all_results, "younger_aligned", n = 5, pcs = setdiff(expected_pcs, c("PC1", "PC2")))
), fill = TRUE)
p <- make_top_term_plot(secondary, "Secondary PC GSEA terms", "PC5 and PC8: additional active PCs beyond the first two variance axes.")
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
    scale_fill_gradient(low = "white", high = "#4E79A7", breaks = c(0, 0.5, 1)) +
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
    "This analysis interprets the final deployment PCA/glmnet ovarian clock: `training_scope = locked_count_tuning_all_136`, selector `",
    pca_selector,
    "`, refit on all ",
    length(artifact$training_keys),
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
  "`rank_score = pc_loading x pc_coefficient`.",
  "",
  "Positive NES means enrichment on the older-predictive side of the final PCA clock. Negative NES means enrichment on the younger-predictive side. Raw PC loading sign alone was not used as the biological direction.",
  "",
  "An additional integrated PCA-clock GSEA ranks each of the same 3000 PCA input genes once by the summed contribution `sum(pc_loading x pc_coefficient)` across all selected PCs. This integrated analysis is better suited for model-level pathway direction than the PC-specific heatmap.",
  "",
  "## Gene Sets",
  "",
  paste0("Gene universe: ", length(gene_universe), " final selected PCA genes."),
  paste0("MSigDB source: ", paste(unique(msig_selected$gene_set_source), collapse = "; "), "."),
  paste0("MSigDB version(s): ", paste(sort(unique(msig_selected$msigdb_version)), collapse = "; "), "."),
  paste0("Collections used: ", paste(sort(unique(msig_selected$collection_label)), collapse = "; "), "."),
  paste0("Primary size filter: ", min_size, "-", max_size, " genes within the 3000-gene universe. Sensitivity counts were also written for max size ", max_size_sensitivity, "."),
  "",
  "## Permutations",
  "",
  paste0("First-pass weighted-KS GSEA used exponent p = 1, seed ", seed, ", and ", n_perm, " random gene-set permutations per PC and set size."),
  "With 500 first-pass permutations, the minimum empirical P value is about 0.002, so P values and BH q values have coarse Monte Carlo resolution.",
  if (run_refined) paste0("Refined permutations were run for candidate terms with ", n_perm_refined, " permutations.") else "Refined permutations were not run; outputs are marked as first-pass. As in the reference workflow, the empirical tail count is divided by all permutations, while NES normalization uses same-sign null scores; this is a custom weighted-KS analysis, not fgsea.",
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
  "These terms use one aggregate rank vector across PC1, PC2, PC5 and PC8 rather than separate PC-specific rank vectors.",
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
  "- PC1/PC2: first two variance axes; their age correlations are measured in the training cohort.",
  "- PC3/PC4/PC6/PC7/PC9/PC10: zero clock coefficients; retained in PCA diagnostics and excluded from coefficient-oriented GSEA because their rank vectors are all zero.",
  "- PC5/PC8: additional active axes; variance and age correlations are reported without importing biological interpretations from the prior clock.",
  "",
  "## Caution",
  "",
  "PCA loadings are dense and pathway results describe model-oriented gene-set enrichment within the selected 3000-gene clock universe. They are not direct causal gene effects."
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

capture.output(sessionInfo(),file=file.path(root,"provenance/sessionInfo_gsea.txt"))
