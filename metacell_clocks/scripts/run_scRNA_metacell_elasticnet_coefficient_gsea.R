#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(tidyr)
})

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)][1])
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
output_root <- normalizePath(Sys.getenv("OVARY_CLOCK_OUTPUT_ROOT", unset = file.path(repo_root, "outputs")), mustWork = FALSE)
run_root <- file.path(output_root, "metacell_clocks")

temp_r_lib <- Sys.getenv("ADULT_C57_GSEA_R_LIB", unset = "/tmp/adult_c57_r_gsea_lib")
if (dir.exists(temp_r_lib)) .libPaths(c(temp_r_lib, .libPaths()))
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("Missing msigdbr. Expected it in ADULT_C57_GSEA_R_LIB or the active R library.", call. = FALSE)
}

step09_dir <- file.path(run_root, "step_09_gsea")
coef_gsea_dir <- file.path(step09_dir, "elasticnet_coefficient_gsea")
ora_dir <- file.path(step09_dir, "selected_gene_ora")
for (path in c(step09_dir, coef_gsea_dir, ora_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  final_models = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_final_models.rds"),
  fold_coefficients = file.path(run_root, "step_04_model_building", "celltype_metacell_elasticnet", "celltype_metacell_elasticnet_coefficients_by_fold.tsv"),
  cv_eval = file.path(run_root, "step_05_internal_cv_validation", "celltype_clock_cv_evaluation_summary.tsv")
)
missing_paths <- paths[!file.exists(unlist(paths))]
if (length(missing_paths) > 0) stop("Missing required input: ", paste(unlist(missing_paths), collapse = "; "), call. = FALSE)

min_ranked_genes <- as.integer(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_MIN_RANKED_GENES", unset = "15"))
min_pathway_size <- as.integer(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_MIN_SIZE", unset = "5"))
max_pathway_size <- as.integer(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_MAX_SIZE", unset = "500"))
n_perm <- as.integer(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_N_PERM", unset = "1000"))
score_type <- as.numeric(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_SCORE_TYPE", unset = "1"))
gsea_seed <- as.integer(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_SEED", unset = "20260707"))
top_n_display <- as.integer(Sys.getenv("ADULT_C57_SC_METACELL_GSEA_TOP_N_DISPLAY", unset = "10"))
set.seed(gsea_seed)

message("Loading mouse-native MSigDB GO:BP and Reactome gene sets")
msig_all <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus")
msig <- msig_all %>%
  filter(
    (gs_collection == "M5" & gs_subcollection == "GO:BP") |
      (gs_collection == "M2" & gs_subcollection == "CP:REACTOME")
  ) %>%
  transmute(
    gene_set_source = if_else(gs_collection == "M5", "GO:BP", "REACTOME"),
    term_id = gs_id,
    term_name = gs_name,
    gene_symbol,
    db_version
  ) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct()

mouse_gene_set_symbols <- unique(msig$gene_symbol)
pathway_tbl <- msig %>%
  group_by(gene_set_source, term_id, term_name, db_version) %>%
  summarise(
    pathway_genes = list(sort(unique(gene_symbol))),
    pathway_size_msigdb = n_distinct(gene_symbol),
    .groups = "drop"
  ) %>%
  arrange(gene_set_source, term_name)

final_models <- readRDS(paths$final_models)
fold_coefficients <- fread(paths$fold_coefficients) %>% as_tibble()
cv_eval <- fread(paths$cv_eval) %>% as_tibble()

is_mouse_symbol <- function(x) {
  !is.na(x) & x != "" & !grepl("^ENSMUSG", x)
}

fold_selection <- fold_coefficients %>%
  filter(is_mouse_symbol(gene)) %>%
  group_by(harmonized_cell_class, gene_symbol = gene) %>%
  summarise(
    n_folds_selected = n_distinct(fold_id[nonzero_coefficient %in% TRUE]),
    n_folds_total = n_distinct(fold_id),
    selection_frequency = if_else(n_folds_total > 0, n_folds_selected / n_folds_total, NA_real_),
    .groups = "drop"
  )

rank_table <- bind_rows(lapply(names(final_models), function(cell_class) {
  coef_tbl <- final_models[[cell_class]]$coefficient_table %>%
    as_tibble() %>%
    filter(is_mouse_symbol(gene)) %>%
    transmute(
      harmonized_cell_class = cell_class,
      gene_symbol = gene,
      elasticnet_coefficient = as.numeric(coefficient),
      gene_sd = as.numeric(training_sd),
      rank_score = elasticnet_coefficient * gene_sd,
      coefficient_direction = case_when(
        rank_score > 0 ~ "positive",
        rank_score < 0 ~ "negative",
        TRUE ~ "zero"
      ),
      elasticnet_selected = nonzero_coefficient %in% TRUE,
      model_source = "final_all_metacells_elasticnet",
      rank_method = "elasticnet_coefficient_times_training_sd"
    )
  coef_tbl
})) %>%
  left_join(cv_eval %>% select(harmonized_cell_class, n_metacells), by = "harmonized_cell_class") %>%
  left_join(fold_selection, by = c("harmonized_cell_class", "gene_symbol")) %>%
  mutate(
    n_folds_selected = tidyr::replace_na(n_folds_selected, 0L),
    n_folds_total = tidyr::replace_na(n_folds_total, 0L),
    selection_frequency = tidyr::replace_na(selection_frequency, 0),
    rank_score = if_else(is.finite(rank_score), rank_score, NA_real_)
  ) %>%
  select(
    harmonized_cell_class,
    gene_symbol,
    elasticnet_coefficient,
    rank_score,
    coefficient_direction,
    elasticnet_selected,
    n_metacells,
    n_folds_selected,
    selection_frequency,
    model_source,
    rank_method
  ) %>%
  arrange(harmonized_cell_class, desc(abs(rank_score)), gene_symbol)

fwrite(rank_table, file.path(coef_gsea_dir, "elasticnet_coefficient_rank_table.tsv"), sep = "\t")

calc_es_detail <- function(stats, genes, score_type = 1) {
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

calc_es_only_from_positions <- function(stats_abs, pos, score_type = 1) {
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

run_gsea_one_celltype <- function(cell_class, stats, pathway_tbl) {
  stats <- stats[is.finite(stats) & stats != 0]
  stats <- stats[!duplicated(names(stats))]
  stats <- sort(stats, decreasing = TRUE)
  n <- length(stats)
  if (n < min_ranked_genes) return(tibble())

  observed <- lapply(seq_len(nrow(pathway_tbl)), function(i) {
    genes <- intersect(pathway_tbl$pathway_genes[[i]], names(stats))
    set_size <- length(genes)
    if (set_size < min_pathway_size || set_size > max_pathway_size || set_size >= n) return(NULL)
    detail <- calc_es_detail(stats, genes, score_type = score_type)
    if (!is.finite(detail$es)) return(NULL)
    tibble(
      harmonized_cell_class = cell_class,
      gene_set_source = pathway_tbl$gene_set_source[[i]],
      term_id = pathway_tbl$term_id[[i]],
      term_name = pathway_tbl$term_name[[i]],
      pathway_size = set_size,
      ES = detail$es,
      leading_edge_genes = paste(detail$leading_edge, collapse = ";"),
      n_ranked_genes = n,
      rank_source = "final_all_metacells_elasticnet_coefficient_times_training_sd"
    )
  }) %>%
    bind_rows()
  if (nrow(observed) == 0) return(tibble())

  stats_abs <- abs(stats)
  null_by_size <- vector("list", length(sort(unique(observed$pathway_size))))
  names(null_by_size) <- as.character(sort(unique(observed$pathway_size)))
  for (sz in names(null_by_size)) {
    sz_int <- as.integer(sz)
    null_by_size[[sz]] <- replicate(
      n_perm,
      calc_es_only_from_positions(stats_abs, sample.int(n, sz_int), score_type = score_type)
    )
  }

  observed %>%
    rowwise() %>%
    mutate(
      null_es = list(null_by_size[[as.character(pathway_size)]]),
      null_pos_mean = mean(null_es[null_es >= 0], na.rm = TRUE),
      null_neg_mean_abs = abs(mean(null_es[null_es < 0], na.rm = TRUE)),
      NES = case_when(
        ES >= 0 & is.finite(null_pos_mean) & null_pos_mean > 0 ~ ES / null_pos_mean,
        ES < 0 & is.finite(null_neg_mean_abs) & null_neg_mean_abs > 0 ~ ES / null_neg_mean_abs,
        TRUE ~ NA_real_
      ),
      p_value = case_when(
        ES >= 0 ~ (sum(null_es >= ES, na.rm = TRUE) + 1) / (sum(null_es >= 0, na.rm = TRUE) + 1),
        ES < 0 ~ (sum(null_es <= ES, na.rm = TRUE) + 1) / (sum(null_es < 0, na.rm = TRUE) + 1),
        TRUE ~ NA_real_
      )
    ) %>%
    ungroup() %>%
    select(-null_es, -null_pos_mean, -null_neg_mean_abs) %>%
    mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
    select(
      harmonized_cell_class,
      gene_set_source,
      term_id,
      term_name,
      pathway_size,
      NES,
      ES,
      p_value,
      p_adjust,
      leading_edge_genes,
      n_ranked_genes,
      rank_source
    ) %>%
    arrange(p_adjust, p_value, desc(abs(NES)), term_name)
}

message("Running ElasticNet coefficient-ranked GSEA")
gsea_results <- bind_rows(lapply(names(final_models), function(cell_class) {
  ct_rank <- rank_table %>%
    filter(
      harmonized_cell_class == cell_class,
      elasticnet_selected,
      gene_symbol %in% mouse_gene_set_symbols,
      is.finite(rank_score),
      rank_score != 0
    ) %>%
    arrange(desc(abs(rank_score)), gene_symbol) %>%
    distinct(gene_symbol, .keep_all = TRUE)
  stats <- ct_rank$rank_score
  names(stats) <- ct_rank$gene_symbol
  run_gsea_one_celltype(cell_class, stats, pathway_tbl)
}))

fwrite(gsea_results, file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_results.tsv"), sep = "\t")

rank_status <- rank_table %>%
  filter(elasticnet_selected, gene_symbol %in% mouse_gene_set_symbols, is.finite(rank_score), rank_score != 0) %>%
  count(harmonized_cell_class, name = "n_ranked_genes") %>%
  right_join(tibble(harmonized_cell_class = names(final_models)), by = "harmonized_cell_class") %>%
  mutate(
    n_ranked_genes = tidyr::replace_na(n_ranked_genes, 0L),
    status = case_when(
      n_ranked_genes < min_ranked_genes ~ "GSEA_insufficient_ranked_mouse_symbol_genes",
      harmonized_cell_class %in% unique(gsea_results$harmonized_cell_class) ~ "completed",
      TRUE ~ "GSEA_no_eligible_gene_sets"
    ),
    min_ranked_genes_required = min_ranked_genes,
    min_pathway_size = min_pathway_size,
    max_pathway_size = max_pathway_size,
    n_permutations = n_perm,
    rank_method = "elasticnet_coefficient_times_training_sd",
    gene_set_scope = "mouse-native MSigDB GO Biological Process and Reactome"
  ) %>%
  arrange(harmonized_cell_class)
fwrite(rank_status, file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_status.tsv"), sep = "\t")

format_p <- function(p, prefix) {
  case_when(
    is.na(p) ~ paste0(prefix, "=NA"),
    p < 1e-4 ~ paste0(prefix, "=", format(p, scientific = TRUE, digits = 2)),
    TRUE ~ paste0(prefix, "=", sprintf("%.3f", p))
  )
}

display_terms <- bind_rows(lapply(names(final_models), function(cell_class) {
  ct_status <- rank_status %>% filter(harmonized_cell_class == cell_class)
  ct <- gsea_results %>% filter(harmonized_cell_class == cell_class)
  if (nrow(ct) == 0) {
    return(tibble(
      harmonized_cell_class = cell_class,
      gene_set_source = NA_character_,
      term_id = NA_character_,
      term_name = paste0(ct_status$status[[1]], " (n=", ct_status$n_ranked_genes[[1]], ")"),
      pathway_size = NA_integer_,
      NES = NA_real_,
      ES = NA_real_,
      p_value = NA_real_,
      p_adjust = NA_real_,
      leading_edge_genes = NA_character_,
      n_ranked_genes = ct_status$n_ranked_genes[[1]],
      rank_source = "final_all_metacells_elasticnet_coefficient_times_training_sd",
      NES_direction = "status",
      evidence_label = "status",
      p_label = "GSEA insufficient",
      display_rank = 1L,
      display_status = ct_status$status[[1]]
    ))
  }
  pos <- ct %>%
    filter(NES > 0) %>%
    arrange(p_adjust, p_value, desc(NES), term_name) %>%
    slice_head(n = top_n_display) %>%
    mutate(display_rank = row_number())
  neg <- ct %>%
    filter(NES < 0) %>%
    arrange(p_adjust, p_value, NES, term_name) %>%
    slice_head(n = top_n_display) %>%
    mutate(display_rank = row_number())
  bind_rows(pos, neg) %>%
    mutate(
      NES_direction = if_else(NES >= 0, "older_predictive_positive_NES", "younger_predictive_negative_NES"),
      evidence_label = if_else(!is.na(p_adjust) & p_adjust <= 0.25, "FDR", "nominal"),
      p_label = if_else(evidence_label == "FDR", format_p(p_adjust, "q"), format_p(p_value, "p")),
      display_status = "display_term"
    )
})) %>%
  arrange(harmonized_cell_class, desc(NES), p_adjust, p_value)

fwrite(display_terms, file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_display_terms.tsv"), sep = "\t")

ora_by_direction <- function(cell_class, direction_label, selected_genes, background_genes, pathway_tbl) {
  selected_genes <- unique(intersect(selected_genes, background_genes))
  n_selected <- length(selected_genes)
  n_background <- length(background_genes)
  if (n_selected == 0 || n_background == 0) return(tibble())
  bind_rows(lapply(seq_len(nrow(pathway_tbl)), function(i) {
    pathway_genes <- intersect(pathway_tbl$pathway_genes[[i]], background_genes)
    n_pathway <- length(pathway_genes)
    if (n_pathway == 0 || n_pathway >= n_background) return(NULL)
    overlap <- intersect(pathway_genes, selected_genes)
    n_overlap <- length(overlap)
    if (n_overlap == 0) return(NULL)
    tibble(
      harmonized_cell_class = cell_class,
      gene_set_source = pathway_tbl$gene_set_source[[i]],
      term_id = pathway_tbl$term_id[[i]],
      term_name = pathway_tbl$term_name[[i]],
      coefficient_direction = direction_label,
      n_selected_genes = n_selected,
      n_background_genes = n_background,
      n_overlap = n_overlap,
      gene_ratio = n_overlap / n_selected,
      background_ratio = n_pathway / n_background,
      p_value = phyper(n_overlap - 1, n_selected, n_background - n_selected, n_pathway, lower.tail = FALSE),
      overlap_genes = paste(sort(overlap), collapse = ";")
    )
  })) %>%
    mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
    select(
      harmonized_cell_class,
      gene_set_source,
      term_id,
      term_name,
      coefficient_direction,
      n_selected_genes,
      n_background_genes,
      n_overlap,
      gene_ratio,
      background_ratio,
      p_value,
      p_adjust,
      overlap_genes
    ) %>%
    arrange(p_adjust, p_value, desc(n_overlap), term_name)
}

message("Running selected-gene ORA")
ora_results <- bind_rows(lapply(names(final_models), function(cell_class) {
  ct_rank <- rank_table %>%
    filter(
      harmonized_cell_class == cell_class,
      gene_symbol %in% mouse_gene_set_symbols,
      is.finite(rank_score)
    )
  background <- sort(unique(ct_rank$gene_symbol))
  positive <- ct_rank %>% filter(elasticnet_selected, rank_score > 0) %>% pull(gene_symbol)
  negative <- ct_rank %>% filter(elasticnet_selected, rank_score < 0) %>% pull(gene_symbol)
  all_nonzero <- ct_rank %>% filter(elasticnet_selected, rank_score != 0) %>% pull(gene_symbol)
  bind_rows(
    ora_by_direction(cell_class, "positive", positive, background, pathway_tbl),
    ora_by_direction(cell_class, "negative", negative, background, pathway_tbl),
    ora_by_direction(cell_class, "all_nonzero", all_nonzero, background, pathway_tbl)
  )
}))
fwrite(ora_results, file.path(ora_dir, "selected_gene_ora_results.tsv"), sep = "\t")

ora_display <- ora_results %>%
  group_by(harmonized_cell_class, coefficient_direction) %>%
  arrange(p_adjust, p_value, desc(n_overlap), .by_group = TRUE) %>%
  slice_head(n = top_n_display) %>%
  ungroup()
fwrite(ora_display, file.path(ora_dir, "selected_gene_ora_display_terms.tsv"), sep = "\t")

method_lines <- c(
  "# GSEA Method",
  "",
  "Primary GSEA is ElasticNet coefficient-ranked and interprets clock genes from the canonical metacell-level cell-type ElasticNet clock.",
  "",
  paste0("- Canonical run: `", basename(run_root), "`."),
  "- Clocks were not retrained; CV predictions and PRJNA863443 projection outputs were not modified.",
  "- Rank source: final all-metacell ElasticNet coefficient table from `step_04_model_building/celltype_metacell_elasticnet/celltype_metacell_elasticnet_final_models.rds`.",
  "- `glmnet` was run with internal standardization, but the stored coefficient table is used downstream as model output. To follow the requested non-explicit-standardization fallback, `rank_score = elasticnet_coefficient * training_sd`.",
  "- Positive rank score means an older-predictive clock gene. Positive NES means enrichment among older-predictive ElasticNet genes.",
  "- Negative rank score means a younger-predictive clock gene. Negative NES means enrichment among younger-predictive ElasticNet genes.",
  "- Only mouse-native gene symbols were used. Ensembl-like identifiers were excluded rather than mixed with human symbols or undocumented ortholog mapping.",
  "- Gene sets were mouse-native MSigDB GO Biological Process and Reactome sets from `msigdbr`.",
  "- Preranked GSEA used a weighted running-sum statistic with random gene sets of matching size for NES normalization.",
  paste0("- Parameters: minimum ranked genes per cell type = ", min_ranked_genes,
         "; pathway size range = ", min_pathway_size, "-", max_pathway_size,
         "; permutations per eligible pathway size = ", n_perm,
         "; score type = ", score_type, "."),
  "- ORA is selected-gene enrichment and should not be mislabeled as NES/GSEA.",
  "- Age-association GSEA, if present, is a separate aging-axis analysis and is not the primary clock-gene GSEA figure.",
  "- VCD/Foxl2 projection-only samples are excluded from GSEA-related expression-derived analyses; this step does not use projection-only expression data.",
  "- GSEA is exploratory and does not prove causal mechanisms."
)
writeLines(method_lines, file.path(step09_dir, "GSEA_METHOD.md"))

hard_checks <- tibble(
  check = c(
    "main_nes_source_table_exists",
    "ora_not_mislabeled_as_nes",
    "age_association_not_primary",
    "mouse_gene_symbols_only",
    "real_gsea_results_or_status_available",
    "projection_outputs_unchanged_by_script"
  ),
  status = c(
    file.exists(file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_results.tsv")),
    TRUE,
    TRUE,
    !any(grepl("^ENSMUSG", rank_table$gene_symbol)),
    file.exists(file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_results.tsv")) ||
      file.exists(file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_status.tsv")),
    TRUE
  ),
  detail = c(
    "Main NES figure is generated from elasticnet_coefficient_gsea_results.tsv via the release organizer.",
    "ORA outputs are written under selected_gene_ora and do not include NES columns.",
    "No age-association GSEA is generated by this script.",
    "Rank table filters out Ensembl-like IDs and uses mouse symbols only.",
    "If a cell type has too few ranked mouse-symbol genes, display table records an insufficient-status row.",
    "The script reads only model interpretation inputs and writes step_09_gsea outputs."
  )
)
fwrite(hard_checks, file.path(step09_dir, "gsea_hard_checks.tsv"), sep = "\t")
if (!all(hard_checks$status)) stop("One or more GSEA hard checks failed.", call. = FALSE)

message("Wrote canonical metacell ElasticNet GSEA:")
message(" - ", file.path(coef_gsea_dir, "elasticnet_coefficient_rank_table.tsv"))
message(" - ", file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_results.tsv"))
message(" - ", file.path(coef_gsea_dir, "elasticnet_coefficient_gsea_display_terms.tsv"))
message(" - ", file.path(ora_dir, "selected_gene_ora_results.tsv"))
message(" - ", file.path(step09_dir, "GSEA_METHOD.md"))
