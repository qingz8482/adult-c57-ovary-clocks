#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else getwd()
if (!dir.exists(script_dir) || !file.exists(file.path(script_dir, "03_train_bulk_pseudobulk_clocks.R"))) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
data_input_dir <- normalizePath(file.path(project_dir, "..", "data_input"), mustWork = TRUE)

metadata_path <- file.path(project_dir, "metadata", "bulk_pseudobulk_sample_metadata.csv")
stopifnot(file.exists(metadata_path))

out_dir <- file.path(project_dir, "results", "alignment_audit", "bulk_group_level_alignment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

target_datasets <- c("GSE273724", "GSE279446", "GSE279573")
metadata_candidate_cols <- c(
  "sample_key", "sample_id", "accession", "run", "SRR", "group",
  "condition", "treatment", "age", "age_week"
)

sanitize_id <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("[^a-z0-9]+", "", x)
}

na_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

bool_chr <- function(x) {
  ifelse(isTRUE(x), "TRUE", "FALSE")
}

read_matrix_header <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  names(fread(path, nrows = 0))
}

read_gene_ids <- function(path) {
  if (is.na(path) || !file.exists(path)) return(character())
  dt <- fread(path, select = 1, showProgress = FALSE)
  unique(trimws(as.character(dt[[1]])))
}

detect_expression_scale <- function(path) {
  if (is.na(path) || !file.exists(path)) return("unknown")
  dt <- fread(path, nrows = 200, showProgress = FALSE)
  if (ncol(dt) <= 1) return("unknown")
  vals <- suppressWarnings(as.numeric(unlist(dt[, -1, with = FALSE], use.names = FALSE)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return("unknown")
  integer_like <- mean(abs(vals - round(vals)) < 1e-8) > 0.99
  if (integer_like && max(vals, na.rm = TRUE) > 50) return("raw_counts")
  if (min(vals, na.rm = TRUE) >= 0 && max(vals, na.rm = TRUE) < 40) return("logcpm_or_normalized")
  "unknown"
}

detect_gene_style <- function(genes) {
  genes <- genes[!is.na(genes) & genes != ""]
  if (length(genes) == 0) return("unknown")
  ens_version <- mean(grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", genes)) > 0.5
  ens_plain <- mean(grepl("^ENS[A-Z]*G[0-9]+$", genes)) > 0.5
  mouse_symbol <- mean(grepl("^[A-Z][A-Za-z0-9._-]*$", genes) & !grepl("^ENS", genes)) > 0.5
  if (ens_version) return("versioned_Ensembl")
  if (ens_plain) return("Ensembl")
  if (mouse_symbol) return("gene_symbol")
  "unknown"
}

clean_gene_ids_for_overlap <- function(genes) {
  genes <- trimws(as.character(genes))
  if (mean(grepl("^ENS[A-Z]*G[0-9]+\\.[0-9]+$", genes), na.rm = TRUE) > 0.5) {
    genes <- sub("\\.[0-9]+$", "", genes)
  }
  unique(genes[!is.na(genes) & genes != ""])
}

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
    maps[[length(maps) + 1L]] <- dt[, .(gene_id = na_chr(V1), gene_symbol = na_chr(V2))]
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
  mapping <- mapping[, .N, by = .(gene_id, gene_symbol)]
  setorder(mapping, gene_id, -N, gene_symbol)
  mapping <- mapping[, .SD[1], by = gene_id][, .(gene_id, gene_symbol)]
  gene_id_symbol_map_cache <<- mapping
  mapping
}

map_gene_ids_for_overlap <- function(genes) {
  clean <- clean_gene_ids_for_overlap(genes)
  if (!any(grepl("^ENSMUSG", clean))) return(clean)
  mapping <- load_gene_id_symbol_map()
  if (nrow(mapping) == 0) return(clean)
  symbol_lookup <- mapping$gene_symbol
  names(symbol_lookup) <- mapping$gene_id
  mapped <- unname(symbol_lookup[clean])
  clean[!is.na(mapped) & mapped != ""] <- mapped[!is.na(mapped) & mapped != ""]
  unique(clean[!is.na(clean) & clean != ""])
}

load_optional_genes <- function(paths) {
  genes <- character()
  for (path in paths[file.exists(paths)]) {
    dt <- fread(path, showProgress = FALSE)
    gene_col <- intersect(c("gene_id", "gene", "feature", names(dt)[1]), names(dt))[1]
    genes <- c(genes, as.character(dt[[gene_col]]))
  }
  unique(trimws(genes[!is.na(genes) & genes != ""]))
}

common_gene_paths <- c(
  file.path(project_dir, "results", "integrated_bulk_pseudobulk_technical_gene_universe.tsv")
)
common_genes <- load_optional_genes(common_gene_paths)

selected_gene_paths <- list.files(
  file.path(project_dir, "results", "selected_genes"),
  pattern = "selected_genes\\.tsv$",
  full.names = TRUE
)
selected_genes <- load_optional_genes(selected_gene_paths)

meta <- fread(metadata_path)
meta <- meta[modality == "bulk_rnaseq" & dataset_folder %chin% target_datasets]
if (nrow(meta) == 0) stop("No target bulk datasets were found in metadata.", call. = FALSE)

external_mapping_paths <- c(
  file.path(project_dir, "results", "external_projection", "alignment_audit", "manual_projection_sample_mapping_APPROVED_AUTO.tsv"),
  file.path(project_dir, "results", "external_projection", "alignment_audit", "manual_projection_sample_mapping_template.tsv")
)
external_mapping_paths <- external_mapping_paths[file.exists(external_mapping_paths)]
external_mapping <- rbindlist(lapply(external_mapping_paths, function(path) {
  dt <- fread(path, sep = "\t")
  dt[, manual_mapping_source := "external_projection_mapping_candidate"]
  dt[, manual_mapping_source_file := path]
  dt
}), fill = TRUE)
if (nrow(external_mapping) > 0 && "dataset_folder" %in% names(external_mapping)) {
  external_mapping <- external_mapping[dataset_folder %chin% target_datasets]
}

audit_rows <- list()
template_rows <- list()

metadata_label <- function(row) {
  vals <- c(row$sample_title %||% NA_character_, row$treatment %||% NA_character_, row$sample_id %||% NA_character_)
  vals <- vals[!is.na(vals) & vals != ""]
  if (length(vals) == 0) NA_character_ else vals[[1]]
}

candidate_expr_for_row <- function(row, expr_cols) {
  accession_norm <- sanitize_id(row$accession %||% "")
  expr_norm <- sanitize_id(expr_cols)
  if (grepl("wtovary", accession_norm)) {
    return(expr_cols[grepl("(^|zeb1)wtovary", expr_norm)])
  }
  if (grepl("zeb1200h", accession_norm)) {
    return(expr_cols[grepl("zeb1200h", expr_norm)])
  }
  if (grepl("zeb1200m", accession_norm)) {
    return(expr_cols[grepl("zeb1200m", expr_norm)])
  }
  raw_terms <- c(row$sample_id, row$accession, row$treatment, row$sample_title)
  tokens <- unique(unlist(strsplit(tolower(paste(raw_terms, collapse = " ")), "[^a-z0-9]+")))
  tokens <- tokens[!is.na(tokens) & (nchar(tokens) >= 3 | tokens %chin% c("wt"))]
  tokens <- tokens[!tokens %chin% c("gse", "sample", "ovary", "subset", "week", "mice", "mouse")]
  hits <- character()
  for (token in tokens) {
    token_norm <- sanitize_id(token)
    if (nchar(token_norm) >= 3) {
      hits <- union(hits, expr_cols[grepl(token_norm, expr_norm, fixed = TRUE)])
    }
  }
  hits
}

append_template_row <- function(dataset_folder, meta_row, expression_column, reason,
                                source = "audit_candidate", aggregation_id = NA_character_,
                                aggregation_method = "unsafe_do_not_use",
                                sample_unit_type = "unsafe_unknown",
                                include = FALSE, output_sample_key = NA_character_,
                                note = NA_character_) {
  data.table(
    dataset_folder = dataset_folder,
    metadata_row_id = meta_row$sample_id,
    metadata_sample_key = meta_row$sample_key,
    metadata_sample_id = meta_row$sample_id,
    metadata_group_label = metadata_label(meta_row),
    metadata_age_week = meta_row$age_week,
    metadata_condition = meta_row$sample_context %||% meta_row$treatment %||% NA_character_,
    metadata_control_like_clock_candidate = bool_chr(isTRUE(meta_row$control_like_clock_candidate)),
    expression_column = expression_column,
    expression_column_sanitized = sanitize_id(expression_column),
    proposed_match_reason = reason,
    sample_unit_type = sample_unit_type,
    aggregation_id = aggregation_id %||% meta_row$sample_key,
    aggregation_method = aggregation_method,
    output_sample_key = output_sample_key %||% paste(meta_row$sample_key, sanitize_id(expression_column), sep = "__"),
    include_in_primary_clock_TRUE_FALSE = bool_chr(include),
    manual_approval_TRUE_FALSE = "FALSE",
    manual_note = note %||% "",
    manual_mapping_source = source
  )
}

for (ds in target_datasets) {
  meta_ds <- meta[dataset_folder == ds]
  if (nrow(meta_ds) == 0) next
  log_path <- unique(meta_ds$expression_logcpm_path)[1]
  count_path <- unique(meta_ds$expression_counts_path)[1]
  header <- read_matrix_header(log_path)
  expr_cols <- header[-1]
  genes <- read_gene_ids(log_path)
  clean_genes <- map_gene_ids_for_overlap(genes)
  common_overlap <- length(intersect(clean_genes, common_genes))
  selected_overlap <- length(intersect(clean_genes, selected_genes))
  exact_values <- unique(unlist(meta_ds[, intersect(c("sample_key", "sample_id", "accession"), names(meta_ds)), with = FALSE], use.names = FALSE))
  exact_values <- exact_values[!is.na(exact_values) & exact_values != ""]
  exact_match_count <- sum(exact_values %chin% expr_cols)
  sanitized_match_count <- sum(sanitize_id(exact_values) %chin% sanitize_id(expr_cols))
  group_candidates <- lapply(seq_len(nrow(meta_ds)), function(i) {
    row <- meta_ds[i]
    hits <- candidate_expr_for_row(row, expr_cols)
    paste(row$sample_key, paste(hits, collapse = ";"), sep = "=")
  })
  includeability <- "needs_manual_group_mapping"
  if (nrow(meta_ds) == length(expr_cols) && exact_match_count == nrow(meta_ds)) {
    includeability <- "ready_exact_1to1"
  }
  if (length(clean_genes) < 5000 || (length(common_genes) > 0 && common_overlap < 5000)) {
    includeability <- "gene_id_problem"
  }
  if (ds == "GSE279446" && (length(common_genes) > 0 && common_overlap < 5000)) {
    includeability <- "gene_id_problem"
  }
  audit_rows[[length(audit_rows) + 1L]] <- data.table(
    dataset_folder = ds,
    n_curated_metadata_rows = nrow(meta_ds),
    n_expression_columns = length(expr_cols),
    metadata_candidate_id_columns = paste(intersect(metadata_candidate_cols, names(meta_ds)), collapse = ";"),
    expression_column_names = paste(expr_cols, collapse = ";"),
    exact_match_count = exact_match_count,
    sanitized_match_count = sanitized_match_count,
    group_label_match_candidates = paste(unlist(group_candidates), collapse = " | "),
    expression_appears = detect_expression_scale(count_path),
    log_expression_appears = detect_expression_scale(log_path),
    gene_id_style = detect_gene_style(genes),
    n_genes_after_cleaning = length(clean_genes),
    overlap_with_integrated_gene_universe = common_overlap,
    overlap_with_old_deployment_selected_genes = selected_overlap,
    selected_gene_coverage = if (length(selected_genes) > 0) selected_overlap / length(selected_genes) else NA_real_,
    includeability_status = includeability
  )

  ext_ds <- if (nrow(external_mapping) > 0) external_mapping[dataset_folder == ds] else data.table()
  if (nrow(ext_ds) > 0 && all(c("metadata_sample_key", "candidate_expression_column") %in% names(ext_ds))) {
    for (i in seq_len(nrow(ext_ds))) {
      map_row <- ext_ds[i]
      meta_row <- meta_ds[sample_key == map_row$metadata_sample_key]
      if (nrow(meta_row) == 0) next
      cols <- trimws(unlist(strsplit(as.character(map_row$candidate_expression_column), ";", fixed = TRUE)))
      cols <- cols[cols != ""]
      for (expr_col in cols) {
        include <- isTRUE(meta_row$control_like_clock_candidate) && includeability != "gene_id_problem"
        template_rows[[length(template_rows) + 1L]] <- append_template_row(
          ds, meta_row, expr_col,
          reason = paste("external projection candidate:", map_row$match_tier %||% map_row$match_status %||% "manual"),
          source = "external_projection_mapping_candidate",
          aggregation_id = meta_row$sample_key,
          aggregation_method = "mean_logcpm",
          sample_unit_type = "biological_replicate",
          include = include,
          output_sample_key = paste0(meta_row$sample_key, "__group_mean"),
          note = map_row$manual_note %||% "Review before benchmark training."
        )
      }
    }
  } else {
    for (i in seq_len(nrow(meta_ds))) {
      row <- meta_ds[i]
      hits <- candidate_expr_for_row(row, expr_cols)
      if (length(hits) == 0) {
        template_rows[[length(template_rows) + 1L]] <- append_template_row(
          ds, row, NA_character_, "no safe automatic expression-column candidate"
        )
      } else {
        for (expr_col in hits) {
          template_rows[[length(template_rows) + 1L]] <- append_template_row(
            ds, row, expr_col, "group-label candidate; manual review required"
          )
        }
      }
    }
  }
}

audit <- rbindlist(audit_rows, fill = TRUE)
template <- rbindlist(template_rows, fill = TRUE)
setorder(audit, dataset_folder)
setorder(template, dataset_folder, metadata_sample_key, expression_column)

audit_path <- file.path(out_dir, "bulk_group_alignment_audit.tsv")
template_path <- file.path(out_dir, "bulk_group_manual_mapping_template.tsv")
approved_path <- file.path(out_dir, "bulk_group_manual_mapping_APPROVED.tsv")
fwrite(audit, audit_path, sep = "\t")
fwrite(template, template_path, sep = "\t")

approved <- template[
  manual_approval_TRUE_FALSE == "TRUE" &
    include_in_primary_clock_TRUE_FALSE == "TRUE" &
    !is.na(expression_column) &
    expression_column != ""
]
existing_approved <- if (file.exists(approved_path)) fread(approved_path, sep = "\t") else data.table()
existing_ready <- FALSE
if (nrow(existing_approved) > 0 && all(c(
  "dataset_folder", "include_in_primary_clock_TRUE_FALSE", "manual_approval_TRUE_FALSE", "expression_column"
) %in% names(existing_approved))) {
  existing_approved[, include_in_primary_clock_TRUE_FALSE := tolower(as.character(include_in_primary_clock_TRUE_FALSE)) %in% c("true", "t", "1")]
  existing_approved[, manual_approval_TRUE_FALSE := tolower(as.character(manual_approval_TRUE_FALSE)) %in% c("true", "t", "1")]
  existing_ready <- all(target_datasets %chin% existing_approved[
    include_in_primary_clock_TRUE_FALSE == TRUE &
      manual_approval_TRUE_FALSE == TRUE &
      !is.na(expression_column) &
      expression_column != "",
    dataset_folder
  ])
}

if (existing_ready) {
  message("Existing approved mapping is present and was preserved: ", approved_path)
} else if (length(target_datasets) > 0 && all(target_datasets %chin% approved$dataset_folder)) {
  fwrite(template, approved_path, sep = "\t")
  message("All target datasets have approved mapping rows; wrote: ", approved_path)
} else {
  missing <- setdiff(target_datasets, unique(approved$dataset_folder))
  message("Wrote audit: ", audit_path)
  message("Wrote manual mapping template: ", template_path)
  message(
    "Manual review is required before benchmark training can include: ",
    paste(missing, collapse = ", ")
  )
}
