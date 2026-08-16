#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(glmnet)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260728)

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)][1])
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
data_root <- normalizePath(Sys.getenv("OVARY_CLOCK_DATA_ROOT", unset = file.path(repo_root, "data")), mustWork = FALSE)
output_root <- normalizePath(Sys.getenv("OVARY_CLOCK_OUTPUT_ROOT", unset = file.path(repo_root, "outputs")), mustWork = FALSE)
analysis_dir <- data_root
run_root <- file.path(output_root, "metacell_clocks")
model_path <- file.path(
  run_root, "step_04_model_building", "celltype_metacell_elasticnet",
  "celltype_metacell_elasticnet_final_models.rds"
)
cv_path <- file.path(
  run_root, "step_05_internal_cv_validation",
  "celltype_clock_cv_evaluation_summary.tsv"
)
sample_metadata_path <- file.path(
  data_root, "scRNA_dataset", "GSE235021", "metadata",
  "GSE235021_scRNA_sample_metadata.tsv"
)
output_dir <- file.path(
  run_root, "step_08_figures_tables_release", "figures",
  "05_stress_test_projection", "GSE235021"
)
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_path <- file.path(log_dir, "GSE235021_metacell_projection.log")
log_connection <- file(log_path, open = "wt")
sink(log_connection, type = "output", split = TRUE)
sink(log_connection, type = "message")
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_connection), silent = TRUE)
}, add = TRUE)

message("GSE235021 frozen metacell-clock projection started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Run root: ", run_root)
message("Output directory: ", output_dir)

required_paths <- c(model_path, cv_path, sample_metadata_path)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = "; "), call. = FALSE)
}

confidence_cutoff <- 0.70
target_cells_per_metacell <- 50L
minimum_cells_per_metacell <- 25L
minimum_features <- 200L
maximum_mitochondrial_percent <- 20
dataset_folder <- "GSE235021"
dataset_id <- "2026_Tu_GSE235021_whole_grain_rice_ovarian_aging"
dataset_name <- "Tu 2026 whole-grain rice ovarian aging"

marker_sets <- list(
  Granulosa = c("Amh", "Fshr", "Foxl2", "Cyp19a1", "Inha", "Inhbb"),
  Theca_Stroma = c("Col1a1", "Col1a2", "Dcn", "Lum", "Lhcgr", "Cyp17a1", "Star"),
  Endothelial = c("Pecam1", "Cdh5", "Kdr", "Vwf", "Eng"),
  Epithelial = c("Epcam", "Krt8", "Krt18", "Krt19"),
  Immune_Myeloid = c("Ptprc", "Lyz2", "Csf1r", "C1qa", "C1qb", "Adgre1"),
  T_NK = c("Cd3d", "Cd3e", "Trac", "Nkg7", "Klrb1c", "Gzma"),
  B_Plasma = c("Ms4a1", "Cd79a", "Cd79b", "Jchain", "Mzb1", "Xbp1"),
  Pericyte_Smooth_Muscle = c("Rgs5", "Pdgfrb", "Acta2", "Tagln", "Myh11"),
  Oocyte_Germline = c("Ddx4", "Gdf9", "Zp1", "Zp2", "Zp3", "Nobox"),
  Erythroid = c("Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2")
)

class_palette <- c(
  B_Plasma = "#1B9E77",
  Endothelial = "#00A6D6",
  Epithelial = "#7CAE00",
  Granulosa = "#F05AA6",
  Immune_Myeloid = "#D95F02",
  Pericyte_Smooth_Muscle = "#6A3D9A",
  T_NK = "#0072B2",
  Theca_Stroma = "#C77CFF"
)
diet_palette <- c(
  polished_rice = "#9E9E9E",
  black_rice = "#3B4CC0",
  white_brown_rice = "#B85C38"
)
diet_labels <- c(
  polished_rice = "Polished rice",
  black_rice = "Black rice",
  white_brown_rice = "White/brown rice"
)

theme_projection <- function(base_size = 8.5) {
  theme_bw(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.22, colour = "grey92"),
      strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.35),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.key = element_blank()
    )
}

save_projection_plot <- function(plot, stem, width, height, plot_data) {
  png_path <- file.path(output_dir, paste0(stem, ".png"))
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  data_path <- file.path(table_dir, paste0(stem, "_plot_data.tsv"))
  ggsave(png_path, plot, width = width, height = height, dpi = 300, limitsize = FALSE)
  ggsave(pdf_path, plot, width = width, height = height, device = grDevices::cairo_pdf, limitsize = FALSE)
  fwrite(as.data.table(plot_data), data_path, sep = "\t", na = "")
  list(png = png_path, pdf = pdf_path, data = data_path)
}

get_assay_matrix <- function(object, layer_or_slot) {
  tryCatch(
    GetAssayData(object = object, assay = "RNA", layer = layer_or_slot),
    error = function(err) GetAssayData(object = object, assay = "RNA", slot = layer_or_slot)
  )
}

logcpm_from_counts <- function(count_mat) {
  library_sizes <- Matrix::colSums(count_mat)
  library_sizes[library_sizes == 0] <- NA_real_
  cpm <- t(t(count_mat) / library_sizes) * 1e6
  cpm[is.na(cpm)] <- 0
  as.matrix(log1p(cpm))
}

align_model_features <- function(log_mat, features, feature_alias_table) {
  aligned <- matrix(
    0,
    nrow = length(features),
    ncol = ncol(log_mat),
    dimnames = list(features, colnames(log_mat))
  )
  aliases <- feature_alias_table %>%
    filter(model_feature %in% features) %>%
    arrange(match(model_feature, features))
  present <- !is.na(aliases$gse235021_gene_symbol) &
    aliases$gse235021_gene_symbol %in% rownames(log_mat)
  if (any(present)) {
    aligned[aliases$model_feature[present], ] <-
      log_mat[aliases$gse235021_gene_symbol[present], , drop = FALSE]
  }
  list(
    matrix = aligned,
    shared = aliases$model_feature[present],
    aliases = aliases
  )
}

build_model_feature_alias_table <- function(models, gse_feature_path) {
  gse_features <- fread(
    cmd = paste("gzip -cd", shQuote(gse_feature_path)),
    header = FALSE,
    col.names = c("ensembl_gene_id", "gse235021_gene_symbol", "feature_type")
  ) %>%
    filter(feature_type == "Gene Expression") %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)

  old_symbol_map_path <- file.path(
    data_root, "bulkseq_dataset", "GSE329024", "input",
    "GSE329024_gene_annotation.csv"
  )
  old_symbol_map <- fread(old_symbol_map_path) %>%
    transmute(
      reference_symbol = as.character(gene_id_or_symbol),
      ensembl_gene_id = as.character(gene_id)
    )
  prjna_feature_paths <- Sys.glob(file.path(
    data_root, "scRNA_dataset", "PRJNA863443",
    "cellranger_counts", "*", "outs", "filtered_feature_bc_matrix", "features.tsv.gz"
  ))
  if (length(prjna_feature_paths) == 0L) {
    stop("No PRJNA863443 feature map is available for frozen-model gene aliases.", call. = FALSE)
  }
  current_symbol_map <- fread(
    cmd = paste("gzip -cd", shQuote(prjna_feature_paths[[1]])),
    header = FALSE,
    col.names = c("ensembl_gene_id", "reference_symbol", "feature_type")
  ) %>%
    filter(feature_type == "Gene Expression") %>%
    select(reference_symbol, ensembl_gene_id)
  reference_map <- bind_rows(old_symbol_map, current_symbol_map) %>%
    filter(
      grepl("^ENSMUSG", ensembl_gene_id),
      !is.na(reference_symbol),
      reference_symbol != ""
    ) %>%
    distinct(reference_symbol, .keep_all = TRUE)
  symbol_to_ensembl <- setNames(reference_map$ensembl_gene_id, reference_map$reference_symbol)
  ensembl_to_gse_symbol <- setNames(gse_features$gse235021_gene_symbol, gse_features$ensembl_gene_id)
  gse_symbols <- unique(gse_features$gse235021_gene_symbol)

  model_features <- unique(unlist(lapply(models, `[[`, "features"), use.names = FALSE))
  bind_rows(lapply(model_features, function(model_feature) {
    syntactic_alias <- gsub("\\.", "-", model_feature)
    syntactic_alias <- sub("^X(?=[0-9])", "", syntactic_alias, perl = TRUE)
    stable_id <- if (grepl("^ENSMUSG", model_feature)) {
      model_feature
    } else {
      unname(symbol_to_ensembl[model_feature])
    }
    if (is.na(stable_id)) {
      stable_id <- unname(symbol_to_ensembl[syntactic_alias])
    }
    stable_symbol <- if (!is.na(stable_id)) {
      unname(ensembl_to_gse_symbol[stable_id])
    } else {
      NA_character_
    }
    if (model_feature %in% gse_symbols) {
      resolved_symbol <- model_feature
      resolution_method <- "direct_symbol"
    } else if (syntactic_alias %in% gse_symbols) {
      resolved_symbol <- syntactic_alias
      resolution_method <- "r_syntactic_name_reversal"
    } else if (!is.null(stable_symbol) && !is.na(stable_symbol) && stable_symbol %in% gse_symbols) {
      resolved_symbol <- stable_symbol
      resolution_method <- "stable_ensembl_alias"
    } else {
      resolved_symbol <- NA_character_
      resolution_method <- "unresolved_zero_filled"
    }
    tibble(
      model_feature = model_feature,
      syntactic_alias = syntactic_alias,
      ensembl_gene_id = stable_id,
      gse235021_gene_symbol = resolved_symbol,
      resolution_method = resolution_method,
      present_in_GSE235021 = !is.na(resolved_symbol)
    )
  }))
}

label_clusters <- function(object, sample_id) {
  expression <- get_assay_matrix(object, "data")
  clusters <- sort(unique(as.character(Idents(object))))
  score_rows <- list()
  label_rows <- list()
  for (cluster_id in clusters) {
    cells <- colnames(object)[as.character(Idents(object)) == cluster_id]
    scores <- vapply(marker_sets, function(markers) {
      features <- intersect(markers, rownames(expression))
      if (length(features) == 0L || length(cells) == 0L) return(0)
      mean(Matrix::rowMeans(expression[features, cells, drop = FALSE]))
    }, numeric(1))
    present <- vapply(marker_sets, function(markers) {
      length(intersect(markers, rownames(expression)))
    }, integer(1))
    score_rows[[cluster_id]] <- tibble(
      sample_id = sample_id,
      seurat_cluster = cluster_id,
      candidate_cell_type = names(scores),
      marker_score = as.numeric(scores),
      n_markers_present = as.integer(present)
    )
    order_index <- order(-scores, -present, names(scores))
    top_index <- order_index[[1]]
    second_index <- order_index[[2]]
    assigned <- if (present[[top_index]] >= 2L && scores[[top_index]] > 0.05) {
      names(scores)[[top_index]]
    } else {
      "Unknown"
    }
    label_rows[[cluster_id]] <- tibble(
      sample_id = sample_id,
      seurat_cluster = cluster_id,
      harmonized_cell_class = assigned,
      label_assignment_confidence = scores[[top_index]] - scores[[second_index]],
      top_marker_score = scores[[top_index]],
      second_marker_score = scores[[second_index]],
      n_markers_present = present[[top_index]],
      annotation_method = "per_library_seurat_cluster_canonical_marker_signature_scores"
    )
  }
  list(scores = bind_rows(score_rows), labels = bind_rows(label_rows))
}

sample_metadata <- fread(sample_metadata_path, na.strings = c("", "NA")) %>%
  as_tibble() %>%
  mutate(
    sample_id = as.character(sample_id),
    sample_title = as.character(sample_title),
    diet = as.character(diet),
    age_month = as.numeric(age_month),
    age_week = as.numeric(age_week),
    projection_group = if_else(diet == "polished_rice", "matched_control", "diet_intervention")
  ) %>%
  arrange(age_month, diet, replicate, sample_id)

if (nrow(sample_metadata) != 17L || n_distinct(sample_metadata$sample_id) != 17L) {
  stop("Expected 17 unique GSE235021 scRNA libraries.", call. = FALSE)
}
if (!all(file.exists(sample_metadata$seurat_rds))) {
  stop("Missing Seurat object(s): ", paste(sample_metadata$seurat_rds[!file.exists(sample_metadata$seurat_rds)], collapse = "; "), call. = FALSE)
}
if (!all(c(6, 9) %in% sample_metadata$age_month)) {
  stop("Both 6- and 9-month strata are required.", call. = FALSE)
}

models <- readRDS(model_path)
model_classes <- names(models)
if (!identical(sort(model_classes), sort(names(class_palette)))) {
  stop("Frozen model classes do not match the expected eight-class vocabulary.", call. = FALSE)
}
model_md5_before <- unname(tools::md5sum(model_path))
cv_metrics <- fread(cv_path) %>% as_tibble()

cluster_scores <- list()
cluster_labels <- list()
cell_annotations <- list()
metacell_metadata <- list()
metacell_counts <- list()
qc_rows <- list()

preprocess_cache_paths <- c(
  file.path(table_dir, "GSE235021_library_qc.tsv"),
  file.path(table_dir, "GSE235021_cluster_marker_signature_scores.tsv"),
  file.path(table_dir, "GSE235021_cluster_cell_type_labels.tsv"),
  file.path(table_dir, "GSE235021_harmonized_cell_metadata.tsv.gz"),
  file.path(table_dir, "GSE235021_metacell_metadata.tsv"),
  file.path(table_dir, "GSE235021_metacell_raw_counts.rds")
)
reuse_preprocessing <- all(file.exists(preprocess_cache_paths))
if (reuse_preprocessing) {
  message("Reusing validated GSE235021 annotation/metacell preprocessing cache")
  cluster_score_table <- fread(preprocess_cache_paths[[2]]) %>% as_tibble()
  cluster_label_table <- fread(preprocess_cache_paths[[3]]) %>% as_tibble()
  cell_annotation_table <- fread(
    cmd = paste("gzip -cd", shQuote(preprocess_cache_paths[[4]]))
  ) %>% as_tibble()
  metacell_meta <- fread(preprocess_cache_paths[[5]]) %>% as_tibble()
  metacell_count_matrix <- readRDS(preprocess_cache_paths[[6]])
  qc_rows <- list(cache = fread(preprocess_cache_paths[[1]]) %>% as_tibble())
  stopifnot(identical(metacell_meta$metacell_id, colnames(metacell_count_matrix)))
} else {
  message("Processing 17 libraries independently with the established external-cohort annotation policy")
  for (index in seq_len(nrow(sample_metadata))) {
  sample_row <- sample_metadata[index, ]
  sample_id <- sample_row$sample_id[[1]]
  message(sprintf("[%02d/17] %s (%s)", index, sample_id, sample_row$sample_title[[1]]))
  object <- readRDS(sample_row$seurat_rds[[1]])
  if ("Assay5" %in% class(object[["RNA"]])) {
    object <- JoinLayers(object, assay = "RNA")
  }
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-|^MT-")
  cells_before_qc <- ncol(object)
  object <- subset(
    object,
    subset = nFeature_RNA >= minimum_features & percent.mt <= maximum_mitochondrial_percent
  )
  cells_after_qc <- ncol(object)
  if (cells_after_qc < minimum_cells_per_metacell) {
    stop("Too few cells after QC for ", sample_id, call. = FALSE)
  }

  object <- NormalizeData(object, verbose = FALSE)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  object <- ScaleData(object, features = VariableFeatures(object), verbose = FALSE)
  object <- RunPCA(object, features = VariableFeatures(object), npcs = 30, verbose = FALSE)
  object <- FindNeighbors(object, dims = 1:20, verbose = FALSE)
  object <- FindClusters(object, resolution = 0.5, verbose = FALSE)

  annotation <- label_clusters(object, sample_id)
  cluster_scores[[sample_id]] <- annotation$scores
  cluster_labels[[sample_id]] <- annotation$labels
  label_map <- setNames(annotation$labels$harmonized_cell_class, annotation$labels$seurat_cluster)
  confidence_map <- setNames(annotation$labels$label_assignment_confidence, annotation$labels$seurat_cluster)
  cell_cluster <- as.character(Idents(object))
  cell_table <- as.data.table(object@meta.data, keep.rownames = "cell_id") %>%
    as_tibble() %>%
    transmute(
      cell_id,
      sample_id = sample_id,
      sample_title = sample_row$sample_title[[1]],
      age_month = sample_row$age_month[[1]],
      age_week = sample_row$age_week[[1]],
      diet = sample_row$diet[[1]],
      replicate = sample_row$replicate[[1]],
      seurat_cluster = cell_cluster,
      harmonized_cell_class = unname(label_map[cell_cluster]),
      label_assignment_confidence = unname(confidence_map[cell_cluster]),
      annotation_method = "per_library_seurat_cluster_canonical_marker_signature_scores",
      passes_confidence_policy = harmonized_cell_class %in% model_classes &
        is.finite(label_assignment_confidence) &
        label_assignment_confidence >= confidence_cutoff,
      metacell_id = NA_character_
    )

  retained <- cell_table %>%
    filter(passes_confidence_policy) %>%
    arrange(harmonized_cell_class, cell_id)
  metacell_rows <- list()
  member_rows <- list()
  for (cell_class in model_classes) {
    class_cells <- retained$cell_id[retained$harmonized_cell_class == cell_class]
    if (length(class_cells) < minimum_cells_per_metacell) next
    chunk <- ceiling(seq_along(class_cells) / target_cells_per_metacell)
    valid_chunks <- as.integer(names(table(chunk))[table(chunk) >= minimum_cells_per_metacell])
    for (chunk_id in valid_chunks) {
      members <- class_cells[chunk == chunk_id]
      metacell_id <- sprintf(
        "%s__%s__%s__mc%04d",
        dataset_folder, sample_id, cell_class, chunk_id
      )
      member_rows[[metacell_id]] <- tibble(
        metacell_id = metacell_id,
        cell_id = members
      )
      metacell_rows[[metacell_id]] <- tibble(
        metacell_id = metacell_id,
        sample_key = paste(dataset_folder, sample_id, sep = "__"),
        sample_id = sample_id,
        sample_title = sample_row$sample_title[[1]],
        dataset_folder = dataset_folder,
        dataset_id = dataset_id,
        dataset_name = dataset_name,
        model_role = "projection_only",
        projection_role = "diet_stress_test",
        projection_group = sample_row$projection_group[[1]],
        diet = sample_row$diet[[1]],
        replicate = sample_row$replicate[[1]],
        age_month = sample_row$age_month[[1]],
        age_weeks = sample_row$age_week[[1]],
        reproductive_age_bin = sample_row$age_bin[[1]],
        harmonized_cell_class = cell_class,
        label_assignment_source = "per_library_seurat_cluster_canonical_marker_signature_scores",
        label_assignment_confidence = mean(retained$label_assignment_confidence[retained$cell_id %in% members]),
        min_label_assignment_confidence = min(retained$label_assignment_confidence[retained$cell_id %in% members]),
        curation_status_summary = paste0("probable=", length(members)),
        n_confirmed_cells = 0L,
        n_probable_cells = length(members),
        n_cells = length(members)
      )
      cell_table$metacell_id[cell_table$cell_id %in% members] <- metacell_id
    }
  }

  sample_metacell_metadata <- bind_rows(metacell_rows)
  sample_members <- bind_rows(member_rows)
  if (nrow(sample_metacell_metadata) == 0L) {
    stop("No eligible metacells produced for ", sample_id, call. = FALSE)
  }
  raw_counts <- get_assay_matrix(object, "counts")
  group_factor <- factor(sample_members$metacell_id, levels = sample_metacell_metadata$metacell_id)
  design <- sparse.model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)
  sample_metacell_counts <- raw_counts[, sample_members$cell_id, drop = FALSE] %*% design
  sample_metacell_counts <- as(sample_metacell_counts, "dgCMatrix")

  qc_rows[[sample_id]] <- tibble(
    sample_id = sample_id,
    sample_title = sample_row$sample_title[[1]],
    age_month = sample_row$age_month[[1]],
    diet = sample_row$diet[[1]],
    cells_before_qc = cells_before_qc,
    cells_after_qc = cells_after_qc,
    cells_passing_annotation_confidence = sum(cell_table$passes_confidence_policy),
    cells_in_metacells = sum(!is.na(cell_table$metacell_id)),
    n_metacells = nrow(sample_metacell_metadata),
    n_projected_cell_classes = n_distinct(sample_metacell_metadata$harmonized_cell_class)
  )
  cell_annotations[[sample_id]] <- cell_table
  metacell_metadata[[sample_id]] <- sample_metacell_metadata
  metacell_counts[[sample_id]] <- sample_metacell_counts
  rm(object, raw_counts, sample_metacell_counts)
  gc()
  }

  cluster_score_table <- bind_rows(cluster_scores)
  cluster_label_table <- bind_rows(cluster_labels)
  cell_annotation_table <- bind_rows(cell_annotations)
  metacell_meta <- bind_rows(metacell_metadata)
  if (!all(vapply(metacell_counts, function(x) identical(rownames(x), rownames(metacell_counts[[1]])), logical(1)))) {
    stop("GSE235021 Seurat objects do not share an identical gene namespace.", call. = FALSE)
  }
  metacell_count_matrix <- do.call(cbind, metacell_counts)
  stopifnot(identical(metacell_meta$metacell_id, colnames(metacell_count_matrix)))
}

message("Projecting metacells through the eight frozen final models")
gse_feature_path <- file.path(
  sample_metadata$count_matrix_dir[[1]], "features.tsv.gz"
)
feature_alias_table <- build_model_feature_alias_table(models, gse_feature_path)
prediction_rows <- list()
coverage_rows <- list()
for (cell_class in model_classes) {
  class_meta <- metacell_meta %>% filter(harmonized_cell_class == cell_class)
  if (nrow(class_meta) == 0L) next
  class_counts <- metacell_count_matrix[, class_meta$metacell_id, drop = FALSE]
  class_logcpm <- logcpm_from_counts(class_counts)
  feature_alignment <- align_model_features(
    class_logcpm, models[[cell_class]]$features, feature_alias_table
  )
  predictions <- as.numeric(predict(models[[cell_class]]$fit, newx = t(feature_alignment$matrix)))
  prediction_rows[[cell_class]] <- class_meta %>%
    mutate(
      split_id = "projection_only",
      fold_id = "final_all_metacells",
      clock_type = "celltype_metacell_elasticnet_clock",
      model_evaluation_role = "external_diet_projection_stress_test",
      predicted_age_weeks_raw = predictions,
      predicted_age_weeks = predictions,
      chronological_delta_age_weeks = predicted_age_weeks - age_weeks,
      n_model_features = length(models[[cell_class]]$features),
      n_model_features_present = length(feature_alignment$shared),
      model_feature_coverage = n_model_features_present / n_model_features,
      best_alpha = models[[cell_class]]$best_alpha,
      best_lambda = models[[cell_class]]$best_lambda,
      inner_cv_mae_weeks = models[[cell_class]]$cv_mae_weeks
    )
  coverage_rows[[cell_class]] <- tibble(
    harmonized_cell_class = cell_class,
    n_model_features = length(models[[cell_class]]$features),
    n_model_features_present = length(feature_alignment$shared),
    model_feature_coverage = length(feature_alignment$shared) / length(models[[cell_class]]$features)
  )
}
predictions <- bind_rows(prediction_rows)
feature_coverage <- bind_rows(coverage_rows)

sample_celltype_summary <- predictions %>%
  group_by(
    sample_key, sample_id, sample_title, age_month, age_weeks, diet, replicate,
    projection_group, harmonized_cell_class
  ) %>%
  summarise(
    n_projected_metacells = n(),
    n_cells = sum(n_cells),
    median_predicted_age_weeks = median(predicted_age_weeks),
    mean_predicted_age_weeks = mean(predicted_age_weeks),
    sd_predicted_age_weeks = sd(predicted_age_weeks),
    median_chronological_delta_age_weeks = median(chronological_delta_age_weeks),
    .groups = "drop"
  )

control_baseline <- predictions %>%
  filter(diet == "polished_rice") %>%
  group_by(age_month, harmonized_cell_class) %>%
  summarise(
    matched_control_mean_predicted_age_weeks = mean(predicted_age_weeks),
    matched_control_median_predicted_age_weeks = median(predicted_age_weeks),
    n_control_samples = n_distinct(sample_id),
    n_control_metacells = n(),
    .groups = "drop"
  )
predictions_delta <- predictions %>%
  left_join(control_baseline, by = c("age_month", "harmonized_cell_class")) %>%
  mutate(
    delta_predicted_age_weeks = predicted_age_weeks - matched_control_mean_predicted_age_weeks
  )
sample_delta <- predictions_delta %>%
  group_by(
    sample_key, sample_id, sample_title, age_month, age_weeks, diet, replicate,
    projection_group, harmonized_cell_class, n_control_samples, n_control_metacells
  ) %>%
  summarise(
    n_projected_metacells = n(),
    median_predicted_age_weeks = median(predicted_age_weeks),
    median_delta_predicted_age_weeks = median(delta_predicted_age_weeks),
    mean_delta_predicted_age_weeks = mean(delta_predicted_age_weeks),
    .groups = "drop"
  )

contrast_definition <- tribble(
  ~contrast_id, ~age_month, ~case_diet, ~control_diet, ~contrast_label,
  "GSE235021_black_rice_6m_vs_polished", 6, "black_rice", "polished_rice", "6 mo: black vs polished",
  "GSE235021_white_brown_rice_6m_vs_polished", 6, "white_brown_rice", "polished_rice", "6 mo: white/brown vs polished",
  "GSE235021_black_rice_9m_vs_polished", 9, "black_rice", "polished_rice", "9 mo: black vs polished",
  "GSE235021_white_brown_rice_9m_vs_polished", 9, "white_brown_rice", "polished_rice", "9 mo: white/brown vs polished"
)
contrast_metacells <- bind_rows(lapply(seq_len(nrow(contrast_definition)), function(index) {
  definition <- contrast_definition[index, ]
  predictions_delta %>%
    filter(
      age_month == definition$age_month[[1]],
      diet %in% c(definition$case_diet[[1]], definition$control_diet[[1]])
    ) %>%
    mutate(
      contrast_id = definition$contrast_id[[1]],
      contrast_label = definition$contrast_label[[1]],
      contrast_group = if_else(diet == definition$case_diet[[1]], "case", "matched_control")
    ) %>%
    filter(is.finite(delta_predicted_age_weeks))
}))
contrast_samples <- bind_rows(lapply(seq_len(nrow(contrast_definition)), function(index) {
  definition <- contrast_definition[index, ]
  sample_delta %>%
    filter(
      age_month == definition$age_month[[1]],
      diet %in% c(definition$case_diet[[1]], definition$control_diet[[1]])
    ) %>%
    mutate(
      contrast_id = definition$contrast_id[[1]],
      contrast_label = definition$contrast_label[[1]],
      contrast_group = if_else(diet == definition$case_diet[[1]], "case", "matched_control")
    ) %>%
    filter(is.finite(median_delta_predicted_age_weeks))
}))

contrast_summary <- contrast_samples %>%
  group_by(contrast_id, contrast_label, contrast_group, diet, harmonized_cell_class) %>%
  summarise(
    n_samples = n_distinct(sample_id),
    n_metacells = sum(n_projected_metacells),
    median_sample_delta_weeks = median(median_delta_predicted_age_weeks),
    mean_sample_delta_weeks = mean(median_delta_predicted_age_weeks),
    sd_sample_delta_weeks = sd(median_delta_predicted_age_weeks),
    .groups = "drop"
  )
contrast_effect <- contrast_summary %>%
  select(contrast_id, contrast_label, contrast_group, harmonized_cell_class, median_sample_delta_weeks) %>%
  pivot_wider(names_from = contrast_group, values_from = median_sample_delta_weeks) %>%
  mutate(case_minus_control_median_delta_weeks = case - matched_control)

sample_level_tests <- contrast_samples %>%
  group_by(contrast_id, contrast_label, harmonized_cell_class) %>%
  group_modify(~{
    case_values <- .x$median_delta_predicted_age_weeks[.x$contrast_group == "case"]
    control_values <- .x$median_delta_predicted_age_weeks[.x$contrast_group == "matched_control"]
    test <- if (length(case_values) >= 2L && length(control_values) >= 2L) {
      tryCatch(t.test(case_values, control_values, var.equal = FALSE), error = function(e) NULL)
    } else NULL
    tibble(
      n_case_samples = length(case_values),
      n_control_samples = length(control_values),
      mean_case_delta_weeks = mean(case_values),
      mean_control_delta_weeks = mean(control_values),
      case_minus_control_mean_delta_weeks = mean(case_values) - mean(control_values),
      welch_t_p_value = if (is.null(test)) NA_real_ else test$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(welch_t_p_value_bh = p.adjust(welch_t_p_value, method = "BH"))

six_clock_classes <- setdiff(model_classes, c("B_Plasma", "T_NK"))
metacell_level_tests_six_clocks <- contrast_metacells %>%
  filter(harmonized_cell_class %in% six_clock_classes) %>%
  group_by(contrast_id, contrast_label, harmonized_cell_class) %>%
  group_modify(~{
    case_values <- .x$delta_predicted_age_weeks[.x$contrast_group == "case"]
    control_values <- .x$delta_predicted_age_weeks[.x$contrast_group == "matched_control"]
    test <- if (length(case_values) >= 2L && length(control_values) >= 2L) {
      tryCatch(t.test(case_values, control_values, var.equal = FALSE), error = function(e) NULL)
    } else NULL
    tibble(
      n_case_metacells = length(case_values),
      n_control_metacells = length(control_values),
      mean_case_delta_weeks = mean(case_values),
      mean_control_delta_weeks = mean(control_values),
      case_minus_control_mean_delta_weeks =
        mean(case_values) - mean(control_values),
      welch_t_p_value = if (is.null(test)) NA_real_ else test$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    welch_t_p_value_bh = p.adjust(welch_t_p_value, method = "BH"),
    inference_unit = "projected_metacell",
    interpretation_scope = paste(
      "descriptive_metacell_level_test;",
      "not_biological_library_level_inference"
    )
  )

message("Writing projection tables")
fwrite(as.data.table(sample_metadata), file.path(table_dir, "GSE235021_sample_metadata.tsv"), sep = "\t", na = "")
fwrite(as.data.table(bind_rows(qc_rows)), file.path(table_dir, "GSE235021_library_qc.tsv"), sep = "\t", na = "")
fwrite(as.data.table(cluster_score_table), file.path(table_dir, "GSE235021_cluster_marker_signature_scores.tsv"), sep = "\t", na = "")
fwrite(as.data.table(cluster_label_table), file.path(table_dir, "GSE235021_cluster_cell_type_labels.tsv"), sep = "\t", na = "")
fwrite(as.data.table(cell_annotation_table), file.path(table_dir, "GSE235021_harmonized_cell_metadata.tsv.gz"), sep = "\t", na = "", compress = "gzip")
fwrite(as.data.table(metacell_meta), file.path(table_dir, "GSE235021_metacell_metadata.tsv"), sep = "\t", na = "")
saveRDS(metacell_count_matrix, file.path(table_dir, "GSE235021_metacell_raw_counts.rds"))
fwrite(as.data.table(feature_alias_table), file.path(table_dir, "GSE235021_model_feature_alias_map.tsv"), sep = "\t", na = "")
fwrite(as.data.table(feature_coverage), file.path(table_dir, "GSE235021_model_feature_coverage.tsv"), sep = "\t", na = "")
fwrite(as.data.table(predictions), file.path(table_dir, "GSE235021_metacell_clock_predictions.tsv"), sep = "\t", na = "")
fwrite(as.data.table(sample_celltype_summary), file.path(table_dir, "GSE235021_sample_celltype_prediction_summary.tsv"), sep = "\t", na = "")
fwrite(as.data.table(control_baseline), file.path(table_dir, "GSE235021_matched_control_baseline.tsv"), sep = "\t", na = "")
fwrite(as.data.table(predictions_delta), file.path(table_dir, "GSE235021_metacell_delta_to_matched_controls.tsv"), sep = "\t", na = "")
fwrite(as.data.table(sample_delta), file.path(table_dir, "GSE235021_sample_celltype_delta_to_matched_controls.tsv"), sep = "\t", na = "")
fwrite(as.data.table(contrast_definition), file.path(table_dir, "GSE235021_contrast_definition.tsv"), sep = "\t", na = "")
fwrite(as.data.table(contrast_metacells), file.path(table_dir, "GSE235021_contrast_metacell_plot_source.tsv"), sep = "\t", na = "")
fwrite(as.data.table(contrast_samples), file.path(table_dir, "GSE235021_contrast_sample_plot_source.tsv"), sep = "\t", na = "")
fwrite(as.data.table(contrast_summary), file.path(table_dir, "GSE235021_contrast_summary_by_celltype.tsv"), sep = "\t", na = "")
fwrite(as.data.table(contrast_effect), file.path(table_dir, "GSE235021_contrast_effect_by_celltype.tsv"), sep = "\t", na = "")
fwrite(as.data.table(sample_level_tests), file.path(table_dir, "GSE235021_sample_level_welch_t_tests.tsv"), sep = "\t", na = "")
fwrite(
  as.data.table(metacell_level_tests_six_clocks),
  file.path(table_dir, "GSE235021_metacell_level_welch_t_tests_six_clocks.tsv"),
  sep = "\t",
  na = ""
)

message("Rendering R figure release")
sample_order <- sample_metadata %>%
  arrange(age_month, factor(diet, levels = c("polished_rice", "black_rice", "white_brown_rice")), replicate) %>%
  pull(sample_title)
availability <- sample_celltype_summary %>%
  complete(
    sample_title = sample_order,
    harmonized_cell_class = model_classes,
    fill = list(n_projected_metacells = 0L, n_cells = 0L)
  ) %>%
  left_join(sample_metadata %>% select(sample_title, age_month, diet), by = "sample_title") %>%
  mutate(
    sample_title = factor(sample_title, levels = rev(sample_order)),
    harmonized_cell_class = factor(harmonized_cell_class, levels = model_classes)
  )
p_availability <- ggplot(availability, aes(harmonized_cell_class, sample_title, fill = n_projected_metacells)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = if_else(n_projected_metacells > 0, as.character(n_projected_metacells), "")), size = 2) +
  scale_fill_gradient(low = "grey96", high = "#3B4CC0") +
  labs(
    title = "GSE235021 metacells available to frozen cell-type clocks",
    subtitle = "All 17 whole-ovary droplet libraries; tiles report projected 25–50-cell metacells",
    x = "Frozen metacell clock", y = "Library", fill = "Metacells"
  ) +
  theme_projection(8) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
plot_availability <- save_projection_plot(
  p_availability, "GSE235021_metacell_availability_by_library", 10.5, 7.3, availability
)

availability_summary <- availability %>%
  mutate(
    sample_title = as.character(sample_title),
    harmonized_cell_class = as.character(harmonized_cell_class),
    diet = factor(coalesce(diet.y, diet.x), levels = names(diet_labels)),
    metacell_available = n_projected_metacells > 0
  ) %>%
  group_by(diet, harmonized_cell_class) %>%
  summarise(
    n_projected_metacells = sum(n_projected_metacells, na.rm = TRUE),
    n_projected_cells = sum(n_cells, na.rm = TRUE),
    n_libraries = n_distinct(sample_title),
    n_libraries_with_metacells = n_distinct(sample_title[metacell_available]),
    .groups = "drop"
  ) %>%
  group_by(harmonized_cell_class) %>%
  mutate(
    total_projected_metacells = sum(n_projected_metacells, na.rm = TRUE),
    n_libraries_with_metacells_total = sum(n_libraries_with_metacells),
    n_libraries_total = sum(n_libraries),
    coverage_label = sprintf(
      "%s/%s libraries",
      n_libraries_with_metacells_total,
      n_libraries_total
    )
  ) %>%
  ungroup() %>%
  mutate(
    harmonized_cell_class = factor(
      harmonized_cell_class,
      levels = unique(
        harmonized_cell_class[order(total_projected_metacells)]
      )
    )
  ) %>%
  mutate(
    label_x = total_projected_metacells +
      max(total_projected_metacells, na.rm = TRUE) * 0.02,
    diet = factor(diet, levels = names(diet_labels), labels = unname(diet_labels))
  ) %>%
  arrange(harmonized_cell_class, diet)

availability_summary_labels <- availability_summary %>%
  distinct(
    harmonized_cell_class,
    total_projected_metacells,
    n_libraries_with_metacells_total,
    n_libraries_total,
    coverage_label,
    label_x
  )
availability_summary_xlim <-
  max(availability_summary$total_projected_metacells, na.rm = TRUE) * 1.36

p_availability_summary <- ggplot(
  availability_summary,
  aes(n_projected_metacells, harmonized_cell_class, fill = diet)
) +
  geom_col(width = 0.72) +
  geom_text(
    data = availability_summary_labels,
    aes(x = label_x, y = harmonized_cell_class, label = coverage_label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.65,
    colour = "grey25"
  ) +
  scale_fill_manual(
    values = setNames(diet_palette, unname(diet_labels)),
    breaks = unname(diet_labels),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, availability_summary_xlim),
    expand = expansion(mult = c(0, 0))
  ) +
  theme_projection() +
  labs(
    x = "Projected metacells across GSE235021 libraries",
    y = "Cell type",
    fill = "Diet",
    title = "GSE235021 metacell availability summary",
    subtitle = paste(
      "Summed projected metacells by rice-diet group;",
      "labels show library-cell-type coverage."
    )
  )
plot_availability_summary <- save_projection_plot(
  p_availability_summary,
  "GSE235021_metacell_availability_summary",
  7.4,
  4.9,
  availability_summary
)

sample_prediction_plot_data <- sample_celltype_summary %>%
  mutate(
    diet = factor(diet, levels = names(diet_labels), labels = unname(diet_labels)),
    age_label = paste0(age_month, " months"),
    harmonized_cell_class = factor(harmonized_cell_class, levels = model_classes)
  )
p_prediction <- ggplot(
  sample_prediction_plot_data,
  aes(diet, median_predicted_age_weeks, colour = diet)
) +
  geom_hline(
    data = distinct(sample_prediction_plot_data, age_label, age_weeks),
    aes(yintercept = age_weeks),
    colour = "grey55", linetype = 2, linewidth = 0.35
  ) +
  geom_point(position = position_jitter(width = 0.10, height = 0), size = 1.4, alpha = 0.85) +
  facet_grid(harmonized_cell_class ~ age_label, scales = "free_y") +
  scale_colour_manual(values = setNames(diet_palette, unname(diet_labels)), guide = "none") +
  labs(
    title = "GSE235021 sample-level median metacell-clock predictions",
    subtitle = "Each point is one biological library × cell-type clock; dashed line is chronological age",
    x = NULL, y = "Median predicted age (weeks)"
  ) +
  theme_projection(7.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
plot_prediction <- save_projection_plot(
  p_prediction, "GSE235021_sample_celltype_predicted_age", 10.8, 12.5, sample_prediction_plot_data
)

delta_plot_data <- contrast_metacells %>%
  mutate(
    contrast_label = factor(contrast_label, levels = contrast_definition$contrast_label),
    contrast_group = factor(contrast_group, levels = c("matched_control", "case")),
    harmonized_cell_class = factor(harmonized_cell_class, levels = model_classes)
  )
p_delta <- ggplot(
  delta_plot_data,
  aes(contrast_label, delta_predicted_age_weeks, fill = contrast_group, colour = contrast_group)
) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = 2, linewidth = 0.35) +
  geom_boxplot(
    position = position_dodge(width = 0.70), width = 0.58,
    outlier.shape = NA, alpha = 0.65, linewidth = 0.35
  ) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.10, dodge.width = 0.70),
    size = 0.35, alpha = 0.35
  ) +
  facet_wrap(~harmonized_cell_class, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(matched_control = "#BDBDBD", case = "#3B4CC0")) +
  scale_colour_manual(values = c(matched_control = "#636363", case = "#253494")) +
  labs(
    title = "GSE235021 metacell delta-age distributions within matched diet contrasts",
    subtitle = "Each point is one projected metacell; controls are pooled only within age and cell-type clock",
    x = NULL, y = "Predicted-age delta from matched-control mean (weeks)",
    fill = "Projection group", colour = "Projection group"
  ) +
  theme_projection(7.5) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
plot_delta <- save_projection_plot(
  p_delta, "GSE235021_metacell_delta_distribution_by_celltype", 11.5, 12, delta_plot_data
)

six_clock_delta_data <- predictions_delta %>%
  filter(harmonized_cell_class %in% six_clock_classes) %>%
  mutate(
    age_label = factor(
      paste0(age_month, " months"),
      levels = c("6 months", "9 months")
    ),
    diet = factor(diet, levels = names(diet_labels), labels = unname(diet_labels)),
    harmonized_cell_class = factor(
      harmonized_cell_class,
      levels = six_clock_classes
    )
  )

six_clock_panel_ranges <- six_clock_delta_data %>%
  group_by(age_month, age_label, harmonized_cell_class) %>%
  summarise(
    panel_min = min(delta_predicted_age_weeks, na.rm = TRUE),
    panel_max = max(delta_predicted_age_weeks, na.rm = TRUE),
    panel_span = max(panel_max - panel_min, 1),
    .groups = "drop"
  )

format_test_value <- function(value) {
  case_when(
    is.na(value) ~ "NE",
    value < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", value)
  )
}

six_clock_test_annotations <- metacell_level_tests_six_clocks %>%
  left_join(
    contrast_definition %>% select(contrast_id, age_month, case_diet),
    by = "contrast_id"
  ) %>%
  mutate(
    age_label = factor(
      paste0(age_month, " months"),
      levels = c("6 months", "9 months")
    ),
    case_diet_label = factor(
      case_diet,
      levels = names(diet_labels),
      labels = unname(diet_labels)
    ),
    harmonized_cell_class = factor(
      harmonized_cell_class,
      levels = six_clock_classes
    ),
    control_x = 1,
    case_x = match(case_diet, names(diet_labels)),
    test_order = if_else(case_diet == "black_rice", 1L, 2L),
    test_label = paste0(
      "Welch p=", format_test_value(welch_t_p_value),
      "; BH q=", format_test_value(welch_t_p_value_bh)
    )
  ) %>%
  left_join(
    six_clock_panel_ranges,
    by = c("age_month", "age_label", "harmonized_cell_class")
  ) %>%
  mutate(
    bracket_y = panel_max + panel_span * if_else(test_order == 1L, 0.18, 0.48),
    bracket_tip_y = bracket_y - panel_span * 0.045,
    label_y = bracket_y + panel_span * 0.035,
    label_x = (control_x + case_x) / 2
  )

p_six_clock_delta_age_tests <- ggplot(
  six_clock_delta_data,
  aes(diet, delta_predicted_age_weeks, fill = diet)
) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = 2, linewidth = 0.35) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    alpha = 0.58,
    linewidth = 0.32,
    colour = "grey25"
  ) +
  geom_point(
    position = position_jitter(width = 0.08, height = 0, seed = 20260728),
    shape = 21,
    size = 0.48,
    stroke = 0.16,
    colour = "grey15",
    alpha = 0.48
  ) +
  geom_segment(
    data = six_clock_test_annotations,
    aes(
      x = control_x,
      xend = case_x,
      y = bracket_y,
      yend = bracket_y,
      colour = case_diet_label
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data = six_clock_test_annotations,
    aes(
      x = control_x,
      xend = control_x,
      y = bracket_tip_y,
      yend = bracket_y,
      colour = case_diet_label
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_segment(
    data = six_clock_test_annotations,
    aes(
      x = case_x,
      xend = case_x,
      y = bracket_tip_y,
      yend = bracket_y,
      colour = case_diet_label
    ),
    inherit.aes = FALSE,
    linewidth = 0.32
  ) +
  geom_text(
    data = six_clock_test_annotations,
    aes(
      x = label_x,
      y = label_y,
      label = test_label,
      colour = case_diet_label
    ),
    inherit.aes = FALSE,
    size = 1.85,
    vjust = 0,
    show.legend = FALSE
  ) +
  facet_grid(harmonized_cell_class ~ age_label, scales = "free_y") +
  scale_fill_manual(values = setNames(diet_palette, unname(diet_labels)), guide = "none") +
  scale_colour_manual(values = setNames(diet_palette, unname(diet_labels)), guide = "none") +
  scale_x_discrete(labels = c(
    "Polished rice" = "Polished",
    "Black rice" = "Black",
    "White/brown rice" = "White/\nbrown"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    title = "GSE235021 age-stratified diet effects across six metacell clocks",
    subtitle = paste(
      "Each point is one projected metacell;",
      "Welch tests are descriptive metacell-level matched-control comparisons."
    ),
    x = NULL,
    y = "Metacell predicted-age delta from matched-control mean (weeks)"
  ) +
  theme_projection(7.5) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1),
    panel.spacing.x = grid::unit(0.55, "lines"),
    panel.spacing.y = grid::unit(0.45, "lines"),
    strip.text.y = element_text(angle = 0)
  )

six_clock_delta_test_plot_source <- bind_rows(
  six_clock_delta_data %>%
    mutate(record_type = "projected_metacell"),
  six_clock_test_annotations %>%
    mutate(record_type = "metacell_level_welch_test_annotation")
)
plot_six_clock_delta_age_tests <- save_projection_plot(
  p_six_clock_delta_age_tests,
  "GSE235021_sample_delta_by_age_six_clocks_welch_tests",
  10.8,
  12.2,
  six_clock_delta_test_plot_source
)

six_clock_quality_labels <- cv_metrics %>%
  filter(harmonized_cell_class %in% six_clock_classes) %>%
  mutate(
    celltype_quality_label = sprintf(
      "%s\nMAE=%.1fw; R2=%.2f; slope=%.2f",
      harmonized_cell_class,
      mae_weeks,
      r2,
      slope
    )
  ) %>%
  select(harmonized_cell_class, celltype_quality_label)
six_clock_quality_label_levels <- six_clock_quality_labels %>%
  arrange(match(harmonized_cell_class, six_clock_classes)) %>%
  pull(celltype_quality_label)

six_clock_horizontal_data <- predictions_delta %>%
  filter(harmonized_cell_class %in% six_clock_classes) %>%
  left_join(six_clock_quality_labels, by = "harmonized_cell_class") %>%
  mutate(
    diet_code = diet,
    diet = factor(diet, levels = names(diet_labels), labels = unname(diet_labels)),
    celltype_quality_label = factor(
      celltype_quality_label,
      levels = six_clock_quality_label_levels
    )
  )

plot_six_clock_horizontal_age <- function(age_value) {
  age_data <- six_clock_horizontal_data %>%
    filter(age_month == age_value) %>%
    mutate(
      clock_label = factor(
        harmonized_cell_class,
        levels = six_clock_classes
      )
    )
  age_tests <- six_clock_test_annotations %>%
    filter(age_month == age_value) %>%
    transmute(
      harmonized_cell_class = as.character(harmonized_cell_class),
      diet_code = case_diet,
      welch_t_p_value,
      welch_t_p_value_bh
    )
  age_ranges <- age_data %>%
    group_by(harmonized_cell_class, clock_label) %>%
    summarise(
      panel_min = min(delta_predicted_age_weeks, na.rm = TRUE),
      panel_max = max(delta_predicted_age_weeks, na.rm = TRUE),
      panel_span = max(panel_max - panel_min, 1),
      .groups = "drop"
    )
  age_labels <- age_data %>%
    group_by(
      harmonized_cell_class,
      clock_label,
      diet_code,
      diet
    ) %>%
    summarise(
      n_metacells = n(),
      .groups = "drop"
    ) %>%
    left_join(age_ranges, by = c("harmonized_cell_class", "clock_label")) %>%
    mutate(
      n_label = paste0("n=", n_metacells),
      n_label_y = panel_max + panel_span * 0.07
    )
  age_brackets <- age_tests %>%
    left_join(age_ranges, by = "harmonized_cell_class") %>%
    mutate(
      clock_label = factor(
        harmonized_cell_class,
        levels = six_clock_classes
      ),
      clock_index = match(
        harmonized_cell_class,
        six_clock_classes
      ),
      case_diet_label = factor(
        diet_code,
        levels = names(diet_labels),
        labels = unname(diet_labels)
      ),
      control_x = clock_index - 0.24,
      case_x = clock_index + if_else(diet_code == "black_rice", 0, 0.24),
      bracket_y = panel_max +
        panel_span * if_else(diet_code == "black_rice", 0.24, 0.52),
      bracket_tick_y = bracket_y - panel_span * 0.045,
      bracket_label_x = (control_x + case_x) / 2,
      bracket_label_y = bracket_y + panel_span * 0.035,
      significance_label = case_when(
        welch_t_p_value_bh < 0.001 ~ "***",
        welch_t_p_value_bh < 0.01 ~ "**",
        welch_t_p_value_bh < 0.05 ~ "*",
        TRUE ~ paste0("p=", format_test_value(welch_t_p_value))
      )
    )

  age_plot <- ggplot(
    age_data,
    aes(
      clock_label,
      delta_predicted_age_weeks,
      fill = diet,
      group = interaction(clock_label, diet, drop = TRUE)
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      colour = "grey55",
      linewidth = 0.35
    ) +
    geom_boxplot(
      width = 0.44,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.30,
      alpha = 0.82,
      position = position_dodge2(width = 0.72, preserve = "single")
    ) +
    geom_point(
      aes(colour = diet),
      position = position_jitterdodge(
        jitter.width = 0.08,
        jitter.height = 0,
        dodge.width = 0.72,
        seed = 20260728 + age_value
      ),
      shape = 21,
      size = 0.48,
      alpha = 0.50,
      stroke = 0.16
    ) +
    geom_segment(
      data = age_brackets,
      aes(
        x = control_x,
        xend = case_x,
        y = bracket_y,
        yend = bracket_y,
        colour = case_diet_label
      ),
      inherit.aes = FALSE,
      linewidth = 0.34
    ) +
    geom_segment(
      data = age_brackets,
      aes(
        x = control_x,
        xend = control_x,
        y = bracket_tick_y,
        yend = bracket_y,
        colour = case_diet_label
      ),
      inherit.aes = FALSE,
      linewidth = 0.34
    ) +
    geom_segment(
      data = age_brackets,
      aes(
        x = case_x,
        xend = case_x,
        y = bracket_tick_y,
        yend = bracket_y,
        colour = case_diet_label
      ),
      inherit.aes = FALSE,
      linewidth = 0.34
    ) +
    geom_text(
      data = age_brackets,
      aes(
        x = bracket_label_x,
        y = bracket_label_y,
        label = significance_label,
        colour = case_diet_label
      ),
      inherit.aes = FALSE,
      size = 2.3,
      hjust = 0,
      fontface = "bold",
      show.legend = FALSE
    ) +
    geom_text(
      data = age_labels,
      aes(
        x = clock_label,
        y = n_label_y,
        label = n_label,
        colour = diet,
        group = diet
      ),
      inherit.aes = FALSE,
      position = position_dodge(width = 0.72),
      size = 1.85,
      vjust = 0,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = setNames(diet_palette, unname(diet_labels)),
      drop = FALSE
    ) +
    scale_colour_manual(
      values = setNames(diet_palette, unname(diet_labels)),
      drop = FALSE
    ) +
    scale_x_discrete(
      labels = c(
        Endothelial = "Endothelial",
        Epithelial = "Epithelial",
        Granulosa = "Granulosa",
        Immune_Myeloid = "Immune\nMyeloid",
        Pericyte_Smooth_Muscle = "Pericyte/\nSmooth muscle",
        Theca_Stroma = "Theca/\nStroma"
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.07, 0.10))
    ) +
    theme_projection(8.2) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(
      x = "Cell-type metacell clock",
      y = "Metacell predicted-age delta from matched-control mean (weeks)",
      fill = "Diet",
      colour = "Diet",
      title = paste0(
        "GSE235021 ", age_value,
        "-month diet effects across six metacell clocks"
      ),
      subtitle = paste(
        "Each point is one projected metacell;",
        "stars use BH q (*<0.05, **<0.01, ***<0.001);",
        "non-significant labels show Welch p."
      )
    )
  age_plot_source <- bind_rows(
    age_data %>% mutate(record_type = "projected_metacell"),
    age_labels %>%
      mutate(record_type = "metacell_level_n_annotation"),
    age_brackets %>%
      mutate(record_type = "metacell_level_significance_bracket")
  )
  stem <- paste0(
    "GSE235021_sample_delta_six_clocks_",
    age_value,
    "m_horizontal_welch_tests"
  )
  save_projection_plot(age_plot, stem, 9.2, 6.6, age_plot_source)
}

plot_six_clock_delta_6m_horizontal <- plot_six_clock_horizontal_age(6)
plot_six_clock_delta_9m_horizontal <- plot_six_clock_horizontal_age(9)

heatmap_data <- tidyr::expand_grid(
  contrast_id = contrast_definition$contrast_id,
  harmonized_cell_class = model_classes
) %>%
  left_join(
    contrast_definition %>% select(contrast_id, contrast_label),
    by = "contrast_id"
  ) %>%
  left_join(
    contrast_effect %>%
      select(
        contrast_id, harmonized_cell_class,
        case_minus_control_median_delta_weeks
      ),
    by = c("contrast_id", "harmonized_cell_class")
  ) %>%
  mutate(
    contrast_label = factor(contrast_label, levels = contrast_definition$contrast_label),
    harmonized_cell_class = factor(harmonized_cell_class, levels = rev(model_classes)),
    effect_estimable = is.finite(case_minus_control_median_delta_weeks),
    effect_label = if_else(
      effect_estimable,
      sprintf("%.1f", case_minus_control_median_delta_weeks),
      "NE"
    )
  )
heat_limit <- max(abs(heatmap_data$case_minus_control_median_delta_weeks), na.rm = TRUE)
p_heatmap <- ggplot(
  heatmap_data,
  aes(contrast_label, harmonized_cell_class, fill = case_minus_control_median_delta_weeks)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = effect_label), size = 2.25) +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "white", high = "#B40426",
    midpoint = 0, limits = c(-heat_limit, heat_limit), na.value = "grey85"
  ) +
  labs(
    title = "GSE235021 diet effects across frozen metacell clocks",
    subtitle = "Case minus matched-control median sample-level delta; NE indicates unavailable case or matched-control metacells",
    x = NULL, y = "Frozen metacell clock", fill = "Effect (weeks)"
  ) +
  theme_projection(8) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
plot_heatmap <- save_projection_plot(
  p_heatmap, "GSE235021_diet_effect_heatmap_by_celltype", 9.2, 5.8, heatmap_data
)

model_md5_after <- unname(tools::md5sum(model_path))
output_files <- unlist(
  list(
    plot_availability,
    plot_availability_summary,
    plot_prediction,
    plot_delta,
    plot_six_clock_delta_age_tests,
    plot_six_clock_delta_6m_horizontal,
    plot_six_clock_delta_9m_horizontal,
    plot_heatmap
  ),
  use.names = FALSE
)
graphic_output_files <- output_files[grepl("\\.(pdf|png)$", output_files)]
validation <- tribble(
  ~check, ~expected, ~observed, ~status,
  "unique_input_libraries", "17", as.character(n_distinct(sample_metadata$sample_id)), if_else(n_distinct(sample_metadata$sample_id) == 17L, "pass", "fail"),
  "all_libraries_have_metacells", "17", as.character(n_distinct(metacell_meta$sample_id)), if_else(n_distinct(metacell_meta$sample_id) == 17L, "pass", "fail"),
  "all_frozen_clock_classes_projected", "8", as.character(n_distinct(predictions$harmonized_cell_class)), if_else(n_distinct(predictions$harmonized_cell_class) == 8L, "pass", "fail"),
  "finite_predictions", "all", paste0(sum(is.finite(predictions$predicted_age_weeks)), "/", nrow(predictions)), if_else(all(is.finite(predictions$predicted_age_weeks)), "pass", "fail"),
  "minimum_model_feature_coverage_after_stable_id_aliasing", ">=0.80", sprintf("%.4f", min(feature_coverage$model_feature_coverage)), if_else(min(feature_coverage$model_feature_coverage) >= 0.80, "pass", "fail"),
  "unique_metacell_ids", "all", paste0(n_distinct(predictions$metacell_id), "/", nrow(predictions)), if_else(!anyDuplicated(predictions$metacell_id), "pass", "fail"),
  "matched_control_age_strata", "6;9 months", paste(sort(unique(control_baseline$age_month)), collapse = ";"), if_else(all(c(6, 9) %in% control_baseline$age_month), "pass", "fail"),
  "finite_matched_delta_rows", "all inferential rows", paste0(sum(is.finite(contrast_metacells$delta_predicted_age_weeks)), "/", nrow(contrast_metacells)), if_else(all(is.finite(contrast_metacells$delta_predicted_age_weeks)), "pass", "fail"),
  "prespecified_diet_contrasts", "4", as.character(n_distinct(contrast_metacells$contrast_id)), if_else(n_distinct(contrast_metacells$contrast_id) == 4L, "pass", "fail"),
  "complete_contrast_by_clock_audit_grid", "32 cells", paste0(sum(heatmap_data$effect_estimable), " estimable; ", sum(!heatmap_data$effect_estimable), " NE"), if_else(nrow(heatmap_data) == 32L, "pass", "fail"),
  "six_clock_metacell_level_welch_tests", "24 finite tests", paste0(sum(is.finite(metacell_level_tests_six_clocks$welch_t_p_value)), "/", nrow(metacell_level_tests_six_clocks)), if_else(nrow(metacell_level_tests_six_clocks) == 24L && all(is.finite(metacell_level_tests_six_clocks$welch_t_p_value)), "pass", "fail"),
  "frozen_model_checksum_unchanged", model_md5_before, model_md5_after, if_else(identical(model_md5_before, model_md5_after), "pass", "fail"),
  "figure_exports_nonempty", "16 nonempty PDF/PNG files", as.character(sum(file.exists(graphic_output_files) & file.info(graphic_output_files)$size > 0)), if_else(length(graphic_output_files) == 16L && all(file.exists(graphic_output_files) & file.info(graphic_output_files)$size > 0), "pass", "fail")
)
fwrite(as.data.table(validation), file.path(table_dir, "GSE235021_projection_validation.tsv"), sep = "\t", na = "")
if (any(validation$status != "pass")) {
  stop("Projection validation failed: ", paste(validation$check[validation$status != "pass"], collapse = "; "), call. = FALSE)
}

non_estimable_rows <- heatmap_data %>%
  filter(!effect_estimable) %>%
  transmute(
    item = paste0(
      "- `", as.character(contrast_label), "` × `",
      as.character(harmonized_cell_class),
      "`: not estimable because a case or age-matched polished-rice metacell group is unavailable."
    )
  ) %>%
  pull(item)
audit_lines <- c(
  "# GSE235021 Frozen Metacell-Clock Projection Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "- Projection-only external cohort; frozen July 2026 metacell clocks were not refit.",
  "- 17/17 deposited whole-ovary droplet scRNA-seq libraries were processed.",
  paste0("- ", nrow(predictions), " metacells received finite predictions across all eight clock classes."),
  paste0("- ", sum(bind_rows(qc_rows)$cells_in_metacells), " QC- and confidence-filtered cells contributed to metacells."),
  "",
  "## Gene alignment",
  "",
  "- Frozen model symbols were aligned to GSE235021 through direct symbols, reversal of R-sanitized names and stable Ensembl identifiers.",
  paste0(
    "- Per-clock model-feature coverage ranged from ",
    sprintf("%.1f%%", 100 * min(feature_coverage$model_feature_coverage)),
    " to ",
    sprintf("%.1f%%", 100 * max(feature_coverage$model_feature_coverage)),
    "."
  ),
  "- Unresolved deprecated/noncoding features were explicitly zero-filled and retained in the alias audit table.",
  "",
  "## Matched contrasts",
  "",
  "- Black-rice and white/brown-rice libraries were compared with polished-rice controls separately at 6 and 9 months.",
  "- Control baselines were defined within dataset, age stratum and cell-type clock.",
  paste0(
    "- Estimable contrast × clock combinations: ",
    sum(heatmap_data$effect_estimable), "/32."
  ),
  "",
  "Non-estimable combinations:",
  "",
  non_estimable_rows,
  "",
  "## Validation",
  "",
  paste0("- Checks passed: ", sum(validation$status == "pass"), "/", nrow(validation), "."),
  paste0("- Frozen-model checksum remained ", model_md5_after, "."),
  paste(
    "- The three six-clock figures requested here and their Welch annotations",
    "use projected metacells as the analysis unit; these are descriptive",
    "subsample-level comparisons and not biological-library-level inference."
  ),
  "- The separate sample-level Welch table retains biological-library medians as the inferential unit."
)
audit_path <- file.path(output_dir, "GSE235021_PROJECTION_AUDIT.md")
writeLines(audit_lines, audit_path)

manifest <- tibble(
  artifact = c(
    "sample_metadata", "library_qc", "cluster_marker_signature_scores",
    "cluster_cell_type_labels", "harmonized_cell_metadata", "metacell_metadata",
    "metacell_raw_counts", "model_feature_alias_map", "model_feature_coverage",
    "metacell_clock_predictions",
    "sample_celltype_prediction_summary", "matched_control_baseline",
    "metacell_delta_to_matched_controls", "sample_celltype_delta_to_matched_controls",
    "contrast_definition", "contrast_summary_by_celltype", "contrast_effect_by_celltype",
    "sample_level_welch_t_tests",
    "metacell_level_welch_t_tests_six_clocks",
    "projection_validation", "projection_audit",
    "metacell_availability_figure", "metacell_availability_summary_figure",
    "sample_predicted_age_figure",
    "metacell_delta_distribution_figure",
    "sample_delta_by_age_six_clocks_welch_tests_figure",
    "sample_delta_six_clocks_6m_horizontal_welch_tests_figure",
    "sample_delta_six_clocks_9m_horizontal_welch_tests_figure",
    "diet_effect_heatmap"
  ),
  path = c(
    file.path(table_dir, "GSE235021_sample_metadata.tsv"),
    file.path(table_dir, "GSE235021_library_qc.tsv"),
    file.path(table_dir, "GSE235021_cluster_marker_signature_scores.tsv"),
    file.path(table_dir, "GSE235021_cluster_cell_type_labels.tsv"),
    file.path(table_dir, "GSE235021_harmonized_cell_metadata.tsv.gz"),
    file.path(table_dir, "GSE235021_metacell_metadata.tsv"),
    file.path(table_dir, "GSE235021_metacell_raw_counts.rds"),
    file.path(table_dir, "GSE235021_model_feature_alias_map.tsv"),
    file.path(table_dir, "GSE235021_model_feature_coverage.tsv"),
    file.path(table_dir, "GSE235021_metacell_clock_predictions.tsv"),
    file.path(table_dir, "GSE235021_sample_celltype_prediction_summary.tsv"),
    file.path(table_dir, "GSE235021_matched_control_baseline.tsv"),
    file.path(table_dir, "GSE235021_metacell_delta_to_matched_controls.tsv"),
    file.path(table_dir, "GSE235021_sample_celltype_delta_to_matched_controls.tsv"),
    file.path(table_dir, "GSE235021_contrast_definition.tsv"),
    file.path(table_dir, "GSE235021_contrast_summary_by_celltype.tsv"),
    file.path(table_dir, "GSE235021_contrast_effect_by_celltype.tsv"),
    file.path(table_dir, "GSE235021_sample_level_welch_t_tests.tsv"),
    file.path(table_dir, "GSE235021_metacell_level_welch_t_tests_six_clocks.tsv"),
    file.path(table_dir, "GSE235021_projection_validation.tsv"),
    audit_path,
    plot_availability$pdf, plot_availability_summary$pdf,
    plot_prediction$pdf, plot_delta$pdf,
    plot_six_clock_delta_age_tests$pdf,
    plot_six_clock_delta_6m_horizontal$pdf,
    plot_six_clock_delta_9m_horizontal$pdf,
    plot_heatmap$pdf
  ),
  md5 = as.character(tools::md5sum(path)),
  bytes = file.info(path)$size
)
fwrite(as.data.table(manifest), file.path(output_dir, "GSE235021_projection_manifest.tsv"), sep = "\t", na = "")

readme <- c(
  "# GSE235021 frozen metacell-clock projection",
  "",
  "This folder contains projection-only results for all 17 whole-ovary droplet scRNA-seq libraries from GSE235021.",
  "",
  "## Frozen clocks",
  "",
  "- Eight cell-type-specific ElasticNet metacell clocks were loaded from the canonical July 2026 release.",
  "- The frozen models, training folds, coefficients and internal-CV results were not refit or modified.",
  "- Predictions and ages are reported in weeks.",
  "",
  "## External-cell processing",
  "",
  "- Each library was normalized and clustered independently in Seurat.",
  "- Broad ovarian cell classes were assigned with the same canonical marker-signature policy used for the existing external PRJNA863443 cohort.",
  sprintf("- Retained labels require confidence >= %.2f and membership in one of the eight frozen clock classes.", confidence_cutoff),
  sprintf("- Metacells were constructed within library and cell class with target %d cells and minimum %d cells.", target_cells_per_metacell, minimum_cells_per_metacell),
  "",
  "## Matched controls",
  "",
  "- Polished-rice libraries are controls.",
  "- Controls are pooled only within GSE235021, chronological-age stratum and cell-type clock.",
  "- Four prespecified contrasts are retained: black or white/brown rice versus polished rice at 6 or 9 months.",
  "- The 9-month control group contains two libraries because SCP-9M-2 is absent from the deposited whole-ovary scRNA-seq series.",
  paste(
    "- The three six-clock Welch-test figures use individual projected",
    "metacells as the analysis unit, as requested; these are descriptive",
    "subsample-level comparisons rather than biological-library-level inference."
  ),
  "- The separate sample-level Welch table uses one median per library and cell-type clock.",
  "",
  "## Outputs",
  "",
  "- `GSE235021_metacell_availability_by_library.*`",
  "- `GSE235021_metacell_availability_summary.*`",
  "- `GSE235021_sample_celltype_predicted_age.*`",
  "- `GSE235021_metacell_delta_distribution_by_celltype.*`",
  "- `GSE235021_sample_delta_by_age_six_clocks_welch_tests.*`",
  "- `GSE235021_sample_delta_six_clocks_6m_horizontal_welch_tests.*`",
  "- `GSE235021_sample_delta_six_clocks_9m_horizontal_welch_tests.*`",
  "- `GSE235021_diet_effect_heatmap_by_celltype.*`",
  "- `tables/`: source tables, predictions, matched-control effects, tests and validation.",
  "- `logs/`: reproducible run log.",
  "- `GSE235021_PROJECTION_AUDIT.md`: projection scope, gene-alignment and estimability audit.",
  "- `GSE235021_projection_manifest.tsv`: checksums and sizes for the release artifacts.",
  "",
  "## Validation",
  "",
  paste0("- Libraries represented: ", n_distinct(predictions$sample_id), "/17."),
  paste0("- Projected metacells: ", nrow(predictions), "."),
  paste0("- Frozen clock classes represented: ", n_distinct(predictions$harmonized_cell_class), "/8."),
  paste0("- Estimable contrast × clock combinations: ", sum(heatmap_data$effect_estimable), "/32; unavailable combinations are marked NE."),
  paste0("- Minimum frozen-model feature coverage: ", sprintf("%.2f%%", 100 * min(feature_coverage$model_feature_coverage)), "."),
  paste0("- Validation checks passed: ", sum(validation$status == "pass"), "/", nrow(validation), ".")
)
writeLines(readme, file.path(output_dir, "README.md"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

figure_manifest_path <- file.path(run_root, "step_08_figures_tables_release", "figure_manifest.tsv")
figure_manifest <- fread(figure_manifest_path, na.strings = c("", "NA")) %>% as_tibble()
gse_entries <- tribble(
  ~figure_id, ~figure_title, ~stem,
  "GSE235021_metacell_availability_by_library", "GSE235021 projected metacell availability by library", "GSE235021_metacell_availability_by_library",
  "GSE235021_metacell_availability_summary", "GSE235021 projected metacell availability summary", "GSE235021_metacell_availability_summary",
  "GSE235021_sample_celltype_predicted_age", "GSE235021 sample-level median predicted age by cell-type clock", "GSE235021_sample_celltype_predicted_age",
  "GSE235021_metacell_delta_distribution_by_celltype", "GSE235021 metacell delta-age distributions by matched diet contrast", "GSE235021_metacell_delta_distribution_by_celltype",
  "GSE235021_sample_delta_by_age_six_clocks_welch_tests", "GSE235021 age-stratified metacell delta across six clocks with metacell-level Welch tests", "GSE235021_sample_delta_by_age_six_clocks_welch_tests",
  "GSE235021_sample_delta_six_clocks_6m_horizontal_welch_tests", "GSE235021 6-month metacell delta across six clocks in horizontal matched-comparison format", "GSE235021_sample_delta_six_clocks_6m_horizontal_welch_tests",
  "GSE235021_sample_delta_six_clocks_9m_horizontal_welch_tests", "GSE235021 9-month metacell delta across six clocks in horizontal matched-comparison format", "GSE235021_sample_delta_six_clocks_9m_horizontal_welch_tests",
  "GSE235021_diet_effect_heatmap_by_celltype", "GSE235021 matched diet effects across frozen metacell clocks", "GSE235021_diet_effect_heatmap_by_celltype"
) %>%
  transmute(
    figure_id, figure_title,
    figure_path_png = file.path("figures", "05_stress_test_projection", "GSE235021", paste0(stem, ".png")),
    figure_path_pdf = file.path("figures", "05_stress_test_projection", "GSE235021", paste0(stem, ".pdf")),
    plot_data_path = file.path("figures", "05_stress_test_projection", "GSE235021", "tables", paste0(stem, "_plot_data.tsv")),
    source_table_path = file.path(table_dir, "GSE235021_metacell_clock_predictions.tsv"),
    section = "05_stress_test_projection",
    status = "complete_external_GSE235021"
  )
figure_manifest <- bind_rows(
  figure_manifest %>% filter(!figure_id %in% gse_entries$figure_id),
  gse_entries
) %>%
  arrange(section, figure_id)
fwrite(as.data.table(figure_manifest), figure_manifest_path, sep = "\t", na = "")

message("GSE235021 frozen metacell-clock projection completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
message("Libraries: ", n_distinct(predictions$sample_id))
message("Projected metacells: ", nrow(predictions))
message("Cell-type clocks: ", n_distinct(predictions$harmonized_cell_class))
message("Validation: ", sum(validation$status == "pass"), "/", nrow(validation), " passed")
