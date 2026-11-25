#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
script_path <- gsub("~\\+~", " ", script_path)
step_dir <- if (length(script_path) == 1) dirname(normalizePath(script_path)) else normalizePath(".")
out_dir <- file.path(step_dir, "somatic_marker_expression_function_panel")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(step_dir, "filtered_metacell_counts.rds")
metadata_path <- file.path(step_dir, "filtered_metacell_metadata.tsv")
stopifnot(file.exists(counts_path), file.exists(metadata_path))

celltype_order <- c(
  "Granulosa", "Theca_Stroma", "Endothelial", "Pericyte_Smooth_Muscle",
  "Epithelial", "Immune_Myeloid", "T_NK", "B_Plasma"
)
celltype_labels <- c(
  Granulosa = "Granulosa",
  Theca_Stroma = "Theca / stroma",
  Endothelial = "Endothelial",
  Pericyte_Smooth_Muscle = "Pericyte / smooth muscle",
  Epithelial = "Epithelial",
  Immune_Myeloid = "Immune / myeloid",
  T_NK = "T / NK",
  B_Plasma = "B / plasma"
)
celltype_colors <- c(
  Granulosa = "#D98972",
  Theca_Stroma = "#C9A457",
  Endothelial = "#4F9D8D",
  Pericyte_Smooth_Muscle = "#6D8F73",
  Epithelial = "#7E9FC5",
  Immune_Myeloid = "#8C79B5",
  T_NK = "#B477A4",
  B_Plasma = "#B76B75"
)

top_n <- 50L
min_samples_per_group <- 3L
min_cpm <- 1
min_positive_datasets <- 2L
top_terms <- 3L

message("Reading filtered metacell inputs")
counts <- readRDS(counts_path)
metadata <- fread(metadata_path)
stopifnot(inherits(counts, "Matrix"), ncol(counts) == nrow(metadata))
if (!identical(colnames(counts), metadata$metacell_id)) {
  idx <- match(colnames(counts), metadata$metacell_id)
  if (anyNA(idx)) stop("Metacell count columns do not match metadata IDs", call. = FALSE)
  metadata <- metadata[idx]
}

keep <- metadata$model_role == "clock_candidate" &
  metadata$include_clock_candidate %in% TRUE &
  metadata$cv_role == "metacell_cv_candidate" &
  metadata$harmonized_cell_class %in% celltype_order
metadata <- metadata[keep]
counts <- counts[, keep, drop = FALSE]
metadata[, harmonized_cell_class := factor(harmonized_cell_class, levels = celltype_order)]

message("Aggregating metacells to biological sample-by-cell-type pseudobulks")
aggregate_key <- interaction(
  metadata$dataset_id,
  metadata$sample_key,
  metadata$harmonized_cell_class,
  drop = TRUE,
  lex.order = TRUE
)
membership <- sparse.model.matrix(~ 0 + aggregate_key)
sample_counts <- counts %*% membership
aggregate_levels <- levels(aggregate_key)
colnames(sample_counts) <- aggregate_levels

sample_meta <- metadata[, .(
  dataset_id = first(dataset_id),
  dataset_folder = first(dataset_folder),
  sample_key = first(sample_key),
  sample_id = first(sample_id),
  harmonized_cell_class = as.character(first(harmonized_cell_class)),
  n_metacells = .N,
  n_cells = sum(n_cells),
  age_weeks = first(age_weeks)
), by = aggregate_key]
sample_meta <- sample_meta[match(aggregate_levels, aggregate_key)]
sample_meta[, aggregate_key := as.character(aggregate_key)]

fwrite(sample_meta, file.path(out_dir, "sample_celltype_pseudobulk_metadata.tsv"), sep = "\t")
saveRDS(sample_counts, file.path(out_dir, "sample_celltype_pseudobulk_counts.rds"), compress = "xz")

message("Running a dataset-adjusted voom/limma model with one-versus-rest contrasts")
dge_all <- DGEList(counts = sample_counts)
dge_all <- calcNormFactors(dge_all)
logcpm_all <- cpm(dge_all, log = TRUE, prior.count = 2)

celltype_factor <- factor(sample_meta$harmonized_cell_class, levels = celltype_order)
dataset_factor <- factor(sample_meta$dataset_id)
design <- model.matrix(~ 0 + celltype_factor + dataset_factor)
colnames(design) <- make.names(colnames(design))
keep_gene <- filterByExpr(dge_all, group = celltype_factor, min.count = 5) &
  rowSums(cpm(dge_all) >= min_cpm) >= min_samples_per_group
dge_fit <- dge_all[keep_gene, , keep.lib.sizes = FALSE]
dge_fit <- calcNormFactors(dge_fit)
voom_fit <- voom(dge_fit, design, plot = FALSE)
linear_fit <- lmFit(voom_fit, design)

run_marker_test <- function(target_celltype) {
  target <- sample_meta$harmonized_cell_class == target_celltype
  dataset_target <- tapply(target, sample_meta$dataset_id, function(x) length(unique(x)) > 1)
  informative_datasets <- names(dataset_target)[dataset_target]
  usable <- sample_meta$dataset_id %in% informative_datasets
  target_usable <- target[usable]
  if (sum(target_usable) < min_samples_per_group || sum(!target_usable) < min_samples_per_group) {
    stop("Insufficient samples for marker test: ", target_celltype, call. = FALSE)
  }

  target_coef <- make.names(paste0("celltype_factor", target_celltype))
  other_coefs <- make.names(paste0("celltype_factor", setdiff(celltype_order, target_celltype)))
  contrast_string <- paste0(target_coef, "-(`", paste(other_coefs, collapse = "+"), "`)/", length(other_coefs))
  contrast_string <- gsub("`", "", contrast_string, fixed = TRUE)
  contrast <- makeContrasts(contrasts = contrast_string, levels = design)
  fit <- contrasts.fit(linear_fit, contrast)
  fit <- eBayes(fit, trend = FALSE, robust = FALSE)
  result <- topTable(fit, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("gene") %>%
    as_tibble() %>%
    transmute(gene, logFC, logCPM = AveExpr, PValue = P.Value, FDR = adj.P.Val)

  per_dataset <- lapply(informative_datasets, function(dataset_name) {
    idx <- sample_meta$dataset_id == dataset_name
    target_idx <- idx & target
    other_idx <- idx & !target
    if (!any(target_idx) || !any(other_idx)) return(NULL)
    target_mean <- rowMeans(logcpm_all[, target_idx, drop = FALSE])
    other_mean <- rowMeans(logcpm_all[, other_idx, drop = FALSE])
    tibble(
      gene = rownames(logcpm_all),
      dataset_id = dataset_name,
      dataset_logFC = target_mean - other_mean,
      positive = dataset_logFC > 0
    )
  }) %>% bind_rows()

  consistency <- per_dataset %>%
    group_by(gene) %>%
    summarise(
      n_informative_datasets = n(),
      n_positive_datasets = sum(positive),
      median_dataset_logFC = median(dataset_logFC),
      min_dataset_logFC = min(dataset_logFC),
      .groups = "drop"
    )

  target_detection <- rowMeans(cpm(dge_all)[, target, drop = FALSE] >= min_cpm)
  other_detection <- rowMeans(cpm(dge_all)[, !target, drop = FALSE] >= min_cpm)

  result %>%
    left_join(consistency, by = "gene") %>%
    mutate(
      harmonized_cell_class = target_celltype,
      target_sample_detection = target_detection[gene],
      other_sample_detection = other_detection[gene],
      detection_difference = target_sample_detection - other_sample_detection,
      specificity_score = pmax(logFC, 0) * pmax(detection_difference, 0) *
        (n_positive_datasets / pmax(n_informative_datasets, 1)),
      marker_eligible = logFC > 0 & FDR < 0.05 & detection_difference > 0 &
        n_positive_datasets >= pmin(min_positive_datasets, n_informative_datasets)
    ) %>%
    arrange(desc(marker_eligible), desc(specificity_score), FDR, desc(logFC), gene)
}

marker_results <- bind_rows(lapply(celltype_order, run_marker_test))
fwrite(marker_results, file.path(out_dir, "somatic_celltype_marker_tests.tsv"), sep = "\t")

message("Assigning unique top-50 markers to each somatic cell type")
candidates <- marker_results %>%
  filter(marker_eligible) %>%
  group_by(harmonized_cell_class) %>%
  mutate(within_celltype_rank = row_number()) %>%
  ungroup() %>%
  arrange(within_celltype_rank, desc(specificity_score), FDR, gene)

selected <- vector("list", length(celltype_order))
names(selected) <- celltype_order
used_genes <- character()
for (rank_round in seq_len(max(candidates$within_celltype_rank))) {
  for (celltype in celltype_order) {
    if (length(selected[[celltype]]) >= top_n) next
    row <- candidates %>%
      filter(harmonized_cell_class == celltype, within_celltype_rank == rank_round, !gene %in% used_genes)
    if (nrow(row) == 0) next
    selected[[celltype]] <- c(selected[[celltype]], row$gene[[1]])
    used_genes <- c(used_genes, row$gene[[1]])
  }
}
if (any(lengths(selected) < top_n)) {
  short <- paste(names(selected)[lengths(selected) < top_n], lengths(selected)[lengths(selected) < top_n], sep = "=", collapse = "; ")
  stop("Could not assign 50 unique eligible markers per cell type: ", short, call. = FALSE)
}

top_markers <- bind_rows(lapply(celltype_order, function(celltype) {
  marker_results %>%
    filter(harmonized_cell_class == celltype, gene %in% selected[[celltype]]) %>%
    mutate(marker_rank = match(gene, selected[[celltype]])) %>%
    arrange(marker_rank)
}))
fwrite(top_markers, file.path(out_dir, "top50_unique_somatic_markers_by_celltype.tsv"), sep = "\t")

message("Preparing cell-type mean expression and row-scaled heatmap data")
mean_logcpm <- sapply(celltype_order, function(celltype) {
  rowMeans(logcpm_all[, sample_meta$harmonized_cell_class == celltype, drop = FALSE])
})
rownames(mean_logcpm) <- rownames(logcpm_all)
ordered_genes <- unlist(selected[celltype_order], use.names = FALSE)
heat_matrix <- mean_logcpm[ordered_genes, celltype_order, drop = FALSE]
heat_z <- t(scale(t(heat_matrix)))
heat_z[!is.finite(heat_z)] <- 0
heat_z <- pmax(-2.5, pmin(2.5, heat_z))
heat_z <- matrix(
  as.numeric(heat_z),
  nrow = length(ordered_genes),
  ncol = length(celltype_order),
  dimnames = list(ordered_genes, celltype_order)
)

heat_data <- tibble(
  gene = rep(rownames(heat_z), times = ncol(heat_z)),
  expression_celltype = rep(colnames(heat_z), each = nrow(heat_z)),
  z_score = as.vector(heat_z)
)
marker_map <- top_markers %>% select(gene, marker_celltype = harmonized_cell_class, marker_rank)
heat_data <- heat_data %>%
  left_join(marker_map, by = "gene") %>%
  mutate(
    marker_celltype = factor(marker_celltype, levels = rev(celltype_order)),
    expression_celltype = factor(expression_celltype, levels = celltype_order),
    gene_order = factor(gene, levels = rev(ordered_genes)),
    display_gene = if_else(marker_rank <= 6, gene, "")
  )
fwrite(heat_data, file.path(out_dir, "marker_heatmap_source_data.tsv"), sep = "\t")

message("Loading mouse-native MSigDB pathways and running ORA")
temp_r_lib <- Sys.getenv("ADULT_C57_GSEA_R_LIB", unset = "/tmp/adult_c57_r_gsea_lib")
if (dir.exists(temp_r_lib)) .libPaths(c(temp_r_lib, .libPaths()))
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("Missing msigdbr in ADULT_C57_GSEA_R_LIB or the active R library", call. = FALSE)
}
msig <- msigdbr::msigdbr(db_species = "MM", species = "Mus musculus") %>%
  filter(
    (gs_collection == "M5" & gs_subcollection == "GO:BP") |
      (gs_collection == "M2" & gs_subcollection == "CP:REACTOME")
  ) %>%
  transmute(
    gene_set_source = if_else(gs_collection == "M5", "GO:BP", "REACTOME"),
    term_id = gs_id,
    term_name = str_replace_all(gs_name, "^(GOBP_|REACTOME_)", ""),
    gene_symbol
  ) %>%
  distinct()

background <- intersect(unique(marker_results$gene[marker_results$marker_eligible | marker_results$FDR < 0.05]), unique(msig$gene_symbol))
pathways <- msig %>%
  filter(gene_symbol %in% background) %>%
  group_by(gene_set_source, term_id, term_name) %>%
  summarise(pathway_genes = list(sort(unique(gene_symbol))), pathway_size = n_distinct(gene_symbol), .groups = "drop") %>%
  filter(pathway_size >= 10, pathway_size <= 500)

ora_results <- bind_rows(lapply(celltype_order, function(celltype) {
  query <- intersect(selected[[celltype]], background)
  bind_rows(lapply(seq_len(nrow(pathways)), function(i) {
    overlap <- intersect(query, pathways$pathway_genes[[i]])
    p_value <- phyper(length(overlap) - 1, pathways$pathway_size[[i]], length(background) - pathways$pathway_size[[i]], length(query), lower.tail = FALSE)
    tibble(
      harmonized_cell_class = celltype,
      gene_set_source = pathways$gene_set_source[[i]],
      term_id = pathways$term_id[[i]],
      term_name = pathways$term_name[[i]],
      pathway_size = pathways$pathway_size[[i]],
      query_size = length(query),
      overlap_size = length(overlap),
      overlap_genes = paste(sort(overlap), collapse = ";"),
      p_value = p_value
    )
  })) %>% mutate(p_adjust = p.adjust(p_value, method = "BH"))
})) %>%
  mutate(enrichment_score = -log10(pmax(p_adjust, .Machine$double.xmin))) %>%
  arrange(harmonized_cell_class, p_adjust, desc(overlap_size), term_name)
fwrite(ora_results, file.path(out_dir, "top50_marker_msigdb_ora_results.tsv"), sep = "\t")

select_nonredundant_terms <- function(tbl, n_keep = 3L) {
  tbl <- tbl %>% filter(p_adjust < 0.05, overlap_size >= 2)
  if (nrow(tbl) == 0) return(tbl)
  kept <- integer()
  kept_sets <- list()
  for (i in seq_len(nrow(tbl))) {
    genes <- str_split(tbl$overlap_genes[[i]], ";", simplify = FALSE)[[1]]
    redundant <- any(vapply(kept_sets, function(existing) {
      length(intersect(genes, existing)) / max(1, length(union(genes, existing))) >= 0.5
    }, logical(1)))
    if (!redundant) {
      kept <- c(kept, i)
      kept_sets[[length(kept_sets) + 1L]] <- genes
    }
    if (length(kept) >= n_keep) break
  }
  tbl[kept, , drop = FALSE]
}

select_display_terms <- function(tbl, n_keep = 3L) {
  selected_terms <- select_nonredundant_terms(tbl, n_keep)
  if (nrow(selected_terms) > 0) return(selected_terms)
  tbl %>%
    slice(1) %>%
    mutate(
      gene_set_source = "none",
      term_id = "none",
      term_name = "NO_FDR_SIGNIFICANT_TERM",
      pathway_size = 0L,
      query_size = top_n,
      overlap_size = 0L,
      overlap_genes = "",
      p_value = 1,
      p_adjust = 1,
      enrichment_score = 0
    )
}

display_terms <- ora_results %>%
  group_by(harmonized_cell_class) %>%
  group_modify(~ select_display_terms(.x, top_terms)) %>%
  ungroup() %>%
  mutate(
    harmonized_cell_class = factor(harmonized_cell_class, levels = rev(celltype_order)),
    significant_fdr_0.05 = p_adjust < 0.05,
    term_label = str_to_sentence(str_replace_all(term_name, "_", " ")),
    term_plot_order = interaction(harmonized_cell_class, row_number(), lex.order = TRUE)
  )
if (nrow(display_terms) == 0) stop("No significant representative marker terms were found", call. = FALSE)
fwrite(display_terms, file.path(out_dir, "top50_marker_representative_terms.tsv"), sep = "\t")

message("Drawing marker heatmap and representative-term bars")
gene_labels <- setNames(ifelse(top_markers$marker_rank <= 6, top_markers$gene, ""), top_markers$gene)
p_heat <- ggplot(heat_data, aes(expression_celltype, gene_order, fill = z_score)) +
  geom_tile() +
  facet_grid(marker_celltype ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#2C6AA0", mid = "#F7F7F7", high = "#B53A3A", midpoint = 0, limits = c(-2.5, 2.5), name = "Row z-score") +
  scale_x_discrete(labels = celltype_labels) +
  scale_y_discrete(labels = gene_labels, expand = c(0, 0)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 7, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(size = 4.6, colour = "black"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 6.5),
    strip.background = element_rect(fill = "white", colour = NA),
    panel.spacing.y = unit(0.7, "mm"),
    legend.position = "bottom",
    plot.margin = margin(2, 1, 2, 2)
  )

display_terms <- display_terms %>%
  group_by(harmonized_cell_class) %>%
  arrange(enrichment_score, .by_group = TRUE) %>%
  mutate(
    term_plot_label = str_wrap(term_label, width = 31),
    term_factor = factor(term_plot_label, levels = unique(term_plot_label))
  ) %>%
  ungroup()
p_terms <- ggplot(display_terms, aes(enrichment_score, term_factor, fill = as.character(harmonized_cell_class), alpha = significant_fdr_0.05)) +
  geom_col(width = 0.68) +
  facet_grid(harmonized_cell_class ~ ., scales = "free_y", space = "fixed") +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.45), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = expression(-log[10](adjusted~italic(P))), y = NULL) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.text.y = element_text(size = 5.2, colour = "black"),
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    panel.spacing.y = unit(0.7, "mm"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(2, 2, 2, 1)
  )

combined <- p_heat + p_terms + plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = list("F")) &
  theme(plot.tag = element_text(face = "bold", size = 9))

width_in <- 11.2
height_in <- 13.5
ggsave(file.path(out_dir, "somatic_top50_marker_heatmap_and_terms.pdf"), combined, width = width_in, height = height_in, device = cairo_pdf)
ggsave(file.path(out_dir, "somatic_top50_marker_heatmap_and_terms.svg"), combined, width = width_in, height = height_in, device = grDevices::svg)
grDevices::png(file.path(out_dir, "somatic_top50_marker_heatmap_and_terms.png"), width = width_in, height = height_in, units = "in", res = 300, type = "cairo", bg = "white")
print(combined)
dev.off()
grDevices::tiff(file.path(out_dir, "somatic_top50_marker_heatmap_and_terms.tiff"), width = width_in, height = height_in, units = "in", res = 600, compression = "lzw", type = "cairo", bg = "white")
print(combined)
dev.off()

audit <- tibble(
  check = c(
    "input_metacells", "analysis_metacells", "analysis_sample_celltype_pseudobulks",
    "somatic_celltypes", "unique_top_markers", "markers_per_celltype",
    "all_markers_positive_logfc", "all_markers_fdr_below_0.05",
    "all_markers_supported_by_required_datasets", "representative_terms"
  ),
  value = c(
    ncol(readRDS(counts_path)), ncol(counts), ncol(sample_counts), length(celltype_order),
    n_distinct(top_markers$gene), paste(table(top_markers$harmonized_cell_class), collapse = ";"),
    all(top_markers$logFC > 0), all(top_markers$FDR < 0.05),
    all(top_markers$n_positive_datasets >= pmin(min_positive_datasets, top_markers$n_informative_datasets)),
    nrow(display_terms)
  ),
  status = c(rep("information", 6), rep("pass", 3), "information")
)
fwrite(audit, file.path(out_dir, "marker_panel_hard_checks.tsv"), sep = "\t")

session <- capture.output(sessionInfo())
writeLines(session, file.path(out_dir, "session_info.txt"))
message("Completed: ", out_dir)
