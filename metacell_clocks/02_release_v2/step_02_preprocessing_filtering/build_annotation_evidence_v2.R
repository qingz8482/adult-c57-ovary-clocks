#!/usr/bin/env Rscript
# Annotation-evidence figures for the v2 labels: quantitative alternatives to UMAP for showing that cell
# classes are correct and consistent across datasets.
#   A  canonical marker dot plot (pooled across datasets) + dataset-resolved signature heatmap
#   B  benchmark against author labels (Morris <-> Isola), row-normalised confusion matrices
#   C  cross-dataset cluster similarity heatmap ordered by final label
#   D  evidence source and confidence by dataset
#   E  v1 -> v2 relabelling matrix
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(Matrix) })
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
step_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce/step_02_preprocessing_filtering")
out_dir <- file.path(step_dir, "figures", "annotation_evidence_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

classes <- c("Granulosa", "Theca_Stroma", "Endothelial", "Pericyte_Smooth_Muscle", "Epithelial",
             "Immune_Myeloid", "T_NK", "B_Plasma", "Oocyte_Germline")
class_colors <- c(Granulosa = "#E7298A", Theca_Stroma = "#B07AE8", Endothelial = "#1FA5C9",
                  Pericyte_Smooth_Muscle = "#6A3D9A", Epithelial = "#7FB800", Immune_Myeloid = "#E06A00",
                  T_NK = "#0067A5", B_Plasma = "#1B9E77", Oocyte_Germline = "#F2C14E", Erythroid = "#C0392B",
                  Unknown = "grey70")
ds_short <- function(x) {
  x <- sub("^(20[0-9]{2})_[A-Za-z]+_([A-Za-z]+)_.*$", "\\2 \\1", x)
  sub("^(20[0-9]{2})_Kim_PRJNA863443.*$", "Kim 2025", x)
}
save_panel <- function(stem, p, w, h) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 200, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf, limitsize = FALSE)
}
theme_ev <- theme_bw(base_size = 10) + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92"))

cons <- fread_gz(file.path(V2, "04_consensus", "harmonized_annotation_consensus_v2.tsv.gz"), colClasses = list(character = "standard_cluster"))
cons[, eligible := curation_status %in% c("confirmed", "probable") & annotation_confidence >= 0.70 &
       !annotation_label_broad %in% c("Unknown", "Erythroid")]
panel <- panel_v2()
sig_order <- unique(panel$signature)
genes <- unique(panel$gene)

# ---- A: marker expression by class, per dataset ------------------------------------------------------------
expr_path <- file.path(out_dir, "A_marker_expression_by_dataset_class.tsv")
if (!file.exists(expr_path)) {
  rows <- list()
  for (ds in DATASETS) {
    message("A: ", ds)
    o <- load_query(ds); o <- NormalizeData(o, verbose = FALSE)
    X <- LayerData(o, layer = "data")
    lab <- cons[dataset_id == ds & eligible == TRUE, .(cell_id, annotation_label_broad)]
    lab <- lab[cell_id %in% colnames(X)]
    g <- intersect(genes, rownames(X))
    for (k in intersect(classes, unique(lab$annotation_label_broad))) {
      cc <- lab[annotation_label_broad == k, cell_id]
      if (length(cc) < 20) next
      M <- X[g, cc, drop = FALSE]
      rows[[length(rows) + 1]] <- data.table(dataset_id = ds, cell_class = k, n_cells = length(cc), gene = g,
                                             mean_expr = Matrix::rowMeans(M), pct_expr = Matrix::rowMeans(M > 0))
    }
    rm(o, X); gc(verbose = FALSE)
  }
  fwrite(rbindlist(rows), expr_path, sep = "\t")
}
ex <- fread(expr_path)
ex[, z := { s <- sd(mean_expr); if (is.na(s) || s == 0) 0 else (mean_expr - mean(mean_expr)) / s }, by = .(dataset_id, gene)]
pool <- ex[, .(z = mean(z), pct = mean(pct_expr), n_datasets = .N), by = .(cell_class, gene)]
gene_sig <- unique(panel[, .(gene, signature)])[!duplicated(gene)]
pool <- merge(pool, gene_sig, by = "gene")
pool[, gene := factor(gene, levels = unique(panel$gene))]
pool[, signature := factor(signature, levels = sig_order)]
pool[, cell_class := factor(cell_class, levels = rev(classes))]
fwrite(pool, file.path(out_dir, "A_marker_dotplot_pooled.tsv"), sep = "\t")
pA <- ggplot(pool, aes(gene, cell_class, size = pct, colour = pmin(pmax(z, -1.5), 2.5))) + geom_point() +
  facet_grid(~signature, scales = "free_x", space = "free_x") +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey92", high = "#B40426", midpoint = 0, name = "Mean z\n(within dataset,\naveraged)") +
  scale_size_area(max_size = 4, name = "Fraction\nexpressing", labels = scales::percent) + theme_ev +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7, face = "italic"),
        strip.text.x = element_text(angle = 90, size = 7), panel.spacing.x = unit(1, "pt")) +
  labs(x = NULL, y = NULL, title = "A  Canonical marker expression in v2 cell classes (confident cells, 8 datasets pooled)",
       subtitle = "Each gene z-scored across classes within a dataset, then averaged; dot size = mean fraction of expressing cells")
save_panel("A_marker_dotplot_pooled", pA, 20, 6)

sig_ds <- merge(ex, gene_sig, by = "gene")[, .(score = mean(z)), by = .(dataset_id, cell_class, signature)]
sig_ds[, row := paste(ds_short(dataset_id), sep = "")]
sig_ds[, cell_class := factor(cell_class, levels = classes)]
sig_ds[, signature := factor(signature, levels = sig_order)]
fwrite(sig_ds, file.path(out_dir, "A_signature_score_by_dataset_class.tsv"), sep = "\t")
pA2 <- ggplot(sig_ds, aes(signature, row, fill = pmin(pmax(score, -1.5), 2.5))) + geom_tile(colour = "white") +
  facet_grid(cell_class ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "grey95", high = "#B40426", midpoint = 0, name = "Signature\nscore (z)") +
  theme_ev + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text.y.left = element_text(angle = 0),
                   strip.placement = "outside", axis.text.y = element_text(size = 6)) +
  labs(x = "Marker signature", y = NULL, title = "A2  Signature scores per dataset and v2 class",
       subtitle = "A correct and consistent annotation shows the matching signature high in every dataset")
save_panel("A2_signature_score_by_dataset_class", pA2, 10, 14)

# ---- B: benchmark confusion --------------------------------------------------------------------------------
cf <- fread(file.path(V2, "03_benchmark", "benchmark_confusion.tsv"))
cf[, frac := N / sum(N), by = .(query, method, truth)]
cf[, panel := paste0(ds_short(query), " author labels vs ", fifelse(method == "reference", "reference mapping", "marker scoring"),
                     fifelse(method == "reference", fifelse(grepl("Morris", query), " (Isola reference)", " (Morris reference)"), ""))]
lv <- c(classes, "Unknown")
cf[, truth := factor(truth, levels = rev(lv))][, pred := factor(pred, levels = lv)]
pB <- ggplot(cf, aes(pred, truth, fill = frac)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(frac >= 0.02, sprintf("%.2f", frac), "")), size = 2.4) +
  facet_wrap(~panel, ncol = 2) + scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), name = "Row fraction") +
  theme_ev + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Predicted class", y = "Author class", title = "B  Benchmark against author labels (each atlas predicted without its own labels)",
       subtitle = "Ambiguous Isola clusters (Stroma B, Stroma C, Epithelium B) excluded from truth")
save_panel("B_benchmark_confusion_vs_author_labels", pB, 13, 11)

# ---- C: cross-dataset cluster similarity ---------------------------------------------------------------------
prof <- readRDS(file.path(V2, "04_consensus", "cluster_profiles_refgenes.rds"))
P <- do.call(cbind, lapply(prof, function(M) M - rowMeans(M)))
P <- P[stats::complete.cases(P), ]
C <- cor(P, method = "spearman")
cl <- fread(file.path(V2, "04_consensus", "cluster_consensus_v2.tsv"), colClasses = list(character = "cluster"))
cl[, key := paste(dataset_id, cluster, sep = "::")]
cl <- cl[key %in% colnames(C) & n_cells >= 50]
cl[, cluster_label := factor(cluster_label, levels = c(classes, "Erythroid", "Unknown"))]
setorder(cl, cluster_label, dataset_id)
# within each label, order by hierarchical clustering
ord <- unlist(lapply(split(cl$key, cl$cluster_label, drop = TRUE), function(k)
  if (length(k) > 2) k[hclust(as.dist(1 - C[k, k]))$order] else k))
cd <- as.data.table(as.table(C[ord, ord]))
setnames(cd, c("a", "b", "rho"))
cd[, a := factor(a, levels = ord)][, b := factor(b, levels = rev(ord))]
ann <- cl[match(ord, key), .(key = factor(key, levels = ord), cluster_label)]
pC <- ggplot(cd, aes(a, b, fill = rho)) + geom_raster() +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-1, 1), name = "Spearman") +
  theme_void(base_size = 10) + theme(legend.position = "right") + coord_fixed() +
  labs(title = "C  Cross-dataset similarity of all clusters (>=50 cells), ordered by v2 label",
       subtitle = "Dataset-centred cluster mean profiles on reference variable genes; blocks on the diagonal = same state across datasets")
bar <- ggplot(ann, aes(key, 1, fill = cluster_label)) + geom_tile() + scale_fill_manual(values = class_colors, name = "v2 label") +
  theme_void() + theme(legend.position = "bottom")
save_panel("C_cross_dataset_cluster_similarity", pC / bar + plot_layout(heights = c(30, 1)), 12, 13.5)
fwrite(fread(file.path(V2, "04_consensus", "cross_dataset_consistency_summary.tsv")), file.path(out_dir, "C_consistency_summary.tsv"), sep = "\t")

# ---- D: evidence source and confidence by dataset ----------------------------------------------------------
cons[, evidence := fcase(
  annotation_label_broad == "Erythroid", "Erythroid (excluded)",
  annotation_label_broad == "Unknown", "Unknown (excluded)",
  !eligible, "Below confidence 0.70 (excluded)",
  annotation_source == "author_metadata", "Author label, supported",
  annotation_source == "author_metadata_unsupported", "Author label only",
  annotation_source == "reference+marker", "Reference + marker agree",
  annotation_source == "curated", "Manual review",
  default = "One automated method")]
ev_levels <- c("Author label, supported", "Author label only", "Reference + marker agree", "One automated method",
               "Manual review", "Below confidence 0.70 (excluded)", "Unknown (excluded)", "Erythroid (excluded)")
ev_cols <- c("#1B7837", "#A6DBA0", "#2166AC", "#92C5DE", "#F4A582", "grey60", "grey80", "#C0392B")
ev <- cons[, .N, by = .(dataset = ds_short(dataset_id), evidence)]
ev[, frac := N / sum(N), by = dataset][, evidence := factor(evidence, levels = ev_levels)]
fwrite(ev, file.path(out_dir, "D_evidence_by_dataset.tsv"), sep = "\t")
pD1 <- ggplot(ev, aes(dataset, frac, fill = evidence)) + geom_col(width = 0.8) + coord_flip() +
  scale_fill_manual(values = setNames(ev_cols, ev_levels), name = NULL) + scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme_ev + labs(x = NULL, y = "Fraction of cells", title = "D  Evidence supporting each label")
pD2 <- ggplot(cons[eligible == TRUE], aes(ds_short(dataset_id), annotation_confidence)) +
  geom_violin(fill = "grey85", colour = "grey40", scale = "width") + coord_flip() + theme_ev +
  labs(x = NULL, y = "Calibrated label confidence (clock-eligible cells)", title = "Confidence of retained labels")
save_panel("D_evidence_and_confidence_by_dataset", pD1 | pD2, 15, 5)

# ---- E: v1 -> v2 relabelling ---------------------------------------------------------------------------------
v1 <- fread(file.path(V2, "..", "v2", "00_input_seurat", "harmonized_annotation_consensus.tsv"))
v1 <- v1[, .(cell_id, v1_label = annotation_label_broad,
             v1_eligible = curation_status %in% c("confirmed", "probable") & annotation_confidence >= 0.70 &
               !annotation_label_broad %in% c("Unknown", "Erythroid"))]
x <- merge(cons[, .(cell_id, v2_label = annotation_label_broad, v2_eligible = eligible)], v1, by = "cell_id")
x[, v1c := fifelse(v1_eligible, v1_label, "Excluded")][, v2c := fifelse(v2_eligible, v2_label, "Excluded")]
tab <- x[, .N, by = .(v1c, v2c)][, frac := N / sum(N), by = v1c]
fwrite(tab, file.path(out_dir, "E_v1_to_v2_relabelling.tsv"), sep = "\t")
lv2 <- c(classes, "Excluded")
tab[, v1c := factor(v1c, levels = rev(lv2))][, v2c := factor(v2c, levels = lv2)]
pE <- ggplot(tab, aes(v2c, v1c, fill = frac)) + geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(N >= 100, format(N, big.mark = ","), "")), size = 2.3) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), name = "Row fraction") + theme_ev +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "v2 class", y = "v1 class", title = "E  Relabelling from v1 to v2 (7 datasets with v1 consensus labels)",
       subtitle = "Counts shown where >= 100 cells; 'Excluded' = Unknown, Erythroid or below the clock confidence policy")
save_panel("E_v1_to_v2_relabelling", pE, 10, 8)

writeLines(c("# Annotation evidence (v2 labels)", "",
  "Built by `../../build_annotation_evidence_v2.R`. These panels replace UMAP-only presentation of the annotation:",
  "", "- A/A2: canonical marker expression by class (pooled and per dataset): shows *what* defines each class and that it holds in every dataset.",
  "- B: confusion matrices against author labels: quantitative accuracy of each automated method.",
  "- C: cross-dataset cluster similarity ordered by label: same states receive the same label across studies.",
  "- D: evidence source and calibrated confidence per dataset: how each label was supported.",
  "- E: v1 to v2 relabelling: what changed and where (e.g. luteal-like Theca_Stroma -> Granulosa, false Oocyte calls removed).",
  "", "All plotted values are written as TSV next to each figure."), file.path(out_dir, "README.md"))
message("Annotation evidence figures written to ", out_dir)
