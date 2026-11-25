#!/usr/bin/env Rscript
# Combined annotation-evidence figure (panels a-d) from the tables written by build_annotation_evidence_v2.R:
#   a marker dot plot (pooled)  b signature scores per dataset x class
#   c benchmark against author labels  d cross-dataset cluster similarity
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
v2_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
ev_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce/step_02_preprocessing_filtering/figures/annotation_evidence_v2")

classes <- c("Granulosa", "Theca_Stroma", "Endothelial", "Pericyte_Smooth_Muscle", "Epithelial",
             "Immune_Myeloid", "T_NK", "B_Plasma", "Oocyte_Germline")
class_lab <- c(Granulosa = "Granulosa", Theca_Stroma = "Theca/stroma", Endothelial = "Endothelial",
               Pericyte_Smooth_Muscle = "Pericyte/SMC", Epithelial = "Epithelial", Immune_Myeloid = "Myeloid",
               T_NK = "T/NK", B_Plasma = "B/plasma", Oocyte_Germline = "Oocyte", Erythroid = "Erythroid", Unknown = "Unknown")
class_colors <- c(Granulosa = "#E7298A", Theca_Stroma = "#B07AE8", Endothelial = "#1FA5C9",
                  Pericyte_Smooth_Muscle = "#6A3D9A", Epithelial = "#7FB800", Immune_Myeloid = "#E06A00",
                  T_NK = "#0067A5", B_Plasma = "#1B9E77", Oocyte_Germline = "#F2C14E", Erythroid = "#C0392B", Unknown = "grey70")
sig_lab <- c(granulosa = "Granulosa", luteal = "Luteal", theca = "Theca", pericyte = "Pericyte", smooth_muscle = "SMC",
             stroma = "Stroma", endothelial = "Endothelial", epithelial = "Epithelial", oocyte = "Oocyte",
             myeloid = "Myeloid", neutrophil = "Neutrophil", t_nk = "T/NK", b_plasma = "B/plasma", erythroid = "Erythroid")
ds_short <- function(x) {
  x <- sub("^(20[0-9]{2})_[A-Za-z]+_([A-Za-z]+)_.*$", "\\2 \\1", x)
  sub("^(20[0-9]{2})_Kim_PRJNA863443.*$", "Kim 2025", x)
}
base <- theme_bw(base_size = 9) + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
                                        plot.title = element_text(face = "bold", size = 10))
div_fill <- function(name, lim) scale_fill_gradient2(low = "#3B4CC0", mid = "grey96", high = "#B40426", midpoint = 0,
                                                    limits = lim, oob = scales::squish, name = name)

# ---- a: marker dot plot ----------------------------------------------------------------------------------------
pool <- fread(file.path(ev_dir, "A_marker_dotplot_pooled.tsv"))
panel <- fread(file.path(v2_dir, "01_marker_panel_v2.tsv"))
pool[, gene := factor(gene, levels = unique(panel$gene))]
pool[, signature := factor(sig_lab[signature], levels = sig_lab[unique(panel$signature)])]
pool[, cell_class := factor(class_lab[cell_class], levels = rev(class_lab[classes]))]
pa <- ggplot(pool, aes(gene, cell_class, size = pct, colour = z)) + geom_point() +
  facet_grid(~signature, scales = "free_x", space = "free_x") +
  scale_colour_gradient2(low = "#3B4CC0", mid = "grey90", high = "#B40426", midpoint = 0, limits = c(-1.5, 2.5),
                         oob = scales::squish, name = "Mean z") +
  scale_size_area(max_size = 3.2, name = "Fraction\nexpressing", labels = scales::percent) + base +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.5, face = "italic"),
        strip.text.x = element_text(size = 7), panel.spacing.x = unit(1.5, "pt"), legend.key.height = unit(8, "pt")) +
  labs(x = NULL, y = NULL, title = "Canonical marker expression in v2 cell classes (8 datasets pooled)")

# ---- b: signature score per dataset x class ----------------------------------------------------------------------
sd_ <- fread(file.path(ev_dir, "A_signature_score_by_dataset_class.tsv"))
sd_[, signature := factor(sig_lab[signature], levels = sig_lab[unique(panel$signature)])]
sd_[, cell_class := factor(class_lab[cell_class], levels = class_lab[classes])]
sd_[, dataset := ds_short(dataset_id)]
pb <- ggplot(sd_, aes(signature, dataset, fill = score)) + geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(cell_class ~ ., scales = "free_y", space = "free_y", switch = "y") + div_fill("Signature\nscore (z)", c(-1.5, 2.5)) +
  base + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5), axis.text.y = element_text(size = 5.5),
               strip.text.y.left = element_text(angle = 0, size = 7), strip.placement = "outside",
               panel.spacing.y = unit(1, "pt"), legend.key.height = unit(8, "pt")) +
  labs(x = NULL, y = NULL, title = "Signature scores per dataset and class")

# ---- c: benchmark ------------------------------------------------------------------------------------------------
cf <- fread(file.path(v2_dir, "03_benchmark", "benchmark_confusion.tsv"))
cf[, frac := N / sum(N), by = .(query, method, truth)]
cf[, q := fifelse(grepl("Morris", query), "Morris", "Isola")]
cf[, other := fifelse(q == "Morris", "Isola", "Morris")]
cf[, panel := fifelse(method == "reference", paste0(q, " cells ← ", other, " reference"),
                      paste0(q, " cells ← marker scoring"))]
cf[, panel := factor(panel, levels = c("Morris cells ← Isola reference", "Isola cells ← Morris reference",
                                       "Morris cells ← marker scoring", "Isola cells ← marker scoring"))]
lv <- c(classes, "Unknown")
cf[, truth := factor(class_lab[truth], levels = rev(class_lab[lv]))][, pred := factor(class_lab[pred], levels = class_lab[lv])]
pc <- ggplot(cf, aes(pred, truth, fill = frac)) + geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(frac >= 0.05, sprintf("%.2f", frac), "")), size = 1.9) +
  facet_wrap(~panel, ncol = 2) + scale_fill_gradient(low = "white", high = "#2166AC", limits = c(0, 1), name = "Row\nfraction") +
  base + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5), axis.text.y = element_text(size = 6.5),
               strip.text = element_text(size = 7.5), legend.key.height = unit(8, "pt"), aspect.ratio = 1) +
  labs(x = "Predicted class", y = "Author class", title = "Accuracy against author labels (each atlas predicted without its own labels)",
       caption = paste("Isola reference lacks a pericyte class (ambiguous 'Stroma C' excluded), so Morris pericyte/SMC map to theca/stroma;",
                       "\nthe Morris reference pulls some Isola endothelial/epithelial/oocyte cells towards luteal granulosa."))

# ---- d: cross-dataset cluster similarity ------------------------------------------------------------------------
prof <- readRDS(file.path(v2_dir, "04_consensus", "cluster_profiles_refgenes.rds"))
P <- do.call(cbind, lapply(prof, function(M) M - rowMeans(M)))
P <- P[stats::complete.cases(P), ]
C <- cor(P, method = "spearman")
cl <- fread(file.path(v2_dir, "04_consensus", "cluster_consensus_v2.tsv"), colClasses = list(character = "cluster"))
cl[, key := paste(dataset_id, cluster, sep = "::")]
cl <- cl[key %in% colnames(C) & n_cells >= 50]
cl[, cluster_label := factor(cluster_label, levels = c(classes, "Erythroid", "Unknown"))]
ord <- unlist(lapply(split(cl$key, cl$cluster_label, drop = TRUE), function(k)
  if (length(k) > 2) k[hclust(as.dist(1 - C[k, k]))$order] else k))
cd <- as.data.table(as.table(C[ord, ord])); setnames(cd, c("a", "b", "rho"))
cd[, a := factor(a, levels = ord)][, b := factor(b, levels = rev(ord))]
ann <- cl[match(ord, key), .(key = factor(key, levels = ord), lab = factor(class_lab[as.character(cluster_label)], levels = class_lab))]
# single plot: square heatmap with the label bar drawn directly underneath (constant-fill rects avoid a 2nd fill scale)
n <- length(ord)
cd[, x := as.integer(factor(a, levels = ord))][, y := n + 1L - as.integer(factor(as.character(b), levels = ord))]
ann[, x := seq_len(.N)]
runs <- ann[, .(xmin = min(x) - 0.5, xmax = max(x) + 0.5), by = lab]
bar_h <- n * 0.035; gap <- n * 0.01
pd <- ggplot() + geom_raster(data = cd, aes(x, y, fill = rho)) + div_fill("Spearman", c(-1, 1))
for (k in seq_len(nrow(runs))) {
  pd <- pd + annotate("rect", xmin = runs$xmin[k], xmax = runs$xmax[k], ymin = 0.5 - gap - bar_h, ymax = 0.5 - gap,
                      fill = setNames(class_colors, class_lab[names(class_colors)])[as.character(runs$lab[k])])
}
pd <- pd + geom_point(data = ann, aes(x, 1, colour = lab), alpha = 0) +
  scale_colour_manual(values = setNames(class_colors, class_lab[names(class_colors)]), name = "v2 label", drop = TRUE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 4), nrow = 2, order = 2),
         fill = guide_colourbar(order = 1, barwidth = unit(70, "pt"), barheight = unit(8, "pt"), title.vjust = 0.8)) +
  coord_fixed(expand = FALSE, ylim = c(0.5 - gap - bar_h, n + 0.5)) + theme_void(base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10, margin = margin(l = 22, b = 4)),
        legend.position = "bottom", legend.box = "vertical", legend.key.size = unit(8, "pt")) +
  labs(title = sprintf("Cross-dataset similarity of %d clusters, ordered by v2 label", n))

# ---- assemble (explicit grid placement: c and d square and bottom-aligned) ---------------------------------------
W <- 21; H <- 16
draw_fig <- function() {
  grid::grid.newpage()
  vp <- function(x, y, w, h) grid::viewport(x = grid::unit(x, "in"), y = grid::unit(y, "in"), width = grid::unit(w, "in"),
                                            height = grid::unit(h, "in"), just = c("left", "bottom"))
  tag <- function(lbl, x, y) grid::grid.text(lbl, x = grid::unit(x, "in"), y = grid::unit(y, "in"), just = c("left", "top"),
                                            gp = grid::gpar(fontsize = 15, fontface = "bold"))
  grid::grid.text("Evidence supporting the v2 cell-class annotation", x = grid::unit(0.2, "in"), y = grid::unit(H - 0.15, "in"),
                  just = c("left", "top"), gp = grid::gpar(fontsize = 15, fontface = "bold"))
  # a: full width, top
  print(pa, vp = vp(0.35, H - 5.55, W - 0.45, 5.0)); tag("a", 0.15, H - 0.5)
  # bottom row, common baseline y0
  y0 <- 0.15
  print(pb, vp = vp(0.35, y0, 5.0, 9.9)); tag("b", 0.15, y0 + 9.95)
  sq <- 8.0   # side of the square blocks for c and d
  print(pc, vp = vp(5.6, y0, sq + 0.9, sq + 0.55)); tag("c", 5.5, y0 + sq + 0.6)
  print(pd, vp = vp(14.4, y0, W - 14.5, sq + 0.55)); tag("d", 14.3, y0 + sq + 0.6)
}
grDevices::png(file.path(ev_dir, "annotation_evidence_combined.png"), width = W, height = H, units = "in", res = 220)
draw_fig(); invisible(grDevices::dev.off())
grDevices::cairo_pdf(file.path(ev_dir, "annotation_evidence_combined.pdf"), width = W, height = H)
draw_fig(); invisible(grDevices::dev.off())
message("written: ", file.path(ev_dir, "annotation_evidence_combined.{png,pdf}"))
