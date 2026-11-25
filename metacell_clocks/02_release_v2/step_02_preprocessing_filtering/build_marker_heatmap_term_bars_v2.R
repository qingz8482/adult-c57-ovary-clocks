#!/usr/bin/env Rscript
# Compact marker heatmap + one representative term per class (layout after the user's reference panel):
# left, block-diagonal heatmap of the 50 unique markers per class (row z-score of mean TMM logCPM);
# right, one bar per class aligned with its block, -log10 nominal p, term written inside the bar.
# Inputs: somatic_marker_expression_function_panel/ (built by build_somatic_marker_expression_function_panel.R).
# Displayed terms are set in display_terms.tsv (edit to swap); terms not FDR-significant are drawn faded and marked.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
step_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce/step_02_preprocessing_filtering")
in_dir <- file.path(step_dir, "somatic_marker_expression_function_panel")
out_dir <- in_dir

classes <- c("Granulosa", "Theca_Stroma", "Epithelial", "Endothelial", "Pericyte_Smooth_Muscle",
             "Immune_Myeloid", "T_NK", "B_Plasma")
short <- c(Granulosa = "GC", Theca_Stroma = "TC/SC", Epithelial = "EPI", Endothelial = "EC",
           Pericyte_Smooth_Muscle = "PC/SMC", Immune_Myeloid = "MYE", T_NK = "T/NK", B_Plasma = "B")
cols <- c(Granulosa = "#2E8B3E", Theca_Stroma = "#3DAA85", Epithelial = "#E3C238", Endothelial = "#E8857F",
          Pericyte_Smooth_Muscle = "#5E5AA0", Immune_Myeloid = "#E0592A", T_NK = "#B3282D", B_Plasma = "#2A72B5")

# ---- displayed marker genes (written inside the bars) ------------------------------------------------------
# default: top 5 marker genes per class by marker rank, skipping unannotated predicted genes; edit display_genes.tsv to change
mk <- fread(file.path(in_dir, "top50_unique_somatic_markers_by_celltype.tsv"))
gene_path <- file.path(out_dir, "display_genes.tsv")
if (!file.exists(gene_path)) {
  def <- mk[!grepl("^Gm[0-9]|^Olfr|Rik$|^LOC", gene)][order(marker_rank), .(genes = paste(head(gene, 5), collapse = ", ")),
                                                        by = harmonized_cell_class]
  fwrite(def[match(classes, harmonized_cell_class)], gene_path, sep = "\t")
}
terms <- fread(gene_path)
# bar length: median log2 fold change of the class's 50 markers versus the other classes
strength <- mk[, .(median_logfc = median(logFC), n_markers = .N), by = harmonized_cell_class]
terms <- merge(terms, strength, by = "harmonized_cell_class")
terms[, label := genes]

# ---- heatmap data (blocks in class order, genes by marker rank) ----------------------------------------------
heat <- fread(file.path(in_dir, "marker_heatmap_source_data.tsv"))
heat <- heat[marker_celltype %in% classes & expression_celltype %in% classes]
gene_order <- unique(heat[order(match(marker_celltype, classes), marker_rank), .(gene, marker_celltype)])
gene_order[, row := rev(seq_len(.N))]                       # first block at the top
heat <- merge(heat, gene_order[, .(gene, row)], by = "gene")
heat[, col := match(expression_celltype, classes)]
blocks <- gene_order[, .(ymin = min(row) - 0.5, ymax = max(row) + 0.5, center = mean(row)), by = marker_celltype]
terms <- merge(terms, blocks, by.x = "harmonized_cell_class", by.y = "marker_celltype")
n_rows <- nrow(gene_order)

# geom_tile (vector) instead of geom_raster so PDF viewers do not interpolate/blur; white separators between blocks
ph <- ggplot(heat, aes(col, row, fill = pmin(pmax(z_score, -2), 3))) + geom_tile(width = 1, height = 1, colour = NA) +
  geom_hline(yintercept = blocks$ymin[blocks$ymin > 0.5], colour = "white", linewidth = 1.1) +
  geom_vline(xintercept = seq(1.5, length(classes) - 0.5, by = 1), colour = "white", linewidth = 1.1) +
  scale_fill_gradientn(colours = c("#2C7BB6", "#91BFDB", "#F7F7F7", "#FC8D59", "#D7301F"),
                       values = scales::rescale(c(-2, -1, 0, 1.2, 3)), limits = c(-2, 3), name = "Row z-score",
                       breaks = -2:3, guide = guide_colourbar(direction = "horizontal", title.position = "left",
                                                              barwidth = unit(110, "pt"), barheight = unit(9, "pt"))) +
  scale_x_continuous(expand = c(0, 0), limits = c(0.5, length(classes) + 0.5)) + scale_y_continuous(expand = c(0, 0), limits = c(0.5, n_rows + 0.5)) +
  theme_void(base_size = 11) + theme(legend.position = "bottom", legend.title = element_text(size = 12))

pdots <- ggplot(data.table(col = seq_along(classes), cls = classes), aes(col, 1, colour = cls)) +
  geom_point(size = 6) + scale_colour_manual(values = cols, guide = "none") +
  scale_y_continuous(limits = c(0.6, 1.4), expand = c(0, 0)) + coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = seq_along(classes), labels = short[classes], limits = c(0.5, length(classes) + 0.5), expand = c(0, 0)) +
  theme_void(base_size = 12) + theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 13, margin = margin(t = 2)), plot.margin = margin(t = 4))

xmax <- ceiling(max(terms$median_logfc) * 1.15)
pbar <- ggplot(terms) +
  geom_rect(aes(xmin = 0, xmax = median_logfc, ymin = center - 18, ymax = center + 18, fill = harmonized_cell_class)) +
  geom_text(aes(x = xmax * 0.012, y = center, label = label), hjust = 0, size = 4.4, fontface = "italic") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0), name = expression("Median "*log[2]*"FC of top-50 markers vs other classes")) +
  scale_y_continuous(limits = c(0.5, n_rows + 0.5), expand = c(0, 0), breaks = terms$center) +
  theme_classic(base_size = 12) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.line.y = element_line(),
        axis.ticks.y = element_line(), plot.margin = margin(l = 8, r = 10))

foot <- "Genes shown: top-ranked unique markers of each class (from the 50 per class in the heatmap)"
fig <- (ph + pbar + pdots + guide_area() + plot_layout(ncol = 2, widths = c(1, 1.05), heights = c(1, 0.06), guides = "collect")) +
  plot_annotation(caption = foot, theme = theme(plot.caption = element_text(size = 10, hjust = 1)))
ggsave(file.path(out_dir, "marker_heatmap_with_gene_bars.png"), fig, width = 11, height = 6.6, dpi = 400)
ggsave(file.path(out_dir, "marker_heatmap_with_gene_bars.pdf"), fig, width = 11, height = 6.6, device = grDevices::cairo_pdf)
fwrite(terms[order(match(harmonized_cell_class, classes)), .(harmonized_cell_class, genes, median_logfc, n_markers)],
       file.path(out_dir, "marker_heatmap_with_gene_bars_genes.tsv"), sep = "\t")
message("written marker_heatmap_with_gene_bars.{png,pdf}")
