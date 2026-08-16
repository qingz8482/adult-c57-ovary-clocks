#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
results_dir <- file.path(project_dir, "results")
selected_dir <- file.path(results_dir, "selected_genes")

source_path <- file.path(
  selected_dir,
  "integrated_bulk_pseudobulk_elasticnet_hvg_nonzero_gene_age_correlations.tsv"
)
output_stem <- file.path(
  selected_dir,
  "integrated_bulk_pseudobulk_elasticnet_hvg_nonzero_gene_age_correlation_dotplot"
)

if (!file.exists(source_path)) stop("Correlation source table is missing: ", source_path, call. = FALSE)

dt <- fread(source_path)
required <- c(
  "gene", "calibrated_gene_coefficient", "model_direction", "n_samples",
  "pearson_r", "pearson_fdr_bh"
)
missing <- setdiff(required, names(dt))
if (length(missing) > 0) {
  stop("Correlation source table is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}
if (nrow(dt) != 66L || uniqueN(dt$gene) != 66L) {
  stop("Expected 66 unique non-zero deployment genes; found ", nrow(dt), " rows and ", uniqueN(dt$gene), " genes.", call. = FALSE)
}
if (any(!is.finite(dt$pearson_r))) stop("One or more Pearson correlations are non-finite.", call. = FALSE)

dt[, significance := factor(
  fifelse(pearson_fdr_bh < 0.05, "BH FDR < 0.05", "BH FDR >= 0.05"),
  levels = c("BH FDR < 0.05", "BH FDR >= 0.05")
)]
dt[, model_direction := factor(
  model_direction,
  levels = c("older_predictive", "younger_predictive"),
  labels = c("Positive coefficient", "Negative coefficient")
)]
setorder(dt, pearson_r, gene)
dt[, gene := factor(gene, levels = gene)]

direction_colors <- c(
  "Positive coefficient" = "#B64342",
  "Negative coefficient" = "#3775BA"
)
significance_shapes <- c(
  "BH FDR < 0.05" = 16,
  "BH FDR >= 0.05" = 1
)

p <- ggplot(dt, aes(y = gene, x = pearson_r)) +
  geom_segment(
    aes(x = 0, xend = pearson_r, yend = gene),
    linewidth = 0.25,
    colour = "#C9C9C9"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "#5F5F5F") +
  geom_point(
    aes(colour = model_direction, shape = significance),
    size = 2.0,
    stroke = 0.65
  ) +
  scale_colour_manual(values = direction_colors, drop = FALSE) +
  scale_shape_manual(values = significance_shapes, drop = FALSE) +
  scale_x_continuous(
    limits = c(-0.65, 0.65),
    breaks = seq(-0.6, 0.6, by = 0.2),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_discrete(expand = expansion(add = 0.55)) +
  labs(
    title = "ElasticNet genes show distinct marginal associations with chronological age",
    subtitle = "Pearson correlations across 117 control-like bulk RNA-seq and scRNA pseudobulk samples",
    x = "Pearson correlation with chronological age (r)",
    y = NULL,
    colour = "ElasticNet direction",
    shape = "Age association"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 16, size = 2.2)),
    shape = guide_legend(order = 2, override.aes = list(colour = "#333333", size = 2.2))
  ) +
  theme_classic(base_size = 6.5, base_family = "sans") +
  theme(
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_line(linewidth = 0.35, colour = "black"),
    axis.text.y = element_text(size = 5.5, colour = "#222222", margin = margin(r = 2)),
    axis.text.x = element_text(size = 6.2, colour = "#222222"),
    axis.title.x = element_text(size = 6.8, margin = margin(t = 5)),
    plot.title = element_text(size = 8, face = "bold", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 6.5, colour = "#4D4D4D", margin = margin(b = 5)),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "left",
    legend.title = element_text(size = 6.2, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.width = grid::unit(4, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    legend.spacing.y = grid::unit(0.5, "mm"),
    plot.margin = margin(5, 5, 4, 4)
  )

width_mm <- 140
height_mm <- 190
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

grDevices::svg(
  paste0(output_stem, ".svg"),
  width = width_in,
  height = height_in,
  family = "sans",
  onefile = TRUE
)
print(p)
dev.off()

grDevices::cairo_pdf(
  paste0(output_stem, ".pdf"),
  width = width_in,
  height = height_in,
  family = "sans"
)
print(p)
dev.off()

ragg::agg_png(
  paste0(output_stem, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300,
  background = "white"
)
print(p)
dev.off()

ragg::agg_tiff(
  paste0(output_stem, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  background = "white",
  compression = "lzw"
)
print(p)
dev.off()

message("Source: ", source_path)
message("Output stem: ", output_stem)
message("Genes plotted: ", nrow(dt))
