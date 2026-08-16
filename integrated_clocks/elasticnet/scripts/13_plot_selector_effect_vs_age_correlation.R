#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- file_arg[grepl("^--file=", file_arg)]
script_file <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
project_dir <- normalizePath(file.path(script_dir, "..", "..", "shared"), mustWork = TRUE)
selected_dir <- file.path(project_dir, "results", "selected_genes")

elastic_path <- file.path(
  selected_dir,
  "integrated_bulk_pseudobulk_elasticnet_hvg_nonzero_gene_age_correlations.tsv"
)
version_b_path <- file.path(
  selected_dir,
  "integrated_bulk_pseudobulk_pca_glmnet_versionB_selected_gene_age_correlations.tsv"
)

for (path in c(elastic_path, version_b_path)) {
  if (!file.exists(path)) stop("Required source table is missing: ", path, call. = FALSE)
}

elastic <- fread(elastic_path)
version_b <- fread(version_b_path)

elastic_required <- c(
  "gene", "raw_glmnet_gene_coefficient", "pearson_r",
  "pearson_direction_matches_coefficient", "n_samples"
)
version_b_required <- c(
  "gene", "age_beta", "pearson_r", "pearson_direction_matches_age_beta",
  "feature_score", "n_samples"
)
elastic_missing <- setdiff(elastic_required, names(elastic))
version_b_missing <- setdiff(version_b_required, names(version_b))
if (length(elastic_missing) > 0) {
  stop("ElasticNet table is missing columns: ", paste(elastic_missing, collapse = ", "), call. = FALSE)
}
if (length(version_b_missing) > 0) {
  stop("Version B table is missing columns: ", paste(version_b_missing, collapse = ", "), call. = FALSE)
}
if (nrow(elastic) != 66L || uniqueN(elastic$gene) != 66L) {
  stop("Expected 66 unique non-zero ElasticNet genes.", call. = FALSE)
}
if (nrow(version_b) != 2000L || uniqueN(version_b$gene) != 2000L) {
  stop("Expected 2,000 unique Version B genes.", call. = FALSE)
}
if (any(elastic$n_samples != 117L) || any(version_b$n_samples != 117L)) {
  stop("All plotted correlations must use the 117-sample deployment reference.", call. = FALSE)
}

elastic[, concordance := factor(
  fifelse(pearson_direction_matches_coefficient, "Concordant", "Discordant"),
  levels = c("Concordant", "Discordant")
)]
version_b[, concordance := factor(
  fifelse(pearson_direction_matches_age_beta, "Concordant", "Discordant"),
  levels = c("Concordant", "Discordant")
)]

elastic_labels <- elastic[order(-abs(raw_glmnet_gene_coefficient))][1:10]
version_b_labels <- version_b[order(-feature_score)][1:10]

elastic_relation_r <- cor(
  elastic$raw_glmnet_gene_coefficient,
  elastic$pearson_r,
  method = "pearson"
)
version_b_relation_r <- cor(version_b$age_beta, version_b$pearson_r, method = "pearson")

concordance_colors <- c(Concordant = "#3775BA", Discordant = "#D95F02")
concordance_shapes <- c(Concordant = 16, Discordant = 17)

theme_scatter <- theme_classic(base_size = 7, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(size = 6.4, colour = "#222222"),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 7.5, face = "bold", margin = margin(b = 3)),
    legend.position = "top",
    legend.title = element_text(size = 6.4, face = "bold"),
    legend.text = element_text(size = 6.2),
    legend.key.width = grid::unit(4, "mm"),
    legend.key.height = grid::unit(3, "mm"),
    panel.grid = element_blank(),
    plot.margin = margin(4, 5, 4, 4)
  )

p_elastic <- ggplot(
  elastic,
  aes(x = raw_glmnet_gene_coefficient, y = pearson_r)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55, colour = "#555555") +
  geom_point(
    aes(colour = concordance, shape = concordance),
    size = 2.0,
    alpha = 0.85,
    stroke = 0.5
  ) +
  geom_text_repel(
    data = elastic_labels,
    aes(label = gene),
    size = 2.0,
    family = "sans",
    colour = "#222222",
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 101,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = sprintf("Across-gene r = %.2f", elastic_relation_r),
    hjust = -0.05,
    vjust = 1.4,
    size = 2.2,
    family = "sans"
  ) +
  scale_colour_manual(values = concordance_colors, drop = FALSE) +
  scale_shape_manual(values = concordance_shapes, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(-4, 4, by = 2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(-0.6, 0.6, by = 0.3),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(xlim = c(-5.4, 5.4), ylim = c(-0.70, 0.70), clip = "off") +
  labs(
    title = "ElasticNet non-zero HVGs (n = 66)",
    x = "Raw ElasticNet coefficient (weeks per SD logCPM)",
    y = "Pearson correlation with age (r)",
    colour = "Effect-correlation direction",
    shape = "Effect-correlation direction"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 2.2, alpha = 1)),
    shape = "none"
  ) +
  theme_scatter

p_version_b <- ggplot(version_b, aes(x = age_beta, y = pearson_r)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "#777777") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.55, colour = "#555555") +
  geom_point(
    aes(colour = concordance, shape = concordance),
    size = 0.9,
    alpha = 0.50,
    stroke = 0.35
  ) +
  geom_text_repel(
    data = version_b_labels,
    aes(label = gene),
    size = 2.0,
    family = "sans",
    colour = "#222222",
    box.padding = 0.25,
    point.padding = 0.10,
    min.segment.length = 0,
    segment.size = 0.25,
    max.overlaps = Inf,
    seed = 102,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = sprintf("Across-gene r = %.2f", version_b_relation_r),
    hjust = -0.05,
    vjust = 1.4,
    size = 2.2,
    family = "sans"
  ) +
  scale_colour_manual(values = concordance_colors, drop = FALSE) +
  scale_shape_manual(values = concordance_shapes, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(-0.02, 0.02, by = 0.01),
    labels = function(x) sprintf("%.2f", x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq(-0.6, 0.6, by = 0.3),
    labels = function(x) sprintf("%.1f", x),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(xlim = c(-0.030, 0.027), ylim = c(-0.70, 0.70), clip = "off") +
  labs(
    title = "PCA/glmnet Version B genes (n = 2,000)",
    x = "Version B age beta (logCPM per week)",
    y = "Pearson correlation with age (r)",
    colour = "Effect-correlation direction",
    shape = "Effect-correlation direction"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 2.2, alpha = 1)),
    shape = "none"
  ) +
  theme_scatter

legend_data <- data.table(
  x = c(0.95, 2.05),
  y = 1,
  concordance = factor(c("Concordant", "Discordant"), levels = c("Concordant", "Discordant"))
)
legend_strip <- ggplot(legend_data, aes(x = x, y = y)) +
  annotate(
    "text", x = 0.05, y = 1, label = "Effect-correlation direction",
    hjust = 0, size = 2.3, fontface = "bold", family = "sans"
  ) +
  geom_point(aes(colour = concordance, shape = concordance), size = 2.2) +
  geom_text(aes(label = concordance), hjust = -0.2, size = 2.2, family = "sans") +
  scale_colour_manual(values = concordance_colors, drop = FALSE) +
  scale_shape_manual(values = concordance_shapes, drop = FALSE) +
  coord_cartesian(xlim = c(0, 2.8), ylim = c(0.8, 1.2), clip = "off") +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 4, 0, 4))

p_elastic_combined <- p_elastic +
  labs(tag = "a") +
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 8, face = "bold", family = "sans")
  )
p_version_b_combined <- p_version_b +
  labs(tag = "b") +
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 8, face = "bold", family = "sans")
  )
combined <- legend_strip / (p_elastic_combined + p_version_b_combined) +
  plot_layout(heights = c(0.10, 1), widths = 1)

save_plot_bundle <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  grDevices::svg(
    paste0(stem, ".svg"),
    width = width_in,
    height = height_in,
    family = "sans",
    onefile = TRUE
  )
  print(plot)
  dev.off()

  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "sans"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(stem, ".png"),
    width = width_in,
    height = height_in,
    units = "in",
    res = 300,
    background = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"),
    width = width_in,
    height = height_in,
    units = "in",
    res = 600,
    background = "white",
    compression = "lzw"
  )
  print(plot)
  dev.off()
}

elastic_stem <- file.path(selected_dir, "elasticnet_hvg_coefficient_vs_age_correlation_scatter")
version_b_stem <- file.path(selected_dir, "pca_glmnet_versionB_age_beta_vs_age_correlation_scatter")
combined_stem <- file.path(selected_dir, "selector_effect_vs_age_correlation_scatter_combined")

save_plot_bundle(p_elastic, elastic_stem, width_mm = 89, height_mm = 82)
save_plot_bundle(p_version_b, version_b_stem, width_mm = 89, height_mm = 82)
save_plot_bundle(combined, combined_stem, width_mm = 183, height_mm = 88)

message("ElasticNet source rows: ", nrow(elastic))
message("Version B source rows: ", nrow(version_b))
message("ElasticNet coefficient-correlation r: ", signif(elastic_relation_r, 4))
message("Version B beta-correlation r: ", signif(version_b_relation_r, 4))
message("Combined output stem: ", combined_stem)
