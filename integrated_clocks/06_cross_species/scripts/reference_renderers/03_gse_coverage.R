coverage_path<-file.path(gcompat,'projection/GSE202601_model_gene_coverage.tsv');gtex_coverage_path<-file.path(figdir,'source_data/panel_B_orthology_coverage.tsv');plot_dir<-figdir;source_dir<-file.path(figdir,'source_data/GSE202601')
x <- fread(coverage_path)
required_columns <- c(
  "model_family", "high_confidence_one_to_one_mapping",
  "human_gene_detected_after_filter", "absolute_expression_space_weight"
)
stopifnot(
  all(required_columns %in% names(x)),
  nrow(x) == 6000L,
  setequal(unique(x$model_family), c("elasticnet", "pca_glmnet"))
)

gse_coverage <- x[, .(
  cohort = "GSE202601",
  model_label = fifelse(
    model_family[[1]] == "elasticnet",
    "ElasticNet/HVG",
    "PCA/glmnet"
  ),
  n_total = .N,
  n_one2one = sum(high_confidence_one_to_one_mapping),
  n_human_detected = sum(human_gene_detected_after_filter),
  raw_gene_coverage = mean(high_confidence_one_to_one_mapping),
  human_detected_coverage = mean(human_gene_detected_after_filter),
  absolute_model_weight_coverage = {
    total_weight <- sum(absolute_expression_space_weight)
    if (total_weight > 0) {
      sum(absolute_expression_space_weight[human_gene_detected_after_filter]) /
        total_weight
    } else {
      NA_real_
    }
  }
), by = model_family]

gtex_source <- fread(gtex_coverage_path)
gtex_coverage <- gtex_source[
  gene_set %chin% c("elasticnet_hvg_3000", "pca_versionB_3000")
]
gtex_coverage[, `:=`(
  cohort = "GTEx",
  model_family = fifelse(
    gene_set == "elasticnet_hvg_3000", "elasticnet", "pca_glmnet"
  ),
  model_label = fifelse(
    gene_set == "elasticnet_hvg_3000", "ElasticNet/HVG", "PCA/glmnet"
  )
)]
gtex_coverage <- gtex_coverage[, .(
  cohort, model_family, model_label, n_total, n_one2one,
  n_human_detected, raw_gene_coverage, human_detected_coverage,
  absolute_model_weight_coverage
)]

coverage <- rbindlist(list(gtex_coverage, gse_coverage), use.names = TRUE)
coverage[, row_label := paste(cohort, model_label, sep = " · ")]
row_order <- c(
  "GTEx · ElasticNet/HVG", "GTEx · PCA/glmnet",
  "GSE202601 · ElasticNet/HVG", "GSE202601 · PCA/glmnet"
)
coverage[, row_label := factor(row_label, levels = rev(row_order))]

coverage_long <- melt(
  coverage,
  id.vars = c(
    "cohort", "model_family", "model_label", "row_label", "n_total",
    "n_one2one", "n_human_detected"
  ),
  measure.vars = c(
    "raw_gene_coverage", "human_detected_coverage",
    "absolute_model_weight_coverage"
  ),
  variable.name = "metric",
  value.name = "coverage"
)
coverage_long[, metric := factor(
  metric,
  levels = c(
    "raw_gene_coverage", "human_detected_coverage",
    "absolute_model_weight_coverage"
  ),
  labels = c(
    "High-confidence 1:1 orthologue",
    "Detected in cohort",
    "Absolute model-weight coverage"
  )
)]

palette <- c(
  "High-confidence 1:1 orthologue" = "#777777",
  "Detected in cohort" = "#4E9B88",
  "Absolute model-weight coverage" = "#222222"
)
shapes <- c(
  "High-confidence 1:1 orthologue" = 16,
  "Detected in cohort" = 17,
  "Absolute model-weight coverage" = 15
)

p <- ggplot(
  coverage_long,
  aes(x = coverage, y = row_label, shape = metric, colour = metric)
) +
  geom_vline(
    xintercept = c(0.6, 0.8),
    linetype = c(3, 2),
    colour = "#AAAAAA",
    linewidth = 0.35
  ) +
  geom_hline(yintercept = 2.5, colour = "#D0D0D0", linewidth = 0.35) +
  geom_point(size = 2) +
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_shape_manual(values = shapes, drop = FALSE) +
  labs(
    x = "Coverage",
    y = NULL,
    colour = NULL,
    shape = NULL,
    title = "Human model-gene coverage"
  ) +
  theme_classic(base_size = 7.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 8.5),
    legend.position = "bottom",
    legend.key.height = grid::unit(3, "mm"),
    legend.key.width = grid::unit(4, "mm"),
    legend.text = element_text(size = 6.2),
    legend.box = "vertical",
    plot.margin = margin(4, 4, 4, 4)
  ) +
  guides(
    colour = guide_legend(nrow = 3, byrow = TRUE),
    shape = guide_legend(nrow = 3, byrow = TRUE)
  )

stem <- file.path(plot_dir, "03_GTEx_GSE202601_gene_coverage")
width_mm <- 90
height_mm <- 85
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

grDevices::svg(
  paste0(stem, ".svg"), width = width_in, height = height_in,
  bg = "white", family = "sans"
)
print(p)
dev.off()

grDevices::cairo_pdf(
  paste0(stem, ".pdf"), width = width_in, height = height_in,
  bg = "white", family = "sans"
)
print(p)
dev.off()

ragg::agg_png(
  paste0(stem, ".png"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white"
)
print(p)
dev.off()

ragg::agg_tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p)
dev.off()

fwrite(
  coverage,
  file.path(source_dir, "GTEx_GSE202601_gene_coverage_source_data.tsv"),
  sep = "\t"
)


z<-copy(coverage_long[cohort=='GSE202601']);z[,row_label:=factor(ifelse(model_family=='elasticnet','ElasticNet/HVG (3,000 genes)','PCA/glmnet (3,000 genes)'),levels=c('PCA/glmnet (3,000 genes)','ElasticNet/HVG (3,000 genes)'))]
z[,metric:=factor(as.character(metric),levels=levels(coverage_long$metric),labels=c('High-confidence 1:1 orthologue','Detected in GSE202601','Absolute model-weight coverage'))]
p2<-p %+% z;p2$layers<-p2$layers[-2];p2<-p2+labs(title='GSE202601 orthology and detectability')+scale_colour_manual(values=setNames(unname(palette),levels(z$metric)))+scale_shape_manual(values=setNames(unname(shapes),levels(z$metric)))
save_bundle(p2,file.path(figdir,'03_GSE202601_gene_coverage'),90,80)
write_tsv(z,file.path(source_dir,'GSE202601_gene_coverage_source_data.tsv'))
