gtex_path<-file.path(figdir,'source_data/panel_F_model_ood.tsv');gse_path<-file.path(gcompat,'projection/GSE202601_two_integrated_clock_projections.tsv');plot_dir<-figdir;source_dir<-file.path(figdir,'source_data/GSE202601')
gtex <- fread(gtex_path)
stopifnot(
  nrow(gtex) == 193L,
  uniqueN(gtex$donor_id) == 193L
)
gtex_plot <- gtex[, .(
  cohort = "GTEx",
  sample_id = donor_id,
  age_label = age_bracket,
  elasticnet_score = cross_species_mouse_ovarian_clock_score_elasticnet,
  pca_glmnet_score = cross_species_mouse_ovarian_clock_score_pca_glmnet,
  elasticnet_ood_distance = projection_ood_distance_elasticnet,
  pca_glmnet_ood_distance = projection_ood_distance_pca_glmnet
)]

gse_long <- fread(gse_path)
stopifnot(
  nrow(gse_long) == 16L,
  uniqueN(gse_long$projection_sample_id) == 8L,
  setequal(unique(gse_long$model_family), c("elasticnet", "pca_glmnet"))
)
gse_wide <- dcast(
  gse_long,
  projection_sample_id + sample_id + age_year + age_group ~ model_family,
  value.var = c(
    "cross_species_mouse_ovarian_clock_score", "projection_ood_distance"
  )
)
gse_plot <- gse_wide[, .(
  cohort = "GSE202601",
  sample_id = projection_sample_id,
  age_label = sprintf("%s years; %s", age_year, age_group),
  elasticnet_score = cross_species_mouse_ovarian_clock_score_elasticnet,
  pca_glmnet_score = cross_species_mouse_ovarian_clock_score_pca_glmnet,
  elasticnet_ood_distance = projection_ood_distance_elasticnet,
  pca_glmnet_ood_distance = projection_ood_distance_pca_glmnet
)]

plot_data <- rbindlist(list(gtex_plot, gse_plot), use.names = TRUE)
plot_data[, cohort := factor(cohort, levels = c("GTEx", "GSE202601"))]
stopifnot(
  all(is.finite(plot_data$elasticnet_score)),
  all(is.finite(plot_data$pca_glmnet_score)),
  all(is.finite(plot_data$elasticnet_ood_distance))
)

statistics <- plot_data[, {
  pearson <- cor.test(elasticnet_score, pca_glmnet_score, method = "pearson")
  spearman <- suppressWarnings(cor.test(
    elasticnet_score, pca_glmnet_score, method = "spearman", exact = FALSE
  ))
  .(
    n_samples = .N,
    pearson_r = unname(pearson$estimate),
    pearson_ci_lower = pearson$conf.int[[1]],
    pearson_ci_upper = pearson$conf.int[[2]],
    pearson_p = pearson$p.value,
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value
  )
}, by = cohort]

pooled_pearson <- cor.test(
  plot_data$elasticnet_score,
  plot_data$pca_glmnet_score,
  method = "pearson"
)
pooled_spearman <- suppressWarnings(cor.test(
  plot_data$elasticnet_score,
  plot_data$pca_glmnet_score,
  method = "spearman",
  exact = FALSE
))
pooled_statistics <- data.table(
  cohort = "All samples",
  n_samples = nrow(plot_data),
  pearson_r = unname(pooled_pearson$estimate),
  pearson_ci_lower = pooled_pearson$conf.int[[1]],
  pearson_ci_upper = pooled_pearson$conf.int[[2]],
  pearson_p = pooled_pearson$p.value,
  spearman_rho = unname(pooled_spearman$estimate),
  spearman_p = pooled_spearman$p.value
)
statistics[, cohort := as.character(cohort)]
statistics <- rbindlist(list(statistics, pooled_statistics), use.names = TRUE)
statistics[, annotation := sprintf(
  "%s: Pearson r = %.02f, P = %.3f (n = %d)",
  cohort, pearson_r, pearson_p, n_samples
)]

score_range <- range(c(plot_data$elasticnet_score, plot_data$pca_glmnet_score))
score_padding <- max(diff(score_range) * 0.06, 0.15)
score_limits <- score_range + c(-score_padding, score_padding)
annotation_data <- data.table(
  annotation_x = score_limits[[1]] + 0.04 * diff(score_limits),
  annotation_y = score_limits[[2]] - 0.04 * diff(score_limits),
  annotation = statistics[cohort == "All samples", annotation]
)

p <- ggplot(
  plot_data,
  aes(
    x = elasticnet_score,
    y = pca_glmnet_score,
    colour = elasticnet_ood_distance,
    shape = cohort,
    size = cohort
  )
) +
  geom_hline(yintercept = 0, colour = "#CCCCCC", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#CCCCCC", linewidth = 0.35) +
  geom_point(alpha = 0.72, stroke = 0.45) +
  geom_text(
    data = annotation_data,
    aes(x = annotation_x, y = annotation_y, label = annotation),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.35,
    family = "sans",
    lineheight = 1.08,
    colour = "#222222"
  ) +
  scale_x_continuous(limits = score_limits) +
  scale_y_continuous(limits = score_limits) +
  scale_shape_manual(
    values = c("GTEx" = 16, "GSE202601" = 17),
    name = "Cohort"
  ) +
  scale_size_manual(
    values = c("GTEx" = 1.05, "GSE202601" = 2.0),
    guide = "none"
  ) +
  scale_colour_viridis_c(
    option = "C",
    name = "ElasticNet\nOOD distance"
  ) +
  coord_equal() +
  guides(
    colour = guide_colourbar(barheight = grid::unit(22, "mm")),
    shape = guide_legend(override.aes = list(size = 1.8, colour = "#555555"))
  ) +
  labs(
    x = "ElasticNet standardized\nmouse-clock score",
    y = "PCA/glmnet standardized\nmouse-clock score",
    title = "Cross-model agreement and OOD"
  ) +
  theme_classic(base_size = 7.5, base_family = "sans") +
  theme(
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.35),
    plot.title = element_text(face = "bold", size = 8.5),
    legend.title = element_text(size = 6.5),
    legend.text = element_text(size = 6.2),
    legend.key.height = grid::unit(4, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )

stem <- file.path(plot_dir, "04_GTEx_GSE202601_cross_model_agreement")
width_mm <- 100
height_mm <- 78
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

source_data_path <- file.path(
  source_dir, "GTEx_GSE202601_cross_model_agreement_source_data.tsv"
)
statistics_path <- file.path(
  source_dir, "GTEx_GSE202601_cross_model_agreement_statistics.tsv"
)
fwrite(plot_data, source_data_path, sep = "\t")
fwrite(
  statistics,
  statistics_path,
  sep = "\t"
)

