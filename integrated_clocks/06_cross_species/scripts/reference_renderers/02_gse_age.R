projection_path<-file.path(gcompat,'projection/GSE202601_two_integrated_clock_projections.tsv');plot_dir<-figdir;source_dir<-file.path(figdir,'source_data/GSE202601');dir.create(source_dir,recursive=TRUE,showWarnings=FALSE);set.seed(202601)
x <- fread(projection_path)
stopifnot(
  nrow(x) == 16L,
  uniqueN(x$projection_sample_id) == 8L,
  setequal(unique(x$model_family), c("elasticnet", "pca_glmnet")),
  all(is.finite(x$age_year)),
  all(is.finite(x$mouse_week_equivalent_output_audit_only)),
  all(is.finite(x$cross_species_mouse_ovarian_clock_score))
)

x[, model_label := factor(
  model_family,
  levels = c("elasticnet", "pca_glmnet"),
  labels = c("ElasticNet (HVG)", "PCA/glmnet (Version B)")
)]
x[, age_group_label := factor(
  age_group,
  levels = c("young", "reproductively_aged"),
  labels = c("Young", "Reproductively aged")
)]

make_stats <- function(y_var, slope_name) {
  out <- rbindlist(lapply(levels(x$model_label), function(label) {
    d <- x[model_label == label]
    fit <- lm(d[[y_var]] ~ d$age_year)
    co <- summary(fit)$coefficients["d$age_year", ]
    pearson_test <- cor.test(d$age_year, d[[y_var]], method = "pearson")
    spearman_test <- cor.test(
      d$age_year,
      d[[y_var]],
      method = "spearman",
      exact = TRUE
    )
    influence <- data.table(
      sample_id = d$sample_id,
      externally_studentized_residual = rstudent(fit),
      cooks_distance = cooks.distance(fit),
      leverage = hatvalues(fit)
    )
    influential <- influence[
      abs(externally_studentized_residual) > 3 |
        cooks_distance > 4 / nrow(d)
    ]
    data.table(
      model_label = factor(label, levels = levels(x$model_label)),
      n_donors = nrow(d),
      spearman_rho = unname(spearman_test$estimate),
      spearman_exact_p = spearman_test$p.value,
      pearson_r_sensitivity = unname(pearson_test$estimate),
      pearson_95ci_lower = pearson_test$conf.int[[1]],
      pearson_95ci_upper = pearson_test$conf.int[[2]],
      pearson_p_sensitivity = pearson_test$p.value,
      slope = unname(co[["Estimate"]]),
      slope_standard_error = unname(co[["Std. Error"]]),
      slope_p_sensitivity = unname(co[["Pr(>|t|)"]]),
      r_squared = summary(fit)$r.squared,
      age_shapiro_p = shapiro.test(d$age_year)$p.value,
      outcome_shapiro_p = shapiro.test(d[[y_var]])$p.value,
      residual_shapiro_p = shapiro.test(resid(fit))$p.value,
      max_abs_externally_studentized_residual = max(
        abs(influence$externally_studentized_residual)
      ),
      max_cooks_distance = max(influence$cooks_distance),
      max_leverage = max(influence$leverage),
      influential_sample_count = nrow(influential),
      influential_samples = if (nrow(influential)) {
        paste(influential$sample_id, collapse = ";")
      } else {
        NA_character_
      },
      annotation_x = min(d$age_year) + 0.6,
      annotation_y = max(d[[y_var]]),
      annotation = sprintf(
        "n = %d\nSpearman ρ = %.2f\nExact P = %.3f",
        nrow(d),
        unname(spearman_test$estimate),
        spearman_test$p.value
      )
    )
  }))
  setnames(out, "slope", slope_name)
  out
}

stats_mouse_week <- make_stats(
  "mouse_week_equivalent_output_audit_only",
  "slope_mouse_weeks_per_human_year"
)
stats_standardized <- make_stats(
  "cross_species_mouse_ovarian_clock_score",
  "slope_standardized_score_per_human_year"
)

palette <- c(
  "Young" = "#3E78B2",
  "Reproductively aged" = "#D98C3F"
)

make_plot <- function(
  y_var,
  stats,
  title,
  subtitle,
  y_label,
  free_y,
  repeat_y_axes = FALSE
) {
  ggplot(
    x,
    aes(
      x = age_year,
      y = .data[[y_var]],
      colour = age_group_label,
      shape = age_group_label
    )
  ) +
    geom_smooth(
      data = x,
      aes(
        x = age_year,
        y = .data[[y_var]],
        group = model_label
      ),
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      inherit.aes = FALSE,
      colour = "#4A4A4A",
      fill = "#BFC3C7",
      alpha = 0.28,
      linewidth = 0.55,
      show.legend = FALSE
    ) +
    geom_point(size = 2.4, stroke = 0.55, alpha = 0.95) +
    geom_text_repel(
      aes(label = sample_id),
      size = 2.25,
      family = "sans",
      min.segment.length = 0,
      segment.colour = "#999999",
      segment.size = 0.25,
      box.padding = 0.25,
      point.padding = 0.18,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    geom_text(
      data = stats,
      aes(x = annotation_x, y = annotation_y, label = annotation),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.35,
      family = "sans",
      lineheight = 1.12,
      colour = "#222222"
    ) +
    facet_wrap(
      ~model_label,
      nrow = 1,
      scales = if (free_y) "free_y" else "fixed",
      axes = if (repeat_y_axes) "all_y" else "margins",
      axis.labels = if (repeat_y_axes) "all_y" else "all"
    ) +
    scale_colour_manual(values = palette, drop = FALSE) +
    scale_shape_manual(
      values = c("Young" = 16, "Reproductively aged" = 17),
      drop = FALSE
    ) +
    scale_x_continuous(
      breaks = c(20, 30, 40, 50, 60),
      limits = c(20.5, 57),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Human chronological age (years)",
      y = y_label,
      colour = NULL,
      shape = NULL,
      caption = paste(
        "Donors, n = 8. Lines and shading show OLS fits and 95% confidence",
        "intervals. Spearman rank tests are exact."
      )
    ) +
    theme_classic(base_size = 7.2, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.38, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = 7.2, colour = "black"),
      axis.text = element_text(size = 6.7, colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = 7.4, face = "bold", colour = "black"),
      plot.title = element_text(size = 8.4, face = "bold", colour = "black"),
      plot.subtitle = element_text(
        size = 6.8,
        colour = "#333333",
        margin = margin(b = 5)
      ),
      plot.caption = element_text(size = 5.7, colour = "#4A4A4A", hjust = 0),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 6.5),
      legend.key.width = grid::unit(4.5, "mm"),
      panel.spacing = grid::unit(9, "mm"),
      plot.margin = margin(6, 7, 5, 6)
    )
}

p_mouse_week <- make_plot(
  "mouse_week_equivalent_output_audit_only",
  stats_mouse_week,
  "Frozen mouse ovarian clock output versus human chronological age",
  paste(
    "GSE202601 whole-ovary snRNA-seq pseudobulk; mouse-week equivalents are",
    "cross-species audit values, not predicted human age"
  ),
  "Frozen clock output (mouse-week equivalent; audit only)",
  free_y = TRUE,
  repeat_y_axes = FALSE
)

p_standardized <- make_plot(
  "cross_species_mouse_ovarian_clock_score",
  stats_standardized,
  "Standardized frozen mouse ovarian clock score versus human chronological age",
  paste(
    "GSE202601 whole-ovary snRNA-seq pseudobulk; standardized scores are",
    "cross-species audit values, not predicted human age"
  ),
  "Within-cohort standardized frozen mouse-clock score (z)",
  free_y = FALSE,
  repeat_y_axes = TRUE
)

source_data <- x[, .(
  projection_sample_id,
  sample_id,
  gsm,
  age_year,
  age_group,
  model_family,
  model_label,
  mouse_week_equivalent_output_audit_only,
  cross_species_mouse_ovarian_clock_score,
  raw_linear_predictor,
  mapped_gene_count,
  raw_gene_coverage,
  absolute_model_weight_coverage,
  projection_qc_class
)]
fwrite(
  source_data,
  file.path(source_dir, "predicted_vs_chronological_age_source_data.tsv"),
  sep = "\t"
)
fwrite(
  stats_mouse_week[
    , !c("annotation_x", "annotation_y", "annotation"),
    with = FALSE
  ],
  file.path(source_dir, "predicted_vs_chronological_age_statistics.tsv"),
  sep = "\t"
)
fwrite(
  source_data,
  file.path(
    source_dir,
    "standardized_score_vs_chronological_age_source_data.tsv"
  ),
  sep = "\t"
)
fwrite(
  stats_standardized[
    , !c("annotation_x", "annotation_y", "annotation"),
    with = FALSE
  ],
  file.path(
    source_dir,
    "standardized_score_vs_chronological_age_statistics.tsv"
  ),
  sep = "\t"
)

width_mm <- 178
height_mm <- 96
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

export_plot <- function(plot, stem) {
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"),
    width = width_in,
    height = height_in,
    family = "sans"
  )
  print(plot)
  dev.off()

  grDevices::svg(
    paste0(stem, ".svg"),
    width = width_in,
    height = height_in,
    family = "sans",
    onefile = TRUE
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
    compression = "lzw",
    background = "white"
  )
  print(plot)
  dev.off()
}

stem_mouse_week <- file.path(plot_dir, "01_predicted_vs_chronological_age")
stem_standardized <- file.path(
  plot_dir,
  "02_standardized_clock_score_vs_chronological_age"
)
export_plot(p_mouse_week, stem_mouse_week)
export_plot(p_standardized, stem_standardized)

make_qa <- function(stem, plot, stats, y_pattern, source_path, stats_path) {
  data.table(
    check = c(
      "sixteen_model_sample_rows",
      "eight_unique_donors",
      "two_model_panels",
      "all_values_finite",
      "spearman_exact_primary",
      "pearson_sensitivity_retained",
      "influence_diagnostics_retained",
      "chronological_age_x_axis_retained",
      "requested_y_axis_explicit",
      "figure2_y_axes_repeated_when_requested",
      "source_data_written",
      "statistics_written",
      "pdf_written",
      "svg_written",
      "png_written",
      "tiff_written"
    ),
    pass = c(
      nrow(x) == 16L,
      uniqueN(x$projection_sample_id) == 8L,
      uniqueN(x$model_family) == 2L,
      all(is.finite(x$age_year)),
      all(grepl("Spearman", stats$annotation)) &&
        all(grepl("Exact P", stats$annotation)),
      all(is.finite(stats$pearson_r_sensitivity)) &&
        all(is.finite(stats$pearson_95ci_lower)) &&
        all(is.finite(stats$pearson_95ci_upper)),
      all(is.finite(stats$max_cooks_distance)) &&
        all(is.finite(stats$max_abs_externally_studentized_residual)),
      identical(plot$labels$x, "Human chronological age (years)"),
      grepl(y_pattern, plot$labels$y),
      if (grepl("standardized", y_pattern)) {
        identical(plot$facet$params$draw_axes$y, TRUE) &&
          identical(plot$facet$params$axis_labels$y, TRUE)
      } else {
        TRUE
      },
      file.exists(source_path),
      file.exists(stats_path),
      file.exists(paste0(stem, ".pdf")),
      file.exists(paste0(stem, ".svg")),
      file.exists(paste0(stem, ".png")),
      file.exists(paste0(stem, ".tiff"))
    )
  )
}

qa_mouse_week <- make_qa(
  stem_mouse_week,
  p_mouse_week,
  stats_mouse_week,
  "mouse-week equivalent",
  file.path(source_dir, "predicted_vs_chronological_age_source_data.tsv"),
  file.path(source_dir, "predicted_vs_chronological_age_statistics.tsv")
)
qa_standardized <- make_qa(
  stem_standardized,
  p_standardized,
  stats_standardized,
  "standardized frozen mouse-clock score",
  file.path(source_dir, "standardized_score_vs_chronological_age_source_data.tsv"),
  file.path(source_dir, "standardized_score_vs_chronological_age_statistics.tsv")
)
fwrite(
  qa_mouse_week,
  file.path(plot_dir, "01_predicted_vs_chronological_age_QA.tsv"),
  sep = "\t"
)
fwrite(
  qa_standardized,
  file.path(plot_dir, "02_standardized_clock_score_vs_chronological_age_QA.tsv"),
  sep = "\t"
)
if (!all(qa_mouse_week$pass) || !all(qa_standardized$pass)) {
  stop("Parallel figure QA failed.")
}

message("PASS: parallel mouse-week and standardized-score figures exported.")
