#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_arg <- script_arg[grepl("^--file=", script_arg)]
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[[1]]) else NA_character_
script_dir <- if (!is.na(script_file)) dirname(normalizePath(script_file)) else file.path(getwd(), "scripts")
if (!dir.exists(script_dir)) {
  script_dir <- file.path(getwd(), "adult_c57_ovary_clocks", "bulk-pseudobulk_clock", "scripts")
}

project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
results_dir <- normalizePath(Sys.getenv("RESULTS_ROOT", unset = file.path(project_dir, "results")), mustWork = TRUE)
projection_dir <- file.path(results_dir, "external_projection")
simplified_dir <- file.path(projection_dir, "simplified_classes")
diagnostics_dir <- file.path(projection_dir, "diagnostics")
output_subdir <- Sys.getenv("RAW_DELTA_OUTPUT_SUBDIR", unset = "raw delta_age by dataset_with t test")
label_style <- Sys.getenv("RAW_DELTA_LABEL_STYLE", unset = "exact_pvalue")
plot_test <- Sys.getenv("RAW_DELTA_PLOT_TEST", unset = "welch_t")
if (!plot_test %in% c("welch_t", "permutation")) {
  stop("RAW_DELTA_PLOT_TEST must be either 'welch_t' or 'permutation'.", call. = FALSE)
}
out_dir <- file.path(diagnostics_dir, output_subdir)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nchar(x) == 0, "missing", x)
}

theme_projection <- function() {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10.5),
      axis.text.x = element_text(angle = 35, hjust = 1, margin = margin(t = 4)),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(10, 20, 18, 20)
    )
}

save_plot <- function(p, prefix, width = 10, height = 6.8) {
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  ggsave(pdf_path, p, width = width, height = height, limitsize = FALSE)
  ggsave(png_path, p, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  if (identical(prefix, "04_raw_delta_age_by_dataset_treatment_GSE329024")) {
    ggsave(file.path(diagnostics_dir, paste0(prefix, ".pdf")), p, width = 10, height = 6.5, limitsize = FALSE)
    ggsave(file.path(diagnostics_dir, paste0(prefix, ".png")), p, width = 10, height = 6.5, dpi = 300, bg = "white", limitsize = FALSE)
  }
  invisible(c(pdf_path, png_path))
}

truthy <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

read_predictions <- function() {
  candidates <- c(
    file.path(simplified_dir, "external_projection_predictions_simplified.tsv"),
    file.path(projection_dir, "external_projection_predictions.tsv")
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) stop("No external projection predictions table found.", call. = FALSE)
  dt <- fread(path, sep = "\t", na.strings = c("", "NA", "NaN"))
  dt[, source_predictions_path := path]
  dt
}

format_p <- function(p) {
  fifelse(
    is.na(p), NA_character_,
    fifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.4f", p))
  )
}

p_to_star <- function(p) {
  fifelse(
    is.na(p), "NA",
    fifelse(p < 1e-4, "****",
      fifelse(p < 1e-3, "***",
        fifelse(p < 1e-2, "**",
          fifelse(p < 0.05, "*", "ns")
        )
      )
    )
  )
}

annotation_label <- function(p, star, style = "stars") {
  if (identical(style, "exact_pvalue")) {
    return(paste0("p=", format_p(p)))
  }
  if (identical(style, "ns_pvalue")) {
    return(fifelse(star == "ns", paste0("p=", format_p(p)), star))
  }
  star
}

mean_diff <- function(condition, control) {
  mean(condition, na.rm = TRUE) - mean(control, na.rm = TRUE)
}

permutation_p <- function(condition, control, max_exact = 20000L, n_perm = 10000L) {
  condition <- condition[is.finite(condition)]
  control <- control[is.finite(control)]
  if (length(condition) == 0 || length(control) == 0) return(list(p = NA_real_, method = "not_tested"))

  observed <- mean_diff(condition, control)
  values <- c(condition, control)
  n_condition <- length(condition)
  n_total <- length(values)
  n_combinations <- choose(n_total, n_condition)

  if (is.finite(n_combinations) && n_combinations <= max_exact) {
    idx <- combn(n_total, n_condition)
    diffs <- apply(idx, 2, function(ii) {
      mean(values[ii]) - mean(values[-ii])
    })
    p <- mean(abs(diffs) >= abs(observed) - .Machine$double.eps)
    return(list(p = p, method = paste0("exact_permutation_", n_combinations, "_partitions")))
  }

  diffs <- replicate(n_perm, {
    ii <- sample.int(n_total, n_condition, replace = FALSE)
    mean(values[ii]) - mean(values[-ii])
  })
  p <- (sum(abs(diffs) >= abs(observed) - .Machine$double.eps) + 1) / (n_perm + 1)
  list(p = p, method = paste0("random_permutation_", n_perm, "_draws"))
}

bootstrap_ci <- function(condition, control, n_boot = 10000L) {
  condition <- condition[is.finite(condition)]
  control <- control[is.finite(control)]
  if (length(condition) == 0 || length(control) == 0) return(c(NA_real_, NA_real_))
  diffs <- replicate(n_boot, {
    mean(sample(condition, length(condition), replace = TRUE)) -
      mean(sample(control, length(control), replace = TRUE))
  })
  as.numeric(stats::quantile(diffs, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

welch_p <- function(condition, control) {
  condition <- condition[is.finite(condition)]
  control <- control[is.finite(control)]
  if (length(condition) < 2 || length(control) < 2) return(NA_real_)
  out <- tryCatch(stats::t.test(condition, control, var.equal = FALSE)$p.value, error = function(e) NA_real_)
  as.numeric(out)
}

contrast_map <- data.table::rbindlist(list(
  data.table(dataset_folder = "CRA017492", contrast_label = "NMN vs aged water", condition_group = "aged NMN", control_group = "aged water control"),
  data.table(dataset_folder = "GSE194196", contrast_label = "Bnc1tr/tr vs WT", condition_group = "Bnc1tr/tr", control_group = "reference control"),
  data.table(dataset_folder = "GSE240547", contrast_label = "HFHS vs control diet", condition_group = "HFHS diet", control_group = "control diet"),
  data.table(dataset_folder = "GSE252791", contrast_label = "Rnf111 het KO vs WT", condition_group = "Rnf111 heterozygous KO", control_group = "WT control"),
  data.table(dataset_folder = "GSE273724", contrast_label = "adult HEP14 vs vehicle", condition_group = "adult 27-28 week plus HEP14", control_group = "adult 27-28 week vehicle"),
  data.table(dataset_folder = "GSE273724", contrast_label = "aged HEP14 vs vehicle", condition_group = "aged 73-74 week plus HEP14", control_group = "aged 73-74 week vehicle"),
  data.table(dataset_folder = "GSE279446", contrast_label = "young Zeb1_200H vs young WT", condition_group = "young Zeb1_200H", control_group = "young WT"),
  data.table(dataset_folder = "GSE279446", contrast_label = "young Zeb1_200M vs young WT", condition_group = "young Zeb1_200M", control_group = "young WT"),
  data.table(dataset_folder = "GSE279446", contrast_label = "adult Zeb1_200H vs adult WT", condition_group = "adult Zeb1_200H", control_group = "adult WT"),
  data.table(dataset_folder = "GSE279446", contrast_label = "adult Zeb1_200M vs adult WT", condition_group = "adult Zeb1_200M", control_group = "adult WT"),
  data.table(dataset_folder = "GSE279573", contrast_label = "KO chow infertile vs WT chow fertile", condition_group = "KO chow infertile", control_group = "WT chow fertile"),
  data.table(dataset_folder = "GSE279573", contrast_label = "KO rescue vs WT chow fertile", condition_group = "KO vegetable-oil-free rescue", control_group = "WT chow fertile"),
  data.table(dataset_folder = "GSE279573", contrast_label = "KO chow infertile vs KO rescue", condition_group = "KO chow infertile", control_group = "KO vegetable-oil-free rescue"),
  data.table(dataset_folder = "GSE329024", contrast_label = "2M-8W CUMS vs matched control", condition_group = "2M-8W CUMS stress", control_group = "2M-8W control"),
  data.table(dataset_folder = "GSE329024", contrast_label = "2M-12W CUMS vs matched control", condition_group = "2M-12W CUMS stress", control_group = "2M-12W control"),
  data.table(dataset_folder = "GSE329024", contrast_label = "6M-8W CUMS vs matched control", condition_group = "6M-8W CUMS stress", control_group = "6M-8W control"),
  data.table(dataset_folder = "GSE329024", contrast_label = "6M-12W CUMS vs matched control", condition_group = "6M-12W CUMS stress", control_group = "6M-12W control"),
  data.table(dataset_folder = "PRJNA1076509", contrast_label = "estropausal donor FMT vs young donor FMT", condition_group = "estropausal female donor FMT", control_group = "young female donor FMT"),
  data.table(dataset_folder = "PRJNA863443", contrast_label = "VCD 3m 30d vs CTL", condition_group = "VCD_3m_30d", control_group = "CTL_3m_30d"),
  data.table(dataset_folder = "PRJNA863443", contrast_label = "VCD 3m 90d vs CTL", condition_group = "VCD_3m_90d", control_group = "CTL_3m_90d"),
  data.table(dataset_folder = "PRJNA863443", contrast_label = "VCD 10m 30d vs CTL", condition_group = "VCD_10m_30d", control_group = "CTL_10m_30d"),
  data.table(dataset_folder = "PRJNA863443", contrast_label = "VCD 10m 90d vs CTL", condition_group = "VCD_10m_90d", control_group = "CTL_10m_90d"),
  data.table(dataset_folder = "PRJNA863443", contrast_label = "Foxl2+/- young vs Foxl2+/+ young", condition_group = "Foxl2+/-_young", control_group = "Foxl2+/+_young"),
  data.table(dataset_folder = "PRJNA863443", contrast_label = "Foxl2+/- mid-age vs Foxl2+/+ mid-age", condition_group = "Foxl2+/-_mid_age", control_group = "Foxl2+/+_mid_age")
))

dataset_plot_group_levels <- list(
  GSE279446 = c(
    "young WT", "young Zeb1_200H", "young Zeb1_200M",
    "adult WT", "adult Zeb1_200H", "adult Zeb1_200M"
  ),
  GSE279573 = c(
    "WT chow fertile",
    "KO chow infertile",
    "KO vegetable-oil-free rescue"
  ),
  GSE329024 = c(
    "2M-8W control",
    "2M-8W CUMS stress",
    "2M-12W control",
    "2M-12W CUMS stress",
    "6M-8W control",
    "6M-8W CUMS stress",
    "6M-12W control",
    "6M-12W CUMS stress"
  )
)

set.seed(20260710)

pred <- read_predictions()
required_cols <- c("sample_key", "dataset_folder", "model_family", "delta_age_week")
missing_cols <- setdiff(required_cols, names(pred))
if (length(missing_cols) > 0) {
  stop("Missing required columns in predictions table: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

if ("plot_include_primary" %in% names(pred)) {
  pred <- pred[truthy(plot_include_primary)]
}

pred[, delta_age_week := as.numeric(delta_age_week)]
if ("projection_group" %in% names(pred)) {
  pred[, plot_group := as.character(projection_group)]
} else if ("treatment" %in% names(pred)) {
  pred[, plot_group := as.character(treatment)]
} else {
  stop("Predictions table needs projection_group or treatment to define plot groups.", call. = FALSE)
}
pred[is.na(plot_group) | plot_group == "", plot_group := as.character(treatment)]
pred[is.na(plot_group) | plot_group == "", plot_group := as.character(projection_class)]
pred[, projection_class_label := as.character(projection_class)]
pred[projection_class_label == "stress", projection_class_label := "stress_model"]

plot_dt <- pred[dataset_folder %in% unique(contrast_map$dataset_folder) & is.finite(delta_age_week)]
plot_dt <- unique(plot_dt, by = c("dataset_folder", "model_family", "sample_key", "plot_group", "delta_age_week"))

stats_rows <- vector("list", nrow(contrast_map) * max(1L, uniqueN(plot_dt$model_family)))
ii <- 0L
for (rr in seq_len(nrow(contrast_map))) {
  cm <- contrast_map[rr]
  for (mf in sort(unique(plot_dt$model_family))) {
    d <- plot_dt[
      dataset_folder == cm$dataset_folder &
        model_family == mf &
        plot_group %in% c(cm$condition_group, cm$control_group)
    ]
    condition <- d[plot_group == cm$condition_group, delta_age_week]
    control <- d[plot_group == cm$control_group, delta_age_week]
    perm <- permutation_p(condition, control)
    ci <- bootstrap_ci(condition, control)
    ii <- ii + 1L
    stats_rows[[ii]] <- data.table(
      dataset_folder = cm$dataset_folder,
      model_family = mf,
      contrast_label = cm$contrast_label,
      condition_group = cm$condition_group,
      control_group = cm$control_group,
      n_condition = length(condition),
      n_control = length(control),
      mean_condition_delta_age_week = mean(condition, na.rm = TRUE),
      mean_control_delta_age_week = mean(control, na.rm = TRUE),
      mean_difference_week = mean_diff(condition, control),
      bootstrap_ci_low_week = ci[[1]],
      bootstrap_ci_high_week = ci[[2]],
      p_permutation = perm$p,
      permutation_method = perm$method,
      p_welch_t = welch_p(condition, control)
    )
  }
}
test_results <- rbindlist(stats_rows[seq_len(ii)], fill = TRUE)
test_results[, p_permutation_fdr_model := p.adjust(p_permutation, method = "BH"), by = model_family]
test_results[, p_welch_t_fdr_model := p.adjust(p_welch_t, method = "BH"), by = model_family]
test_results[, p_star := p_to_star(p_permutation)]
test_results[, p_star_fdr_model := p_to_star(p_permutation_fdr_model)]
test_results[, `:=`(
  plot_test = plot_test,
  p_plot_value = if (plot_test == "welch_t") p_welch_t else p_permutation
)]
test_results[, p_plot_star := p_to_star(p_plot_value)]
test_results[, p_plot_label := annotation_label(p_plot_value, p_plot_star, label_style)]
test_results[, `:=`(
  p_permutation_label = format_p(p_permutation),
  p_permutation_fdr_model_label = format_p(p_permutation_fdr_model),
  p_welch_t_label = format_p(p_welch_t),
  p_welch_t_fdr_model_label = format_p(p_welch_t_fdr_model)
)]
setorder(test_results, dataset_folder, model_family, contrast_label)

group_summary <- plot_dt[, .(
  n_samples = uniqueN(sample_key),
  mean_delta_age_week = mean(delta_age_week, na.rm = TRUE),
  median_delta_age_week = median(delta_age_week, na.rm = TRUE),
  sd_delta_age_week = sd(delta_age_week, na.rm = TRUE)
), by = .(dataset_folder, model_family, plot_group)]
setorder(group_summary, dataset_folder, model_family, plot_group)

fwrite(contrast_map, file.path(out_dir, "raw_delta_age_condition_control_contrast_map.tsv"), sep = "\t", na = "NA")
fwrite(group_summary, file.path(out_dir, "raw_delta_age_group_summary_by_dataset_model.tsv"), sep = "\t", na = "NA")
fwrite(test_results, file.path(out_dir, "raw_delta_age_condition_vs_control_tests_by_dataset_model.tsv"), sep = "\t", na = "NA")

for (ds in sort(unique(plot_dt$dataset_folder))) {
  dt <- copy(plot_dt[dataset_folder == ds])
  ds_map <- contrast_map[dataset_folder == ds]
  if (ds %in% names(dataset_plot_group_levels)) {
    ordered_levels <- unique(c(dataset_plot_group_levels[[ds]], sort(unique(dt$plot_group))))
  } else if (ds == "GSE273724") {
    paired_contrast_levels <- as.vector(t(as.matrix(ds_map[, .(control_group, condition_group)])))
    ordered_levels <- unique(c(paired_contrast_levels, sort(unique(dt$plot_group))))
  } else {
    ordered_levels <- unique(c(ds_map$control_group, ds_map$condition_group, sort(unique(dt$plot_group))))
  }
  ordered_levels <- ordered_levels[ordered_levels %in% dt$plot_group]
  dt[, plot_group := factor(plot_group, levels = ordered_levels)]

  ann <- copy(test_results[dataset_folder == ds & n_condition > 0 & n_control > 0])
  ann <- ann[condition_group %in% ordered_levels & control_group %in% ordered_levels]
  if (nrow(ann) > 0) {
    ranges <- dt[, .(
      y_min = min(delta_age_week, na.rm = TRUE),
      y_max = max(delta_age_week, na.rm = TRUE)
    ), by = model_family]
    ranges[, y_span := pmax(y_max - y_min, 1)]
    ann <- merge(ann, ranges, by = "model_family", all.x = TRUE)
    ann[, contrast_index := seq_len(.N), by = model_family]
    ann[, `:=`(
      x_control = match(control_group, ordered_levels),
      x_condition = match(condition_group, ordered_levels),
      y = y_max + y_span * (0.12 + 0.13 * (contrast_index - 1)),
      y_tip_low = y_max + y_span * (0.08 + 0.13 * (contrast_index - 1)),
      y_label = y_max + y_span * (0.16 + 0.13 * (contrast_index - 1)),
      p_label = p_plot_label
    )]
    ann[, `:=`(
      xmin = pmin(x_control, x_condition),
      xmax = pmax(x_control, x_condition)
    )]
    ann[, y_plot_max := max(y_label, na.rm = TRUE), by = model_family]
  }

  p <- ggplot(dt, aes(plot_group, delta_age_week, color = projection_class_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45") +
    geom_boxplot(outlier.shape = NA, alpha = 0.12) +
    geom_jitter(width = 0.16, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~model_family, scales = "free_x") +
    labs(
      title = paste0("Raw delta age by projection group: ", ds),
      subtitle = if (plot_test == "welch_t") {
        "Unadjusted two-sided Welch t-tests; design-matched contrasts evaluated separately within each model family."
      } else if (identical(label_style, "ns_pvalue")) {
        "Labels show two-sided permutation p-values computed separately within each model family; significant contrasts use stars."
      } else {
        "Stars show two-sided permutation p-values computed separately within each model family."
      },
      x = "Projection group",
      y = "Delta age (projected - chronological weeks)",
      color = "Projection class"
    ) +
    theme_projection() +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    coord_cartesian(clip = "off")

  if (exists("ann") && nrow(ann) > 0) {
    p <- p +
      geom_segment(
        data = ann,
        aes(x = xmin, xend = xmax, y = y, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = ann,
        aes(x = xmin, xend = xmin, y = y_tip_low, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = ann,
        aes(x = xmax, xend = xmax, y = y_tip_low, yend = y),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_text(
        data = ann,
        aes(x = (xmin + xmax) / 2, y = y_label, label = p_label),
        inherit.aes = FALSE,
        size = if (identical(label_style, "ns_pvalue")) 3.1 else 3.6,
        fontface = "bold"
      ) +
      expand_limits(y = max(ann$y_label, na.rm = TRUE))
  }

  width <- max(10, min(15, 1.05 * length(ordered_levels) + 5))
  save_plot(p, paste0("04_raw_delta_age_by_dataset_treatment_", safe_name(ds)), width = width, height = 6.8)
}

message("Wrote plots and tables to: ", out_dir)
