#!/usr/bin/env Rscript
# Granulosa clock with vs without luteal-state cells (v2 main CV vs v2 no-luteal sensitivity CV).
# Reads existing out-of-fold predictions only; nothing is refitted.
#   main:       06_clock_run/cv            (Granulosa metacells include luteal-state cells)
#   no luteal:  06_clock_run/cv_no_luteal  (luteal-state Granulosa cells removed before metacell construction;
#                                           V2_EXCLUDE_GRANULOSA_LUTEAL=1, same CV protocol)
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

V2 <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
run_dir <- file.path(V2, "06_clock_run")
out_dir <- file.path(run_dir, "comparison", "luteal_sensitivity")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

oof_path <- function(run) file.path(run_dir, run, "step_08_figures_tables_release", "tables",
                                    "celltype_metacell_elasticnet_oof_predictions.tsv")
runs <- c(with_luteal = "cv", without_luteal = "cv_no_luteal")
lab <- c(with_luteal = "With luteal-state cells", without_luteal = "Luteal-state cells excluded")
col <- c(with_luteal = "#2a78d6", without_luteal = "#e1891a")

oof <- rbindlist(lapply(names(runs), function(k) {
  d <- fread(oof_path(runs[[k]]))[harmonized_cell_class == "Granulosa"]
  d[, run := k]
}))
oof[, run := factor(run, levels = names(runs))]
oof[, abs_err := abs(predicted_age_weeks - age_weeks)]

# dataset short labels
short <- function(x) {
  y <- sub("^(\\d{4})_[^_]+_([^_]+)_.*$", "\\2 \\1", x)
  ifelse(grepl("PRJNA863443", x), "Kim 2025", y)
}
oof[, dataset := short(dataset_id)]
bins <- c(young_adult = "Young adult", mature_adult = "Mature adult",
          advanced_maternal_age = "Advanced maternal age", post_reproductive_ageing = "Post-reproductive")
oof[, age_bin := factor(bins[reproductive_age_bin], levels = bins)]

# ---- metrics (same definitions as the release summary)
met <- oof[, .(n_metacells = .N, n_samples = uniqueN(sample_key),
               mae_weeks = mean(abs_err), rmse_weeks = sqrt(mean((predicted_age_weeks - age_weeks)^2)),
               pearson_r = cor(predicted_age_weeks, age_weeks),
               slope = coef(lm(predicted_age_weeks ~ age_weeks))[2],
               r2 = summary(lm(age_weeks ~ predicted_age_weeks))$r.squared), by = run]
fwrite(met, file.path(out_dir, "granulosa_luteal_sensitivity_metrics.tsv"), sep = "\t")

# cross-check against the stored CV summaries
for (k in names(runs)) {
  s <- fread(file.path(run_dir, runs[[k]], "step_08_figures_tables_release", "tables",
                       "celltype_metacell_elasticnet_cv_evaluation_summary.tsv"))[harmonized_cell_class == "Granulosa"]
  stopifnot(abs(s$mae_weeks - met[run == k, mae_weeks]) < 1e-8, s$n_metacells == met[run == k, n_metacells])
}

per_ds <- oof[, .(mae_weeks = mean(abs_err), n = .N), by = .(run, dataset)]
fwrite(per_ds, file.path(out_dir, "granulosa_luteal_sensitivity_mae_by_dataset.tsv"), sep = "\t")

# luteal fraction per Granulosa metacell (main build)
md <- fread(file.path(run_dir, "build", "step_02_preprocessing_filtering", "filtered_metacell_metadata.tsv"))
md <- md[harmonized_cell_class == "Granulosa"]
main <- merge(oof[run == "with_luteal"], md[, .(metacell_id, frac_luteal)], by = "metacell_id")
stopifnot(nrow(main) == met[run == "with_luteal", n_metacells])
main[, lut_bin := cut(frac_luteal, c(-Inf, 0, 0.25, 0.5, Inf), labels = c("0", "0–25%", "25–50%", ">50%"))]
fwrite(main[, .(n = .N, mae_weeks = mean(abs_err), median_abs_err = median(abs_err)), by = lut_bin][order(lut_bin)],
       file.path(out_dir, "granulosa_main_error_by_luteal_fraction.tsv"), sep = "\t")

FS <- 12  # all text: 12 pt Arial
theme_set(theme_classic(base_size = FS, base_family = "Arial") +
            theme(panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
                  strip.background = element_rect(fill = "grey95", colour = NA),
                  strip.text = element_text(face = "bold"), legend.position = "bottom",
                  plot.tag = element_text(face = "bold", size = FS),
                  text = element_text(size = FS, family = "Arial"),
                  axis.text = element_text(size = FS, colour = "black"),
                  strip.text.x = element_text(size = FS, face = "bold"),
                  legend.text = element_text(size = FS)))

# a. luteal fraction by dataset (context: why the sensitivity matters)
pa <- ggplot(main, aes(x = reorder(short(dataset_id), frac_luteal, median), y = 100 * frac_luteal)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, colour = "grey40", linewidth = 0.3) +
  geom_jitter(aes(colour = age_bin), width = 0.15, size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = c("#2a78d6", "#1baf7a", "#eda100", "#e87ba4"), name = NULL) +
  coord_flip() + labs(x = NULL, y = "Luteal-state cells per Granulosa metacell (%)") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2, alpha = 1)))

# b. out-of-fold predicted vs observed, both runs
lbl <- met[, .(run, txt = sprintf("n = %s\nMAE = %.2f wk\nr = %.2f\nslope = %.2f",
                                   format(n_metacells, big.mark = ","), mae_weeks, pearson_r, slope))]
oof[, run_lab := factor(lab[as.character(run)], levels = lab)]
lbl[, run_lab := factor(lab[as.character(run)], levels = lab)]
lim <- range(c(oof$age_weeks, oof$predicted_age_weeks))
pb <- ggplot(oof, aes(age_weeks, predicted_age_weeks, colour = run)) +
  geom_abline(linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_point(size = 0.6, alpha = 0.35, position = position_jitter(width = 0.6, height = 0, seed = 1)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
  geom_text(data = lbl, aes(x = lim[2], y = lim[1], label = txt), hjust = 1, vjust = 0, size = FS / .pt, family = "Arial",
            colour = "black", inherit.aes = FALSE) +
  facet_wrap(~run_lab) + scale_colour_manual(values = col, guide = "none") +
  coord_equal(xlim = lim, ylim = lim) +
  labs(x = "Chronological age (weeks)", y = "Out-of-fold predicted age (weeks)")

# c. per-dataset MAE, paired
per_ds[, run_lab := factor(lab[as.character(run)], levels = lab)]
ord <- per_ds[run == "with_luteal"][order(mae_weeks), dataset]
per_ds[, dataset := factor(dataset, levels = ord)]
pc <- ggplot(per_ds, aes(mae_weeks, dataset)) +
  geom_line(aes(group = dataset), colour = "grey70", linewidth = 0.5) +
  geom_point(aes(colour = run_lab), size = 2.2) +
  scale_colour_manual(values = setNames(col, lab), name = NULL) +
  labs(x = "Out-of-fold MAE by dataset (weeks)", y = NULL)

# d. main clock error vs luteal content of the metacell
nlab <- main[, .N, by = lut_bin]
pd <- ggplot(main, aes(lut_bin, abs_err)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, fill = "#2a78d6", alpha = 0.25, colour = "#2a78d6", linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.25, colour = "#2a78d6") +
  geom_text(data = nlab, aes(lut_bin, y = max(main$abs_err) * 1.05, label = paste0("n=", N)), size = FS / .pt, family = "Arial") +
  labs(x = "Luteal-state cells in metacell", y = "Absolute error, main clock (weeks)")

p <- (pa | pb) / (pc | pd) + plot_layout(widths = c(1, 1.6)) + plot_annotation(tag_levels = "a")
ggsave(file.path(out_dir, "granulosa_luteal_sensitivity.pdf"), p, width = 11, height = 10, device = cairo_pdf)
ggsave(file.path(out_dir, "granulosa_luteal_sensitivity.png"), p, width = 11, height = 10, dpi = 300)
print(met)
