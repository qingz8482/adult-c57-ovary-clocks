#!/usr/bin/env Rscript

# Cell-class composition of every stress-test library (v2 labels):
#   PRJNA863443 VCD (vehicle vs VCD by injection age / post-injection time) and Foxl2 (WT vs +/- by age),
#   GSE235021 rice diet (black, white/brown vs polished at 6 and 9 months).
# Cells counted with the clock cell policy (curation confirmed/probable, label confidence >= 0.70),
# restricted to the 8 clock classes; percentages are within those 8 classes (sum to 100 per library).
# Delta = case-library percentage - mean percentage of its matched control libraries (percentage points).
# Outputs in step_08 .../05_stress_test_projection/cell composition.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

review_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review")
run_root <- file.path(review_root, "reproduce")
stress_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection")
out_dir <- file.path(stress_dir, "cell composition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

classes <- c("Theca_Stroma", "Granulosa", "Endothelial", "Epithelial", "T_NK",
             "Pericyte_Smooth_Muscle", "Immune_Myeloid", "B_Plasma")
class_palette <- c(B_Plasma = "#1B9E77", Endothelial = "#00A6D6", Epithelial = "#7CAE00", Granulosa = "#F05AA6",
                   Immune_Myeloid = "#D95F02", Pericyte_Smooth_Muscle = "#6A3D9A", T_NK = "#0072B2",
                   Theca_Stroma = "#C77CFF")
min_confidence <- 0.70

# ---------------------------------------------------------------- PRJNA863443 (v2 consensus labels)
kim <- fread(cmd = paste("zcat", shQuote(file.path(review_root, "v2/04_consensus/harmonized_annotation_consensus_v2.tsv.gz"))),
             select = c("dataset_id", "sample_id", "annotation_label_broad", "annotation_confidence", "curation_status"))
kim <- kim[dataset_id == "2025_Kim_PRJNA863443_mouse_menopause_models" & grepl("^(CTL|VCD|Foxl2)_", sample_id)]
kim[, passes := curation_status %in% c("confirmed", "probable") & annotation_confidence >= min_confidence]
kim_cells <- kim[, .(cohort = "PRJNA863443", sample_id, sample_label = sample_id,
                     harmonized_cell_class = annotation_label_broad, passes)]
kim_meta <- unique(kim[, .(sample_id)])
kim_meta[, `:=`(
  family = fifelse(grepl("^Foxl2", sample_id), "Foxl2", "VCD"),
  role = fifelse(grepl("^(CTL|Foxl2_wt)", sample_id), "control", "case"),
  contrast = fcase(
    grepl("^(CTL|VCD)_", sample_id), paste0("VCD ", sub("^(CTL|VCD)_(\\w+)_(\\w+)_\\d+$", "\\2 \\3", sample_id)),
    grepl("_young_", sample_id), "Foxl2 4 months",
    grepl("_old_", sample_id), "Foxl2 9 months"),
  group = fcase(grepl("^CTL", sample_id), "Vehicle control", grepl("^VCD", sample_id), "VCD",
                grepl("^Foxl2_wt", sample_id), "Foxl2 WT", grepl("^Foxl2_het", sample_id), "Foxl2+/-")
)]

# ---------------------------------------------------------------- GSE235021 (reproduce v2 labels)
gse <- fread(cmd = paste("zcat", shQuote(file.path(stress_dir, "GSE235021/tables/GSE235021_harmonized_cell_metadata.tsv.gz"))),
             select = c("sample_id", "sample_title", "age_month", "diet", "harmonized_cell_class", "passes_confidence_policy"))
gse_cells <- gse[, .(cohort = "GSE235021", sample_id, sample_label = sample_title, harmonized_cell_class,
                     passes = passes_confidence_policy)]
diet_lab <- c(polished_rice = "Polished rice", black_rice = "Black rice", white_brown_rice = "White/brown rice")
gse_meta <- unique(gse[, .(sample_id, age_month, diet)])
# Each diet library compared with the polished-rice libraries of the same age; polished libraries are the controls
gse_meta <- rbind(
  gse_meta[diet != "polished_rice", .(sample_id, family = "Rice diet", role = "case",
                                      contrast = paste0("Rice ", age_month, " mo: ", sub(" rice", "", diet_lab[diet]), " vs polished"),
                                      group = unname(diet_lab[diet]))],
  rbindlist(lapply(c("Black rice", "White/brown rice"), function(cd) {
    gse_meta[diet == "polished_rice", .(sample_id, family = "Rice diet", role = "control",
                                        contrast = paste0("Rice ", age_month, " mo: ", sub(" rice", "", cd), " vs polished"),
                                        group = "Polished rice")]
  }))
)

# ---------------------------------------------------------------- composition per library
cells <- rbind(kim_cells, gse_cells)
cells <- cells[(passes) & harmonized_cell_class %in% classes]
counts <- dcast(cells[, .N, by = .(cohort, sample_id, sample_label, harmonized_cell_class)],
                cohort + sample_id + sample_label ~ harmonized_cell_class, value.var = "N", fill = 0)
comp <- melt(counts, id.vars = c("cohort", "sample_id", "sample_label"), variable.name = "harmonized_cell_class",
             value.name = "n_cells")
comp[, harmonized_cell_class := as.character(harmonized_cell_class)]
comp[, n_cells_8_classes := sum(n_cells), by = sample_id]
comp[, percent := 100 * n_cells / n_cells_8_classes]
totals <- rbind(kim_cells, gse_cells)[, .(n_cells_all_labelled = .N, n_cells_pass_policy = sum(passes)), by = sample_id]
comp <- merge(comp, totals, by = "sample_id")
stopifnot(all(abs(comp[, sum(percent), by = sample_id]$V1 - 100) < 1e-8))

meta <- rbind(kim_meta, gse_meta)
contrast_levels <- c(paste("VCD", c("3m 30d", "3m 90d", "10m 30d", "10m 90d")), "Foxl2 4 months", "Foxl2 9 months",
                     paste0("Rice ", rep(c(6, 9), each = 2), " mo: ", c("Black", "White/brown"), " vs polished"))
stopifnot(all(meta$contrast %in% contrast_levels))
comp_long <- merge(comp, meta, by = "sample_id", allow.cartesian = TRUE)

fwrite(comp[order(cohort, sample_id, match(harmonized_cell_class, classes))],
       file.path(out_dir, "stress_test_library_cell_class_composition.tsv"), sep = "\t")

# ---------------------------------------------------------------- delta vs matched controls
ctrl_mean <- comp_long[role == "control", .(control_percent = mean(percent), n_control_libraries = uniqueN(sample_id)),
                       by = .(contrast, harmonized_cell_class)]
case <- merge(comp_long[role == "case"], ctrl_mean, by = c("contrast", "harmonized_cell_class"))
case[, delta_percent_points := percent - control_percent]
delta_summary <- case[, .(case_percent = mean(percent), control_percent = control_percent[1],
                          mean_delta_percent_points = mean(delta_percent_points),
                          sd_delta_percent_points = sd(delta_percent_points),
                          n_case_libraries = uniqueN(sample_id), n_control_libraries = n_control_libraries[1]),
                      by = .(family, contrast, harmonized_cell_class)]
fwrite(case[order(match(contrast, contrast_levels), match(harmonized_cell_class, classes), sample_id),
            .(family, contrast, group, sample_id, sample_label, harmonized_cell_class, n_cells, n_cells_8_classes,
              percent, control_percent, n_control_libraries, delta_percent_points)],
       file.path(out_dir, "stress_test_cell_class_delta_by_case_library.tsv"), sep = "\t")
fwrite(delta_summary[order(match(contrast, contrast_levels), match(harmonized_cell_class, classes))],
       file.path(out_dir, "stress_test_cell_class_delta_summary.tsv"), sep = "\t")

# ---------------------------------------------------------------- figures
class_lab <- c(Theca_Stroma = "Theca / Stroma", Granulosa = "Granulosa", Endothelial = "Endothelial",
               Epithelial = "Epithelial", T_NK = "T / NK", Pericyte_Smooth_Muscle = "Smooth Muscle",
               Immune_Myeloid = "Immune / Myeloid", B_Plasma = "B Plasma")
theme_comp <- function(base_size = 8) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.25),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
          legend.key = element_blank())
}

# (1) Delta barplot: one panel per matched contrast
for (d in list(delta_summary, case)) {
  d[, `:=`(contrast = factor(contrast, levels = contrast_levels),
           harmonized_cell_class = factor(harmonized_cell_class, levels = classes))]
}
n_lab <- delta_summary[, .(lab = sprintf("n = %d vs %d", n_case_libraries[1], n_control_libraries[1])), by = contrast]
p_delta <- ggplot(delta_summary, aes(harmonized_cell_class, mean_delta_percent_points, fill = harmonized_cell_class)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.2, alpha = 0.9) +
  geom_point(data = case, aes(harmonized_cell_class, delta_percent_points), inherit.aes = FALSE,
             shape = 21, fill = "white", colour = "black", size = 1.1, stroke = 0.3,
             position = position_jitter(width = 0.12, height = 0, seed = 1)) +
  geom_text(data = n_lab, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.4, size = 2.2, colour = "grey30") +
  facet_wrap(~ contrast, ncol = 4, drop = FALSE) +
  scale_fill_manual(values = class_palette, labels = class_lab, guide = "none") +
  scale_x_discrete(labels = class_lab) +
  labs(title = "Cell-class composition change vs matched controls",
       subtitle = paste("Bars: mean of case libraries; points: each case library. Percent of cells in the 8 clock classes",
                        "(clock cell policy: confirmed/probable, confidence >= 0.70); delta = case % - mean matched-control %."),
       x = NULL, y = "Delta cell-class percentage (percentage points)") +
  theme_comp() +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1), plot.subtitle = element_text(size = 7))
ggsave(file.path(out_dir, "stress_test_cell_class_delta_barplot.pdf"), p_delta, width = 11, height = 8.5, device = cairo_pdf)
ggsave(file.path(out_dir, "stress_test_cell_class_delta_barplot.png"), p_delta, width = 11, height = 8.5, dpi = 300, bg = "white")

# (2) Composition of every library (stacked), controls first within each contrast
lib <- unique(comp_long[, .(sample_id, sample_label, family, role, group)])
lib <- lib[, .SD[1], by = sample_id]   # polished-rice controls appear in two contrasts; draw each library once
lib[, panel := fcase(family == "VCD", "PRJNA863443 VCD", family == "Foxl2", "PRJNA863443 Foxl2",
                     default = "GSE235021 rice diet")]
lib_order <- c(
  sort(kim_meta[family == "VCD"][order(sub("^(CTL|VCD)_", "", sample_id), role != "control"), sample_id]),
  kim_meta[family == "Foxl2"][order(!grepl("young", sample_id), role != "control", sample_id), sample_id],
  gse[, .(sample_id, sample_title, age_month, diet)][, unique(.SD)][order(age_month, match(diet, names(diet_lab)), sample_title), sample_id]
)
lib_order <- unique(lib_order)
vcd_order <- kim_meta[family == "VCD"][, key := sub("_\\d+$", "", sub("^(CTL|VCD)_", "", sample_id))][
  order(match(key, c("3m_30d", "3m_90d", "10m_30d", "10m_90d")), role != "control", sample_id), sample_id]
lib_order <- c(vcd_order, setdiff(lib_order, vcd_order))
stack <- merge(comp, lib[, .(sample_id, panel)], by = "sample_id")
stack[, `:=`(sample_label = factor(sample_label, levels = unique(comp[match(lib_order, sample_id), sample_label])),
             harmonized_cell_class = factor(harmonized_cell_class, levels = rev(classes)),
             panel = factor(panel, levels = c("PRJNA863443 VCD", "PRJNA863443 Foxl2", "GSE235021 rice diet")))]
p_stack <- ggplot(stack, aes(sample_label, percent, fill = harmonized_cell_class)) +
  geom_col(width = 0.85, colour = "white", linewidth = 0.1) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = class_palette, labels = class_lab, breaks = classes, name = "Cell class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  labs(title = "Cell-class composition of each stress-test library",
       subtitle = "Percent of cells in the 8 clock classes (clock cell policy: confirmed/probable, confidence >= 0.70).",
       x = NULL, y = "Cells (%)") +
  theme_comp() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6), legend.position = "right")
ggsave(file.path(out_dir, "stress_test_library_cell_class_composition.pdf"), p_stack, width = 12, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "stress_test_library_cell_class_composition.png"), p_stack, width = 12, height = 5, dpi = 300, bg = "white")

message("Libraries: ", uniqueN(comp$sample_id), " (PRJNA863443 ", uniqueN(kim_cells$sample_id), ", GSE235021 ",
        uniqueN(gse_cells$sample_id), "); cells counted: ", nrow(cells))
message("Wrote ", out_dir)
