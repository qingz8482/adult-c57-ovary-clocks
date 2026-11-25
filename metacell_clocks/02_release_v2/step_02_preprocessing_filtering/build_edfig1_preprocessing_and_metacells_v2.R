#!/usr/bin/env Rscript
# ED Fig. 1 (preprocessing and metacells) rebuilt on the v2 cell labels.
# Replaces the panels produced in the 20260707 run by update_prjna863443_scRNA_metadata_and_edfig1.R and
# regenerate_confident_celltype_build_plan_composition.R (which read the superseded v1 harmonized metadata).
# Panels: a metacell PCA; b per-dataset UMAP by class; c UMAP by dataset x age bin; d metacell composition by
# sample; e confident-cell composition by sample; f UMAP by age bin. One PNG + one PDF per panel.
suppressPackageStartupMessages({
  library(data.table); library(Matrix); library(ggplot2); library(patchwork)
})
options(stringsAsFactors = FALSE)

step_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce/step_02_preprocessing_filtering")
run_root <- dirname(step_dir)
v2_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
umap_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/data_input/standardized_seurat/05_cell_type_features/umap")
out_dir <- file.path(step_dir, "figures", "ED_Fig1_preprocessing_and_metacells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260707)

class_levels <- c("Granulosa", "Theca_Stroma", "Endothelial", "Pericyte_Smooth_Muscle", "Epithelial",
                  "Immune_Myeloid", "T_NK", "B_Plasma", "Oocyte_Germline", "Erythroid", "Not used")
class_colors <- c(Granulosa = "#E7298A", Theca_Stroma = "#B07AE8", Endothelial = "#1FA5C9",
                  Pericyte_Smooth_Muscle = "#6A3D9A", Epithelial = "#7FB800", Immune_Myeloid = "#E06A00",
                  T_NK = "#0067A5", B_Plasma = "#1B9E77", Oocyte_Germline = "#F2C14E", Erythroid = "#C0392B",
                  `Not used` = "grey80")
age_levels <- c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
age_labels <- c(young_adult = "Young adult", mature_adult = "Mature adult",
                advanced_maternal_age = "Advanced maternal age", post_reproductive_ageing = "Post-reproductive")
age_colors <- c(`Young adult` = "#4DAF4A", `Mature adult` = "#377EB8",
                `Advanced maternal age` = "#FF7F00", `Post-reproductive` = "#984EA3")
ds_short <- function(x) {
  x <- sub("^(20[0-9]{2})_[A-Za-z]+_([A-Za-z]+)_.*$", "\\2 \\1", x)
  sub("^(20[0-9]{2})_Kim_PRJNA863443.*$", "Kim 2025", x)
}
save_panel <- function(stem, p, w, h) {
  ggsave(file.path(out_dir, paste0(stem, ".png")), p, width = w, height = h, dpi = 200, limitsize = FALSE)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf, limitsize = FALSE)
}
theme_ed <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))

# ---- inputs -------------------------------------------------------------------------------------------------
cons <- fread(cmd = paste("zcat", shQuote(file.path(v2_dir, "04_consensus", "harmonized_annotation_consensus_v2.tsv.gz"))))
cons[, eligible := curation_status %in% c("confirmed", "probable") & annotation_confidence >= 0.70 &
       !annotation_label_broad %in% c("Unknown", "Erythroid")]
cons[, plot_class := fifelse(eligible, annotation_label_broad, "Not used")]
roles <- fread(file.path(run_root, "step_01_metadata", "scRNA_sample_role_metadata.tsv"))
roles <- roles[, .(dataset_id, sample_id = as.character(sample_id), model_role, reproductive_age_bin, age_weeks)]
cons[, sample_id := as.character(sample_id)]
cons <- merge(cons, roles, by = c("dataset_id", "sample_id"), all.x = TRUE)
cons[, age_bin := factor(age_labels[reproductive_age_bin], levels = age_labels)]
meta <- fread(file.path(step_dir, "filtered_metacell_metadata.tsv"))
counts <- readRDS(file.path(step_dir, "filtered_metacell_counts.rds"))
stopifnot(identical(colnames(counts), meta$metacell_id))

# ---- a: metacell PCA ----------------------------------------------------------------------------------------
lib <- Matrix::colSums(counts)
logcpm <- log1p(t(t(counts) / lib) * 1e6)
keep_g <- Matrix::rowSums(counts > 0) >= 0.1 * ncol(counts)
v <- apply(as.matrix(logcpm[keep_g, ]), 1, var)
hv <- names(sort(v, decreasing = TRUE))[seq_len(min(2000, length(v)))]
pc <- prcomp(t(as.matrix(logcpm[hv, ])), center = TRUE, scale. = TRUE, rank. = 5)
ve <- round(100 * pc$sdev[1:2]^2 / sum(pc$sdev^2), 1)
pca_dt <- data.table(meta[, .(metacell_id, dataset_id, harmonized_cell_class, model_role)], pc$x[, 1:2])
pca_dt[, dataset := ds_short(dataset_id)]
fwrite(pca_dt, file.path(out_dir, "edfig1a_metacell_pca_coordinates.tsv"), sep = "\t")
pa1 <- ggplot(pca_dt[sample(.N)], aes(PC1, PC2, colour = harmonized_cell_class)) + geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = class_colors, name = "Cell class") + theme_ed +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), title = "Coloured by cell class") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
pa2 <- ggplot(pca_dt[sample(.N)], aes(PC1, PC2, colour = dataset)) + geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_brewer(palette = "Dark2", name = "Dataset") + theme_ed +
  labs(x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)"), title = "Coloured by dataset") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
save_panel("edfig1a_metacell_exploratory_pca", (pa1 | pa2) +
             plot_annotation(title = sprintf("Metacell PCA (v2 labels; %s metacells, top %d variable genes)",
                                             format(nrow(pca_dt), big.mark = ","), length(hv))), 15, 6)

# ---- UMAP table ---------------------------------------------------------------------------------------------
umap <- rbindlist(lapply(unique(cons$dataset_id), function(ds) {
  u <- fread(cmd = paste("zcat", shQuote(file.path(umap_dir, paste0(ds, ".processed_umap.tsv.gz")))))
  setnames(u, tolower(names(u)))
  u[, .(cell_id, umap_1, umap_2)]
}))
um <- merge(cons[, .(cell_id, dataset_id, plot_class, age_bin, model_role)], umap, by = "cell_id")
um[, dataset := ds_short(dataset_id)]
um[, plot_class := factor(plot_class, levels = class_levels)]
cat("UMAP cells matched:", nrow(um), "of", nrow(cons), "\n")
um_plot <- um[sample(.N)][order(plot_class == "Not used", decreasing = TRUE)]  # grey drawn first

pb <- ggplot(um_plot, aes(umap_1, umap_2, colour = plot_class)) + geom_point(size = 0.25, alpha = 0.8, stroke = 0) +
  facet_wrap(~dataset, scales = "free", ncol = 4) + scale_colour_manual(values = class_colors, name = "v2 cell class", drop = TRUE) +
  theme_ed + labs(x = "UMAP 1", y = "UMAP 2", title = "v2 cell classes on standardized per-dataset UMAPs",
                  subtitle = "Grey: Unknown, Erythroid or below the clock confidence policy") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) + theme(legend.position = "bottom")
save_panel("edfig1b_v2_class_umap_by_dataset", pb, 16, 9.5)

umc <- um_plot[!is.na(age_bin) & model_role == "clock_candidate"]
pc_ <- ggplot(umc, aes(umap_1, umap_2, colour = plot_class)) + geom_point(size = 0.25, alpha = 0.8, stroke = 0) +
  facet_grid(dataset ~ age_bin, scales = "free") + scale_colour_manual(values = class_colors, name = "v2 cell class", drop = TRUE) +
  theme_ed + labs(x = "UMAP 1", y = "UMAP 2", title = "v2 cell classes by dataset and reproductive age bin (natural-ageing samples)") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), nrow = 2)) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0))
save_panel("edfig1c_v2_class_umap_by_dataset_age", pc_, 14, 20)

pf <- ggplot(um_plot[!is.na(age_bin)], aes(umap_1, umap_2, colour = age_bin)) + geom_point(size = 0.25, alpha = 0.8, stroke = 0) +
  facet_wrap(~dataset, scales = "free", ncol = 4) + scale_colour_manual(values = age_colors, name = "Age bin") +
  theme_ed + labs(x = "UMAP 1", y = "UMAP 2", title = "Per-dataset UMAPs coloured by reproductive age bin") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) + theme(legend.position = "bottom")
save_panel("edfig1f_umap_by_four_age_bins", pf, 16, 9.5)

# ---- d/e: composition by sample ------------------------------------------------------------------------------
comp_plot <- function(dt, value, ylab, title) {
  dt[, dataset := factor(ds_short(dataset_id))]
  ggplot(dt, aes(sample_id, frac, fill = cls)) + geom_col(width = 0.85) +
    facet_grid(~dataset, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = class_colors, name = "v2 cell class", drop = TRUE) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) + theme_ed +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6), legend.position = "bottom",
          strip.text = element_text(size = 7)) +
    labs(x = NULL, y = ylab, title = title)
}
md <- meta[, .(n = .N), by = .(dataset_id, sample_id = as.character(sample_id), cls = harmonized_cell_class)]
md[, frac := n / sum(n), by = .(dataset_id, sample_id)]
fwrite(md, file.path(out_dir, "edfig1d_metacell_composition_by_sample.tsv"), sep = "\t")
save_panel("edfig1d_v2_metacell_composition_by_sample",
           comp_plot(md, "frac", "Metacell fraction", "v2 metacell composition by sample (clock-candidate and projection-only samples)"), 22, 7)
ce <- cons[eligible == TRUE, .(n = .N), by = .(dataset_id, sample_id, cls = annotation_label_broad)]
ce[, frac := n / sum(n), by = .(dataset_id, sample_id)]
fwrite(ce, file.path(out_dir, "edfig1e_confident_cell_composition_by_sample.tsv"), sep = "\t")
save_panel("edfig1e_v2_confident_cell_composition_by_sample",
           comp_plot(ce, "frac", "Cell fraction", "v2 confident-cell composition by sample"), 22, 7)

writeLines(c(
  "# ED Fig. 1 — preprocessing and metacells (v2 labels)", "",
  "Built by `../../build_edfig1_preprocessing_and_metacells_v2.R` from the v2 consensus",
  "(`v2/04_consensus/harmonized_annotation_consensus_v2.tsv.gz`), the v2 metacells in this step folder and the",
  "standardized per-dataset UMAP sidecars. Grey cells in UMAPs are Unknown, Erythroid or below the clock policy",
  "(confirmed/probable, confidence >= 0.70).", "",
  "| Panel | File |", "|---|---|",
  "| a | `edfig1a_metacell_exploratory_pca` (+ coordinates TSV) |",
  "| b | `edfig1b_v2_class_umap_by_dataset` |",
  "| c | `edfig1c_v2_class_umap_by_dataset_age` |",
  "| d | `edfig1d_v2_metacell_composition_by_sample` (+ TSV) |",
  "| e | `edfig1e_v2_confident_cell_composition_by_sample` (+ TSV) |",
  "| f | `edfig1f_umap_by_four_age_bins` |", "",
  "Differences from 20260707: v2 labels; one PNG and one PDF per panel (no SVG/TIFF copies); panel a uses all v2",
  "metacells (clock-candidate and projection-only) instead of the superseded somatic pseudometacells."),
  file.path(out_dir, "README.md"))
message("ED Fig. 1 v2 written to ", out_dir)
