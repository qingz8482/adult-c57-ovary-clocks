#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grepl("^--file=", file_arg)][1])
script_dir <- dirname(normalizePath(file_arg, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
clock_dir <- normalizePath(
  Sys.getenv("OVARY_INTEGRATED_CLOCK_ROOT", unset = file.path(repo_root, "integrated_clocks", "shared")),
  mustWork = FALSE
)

projection_dir <- file.path(clock_dir, "results", "external_projection")
output_dir <- file.path(projection_dir, "GSE235021")
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")
log_dir <- file.path(output_dir, "logs")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pca_path <- file.path(projection_dir, "pca_projection_scores.tsv")
metadata_path <- file.path(clock_dir, "methods_tables", "04_external_projection_samples.tsv")
for (path in c(pca_path, metadata_path)) {
  if (!file.exists(path)) stop("Missing required input: ", path, call. = FALSE)
}

log_path <- file.path(log_dir, "GSE235021_pca_pc_contribution_heatmap.log")
log_lines <- character()
log_msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "  ", paste0(..., collapse = ""))
  log_lines <<- c(log_lines, line)
  message(line)
}

log_msg("Reading frozen PCA-clock projection scores and GSE235021 metadata.")
pca_scores <- fread(pca_path, sep = "\t", na.strings = c("", "NA", "NaN"))
metadata <- fread(metadata_path, sep = "\t", na.strings = c("", "NA", "NaN"))

required_meta <- c("sample_key", "dataset_folder", "modality", "age_week", "treatment", "sample_title")
missing_meta <- setdiff(required_meta, names(metadata))
if (length(missing_meta) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

pc_contribution_cols <- grep("^PC[0-9]+_contribution_week$", names(pca_scores), value = TRUE)
if (length(pc_contribution_cols) == 0L) {
  stop("No frozen PCA-clock PC contribution columns were found.", call. = FALSE)
}
pc_number <- as.integer(sub("^PC([0-9]+)_.*$", "\\1", pc_contribution_cols))
pc_contribution_cols <- pc_contribution_cols[order(pc_number)]

gse_meta <- unique(metadata[dataset_folder == "GSE235021", ..required_meta])
if (nrow(gse_meta) != 35L || uniqueN(gse_meta$sample_key) != 35L) {
  stop(
    "Expected 35 unique GSE235021 samples (18 bulk + 17 pseudobulk); observed ",
    nrow(gse_meta), " rows and ", uniqueN(gse_meta$sample_key), " unique keys.",
    call. = FALSE
  )
}

gse_scores <- pca_scores[sample_key %chin% gse_meta$sample_key]
if (nrow(gse_scores) != 35L || uniqueN(gse_scores$sample_key) != 35L) {
  stop(
    "Expected one frozen PCA-clock score row for each of 35 GSE235021 samples; observed ",
    nrow(gse_scores), " rows and ", uniqueN(gse_scores$sample_key), " unique keys.",
    call. = FALSE
  )
}

gse_meta_join <- gse_meta[, setdiff(names(gse_meta), "dataset_folder"), with = FALSE]
gse_scores <- merge(gse_scores, gse_meta_join, by = "sample_key", all.x = TRUE, all.y = FALSE)
if (anyNA(gse_scores$modality) || anyNA(gse_scores$treatment) || anyNA(gse_scores$age_week)) {
  stop("GSE235021 PCA scores did not join completely to modality, diet, and age metadata.", call. = FALSE)
}

expected_modalities <- c("bulk_rnaseq", "scrna_pseudobulk")
modality_counts <- gse_scores[, .(n_samples = uniqueN(sample_key)), by = modality]
if (!setequal(modality_counts$modality, expected_modalities) ||
    modality_counts[modality == "bulk_rnaseq", n_samples] != 18L ||
    modality_counts[modality == "scrna_pseudobulk", n_samples] != 17L) {
  stop("Expected 18 bulk RNA-seq and 17 scRNA-derived pseudobulk samples.", call. = FALSE)
}

gse_scores[, age_month := fifelse(age_week < 32, 6L, 9L)]
gse_scores[, diet_label := fcase(
  treatment == "polished rice", "Polished rice",
  treatment == "black rice", "Black rice",
  treatment == "white brown rice", "White/brown rice",
  default = NA_character_
)]
if (anyNA(gse_scores$diet_label)) {
  stop("Unexpected GSE235021 diet label: ", paste(unique(gse_scores[is.na(diet_label), treatment]), collapse = ", "), call. = FALSE)
}
gse_scores[, modality_label := fcase(
  modality == "bulk_rnaseq", "Bulk RNA-seq",
  modality == "scrna_pseudobulk", "scRNA pseudobulk",
  default = modality
)]

gse_scores[, display_group := paste0(age_month, " months · ", diet_label, " · ", modality_label)]
diet_order <- c("Polished rice", "Black rice", "White/brown rice")
modality_order <- c("Bulk RNA-seq", "scRNA pseudobulk")
display_levels <- CJ(
  modality_label = modality_order,
  age_month = c(6L, 9L),
  diet_label = diet_order,
  sorted = FALSE
)[, paste0(age_month, " months · ", diet_label, " · ", modality_label)]
gse_scores[, display_group := factor(display_group, levels = rev(display_levels))]

id_vars <- c(
  "sample_key", "model_family", "selector", "dataset_folder", "modality",
  "modality_label", "age_week", "age_month", "treatment", "diet_label",
  "sample_title", "display_group"
)
contribution_long <- melt(
  gse_scores,
  id.vars = id_vars,
  measure.vars = pc_contribution_cols,
  variable.name = "pc",
  value.name = "pc_contribution_week"
)
contribution_long[, pc := sub("_contribution_week$", "", pc)]
contribution_long[, pc_number := as.integer(sub("^PC", "", pc))]
setorder(contribution_long, sample_key, pc_number)

if (any(!is.finite(contribution_long$pc_contribution_week))) {
  stop("Non-finite PC contributions are present.", call. = FALSE)
}

group_summary <- contribution_long[, .(
  mean_pc_contribution_week = mean(pc_contribution_week),
  sd_pc_contribution_week = sd(pc_contribution_week),
  n_samples = uniqueN(sample_key)
), by = .(age_month, diet_label, modality, modality_label, display_group, pc, pc_number)]
group_summary[, pc := factor(pc, levels = paste0("PC", sort(unique(pc_number))))]
setorder(group_summary, display_group, pc_number)

fwrite(
  gse_scores[order(age_month, diet_label, modality, sample_key)],
  file.path(output_dir, "gse235021_pca_scores.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  contribution_long,
  file.path(output_dir, "gse235021_pca_pc_contributions.tsv"),
  sep = "\t",
  na = "NA"
)
fwrite(
  group_summary,
  file.path(table_dir, "gse235021_pca_pc_contribution_group_summary.tsv"),
  sep = "\t",
  na = "NA"
)

theme_gse <- function() {
  theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      strip.text = element_text(face = "bold")
    )
}

fill_limits <- range(group_summary$mean_pc_contribution_week)

p <- ggplot(group_summary, aes(pc, display_group, fill = mean_pc_contribution_week)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = fill_limits
  ) +
  labs(
    x = "Frozen selected PC",
    y = "Age, diet, and modality",
    fill = "Weeks"
  ) +
  theme_gse()

stem <- "08_gse235021_pca_pc_contribution_heatmap_modality_blocked"
pdf_path <- file.path(plot_dir, paste0(stem, ".pdf"))
png_path <- file.path(plot_dir, paste0(stem, ".png"))
plot_width_in <- 9
# Geometry measured from the original 9 x 7 inch PRJNA863443 ggplot:
# 5.332861 inch panel / 14 rows = 0.3809186 inch (9.675334 mm) per tile.
# With unchanged axes and legend, 12 GSE235021 rows require a 6.238163 inch page.
plot_height_in <- 6.238163
reference_panel_height_in <- 5.332861
reference_n_rows <- 14L
reference_tile_height_in <- reference_panel_height_in / reference_n_rows
gse_panel_height_in <- reference_tile_height_in * 12L
gse_tile_height_in <- gse_panel_height_in / 12L

ggsave(pdf_path, p, width = plot_width_in, height = plot_height_in, limitsize = FALSE)
ggsave(png_path, p, width = plot_width_in, height = plot_height_in, dpi = 300, limitsize = FALSE)

bulk_summary <- copy(group_summary[modality == "bulk_rnaseq"])
bulk_display_levels <- CJ(
  age_month = c(6L, 9L),
  diet_label = diet_order,
  sorted = FALSE
)[, paste0(age_month, " months · ", diet_label)]
bulk_summary[, bulk_display_group := factor(
  paste0(age_month, " months · ", diet_label),
  levels = rev(bulk_display_levels)
)]

p_bulk <- ggplot(bulk_summary, aes(pc, bulk_display_group, fill = mean_pc_contribution_week)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#3B6EA8",
    mid = "white",
    high = "#C44E52",
    midpoint = 0,
    limits = fill_limits
  ) +
  labs(
    x = "Frozen selected PC",
    y = "Age and diet",
    fill = "Weeks"
  ) +
  theme_gse()

bulk_stem <- "09_gse235021_pca_pc_contribution_heatmap_bulk_only"
bulk_pdf_path <- file.path(plot_dir, paste0(bulk_stem, ".pdf"))
bulk_png_path <- file.path(plot_dir, paste0(bulk_stem, ".png"))
bulk_n_rows <- uniqueN(bulk_summary$bulk_display_group)
fixed_plot_overhead_in <- plot_height_in - gse_panel_height_in
bulk_panel_height_in <- reference_tile_height_in * bulk_n_rows
bulk_plot_height_in <- fixed_plot_overhead_in + bulk_panel_height_in
bulk_tile_height_in <- bulk_panel_height_in / bulk_n_rows

ggsave(
  bulk_pdf_path,
  p_bulk,
  width = plot_width_in,
  height = bulk_plot_height_in,
  limitsize = FALSE
)
ggsave(
  bulk_png_path,
  p_bulk,
  width = plot_width_in,
  height = bulk_plot_height_in,
  dpi = 300,
  limitsize = FALSE
)

geometry <- data.table(
  figure = c("PRJNA863443 reference", "GSE235021 matched", "GSE235021 bulk only"),
  page_width_in = c(9, plot_width_in, plot_width_in),
  page_height_in = c(7, plot_height_in, bulk_plot_height_in),
  panel_height_in = c(
    reference_panel_height_in,
    gse_panel_height_in,
    bulk_panel_height_in
  ),
  n_heatmap_rows = c(reference_n_rows, 12L, bulk_n_rows),
  tile_height_in = c(
    reference_tile_height_in,
    gse_tile_height_in,
    bulk_tile_height_in
  ),
  tile_height_mm = c(
    reference_tile_height_in,
    gse_tile_height_in,
    bulk_tile_height_in
  ) * 25.4
)
fwrite(
  geometry,
  file.path(table_dir, "gse235021_pca_pc_contribution_geometry.tsv"),
  sep = "\t"
)

validation <- rbindlist(list(
  data.table(
    check = "unique_GSE235021_samples",
    expected = "35",
    observed = as.character(uniqueN(gse_scores$sample_key)),
    status = ifelse(uniqueN(gse_scores$sample_key) == 35L, "pass", "fail")
  ),
  data.table(
    check = "bulk_RNAseq_samples",
    expected = "18",
    observed = as.character(modality_counts[modality == "bulk_rnaseq", n_samples]),
    status = ifelse(modality_counts[modality == "bulk_rnaseq", n_samples] == 18L, "pass", "fail")
  ),
  data.table(
    check = "scRNA_pseudobulk_samples",
    expected = "17",
    observed = as.character(modality_counts[modality == "scrna_pseudobulk", n_samples]),
    status = ifelse(modality_counts[modality == "scrna_pseudobulk", n_samples] == 17L, "pass", "fail")
  ),
  data.table(
    check = "finite_PC_contributions",
    expected = as.character(nrow(contribution_long)),
    observed = as.character(sum(is.finite(contribution_long$pc_contribution_week))),
    status = ifelse(all(is.finite(contribution_long$pc_contribution_week)), "pass", "fail")
  ),
  data.table(
    check = "complete_age_diet_modality_groups",
    expected = "12",
    observed = as.character(uniqueN(group_summary$display_group)),
    status = ifelse(uniqueN(group_summary$display_group) == 12L, "pass", "fail")
  ),
  data.table(
    check = "tile_height_matches_PRJNA863443",
    expected = sprintf("%.6f mm", reference_tile_height_in * 25.4),
    observed = sprintf("%.6f mm", gse_tile_height_in * 25.4),
    status = ifelse(abs(gse_tile_height_in - reference_tile_height_in) < 1e-10, "pass", "fail")
  ),
  data.table(
    check = "bulk_only_groups",
    expected = "6",
    observed = as.character(bulk_n_rows),
    status = ifelse(
      bulk_n_rows == 6L &&
        uniqueN(bulk_summary$modality) == 1L &&
        unique(bulk_summary$modality) == "bulk_rnaseq",
      "pass",
      "fail"
    )
  ),
  data.table(
    check = "bulk_tile_height_matches_PRJNA863443",
    expected = sprintf("%.6f mm", reference_tile_height_in * 25.4),
    observed = sprintf("%.6f mm", bulk_tile_height_in * 25.4),
    status = ifelse(abs(bulk_tile_height_in - reference_tile_height_in) < 1e-10, "pass", "fail")
  ),
  data.table(
    check = "nonempty_PDF_PNG_exports",
    expected = "4",
    observed = as.character(sum(
      file.exists(c(pdf_path, png_path, bulk_pdf_path, bulk_png_path)) &
        file.info(c(pdf_path, png_path, bulk_pdf_path, bulk_png_path))$size > 0
    )),
    status = ifelse(
      all(file.exists(c(pdf_path, png_path, bulk_pdf_path, bulk_png_path))) &&
        all(file.info(c(pdf_path, png_path, bulk_pdf_path, bulk_png_path))$size > 0),
      "pass",
      "fail"
    )
  )
))
fwrite(
  validation,
  file.path(table_dir, "gse235021_pca_pc_contribution_validation.tsv"),
  sep = "\t"
)
if (any(validation$status != "pass")) {
  stop("One or more GSE235021 PCA contribution validation checks failed.", call. = FALSE)
}

readme <- c(
  "# GSE235021 frozen PCA-clock PC contribution heatmap",
  "",
  "This folder mirrors the PRJNA863443 PC-contribution heatmap release for GSE235021.",
  "The frozen PCA/glmnet projection scores are not refitted. Sample-level calibrated",
  "PC contributions are summarized within age × rice diet × modality. Separate",
  "modality-blocked and bulk-only heatmaps are provided.",
  "",
  "## Cohort",
  "",
  "- 35 samples: 18 bulk RNA-seq and 17 scRNA-derived pseudobulk.",
  "- Ages: 6 and 9 months.",
  "- Diets: polished rice, black rice, and white/brown rice.",
  "- Seven selected frozen PCs are shown.",
  "",
  "## Interpretation",
  "",
  "Each tile is the mean calibrated contribution of one frozen PCA-clock component",
  "to predicted age, in weeks. The figure is descriptive; assay modality is partly",
  "confounded with library preparation and should not be interpreted as a biological",
  "diet effect.",
  "",
  "## Files",
  "",
  "- `plots/07_gse235021_pca_pc_contribution_heatmap.pdf` and `.png`: original interleaved-modality figure pair.",
  paste0("- `plots/", stem, ".pdf` and `.png`: bulk-first, scRNA-pseudobulk-second figure pair."),
  paste0("- `plots/", bulk_stem, ".pdf` and `.png`: bulk RNA-seq-only figure pair."),
  "- `gse235021_pca_scores.tsv`: GSE235021 frozen PCA projection scores.",
  "- `gse235021_pca_pc_contributions.tsv`: sample-level long contribution table.",
  "- `tables/gse235021_pca_pc_contribution_group_summary.tsv`: heatmap source data.",
  "- `tables/gse235021_pca_pc_contribution_geometry.tsv`: reference and matched tile dimensions.",
  "- `tables/gse235021_pca_pc_contribution_validation.tsv`: automated validation."
)
writeLines(readme, file.path(output_dir, "README.md"))

log_msg(
  "Completed GSE235021 PCA contribution heatmap: ",
  uniqueN(gse_scores$sample_key), " samples; ",
  modality_counts[modality == "bulk_rnaseq", n_samples], " bulk; ",
  modality_counts[modality == "scrna_pseudobulk", n_samples], " pseudobulk; ",
  length(pc_contribution_cols), " PCs; validation ",
  sum(validation$status == "pass"), "/", nrow(validation), " passed."
)
writeLines(log_lines, log_path)
