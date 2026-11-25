#!/usr/bin/env Rscript
# Redraw the 8-cell-type ElasticNet coefficient-ranked GSEA NES panels sized for
# Supplementary Fig. 5 panel e (artboard 1570.54 x 996.42 pt; panel e area below panels a-d).
# All text Arial 11 pt; every term labelled with nominal GSEA p (no FDR q mixing). Reads the plot data written by organize_scRNA_metacell_cv_figure_package.R.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

release_root <- "step_08_figures_tables_release"
in_tsv <- file.path(release_root, "plot_data", "04_gsea",
                    "elasticnet_coefficient_gsea_nes_panels_by_celltype_plot_data.tsv")
out_stem <- file.path(release_root, "figures", "04_gsea",
                      "elasticnet_coefficient_gsea_nes_panels_by_celltype_SF5_artboard")

font_family <- "Arial"
font_pt <- 11
width_pt <- 1560
height_pt <- 690

# Display names matching SF5 panel d
celltype_labels <- c(
  B_Plasma = "B Plasma",
  Endothelial = "Endothelial",
  Epithelial = "Epithelial",
  Granulosa = "Granulosa",
  Immune_Myeloid = "Immune Myeloid",
  Pericyte_Smooth_Muscle = "Smooth Muscle",
  T_NK = "T / NK",
  Theca_Stroma = "Theca / Stroma"
)

# Abbreviations so each term fits on one line at 11 pt
abbreviate_term <- function(x) {
  replacements <- c(
    "^Adaptive Immune Response Based on Somatic Recombination of Immune Receptors Built from Immunoglobulin Superfamily Domains$" =
      "Adaptive Immune Resp. (Somatic Recomb.)",
    "Negative Regulation" = "Neg. Reg.",
    "Positive Regulation" = "Pos. Reg.",
    "^Regulation" = "Reg.",
    "Multicellular Organismal" = "Multicellular Org.",
    "Cell Population Proliferation" = "Cell Proliferation",
    "Anatomical Structure Formation Involved in Morphogenesis" = "Anatomical Structure Formation in Morphogenesis",
    "Gtpases" = "GTPases",
    "Rhobtb3" = "RHOBTB3",
    "^t Cell" = "T Cell",
    "^b Cell" = "B Cell"
  )
  for (pattern in names(replacements)) x <- gsub(pattern, replacements[[pattern]], x)
  x
}

gsea <- fread(in_tsv) %>%
  as_tibble() %>%
  mutate(
    harmonized_cell_class = factor(harmonized_cell_class, levels = unique(harmonized_cell_class)),
    term_display = if_else(is_padding_row, "", abbreviate_term(term_display)),
    # Nominal p for all terms (same format as run_scRNA_metacell_elasticnet_coefficient_gsea.R)
    p_label = if_else(p_value < 1e-4, paste0("p=", format(p_value, scientific = TRUE, digits = 2)),
                      paste0("p=", sprintf("%.3f", p_value))),
    term_plot_id = factor(term_plot_id, levels = unique(term_plot_id))
  )
term_labels <- setNames(gsea$term_display, as.character(gsea$term_plot_id))
bars <- gsea %>% filter(!is_padding_row, !is_status_row)

p <- ggplot(gsea, aes(NES, term_plot_id)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_col(data = bars, aes(fill = NES_direction_display), width = 0.72) +
  # q-value printed on the empty side of zero, opposite each bar
  geom_text(
    data = bars,
    aes(x = if_else(NES >= 0, -0.06, 0.06), label = p_label, hjust = if_else(NES >= 0, 1, 0)),
    size = font_pt / .pt, family = font_family, color = "grey20"
  ) +
  facet_wrap(~ harmonized_cell_class, scales = "free_y", ncol = 4, drop = FALSE,
             labeller = as_labeller(celltype_labels)) +
  scale_y_discrete(labels = term_labels) +
  scale_x_continuous(limits = c(-2.15, 2.15), breaks = -2:2, expand = c(0, 0)) +
  scale_fill_manual(values = c(
    "Positive NES = older-predictive clock-gene programs" = "#B2182B",
    "Negative NES = younger-predictive clock-gene programs" = "#2166AC"
  )) +
  theme_bw(base_size = font_pt, base_family = font_family) +
  theme(
    text = element_text(size = font_pt, family = font_family, color = "black"),
    axis.text = element_text(size = font_pt, color = "black"),
    axis.title = element_text(size = font_pt),
    strip.text = element_text(size = font_pt, face = "bold", margin = margin(2, 0, 2, 0)),
    strip.background = element_rect(fill = "grey95", color = "grey70", linewidth = 0.25),
    legend.text = element_text(size = font_pt),
    legend.title = element_text(size = font_pt, face = "bold"),
    legend.position = "bottom",
    legend.key.size = unit(10, "pt"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(4, "pt"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
    panel.spacing.x = unit(6, "pt"),
    panel.spacing.y = unit(6, "pt"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(margin = margin(r = 2)),
    plot.margin = margin(2, 4, 2, 2)
  ) +
  labs(x = "NES", y = NULL, fill = "Coefficient-ranked enrichment direction")

ggsave(paste0(out_stem, ".pdf"), p, width = width_pt / 72, height = height_pt / 72,
       units = "in", device = cairo_pdf)
ggsave(paste0(out_stem, ".png"), p, width = width_pt / 72, height = height_pt / 72,
       units = "in", dpi = 300, device = ragg::agg_png)
message("Wrote ", out_stem, ".{pdf,png}")
