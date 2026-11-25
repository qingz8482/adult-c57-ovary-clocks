#!/usr/bin/env Rscript
# Horizontal (transposed) view of A2: signature score per dataset x v2 class.
# Classes are column facets, datasets on x within each class, signatures on y.
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
v2_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
ev_dir <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce/step_02_preprocessing_filtering/figures_cell annotation/annotation_evidence_v2")

classes <- c("Granulosa", "Theca_Stroma", "Endothelial", "Pericyte_Smooth_Muscle", "Epithelial",
             "Immune_Myeloid", "T_NK", "B_Plasma", "Oocyte_Germline")
class_lab <- c(Granulosa = "Granulosa", Theca_Stroma = "Theca/stroma", Endothelial = "Endothelial",
               Pericyte_Smooth_Muscle = "Pericyte/SMC", Epithelial = "Epithelial", Immune_Myeloid = "Myeloid",
               T_NK = "T/NK", B_Plasma = "B/plasma", Oocyte_Germline = "Oocyte")
sig_lab <- c(granulosa = "Granulosa", luteal = "Luteal", theca = "Theca", pericyte = "Pericyte", smooth_muscle = "SMC",
             stroma = "Stroma", endothelial = "Endothelial", epithelial = "Epithelial", oocyte = "Oocyte",
             myeloid = "Myeloid", neutrophil = "Neutrophil", t_nk = "T/NK", b_plasma = "B/plasma", erythroid = "Erythroid")
ds_short <- function(x) {
  x <- sub("^(20[0-9]{2})_[A-Za-z]+_([A-Za-z]+)_.*$", "\\2 \\1", x)
  sub("^(20[0-9]{2})_Kim_PRJNA863443.*$", "Kim 2025", x)
}

panel <- fread(file.path(v2_dir, "01_marker_panel_v2.tsv"))
sd_ <- fread(file.path(ev_dir, "A_signature_score_by_dataset_class.tsv"))
sd_[, signature := factor(sig_lab[signature], levels = rev(sig_lab[unique(panel$signature)]))]
sd_[, cell_class := factor(class_lab[cell_class], levels = class_lab[classes])]
sd_[, dataset := ds_short(dataset_id)]

p <- ggplot(sd_, aes(dataset, signature, fill = score)) + geom_tile(colour = "white", linewidth = 0.2) +
  facet_grid(. ~ cell_class, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "grey96", high = "#B40426", midpoint = 0,
                       limits = c(-1.5, 2.5), oob = scales::squish, name = "Signature\nscore (z)") +
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(text = element_text(family = "Arial", size = 12, colour = "black"),
        axis.text = element_text(family = "Arial", size = 12, colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1),
        plot.title = element_text(family = "Arial", size = 12, face = "bold"),
        strip.text.x = element_text(family = "Arial", size = 12),
        legend.title = element_text(family = "Arial", size = 12), legend.text = element_text(family = "Arial", size = 12),
        panel.grid = element_blank(), strip.background = element_rect(fill = "grey92", colour = NA),
        panel.spacing.x = unit(3, "pt"), legend.key.height = unit(12, "pt"),
        plot.margin = margin(10, 10, 10, 10, "pt")) +
  labs(x = NULL, y = "Marker signature", title = "Signature scores per dataset and class")

out <- file.path(ev_dir, "A2_signature_score_by_dataset_class_horizontal")
# width matches the Illustrator artboard (1570.5369 pt); height kept short so tiles are wider than tall
w_in <- 1570.5369 / 72; h_in <- 4.6
ggsave(paste0(out, ".pdf"), p, width = w_in, height = h_in, device = grDevices::cairo_pdf)
ggsave(paste0(out, ".png"), p, width = w_in, height = h_in, dpi = 150)
message("written: ", out, ".{pdf,png}")
