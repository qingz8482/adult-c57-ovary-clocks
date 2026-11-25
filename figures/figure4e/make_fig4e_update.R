suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)})

args <- commandArgs(TRUE)
in_tsv <- args[1]; out_prefix <- args[2]

d <- read_tsv(in_tsv, show_col_types = FALSE)

# panel order / names as in original Fig 4e
ct_map <- c(B_Plasma = "B cells", Endothelial = "Endothelial", Epithelial = "Epithelial",
            Theca_Stroma = "Theca / Stroma", Granulosa = "Granulosa", Immune_Myeloid = "Immune Myeloid",
            Pericyte_Smooth_Muscle = "Smooth Muscle", T_NK = "T / NK")
ds_map <- c(
  "2022_eLife_Morris_cycling_murine_ovary_atlas" = "2022 Morris",
  "2024_Cell_Winkler_FRT_scRNA_spatial" = "2024 Winkler",
  "2024_LifeMedicine_Jiang_whole_life_mouse_ovary" = "2024 Jiang",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "2024 Isola",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary" = "2025 Zhang",
  "2025_CommBio_Zhao_developing_ovary_scRNA_spatial" = "2025 Zhao",
  "2026_CommBio_Galligos_ovarian_inflammaging" = "2026 Galligos",
  "2025_Kim_PRJNA863443_mouse_menopause_models" = "2025 Kim")
ds_pal <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#7F3C8D", "#11A579", "#E69F00", "#56B4E9"), ds_map)
age_levels <- c("young_adult", "mature_adult", "advanced_maternal_age", "post_reproductive_ageing")
age_labels <- c("Young adult", "Mature adult", "Advanced maternal age", "Post−reproductive")
age_shapes <- c(16, 17, 15, 18)

d <- d %>% mutate(
  cell = factor(ct_map[harmonized_cell_class], levels = ct_map),
  dataset = factor(ds_map[dataset_id], levels = ds_map),
  age_bin = factor(reproductive_age_bin, levels = age_levels, labels = age_labels))

# stats from the release table's facet labels (n, MAE, slope, R2 from step_05 CV evaluation)
lab <- d %>% group_by(cell, facet_label_tsv) %>%
  summarise(x = min(age_weeks) - 0.01 * diff(range(age_weeks)),
            y = max(predicted_age_weeks) + 0.025 * diff(range(predicted_age_weeks)), .groups = "drop") %>%
  mutate(label = sub("^[^;]+; ", "", facet_label_tsv),
         label = gsub(", ", "\n", gsub("; ", "\n", label)))
print(lab %>% select(cell, label))

pt <- ggplot2::.pt
p <- ggplot(d, aes(age_weeks, predicted_age_weeks, colour = dataset, shape = age_bin)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.5) +
  geom_point(alpha = 0.55, size = 2.0, stroke = 0) +
  geom_text(data = lab, aes(x = x, y = y, label = label), inherit.aes = FALSE,
            hjust = 0, vjust = 1, family = "Arial", fontface = "italic",
            size = 12 / pt, lineheight = 0.92) +
  # extend B cells y range below 0 so all clocks share the same y-axis floor
  geom_blank(data = data.frame(cell = factor("B cells", levels = ct_map), predicted_age_weeks = -10),
             aes(y = predicted_age_weeks), inherit.aes = FALSE) +
  facet_wrap(~ cell, nrow = 2, scales = "free") +
  scale_x_continuous(breaks = c(25, 50, 75)) +
  scale_y_continuous(breaks = scales::breaks_width(25), expand = expansion(mult = c(0.05, 0.05))) +
  scale_colour_manual(values = ds_pal, drop = FALSE) +
  scale_shape_manual(values = age_shapes, drop = FALSE) +
  guides(colour = guide_legend(title = "Dataset", nrow = 1, title.position = "top", order = 1,
                               override.aes = list(shape = 16, size = 1, alpha = 0.55)),
         shape = guide_legend(title = "Age bin", nrow = 1, title.position = "top", order = 2,
                              override.aes = list(colour = "grey35", size = 1, alpha = 0.8))) +
  labs(x = "Observed age (weeks)", y = "Predicted age (weeks)") +
  theme_bw(base_family = "Arial") +
  theme(
    aspect.ratio = 1,
    text = element_text(colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    strip.text = element_text(size = 12, face = "bold", margin = margin(3, 0, 3, 0)),
    strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.8),
    panel.border = element_rect(colour = "black", linewidth = 0.8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.5),
    panel.spacing.x = unit(52, "pt"), panel.spacing.y = unit(14, "pt"),
    legend.position = "bottom", legend.box = "horizontal", legend.justification = "left",
    legend.text = element_text(size = 12),
    legend.key = element_blank(), legend.key.size = unit(8, "pt"),
    legend.key.spacing.x = unit(10, "pt"), legend.spacing.x = unit(100, "pt"),
    legend.title.position = "top", legend.title = element_text(size = 12, face = "bold", margin = margin(b = 1)),
    legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 2)),
    plot.margin = margin(2, 2, 2, 2))

# original panel e occupies ~1250 x 605 pt (+ y-axis title) on the Figure 4 artboard
W <- 1272 / 72; H <- 605 / 72
ggsave(paste0(out_prefix, ".pdf"), p, width = W, height = H, device = cairo_pdf)
ggsave(paste0(out_prefix, ".png"), p, width = W, height = H, dpi = 150)
