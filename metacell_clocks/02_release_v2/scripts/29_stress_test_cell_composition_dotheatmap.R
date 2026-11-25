#!/usr/bin/env Rscript

# Composition-aware view of cell-class changes per stress-test contrast (input: script 28's per-library
# counts in the 8 clock classes).
#   Effect (colour): CLR log2 fold change = case - control difference in centred log2-ratio composition
#     (log2 proportion minus the library's mean log2 proportion), which removes most of the closure
#     effect where one class's loss inflates every other class's share.
#   Test (stars): limma on the same CLR values (lmFit + eBayes(robust = TRUE)), BH within each contrast
#     across the 8 classes, so stars and colour use one composition-aware scale.
#   Sensitivity (table only): propeller (speckle; Phipson et al. 2022) re-implemented with limma, since
#     speckle is not installed: logit-transformed proportions (pseudocount 0.5, as getTransformedProps).
#     It is not closure-aware, so a class can reach significance only because another class shrank.
#   One model per cohort:
#     VCD   ~ 0 + group (condition x injection age x post-injection days) + replicate (pairs CTL/VCD _1 and _2)
#     Foxl2 ~ 0 + group (genotype x age)
#     Rice  ~ 0 + group (diet x age); contrasts vs polished rice of the same age
#   Dot size: mean percentage of the class in the matched controls.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
})

run_root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/reproduce")
comp_dir <- file.path(run_root, "step_08_figures_tables_release/figures/05_stress_test_projection/cell composition")
comp <- fread(file.path(comp_dir, "stress_test_library_cell_class_composition.tsv"))

classes <- c("Theca_Stroma", "Granulosa", "Endothelial", "Epithelial", "T_NK",
             "Pericyte_Smooth_Muscle", "Immune_Myeloid", "B_Plasma")
class_lab <- c(Theca_Stroma = "Theca / Stroma", Granulosa = "Granulosa", Endothelial = "Endothelial",
               Epithelial = "Epithelial", T_NK = "T / NK", Pericyte_Smooth_Muscle = "Smooth Muscle",
               Immune_Myeloid = "Immune / Myeloid", B_Plasma = "B Plasma")

counts <- as.matrix(dcast(comp, sample_label ~ harmonized_cell_class, value.var = "n_cells"), rownames = "sample_label")[, classes]
pseudo <- (counts + 0.5) / rowSums(counts + 0.5)
logit_prop <- t(log(pseudo / (1 - pseudo)))                       # classes x libraries
clr <- t(log2(pseudo) - rowMeans(log2(pseudo)))                   # classes x libraries
pct <- t(100 * counts / rowSums(counts))

lib <- unique(comp[, .(sample_label, cohort)])
lib[, `:=`(
  family = fcase(grepl("^(CTL|VCD)_", sample_label), "VCD", grepl("^Foxl2_", sample_label), "Foxl2", default = "Rice diet"),
  group = fcase(
    grepl("^(CTL|VCD)_", sample_label), sub("_\\d+$", "", sample_label),
    grepl("^Foxl2_", sample_label), sub("_\\d+$", "", sample_label),
    default = sub("-\\d+$", "", sample_label)),
  replicate = sub("^.*[_-](\\d+)$", "\\1", sample_label)
)]

fit_cohort <- function(fam, contrasts, block_replicate = FALSE) {
  L <- lib[family == fam]
  grp <- factor(make.names(L$group))
  design <- if (block_replicate) model.matrix(~ 0 + grp + factor(L$replicate)) else model.matrix(~ 0 + grp)
  colnames(design) <- sub("^grp", "", colnames(design))
  colnames(design) <- make.names(colnames(design))
  cm <- makeContrasts(contrasts = unname(contrasts), levels = design)
  colnames(cm) <- names(contrasts)
  run <- function(y) eBayes(contrasts.fit(lmFit(y[, L$sample_label], design), cm), robust = TRUE)
  f_logit <- run(logit_prop); f_clr <- run(clr)
  rbindlist(lapply(names(contrasts), function(cn) {
    tt <- topTable(f_logit, coef = cn, number = Inf, sort.by = "none")
    tc <- topTable(f_clr, coef = cn, number = Inf, sort.by = "none")[rownames(tt), ]
    parts <- strsplit(contrasts[[cn]], " - ", fixed = TRUE)[[1]]
    case_libs <- L[make.names(group) == parts[1], sample_label]
    ctrl_libs <- L[make.names(group) == parts[2], sample_label]
    data.table(family = fam, contrast = cn, harmonized_cell_class = rownames(tt),
               clr_log2fc = tc$logFC, p_value = tc$P.Value, fdr = p.adjust(tc$P.Value, "BH"),
               propeller_logit_diff = tt$logFC, propeller_p_value = tt$P.Value,
               propeller_fdr = p.adjust(tt$P.Value, "BH"),
               control_percent = rowMeans(pct[rownames(tt), ctrl_libs, drop = FALSE]),
               case_percent = rowMeans(pct[rownames(tt), case_libs, drop = FALSE]),
               n_case_libraries = length(case_libs), n_control_libraries = length(ctrl_libs))
  }))
}

vcd_ct <- c("VCD 3m 30d" = "VCD_3m_30d - CTL_3m_30d", "VCD 3m 90d" = "VCD_3m_90d - CTL_3m_90d",
            "VCD 10m 30d" = "VCD_10m_30d - CTL_10m_30d", "VCD 10m 90d" = "VCD_10m_90d - CTL_10m_90d")
foxl2_ct <- c("Foxl2 4 months" = "Foxl2_het_young - Foxl2_wt_young", "Foxl2 9 months" = "Foxl2_het_old - Foxl2_wt_old")
rice_ct <- c("Rice 6 mo: Black" = "SCB.6M - SCP.6M", "Rice 6 mo: White/brown" = "SCW.6M - SCP.6M",
             "Rice 9 mo: Black" = "SCB.9M - SCP.9M", "Rice 9 mo: White/brown" = "SCW.9M - SCP.9M")
res <- rbind(fit_cohort("VCD", vcd_ct, block_replicate = TRUE), fit_cohort("Foxl2", foxl2_ct),
             fit_cohort("Rice diet", rice_ct))
res[, pct_point_change := case_percent - control_percent]
res[, sig := fcase(fdr < 0.001, "***", fdr < 0.01, "**", fdr < 0.05, "*", default = "")]
fwrite(res[order(match(contrast, c(names(vcd_ct), names(foxl2_ct), names(rice_ct))), match(harmonized_cell_class, classes))],
       file.path(comp_dir, "stress_test_cell_class_composition_propeller_clr.tsv"), sep = "\t")

# ---------------------------------------------------------------- dot-heatmap
res[, `:=`(contrast = factor(contrast, levels = c(names(vcd_ct), names(foxl2_ct), names(rice_ct))),
           harmonized_cell_class = factor(harmonized_cell_class, levels = rev(classes)),
           family = factor(family, levels = c("VCD", "Foxl2", "Rice diet")))]
lim <- 3.5   # B_Plasma in rice 9 mo (0% in controls) reaches |4| through the pseudocount; squished to the end
p <- ggplot(res, aes(contrast, harmonized_cell_class)) +
  geom_point(aes(size = control_percent, fill = clr_log2fc), shape = 21, colour = "grey30", stroke = 0.25) +
  geom_text(aes(label = sig), size = 3, fontface = "bold", vjust = 0.75, nudge_x = 0.33) +
  facet_grid(. ~ family, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-lim, lim), oob = scales::squish,
                       name = "CLR log2 fold change\n(case vs matched control)") +
  scale_size_area(max_size = 7, breaks = c(2, 10, 25, 50), name = "Control abundance (%)") +
  scale_y_discrete(labels = class_lab) +
  labs(title = "Cell-class composition changes by contrast",
       subtitle = paste0("Colour: CLR log2 fold change (composition-aware; capped at ±3.5). Stars: limma robust eBayes on the same CLR values,\n",
                         "BH FDR within contrast (* < 0.05, ** < 0.01, *** < 0.001). VCD model adjusts for replicate (paired _1/_2 libraries).\n",
                         "Libraries per contrast: VCD 2 vs 2; Foxl2 2 vs 2 (4 mo), 3 vs 3 (9 mo); rice 3 vs 3 (6 mo), 3 vs 2 (9 mo)."),
       x = NULL, y = NULL) +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey92"),
        strip.background = element_rect(fill = "grey95", colour = "grey70", linewidth = 0.25),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 7.5), axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
        legend.key = element_blank())
ggsave(file.path(comp_dir, "stress_test_cell_class_composition_dotheatmap.pdf"), p, width = 10, height = 3.8, device = cairo_pdf)
ggsave(file.path(comp_dir, "stress_test_cell_class_composition_dotheatmap.png"), p, width = 10, height = 3.8, dpi = 300, bg = "white")
print(res[fdr < 0.05 | propeller_fdr < 0.05, .(contrast, harmonized_cell_class, clr_log2fc = round(clr_log2fc, 2),
                        pct_change = round(pct_point_change, 1), fdr_clr = signif(fdr, 2),
                        fdr_propeller = signif(propeller_fdr, 2))][order(contrast)])
message("Wrote ", comp_dir)
