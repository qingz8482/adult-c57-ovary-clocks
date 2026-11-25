# Correlate back-projected gene weights of the previous Version B PCA/glmnet clock
# (2,000 genes, 7 PCs, ridge) with the locked rebutal4 Version B PCA/glmnet clock
# (3,000 genes, 10 PCs, lasso). Back-projected weight = PCA loadings x PC coefficients,
# on the standardized-gene scale.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel)
})
root <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks")
out  <- file.path(root, "rebutal/rebutal4/backprojection_PCA")

## Previous clock (deployment_all_control_like)
old_rds <- file.path(root, "bulk-pseudobulk_clock/results/model_artifacts/deployment/pca_glmnet__version_b_age_dataset_penalized/deployment_all_control_like_pca_glmnet_model.rds")
m <- readRDS(old_rds)
pcs  <- m$pc_coefficient_table$pc
rot  <- m$pca$rotation[, pcs, drop = FALSE]
w_old <- as.numeric(rot %*% m$pc_coefficient_table$calibrated_pc_coefficient)
old <- data.table(gene = rownames(rot), w_old = w_old)
stopifnot(nrow(old) == 2000L, identical(old$gene, m$genes))

# Cross-check with the archived back-projection of the same model
arch <- fread(file.path(root, "bulk-pseudobulk_clock/results/external_projection_archive/external_projection_before_added7_20260727_112828/tAge/feature_overlap_analysis/audit/ovary_pca_deployment_rds_backprojection.tsv"))
chk <- merge(old, arch[, .(gene = gene_id, gene_score)], by = "gene")
stopifnot(nrow(chk) == 2000L, max(abs(chk$w_old - chk$gene_score)) < 1e-10)

## Rebutal4 locked clock
L <- file.path(root, "rebutal/rebutal4/lockedmodel/VersionB-PCA-glmnet")
new <- fread(file.path(L, "coefficients.csv"))
load <- as.matrix(fread(file.path(L, "pca_loadings.csv")), rownames = 1)
pcc  <- fread(file.path(L, "pc_coefficients.csv"))
w_chk <- as.numeric(load %*% pcc$coefficient)
stopifnot(identical(rownames(load), new$gene), max(abs(w_chk - new$standardized_coefficient)) < 1e-12)
new <- new[, .(gene, w_new = standardized_coefficient)]
stopifnot(nrow(new) == 3000L)

## Merge and statistics
d <- merge(old, new, by = "gene")
n_shared <- nrow(d)
pr <- cor.test(d$w_old, d$w_new, method = "pearson")
sp <- suppressWarnings(cor.test(d$w_old, d$w_new, method = "spearman"))
sign_agree <- mean(sign(d$w_old) == sign(d$w_new))

top_overlap <- function(k) {
  a <- d[order(-abs(w_old))][1:k]
  b <- d[order(-abs(w_new))][1:k]
  shared <- intersect(a$gene, b$gene)
  data.table(top_k = k, n_overlap = length(shared),
             sign_concordant = sum(sign(d[gene %in% shared, w_old]) == sign(d[gene %in% shared, w_new])),
             expected_by_chance = k * k / n_shared)
}
topk <- rbindlist(lapply(c(50, 100, 200), top_overlap))

summ <- data.table(
  n_previous_genes = nrow(old), n_rebutal4_genes = nrow(new), n_shared = n_shared,
  n_previous_only = nrow(old) - n_shared, n_rebutal4_only = nrow(new) - n_shared,
  previous_pcs = paste(pcs, collapse = ","), rebutal4_active_pcs = paste0("PC", pcc[coefficient != 0, PC], collapse = ","),
  pearson_r = unname(pr$estimate), pearson_p = pr$p.value,
  spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
  sign_agreement = sign_agree)
fwrite(summ, file.path(out, "backprojection_weight_correlation_summary.tsv"), sep = "\t")
fwrite(topk, file.path(out, "backprojection_top_k_overlap.tsv"), sep = "\t")

d[, rank_old := frank(-abs(w_old))][, rank_new := frank(-abs(w_new))]
fwrite(d[order(rank_new)], file.path(out, "backprojection_weights_previous_vs_rebutal4.tsv"), sep = "\t")

## Plot
d[, label := ifelse(rank_old <= 6 | rank_new <= 6, gene, NA_character_)]
d[, concord := ifelse(sign(w_old) == sign(w_new), "Same sign", "Opposite sign")]
fmt_p <- function(p) if (p < 1e-300) "< 1e-300" else paste0("= ", formatC(p, format = "e", digits = 1))
lab <- sprintf("n = %d shared genes\nPearson r = %.2f, P %s\nSpearman ρ = %.2f, P %s\nSign agreement = %.1f%%",
               n_shared, pr$estimate, fmt_p(pr$p.value), sp$estimate, fmt_p(sp$p.value), 100 * sign_agree)

fam <- "Arial"
p <- ggplot(d, aes(w_old, w_new)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_point(aes(colour = concord), size = 1, alpha = 0.6, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.5) +
  geom_text_repel(aes(label = label), family = fam, size = 12 / .pt, na.rm = TRUE,
                  min.segment.length = 0, segment.size = 0.25, max.overlaps = Inf, seed = 1, box.padding = 0.5) +
  annotate("text", x = min(d$w_old), y = max(d$w_new) * 1.12, label = lab, hjust = 0, vjust = 1,
           family = fam, size = 12 / .pt, lineheight = 0.95) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  scale_colour_manual(values = c("Same sign" = "#2166AC", "Opposite sign" = "#B2182B"), name = NULL) +
  labs(x = "Back-projected gene weight, previous clock\n(2,000 genes, PC1–PC7)",
       y = "Back-projected gene weight, revised clock\n(3,000 genes, PC1, PC2, PC5, PC8)") +
  theme_classic(base_size = 12, base_family = fam) +
  theme(text = element_text(size = 12, family = fam, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12), legend.position = "bottom")

ggsave(file.path(out, "backprojection_weight_correlation.pdf"), p, width = 6.5, height = 6.8, device = cairo_pdf)
ggsave(file.path(out, "backprojection_weight_correlation.png"), p, width = 6.5, height = 6.8, dpi = 300, device = ragg::agg_png)
ggsave(file.path(out, "backprojection_weight_correlation.tiff"), p, width = 6.5, height = 6.8, dpi = 300, device = ragg::agg_tiff, compression = "lzw")
print(summ); print(topk)
# Opposite-sign genes: how large are their weights?
d[, pct_old := frank(abs(w_old)) / .N][, pct_new := frank(abs(w_new)) / .N]
opp <- d[concord == "Opposite sign"]
cat("opposite-sign genes:", nrow(opp), "; in top 10% |weight| of either clock:", opp[pct_old > 0.9 | pct_new > 0.9, .N],
    "; median |w| percentile old/new:", round(median(opp$pct_old), 2), round(median(opp$pct_new), 2), "\n")
