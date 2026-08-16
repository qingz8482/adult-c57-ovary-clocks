#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common.R"))

message("[1/4] Loading GTEx metadata, frozen scores, and expression")
x <- load_core_inputs()
gtex <- x$gtex; sa <- x$sa; sp <- x$sp; meta <- x$meta; proj <- x$proj
counts <- gtex$counts_filtered
logcpm <- gtex$logcpm
stopifnot(ncol(counts) == 193L, identical(colnames(counts), meta$sample_id))

## Donor-level expression QC ----------------------------------------------------
message("Computing sample-level expression QC and robust PCA")
lib_size <- colSums(counts)
detected_ge1 <- colSums(counts > 0)
detected_ge10 <- colSums(counts >= 10)
global_median <- rowMedians(logcpm)
pearson_global <- apply(logcpm, 2, cor, y = global_median, method = "pearson")
spearman_global <- apply(logcpm, 2, cor, y = global_median, method = "spearman")

pca <- prcomp(t(logcpm), center = TRUE, scale. = FALSE, rank. = 30)
pcs <- as.data.table(pca$x)
pcs[, sample_id := rownames(pca$x)]
setcolorder(pcs, "sample_id")
pve <- pca$sdev^2 / sum(pca$sdev^2)

robust_dims <- seq_len(min(10L, ncol(pca$x)))
robfit <- cov.rob(pca$x[, robust_dims, drop = FALSE], method = "mcd", quantile.used = floor(0.75 * nrow(pca$x)))
robust_md2 <- mahalanobis(pca$x[, robust_dims, drop = FALSE], center = robfit$center, cov = robfit$cov)
robust_md <- sqrt(robust_md2)
robust_cut <- sqrt(qchisq(0.999, df = length(robust_dims)))

qc <- data.table(
  sample_id = meta$sample_id, donor_id = meta$donor_id,
  raw_library_size = lib_size[meta$sample_id],
  expressed_genes_count_gt0 = detected_ge1[meta$sample_id],
  expressed_genes_count_ge10 = detected_ge10[meta$sample_id],
  sample_global_median_pearson = pearson_global[meta$sample_id],
  sample_global_median_spearman = spearman_global[meta$sample_id],
  robust_pca_mahalanobis = robust_md[meta$sample_id]
)
qc <- merge(qc, pcs, by = "sample_id", all.x = TRUE, sort = FALSE)

## Frozen scores/OOD in wide form ---------------------------------------------
proj_keep <- proj[, .(
  donor_id, sample_id, model_family,
  existing_score = cross_species_mouse_ovarian_clock_score,
  mouse_week_equivalent_output_audit_only,
  projection_ood_distance,
  nearest_mouse_reference_distance,
  fraction_features_abs_z_gt3,
  fraction_features_abs_z_gt5,
  fraction_mapped_features_outside_mouse_range,
  projection_qc_class
)]
proj_wide <- dcast(proj_keep, donor_id + sample_id ~ model_family,
  value.var = setdiff(names(proj_keep), c("donor_id", "sample_id", "model_family")), sep = "_")

## Retain every official public field with stable lower-case names -------------
sa_all <- copy(sa); sp_all <- copy(sp)
setnames(sa_all, names(sa_all), paste0("gtex_sample_", safe_name(names(sa_all))))
setnames(sp_all, names(sp_all), paste0("gtex_donor_", safe_name(names(sp_all))))
sa_all[, sample_id := get("gtex_sample_sampid")]
sp_all[, donor_id := get("gtex_donor_subjid")]

audit <- data.table(
  donor_id = meta$donor_id, sample_id = meta$sample_id,
  tissue = meta$tissue, age_bracket = factor(as.character(meta$age_bracket), levels = age_levels, ordered = TRUE),
  age_bracket_midpoint_sensitivity_only = meta$age_midpoint_sensitivity,
  raw_library_size = lib_size[meta$sample_id],
  expressed_genes_count_gt0 = detected_ge1[meta$sample_id],
  expressed_genes_count_ge10 = detected_ge10[meta$sample_id],
  rin = safe_num(sa$SMRIN), ischemic_time_minutes = safe_num(sa$SMTSISCH),
  mapping_rate = safe_num(sa$SMMAPRT), exonic_read_rate = safe_num(sa$SMEXNCRT),
  duplicate_mapped_read_rate = safe_num(sa$SMDPMPRT),
  intergenic_read_rate = safe_num(sa$SMNTRART), rrna_read_rate = safe_num(sa$SMRRNART),
  expression_profiling_efficiency = safe_num(sa$SMEXPEFF),
  sequencing_center = as.character(sa$SMCENTER), sequencing_batch = as.character(sa$SMGEBTCH),
  library_extraction_batch = as.character(sa$SMNABTCH),
  donor_death_hardy_scale = safe_num(sp$DTHHRDY),
  pathology_notes = as.character(sa$SMPTHNTS)
)
audit <- merge(audit, qc, by = c("donor_id", "sample_id"), all.x = TRUE, sort = FALSE,
  suffixes = c("", "_qc"))
audit <- merge(audit, proj_wide, by = c("donor_id", "sample_id"), all.x = TRUE, sort = FALSE)
audit <- merge(audit, sa_all, by = "sample_id", all.x = TRUE, sort = FALSE)
audit <- merge(audit, sp_all, by = "donor_id", all.x = TRUE, sort = FALSE)
audit <- audit[match(meta$sample_id, sample_id)]

## Direct pathology-text audit; no expression-derived histology ----------------
ptxt <- tolower(audit$pathology_notes)
flag_kw <- function(pattern) grepl(pattern, ptxt, perl = TRUE)
audit[, `:=`(
  histology_annotation_available = !is.na(pathology_notes) & nzchar(pathology_notes),
  histology_atrophy_or_postmenopausal = flag_kw("atroph|post[ -]?menopaus|involution"),
  histology_follicular_material = flag_kw("follic|corpus lute|corpora lute|corpora albican|corpus albican"),
  histology_stromal_predominance = flag_kw("stroma|fibrovascular|fibrous"),
  histology_adipose_predominance = flag_kw("adipos|fat "),
  histology_cystic_or_pathological = flag_kw("cyst|tumou?r|neoplas|carcin|metasta|endometri|hemorrhag|necros|inflamm"),
  histology_cortex_mentioned = flag_kw("cortex|cortical"),
  histology_medulla_mentioned = flag_kw("medulla|medullary"),
  histology_possible_adjacent_tissue = flag_kw("fallopian|tube|uter|omentum|adjacent|skeletal muscle|lymph")
)]

## Predefined technical-QC exclusions (never based on age or clock score) -------
low_cor_cut <- median(audit$sample_global_median_pearson, na.rm = TRUE) - 5 * mad(audit$sample_global_median_pearson, constant = 1.4826, na.rm = TRUE)
lib_lo <- quantile(audit$raw_library_size, 0.25) - 3 * IQR(audit$raw_library_size)
det_lo <- quantile(audit$expressed_genes_count_ge10, 0.25) - 3 * IQR(audit$expressed_genes_count_ge10)
map_lo <- quantile(audit$mapping_rate, 0.25) - 3 * IQR(audit$mapping_rate)
exon_lo <- quantile(audit$exonic_read_rate, 0.25) - 3 * IQR(audit$exonic_read_rate)
audit[, `:=`(
  qc_flag_rin_lt_5_5 = is.na(rin) | rin < 5.5,
  qc_flag_library_extreme_low = raw_library_size < lib_lo,
  qc_flag_detected_gene_extreme_low = expressed_genes_count_ge10 < det_lo,
  qc_flag_global_correlation_low = sample_global_median_pearson < low_cor_cut,
  qc_flag_robust_pca_outlier = robust_pca_mahalanobis > robust_cut,
  qc_flag_mapping_extreme_low = mapping_rate < map_lo,
  qc_flag_exonic_extreme_low = exonic_read_rate < exon_lo
)]
flag_cols <- grep("^qc_flag_", names(audit), value = TRUE)
audit[, qc_exclusion_any := rowSums(.SD, na.rm = TRUE) > 0, .SDcols = flag_cols]
audit[, qc_exclusion_reasons := apply(.SD, 1, function(z) paste(sub("^qc_flag_", "", flag_cols[which(z)]), collapse = ";")), .SDcols = flag_cols]
audit[, primary_all_eligible_samples := TRUE]
audit[, sensitivity_technical_qc_pass := !qc_exclusion_any]

## Missingness indicators for all requested/curated fields ----------------------
missing_basis <- c(
  "rin", "ischemic_time_minutes", "mapping_rate", "exonic_read_rate", "duplicate_mapped_read_rate",
  "intergenic_read_rate", "rrna_read_rate", "expression_profiling_efficiency", "sequencing_center",
  "sequencing_batch", "library_extraction_batch", "donor_death_hardy_scale", "pathology_notes",
  "existing_score_elasticnet", "existing_score_pca_glmnet", "projection_ood_distance_elasticnet",
  "projection_ood_distance_pca_glmnet"
)
for (v in intersect(missing_basis, names(audit))) {
  miss <- is.na(audit[[v]]) | (is.character(audit[[v]]) & audit[[v]] == "")
  audit[, (paste0("missing_", v)) := miss]
}
audit[, n_missing_requested_fields := rowSums(.SD), .SDcols = patterns("^missing_")]

write_tsv(audit, paudit("results/table_gtex_ovary_donor_level_audit.tsv"))
write_tsv(audit[, .(sample_id, donor_id, age_bracket, pathology_notes,
  histology_annotation_available, histology_atrophy_or_postmenopausal, histology_follicular_material,
  histology_stromal_predominance, histology_adipose_predominance, histology_cystic_or_pathological,
  histology_cortex_mentioned, histology_medulla_mentioned, histology_possible_adjacent_tissue)],
  paudit("results/table_gtex_ovary_histology_pathology_audit.tsv"))

## Exact age-bin counts ---------------------------------------------------------
age_counts <- audit[, .(
  n_samples = .N, n_donors = uniqueN(donor_id),
  n_primary = sum(primary_all_eligible_samples), n_qc_sensitivity = sum(sensitivity_technical_qc_pass)
), by = age_bracket]
age_counts[, proportion := n_samples / sum(n_samples)]
write_tsv(age_counts, paudit("results/table_gtex_ovary_age_bin_counts.tsv"))

## QC summaries by age bin -----------------------------------------------------
qc_vars <- c(
  "raw_library_size", "expressed_genes_count_gt0", "expressed_genes_count_ge10", "rin",
  "ischemic_time_minutes", "mapping_rate", "exonic_read_rate", "duplicate_mapped_read_rate",
  "intergenic_read_rate", "rrna_read_rate", "expression_profiling_efficiency",
  "sample_global_median_pearson", "sample_global_median_spearman", "robust_pca_mahalanobis",
  "projection_ood_distance_elasticnet", "projection_ood_distance_pca_glmnet"
)
qc_by_bin <- rbindlist(lapply(intersect(qc_vars, names(audit)), function(v) {
  audit[, {
    z <- safe_num(get(v)); ok <- is.finite(z)
    .(variable = v, n_nonmissing = sum(ok), n_missing = sum(!ok), mean = mean(z, na.rm = TRUE),
      sd = sd(z, na.rm = TRUE), median = median(z, na.rm = TRUE), q25 = quantile(z, .25, na.rm = TRUE),
      q75 = quantile(z, .75, na.rm = TRUE), min = min(z, na.rm = TRUE), max = max(z, na.rm = TRUE))
  }, by = age_bracket]
}), fill = TRUE)
write_tsv(qc_by_bin, paudit("results/table_gtex_ovary_qc_by_age_bin.tsv"))

## Technical-variable balance and missingness across age bins ------------------
message("Testing technical-variable balance across ordered age brackets")
age_factor <- audit$age_bracket; age_order <- as.integer(age_factor)
raw_sample_names <- grep("^gtex_sample_", names(audit), value = TRUE)
exclude_balance <- c("gtex_sample_sampid", "gtex_sample_smts", "gtex_sample_smtsd", "gtex_sample_smpthnts")
technical_vars <- unique(c(qc_vars, setdiff(raw_sample_names, exclude_balance),
  "sequencing_center", "sequencing_batch", "library_extraction_batch", "donor_death_hardy_scale"))
technical_vars <- intersect(technical_vars, names(audit))

balance_rows <- lapply(technical_vars, function(v) {
  raw <- audit[[v]]
  numeric_like <- is.numeric(raw) || (mean(is.finite(safe_num(raw))) >= 0.8 && uniqueN(raw, na.rm = TRUE) > 6)
  missing <- is.na(raw) | (is.character(raw) & raw == "")
  miss_p <- if (length(unique(missing)) > 1L) suppressWarnings(chisq.test(table(age_factor, missing), simulate.p.value = TRUE, B = 2000)$p.value) else NA_real_
  if (numeric_like) {
    z <- safe_num(raw); a <- association_numeric(z, age_factor, age_order)
    data.table(variable = v, variable_type = "continuous", n_nonmissing = sum(is.finite(z)), n_missing = sum(!is.finite(z)),
      n_levels = uniqueN(z[is.finite(z)]), standardized_mean_difference_max_pairwise = max_pairwise_smd(z, age_factor),
      cramers_v = NA_real_, ordered_spearman_rho = a$rho, omnibus_p = a$p_omnibus,
      ordered_association_p = a$p_ordered, missingness_age_association_p = miss_p,
      test = "ANOVA across age bins; Spearman ordered trend")
  } else {
    z <- as.character(raw); ok <- !missing
    tab <- table(age_factor[ok], z[ok]); nl <- ncol(tab)
    tst <- if (nrow(tab) >= 2L && nl >= 2L) suppressWarnings(chisq.test(tab, simulate.p.value = nl <= 20L, B = 5000)) else NULL
    data.table(variable = v, variable_type = "categorical", n_nonmissing = sum(ok), n_missing = sum(!ok),
      n_levels = uniqueN(z[ok]), standardized_mean_difference_max_pairwise = NA_real_,
      cramers_v = if (!is.null(tst)) cramers_v(tab) else NA_real_, ordered_spearman_rho = NA_real_,
      omnibus_p = if (!is.null(tst)) tst$p.value else NA_real_, ordered_association_p = NA_real_,
      missingness_age_association_p = miss_p,
      test = ifelse(nl <= 20L, "chi-square with Monte Carlo p", "asymptotic chi-square; sparse high-cardinality factor"))
  }
})
balance <- rbindlist(balance_rows, fill = TRUE)
balance[, `:=`(omnibus_fdr = p.adjust(omnibus_p, "BH"), ordered_fdr = p.adjust(ordered_association_p, "BH"),
  missingness_fdr = p.adjust(missingness_age_association_p, "BH"))]
write_tsv(balance, paudit("results/table_gtex_ovary_metadata_balance.tsv"))

## Expression-PC associations --------------------------------------------------
pc_long <- merge(pcs, audit[, .(sample_id, age_bracket, rin, ischemic_time_minutes, sequencing_center,
  sequencing_batch, library_extraction_batch, mapping_rate, exonic_read_rate, duplicate_mapped_read_rate,
  raw_library_size, expressed_genes_count_ge10, donor_death_hardy_scale)], by = "sample_id")
pc_covars <- setdiff(names(pc_long), c("sample_id", paste0("PC", seq_len(ncol(pca$x)))))
pc_assoc <- list()
for (pcn in paste0("PC", 1:10)) for (v in pc_covars) {
  y <- pc_long[[pcn]]; raw <- pc_long[[v]]
  if (is.numeric(raw)) {
    ok <- is.finite(y) & is.finite(raw)
    tst <- if (sum(ok) >= 10L && sd(raw[ok]) > 0) suppressWarnings(cor.test(y[ok], raw[ok], method = "spearman", exact = FALSE)) else NULL
    pc_assoc[[length(pc_assoc) + 1L]] <- data.table(pc = pcn, covariate = v, covariate_type = "continuous",
      n = sum(ok), n_levels = uniqueN(raw[ok]), effect = if (!is.null(tst)) unname(tst$estimate) else NA_real_,
      p = if (!is.null(tst)) tst$p.value else NA_real_, test = "Spearman correlation", status = "tested")
  } else {
    z <- factor(raw); ok <- is.finite(y) & !is.na(z)
    nl <- nlevels(droplevels(z[ok]))
    if (nl >= 2L && nl <= 20L) {
      av <- anova(lm(y[ok] ~ droplevels(z[ok])))
      eta <- av[1, "Sum Sq"] / sum(av[, "Sum Sq"])
      pc_assoc[[length(pc_assoc) + 1L]] <- data.table(pc = pcn, covariate = v, covariate_type = "categorical",
        n = sum(ok), n_levels = nl, effect = eta, p = av[1, "Pr(>F)"], test = "ANOVA partial eta-squared", status = "tested")
    } else {
      pc_assoc[[length(pc_assoc) + 1L]] <- data.table(pc = pcn, covariate = v, covariate_type = "categorical",
        n = sum(ok), n_levels = nl, effect = NA_real_, p = NA_real_, test = "not estimable",
        status = "sparse high-cardinality factor; audited descriptively")
    }
  }
}
pc_assoc <- rbindlist(pc_assoc, fill = TRUE)
pc_assoc[, fdr := p.adjust(p, "BH")]
write_tsv(pc_assoc, paudit("results/table_expression_pc_technical_associations.tsv"))
write_tsv(data.table(pc = paste0("PC", seq_along(pve)), proportion_variance = pve,
  cumulative_variance = cumsum(pve)), paudit("results/table_expression_pca_variance.tsv"))

## Hierarchical clustering -----------------------------------------------------
hc <- hclust(dist(pca$x[, 1:20, drop = FALSE]), method = "ward.D2")
clust <- data.table(sample_id = hc$labels[hc$order], dendrogram_order = seq_along(hc$order))
clust <- merge(clust, audit[, .(sample_id, donor_id, age_bracket, sensitivity_technical_qc_pass)], by = "sample_id", sort = FALSE)
write_tsv(clust, paudit("results/table_hierarchical_sample_clustering.tsv"))
save_base_bundle(paudit("figures/sample_hierarchical_clustering"), function() {
  par(mar = c(2.5, 3.5, 1.5, 0.5), cex = 0.45)
  plot(hc, labels = FALSE, hang = -1, main = "GTEx ovary sample clustering", xlab = "193 donors", sub = "Ward.D2 on expression PCs 1-20", ylab = "Height")
}, width_mm = 183, height_mm = 80)

## PCA plots requested as separate colour encodings ----------------------------
pdat <- merge(pcs[, 1:3], audit[, .(sample_id, age_bracket, rin, ischemic_time_minutes, sequencing_center,
  sequencing_batch, expressed_genes_count_ge10, projection_ood_distance_elasticnet,
  projection_ood_distance_pca_glmnet, sensitivity_technical_qc_pass)], by = "sample_id")
axis_labs <- c(x = sprintf("PC1 (%.1f%%)", 100 * pve[1]), y = sprintf("PC2 (%.1f%%)", 100 * pve[2]))
plot_cont <- function(v, title, low = "#D8E9F3", high = "#B2182B") {
  ggplot(pdat, aes(PC1, PC2, colour = .data[[v]])) + geom_point(size = 1.4, alpha = .9) +
    scale_colour_gradient(low = low, high = high, na.value = "#BDBDBD") +
    labs(x = axis_labs["x"], y = axis_labs["y"], colour = title, title = title)
}
p_age <- ggplot(pdat, aes(PC1, PC2, colour = age_bracket)) + geom_point(size = 1.4, alpha = .9) +
  scale_colour_manual(values = palette_age, drop = FALSE) + labs(x = axis_labs["x"], y = axis_labs["y"], colour = "Age", title = "Age bracket")
p_rin <- plot_cont("rin", "RIN")
p_isch <- plot_cont("ischemic_time_minutes", "Ischemic time")
p_det <- plot_cont("expressed_genes_count_ge10", "Detected genes")
p_ood_el <- plot_cont("projection_ood_distance_elasticnet", "ElasticNet OOD")
p_ood_pc <- plot_cont("projection_ood_distance_pca_glmnet", "PCA/glmnet OOD")
p_center <- ggplot(pdat, aes(PC1, PC2, colour = sequencing_center)) + geom_point(size = 1.4, alpha = .9) +
  scale_colour_brewer(palette = "Dark2") + labs(x = axis_labs["x"], y = axis_labs["y"], colour = "Center", title = "Sequencing center")

save_pub_r(p_age, paudit("figures/pca_by_age_bracket"), 90, 75)
save_pub_r(p_rin, paudit("figures/pca_by_rin"), 90, 75)
save_pub_r(p_isch, paudit("figures/pca_by_ischemic_time"), 90, 75)
save_pub_r(p_center, paudit("figures/pca_by_sequencing_batch"), 90, 75)
save_pub_r(p_det, paudit("figures/pca_by_detected_gene_count"), 90, 75)
save_pub_r(p_ood_el, paudit("figures/pca_by_elasticnet_ood_distance"), 90, 75)
save_pub_r(p_ood_pc, paudit("figures/pca_by_pca_glmnet_ood_distance"), 90, 75)
p_grid <- (p_age | p_rin | p_isch) / (p_center | p_det | p_ood_el) + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8), legend.position = "right")
save_pub_r(p_grid, paudit("figures/expression_pca_qc_panel"), 183, 118)

## Library/detected-gene distributions and outlier summary ---------------------
qclong <- melt(audit[, .(sample_id, age_bracket, raw_library_size, expressed_genes_count_ge10)],
  id.vars = c("sample_id", "age_bracket"), variable.name = "metric", value.name = "value")
p_dist <- ggplot(qclong, aes(age_bracket, value, fill = age_bracket)) +
  geom_violin(scale = "width", linewidth = .25, alpha = .7) + geom_boxplot(width = .18, outlier.shape = NA, linewidth = .25) +
  scale_fill_manual(values = palette_age, drop = FALSE) + facet_wrap(~metric, scales = "free_y", nrow = 1) +
  labs(x = "Age bracket", y = NULL, title = "Library size and detected-gene distributions") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
save_pub_r(p_dist, paudit("figures/library_and_detected_gene_distributions"), 150, 65)

outlier_summary <- audit[, c(list(n = .N, n_any_qc_exclusion = sum(qc_exclusion_any)),
  lapply(.SD, sum, na.rm = TRUE)), by = age_bracket, .SDcols = flag_cols]
write_tsv(outlier_summary, paudit("results/table_predefined_technical_qc_outliers.tsv"))

saveRDS(list(audit = audit, pca = pca, pve = pve, hclust = hc, technical_variables = technical_vars,
  qc_thresholds = list(rin = 5.5, library_low = lib_lo, detected_low = det_lo, pearson_low = low_cor_cut,
    robust_md = robust_cut, mapping_low = map_lo, exonic_low = exon_lo)),
  paudit("derived_data/gtex_age_signal_qc_objects.rds"))
message("Metadata/expression QC complete: ", nrow(audit), " donors; ", sum(audit$qc_exclusion_any), " predefined technical-QC exclusions.")
