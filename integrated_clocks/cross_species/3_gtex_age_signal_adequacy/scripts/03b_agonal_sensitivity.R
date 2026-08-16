#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]))), "common.R"))

message("[3b] Donor-death/agonal sensitivity")
x <- load_core_inputs()
gtex <- x$gtex; meta <- x$meta
qc_obj <- readRDS(paudit("derived_data/gtex_age_signal_qc_objects.rds"))
audit <- qc_obj$audit[match(meta$sample_id, sample_id)]
keep <- is.finite(audit$donor_death_hardy_scale)
dd <- data.frame(
  age_order = as.integer(audit$age_bracket[keep]) - 1L,
  rin_z = zsafe(audit$rin[keep]), ischemia_z = zsafe(audit$ischemic_time_minutes[keep]),
  mapping_z = zsafe(audit$mapping_rate[keep]), exonic_z = zsafe(audit$exonic_read_rate[keep]),
  duplication_z = zsafe(audit$duplicate_mapped_read_rate[keep]),
  sequencing_center = droplevels(factor(audit$sequencing_center[keep])),
  death_hardy = droplevels(factor(audit$donor_death_hardy_scale[keep]))
)
design <- model.matrix(~ age_order + rin_z + ischemia_z + mapping_z + exonic_z + duplication_z + sequencing_center + death_hardy, dd)
q <- qr(design); design <- design[, sort(q$pivot[seq_len(q$rank)]), drop = FALSE]
dge <- DGEList(gtex$counts_filtered[, keep, drop = FALSE])
dge$samples$norm.factors <- gtex$normalization_factors[keep]
fit <- glmQLFit(dge, design, robust = TRUE)
test <- glmQLFTest(fit, coef = match("age_order", colnames(design)))
tab <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "human_gene_id_stable")
tab[, signed_quasi_likelihood_t := sign(logFC) * sqrt(pmax(F, 0))]
gm <- gtex$age_effects[, .(human_gene_id_stable, human_gene_symbol, human_detection_fraction)]
tab <- merge(gm, tab, by = "human_gene_id_stable", sort = FALSE)
tab[, `:=`(analysis = "technical_plus_donor_death_hardy_factor", n_donors = sum(keep),
  model = "edgeR QL raw counts ~ ordered age + prespecified technical covariates + donor-death Hardy factor")]
write_tsv(tab, paudit("results/table_ordered_age_trend_donor_death_sensitivity_full.tsv"))
summary <- tab[, .(
  n_donors = unique(n_donors), total_genes = .N, fdr_lt_0_05 = sum(FDR < .05), fdr_lt_0_10 = sum(FDR < .10),
  fdr_lt_0_05_abs_log2fc_gt_0_10 = sum(FDR < .05 & abs(logFC) > .10),
  fdr_lt_0_05_abs_log2fc_gt_0_25 = sum(FDR < .05 & abs(logFC) > .25),
  fdr_lt_0_05_abs_log2fc_gt_0_50 = sum(FDR < .05 & abs(logFC) > .50)
)]
write_tsv(summary, paudit("results/table_ordered_age_trend_donor_death_sensitivity_summary.tsv"))
message("Agonal sensitivity complete: ", summary$fdr_lt_0_05, " genes at FDR<0.05 after donor-death adjustment.")
