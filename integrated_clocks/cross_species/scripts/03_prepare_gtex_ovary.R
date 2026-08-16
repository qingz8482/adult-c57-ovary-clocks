#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))
suppressPackageStartupMessages({ library(edgeR); library(limma) })

raw_dir <- path_module("derived_data/gtex/raw")
count_path <- file.path(raw_dir, basename(config$gtex$counts_url))
sample_path <- file.path(raw_dir, basename(config$gtex$sample_attributes_url))
subject_path <- file.path(raw_dir, basename(config$gtex$subject_phenotypes_url))

gct <- fread(cmd = paste("gzip -cd", shQuote(count_path)), skip = 2)
if (!all(c("Name", "Description") %in% names(gct))) stop("Unexpected GTEx GCT schema.")
sample_ids <- setdiff(names(gct), c("Name", "Description"))
counts <- as.matrix(gct[, ..sample_ids]); mode(counts) <- "numeric"
gene_id_original <- as.character(gct$Name)
gene_symbol <- as.character(gct$Description)
gene_id <- strip_version(gene_id_original)
if (anyDuplicated(gene_id)) {
  groups <- split(seq_along(gene_id), gene_id)
  counts <- do.call(rbind, lapply(groups, function(i) colSums(counts[i, , drop = FALSE])))
  gene_symbol <- vapply(groups, function(i) gene_symbol[i[1]], character(1))
  gene_id_original <- vapply(groups, function(i) gene_id_original[i[1]], character(1))
  gene_id <- names(groups)
}
rownames(counts) <- gene_id

sa <- fread(sample_path)
sp <- fread(subject_path)
sa <- sa[SAMPID %in% sample_ids]
if (nrow(sa) != length(sample_ids)) stop("GTEx sample metadata does not match ovary count columns.")
sa[, donor_id := sub("^(GTEX-[^-]+).*", "\\1", SAMPID)]
meta <- merge(sa, sp, by.x = "donor_id", by.y = "SUBJID", all.x = TRUE, sort = FALSE)
meta <- meta[match(sample_ids, SAMPID)]
if (anyNA(meta$AGE)) stop("Public age brackets missing for GTEx ovary samples.")
meta[, age_midpoint_sensitivity := age_midpoint(AGE)]
meta[, age_decade := age_midpoint_sensitivity / 10]
meta[, rin := safe_num(SMRIN)]
meta[, ischemic_time := safe_num(SMTSISCH)]
meta[, library_size := colSums(counts)[SAMPID]]

dup <- meta[, .N, by = donor_id][N > 1]
if (nrow(dup)) {
  meta[, quality_rank := frank(-rin, -library_size, ties.method = "first"), by = donor_id]
  keep_ids <- meta[quality_rank == 1, SAMPID]
} else {
  meta[, quality_rank := 1L]
  keep_ids <- meta$SAMPID
}
counts <- counts[, keep_ids, drop = FALSE]
meta <- meta[match(keep_ids, SAMPID)]
if (anyDuplicated(meta$donor_id)) stop("Donor deduplication failed.")

min_n <- ceiling(config$gtex$filter_min_fraction * ncol(counts))
keep_gene <- rowSums(counts >= config$gtex$filter_min_count) >= min_n
counts_f <- counts[keep_gene, , drop = FALSE]
dge <- DGEList(counts_f)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = config$gtex$prior_count)

cov_audit <- rbindlist(list(
  data.table(covariate = "RIN", completeness = mean(is.finite(meta$rin)), variation = sd(meta$rin, na.rm = TRUE),
             included = mean(is.finite(meta$rin)) >= .8 && sd(meta$rin, na.rm = TRUE) > 0,
             reason = "continuous technical quality covariate"),
  data.table(covariate = "ischemic_time", completeness = mean(is.finite(meta$ischemic_time)), variation = sd(meta$ischemic_time, na.rm = TRUE),
             included = mean(is.finite(meta$ischemic_time)) >= .8 && sd(meta$ischemic_time, na.rm = TRUE) > 0,
             reason = "included only with >=80% completeness"),
  data.table(covariate = "extraction_batch", completeness = mean(!is.na(meta$SMNABTCH) & meta$SMNABTCH != ""), variation = uniqueN(meta$SMNABTCH),
             included = FALSE, reason = "too many sparse levels for n=193 ovary donors"),
  data.table(covariate = "sequencing_batch", completeness = mean(!is.na(meta$SMGEBTCH) & meta$SMGEBTCH != ""), variation = uniqueN(meta$SMGEBTCH),
             included = FALSE, reason = "too many sparse levels for n=193 ovary donors")
), fill = TRUE)
write_tsv(cov_audit, path_module("results/qc", "gtex_technical_covariate_audit.tsv"))

analysis_meta <- data.table(sample_id = meta$SAMPID, donor_id = meta$donor_id, tissue = meta$SMTSD,
                            age_bracket = meta$AGE, age_midpoint_sensitivity = meta$age_midpoint_sensitivity,
                            sex = meta$SEX, rin = meta$rin, ischemic_time = meta$ischemic_time,
                            sequencing_batch = meta$SMGEBTCH, extraction_batch = meta$SMNABTCH,
                            library_size = meta$library_size, included_primary = TRUE,
                            exclusion_reason = NA_character_)
analysis_meta[, age_bracket := factor(age_bracket, levels = sort(unique(age_bracket), method = "radix"), ordered = TRUE)]
rownames(analysis_meta) <- analysis_meta$sample_id

design_data <- as.data.frame(analysis_meta)
design_data$age_decade <- design_data$age_midpoint_sensitivity / 10
formula_terms <- "age_decade"
if (cov_audit[covariate == "RIN", included]) { design_data$rin_z <- as.numeric(scale(design_data$rin)); formula_terms <- paste(formula_terms, "+ rin_z") }
if (cov_audit[covariate == "ischemic_time", included]) { design_data$ischemia_z <- as.numeric(scale(design_data$ischemic_time)); formula_terms <- paste(formula_terms, "+ ischemia_z") }
design <- model.matrix(as.formula(paste("~", formula_terms)), design_data)
vo <- voom(dge, design, plot = FALSE)
fit <- eBayes(lmFit(vo, design))
k <- which(colnames(design) == "age_decade")
rho <- apply(logcpm, 1, function(z) suppressWarnings(cor(z, analysis_meta$age_midpoint_sensitivity, method = "spearman")))
age_effects <- data.table(human_gene_id_stable = rownames(logcpm), human_gene_symbol = gene_symbol[match(rownames(logcpm), gene_id)],
                          human_age_beta_per_decade = fit$coefficients[, k], human_age_t = fit$t[, k],
                          human_age_p = fit$p.value[, k], human_age_fdr = p.adjust(fit$p.value[, k], "BH"),
                          human_age_spearman = rho, human_age_monotonicity_flag = sign(rho) == sign(fit$coefficients[, k]),
                          n_donors = ncol(logcpm), n_age_brackets = uniqueN(analysis_meta$age_bracket),
                          technical_model = deparse(as.formula(paste("~", formula_terms))),
                          human_mean_expression = rowMeans(logcpm), human_expression_variance = apply(logcpm, 1, var),
                          human_detection_fraction = rowMeans(counts_f >= config$gtex$filter_min_count))

saveRDS(list(counts_filtered = counts_f, logcpm = logcpm, dge = dge, metadata = analysis_meta,
             age_effects = age_effects, design = design, normalization_factors = dge$samples$norm.factors,
             gene_id_original = gene_id_original, gene_symbol = gene_symbol, filter_min_n = min_n),
        path_module("derived_data/gtex", "gtex_ovary_prepared.rds"))
write_tsv(analysis_meta, path_module("derived_data/gtex", "gtex_ovary_sample_audit.tsv"))
write_tsv(age_effects, path_module("derived_data/gtex", "gtex_gene_age_effects.tsv"))
write_tsv(data.table(sample_id = colnames(logcpm), tmm_factor = dge$samples$norm.factors,
                     library_size = dge$samples$lib.size), path_module("derived_data/gtex", "gtex_tmm_factors.tsv"))
write_tsv(analysis_meta[, .(n = .N), by = age_bracket], path_module("results/tables", "table_gtex_sample_characteristics.tsv"))
message("GTEx ovary prepared: ", nrow(analysis_meta), " donors; ", nrow(logcpm), " filtered genes.")
