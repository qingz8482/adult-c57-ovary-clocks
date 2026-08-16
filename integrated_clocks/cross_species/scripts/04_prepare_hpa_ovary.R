#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))

zip_expr <- path_module("derived_data/hpa/raw", basename(config$hpa$transcript_url))
zip_meta <- path_module("derived_data/hpa/raw", basename(config$hpa$samples_url))
samples <- fread(cmd = paste("unzip -p", shQuote(zip_meta)))
setnames(samples, c("Tissue", "HPA sample id", "Age", "Sex"), c("tissue", "sample_id", "age", "sex"), skip_absent = TRUE)
ov <- samples[tolower(tissue) == "ovary"]
if (nrow(ov) != 7L) stop("Expected seven official HPA ovary samples; found ", nrow(ov))
ov[, sample_id := as.character(sample_id)]
ov[, age := safe_num(age)]

header <- fread(cmd = paste("unzip -p", shQuote(zip_expr)), nrows = 0)
tpm_cols <- paste0("TPM.ovary.", ov$sample_id)
if (!all(tpm_cols %in% names(header))) stop("HPA transcript matrix lacks expected ovary sample columns.")
tx <- fread(cmd = paste("unzip -p", shQuote(zip_expr)), select = c("ensgid", tpm_cols))
for (j in tpm_cols) set(tx, j = j, value = safe_num(tx[[j]]))
gene <- tx[, lapply(.SD, sum, na.rm = TRUE), by = ensgid, .SDcols = tpm_cols]
mat <- as.matrix(gene[, ..tpm_cols]); mode(mat) <- "numeric"; rownames(mat) <- gene$ensgid
colnames(mat) <- ov$sample_id
logexpr <- log2(mat + 1)

rho <- apply(logexpr, 1, function(z) suppressWarnings(cor(z, ov$age, method = "spearman")))
loo <- vapply(seq_len(nrow(logexpr)), function(i) {
  signs <- vapply(seq_len(ncol(logexpr)), function(j) sign(suppressWarnings(cor(logexpr[i, -j], ov$age[-j], method = "spearman"))), numeric(1))
  mean(signs == sign(rho[i]), na.rm = TRUE)
}, numeric(1))
effects <- data.table(human_gene_id_stable = rownames(logexpr), hpa_age_spearman = rho,
                      hpa_leave_one_out_direction_stability = loo, n_hpa_samples = ncol(logexpr), exploratory = TRUE)
saveRDS(list(tpm = mat, log2_tpm_plus1 = logexpr, metadata = ov, age_effects = effects,
             release = config$hpa$release, ensembl_version = config$hpa$ensembl_version),
        path_module("derived_data/hpa", "hpa_ovary_prepared.rds"))
write_tsv(ov, path_module("derived_data/hpa", "hpa_ovary_sample_metadata.tsv"))
write_tsv(effects, path_module("derived_data/hpa", "hpa_gene_age_effects_exploratory.tsv"))
message("HPA exploratory ovary data prepared: seven exact-age samples.")
