#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))

raw_path <- path_module("derived_data/orthology", paste0("ensembl_release_", config$ensembl$release, "_mouse_human_biomart.tsv"))
xml_path <- path_module("derived_data/orthology", "biomart_query.xml")
if (!file.exists(raw_path) || file.info(raw_path)$size == 0) {
  xml <- paste0('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="mmusculus_gene_ensembl" interface="default">',
                '<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>',
                '<Attribute name="hsapiens_homolog_ensembl_gene"/><Attribute name="hsapiens_homolog_associated_gene_name"/>',
                '<Attribute name="hsapiens_homolog_orthology_type"/><Attribute name="hsapiens_homolog_orthology_confidence"/>',
                '<Attribute name="hsapiens_homolog_perc_id"/><Attribute name="hsapiens_homolog_perc_id_r1"/>',
                '</Dataset></Query>')
  writeLines(xml, xml_path)
  rc <- system2("curl", c("-fsSL", "--retry", "3", "--get", "--data-urlencode", shQuote(paste0("query@", xml_path)),
                          shQuote(config$ensembl$mart_url), "-o", shQuote(raw_path)))
  if (rc != 0 || !file.exists(raw_path) || file.info(raw_path)$size == 0) stop("Ensembl BioMart download failed.")
}

x <- fread(raw_path)
if (ncol(x) != 8L) stop("Unexpected BioMart schema.")
setnames(x, c("mouse_gene_id_stable", "mouse_gene_symbol", "human_gene_id_stable", "human_gene_symbol",
              "orthology_type", "orthology_confidence", "mouse_percent_identity", "human_percent_identity"))
x[, `:=`(mouse_gene_id_original = mouse_gene_id_stable,
          ensembl_release = as.integer(config$ensembl$release), mapping_source = "Ensembl BioMart",
          mapping_status = fifelse(human_gene_id_stable == "", "unmapped",
                           fifelse(orthology_type == "ortholog_one2one" & orthology_confidence == 1, "primary_one2one_high_confidence", orthology_type)))]
x <- x[, .(mouse_gene_id_original, mouse_gene_id_stable, mouse_gene_symbol, human_gene_id_stable,
           human_gene_symbol, orthology_type, orthology_confidence, mouse_percent_identity,
           human_percent_identity, ensembl_release, mapping_source, mapping_status)]
write_tsv(x, path_module("derived_data/orthology", "mouse_human_orthology_full.tsv"))
primary <- x[mapping_status == "primary_one2one_high_confidence" & human_gene_id_stable != ""]
primary <- primary[!duplicated(mouse_gene_id_stable) & !duplicated(human_gene_id_stable)]
write_tsv(primary, path_module("derived_data/orthology", "mouse_human_orthology_primary_one2one.tsv"))

master <- fread(path_module("derived_data/mouse_models", "mouse_gene_direction_and_weight_master.tsv"))
audit <- merge(master, x, by = "mouse_gene_symbol", all.x = TRUE, allow.cartesian = TRUE)
audit[is.na(mapping_status), mapping_status := "unmapped"]
sets <- list(integrated_universe = audit$integrated_gene_universe,
             elasticnet_hvg_2000 = audit$elasticnet_hvg_2000,
             pca_versionB_2000 = audit$pca_versionB_2000,
             deployment_shared_artifact = audit$deployment_shared_artifact,
             elasticnet_nonzero = audit$elasticnet_nonzero %in% TRUE)
summary <- rbindlist(lapply(names(sets), function(s) {
  d <- audit[sets[[s]] %in% TRUE]
  data.table(gene_set = s, mapping_status = d$mapping_status)[, .N, by = .(gene_set, mapping_status)]
}))
write_tsv(summary, path_module("results/qc", "orthology_mapping_summary.tsv"))
write_tsv(data.table(ensembl_release = config$ensembl$release, source = "Ensembl BioMart",
                     url = config$ensembl$mart_url, access_date = as.character(Sys.Date()), raw_file = normalizePath(raw_path)),
          path_module("data_manifest", "orthology_source_manifest.tsv"))
message("Official Ensembl orthology mapping complete.")
