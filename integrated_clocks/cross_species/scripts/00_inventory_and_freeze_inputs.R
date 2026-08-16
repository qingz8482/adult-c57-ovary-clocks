#!/usr/bin/env Rscript
source(file.path(dirname(normalizePath(commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1] |> sub("^--file=", "", x = _))), "utilities.R"))

download_one <- function(label, url, dest, release) {
  ensure_dir(dirname(dest))
  status <- "cached"
  if (!file.exists(dest) || file.info(dest)$size == 0) {
    status <- "downloaded"
    rc <- system2("curl", c("-fL", "--retry", "3", "--retry-delay", "2", "-o", shQuote(dest), shQuote(url)))
    if (rc != 0 || !file.exists(dest) || file.info(dest)$size == 0) stop("Download failed: ", url)
  }
  data.table(label, release, url, local_path = normalizePath(dest), file_size = file.info(dest)$size,
             access_date = as.character(Sys.Date()), status)
}

artifacts <- list.files(path_repo(config$mouse$artifact_dir), pattern = "deployment_all_control_like[.]rds$",
                        full.names = TRUE, recursive = FALSE)
artifacts <- artifacts[grepl("integrated_bulk_pseudobulk__", basename(artifacts))]
if (length(artifacts) != 2L) stop("Expected exactly two authoritative top-level deployment artifacts; found ", length(artifacts))

input_paths <- unique(c(
  artifacts,
  path_repo(config$mouse$alignment_file),
  path_repo(config$mouse$authoritative_results, "FROZEN_RESULTS_INSTRUCTIONS.md"),
  path_repo(config$mouse$authoritative_results, "model_artifacts", "model_artifact_manifest.tsv")
))
input_manifest <- data.table(
  input_type = c(rep("mouse_frozen_input", length(input_paths))),
  path = normalizePath(input_paths),
  size_bytes = file.info(input_paths)$size,
  modified_time = as.character(file.info(input_paths)$mtime)
)
write_tsv(input_manifest, path_module("data_manifest", "input_file_manifest.tsv"))

downloads <- rbindlist(list(
  download_one("gtex_ovary_gene_reads", config$gtex$counts_url,
               path_module("derived_data/gtex/raw", basename(config$gtex$counts_url)), config$gtex$release),
  download_one("gtex_sample_attributes", config$gtex$sample_attributes_url,
               path_module("derived_data/gtex/raw", basename(config$gtex$sample_attributes_url)), config$gtex$release),
  download_one("gtex_subject_phenotypes", config$gtex$subject_phenotypes_url,
               path_module("derived_data/gtex/raw", basename(config$gtex$subject_phenotypes_url)), config$gtex$release),
  download_one("hpa_transcript_tissue", config$hpa$transcript_url,
               path_module("derived_data/hpa/raw", basename(config$hpa$transcript_url)), config$hpa$release),
  download_one("hpa_tissue_samples", config$hpa$samples_url,
               path_module("derived_data/hpa/raw", basename(config$hpa$samples_url)), config$hpa$release)
), fill = TRUE)
write_tsv(downloads, path_module("data_manifest", "download_manifest.tsv"))

session <- capture.output(sessionInfo())
writeLines(session, path_module("data_manifest", "software_versions.txt"))

all_files <- c(input_manifest$path, downloads$local_path)
cmd <- c("sha256sum", vapply(all_files, shQuote, character(1)), ">", shQuote(path_module("data_manifest", "checksums.sha256")))
rc <- system(paste(cmd, collapse = " "))
if (rc != 0) stop("sha256sum failed")
message("Inventory and public input acquisition complete.")
