# Step B2: transfer reference labels to one query dataset.
# Usage: Rscript 02b_reference_transfer.R <dataset_id> <ref: combined|morris|isola>
args <- commandArgs(trailingOnly = TRUE)
ds <- args[[1]]
ref_name <- args[[2]]
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
set.seed(20261020)
CHUNK <- 40000
out_dir <- file.path(V2, "02_reference_transfer")

ref <- readRDS(file.path(out_dir, paste0("ref_", ref_name, ".rds")))
q <- load_query(ds)
q <- NormalizeData(q, verbose = FALSE)
features <- intersect(VariableFeatures(ref), rownames(q))
message(ds, " <- ref_", ref_name, ": ", length(features), " shared variable features, ", ncol(q), " cells")

cells <- colnames(q)
chunks <- split(cells, ceiling(seq_along(cells) / CHUNK))
res <- rbindlist(lapply(chunks, function(cc) {
  qq <- subset(q, cells = cc)
  anchors <- FindTransferAnchors(reference = ref, query = qq, dims = 1:30, features = features,
                                 reference.reduction = "pca", normalization.method = "LogNormalize",
                                 verbose = FALSE)
  pred <- TransferData(anchorset = anchors, refdata = ref$ref_label, dims = 1:30, verbose = FALSE)
  sc <- as.matrix(pred[, grep("^prediction.score\\.", colnames(pred)) , drop = FALSE])
  sc <- sc[, colnames(sc) != "prediction.score.max", drop = FALSE]
  ord <- t(apply(sc, 1, sort, decreasing = TRUE))
  # broad-level score: sum of fine-label scores within broad class
  labs <- sub("^prediction.score.", "", colnames(sc))
  labs <- gsub(".", "|", labs, fixed = TRUE)  # Seurat replaces '|' with '.'
  broad_of <- sub("\\|.*$", "", labs)
  broad_sc <- sapply(split(seq_along(broad_of), broad_of), function(ix) rowSums(sc[, ix, drop = FALSE]))
  if (is.null(dim(broad_sc))) broad_sc <- matrix(broad_sc, nrow = 1, dimnames = list(rownames(sc), names(broad_sc)))
  data.table(
    cell_id = rownames(pred), dataset_id = ds, reference = ref_name,
    ref_pred_label = pred$predicted.id,
    ref_pred_score = pred$prediction.score.max,
    ref_pred_margin = ord[, 1] - ord[, 2],
    ref_pred_broad = colnames(broad_sc)[max.col(broad_sc, ties.method = "first")],
    ref_pred_broad_score = apply(broad_sc, 1, max)
  )
}))
res[, ref_pred_fine := sub("^.*\\|", "", ref_pred_label)]
fwrite(res, file.path(out_dir, paste0("transfer_", ds, "__ref_", ref_name, ".tsv.gz")), sep = "\t")
message("done ", ds)
