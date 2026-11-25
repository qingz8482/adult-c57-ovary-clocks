# v2 cell labelling for an external library (GSE235021): same marker panel v2.2, gate, reference
# (ref_combined) and consensus rules as v2/scripts/01, 02b, 04. No manual curation is available for external
# libraries, so cluster-level conflicts become Unknown.
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
V2_DIR <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")

.v2_panel <- fread(file.path(V2_DIR, "01_marker_panel_v2.tsv"))
.v2_cal <- fread(file.path(V2_DIR, "03_benchmark", "calibration_reference_pooled.tsv"))
.v2_bench <- fread(file.path(V2_DIR, "03_benchmark", "benchmark_per_class.tsv"))
.v2_mprec <- .v2_bench[method == "marker", .(precision = mean(precision)), by = class]
.v2_ref <- NULL

v2_ref_prec <- function(score) {
  bin <- as.character(cut(score, c(0, 0.5, 0.7, 0.9, 1.0001), right = FALSE,
                          labels = c("[0,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1]")))
  .v2_cal$precision[match(bin, as.character(.v2_cal$score_bin))]
}
v2_mk_prec <- function(cls) { p <- .v2_mprec$precision[match(cls, .v2_mprec$class)]; ifelse(is.na(p), 0.8, p) }

v2_label_library <- function(object, sample_id) {
  if (is.null(.v2_ref)) .v2_ref <<- readRDS(file.path(V2_DIR, "02_reference_transfer", "ref_combined.rds"))
  X <- GetAssayData(object, assay = "RNA", layer = "data")
  cl <- as.character(Idents(object)); clusters <- sort(unique(cl)); cells <- colnames(X)
  G <- sparseMatrix(i = seq_along(cl), j = match(cl, clusters), x = 1, dims = c(length(cl), length(clusters)))
  n_cl <- colSums(G)
  mean_cl <- sweep(as.matrix(X %*% G), 2, n_cl, "/")
  det_cl <- sweep(as.matrix((X > 0) %*% G), 2, n_cl, "/")
  colnames(mean_cl) <- colnames(det_cl) <- clusters
  z_gene <- t(scale(t(mean_cl))); z_gene[is.na(z_gene)] <- 0
  info <- unique(.v2_panel[, .(signature, broad_class, fine_state)])
  long <- list(); cs <- matrix(0, length(cells), nrow(info), dimnames = list(cells, info$signature))
  for (s in info$signature) {
    g <- intersect(.v2_panel[signature == s, gene], rownames(X)); if (!length(g)) next
    n_det <- colSums(det_cl[g, , drop = FALSE] >= 0.20)
    long[[s]] <- data.table(sample_id = sample_id, seurat_cluster = clusters, signature = s,
                            z_score = colMeans(z_gene[g, , drop = FALSE]),
                            gate_pass = n_det / length(g) >= 0.35 & n_det >= min(3, ceiling(length(g) / 2)))
    cs[, s] <- Matrix::colMeans(X[g, , drop = FALSE])
  }
  long <- merge(rbindlist(long), info, by = "signature")
  csz <- scale(cs); csz[is.na(csz)] <- 0
  cell_top <- colnames(csz)[max.col(csz, ties.method = "first")]

  # reference transfer (same parameters as 02b)
  feats <- intersect(VariableFeatures(.v2_ref), rownames(object))
  anchors <- FindTransferAnchors(reference = .v2_ref, query = object, dims = 1:30, features = feats,
                                 reference.reduction = "pca", normalization.method = "LogNormalize", verbose = FALSE)
  pred <- TransferData(anchorset = anchors, refdata = .v2_ref$ref_label, dims = 1:30, verbose = FALSE)
  sc <- as.matrix(pred[, setdiff(grep("^prediction.score\\.", colnames(pred), value = TRUE), "prediction.score.max"), drop = FALSE])
  labs <- gsub(".", "|", sub("^prediction.score.", "", colnames(sc)), fixed = TRUE)
  broad_of <- sub("\\|.*$", "", labs)
  bsc <- sapply(split(seq_along(broad_of), broad_of), function(ix) rowSums(sc[, ix, drop = FALSE]))
  x <- data.table(cell_id = cells, seurat_cluster = cl, cell_top_signature = cell_top,
                  ref_pred_broad = colnames(bsc)[max.col(bsc, ties.method = "first")],
                  ref_pred_broad_score = apply(bsc, 1, max), ref_pred_label = pred$predicted.id)
  x[, ref_pred_fine := sub("^.*\\|", "", ref_pred_label)]
  x[, ref_ok := ref_pred_broad_score >= 0.5]

  calls <- long[, {
    gated <- .SD[gate_pass == TRUE & z_score > 0][order(-z_score)]
    .(marker_broad = if (nrow(gated)) gated$broad_class[1] else "Unknown",
      marker_fine = if (nrow(gated)) gated$fine_state[1] else "unknown")
  }, by = seurat_cluster]
  x <- merge(x, calls, by = "seurat_cluster")
  clt <- x[, {
    rb <- ref_pred_broad[ref_ok]; tb <- if (length(rb)) sort(table(rb), decreasing = TRUE) else integer()
    .(ref_major = if (length(tb)) names(tb)[1] else NA_character_,
      ref_major_frac = if (length(tb)) as.numeric(tb[1]) / .N else 0, marker_broad = marker_broad[1])
  }, by = seurat_cluster]
  clt[, case := fcase(marker_broad == "Erythroid", "marker_erythroid",
                      marker_broad != "Unknown" & marker_broad == ref_major & ref_major_frac >= 0.5, "agree",
                      marker_broad == "Unknown" & ref_major_frac >= 0.7, "reference_only",
                      marker_broad != "Unknown" & ref_major_frac < 0.5 & (is.na(ref_major) | ref_major == marker_broad), "marker_only",
                      default = "conflict")]
  clt[, cluster_label := fcase(case %in% c("agree", "marker_only", "marker_erythroid"), marker_broad,
                               case == "reference_only", ref_major, default = "Unknown")]
  x <- merge(x, clt[, .(seurat_cluster, case, cluster_label)], by = "seurat_cluster")
  x[, `:=`(label = cluster_label, status = "probable", conf = 0)]
  veto <- x$label != "Unknown" & x$label != "Erythroid" & x$ref_pred_broad_score >= 0.8 & x$ref_pred_broad != x$label
  support <- !veto & x$ref_ok & x$ref_pred_broad == x$label
  x[veto, `:=`(label = "Unknown")]
  x[support & case == "agree", `:=`(status = "confirmed", conf = pmin(0.99, 1 - (1 - v2_ref_prec(ref_pred_broad_score)) * (1 - v2_mk_prec(label))))]
  x[support & case %in% c("reference_only", "marker_only"), conf := v2_ref_prec(ref_pred_broad_score)]
  x[!veto & !support & case %in% c("agree", "marker_only", "marker_erythroid"), conf := v2_mk_prec(label)]
  x[!veto & !support & case == "reference_only", label := "Unknown"]
  x[label == "Unknown", `:=`(status = "ambiguous", conf = 0)]

  list(
    scores = long[, .(sample_id, seurat_cluster, candidate_cell_type = signature, marker_score = z_score, gate_pass)],
    labels = clt[, .(sample_id = sample_id, seurat_cluster, harmonized_cell_class = cluster_label,
                     consensus_case = case, marker_broad, ref_major, ref_major_frac,
                     annotation_method = "v2_marker_reference_consensus")],
    cells = x[, .(cell_id, label, status, conf, case, ref_pred_broad, ref_pred_broad_score, cell_top_signature)]
  )
}
