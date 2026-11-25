# Step C1/C2: marker panel v2 scoring for one dataset.
# Usage: Rscript 01_marker_v2_scoring.R <dataset_id>
args <- commandArgs(trailingOnly = TRUE)
ds <- args[[1]]
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
set.seed(20261020)

# v2.1 (after benchmark B3, see REVISIONS.md): gate = >=3 genes and >=35% of genes detected; z > 0.
MIN_Z <- 0              # leading signature must be above the dataset's average cluster
GATE_DETECT <- 0.20     # a gene counts as detected if expressed in >= 20% of cluster cells
GATE_FRAC <- 0.35       # >= 35% of the signature's genes must be detected
GATE_MIN_GENES <- 3     # and at least min(3, ceiling(n_genes / 2)) genes detected (v2.2)
CYCLING_TOP_FRAC <- 0.4 # >= 40% of top-15 DE genes proliferation-related -> cycling-dominated

out_dir <- file.path(V2, "01_marker_v2", ds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

o <- load_query(ds)
md <- o@meta.data
cl <- as.character(md$standard_cluster)
cluster_source <- "standard_cluster"
if (all(is.na(cl) | cl %in% c("NA", ""))) {
  message("No standard_cluster: reclustering ", ds)
  o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE)
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- FindNeighbors(o, dims = 1:30, verbose = FALSE)
  o <- FindClusters(o, resolution = 0.5, random.seed = 20261020, verbose = FALSE)
  cl <- as.character(o$seurat_clusters)
  cluster_source <- "v2_recluster_res0.5"
}
X <- LayerData(o, assay = "RNA", layer = "data")
cells <- colnames(X)
clusters <- sort(unique(cl))
G <- sparseMatrix(i = seq_along(cl), j = match(cl, clusters), x = 1,
                  dims = c(length(cl), length(clusters)), dimnames = list(cells, clusters))
n_cl <- colSums(G)
sum_cl <- as.matrix(X %*% G)
mean_cl <- sweep(sum_cl, 2, n_cl, "/")
det_cl <- sweep(as.matrix((X > 0) %*% G), 2, n_cl, "/")

# Signature scores ----------------------------------------------------------
panel <- panel_v2()
sigs <- unique(panel$signature)
sig_info <- unique(panel[, .(signature, broad_class, fine_state)])
z_gene <- t(scale(t(mean_cl)))
z_gene[is.na(z_gene)] <- 0
cl_scores <- list()
cell_scores <- matrix(NA_real_, nrow = length(cells), ncol = length(sigs), dimnames = list(cells, sigs))
for (s in sigs) {
  genes <- intersect(panel[signature == s, gene], rownames(X))
  if (length(genes) == 0) next
  z <- colMeans(z_gene[genes, , drop = FALSE])
  n_det <- colSums(det_cl[genes, , drop = FALSE] >= GATE_DETECT)
  gate_frac <- n_det / length(genes)
  cl_scores[[s]] <- data.table(cluster = clusters, signature = s, z_score = z,
                               abs_mean = colMeans(mean_cl[genes, , drop = FALSE]),
                               gate_frac = gate_frac,
                               gate_pass = gate_frac >= GATE_FRAC & n_det >= min(GATE_MIN_GENES, ceiling(length(genes) / 2)),
                               n_genes_present = length(genes))
  cell_scores[, s] <- Matrix::colMeans(X[genes, , drop = FALSE])
}
cl_long <- rbindlist(cl_scores)
cl_long <- merge(cl_long, sig_info, by = "signature")
fwrite(cl_long, file.path(out_dir, "cluster_signature_scores_long.tsv"), sep = "\t")

# Cell cycle ------------------------------------------------------------------
cc <- cc_genes_mouse()
o <- CellCycleScoring(o, s.features = intersect(cc$s, rownames(o)),
                      g2m.features = intersect(cc$g2m, rownames(o)), set.ident = FALSE)

# Per-cell argmax signature (z per signature across cells) -> cluster purity
cs <- scale(cell_scores)
cs[is.na(cs)] <- 0
cell_top_sig <- colnames(cs)[max.col(cs, ties.method = "first")]
cell_top_broad <- sig_info$broad_class[match(cell_top_sig, sig_info$signature)]

# Top DE genes (mean log-normalised difference vs rest, detection >= 25%) ----
tot_sum <- rowSums(sum_cl)
N <- sum(n_cl)
top_de <- sapply(clusters, function(k) {
  rest <- (tot_sum - sum_cl[, k]) / (N - n_cl[[k]])
  d <- mean_cl[, k] - rest
  d[det_cl[, k] < 0.25] <- -Inf
  names(sort(d, decreasing = TRUE))[1:15]
}, simplify = FALSE)

# Cluster calls -----------------------------------------------------------------
calls <- rbindlist(lapply(clusters, function(k) {
  d <- cl_long[cluster == k][order(-z_score)]
  gated <- d[gate_pass == TRUE & z_score > MIN_Z]
  best <- if (nrow(gated) > 0) gated[1] else NULL
  second_z <- if (!is.null(best)) max(d[signature != best$signature, z_score]) else NA_real_
  de <- top_de[[k]]
  idx <- which(cl == k)
  label <- if (!is.null(best)) best$broad_class else "Unknown"
  data.table(
    dataset_id = ds, cluster = k, cluster_source = cluster_source, n_cells = n_cl[[k]],
    marker_signature = if (!is.null(best)) best$signature else NA_character_,
    marker_broad = label,
    marker_fine = if (!is.null(best)) best$fine_state else "unknown",
    marker_z = if (!is.null(best)) best$z_score else d$z_score[1],
    marker_margin = if (!is.null(best)) best$z_score - second_z else NA_real_,
    top_ungated_signature = d$signature[1],
    n_signatures_gated = nrow(gated),
    cell_purity = mean(cell_top_broad[idx] == label),
    mean_S = mean(o$S.Score[idx]), mean_G2M = mean(o$G2M.Score[idx]),
    frac_cycling_cells = mean(o$Phase[idx] != "G1"),
    prolif_top15_frac = mean(is_prolif_gene(de)),
    cycling_dominated = mean(is_prolif_gene(de)) >= CYCLING_TOP_FRAC,
    top_de_genes = paste(de, collapse = ",")
  )
}))
fwrite(calls, file.path(out_dir, "cluster_marker_calls.tsv"), sep = "\t")

cell_tbl <- data.table(
  cell_id = cells, dataset_id = ds, sample_id = as.character(md[cells, "sample_id"]),
  cluster = cl, cluster_source = cluster_source,
  cell_top_signature = cell_top_sig, cell_top_broad = cell_top_broad,
  S_score = o$S.Score, G2M_score = o$G2M.Score, phase = o$Phase
)
if ("author_cell_type" %in% colnames(md)) cell_tbl[, author_cell_type := as.character(md[cells, "author_cell_type"])]
fwrite(cell_tbl, file.path(out_dir, "cell_marker_table.tsv.gz"), sep = "\t")
message("done ", ds, ": ", length(clusters), " clusters, ", length(cells), " cells")
