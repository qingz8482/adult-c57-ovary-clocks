# Step C5: the same transcriptional state must get the same label in every dataset.
# Cluster mean profiles on the reference variable genes, centred within dataset (removes dataset offset),
# then each cluster's best-correlated cluster in every other dataset is compared on final label.
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
out_dir <- file.path(V2, "04_consensus")
SIM_MIN <- 0.6

genes <- VariableFeatures(readRDS(file.path(V2, "02_reference_transfer", "ref_combined.rds")))
cons <- fread_gz(file.path(out_dir, "harmonized_annotation_consensus_v2.tsv.gz"), colClasses = list(character = "standard_cluster"))
clus <- fread(file.path(out_dir, "cluster_consensus_v2.tsv"), colClasses = list(character = "cluster"))

prof_path <- file.path(out_dir, "cluster_profiles_refgenes.rds")
if (file.exists(prof_path)) {
  profiles <- readRDS(prof_path)
} else {
  profiles <- list()
  for (ds in DATASETS) {
    o <- load_query(ds)
    o <- NormalizeData(o, verbose = FALSE)
    X <- LayerData(o, layer = "data")
    cl <- cons[dataset_id == ds][match(colnames(X), cell_id), standard_cluster]
    g <- intersect(genes, rownames(X))
    G <- sparseMatrix(i = seq_along(cl), j = match(cl, sort(unique(cl))), x = 1)
    M <- as.matrix(X[g, ] %*% G)
    M <- sweep(M, 2, colSums(G), "/")
    colnames(M) <- paste(ds, sort(unique(cl)), sep = "::")
    full <- matrix(NA_real_, length(genes), ncol(M), dimnames = list(genes, colnames(M)))
    full[g, ] <- M
    profiles[[ds]] <- full
    rm(o, X); gc(verbose = FALSE)
  }
  saveRDS(profiles, prof_path)
}

# centre each gene within dataset, keep genes present everywhere
P <- do.call(cbind, lapply(profiles, function(M) M - rowMeans(M)))
P <- P[stats::complete.cases(P), ]
C <- cor(P, method = "spearman")

lab <- clus[, .(key = paste(dataset_id, cluster, sep = "::"), cluster_label, n_cells, cluster_case)]
ds_of <- sub("::.*$", "", colnames(C))
rows <- list()
for (k in colnames(C)) {
  for (d2 in setdiff(unique(ds_of), ds_of[colnames(C) == k])) {
    cand <- colnames(C)[ds_of == d2]
    best <- cand[which.max(C[k, cand])]
    rows[[length(rows) + 1]] <- data.table(cluster_a = k, cluster_b = best, similarity = C[k, best])
  }
}
nn <- rbindlist(rows)
nn <- merge(nn, lab[, .(cluster_a = key, label_a = cluster_label, n_a = n_cells, case_a = cluster_case)], by = "cluster_a")
nn <- merge(nn, lab[, .(cluster_b = key, label_b = cluster_label, n_b = n_cells)], by = "cluster_b")
# reciprocal best hits only
nn[, pair := paste(pmin(cluster_a, cluster_b), pmax(cluster_a, cluster_b))]
recip <- nn[, .N, by = pair][N == 2, pair]
nn[, reciprocal := pair %in% recip]
nn[, consistent := label_a == label_b]
fwrite(nn[order(-similarity)], file.path(out_dir, "cross_dataset_nearest_clusters.tsv"), sep = "\t")
flag <- unique(nn[reciprocal & similarity >= SIM_MIN & !consistent & label_a != "Unknown" & label_b != "Unknown",
                  .(pair, cluster_a, label_a, n_a, cluster_b, label_b, n_b, similarity)], by = "pair")
fwrite(flag[order(-similarity)], file.path(out_dir, "cross_dataset_inconsistencies.tsv"), sep = "\t")
summ <- nn[reciprocal & similarity >= SIM_MIN, .(n_pairs = uniqueN(pair), n_inconsistent = uniqueN(pair[!consistent]))]
fwrite(summ, file.path(out_dir, "cross_dataset_consistency_summary.tsv"), sep = "\t")
print(summ); print(flag[order(-similarity)])
