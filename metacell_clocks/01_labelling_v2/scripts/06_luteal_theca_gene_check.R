# Diagnostic for v2.2: detection / mean of candidate luteal vs theca genes by author label (Morris, Isola).
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
genes <- c("Akr1c18","Prlr","Ptgfr","Plin5","Sfrp4","Cyp11a1","Star","Hsd3b1","Lhcgr","Cyp17a1","Sct","Bpifb4","Mgarp",
           "Aldh1a1","Ephx2","S100a6","Gm2a","Runx2","Cd36","Lgmn","Tinagl1","Prss35","Hsd17b7","Scarb1","Ldlr","Fdx1",
           "Ptgds","Hhip","Mro","Lrp8","Nr4a1","Cyp19a1","Inha","Foxl2","Ptgs2","Car2","Vnn1","Tnfaip6")
out <- list()
for (ds in AUTHOR_DATASETS) {
  o <- load_query(ds); X <- LayerData(o, layer = "data")
  g <- intersect(genes, rownames(X)); lab <- o$author_cell_type
  for (l in unique(lab)) {
    idx <- which(lab == l)
    out[[length(out)+1]] <- data.table(dataset = substr(ds, 1, 18), label = l, n = length(idx), gene = g,
      det = Matrix::rowMeans(X[g, idx, drop = FALSE] > 0), mean = Matrix::rowMeans(X[g, idx, drop = FALSE]))
  }
}
d <- rbindlist(out)
fwrite(d, file.path(V2, "03_benchmark", "luteal_theca_candidate_gene_detection.tsv"), sep = "\t")
keep <- c("GC_CL_Active","GC_CL_Lytic","GC_Estrous","GC_Antral","GC_Mural","M_Steroidogenic Theca","M_Immature Theca","M_Cortical Stroma","Luteal","Theca","Granulosa A","Stroma A","Stroma B")
w <- dcast(d[label %in% keep], gene ~ paste(dataset, label), value.var = "det")
options(width = 250); print(w, digits = 2)
