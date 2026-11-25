# Shared helpers for the v2 label rebuild. Reads only the copies in v2/00_input_seurat.
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

V2 <- file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2")
INPUT_DIR <- file.path(V2, "00_input_seurat")

DATASETS <- c(
  "2022_eLife_Morris_cycling_murine_ovary_atlas",
  "2024_NatureAging_Isola_aging_mouse_ovary_atlas",
  "2024_LifeMedicine_Jiang_whole_life_mouse_ovary",
  "2025_GeroScience_Zhang_multiomic_aging_mouse_ovary",
  "2026_CommBio_Galligos_ovarian_inflammaging",
  "2025_CommBio_Zhao_developing_ovary_scRNA_spatial",
  "2024_Cell_Winkler_FRT_scRNA_spatial",
  "2025_Kim_PRJNA863443_mouse_menopause_models"
)
AUTHOR_DATASETS <- DATASETS[1:2]

kim_features <- function() {
  fm <- read.delim(gzfile(file.path(INPUT_DIR, "kim_features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  setNames(fm$V2, fm$V1)
}

# fread cannot read .gz without R.utils; stream through zcat instead.
fread_gz <- function(path, ...) fread(cmd = paste("zcat", shQuote(path)), ...)

# Repair R-mangled symbols (leading "X" before a digit, "." for "-") against a reference universe.
fix_symbols <- function(x, universe) {
  out <- x
  bad <- !(x %in% universe)
  cand1 <- sub("^X(?=[0-9])", "", x, perl = TRUE)
  cand2 <- gsub(".", "-", x, fixed = TRUE)
  cand3 <- gsub(".", "-", cand1, fixed = TRUE)
  for (cand in list(cand1, cand2, cand3)) {
    hit <- bad & cand %in% universe
    out[hit] <- cand[hit]
    bad <- bad & !hit
  }
  make.unique(out)
}

rebuild_with_names <- function(o, new_names) {
  counts <- LayerData(o, assay = "RNA", layer = "counts")
  rownames(counts) <- new_names
  o2 <- CreateSeuratObject(counts = counts, meta.data = o@meta.data, min.cells = 0, min.features = 0)
  o2
}

# Load a copied object, join split layers, convert to gene symbols, ensure log-normalised data.
load_query <- function(ds) {
  o <- readRDS(file.path(INPUT_DIR, paste0(ds, ".processed.seurat.rds")))
  DefaultAssay(o) <- "RNA"
  if (inherits(o[["RNA"]], "Assay5") && length(Layers(o[["RNA"]], search = "counts")) > 1) {
    o[["RNA"]] <- JoinLayers(o[["RNA"]])
  }
  universe <- unique(kim_features())
  rn <- rownames(o)
  if (grepl("^ENSMUSG", rn[1])) {
    km <- kim_features()
    o <- rebuild_with_names(o, make.unique(unname(km[rn])))
  } else if (mean(rn %in% universe) < 0.98) {
    o <- rebuild_with_names(o, fix_symbols(rn, universe))
  }
  if (!"data" %in% Layers(o[["RNA"]])) {
    o <- NormalizeData(o, verbose = FALSE)
  }
  o
}

panel_v2 <- function() fread(file.path(V2, "01_marker_panel_v2.tsv"))
author_map <- function() fread(file.path(V2, "03_author_label_map.tsv"))

cc_genes_mouse <- function() {
  to_mouse <- function(g) paste0(substr(g, 1, 1), tolower(substr(g, 2, nchar(g))))
  list(s = to_mouse(Seurat::cc.genes.updated.2019$s.genes),
       g2m = to_mouse(Seurat::cc.genes.updated.2019$g2m.genes))
}

# Extra proliferation genes seen in DE lists, used only to flag cycling-dominated clusters.
PROLIF_EXTRA <- c("Mki67", "Top2a", "Pbk", "Cdc20", "Ccnb1", "Ccnb2", "Ube2c", "Birc5", "Cenpf", "Hmmr",
                  "Depdc1a", "Kif20a", "Cdca2", "Troap", "Cdc25c", "Kif18b", "Cep55", "Nek2", "Mxd3",
                  "Pimreg", "Neil3", "Shcbp1", "Ska1", "Nuf2", "Esco2", "Knl1", "Ckap2l", "Kif15", "Bub1",
                  "Bub1b", "Melk", "E2f8", "Ccna2", "Foxm1", "Sapcd2", "Sgo2a", "Oip5", "Kifc1", "Eme1",
                  "Mis18bp1", "Ankle1", "Arhgef39", "Iqgap3", "Pclaf", "Cdk1", "Kif11", "Prc1")

is_prolif_gene <- function(g) {
  cc <- cc_genes_mouse()
  g %in% c(cc$s, cc$g2m, PROLIF_EXTRA) | grepl("^Hist[0-9]|^H2a[c-z]|^H2b[c-z]|^H1f|^H4c|^H3c", g)
}
