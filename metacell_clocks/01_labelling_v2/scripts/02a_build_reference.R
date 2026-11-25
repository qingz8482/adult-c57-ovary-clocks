# Step B2: build reference objects from author-labelled cells (Morris, Isola).
# Outputs v2/02_reference_transfer/ref_{combined,morris,isola}.rds
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
set.seed(20261020)
MAX_PER_LABEL <- 1500
out_dir <- file.path(V2, "02_reference_transfer")
dir.create(out_dir, showWarnings = FALSE)

amap <- author_map()
load_ref <- function(ds) {
  o <- load_query(ds)
  lab <- as.character(o$author_cell_type)
  m <- amap[dataset_id == ds]
  o$ref_broad <- m$broad_class[match(lab, m$author_label)]
  o$ref_fine <- m$fine_state[match(lab, m$author_label)]
  o$ref_label <- paste(o$ref_broad, o$ref_fine, sep = "|")
  o$ref_dataset <- ds
  unmapped <- unique(lab[is.na(o$ref_broad)])
  if (length(unmapped)) message("Unmapped labels in ", ds, ": ", paste(unmapped, collapse = ", "))
  use <- m$use_in_reference[match(lab, m$author_label)]
  o <- subset(o, cells = colnames(o)[!is.na(o$ref_broad) & use %in% TRUE])
  # downsample per label
  keep <- unlist(lapply(split(colnames(o), o$ref_label), function(x) if (length(x) > MAX_PER_LABEL) sample(x, MAX_PER_LABEL) else x))
  subset(o, cells = keep)
}

process_ref <- function(o) {
  o <- NormalizeData(o, verbose = FALSE)
  o <- FindVariableFeatures(o, nfeatures = 3000, verbose = FALSE)
  o <- ScaleData(o, verbose = FALSE)
  RunPCA(o, npcs = 50, verbose = FALSE)
}

morris <- load_ref(AUTHOR_DATASETS[1])
isola <- load_ref(AUTHOR_DATASETS[2])
common <- intersect(rownames(morris), rownames(isola))
message("Common genes Morris/Isola: ", length(common))

ref_m <- process_ref(morris)
ref_i <- process_ref(isola)
saveRDS(ref_m, file.path(out_dir, "ref_morris.rds"))
saveRDS(ref_i, file.path(out_dir, "ref_isola.rds"))

comb <- merge(subset(morris, features = common), subset(isola, features = common))
comb[["RNA"]] <- JoinLayers(comb[["RNA"]])
ref_c <- process_ref(comb)
saveRDS(ref_c, file.path(out_dir, "ref_combined.rds"))

summ <- rbind(
  data.table(ref = "morris", as.data.table(table(ref_label = ref_m$ref_label))),
  data.table(ref = "isola", as.data.table(table(ref_label = ref_i$ref_label))),
  data.table(ref = "combined", as.data.table(table(ref_label = ref_c$ref_label)))
)
fwrite(summ, file.path(out_dir, "reference_label_counts.tsv"), sep = "\t")
message("done")
