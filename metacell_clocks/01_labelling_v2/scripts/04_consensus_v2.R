# Step C3/C4: combine author labels, reference transfer and marker v2 into the v2 consensus.
# Curated decisions in 05_curation/cluster_decisions.tsv override automatic cluster calls.
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
out_dir <- file.path(V2, "04_consensus")
`%||%` <- function(a, b) if (is.null(a)) b else a
dir.create(out_dir, showWarnings = FALSE)

REF_MIN <- 0.5      # reference broad score counted as evidence
REF_STRONG <- 0.8   # confident reference call that can veto a cluster label for a single cell
AGREE_FRAC <- 0.5   # reference majority needed for cluster agreement
REF_ONLY_FRAC <- 0.7

amap <- author_map()
ont <- fread(file.path(V2, "00_ontology_v2.tsv"))
cal <- fread(file.path(V2, "03_benchmark", "calibration_reference_pooled.tsv"))
bench <- fread(file.path(V2, "03_benchmark", "benchmark_per_class.tsv"))
marker_prec <- bench[method == "marker", .(precision = mean(precision)), by = class]
dec_path <- file.path(V2, "05_curation", "cluster_decisions.tsv")
decisions <- if (file.exists(dec_path)) fread(dec_path, colClasses = list(character = "cluster")) else
  data.table(dataset_id = character(), cluster = character(), final_broad = character(),
             final_fine = character(), final_status = character(), rationale = character())

ref_prec <- function(score) {
  bin <- as.character(cut(score, c(0, 0.5, 0.7, 0.9, 1.0001), right = FALSE,
                          labels = c("[0,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1]")))
  cal$precision[match(bin, as.character(cal$score_bin))]
}
mk_prec <- function(cls) {
  p <- marker_prec$precision[match(cls, marker_prec$class)]
  ifelse(is.na(p), 0.8, p)
}
default_fine <- function(b) ont[match(b, broad_class)]$fine_state  # first fine state per class

cluster_tables <- list(); cell_tables <- list()
for (ds in DATASETS) {
  cm <- fread_gz(file.path(V2, "01_marker_v2", ds, "cell_marker_table.tsv.gz"), colClasses = list(character = "cluster"))
  calls <- fread(file.path(V2, "01_marker_v2", ds, "cluster_marker_calls.tsv"), colClasses = list(character = "cluster"))
  ref <- switch(ds, "2022_eLife_Morris_cycling_murine_ovary_atlas" = "isola",
                "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "morris", "combined")
  tr <- fread_gz(file.path(V2, "02_reference_transfer", paste0("transfer_", ds, "__ref_", ref, ".tsv.gz")))
  x <- merge(cm, tr[, .(cell_id, ref_pred_broad, ref_pred_broad_score, ref_pred_fine, ref_pred_label)], by = "cell_id", all.x = TRUE)
  x <- merge(x, calls[, .(cluster, marker_broad, marker_fine, marker_z, marker_margin, cell_purity,
                          cycling_dominated, top_de_genes)], by = "cluster", all.x = TRUE)
  if (!"author_cell_type" %in% names(x)) x[, author_cell_type := NA_character_]
  x[, ref_ok := !is.na(ref_pred_broad_score) & ref_pred_broad_score >= REF_MIN]

  # author labels (Morris, Isola only; Kim's author_cell_type column is our own marker output)
  x[, author_broad := NA_character_][, author_fine := NA_character_][, author_use := FALSE]
  if (ds %in% AUTHOR_DATASETS) {
    m <- amap[dataset_id == ds]
    i <- match(x$author_cell_type, m$author_label)
    x[, `:=`(author_broad = m$broad_class[i], author_fine = m$fine_state[i], author_use = m$use_in_reference[i] %in% TRUE)]
  }

  # cluster-level reference majority
  cl <- x[, {
    rb <- ref_pred_broad[ref_ok]
    tb <- if (length(rb)) sort(table(rb), decreasing = TRUE) else integer()
    .(n_cells = .N,
      ref_major = if (length(tb)) names(tb)[1] else NA_character_,
      ref_major_frac = if (length(tb)) as.numeric(tb[1]) / .N else 0,
      ref_mean_score = mean(ref_pred_broad_score, na.rm = TRUE),
      author_major = if (any(author_use)) names(sort(table(author_broad[author_use]), decreasing = TRUE))[1] else NA_character_,
      author_frac = mean(author_use))
  }, by = cluster]
  cl <- merge(cl, calls[, .(cluster, marker_broad, marker_fine, marker_z, marker_margin, cell_purity,
                            cycling_dominated, top_de_genes)], by = "cluster")
  cl[, dataset_id := ds]
  cl[, auto_case := fcase(
    marker_broad == "Erythroid", "marker_erythroid",
    marker_broad != "Unknown" & marker_broad == ref_major & ref_major_frac >= AGREE_FRAC, "agree",
    marker_broad == "Unknown" & ref_major_frac >= REF_ONLY_FRAC, "reference_only",
    marker_broad != "Unknown" & ref_major_frac < AGREE_FRAC & (is.na(ref_major) | ref_major == marker_broad), "marker_only",
    default = "conflict")]
  cl[, auto_label := fcase(auto_case %in% c("agree", "marker_only", "marker_erythroid"), marker_broad,
                           auto_case == "reference_only", ref_major, default = "Unknown")]
  cl <- merge(cl, decisions[dataset_id == ds, .(cluster, final_broad, final_fine, final_status, rationale)],
              by = "cluster", all.x = TRUE)
  cl[, cluster_label := fifelse(!is.na(final_broad), final_broad, auto_label)]
  cl[, cluster_case := fifelse(!is.na(final_broad), "curated", auto_case)]
  cluster_tables[[ds]] <- cl

  x <- merge(x, cl[, .(cluster, cluster_label, cluster_case, final_fine, final_status, ref_major)], by = "cluster")

  # cell-level decision ------------------------------------------------------
  x[, `:=`(label = "Unknown", status = "ambiguous", conf = 0, source = "unresolved")]
  # author path
  a <- x$author_use
  x[a, label := author_broad]
  x[a, agree_any := (ref_ok & ref_pred_broad == author_broad) | (marker_broad == author_broad)]
  x[a & agree_any == TRUE, `:=`(status = "confirmed", conf = 1, source = "author_metadata")]
  x[a & agree_any == FALSE, `:=`(status = "probable", conf = 0.8, source = "author_metadata_unsupported")]
  # evidence path
  e <- !a & x$cluster_label != "Unknown"
  x[e, label := cluster_label]
  veto <- e & x$ref_pred_broad_score >= REF_STRONG & !is.na(x$ref_pred_broad) &
    x$ref_pred_broad != x$cluster_label & x$cluster_label != "Erythroid" &
    # curated clusters: the overruled reference class (ref_major) does not veto its own cells
    !(x$cluster_case == "curated" & x$ref_pred_broad %in% x$ref_major)
  x[veto, `:=`(label = "Unknown", status = "ambiguous", conf = 0, source = "cell_conflicts_cluster")]
  support <- e & !veto & x$ref_ok & x$ref_pred_broad == x$cluster_label
  x[support & cluster_case %in% c("agree"), `:=`(status = "confirmed",
      conf = pmin(0.99, 1 - (1 - ref_prec(ref_pred_broad_score)) * (1 - mk_prec(cluster_label))), source = "reference+marker")]
  x[support & cluster_case %in% c("reference_only", "marker_only"), `:=`(status = "probable",
      conf = ref_prec(ref_pred_broad_score), source = cluster_case)]
  x[support & cluster_case == "curated", `:=`(status = final_status, conf = fifelse(final_status == "confirmed", 0.95, 0.8), source = "curated")]
  nos <- e & !veto & !support
  x[nos & cluster_case %in% c("agree", "marker_only", "marker_erythroid"), `:=`(status = "probable", conf = mk_prec(cluster_label), source = "marker_cluster")]
  x[nos & cluster_case == "curated", `:=`(status = final_status, conf = 0.75, source = "curated")]
  x[nos & cluster_case == "reference_only", `:=`(label = "Unknown", status = "ambiguous", conf = 0, source = "reference_low_score")]
  x[label == "Unknown", `:=`(status = "ambiguous", conf = 0)]

  # fine state ---------------------------------------------------------------
  x[, fine := NA_character_]
  x[a, fine := author_fine]
  # top fine label must belong to the same broad class (broad call is a sum over fine labels)
  x[is.na(fine) & label == ref_pred_broad & ref_ok & sub("\\|.*$", "", ref_pred_label) == label, fine := ref_pred_fine]
  # v2.2: reference alone over-calls the luteal state (benchmark precision 0.22-0.53);
  # luteal requires marker support (cluster luteal signature or the cell's top signature).
  x[, lut_marker := marker_fine == "luteal" | cell_top_signature == "luteal"]
  x[!a & label == "Granulosa" & fine %in% "luteal" & !(lut_marker %in% TRUE), fine := NA_character_]
  x[!a & label == "Granulosa" & is.na(fine) & marker_fine %in% "luteal", fine := "follicular"]
  x[is.na(fine) & !is.na(final_fine) & label == cluster_label, fine := final_fine]
  x[is.na(fine) & label == marker_broad, fine := marker_fine]
  x[is.na(fine), fine := default_fine(label)]
  x[label == "Unknown", fine := "unknown"]

  cell_tables[[ds]] <- x[, .(dataset_id = ds, cell_id, sample_id, standard_cluster = cluster,
    annotation_label_broad = label, annotation_label_fine = fine, annotation_state = fine,
    annotation_source = source, annotation_confidence = round(conf, 4), curation_status = status,
    author_label = author_cell_type, author_label_broad = author_broad,
    marker_cluster_label = marker_broad, marker_cluster_fine = marker_fine,
    ref_pred_broad, ref_pred_fine, ref_pred_broad_score = round(ref_pred_broad_score, 4),
    cluster_case, phase,
    evidence_notes = paste0("cluster_case=", cluster_case, ";marker=", marker_broad, ";ref=", ref_pred_broad))]
  message(ds, ": ", nrow(x), " cells")
}
cons <- rbindlist(cell_tables, fill = TRUE)
fwrite(cons, file.path(out_dir, "harmonized_annotation_consensus_v2.tsv.gz"), sep = "\t")
clus <- rbindlist(cluster_tables, fill = TRUE)
setcolorder(clus, c("dataset_id", "cluster", "n_cells", "cluster_case", "cluster_label", "marker_broad", "marker_fine",
                    "ref_major", "ref_major_frac", "author_major", "author_frac"))
fwrite(clus, file.path(out_dir, "cluster_consensus_v2.tsv"), sep = "\t")

summ <- cons[, .N, by = .(dataset_id, annotation_source, curation_status)][order(dataset_id, -N)]
fwrite(summ, file.path(out_dir, "consensus_source_summary.tsv"), sep = "\t")
clock_ok <- cons[curation_status %in% c("confirmed", "probable") & annotation_confidence >= 0.7 &
                   !annotation_label_broad %in% c("Unknown", "Erythroid", "Oocyte_Germline")]
fwrite(clock_ok[, .N, by = .(dataset_id, annotation_label_broad)][order(dataset_id, annotation_label_broad)],
       file.path(out_dir, "clock_eligible_cells_by_class.tsv"), sep = "\t")
print(clus[, .N, by = .(cluster_case)])
print(clus[cluster_case == "conflict", .(dataset_id, cluster, n_cells, marker_broad, ref_major, ref_major_frac)])
