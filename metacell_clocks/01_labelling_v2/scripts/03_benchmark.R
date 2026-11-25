# Step B3: benchmark reference transfer and marker v2 against author labels (Morris <-> Isola).
source(file.path(Sys.getenv("OVARY_CLOCK_ROOT"), "adult_c57_ovary_clocks/rebutal/rebutal4/metacell_review/v2/scripts/helpers_v2.R"))
out_dir <- file.path(V2, "03_benchmark")
dir.create(out_dir, showWarnings = FALSE)
amap <- author_map()
MAIN <- c("Granulosa", "Theca_Stroma", "Endothelial", "Immune_Myeloid")
PASS <- 0.85

refs <- c("2022_eLife_Morris_cycling_murine_ovary_atlas" = "isola",
          "2024_NatureAging_Isola_aging_mouse_ovary_atlas" = "morris")

per_class <- list(); curves <- list(); conf <- list(); luteal <- list(); calib <- list()
for (ds in names(refs)) {
  cm <- fread_gz(file.path(V2, "01_marker_v2", ds, "cell_marker_table.tsv.gz"))
  calls <- fread(file.path(V2, "01_marker_v2", ds, "cluster_marker_calls.tsv"), colClasses = list(character = "cluster"))
  tr <- fread_gz(file.path(V2, "02_reference_transfer", paste0("transfer_", ds, "__ref_", refs[[ds]], ".tsv.gz")))
  cm[, cluster := as.character(cluster)]
  m <- amap[dataset_id == ds]
  cm[, author_broad := m$broad_class[match(author_cell_type, m$author_label)]]
  cm[, author_fine := m$fine_state[match(author_cell_type, m$author_label)]]
  cm <- merge(cm, calls[, .(cluster, marker_broad, marker_fine)], by = "cluster", all.x = TRUE)
  cm <- merge(cm, tr[, .(cell_id, ref_pred_broad, ref_pred_broad_score, ref_pred_fine, ref_pred_score)], by = "cell_id")
  cm[, author_use := m$use_in_reference[match(author_cell_type, m$author_label)]]
  cm <- cm[!is.na(author_broad) & author_use %in% TRUE]  # ambiguous author clusters excluded from truth
  ref_classes <- unique(sub("\\|.*$", "", readRDS(file.path(V2, "02_reference_transfer", paste0("ref_", refs[[ds]], ".rds")))$ref_label))

  for (method in c("reference", "marker")) {
    pred <- if (method == "reference") cm$ref_pred_broad else cm$marker_broad
    truth <- cm$author_broad
    conf[[paste(ds, method)]] <- data.table(truth = truth, pred = pred)[, .N, by = .(truth, pred)][
      , `:=`(query = ds, method = method)]
    for (k in sort(unique(truth))) {
      tp <- sum(pred == k & truth == k); fp <- sum(pred == k & truth != k); fn <- sum(pred != k & truth == k)
      per_class[[length(per_class) + 1]] <- data.table(
        query = ds, reference = if (method == "reference") refs[[ds]] else "marker_v2", method = method, class = k,
        n_truth = sum(truth == k), precision = tp / max(tp + fp, 1), recall = tp / max(tp + fn, 1),
        class_in_reference = if (method == "reference") k %in% ref_classes else TRUE)
    }
  }
  # Calibration bins (truth classes present in the reference only)
  cb <- cm[author_broad %in% ref_classes]
  cb[, score_bin := cut(ref_pred_broad_score, c(0, 0.5, 0.7, 0.9, 1.0001), right = FALSE,
                        labels = c("[0,0.5)", "[0.5,0.7)", "[0.7,0.9)", "[0.9,1]"))]
  calib[[ds]] <- cb[, .(n = .N, precision = mean(ref_pred_broad == author_broad)), by = score_bin][, query := ds]
  # Threshold curve for reference transfer (broad score)
  for (t in seq(0.3, 0.95, by = 0.05)) {
    keep <- cm$ref_pred_broad_score >= t
    curves[[length(curves) + 1]] <- data.table(query = ds, threshold = t, coverage = mean(keep),
                                               precision = mean(cm$ref_pred_broad[keep] == cm$author_broad[keep]))
  }
  # Luteal state within Granulosa
  g <- cm[author_broad == "Granulosa"]
  luteal[[ds]] <- data.table(
    query = ds, n_author_luteal = sum(g$author_fine == "luteal"),
    ref_luteal_recall = mean(g[author_fine == "luteal", ref_pred_broad == "Granulosa" & ref_pred_fine == "luteal"]),
    ref_luteal_precision = { p <- cm$ref_pred_broad == "Granulosa" & cm$ref_pred_fine == "luteal"
                             mean(cm$author_fine[p] == "luteal") },
    marker_luteal_recall = mean(g[author_fine == "luteal", marker_broad == "Granulosa" & marker_fine == "luteal"]),
    marker_luteal_precision = { p <- cm$marker_broad == "Granulosa" & cm$marker_fine == "luteal"
                                mean(cm$author_fine[p] == "luteal") },
    both_luteal_recall = { b <- cm$ref_pred_broad == "Granulosa" & cm$ref_pred_fine == "luteal" &
                             (cm$marker_fine == "luteal" | cm$cell_top_signature == "luteal")
                           mean(b[cm$author_fine == "luteal"]) },
    both_luteal_precision = { b <- cm$ref_pred_broad == "Granulosa" & cm$ref_pred_fine == "luteal" &
                                (cm$marker_fine == "luteal" | cm$cell_top_signature == "luteal")
                              mean(cm$author_fine[b] == "luteal") })
}
pc <- rbindlist(per_class)
pc[, f1 := ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)]
pc[, gate_pass := class %in% MAIN & precision >= PASS & recall >= PASS]
fwrite(pc, file.path(out_dir, "benchmark_per_class.tsv"), sep = "\t")
fwrite(rbindlist(curves), file.path(out_dir, "benchmark_reference_threshold_curve.tsv"), sep = "\t")
fwrite(rbindlist(conf), file.path(out_dir, "benchmark_confusion.tsv"), sep = "\t")
fwrite(rbindlist(luteal), file.path(out_dir, "benchmark_luteal_state.tsv"), sep = "\t")
cal <- rbindlist(calib)
cal_pooled <- cal[, .(n = sum(n), precision = sum(n * precision) / sum(n)), by = score_bin][order(score_bin)]
fwrite(cal, file.path(out_dir, "calibration_reference_by_query.tsv"), sep = "\t")
fwrite(cal_pooled, file.path(out_dir, "calibration_reference_pooled.tsv"), sep = "\t")
print(cal_pooled)

gate <- pc[class %in% MAIN & class_in_reference, .(n_main = .N, n_pass = sum(gate_pass)), by = .(query, method)]
fwrite(gate, file.path(out_dir, "gate_B_summary.tsv"), sep = "\t")
print(pc[order(query, method, class)])
print(gate)
print(rbindlist(luteal))
