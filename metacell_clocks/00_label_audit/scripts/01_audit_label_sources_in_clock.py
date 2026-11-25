#!/usr/bin/env python3
"""Audit which annotation sources and flagged clusters reached the
20260707 metacell clocks. Standard library only.

Outputs (tables/):
  01_label_source_by_dataset_clock.tsv   clock-candidate metacells/cells by dataset x label source
  02_cluster_clock_contribution.tsv      per annotation cluster: cells, clock-eligible cells, source
  04_flagged_cluster_clock_impact.tsv    review flags joined to clock-eligible cell counts
  05_flag_impact_by_class.tsv            clock-eligible cells affected, by current class and issue
"""
import csv
import os
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

ROOT = os.path.join(os.environ["OVARY_CLOCK_ROOT"], "adult_c57_ovary_clocks")
RUN = os.path.join(ROOT, "7 scRNAseq_clock_results/clocks/metacell_level_scRNA/"
                   "scRNA_metacell_elasticnet_metacell_cv_20260707")
CONSENSUS = os.path.join(ROOT, "data_input/standardized_seurat/00_metadata/harmonized_annotation_consensus.tsv")
SAMPLE_ROLES = os.path.join(RUN, "step_01_metadata/scRNA_sample_role_metadata.tsv")
METACELLS = os.path.join(RUN, "step_02_preprocessing_filtering/filtered_metacell_metadata.tsv")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tables")
FLAGS = os.path.join(OUT, "03_cluster_review_flags.tsv")

# Same filters as run_scRNA_metacell_elasticnet_rerun.R
STATUSES = {"confirmed", "probable"}
MIN_CONF = 0.70
EXCLUDED = {"Unknown", "Erythroid"}
NOT_MODELLED = {"Oocyte_Germline"}  # filtered metacells exist but class failed CV eligibility


def read_tsv(path):
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield {k: (v.strip('"') if v is not None else v) for k, v in row.items()}


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


# 1. Metacell-level source summary -------------------------------------------
mc = defaultdict(lambda: [0, 0])
total_cells = 0
for r in read_tsv(METACELLS):
    if r["model_role"] != "clock_candidate":
        continue
    k = (r["dataset_id"], r["label_assignment_source"])
    mc[k][0] += 1
    mc[k][1] += int(r["n_cells"])
    total_cells += int(r["n_cells"])
rows = [[d, s, n, c, f"{100 * c / total_cells:.1f}"] for (d, s), (n, c) in sorted(mc.items())]
write_tsv(os.path.join(OUT, "01_label_source_by_dataset_clock.tsv"),
          ["dataset_id", "label_assignment_source", "n_clock_metacells", "n_cells", "pct_of_clock_cells"], rows)

# 2. Cluster-level contribution (consensus datasets only; PRJNA863443 is labelled outside the consensus) --
clock_samples = set()
for r in read_tsv(SAMPLE_ROLES):
    if r["model_role"] == "clock_candidate" and r["age_weeks"] not in ("", "NA") \
            and r["reproductive_age_bin"] != "out_of_scope":
        clock_samples.add((r["dataset_id"], r["sample_id"]))

clu = defaultdict(lambda: {"n": 0, "elig": 0, "label": defaultdict(int), "src": defaultdict(int)})
for r in read_tsv(CONSENSUS):
    k = (r["dataset_id"], r["standard_cluster"])
    c = clu[k]
    c["n"] += 1
    c["label"][r["annotation_label_broad"]] += 1
    c["src"][r["annotation_source"]] += 1
    if ((r["dataset_id"], r["sample_id"]) in clock_samples
            and r["curation_status"] in STATUSES
            and fnum(r["annotation_confidence"]) >= MIN_CONF
            and r["annotation_label_broad"] not in EXCLUDED
            and r["annotation_label_broad"] not in NOT_MODELLED):
        c["elig"] += 1


def top(d):
    return max(d.items(), key=lambda kv: kv[1])[0] if d else "NA"


rows = []
for (d, cl), c in sorted(clu.items()):
    rows.append([d, cl, c["n"], c["elig"], top(c["label"]), top(c["src"]),
                 ";".join(f"{k}={v}" for k, v in sorted(c["src"].items()))])
write_tsv(os.path.join(OUT, "02_cluster_clock_contribution.tsv"),
          ["dataset_id", "standard_cluster", "n_cells", "n_clock_eligible_cells",
           "main_label", "main_source", "source_breakdown"], rows)

# 4/5. Flag impact -----------------------------------------------------------
# Excluded-class clusters (Erythroid/Oocyte) contribute 0 eligible cells by construction;
# their impact is cells LOST, so report cells in clock samples regardless of label.
lost = defaultdict(int)
for r in read_tsv(CONSENSUS):
    if (r["dataset_id"], r["sample_id"]) in clock_samples:
        lost[(r["dataset_id"], r["standard_cluster"])] += 1

flag_rows, by_class = [], defaultdict(lambda: [0, 0])
for f in read_tsv(FLAGS):
    k = (f["dataset_id"], f["standard_cluster"])
    c = clu.get(k, {"n": 0, "elig": 0})
    in_clock_samples = lost.get(k, 0)
    flag_rows.append([f["dataset_id"], f["standard_cluster"], f["current_label"], f["proposed_label"],
                      f["issue_category"], f["severity"], c["n"], in_clock_samples, c["elig"]])
    key = (f["current_label"], f["issue_category"])
    by_class[key][0] += c["elig"]
    by_class[key][1] += in_clock_samples
write_tsv(os.path.join(OUT, "04_flagged_cluster_clock_impact.tsv"),
          ["dataset_id", "standard_cluster", "current_label", "proposed_label", "issue_category", "severity",
           "n_cells_total", "n_cells_in_clock_samples", "n_clock_eligible_cells_current_label"], flag_rows)
write_tsv(os.path.join(OUT, "05_flag_impact_by_class.tsv"),
          ["current_label", "issue_category", "n_clock_eligible_cells", "n_cells_in_clock_samples"],
          [[a, b, v[0], v[1]] for (a, b), v in sorted(by_class.items())])
print("done")
