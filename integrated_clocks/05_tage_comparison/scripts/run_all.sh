#!/usr/bin/env bash
set -euo pipefail
TASK_SCRIPTS=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
TASK_OUT=$(dirname -- "$TASK_SCRIPTS")
TASK_OVARY=$(dirname -- "$(dirname -- "$(dirname -- "$TASK_OUT")")")
TASK_R=${RSCRIPT:-Rscript}
TASK_PY="${TAGE_PYTHON:-$TASK_OVARY/bulk-pseudobulk_clock/results/external_projection/tAge/environment/python_env_locked/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
mkdir -p "$TASK_OUT/logs"
"$TASK_R" "$TASK_SCRIPTS/01_compute_comparisons.R" > "$TASK_OUT/logs/01_compute.log" 2>&1
"$TASK_PY" "$TASK_SCRIPTS/02_extract_tage_features.py" > "$TASK_OUT/logs/02_features.log" 2>&1
"$TASK_R" "$TASK_SCRIPTS/03_figures.R" > "$TASK_OUT/logs/03_figures.log" 2>&1
"$TASK_PY" "$TASK_SCRIPTS/04_validate.py" > "$TASK_OUT/logs/04_validation.log" 2>&1
"$TASK_R" "$TASK_SCRIPTS/05_pathways.R" > "$TASK_OUT/logs/05_pathways.log" 2>&1
"$TASK_PY" "$TASK_SCRIPTS/06_validate_pathways.py" > "$TASK_OUT/logs/06_pathway_validation.log" 2>&1
"$TASK_R" "$TASK_SCRIPTS/07_model_resolved_pathways.R" > "$TASK_OUT/logs/07_model_resolved.log" 2>&1
"$TASK_R" "$TASK_SCRIPTS/08_validate_permutation_engine.R" > "$TASK_OUT/logs/08_permutation_validation.log" 2>&1
"$TASK_PY" "$TASK_SCRIPTS/09_validate_model_resolved.py" > "$TASK_OUT/logs/09_model_resolved_validation.log" 2>&1
"$TASK_PY" "$TASK_SCRIPTS/10_compare_rebutal3_pca.py" > "$TASK_OUT/logs/10_compare_rebutal3_pca.log" 2>&1
