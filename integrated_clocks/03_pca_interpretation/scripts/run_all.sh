#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
rscript="${RSCRIPT:-Rscript}"
"$rscript" "$task_dir/scripts/01_pca_diagnostics.R" > "$task_dir/logs/01_pca_diagnostics.log" 2>&1
"$rscript" "$task_dir/scripts/02_deployment_gsea.R" > "$task_dir/logs/02_deployment_gsea.log" 2>&1
"$rscript" "$task_dir/scripts/03_validate.R" > "$task_dir/logs/03_validate.log" 2>&1
python3 "$task_dir/scripts/07_validate_heatmap.py" > "$task_dir/logs/07_validate_heatmap.log" 2>&1
"$rscript" "$task_dir/scripts/05_visual_qa.R" > "$task_dir/logs/05_visual_qa.log" 2>&1
bash "$task_dir/supple2_PCA/scripts/run_all.sh"
python3 "$task_dir/scripts/04_report.py"
