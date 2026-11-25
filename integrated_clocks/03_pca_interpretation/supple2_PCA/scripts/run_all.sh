#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
rscript="${RSCRIPT:-Rscript}"
"$rscript" "$task_dir/scripts/01_compute.R" > "$task_dir/logs/01_compute.log" 2>&1
"$rscript" "$task_dir/scripts/02_gsea_curves.R" > "$task_dir/logs/02_gsea_curves.log" 2>&1
"$rscript" "$task_dir/scripts/03_plot.R" > "$task_dir/logs/03_plot.log" 2>&1
"$rscript" "$task_dir/scripts/04_validate.R" > "$task_dir/logs/04_validate.log" 2>&1
python3 "$task_dir/scripts/06_validate_f4.py" > "$task_dir/logs/06_validate_f4.log" 2>&1
python3 "$task_dir/scripts/05_report.py"
