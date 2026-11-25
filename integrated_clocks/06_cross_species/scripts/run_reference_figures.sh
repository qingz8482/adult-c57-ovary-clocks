#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 "$task_dir/scripts/07_build_reference_scripts.py"
"${RSCRIPT:-Rscript}" "$task_dir/scripts/04_figures.R" > "$task_dir/logs/reference_figures.log" 2>&1
"${RSCRIPT:-Rscript}" "$task_dir/scripts/05_validate.R" > "$task_dir/logs/05_validate.log" 2>&1
python3 "$task_dir/scripts/08_validate_reference_figures.py"
python3 "$task_dir/scripts/06_report.py"
