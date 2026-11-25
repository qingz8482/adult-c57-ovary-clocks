#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
rscript="${RSCRIPT:-Rscript}"
python3 "$task_dir/../lockedmodel/verify_checksums.py" > "$task_dir/logs/release_checksums.log" 2>&1
python3 "$task_dir/scripts/07_build_reference_scripts.py"
for step in 01_prepare_project 02_statistics 03_nulls 04_figures 05_validate; do
 "$rscript" "$task_dir/scripts/$step.R" > "$task_dir/logs/$step.log" 2>&1
done
python3 "$task_dir/scripts/08_validate_reference_figures.py"
python3 "$task_dir/scripts/06_report.py"
"$rscript" "$task_dir/scripts/09_compare_rebutal3_pca.R" > "$task_dir/logs/09_compare_rebutal3_pca.log" 2>&1
