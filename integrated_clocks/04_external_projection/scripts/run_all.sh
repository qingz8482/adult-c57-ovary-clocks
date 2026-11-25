#!/usr/bin/env bash
set -euo pipefail
task_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
rscript="${RSCRIPT:-Rscript}"
"$rscript" "$task_dir/scripts/01_project_locked.R" > "$task_dir/logs/01_project_locked.log" 2>&1
"$rscript" "$task_dir/scripts/02_statistics_plots.R" > "$task_dir/logs/02_statistics_plots.log" 2>&1
"$rscript" "$task_dir/scripts/03_validate.R" > "$task_dir/logs/03_validate.log" 2>&1
python3 "$task_dir/scripts/12_record_numeric_baseline.py"
"$rscript" "$task_dir/scripts/13_compare_rebutal3_pca.R" > "$task_dir/logs/13_compare_rebutal3_pca.log" 2>&1
"$rscript" "$task_dir/scripts/07_complete_reference_outputs.R" > "$task_dir/logs/07_complete_reference_outputs.log" 2>&1
python3 "$task_dir/scripts/06_validate_plot_revision.py"
"$rscript" "$task_dir/scripts/09_validate_derived_analyses.R" > "$task_dir/logs/09_validate_derived_analyses.log" 2>&1
python3 "$task_dir/scripts/08_audit_reference_completeness.py"
"$rscript" "$task_dir/scripts/11_gse235021_mixed_modality_delta.R" > "$task_dir/logs/11_gse235021_mixed_modality_delta.log" 2>&1
python3 "$task_dir/scripts/10_review_application_files.py"
python3 "$task_dir/scripts/04_report.py"
