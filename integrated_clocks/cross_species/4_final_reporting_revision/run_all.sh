#!/usr/bin/env bash
set -euo pipefail

REVISION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_BIN="${R_BIN:-$(command -v Rscript)}"
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-/tmp/gtex_final_reporting_fontcache}"

cd "$REVISION_DIR"
"$R_BIN" scripts/01_centroid_similarity.R 2>&1 | tee logs/01_centroid_similarity.log
"$R_BIN" scripts/02_hardy_features_patterns.R 2>&1 | tee logs/02_hardy_features_patterns.log
"$R_BIN" scripts/03_sparse_oldest_sensitivity.R 2>&1 | tee logs/03_sparse_oldest_sensitivity.log
"$R_BIN" scripts/04_nested_cv_validation.R 2>&1 | tee logs/04_nested_cv_validation.log
"$R_BIN" scripts/05_reports_figures.R 2>&1 | tee logs/05_reports_figures.log
"$R_BIN" scripts/06_visual_qa_repairs.R 2>&1 | tee logs/06_visual_qa_repairs.log
"$R_BIN" validate_final_revision.R 2>&1 | tee logs/validate_final_revision.log
