#!/usr/bin/env bash
set -euo pipefail

AUDIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${AUDIT_DIR}/../../../.." && pwd)"
R_BIN="${R_BIN:-$(command -v Rscript)}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export XDG_CACHE_HOME="${AUDIT_DIR}/logs/fontcache"
mkdir -p "${AUDIT_DIR}/logs" "${XDG_CACHE_HOME}"
cd "${REPO_DIR}"

for script in \
  01_metadata_expression_qc.R \
  02_age_de_pairwise_permutation.R \
  03_multivariate_pathway_prediction.R \
  03b_agonal_sensitivity.R \
  04_figures_reports.R
do
  "${R_BIN}" --vanilla "${AUDIT_DIR}/scripts/${script}" 2>&1 | tee "${AUDIT_DIR}/logs/${script%.R}.log"
done
"${R_BIN}" --vanilla "${AUDIT_DIR}/validate_outputs.R" 2>&1 | tee "${AUDIT_DIR}/logs/validate_outputs.log"
