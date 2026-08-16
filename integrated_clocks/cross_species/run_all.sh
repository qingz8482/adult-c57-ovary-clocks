#!/usr/bin/env bash
set -euo pipefail

MODULE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${MODULE_DIR}/../../.." && pwd)"
R_BIN="${R_BIN:-$(command -v Rscript)}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export XDG_CACHE_HOME="${MODULE_DIR}/logs/fontcache"
mkdir -p "${XDG_CACHE_HOME}"
cd "${REPO_DIR}"

for script in \
  00_inventory_and_freeze_inputs.R \
  01_extract_mouse_models_and_weights.R \
  02_build_mouse_human_orthology.R \
  03_prepare_gtex_ovary.R \
  04_prepare_hpa_ovary.R \
  05_stageA_age_direction_conservation.R \
  06_stageB_frozen_zero_shot_projection.R \
  07_random_signature_and_permutation_tests.R \
  08_build_publication_panels.R \
  09_render_report.R
do
  "${R_BIN}" --vanilla "${MODULE_DIR}/scripts/${script}" 2>&1 | tee "${MODULE_DIR}/logs/${script%.R}.log"
done
"${R_BIN}" --vanilla "${MODULE_DIR}/validate_outputs.R" 2>&1 | tee "${MODULE_DIR}/logs/validate_outputs.log"
