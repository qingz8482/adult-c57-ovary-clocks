#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
rscript_bin="${RSCRIPT:-Rscript}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 "$repo_dir/integrated_models/verify_checksums.py"
"$rscript_bin" "$repo_dir/tests/validate_metadata.R" "$repo_dir"
"$rscript_bin" "$repo_dir/tests/test_integrated_models.R"
echo "smoke tests passed"
