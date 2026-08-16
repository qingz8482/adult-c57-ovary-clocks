#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"
rscript_bin="${RSCRIPT:-$(command -v Rscript || true)}"
test -n "$rscript_bin" || { echo "Rscript is required" >&2; exit 1; }
"$rscript_bin" tests/validate_metadata.R "$repo_dir"
"$rscript_bin" tests/smoke_models.R
echo "all smoke tests passed"
