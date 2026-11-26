#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

required=(README.md LICENSE CITATION.cff environment.yml renv.lock integrated_models/HVG-ElasticNet-exclude-single-age/model.rds integrated_models/VersionB-PCA-glmnet/model.rds integrated_clocks metacell_clocks metadata/manuscript_tables metadata/pipeline_inputs docs/manuscript_code_map.tsv tests)
for path in "${required[@]}"; do test -e "$path" || { echo "missing required path: $path" >&2; exit 1; }; done

if grep -RInE '/mnt/beegfs|/home/qingz8482' --exclude-dir=.git --exclude=validate_repository.sh .; then
  echo "personal or cluster-specific path found" >&2
  exit 1
fi

if find . -type f -size +50M -print | grep -q .; then
  echo "file larger than 50 MiB found" >&2
  exit 1
fi

if find . -type d \( -name archive -o -name results -o -name outputs -o -name __pycache__ -o -name .venv \) -print | grep -q .; then
  echo "generated or archived directory found" >&2
  exit 1
fi

if grep -RInE '(BEGIN (RSA|OPENSSH|EC) PRIVATE KEY|ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})' --exclude-dir=.git .; then
  echo "possible credential found" >&2
  exit 1
fi

while IFS= read -r -d '' script; do python3 -c 'import ast, pathlib, sys; ast.parse(pathlib.Path(sys.argv[1]).read_text())' "$script"; done < <(find integrated_clocks metacell_clocks integrated_models -name '*.py' -print0)

rscript_bin="${RSCRIPT:-$(command -v Rscript || true)}"
if [[ -n "$rscript_bin" ]]; then
  "$rscript_bin" -e 'files <- list.files(c("integrated_clocks", "metacell_clocks", "integrated_models", "tests"), pattern="[.]R$", recursive=TRUE, full.names=TRUE); invisible(lapply(files, parse))'
  "$rscript_bin" tests/validate_metadata.R "$repo_dir"
  while IFS=$'\t' read -r _ _ analysis; do
    for p in $(printf '%s' "$analysis" | tr ';' '\n' | sed -E 's/^ +//'); do
      test -e "$p" || { echo "manuscript_code_map.tsv path missing: $p" >&2; exit 1; }
    done
  done < <(tail -n +2 docs/manuscript_code_map.tsv)
else
  echo "Rscript unavailable: R parse and metadata checks skipped" >&2
fi

echo "repository validation passed"
