# Environment

The reference environment uses R 4.4 and Python 3.12. `environment.yml` provides the complete user-facing environment, while `renv.lock` records direct R dependencies for R-only restoration.

Core R packages include `data.table`, `Matrix`, `dplyr`, `tidyr`, `tibble`, `stringr`, `ggplot2`, `glmnet`, `msigdbr`, and `Seurat`. Cross-species analyses additionally use `edgeR`, `limma`, and `matrixStats`; plot assembly uses `patchwork`, `ggrepel`, `forcats`, `scales`, and `ragg`.

External requirements:

- sufficient memory for the prepared matrices and Seurat objects;
- access to the public repositories listed in `metadata/dataset_accessions.tsv`;
- acceptance of MSigDB terms for pathway retrieval;
- optional SLURM for cluster-scale execution.

No script requires an Abraxas-specific path. Use `OVARY_CLOCK_DATA_ROOT`, `OVARY_CLOCK_OUTPUT_ROOT`, and, for metacells, `OVARY_METACELL_PREPROCESS_ROOT`.
