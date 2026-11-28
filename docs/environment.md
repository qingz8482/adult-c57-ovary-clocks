# Environment

The locked fits used R 4.4.3 with glmnet 5.0, data.table 1.17.8, Matrix 1.7-5 and jsonlite 2.0.0. The complete session records are in `integrated_clocks/final_models/*/provenance/sessionInfo.txt`. `environment.yml` (conda) and `renv.lock` (R only) give a matching user environment.

Additional packages by component:

| Component | Packages |
|---|---|
| Model development | glmnet, data.table, jsonlite, Matrix |
| PCA/EN GSEA, metacell GSEA | msigdbr (MSigDB 2026.1.Mm) |
| EN non-zero-gene ORA | clusterProfiler, org.Mm.eg.db (run in a separate R environment; `07_en_ora`) |
| Cross-species | limma, edgeR, matrixStats, MASS, Rcpp |
| Metacell labelling and clocks | Seurat (anchor-based label transfer), Matrix, glmnet |
| tAge comparison | tAge v1.0.1 R package and Biobase, for predictions and control subtraction; a Python environment with the tAge model files (`TAGE_PYTHON`), for feature extraction and validation |
| Plotting and tables | ggplot2, patchwork, ggrepel, scales, ragg, systemfonts, png, dplyr, tidyr, tibble, stringr, readr, readxl |

Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` for bit-identical replay of the locked fits; the shell runners do this. The metacell labelling and CV steps need enough memory to load the standardized Seurat objects of eight datasets (about 350,000 cells). They were run on a SLURM cluster.
