# Integrated EN and PCA/EN clocks

Each numbered folder is one analysis step. Within a step, scripts keep their original relative layout (`scripts/`, plus sibling folders such as `supple2_PCA/`), so relative paths in the code resolve when the step folder is placed at its workspace location ([`../docs/workflow.md`](../docs/workflow.md)).

| Step | Folder | What it does |
|---|---|---|
| 1 | `01_input_preparation/` | `prepare_inputs.R` builds the 136-profile matrix. Expression is ln(1 + 10⁶ × count / full pre-filter library total), duplicate mapped identifiers are summed, and identifiers are mapped through the frozen dictionary. It also writes the fold assignments (Stage A plus 30 repeats) and input checksums. Nothing is filtered or batch-corrected here. |
| 2 | `02_model_development/` | Study-grouped nested CV and final fits (`core.R`, `run_models.R`, `run_all.py`, `summarize_results.R`). |
| 3 | `03_pca_interpretation/` | Diagnostics for the locked PCA/EN clock. `scripts/` covers variance, scores and per-PC gene weights, plus deployment-model pre-ranked GSEA (Fig. 2). `supple2_PCA/` covers coefficient paths, conditional leave-one-active-PC-out attribution and the contribution heatmaps (Supplementary Fig. 2). `backprojection/` compares back-projected gene weights between the pre- and post-revision PCA/EN fits. |
| 4 | `04_external_projection/` | Frozen projection of the intervention cohorts. Also matched-control Welch tests with BH correction, phenotype benchmarking, PC-contribution heatmaps, GSE235021 mixed-modality views and the minimum detectable shift (Fig. 3, Fig. 5, Supplementary Fig. 4). |
| 5 | `05_tage_comparison/` | Comparison with the tAge scaled-difference and YuGene-difference clocks: contrast effects, gene-weight concordance and model-resolved pathway permutation tests (Fig. 3, Supplementary Fig. 3). |
| 6 | `06_cross_species/` | Stage A: conservation of gene-level age direction (mouse vs GTEx). Stage B: frozen zero-shot projection to GTEx and GSE202601, with permuted-weight, matched-random-signature and Freedman–Lane nulls (Supplementary Fig. 7). |
| 7 | `07_en_ora/` | One-sided hypergeometric over-representation of the 79 non-zero EN genes, split by coefficient sign, against the 2,952 Entrez-mapped HVG candidate genes (GO, KEGG, Reactome, Hallmark; BH correction). |

## Model development (step 2)

`run_all.py` drives `run_models.R` over the phases `A`, `B`, `C`, `final` and `A_repeated` with checkpointing:

```bash
python3 scripts/run_all.py --workers=8 --phases='A;B,C;final;A_repeated'
```

| Config | Clock | Difference |
|---|---|---|
| `config/production.json` | EN (locked) | `"training_exclusion": "single_age_studies"`: studies with a single distinct age (Morris, CRA017492, GSE252791, Zhao, GSE194196; 26 profiles) are dropped from each training set before filtering, selection, inner CV or fitting. Held-out sets are unchanged. |
| `config/production_rebutal3_reference.json` | PCA/EN (locked) | No training exclusion. All other settings are identical. |

`run_models.R` applies the exclusion only when the config requests it, so the same code reproduces both locked fits. Shared settings:

- **Gene filter:** genes measured in every training profile, with CPM ≥ 1 in ≥ 50% of training profiles and positive variance.
- **Gene count:** N ∈ {500, 1,000, 2,000, 3,000, 5,000}.
- **Penalty:** α ∈ {0.1, 0.25, 0.5, 0.75, 0.9, 1}, over 100 fixed λ values from 10³ to 10⁻⁴.
- **PCA/EN only:** 1–10 PCs.
- **Weights:** study-balanced.
- **Selection:** candidates are scored by inner leave-one-study-out MAE. λ follows the one-SE rule; ties go to the smaller N, then the larger α.

`tests/test_pipeline.R` checks the CV, selection and leakage guards. It runs on the prepared 136-profile matrix.

## Locked models

The final fits are copied unchanged to [`../locked_models/hvg_exclusion_vbpca_20261012`](../locked_models/README.md). Steps 3–7 check that they read this release before reporting.

## Third-party components

The tAge models (v1.0.1, Zenodo record 18763485), MSigDB 2026.1.Mm and KEGG gene sets are not redistributed. `05_tage_comparison/scripts/run_all.sh` takes the tAge Python environment from `TAGE_PYTHON`.
