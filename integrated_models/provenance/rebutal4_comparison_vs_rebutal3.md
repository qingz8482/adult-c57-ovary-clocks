# Exclusion of single-age training studies vs locked rebutal3 models

Completed 2026-10-12. All phases ran: A, B, C, final and A_repeated, with 30 repeats × 4 folds. The comparison uses count-tuned held-out predictions, paired on identical held-out profiles. "Locked" = rebutal3 `count_tuning`; "Excl" = this run. MAE is in weeks.

Reproduction check: rebutal4 main A reproduces rebutal3 `A_exclude_single_age` predictions exactly (max |difference| = 0) for all 3 recipes × 2 modes.

## Held-out MAE

| Clock | A pooled | A study-mean | B NIHMS | C Winkler | 30 repeats pooled, median [IQR] | Excl better (of 30) |
|---|---|---|---|---|---|---|
| HVG–EN locked | 16.10 | 16.17 | 23.18 | 17.65 | 15.96 [15.3–17.1] | |
| HVG–EN excl | **9.90** | **10.07** | **9.33** | 16.91 | **10.98 [10.6–11.3]** | **30/30** |
| Version B–EN locked | 15.00 | 16.79 | 23.39 | 17.12 | 15.74 [14.7–16.3] | |
| Version B–EN excl | 14.13 | 14.27 | **8.18** | 19.81 (worse) | 12.86 [12.2–13.5] | 30/30 |
| Version B–PCA locked | 15.73 | 15.44 | 21.30 | 17.73 | 15.70 [15.2–16.9] | |
| Version B–PCA excl | 15.23 | 14.91 | 16.09 | 19.33 (worse) | 15.47 [15.3–15.7] | 19/30 |

Repeated splits, median MAE by test subset (locked → excl, excl better in n/30):

| Clock | Single-age test studies | Multi-age test studies | scRNA pseudobulk test | Bulk test |
|---|---|---|---|---|
| HVG–EN | 18.8 → 9.9 (30) | 15.5 → 11.2 (30) | 16.4 → 12.6 (30) | 15.7 → 9.4 (30) |
| Version B–EN | 19.6 → 15.1 (28) | 15.0 → 12.5 (26) | 15.9 → 14.9 (22) | 15.5 → 11.2 (30) |
| Version B–PCA | 16.9 → 20.5 (0) | 15.4 → 14.4 (27) | 17.6 → 18.0 (14) | 14.0 → 13.0 (24) |

Selected studies in the repeats (median per-study MAE, locked → excl, excl better in n/30):

| Clock | Winkler | NIHMS | Morris | CRA003645 | Jiang |
|---|---|---|---|---|---|
| HVG–EN | 17.0 → 17.0 (9) | 23.6 → 9.7 (22) | 18.7 → 7.0 (30) | 17.7 → 7.7 (29) | 16.1 → 10.0 (28) |
| Version B–EN | 15.1 → 19.8 (8) | 25.2 → 13.0 (22) | 18.1 → 13.7 (28) | 12.8 → 8.7 (25) | 15.7 → 8.1 (23) |
| Version B–PCA | 18.4 → 19.3 (12) | 22.9 → 16.3 (29) | 19.5 → 22.2 (0) | 12.2 → 12.5 (15) | 15.6 → 12.1 (29) |

## Against the prespecified rules (README)

- **HVG–EN: passes.** The gain holds in every repeated split (30/30), and for every subset: single-age vs multi-age test studies, scRNA vs bulk. B improves sharply. C does not get worse (17.65 → 16.91; Winkler is unchanged in the repeats). Single-age studies improve even though they are never used for training.
- **Version B–EN: partial.** It improves in the repeats and on B, but Winkler gets worse (C 17.1 → 19.8; worse in 22/30 repeats).
- **Version B–PCA: does not pass.** There is no consistent gain, and single-age test studies get worse in 30/30 splits.

## Caveats

- The exclusion rule was chosen after seeing the rebutal3 sensitivity results. Its improvement is a model-development finding, not an independent validation. Independent data is still needed.
- Repeated splits share the same 14 studies; their spread is not a confidence interval.
- The candidate final models (`*/count_tuning/results/final_model/`) were trained on 110 profiles from 9 studies, with minimum training age 8 wk. They are not locked. Their inner MAEs are not comparable with the locked models' (9 vs 14 inner studies).
- Candidate final settings: HVG–EN N = 3,000, α = 0.5, λ = 0.179, 79 nonzero genes. Version B–EN N = 2,000, α = 0.1, λ = 1.26, 155 nonzero genes. Version B–PCA N = 500, 3 PCs, α = 0.1.
