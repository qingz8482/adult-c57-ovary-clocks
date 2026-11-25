# Rebutal4: excluding single-age studies from all training sets

Started 2026-10-12 (Asia/Singapore). See `run_status.json` and `controller.log` for progress. Do not treat directory existence as evidence of completion.

## Question

In the rebutal3 sensitivity `A_exclude_single_age`, the HVG–ElasticNet count-tuned clock dropped from 16.10 to 9.90 weeks pooled MAE. That result came from one partition and was the best of 24 sensitivity results. This run tests whether the exclusion rule holds up under the same evaluation that justified the locked rebutal3 model:

- main A (the same four study-grouped folds);
- the same 30 repeated study partitions (120 outer fits per procedure);
- B (NIHMS held out) and C (Winkler held out), each with training-only tuning;
- a final fit, for a possible replacement model.

It runs all three recipes, in both fixed-2,000 and count-tuning modes.

## The only change from rebutal3 production

`config/production.json` sets `"training_exclusion": "single_age_studies"`. In `scripts/run_models.R`, every outer task's training set drops any study with one distinct `age_week` within that training set. This happens before filtering, gene selection, inner CV or fitting. Held-out sets are unchanged, so the 136 A predictions, 18 B predictions and 29 C predictions are paired with rebutal3 on identical profiles.

The excluded studies are all young: Morris (16 profiles, 7 wk), CRA017492 (3, 8.7 wk), GSE252791 (3, 12 wk), Zhao (1, 12.9 wk) and GSE194196 (3, 16 wk). That is 26 profiles in total.

| Stage | Training profiles / studies |
|---|---|
| A fold 1 / 2 / 3 / 4 | 91/8, 65/7, 93/6, 81/6 |
| B (NIHMS held out) | 92/8 |
| C (Winkler held out) | 81/8 |
| Final | 110/9 (ages 8–83 wk) |
| A_repeated | 5–8 studies per training set |

Everything else is identical: grids, the one-SE rule, study-balanced weights, ln(1+CPM) inputs, seeds, fold assignments and the 30 repeated assignments. The prespecified sensitivities are not rerun. The inputs in `shared_inputs/` are byte copies of rebutal3's. `config/production_rebutal3_reference.json` is the unmodified rebutal3 config.

## Built-in reproduction check

In rebutal3, the `exclude_single_age` variant changed only training membership, not weights or target. Main A here should therefore reproduce rebutal3's `A_exclude_single_age` predictions; for example, HVG count-tuned pooled MAE should be 9.90.

## Interpretation rules (set before seeing results)

- Compare paired against the locked rebutal3 count-tuning results on identical held-out profiles: main A, the median/IQR across the 30 repeats, B and C.
- Consider replacing the locked model only if the gain holds across the repeated partitions and B/C do not get worse. Any replacement is still a model-development choice made after seeing outer results. Independent data is needed to validate it.
- A final model trained here excludes single-age studies. That restricts its training age support (the minimum training age becomes 8 wk). Report this alongside any adoption.

## Run

```bash
python3 scripts/run_all.py --workers=8 --phases='A;B,C;final;A_repeated'
```

The run resumes from checkpoints. Remove a stale `RUNNING.lock` only after confirming the recorded PID has exited.
