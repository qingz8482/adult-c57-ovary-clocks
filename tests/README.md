# Tests

- `validate_repository.sh`: checks required paths, scans for cluster paths, credentials, files over 50 MiB and generated or archive folders, parses every R and Python script, validates the metadata and checks that every path in `docs/manuscript_code_map.tsv` exists.
- `run_smoke_tests.sh`: runs three checks:
  - locked-release checksums;
  - metadata validation (136 profiles and 14 studies; the 26-profile single-age exclusion leaving 110 EN training profiles; Stage A folds of 36/48/20/32 with whole studies held out; 30 repeats; 281 external profiles; 36 contrasts; 7,049 metacells in eight classes with five folds each);
  - `test_locked_models.R`: synthetic input through `predict_frozen.R` must match a manual replay from the exported coefficients, missing-gene imputation must behave as documented, and the exports must match Supplementary Data 3a;

The tests use synthetic data only and do not reproduce biological results. Set `RSCRIPT` if `Rscript` is not on `PATH`.

`integrated_clocks/02_model_development/tests/test_pipeline.R` checks the nested-CV pipeline for leakage, split and selection invariants. It needs the prepared 136-profile matrix, so it is not part of the smoke tests.
