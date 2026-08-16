# Tests

- `validate_repository.sh` checks required structure, metadata accounting, forbidden paths/directories, file sizes, and script syntax.
- `run_smoke_tests.sh` runs metadata checks plus deterministic synthetic ElasticNet, PCA/glmnet, and cell-class metacell fits.
- `smoke_models.R` uses synthetic data only and does not test biological results.

Run both shell scripts from any working directory; they locate the repository from their own paths.
