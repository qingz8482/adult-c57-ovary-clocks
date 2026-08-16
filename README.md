# Adult C57 ovary clocks

Publication code for three mouse ovarian transcriptomic ageing-clock analyses:

1. an integrated bulk RNA-seq plus scRNA pseudobulk ElasticNet clock;
2. an integrated bulk RNA-seq plus scRNA pseudobulk PCA/glmnet clock; and
3. cell-type-specific scRNA-seq metacell ElasticNet clocks.

The frozen integrated benchmark contains 117 model-ready biological samples (83 training, 14 validation, and 20 held-out test samples). The expanded curated inventory contains 411 sample rows: 294 bulk RNA-seq and 117 scRNA pseudobulk profiles. The external-projection collection contains 300 samples, of which 288 are analysis-domain samples and 12 are developmental/out-of-domain placements. It supports 41 prespecified contrasts, including 34 primary adult perturbation contrasts.

## Repository map

- [`integrated_clocks/`](integrated_clocks/README.md): shared preprocessing, fitting, evaluation, projection, and separate ElasticNet and PCA/glmnet interpretation code.
- [`metacell_clocks/`](metacell_clocks/README.md): the active July 2026 metacell-level internal-CV workflow.
- [`metadata/`](metadata/README.md): public annotations, accessions, frozen split assignments, and approved mappings. Expression-file paths were intentionally removed.
- [`docs/`](docs/workflow.md): data access, environment, reproducibility boundary, and manuscript-to-code mapping.
- [`tests/`](tests/README.md): static validation and synthetic-data smoke tests.

## Installation

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate adult-c57-ovary-clocks
```

Alternatively, restore R dependencies with `renv::restore()` using `renv.lock`. See [`docs/environment.md`](docs/environment.md) for external resources and optional packages.

## Data and output locations

Large expression matrices and generated results are not distributed in this repository. Set:

```bash
export OVARY_CLOCK_DATA_ROOT=/path/to/prepared/data
export OVARY_CLOCK_OUTPUT_ROOT=/path/to/new/outputs
```

Expected input layouts and public accession links are documented in [`docs/data_access.md`](docs/data_access.md). Full analyses must be run with the exact public input matrices described there; the included synthetic smoke tests validate software plumbing only.

## Canonical workflow

Run the scripts in the order documented in [`docs/workflow.md`](docs/workflow.md). The integrated scripts preserve the frozen family-specific selector lock: HVG for ElasticNet and Version B for PCA/glmnet. Metacell evaluation holds out metacells, not biological samples.

Validate a clone before use or release:

```bash
bash tests/validate_repository.sh
bash tests/run_smoke_tests.sh
```

## Citation and license

Citation metadata are in [`CITATION.cff`](CITATION.cff). The intended publication snapshot will be tagged `v1.0.0-publication` and archived on Zenodo. Code is released under the [MIT License](LICENSE), subject to institutional confirmation. Public datasets and third-party software retain their original licenses and terms.
