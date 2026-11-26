# Reproducibility boundary

**Included:**

- the analysis code used for the revised manuscript;
- the two locked integrated clocks, with replay checks and checksums;
- the published supplementary tables;
- the fold, contrast and metacell assignments the pipelines read;
- environment specifications;
- tests.

**Excluded:**

- expression inputs and Seurat objects;
- generated result trees, figures and manuscript files;
- gene-set snapshots;
- the third-party tAge software.

The integrated-clock predictions can be reproduced exactly from the locked models given correctly prepared ln(1 + CPM) inputs. `integrated_models/validation/replay_checks.csv` records that each locked object reproduces its 110 or 136 training-profile predictions exactly. The exported coefficients reproduce them to within 10⁻¹³ weeks. Repeating model development needs the prepared 136-profile matrix, which `integrated_clocks/01_input_preparation` builds from the public datasets.

Interpretation limits stated in the manuscript:

- Nested-CV estimates describe the development procedure. The EN single-age exclusion rule was adopted after outer results were seen, so its estimates are development estimates, not independent validation.
- External projections are uncalibrated deployment outputs. A delta-age shift is a descriptive transcriptomic shift relative to matched controls. It is not by itself proof of biological ageing or rejuvenation.
- Cross-species outputs are in mouse weeks and are standardized within each cohort before testing.

The `v0.1-initial-submission` tag holds the earlier code from the original submission (fixed train/validation/test split, v1 metacell labels). It is superseded and is kept only for the review record.
