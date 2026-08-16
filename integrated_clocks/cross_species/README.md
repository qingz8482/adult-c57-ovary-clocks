# Cross-species integrated-clock analysis

This module evaluates mouse-to-human gene-direction conservation and frozen zero-shot transfer, then performs GTEx ovarian age-signal adequacy and final sensitivity analyses.

Run `run_all.sh`, followed by `3_gtex_age_signal_adequacy/run_all.sh` and `4_final_reporting_revision/run_all.sh`. Frozen mouse artifacts and the expression-alignment table are not distributed; place them under `data/integrated_clock_artifacts/` as specified in `config/analysis_config.yaml`.

GTEx and Human Protein Atlas resources are downloaded from the URLs frozen in the configuration. Their respective data-use terms apply. The workflow writes generated data, figures, reports, and logs within this module; those directories are ignored by Git.
