# Metadata

- `sample_annotations.tsv`: 411 public sample annotations (294 bulk RNA-seq and 117 scRNA pseudobulk). Local expression and source-file paths were removed.
- `cohort_manifest.tsv`: modality/dataset-level inventory derived from the sample annotations.
- `dataset_accessions.tsv`: dataset-to-repository/accession lookup.
- `split_assignments.tsv`: frozen integrated-clock sample roles and balance fields.
- `metacell_cv_assignments.tsv`: canonical metacell-level fold assignments.
- `manual_mappings/`: approved sample/group mappings required by the frozen integrated workflow.

`sample_key` is the stable integrated-clock join key. `dataset_folder` is the stable dataset key. Empty values mean unavailable or not applicable; they must not be interpreted as negative evidence.

The 411-row inventory is broader than the frozen 117-sample model-ready matrix. Eligibility flags and `split_role` determine whether a row contributed to model development, external projection, or descriptive placement.
