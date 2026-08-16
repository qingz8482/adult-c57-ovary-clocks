# Data access

All source datasets are publicly deposited. Dataset identifiers and repositories are listed in `metadata/dataset_accessions.tsv`; sample-level accessions are in `metadata/sample_annotations.tsv`.

The repository does not redistribute FASTQ files, count matrices, Seurat objects, pseudobulk matrices, model objects, or complete generated results. Users should retrieve source data from NCBI GEO/SRA/BioProject, GSA, or the repository named in the accession table and follow the corresponding study terms.

## Prepared data layout

`OVARY_CLOCK_DATA_ROOT` should contain the applicable parts of:

```text
bulkseq_dataset/<dataset>/input/
aggregated_pseudobulk/<dataset>/input_genelevel/
scRNA_dataset/<dataset>/metadata/
scRNA_dataset/<dataset>/cellranger_counts/
metacell_preprocessed/step_02_preprocessing_filtering/
```

Integrated inputs are gene-by-sample logCPM matrices accompanied by count matrices used for expression filtering. The documented filter retains a gene when count is at least 10 in at least 20% of samples within a dataset/modality. Metacell inputs are summed raw counts plus metacell annotations.

MSigDB gene sets are obtained at runtime through `msigdbr` and remain subject to MSigDB terms. The external tAge implementation is not redistributed in this repository.
