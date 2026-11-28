# Data access

All source datasets are public. Dataset accessions and repositories are listed in `metadata/dataset_accessions.tsv`. Sample-level accessions are in `metadata/manuscript_tables/SD2a_model_development_profiles.tsv` and `SD2b_external_projection_profiles.tsv`.

This repository does not redistribute FASTQ files, count matrices, Seurat objects, pseudobulk or metacell matrices, or generated results. The two locked integrated clocks are included in `integrated_clocks/final_models/`.

## Inputs expected by the code

All expression inputs are raw gene-level counts; no author-normalized values are used.

- **Bulk RNA-seq.** Raw counts deposited by the authors (featureCounts), or counts reconstructed from FASTQ:
  - Salmon 1.10.3 against the GENCODE M39/GRCm39 transcriptome for paired-end libraries;
  - STAR 2.7.11b for strand-specific QuantSeq libraries.
- **scRNA-seq.**
  - Author Cell Ranger matrices or processed Seurat objects, or FASTQ reprocessed with Cell Ranger 10.0.0 against GRCm39 (Ensembl 110).
  - Cells with < 200 detected genes or > 20% mitochondrial reads are removed.
  - Pseudobulk profiles are raw counts summed per biological sample.
  - Metacells are raw counts summed per dataset × sample × cell class.
- **Human cohorts.**
  - GTEx v11 ovary gene counts with public age brackets, RIN and ischaemic time; the protected-access fields are not used.
  - GSE202601 whole-ovary snRNA-seq, summed per donor.
- **Orthology.** Ensembl BioMart (release 116) one-to-one mouse–human orthologues.

## Gene sets and third-party models

- MSigDB 2026.1.Mm (Hallmark, M2 canonical pathways, M5 GO BP, M7) is loaded through `msigdbr`.
- clusterProfiler/org.Mm.eg.db (GO, KEGG, Reactome) is used for the EN non-zero-gene ORA.
- GO BP, Reactome and KEGG Mouse are used for the tAge pathway comparison.
- Gene-set snapshots are not redistributed, because KEGG and MSigDB have their own terms.
- The tAge models (v1.0.1) come from Zenodo record 18763485.
