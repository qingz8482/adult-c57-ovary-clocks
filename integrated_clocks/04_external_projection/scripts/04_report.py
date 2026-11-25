#!/usr/bin/env python3
"""Package local provenance and a reviewable application report."""
import csv
import hashlib
import json
import shutil
from collections import Counter
from pathlib import Path

out = Path(__file__).resolve().parent.parent
project = out.parent
root = project.parent.parent
dev = project / 'model_development'
release = project / 'lockedmodel'
supp = root / 'bulk-pseudobulk_clock/methods_tables/publication_supplementary_tables'

def read(path):
    with open(out / path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

for name in ['Supplementary_Table_3b_external_projection_samples.tsv', 'Supplementary_Table_4_external_projection_contrasts.tsv', 'Supplementary_Table_5_primary_integrated_clock_projections.tsv']:
    shutil.copy2(supp / name, out / 'inputs' / name)
shutil.copy2(project / 'MODEL_SELECTION.json', out / 'provenance/MODEL_SELECTION.json')
shutil.copy2(release / 'manifest.json', out / 'provenance/locked_release_manifest.json')
shutil.copy2(release / 'checksums.tsv', out / 'provenance/locked_release_checksums.tsv')
samples = read('tables/external_projection_samples.tsv')
stats = read('tables/external_projection_contrast_statistics.tsv')
pred = read('external_projection_predictions.tsv')
coverage = read('tables/external_projection_gene_coverage.tsv')
counts = read('tables/cohort_counts.tsv')
groups = read('tables/all_sample_group_summary.tsv')
for dataset in ['GSE235021', 'PRJNA863443', 'PRJNA1002222']:
    rows = [r for r in samples if r['dataset_folder'] == dataset]
    summaries = [r for r in groups if r['dataset_cohort'].startswith(dataset + '__')]
    for filename, records in [('projection_metadata.tsv', rows), ('group_summary.tsv', summaries)]:
        with open(out / dataset / filename, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=list(records[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(records)
    ncontrasts = len({r['analysis_contrast_id'] for r in stats if r['dataset_folder'] == dataset})
    text = f'''# {dataset}: external projection samples

{len(rows)} sample-level profiles projected through both locked models; {ncontrasts} modality-specific contrasts per model.

- `projection_predictions.tsv`: ages in weeks, delta ages, gene coverage and development-overlap flags.
- `projection_metadata.tsv` and `group_summary.tsv`: sample metadata and all-group summaries.
- `contrast_statistics.tsv`: two-sided Welch tests, exact P, BH q, t, degrees of freedom, effects and 95% confidence intervals.
- `relative_to_controls.tsv`: per-contrast reference-centered delta ages with original reference and target samples retained.
- `pca_scores.tsv` and `pca_pc_contributions.tsv`: frozen PCA scores and additive contributions in weeks.
- `plots/`: delta-age distributions, annotated delta-age comparisons and predicted-versus-chronological-age plots.

See the parent README for model provenance, normalization, statistics and limitations. Both models are applied without refitting. Significance asterisks denote unadjusted P, not BH q.
'''
    if dataset == 'GSE235021':
        text += '\nThe 18 bulk profiles and 17 pseudobulk profiles are analyzed separately. The historical pooled-modality contrasts are split into four contrasts per modality.\n'
        if (out / dataset / 'plots/GSE235021__mixed_modality_all_groups_delta_age.pdf').exists():
            text += '\n[Mixed-modality delta-age plot](plots/GSE235021__mixed_modality_all_groups_delta_age.pdf): circles denote bulk RNA-seq and triangles denote pseudobulk, with side-by-side boxes and modality-specific Welch annotations. The existing tests remain separate by modality. Source data, annotations, counts, caption and validation are in `tables/mixed_modality_delta_age_*`; reproduce with `scripts/11_gse235021_mixed_modality_delta.R`.\n'
    (out / dataset / 'README.md').write_text(text)
models = list(dict.fromkeys(r['model'] for r in stats))
ntrain = sum(r['locked_training_sample'] == 'TRUE' for r in samples)
nshared = sum(r['locked_training_study'] == 'TRUE' for r in samples)
lowcov = [r for r in coverage if float(r['min_gene_coverage']) < .9]
negative = Counter(r['model'] for r in pred if float(r['predicted_age_week']) < 0)
permodel = {m: (sum(r['locked_training_sample'] == 'TRUE' for r in pred if r['model'] == m), sum(r['locked_training_study'] == 'TRUE' for r in pred if r['model'] == m)) for m in models}
lines = [
 '# rebutal4 locked-model application — external projection samples', '',
 'The two locked clocks of release `hvg_exclusion_vbpca_20261012` (`../lockedmodel`) were applied without refitting, recalibration, batch correction or feature selection: **HVG–ElasticNet trained without single-age studies** (110 profiles, 9 studies) and **Version B–PCA/glmnet** (136 profiles, 14 studies; unchanged from rebutal3). Version B–ElasticNet is not part of this release. Figures follow the rebutal3 application layouts, adapted to two clocks.', '',
 'Version B–PCA outputs reproduce the completed rebutal3 application exactly; see `provenance/rebutal3_pca_reuse_check.tsv`. HVG–ElasticNet outputs are new.', '',
 'The projection set contains **281 samples from 20 modality-specific datasets representing 19 study accessions**: 231 bulk RNA-seq samples from 18 datasets and 50 sample-level pseudobulk profiles from two datasets. All 19 GSE154890 profiles are excluded because that study contributed to model development. GSE235021 consists of 18 bulk plus 17 pseudobulk profiles, analysed as separate datasets; PRJNA863443 contributes the other 33 pseudobulk profiles. No samples were dropped for an adverse projection result.', '',
 '## Outputs', '',
 '- [All predictions](external_projection_predictions.tsv): 562 sample-by-model predictions in weeks, delta ages, gene coverage, age-range flags and model-specific development overlap.',
 '- [Contrast statistics](external_projection_contrast_statistics.tsv): 44 modality-specific comparisons × two clocks = 88 two-sided Welch tests, effect estimates, 95% confidence intervals, t, degrees of freedom, P and BH q.',
 '- [Complete reference-folder review](REFERENCE_REVIEW.md) and [file-by-file coverage](provenance/reference_file_coverage.tsv).',
 '- [Diagnostics](diagnostics/): raw and matched-control relative delta age with named groups and Welch annotations; concordance, disagreement, QC/OOD and interpretation grading.',
 '- [Overview plots](plots/), and dataset folders [GSE235021](GSE235021/), [PRJNA863443](PRJNA863443/), [PRJNA1002222](PRJNA1002222/).',
 '- Retired duplicate diagnostics (as in rebutal3) are not generated; see [retired output scope](provenance/retired_diagnostic_outputs.tsv).',
 '- [Source tables](tables/), [input snapshots](inputs/), [provenance and validation](provenance/), [scripts](scripts/) and [logs](logs/). `provenance/port_from_rebutal3.py` records every change from the rebutal3 scripts.', '',
 '## Frozen models and expression preprocessing', '',
 '| Model | Training set | Input genes | Active predictors |', '|---|---|---:|---|',
 '| HVG–ElasticNet (single-age excluded) | 110 profiles / 9 studies | 3,000 | 79 genes |',
 '| Version B–PCA/glmnet | 136 profiles / 14 studies | 3,000 | PC1, PC2, PC5, PC8 (10 PC inputs) |', '',
 'Release manifest, checksums and selection record are copied into provenance. Expression is ln(1 + 10^6 × gene quantity / full available prefilter sample total), with the frozen development identifier dictionary; duplicate symbols are summed on the quantity scale and denominators are computed before clock-gene subsetting. Missing model genes receive the frozen training mean. The input RDS retains the union of genes required by both models. Projection inputs were rebuilt from source because the HVG model uses a different gene set; every recovered dataset was checked against its existing filtered quantity matrix (maximum relative discrepancy < 10^-6).', '',
 '## Statistical interpretation and limitations', '',
 'Delta age = predicted age − chronological age, in weeks. Tests compare target versus reference delta age within the same dataset and modality, using sample profiles as units, with unpaired two-sided Welch t-tests. Effect = target mean − reference mean delta age. Relative delta age subtracts the corresponding observed reference mean. Asterisks denote unadjusted P (* <0.05, ** <0.01, *** <0.001, **** <0.0001); BH q-values across all 88 tests and within each model are supplied. Two additional clock-specific tests restore the GSE279573 oil-free rescue versus WT reference contrast (`tables/additional_reference_contrast_statistics.tsv`; 90-test family in `tables/all_reference_contrast_statistics.tsv`).', '',
 '**Development overlap is model-specific.** ' + '; '.join(f'{m}: {permodel[m][0]} profiles used to fit it, {permodel[m][1]} profiles from its development studies' for m in models) + '. HVG–ElasticNet never trained on the single-age studies CRA017492, GSE194196 and GSE252791, so their profiles are external to it but not to Version B–PCA. Overlapping projections are not independent external validation.', '',
 'Natural-age placement comparisons, neonatal/prepubertal profiles and extrapolations beyond the training age range (HVG 8–83 weeks; PCA 7–83 weeks) are descriptive. Delta-age shifts alone do not establish rejuvenation or accelerated biological ageing. CRA003645 has six developmental profiles and no prespecified contrast. GSE228975 sample mapping follows SRA titles and replicate numbers, as in rebutal3.', '',
 f'{len(lowcov)} dataset-by-model combinations have less than 90% measured input-gene coverage.', '',
 'Negative predicted ages are retained and flagged, without clipping: ' + ', '.join(f'{model}: {negative[model]}' for model in models) + '. These are unconstrained model outputs under external projection, not literal negative ages.', '',
 '## Result counts', '',
 '| Model | Contrasts | Nominal P < 0.05 | BH q < 0.05 across 88 tests |', '|---|---:|---:|---:|',
 ]
for model in models:
    ss = [r for r in stats if r['model'] == model]
    lines.append(f"| {model} | {len(ss)} | {sum(float(r['welch_p']) < .05 for r in ss)} | {sum(float(r['bh_q_all_tests']) < .05 for r in ss)} |")
lines += ['', 'These counts include the explicitly labelled placement and developmental contrasts; they are not a count of demonstrated interventions.', '',
 '## Dataset accounting', '', '| Study | Modality | Samples | Overlap with any locked model |', '|---|---|---:|---:|']
for row in counts:
    lines.append(f"| {row['dataset_folder']} | {row['modality']} | {row['n_samples']} | {row['n_training_samples']} |")
lines += ['', '## Validation and reproduction', '',
 'The release checksum verification passed. Frozen glmnet predictions were independently replayed from exported gene-scale coefficients; PCA scores × coefficients plus the frozen intercept reconstruct PCA predictions. Shared sample expression was compared with the development input RDS. Every Welch statistic was recomputed from the exported plotted samples. The checks and numerical errors are saved in provenance.', '',
 'Run `bash scripts/run_all.sh` from this folder to reproduce using the existing R environment and the recorded local source data. For a figure/derived-analysis refresh, run R script `scripts/07_complete_reference_outputs.R`, Python script `scripts/06_validate_plot_revision.py`, R script `scripts/09_validate_derived_analyses.R`, then Python scripts `scripts/08_audit_reference_completeness.py`, `scripts/10_review_application_files.py` and `scripts/04_report.py`, in that order. Use Rscript for R files and python3 for Python files. The older 05 entry point delegates to 07. The scripts write only into this application folder. Input supplementary tables are provenance snapshots; historical Table 5 numerical results are not reused as current predictions.', '']
(out / 'README.md').write_text('\n'.join(lines))
comprehensive = json.loads((out / 'provenance/comprehensive_validation.json').read_text()) if (out / 'provenance/comprehensive_validation.json').exists() else {}
(out / 'run_status.json').write_text(json.dumps({'status': 'complete' if comprehensive.get('status') == 'PASS' else 'awaiting_completeness_audit', 'comprehensive_review': comprehensive, 'date': '2026-10-12', 'timezone': 'Asia/Singapore', 'samples': 281, 'modality_specific_datasets': 20, 'study_accessions': 19, 'models': 2, 'predictions': len(pred), 'prespecified_modality_specific_contrasts': 44, 'welch_tests': len(stats), 'excluded_study': 'GSE154890'}, indent=2) + '\n')
with open(out / 'provenance/output_checksums.tsv', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['relative_path', 'bytes', 'sha256'])
    for path in sorted(out.rglob('*')):
        if path.is_file() and path.name != 'output_checksums.tsv' and 'logs' not in path.relative_to(out).parts:
            writer.writerow([str(path.relative_to(out)), path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest()])
print('Packaged report, supplementary input snapshots, provenance and output checksums.')
