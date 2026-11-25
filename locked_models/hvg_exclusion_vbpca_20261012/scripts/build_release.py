#!/usr/bin/env python3
"""Assemble the 2026-10-12 locked release by byte-copying saved artifacts; no model is refitted."""
import hashlib, json, shutil
from pathlib import Path

REL = Path(__file__).resolve().parent.parent
REBUTAL = REL.parent.parent
R3 = REBUTAL / 'rebutal3_batcheffects removal_CV_LOO/model_development'
R4 = REBUTAL / 'rebutal4/model_development'
R3_LOCK = R3 / 'locked_models/count_tuning_20261011'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

sources = []
def copy(src, dst):
    dst = REL / dst
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        raise SystemExit(f'Refusing to overwrite locked file: {dst}')
    shutil.copy2(src, dst)
    assert sha(src) == sha(dst)
    sources.append({'source': str(src), 'locked_path': str(dst.relative_to(REL)), 'sha256': sha(dst)})

# HVG-ElasticNet, exclusion of single-age training studies (rebutal4 final fit on 110 profiles).
hvg = 'HVG-ElasticNet-exclude-single-age'
fm = R4 / 'HVG-ElasticNet/count_tuning/results/final_model'
for f in ['model.rds', 'coefficients.csv', 'genes.txt', 'hyperparameters.json', 'intercept.txt', 'scaler.csv']:
    copy(fm / f, f'{hvg}/{f}')
for f in ['config_used.json', 'run_signature.txt', 'sessionInfo.txt']:
    copy(R4 / 'HVG-ElasticNet/count_tuning/results' / f, f'{hvg}/provenance/{f}')
for f in ['count_tuning_candidate_summary.csv', 'inner_split_audit.csv']:
    copy(R4 / 'HVG-ElasticNet/tuning_audit/final_r00_f0' / f, f'{hvg}/provenance/{f}')

# Version B-PCA/glmnet, unchanged from the 2026-10-11 locked release (all 136 profiles).
pca = 'VersionB-PCA-glmnet'
for p in sorted((R3_LOCK / pca).rglob('*')):
    if p.is_file():
        copy(p, f'{pca}/{p.relative_to(R3_LOCK / pca)}')

# Frozen prediction code (byte-identical in rebutal3, its locked release and rebutal4).
for f in ['core.R', 'predict_frozen.R']:
    assert sha(R3_LOCK / 'scripts' / f) == sha(R4 / 'scripts' / f)
    copy(R4 / 'scripts' / f, f'scripts/{f}')

# Shared provenance.
for f in ['sample_metadata.csv', 'source_checksums.tsv', 'normalization_comparison.tsv']:
    copy(R4 / 'shared_inputs' / f, f'provenance/{f}')
copy(R4 / 'config/production.json', 'provenance/rebutal4_production_config.json')
copy(R4 / 'README.md', 'provenance/rebutal4_development_README.md')
copy(R4 / 'COMPARISON_vs_locked_rebutal3.md', 'provenance/rebutal4_comparison_vs_rebutal3.md')
copy(R3 / 'RESULTS_REPORT.md', 'provenance/rebutal3_results_report.md')
copy(R3_LOCK / 'manifest.json', 'provenance/previous_release_count_tuning_20261011_manifest.json')

def hp(path):
    return json.loads((REL / path).read_text())

manifest = {
    'release': 'hvg_exclusion_vbpca_20261012',
    'locked_date': '2026-10-12',
    'timezone': 'Asia/Singapore',
    'selection_authority': 'User explicitly selected HVG-ElasticNet with single-age training-study exclusion (rebutal4) and the previously locked Version B-PCA/glmnet (rebutal3 count_tuning_20261011).',
    'supersedes': 'locked_models/count_tuning_20261011 for deployment of these two clocks; Version B-ElasticNet is not carried into this release. The previous release is unchanged.',
    'selection_scope': 'The exclusion rule was adopted after seeing rebutal3 sensitivity and rebutal4 outer results; this is a model-development decision, not independent validation.',
    'models': [
        {'recipe': 'HVG-ElasticNet', 'artifact': f'{hvg}/model.rds', 'procedure': 'count_tuning',
         'training_rule': 'exclude studies with a single distinct age_week from training',
         'source_dev': str(R4), 'source_model': 'HVG-ElasticNet/count_tuning/results/final_model/model.rds',
         'training_profiles': 110, 'training_studies': 9,
         'excluded_studies': ['2022_eLife_Morris_cycling_murine_ovary_atlas', 'CRA017492', 'GSE252791',
                              '2025_CommBio_Zhao_developing_ovary_scRNA_spatial', 'GSE194196'],
         'input_genes': 3000, 'nonzero_effective_gene_coefficients': 79, 'PC_inputs': 0, 'nonzero_PCs': [],
         'hyperparameters': hp(f'{hvg}/hyperparameters.json')},
        {'recipe': 'VersionB-PCA/glmnet', 'artifact': f'{pca}/model.rds', 'procedure': 'count_tuning',
         'training_rule': 'all studies (rebutal3 production)',
         'source_dev': str(R3), 'source_model': 'VersionB-PCA/glmnet/count_tuning/results/final_model/model.rds',
         'training_profiles': 136, 'training_studies': 14, 'excluded_studies': [],
         'input_genes': 3000, 'nonzero_effective_gene_coefficients': 3000, 'PC_inputs': 10, 'nonzero_PCs': [1, 2, 5, 8],
         'hyperparameters': hp(f'{pca}/hyperparameters.json')},
    ],
    'sources': sources,
}
(REL / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
print(f'Copied {len(sources)} files')
