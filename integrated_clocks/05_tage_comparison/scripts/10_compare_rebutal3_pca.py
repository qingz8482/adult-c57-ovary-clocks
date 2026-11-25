#!/usr/bin/env python3
"""Version B-PCA is unchanged between releases: compare its compare_tAge outputs with rebutal3."""
from pathlib import Path
import numpy as np
import pandas as pd

out = Path(__file__).resolve().parent.parent
r3 = out.parent.parent / 'rebutal3_batcheffects removal_CV_LOO/compare_tAge'
def read(base, name): return pd.read_csv(base / name, sep='\t')
PCA = 'VersionB-PCA/glmnet'
rows = []
def add(check, value, exact):
    rows.append(dict(check=check, value=value, must_be_exact=exact, passed=(not exact) or value < 1e-9))

a, b = read(r3, 'tables/sample_contrast_predictions.tsv'), read(out, 'tables/sample_contrast_predictions.tsv')
m = a[a.model == PCA].merge(b[b.model == PCA], on=['analysis_contrast_id', 'sample_key'], suffixes=('_3', '_4'))
assert len(m) == (a.model == PCA).sum() == (b.model == PCA).sum()
add('PCA sample predictions max |difference| (weeks)', float(np.abs(m.value_weeks_3 - m.value_weeks_4).max()), True)
for model in [PCA, 'tAge scaled-diff', 'tAge YuGene-diff']:
    a, b = read(r3, 'tables/contrast_effects.tsv'), read(out, 'tables/contrast_effects.tsv')
    m = a[a.model == model].merge(b[b.model == model], on='analysis_contrast_id', suffixes=('_3', '_4'))
    add(f'{model} contrast effects max |difference| (weeks)', float(np.abs(m.effect_weeks_3 - m.effect_weeks_4).max()), True)
a, b = read(r3, 'tables/ovary_pca_backprojected_features.tsv'), read(out, 'tables/ovary_pca_backprojected_features.tsv')
m = a.merge(b, on='gene', suffixes=('_3', '_4'))
assert len(m) == len(a) == len(b)
add('PCA backprojected gene weights max |difference|', float(np.abs(m.standardized_coefficient_3 - m.standardized_coefficient_4).max()), True)
a, b = read(r3, 'tables/feature_direction_statistics_all_pca.tsv'), read(out, 'tables/feature_direction_statistics_all_pca.tsv')
m = a[a.comparison.str.startswith('VersionB-PCA')].merge(b, on='comparison', suffixes=('_3', '_4'))
add('PCA-vs-tAge shared-gene n difference', float(np.abs(m.n_3 - m.n_4).max()), True)
add('PCA-vs-tAge gene-weight correlation max |difference|', float(np.abs(m.pearson_3 - m.pearson_4).max()), True)
u3, u4 = read(r3, 'tables/feature_overlap_gene_universe.tsv'), read(out, 'tables/feature_overlap_gene_universe.tsv')
add('shared pathway background genes, rebutal3', len(u3), False)
add('shared pathway background genes, rebutal4 (union)', len(u4), False)
a, b = read(r3, 'tables/pathway_overlap_enrichment.tsv'), read(out, 'tables/pathway_overlap_enrichment.tsv')
key = ['gene_set', 'database', 'pathway', 'method']
m = a[a.gene_set.str.contains('pca')].merge(b[b.gene_set.str.contains('pca')], on=key, suffixes=('_3', '_4'))
lp3, lp4 = -np.log10(m.p_value_3.clip(lower=1e-300)), -np.log10(m.p_value_4.clip(lower=1e-300))
add('PCA pathway tests matched', len(m), False)
add('PCA pathway tests: Spearman of -log10 P (rebutal3 vs rebutal4)', float(pd.Series(lp3).corr(pd.Series(lp4), method='spearman')), False)
add('PCA pathway tests: median |delta -log10 P|', float(np.median(np.abs(lp3 - lp4))), False)
add('PCA pathway tests: FDR<0.05 status changed', int(((m.FDR_3 < .05) != (m.FDR_4 < .05)).sum()), False)
res = pd.DataFrame(rows)
res.to_csv(out / 'provenance/rebutal3_pca_reuse_check.tsv', sep='\t', index=False)
print(res.to_string(index=False))
assert res.passed.all(), res[~res.passed]
print('PASS: Version B-PCA predictions, effects, weights and gene concordance identical to rebutal3; pathway shifts reported.')
