#!/usr/bin/env python3
"""Non-visual extraction of frozen tAge weights; R makes every figure."""
from pathlib import Path
import hashlib
import json
import joblib
import numpy as np
import pandas as pd
import sklearn

out = Path(__file__).resolve().parent.parent
ovary = out.parent.parent.parent
tage = ovary / 'bulk-pseudobulk_clock/results/external_projection/tAge'
universe = pd.read_csv(out / 'tables/current_ovary_gene_universe.tsv', sep='\t').gene.tolist()
universe_set = set(universe)
canonical = {g.upper(): g for g in universe}
assert len(canonical) == len(universe), 'Ambiguous case-normalized ovary gene symbols'
mapping = pd.read_csv(tage / 'audit/tage_official_mouse_gene_mapping.tsv', sep='\t', dtype=str)
counts = mapping.dropna(subset=['mouse_entrez','original_genes']).groupby(['mouse_entrez','original_genes']).size().reset_index(name='support')
counts['current_universe_match'] = counts.original_genes.str.upper().isin(canonical)
lookup, audit = {}, []
for entrez, g in counts.groupby('mouse_entrez', sort=False):
    g = g.copy()
    g = g.sort_values(['current_universe_match','support','original_genes'], ascending=[False,False,True])
    symbol = g.iloc[0].original_genes
    symbol = canonical.get(symbol.upper(), symbol)
    lookup[entrez] = symbol
    audit.append(dict(entrez=entrez,gene=symbol,n_candidates=len(g),current_universe_match=symbol in universe_set))
pd.DataFrame(audit).to_csv(out / 'provenance/tage_entrez_symbol_mapping.tsv',sep='\t',index=False)
records=[]
for variant in ['scaleddiff','yugenediff']:
    path = tage / 'environment/models' / f'EN_Chronoage_Mouse_Multitissue_{variant}.pkl'
    pipeline = joblib.load(path)
    features = np.asarray(pipeline.feature_names_in_, dtype=str)
    selected = features[pipeline.named_steps['featureselection'].get_support()]
    weights = np.asarray(pipeline.named_steps['estimator'].coef_)
    assert len(selected)==len(weights)
    raw = pd.DataFrame(dict(entrez=selected,coefficient=weights))
    raw['gene']=[lookup.get(g,'ENTREZ_'+g) for g in selected]
    raw.to_csv(out / f'tables/tage_{variant}_raw_features.tsv',sep='\t',index=False)
    grouped = raw.groupby('gene',sort=False).agg(coefficient=('coefficient','sum'),entrez=('entrez',lambda z:';'.join(sorted(set(z)))))
    grouped.reset_index().to_csv(out / f'tables/tage_{variant}_features.tsv',sep='\t',index=False)
    records.append(dict(variant=variant,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),n_required=len(features),n_selected=len(selected),n_nonzero=int(np.count_nonzero(weights)),n_symbols=len(grouped),sklearn=sklearn.__version__,steps=str(pipeline.steps)))
(out/'provenance/tage_model_extraction.json').write_text(json.dumps(records,indent=2))
print('PASS: tAge weights extracted directly from both frozen serialized models.')
