#!/usr/bin/env python3
"""Independent validation of pathway statistics; no graphics."""
from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
import pandas as pd
from scipy import stats

out=Path(__file__).resolve().parent.parent
def read(p):return pd.read_csv(out/p,sep='\t')
checks=[]
def check(name,ok,detail=''):
    checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name+': '+str(detail))
e=read('tables/pathway_overlap_enrichment.tsv')
p=read('tables/pathway_program_classification.tsv')
s=read('tables/pathway_program_classification_full_family_sensitivity.tsv')
chosen=read('tables/shared_ageing_pathways_source.tsv')
u=read('tables/feature_overlap_gene_universe.tsv').gene.tolist();uset=set(u)
q=read('tables/pathway_query_gene_membership.tsv')
members=read('tables/pathway_library_membership.tsv')
libraries={(db,pathway):set(g.gene) for (db,pathway),g in members.groupby(['database','pathway'],sort=False)}
queries={nm:set(g.loc[g.in_common_universe,'gene']) for nm,g in q.groupby('gene_set')}
check('all_pathway_genes_in_background',set(members.gene)<=uset)
check('pathway_size_bounds',all(5<=len(g)<=500 for g in libraries.values()))
check('unique_enrichment_keys',not e.duplicated(['gene_set','database','pathway','method']).any())
check('valid_probabilities',e[['p_value','FDR','FDR_full_pathway_family']].ge(0).all().all() and e[['p_value','FDR','FDR_full_pathway_family']].le(1).all().all())
ora=e[e.method=='ORA_Fisher'];err=0
for r in ora.itertuples():
    a=queries[r.gene_set];b=libraries[r.database,r.pathway];k=len(a&b)
    val=stats.hypergeom.sf(k-1,len(u),len(b),len(a))
    err=max(err,abs(val-r.p_value))
    if k!=r.overlap_gene_count or len(a)!=r.query_size or set(str(r.shared_genes).split(';'))!=a&b:raise AssertionError('ORA source membership mismatch')
check('independent_ORA_probabilities_and_membership',err<1e-9,err)
rankmap={
 'ovary_elasticnet_ranked':('ovary_elasticnet_features.tsv','standardized_coefficient'),
 'ovary_pca_backprojected_ranked':('ovary_pca_backprojected_features.tsv','standardized_coefficient'),
 'tage_scaleddiff_ranked':('tage_scaleddiff_features.tsv','coefficient'),
 'tage_yugenediff_ranked':('tage_yugenediff_features.tsv','coefficient')}
positions={g:i for i,g in enumerate(u)}
indices={k:np.array([positions[g] for g in genes]) for k,genes in libraries.items()}
err=0
for nm,(fn,col) in rankmap.items():
    w=read('tables/'+fn).set_index('gene')[col].reindex(u).fillna(0).to_numpy()
    ranks=stats.rankdata(abs(w),method='average');N=len(w)
    for r in e[e.gene_set==nm].itertuples():
        ii=indices[r.database,r.pathway];n=len(ii)
        z=(ranks[ii].sum()-n*(n+1)/2-n*(N-n)/2)/np.sqrt(n*(N-n)*(N+1)/12)
        err=max(err,abs(z-r.NES_or_OR),abs(stats.norm.sf(z)-r.p_value))
check('independent_absolute_rank_tests',err<1e-9,err)
def bh(vals,total=None):
    a=np.asarray(vals);order=np.argsort(a);n=len(a)
    q=np.minimum.accumulate((a[order]*(total or n)/np.arange(1,n+1))[::-1])[::-1]
    res=np.empty(n);res[order]=np.minimum(q,1);return res
err=0
for _,g in e.groupby(['gene_set','database','method'],sort=False):
    full=int(g.n_eligible_pathways.iloc[0]) if g.method.iloc[0]=='ORA_Fisher' else len(g)
    err=max(err,np.max(abs(bh(g.p_value)-g.FDR)),np.max(abs(bh(g.p_value,full)-g.FDR_full_pathway_family)))
check('independent_BH_and_full_family_sensitivity',err<1e-9,err)
for table,col,label in [(p,'FDR','legacy'),(s,'FDR_full_pathway_family','full_family')]:
    recalcs=[]
    for (db,pathway),g in e.groupby(['database','pathway']):
        ov=g.gene_set.str.startswith('ovary_');a=g.loc[ov,col].min();b=g.loc[~ov,col].min()
        cls='conserved' if a<.05 and b<.05 else 'ovary_specific' if a<.05 else 'tage_specific' if b<.05 else 'not_FDR_supported'
        score=-np.log10(max(a,b,1e-300)) if cls=='conserved' else -np.log10(max(a,1e-300)) if cls=='ovary_specific' else -np.log10(max(b,1e-300)) if cls=='tage_specific' else 0
        recalcs.append(dict(database=db,pathway=pathway,recomputed_class=cls,recomputed_score=score))
    joined=table.merge(pd.DataFrame(recalcs),on=['database','pathway'])
    check('classification_and_scores_'+label,len(joined)==len(table) and joined.program_class.eq(joined.recomputed_class).all() and np.max(abs(joined.support_score-joined.recomputed_score))<1e-9)
expected=pd.concat([p[p.program_class==cls].sort_values(['support_score','pathway','database'],ascending=[False,True,True]).head(10) for cls in ['conserved','ovary_specific']])
check('plotted_pathways_ranked_without_manual_selection',set(zip(expected.database,expected.pathway))==set(zip(chosen.database,chosen.pathway)))
for r in read('provenance/pathway_input_integrity.tsv').itertuples():
    check('unchanged_input_'+Path(r.path).name,hashlib.md5(Path(r.path).read_bytes()).hexdigest()==r.md5_before)
replay=read('provenance/pathway_historical_replay.tsv')
check('historical_full_replay',replay.passed.all())
pdf=out/'figures/shared_ageing_pathways.pdf';png=pdf.with_suffix('.png')
txt=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True)
check('readable_figure_exports',pdf.stat().st_size>1000 and png.stat().st_size>1000 and 'Conserved programs' in txt and 'Ovary-supported only' in txt)
pd.DataFrame(checks).to_csv(out/'provenance/pathway_validation_checks.tsv',sep='\t',index=False)
diff=p.merge(s,on=['database','pathway'],suffixes=('_legacy','_full'))
changed=diff[diff.program_class_legacy!=diff.program_class_full]
changed.to_csv(out/'tables/pathway_classification_sensitivity_changes.tsv',sep='\t',index=False)
marker='## Validation and sensitivity results'
report=(out/'PATHWAY_ANALYSIS.md').read_text().split(marker)[0].rstrip()
section=[marker,'',f'All **{len(checks)} pathway checks passed**, including an independent SciPy recalculation of every current ORA probability and rank statistic, both BH families, classifications, score-based plot selection, historical numerical replay, and unchanged source/previous-figure hashes.','',f'Including zero-overlap pathways in the ORA multiplicity family changes **{len(changed)} of {len(p)}** pathway classes. The complete changes are in `tables/pathway_classification_sensitivity_changes.tsv`.','', '| Class | Historical multiplicity rule | Full eligible-pathway sensitivity |','|---|---:|---:|']
for cls in ['conserved','ovary_specific','tage_specific','not_FDR_supported']:
    section.append(f'| {cls} | {int(p.program_class.eq(cls).sum())} | {int(s.program_class.eq(cls).sum())} |')
(out/'PATHWAY_ANALYSIS.md').write_text(report+'\n\n'+'\n'.join(section)+'\n')
status=json.loads((out/'run_status.json').read_text())
status['pathway_analysis']=dict(status='PASS',n_checks=len(checks),background_genes=len(u),enrichment_rows=len(e),program_counts={str(k):int(v) for k,v in p.program_class.value_counts().items()},full_family_classification_changes=len(changed))
if 'shared_ageing_pathways' not in status['figures']:status['figures'].append('shared_ageing_pathways')
(out/'run_status.json').write_text(json.dumps(status,indent=2))
files=[]
for folder in ['tables','figures','scripts','provenance','inputs']:
    for path in sorted((out/folder).rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts:files.append(dict(path=str(path.relative_to(out)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
for fn in ['README.md','PATHWAY_ANALYSIS.md','run_status.json']:
    files.append(dict(path=fn,sha256=hashlib.sha256((out/fn).read_bytes()).hexdigest()))
pd.DataFrame(files).to_csv(out/'output_checksums.tsv',sep='\t',index=False)
print(json.dumps(status['pathway_analysis'],indent=2))
