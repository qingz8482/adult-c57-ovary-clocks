#!/usr/bin/env python3
"""Validate revised pathway comparisons independently; no rendering."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
import pandas as pd
from scipy import stats
out=Path(__file__).resolve().parent.parent
def read(p):return pd.read_csv(out/p,sep='\t')
checks=[]
def check(name,ok,detail=''):
    checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name+': '+str(detail))
st=read('tables/model_resolved_pathway_statistics.tsv')
p=read('tables/model_resolved_pathway_comparison.tsv')
sel=read('tables/model_resolved_selected_pathways.tsv')
w=read('tables/model_resolved_gene_weights.tsv').set_index('gene');models=list(w.columns)
members=read('tables/pathway_library_membership.tsv');members['pathway_id']=members.database+' | '+members.pathway
sets={k:set(g.gene) for k,g in members.groupby('pathway_id')}
check('complete_model_analysis_grid',len(st)==len(p)*8 and not st.duplicated(['pathway_id','model','analysis']).any())
check('effect_bounds',st.rank_biserial.between(-1,1).all())
check('permutation_count',st.permutations.eq(100000).all())
check('finite_monte_carlo_probabilities',np.max(abs(st.p_value-(st.permutation_exceedances+1)/(st.permutations+1)))<1e-12)
def bh(vals):
    v=np.asarray(vals);order=np.argsort(v);q=np.minimum.accumulate((v[order]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1];res=np.empty(len(v));res[order]=np.minimum(q,1);return res
err=0
for _,g in st.groupby(['model','analysis']):err=max(err,np.max(abs(bh(g.p_value)-g.FDR)))
check('BH_all_pathways_per_model_and_analysis',err<1e-10,err)
# Independent SciPy U/effect/asymptotic diagnostics for all representatives and a random audit sample.
rng=np.random.default_rng(20261012);auditids=set(sel.pathway_id)|set(rng.choice(p.pathway_id,size=80,replace=False))
err=0;n=0
for r in st[st.pathway_id.isin(auditids)].itertuples():
    a=w[r.model].to_numpy();mask=w.index.isin(sets[r.pathway_id])
    if r.analysis=='absolute':a=abs(a)
    test=stats.mannwhitneyu(a[mask],a[~mask],alternative='greater' if r.analysis=='absolute' else 'two-sided',method='asymptotic',use_continuity=True)
    eff=2*test.statistic/(mask.sum()*(~mask).sum())-1
    err=max(err,abs(eff-r.rank_biserial),abs(test.pvalue-r.asymptotic_p_value));n+=1
check('independent_scipy_rank_statistics',err<1e-9,f'{n} rows; max error {err}')
e=st[st.analysis=='absolute'].pivot(index='pathway_id',columns='model',values='rank_biserial').reindex(columns=models)
q=st[st.analysis=='absolute'].pivot(index='pathway_id',columns='model',values='FDR').reindex(columns=models)
pp=p.set_index('pathway_id').loc[e.index]
diff=e.iloc[:,:2].mean(axis=1)-e.iloc[:,2:].mean(axis=1)
check('descriptive_group_difference',np.max(abs(diff-pp.difference))<1e-12)
check('no_invented_difference_test',not any('difference_p' in c or 'difference_ci' in c for c in p.columns))
os=((q.iloc[:,:2]<.05)&(e.iloc[:,:2]>0)).sum(axis=1);ts=((q.iloc[:,2:]<.05)&(e.iloc[:,2:]>0)).sum(axis=1)
check('individual_model_support_counts',np.array_equal(os,pp.ovary_model_support) and np.array_equal(ts,pp.tage_model_support))
classes=np.where((os>0)&(ts>0),'Shared support',np.where(os>0,'Ovary support only',np.where(ts>0,'tAge support only','No FDR support')))
check('group_classification_transparent',np.array_equal(classes,pp.support_class))
check('all_four_support',np.array_equal((os+ts)==4,pp.all_four_supported))
for r in sel.itertuples():
    if r.display_group=='Ovary-preferential candidates':check('ovary_preference_'+str(r.row_index),r.difference>0 and r.ovary_model_support>0 and r.tage_model_support==0)
    elif r.display_group=='tAge-preferential candidates':check('tage_preference_'+str(r.row_index),r.difference<0 and r.tage_model_support>0 and r.ovary_model_support==0)
maxj=0
for i,a in enumerate(sel.pathway_id):
    for b in sel.pathway_id.iloc[:i]:maxj=max(maxj,len(sets[a]&sets[b])/len(sets[a]|sets[b]))
check('nonredundant_representatives',maxj<.35,maxj)
selection=read('tables/model_resolved_representative_selection.tsv')
check('selection_audit_matches_figure',set(selection.loc[selection.selected,'pathway_id'])==set(sel.pathway_id))
check('maximum_five_per_group',sel.groupby('display_group').size().le(5).all())
check('representatives_have_at_least_15_genes',sel.pathway_size.ge(15).all())
check('exact_null_engine_validation',read('provenance/permutation_engine_exact_null_check.tsv').absolute_error.max()<.02)
for r in read('provenance/model_resolved_input_integrity.tsv').itertuples():check('preserved_'+Path(r.path).name,hashlib.md5(Path(r.path).read_bytes()).hexdigest()==r.md5_before)
stems=['pathway_model_resolved_comparison','pathway_signed_weight_alignment','pathway_all_programs_comparison']
for stem in stems:
    pdf=out/'figures'/f'{stem}.pdf';txt=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True)
    check('figure_'+stem,len(txt)>200 and pdf.with_suffix('.png').stat().st_size>1000)
pd.DataFrame(checks).to_csv(out/'provenance/model_resolved_validation_checks.tsv',sep='\t',index=False)
coverage=read('tables/model_resolved_gene_coverage.tsv')
legacy=read('tables/pathway_program_classification.tsv')
legacy['ovary_PCA_only']=legacy.ovary_support.fillna('').str.contains('ovary_pca') & ~legacy.ovary_support.fillna('').str.contains('ovary_elasticnet')
legacy[['database','pathway','program_class','ovary_support','tage_support','ovary_PCA_only']].to_csv(out/'tables/historical_pathway_model_attribution.tsv',sep='\t',index=False)
legacy_ov=legacy[legacy.program_class=='ovary_specific'];legacy_shared=legacy[legacy.program_class=='conserved']
counts=p.support_class.value_counts()
all4=int(p.all_four_supported.sum());same=int((p.all_four_signed_same_direction & p.all_four_signed_supported).sum())
rows=['# A more defensible view of shared and preferential pathways','',
 '## Recommended figures','',
 '- [Four-clock importance matrix and direct descriptive contrast](figures/pathway_model_resolved_comparison.pdf): recommended main figure.',
 '- [Signed age-weight alignment](figures/pathway_signed_weight_alignment.pdf): the same pathways, with positive and negative weight-rank effects.',
 '- [All-pathway comparison](figures/pathway_all_programs_comparison.pdf): complete landscape; numbers identify the displayed representatives.','',
 '## Why the original figure is insufficient','',
 'The historical plot mixes ORA and absolute-rank tests, chooses the smallest FDR across multiple analyses, and suppresses disagreement between the two ovarian clocks. Its group-specific score formulas differ. Two separate top-ten lists do not show how the same pathway behaves in the other clocks. A significance/non-significance split is not a test of between-clock difference or tissue specificity. Related GO/Reactome/KEGG terms can also repeat overlapping gene sets.','',
 f'**The hidden model dependence is substantial: {int(legacy_ov.ovary_PCA_only.sum())} of {len(legacy_ov)} historical ovary-only pathways had ovarian support from PCA analyses alone. Likewise, {int(legacy_shared.ovary_PCA_only.sum())} of {len(legacy_shared)} historical conserved pathways had ovarian support only from PCA.** These pooled labels therefore did not establish agreement between HVG-ElasticNet and PCA. Exact attribution is exported in `tables/historical_pathway_model_attribution.tsv`.','',
 '## Revised analysis','',
 f'All four models use the same {len(w):,}-gene background and {len(p):,} pathways. Every pathway is assessed using the same absolute-weight rank effect, **2 × AUC − 1**, where AUC is the probability that a pathway gene has a larger weight than a non-pathway gene (ties count one half). Positive effects indicate enrichment of large absolute weights, zero indicates no rank separation, and negative effects indicate smaller weights. This is model importance relative to each clock\'s own gene background; it is not explained variance or an expression/activity score.','',
 'Ranks are computed within each model, so absolute coefficient units are not compared across models. Genes outside the model have weight zero. Feature transformations, regularization, training cohorts and feature-selection rules still differ; this design cannot separate their effects from tissue biology. Sparsity affects the resulting importance distribution and must be read alongside the number of active genes:','',
 '| Model | Nonzero genes in common background |','|---|---:|']
for r in coverage.itertuples():rows.append(f'| {r.model} | {r.active_genes_in_background:,} |')
rows+=['',
 '**Significance:** 100,000 uniform permutations of gene labels, shared across all models and absolute/signed score columns, preserve each model\'s weight distribution, cross-model gene-weight relationships and pathway sizes/overlaps. One-sided tests assess large absolute weights; two-sided tests assess signed weight-rank shifts. P=(exceedances+1)/(100,001), with minimum attainable P about 0.00001. BH is applied across all 6,592 pathways/libraries per model and analysis type. Ties are handled by average ranks and the empirical null; this avoids relying on the normal approximation with the sparse HVG model. Asymptotic tie-corrected results are retained only as diagnostics.','',
 'These are competitive random-gene null tests of fixed coefficients. Gene exchangeability is an assumption; correlated gene modules, model fitting/selection uncertainty and biological sample variability are not represented. Stars are exploratory support, not biological replication. Near-threshold calls and weak/sparse-model calls can be limited by finite permutation resolution.','',
 '**Difference column:** mean importance of the two ovary clocks minus mean importance of the two tAge variants. This is a descriptive contrast. The two algorithms are not biological replicates; no differential P value or confidence interval is manufactured. A positive value means greater pathway rank separation in ovary models under this metric, not ovarian tissue specificity.','',
 '**Direction figure:** the same rank effect computed from signed coefficients. Red/positive means pathway weights shift toward positive age coefficients relative to non-pathway genes; blue/negative means the reverse. Mixed positive/negative genes can cancel at pathway level. Coefficients are conditional model associations, not measured expression changes, causal effects or pathway activation.','',
 '**Displayed pathways:** up to five per group, restricted to pathways with at least 15 mapped background genes so tiny sets do not dominate effect-size selection. All 5–500-gene pathways remain tested and exported. Shared support requires at least one supported model in each group and is ranked by the smaller of the two group maximum effects. Ovary-preferential candidates require ovary-only support and a positive group-mean difference, ranked by that difference; tAge candidates use the reverse. Greedy Jaccard filtering excludes pathways overlapping any already selected representative by >=0.35. All candidates, ordering and exclusions are exported. This reduces display repetition; it does not alter testing or prove independent biological programs. No pathway is selected manually for its biological name.','',
 '## What the data support','',
 '| Support category | Pathways |','|---|---:|']
for cls in ['Shared support','Ovary support only','tAge support only','No FDR support']:rows.append(f'| {cls} | {int(counts.get(cls,0))} |')
rows += ['',f'**{all4} pathways have significant absolute-weight enrichment in all four clocks.** The broader shared-support category allows support from one ovary model and one tAge variant; it must not be called four-clock consensus.',
 f'Across all pathways, **{same} have significant signed-weight shifts in all four models with the same sign**. This count is independent of absolute-weight support, and still describes coefficient alignment rather than biological conservation.','',
 'Selected pathways make the model dependence visible:','',
 '| Group | Pathway | Ovary models supported / 2 | tAge variants supported / 2 | Ovary minus tAge effect |','|---|---|---:|---:|---:|']
for r in sel.itertuples():rows.append(f'| {r.display_group} | {r.pathway} | {r.ovary_model_support} | {r.tage_model_support} | {r.difference:.3f} |')
rows+=['','## Interpretation for the manuscript','',
 'Prefer “shared pathway weighting” and “ovary-preferential model programs” over “conserved ageing mechanisms” and “ovary-specific pathways”. Stronger biological conservation would require replicated age-associated expression or pathway activity in matching directions. Tissue specificity would require ovarian versus non-ovarian age effects, a direct age-by-tissue comparison, and independent replication. The present locked feature weights alone cannot establish those claims.','',
 '## Sources and reproducibility','',
 'The interpretation of Mann–Whitney rank statistics and tie handling follows the [SciPy reference](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html). The use of shared genes to recognize redundant pathway annotations is motivated by [Merico et al., Enrichment Map (2010)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0013984); the greedy Jaccard display rule here is our specified implementation, not the Enrichment Map algorithm.','',
 f'**{len(checks)} validation checks passed.** Independent SciPy tests verified U/effect statistics and asymptotic diagnostics for {n} rows, including every representative and a random pathway audit. The C++ permutation engine was compared against fully enumerated toy nulls with ties and mixed signs. All FDRs, group differences, support rules, display-overlap constraints and prior-file integrity were verified.','',
 'Run `scripts/07_model_resolved_pathways.R`, `scripts/08_validate_permutation_engine.R`, then `scripts/09_validate_model_resolved.py` using the recorded R/Python environments. They are appended to `scripts/run_all.sh`. The first revised step caches permutations under `cache/`, invalidated by input/source hashes. All graphics are R-generated. Detailed results use the `model_resolved_` prefix under `tables/` and `provenance/`. Previous figures and their original-method analyses are preserved.']
(out/'MODEL_RESOLVED_PATHWAYS.md').write_text('\n'.join(rows)+'\n')
readme=(out/'README.md').read_text();marker='## Recommended model-resolved pathway comparison'
readme=readme.split(marker)[0].rstrip()
(out/'README.md').write_text(readme+'\n\n'+marker+'\n\nUse the [four-clock comparison](figures/pathway_model_resolved_comparison.pdf), [signed-weight view](figures/pathway_signed_weight_alignment.pdf), and [all-pathway overview](figures/pathway_all_programs_comparison.pdf) for interpreting agreement and differences. The historical minimum-FDR plot remains a reproduction reference. [Methods, limitations and findings](MODEL_RESOLVED_PATHWAYS.md).\n')
status=json.loads((out/'run_status.json').read_text());status['model_resolved_pathways']=dict(status='PASS',n_checks=len(checks),pathway_counts={str(k):int(v) for k,v in counts.items()},all_four_absolute_supported=all4,all_four_signed_supported_same_direction=same,permutations=100000)
for stem in stems:
    if stem not in status['figures']:status['figures'].append(stem)
(out/'run_status.json').write_text(json.dumps(status,indent=2))
files=[]
for folder in ['tables','figures','scripts','provenance','inputs']:
    for path in sorted((out/folder).rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts:files.append(dict(path=str(path.relative_to(out)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
for fn in ['README.md','PATHWAY_ANALYSIS.md','MODEL_RESOLVED_PATHWAYS.md','run_status.json']:files.append(dict(path=fn,sha256=hashlib.sha256((out/fn).read_bytes()).hexdigest()))
pd.DataFrame(files).to_csv(out/'output_checksums.tsv',sep='\t',index=False)
print(json.dumps(status['model_resolved_pathways'],indent=2))
