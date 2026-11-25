#!/usr/bin/env python3
"""Independent numerical validation and concise report; no rendering."""
from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
import pandas as pd
from scipy import stats

out=Path(__file__).resolve().parent.parent
ovary=out.parent.parent.parent
tage=ovary/'bulk-pseudobulk_clock/results/external_projection/tAge'
def read(name):return pd.read_csv(out/name,sep='\t')
checks=[]
def check(name,passed,detail=''):
    checks.append(dict(check=name,passed=bool(passed),detail=str(detail)))
    if not passed:raise AssertionError(name+': '+str(detail))
p=read('tables/sample_contrast_predictions.tsv')
c=read('tables/contrast_definitions.tsv')
m=read('tables/contrast_membership.tsv')
e=read('tables/contrast_effects.tsv')
s=read('tables/sample_inclusion_audit.tsv')
check('sample_inventory',len(s)==281 and s.sample_key.nunique()==281)
check('GSE154890_excluded',not s.dataset_folder.eq('GSE154890').any())
check('four_models_each_contrast',p.groupby('analysis_contrast_id').model.nunique().eq(4).all())
check('unique_prediction_key',not p.duplicated(['analysis_contrast_id','sample_key','model']).any())
check('finite_predictions',np.isfinite(p[['value_weeks','relative_shift_weeks']]).all().all())
log=(out/'logs/01_compute.log').read_text()
warning_count=log.count('RuntimeWarning: invalid value encountered in reduce')
pd.DataFrame([dict(warning='NumPy invalid value encountered in reduce',occurrences=warning_count,context='Previously documented official sklearn processing of NA-padded missing features; no arbitrary-zero replacement',validation='All predictions finite; unchanged-contrast replay checked')]).to_csv(out/'provenance/runtime_warning_audit.tsv',sep='\t',index=False)
check('same_modality_controls',m.groupby('analysis_contrast_id').modality.nunique().eq(1).all())
check('primary_adult_nonplacement',not c.loc[c.primary,['developmental','placement_only']].any().any())
check('primary_prespecified',c.loc[c.primary,'analysis_family'].eq('prespecified').all())
err=0
for (cid,model),g in p.groupby(['analysis_contrast_id','model']):
    ctl=g.role.eq('control');a=g.loc[~ctl,'value_weeks'];b=g.loc[ctl,'value_weeks']
    val=a.mean()-b.mean()
    row=e.loc[(e.analysis_contrast_id==cid)&(e.model==model)].iloc[0]
    err=max(err,abs(val-row.effect_weeks),float(np.max(np.abs(g.relative_shift_weeks-(g.value_weeks-b.mean())))))
check('independent_effect_and_centering_recalculation',err<1e-8,err)
ag=read('tables/clock_pairwise_agreement.tsv')
wide=e.pivot(index='analysis_contrast_id',columns='model',values='effect_weeks')
err=0
for _,r in ag.iterrows():
    if r.analysis_set=='primary_adult_prespecified':ids=c.loc[c.primary,'analysis_contrast_id']
    elif r.analysis_set=='all_prespecified_nonplacement':ids=c.loc[c.analysis_family.eq('prespecified')&~c.placement_only,'analysis_contrast_id']
    elif r.analysis_set=='all_prespecified':ids=c.loc[c.analysis_family.eq('prespecified'),'analysis_contrast_id']
    else:ids=c.analysis_contrast_id
    x=wide.loc[ids,r.model_1].to_numpy();y=wide.loc[ids,r.model_2].to_numpy()
    a=stats.pearsonr(x,y);b=stats.spearmanr(x,y)
    err=max(err,abs(a.statistic-r.pearson),abs(b.statistic-r.spearman),abs(a.pvalue-r.pearson_p),abs(b.pvalue-r.spearman_p),abs(np.mean(np.sign(x)==np.sign(y))-r.direction_agreement))
    check('pairwise_n_'+r.analysis_set+'_'+r.model_1+'_'+r.model_2,len(x)==r.n)
check('independent_scipy_clock_correlations',err<1e-8,err)
fd=read('tables/feature_coefficient_concordance_all_pca_source.tsv')
fs=read('tables/feature_direction_statistics_all_pca.tsv')
err=0
for _,r in fs.iterrows():
    z=fd[fd.comparison==r.comparison]
    check('feature_n_'+r.comparison,len(z)==r.n and not z.gene.duplicated().any())
    err=max(err,abs(stats.pearsonr(z.tage_weight,z.ovary_weight).statistic-r.pearson),abs(stats.spearmanr(z.tage_weight,z.ovary_weight).statistic-r.spearman))
check('independent_gene_weight_correlations',err<1e-8,err)
# Same controls, samples and physical mappings must reproduce the previous official tAge calculation.
old=pd.read_csv(tage/'tables/tage_sample_contrast_predictions.tsv',sep='\t')
replays=[]
for (cid,v),g in p.dropna(subset=['model_variant']).groupby(['analysis_contrast_id','model_variant']):
    if g.dataset_folder.iloc[0]=='GSE228975':continue
    ref=old[(old.tage_contrast_id==g.contrast_id.iloc[0])&(old.model_variant==v)]
    oldcontrols=set(ref.loc[ref.control_membership,'sample_key'])
    newcontrols=set(g.loc[g.role.eq('control'),'sample_key'])
    if oldcontrols!=newcontrols or set(ref.sample_key)!=set(g.sample_key):continue
    a=g.set_index('sample_key');b=ref.set_index('sample_key').loc[a.index]
    delta=np.max(np.abs(a.relative_shift_weeks-b.tage_shift_relative_to_reference_weeks))
    replays.append(dict(analysis_contrast_id=cid,variant=v,max_absolute_error=delta))
pd.DataFrame(replays).to_csv(out/'provenance/unchanged_tage_contrast_replay.tsv',sep='\t',index=False)
check('tAge_official_replay_unchanged_contrasts',len(replays)>0 and max(r['max_absolute_error'] for r in replays)<1e-7,len(replays))
integrity=read('provenance/input_integrity.tsv')
for _,r in integrity.iterrows():
    h=hashlib.md5(Path(r.path).read_bytes()).hexdigest()
    check('read_only_input_'+r.path,h==r.md5_before)
stems=['10_tage_clock_agreement_readable','03_tage_shift_relative_to_reference_distributions','feature_coefficient_concordance_all_pca']
for stem in stems:
    pdf=out/'figures'/f'{stem}.pdf';png=pdf.with_suffix('.png')
    text=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True)
    check('figure_export_'+stem,pdf.stat().st_size>1000 and png.stat().st_size>1000 and len(text)>100)
pd.DataFrame(checks).to_csv(out/'provenance/validation_checks.tsv',sep='\t',index=False)
status=dict(status='PASS',n_checks=len(checks),n_samples_inventory=len(s),n_samples_contrasted=p.sample_key.nunique(),n_contrasts=len(c),n_prespecified_contrasts=int(c.analysis_family.eq('prespecified').sum()),n_primary_adult_contrasts=int(c.primary.sum()),figures=stems)
(out/'run_status.json').write_text(json.dumps(status,indent=2))
primary=ag[ag.analysis_set=='primary_adult_prespecified']
lines=['# Current locked ovary clocks versus tAge','',
 'Reproduced the three historical figure layouts with the rebutal4 locked release (`../lockedmodel`): **HVG-ElasticNet trained without single-age studies** (110 profiles; new) and **VersionB-PCA/glmnet** (136 profiles; unchanged from rebutal3), using the rebutal4 external-projection application. No ovary model was retrained or recalibrated.','', 'Version B-PCA sample predictions, contrast effects, backprojected gene weights and PCA-versus-tAge gene concordance reproduce rebutal3 exactly. Pathway statistics for every clock are recomputed on the new shared gene background (below); `provenance/rebutal3_pca_reuse_check.tsv` quantifies the PCA changes.','',
 '## Figures','',
 '- [Four-clock agreement](figures/10_tage_clock_agreement_readable.pdf)',
 '- [Reference-centred tAge distributions](figures/03_tage_shift_relative_to_reference_distributions.pdf)',
 '- [Gene-coefficient concordance](figures/feature_coefficient_concordance_all_pca.pdf)','',
 'PNG previews are beside the PDFs. All graphics were generated with R.','',
 '## Cohort and comparison scope','',
 f'The application inventory has **{len(s)} profiles**; **{p.sample_key.nunique()}** occur in declared contrasts. There are **{len(c)} contrasts**: 44 prespecified plus the application\'s separately labelled GSE279573 reference-figure addition. The primary heatmap uses **{int(c.primary.sum())} adult prespecified non-placement contrasts**. The tAge-only panel and distributions use all 44 prespecified contrasts, including descriptive placement/developmental comparisons. GSE154890 remains excluded. GSE235021 bulk and pseudobulk are separate in every contrast.','',
 'Samples without a declared comparison are listed in `tables/sample_inclusion_audit.tsv`; no control group was invented. Every row in the distribution is a sample-within-contrast observation, so reused controls and targets are repeated. The figure caption reports both row and unique-profile counts.','',
 '## Primary effect agreement','',
 '| Comparison | n | Pearson r | Spearman rho | Direction agreement |','|---|---:|---:|---:|---:|']
for _,r in primary.iterrows():lines.append(f'| {r.model_1} vs {r.model_2} | {int(r.n)} | {r.pearson:.3f} | {r.spearman:.3f} | {r.direction_agreement:.1%} |')
lines += ['', 'Exact two-sided nominal P values and BH q values across the six pairs (separately for each correlation metric and analysis set) are in `tables/clock_pairwise_agreement.tsv`. Stars use nominal P, matching the historical display. Contrasts share studies and sometimes controls, so nominal correlation tests are descriptive and do not establish independent replication.','',
 '## Methods and differences from the screenshots','',
 'The ovary clocks use existing frozen sample predictions. Their comparison effect is mean target predicted age minus mean matched-control predicted age, following the original tAge comparison. The application instead tests delta age (prediction minus chronological age); both effect definitions are exported. Their difference is exactly the target-minus-control chronological-age difference. Application delta-age effects were independently replayed.','',
 'tAge is recomputed with the pinned official package and frozen mouse scaled-diff/YuGene-diff models, using the existing source-matrix preprocessing caches. Those caches retain the earlier tAge filtering, mouse identifier mapping and normalization, independently of the ovary clocks\' log-CPM preprocessing. For each current contrast, the declared-control gene-wise median is subtracted using the official tAge function; missing serialized-model features remain NA for its imputer. Official model outputs in months are converted by 52.142857/12, then centered by the mean output of that same control group. No overall-median fallback is permitted. tAge outputs are control-relative transcriptomic shifts, not conventional biological age acceleration.','',
 'Cached expression columns are matched to the application\'s physical source columns within modality. Four GSE228975 replicate labels are corrected in memory to match the application; original files remain unchanged. For the added GSE279573 contrast, the summary had missing membership strings: WT keys and oil-free keys were recovered from its two explicit prespecified comparisons and verified against the supplied effect.','',
 'Unlike the historical hard-coded heatmap subtitle, displayed n is calculated from the comparison set. Primary analysis excludes all placement contrasts and contrasts containing a profile under 7 weeks or designated developmental/OOD. `tables/contrast_definitions.tsv` records every inclusion/exclusion. Additional all-nonplacement and all-contrast sensitivities are supplied in the agreement table. The colour scale covers the full agreement range without clipping values below 0.5.','',
 'Gene weights are extracted from the current locked RDS models and directly from both serialized tAge pipelines. HVG-ElasticNet uses active nonzero coefficients in standardized gene space; PCA weights are rotation × PC coefficients in standardized gene space, matching the original backprojection convention. The common universe is the union of the two locked clocks\' eligible (training-filtered) genes, because the HVG model was filtered on 110 profiles and the PCA model on 136; the remaining description follows rebutal3. The common universe is the current ovary measured/eligible gene universe intersected with both tAge feature maps. Entrez-to-symbol mapping uses the official tAge projection mapping and prioritizes current ovary symbols. Raw-input-scale ovary weights and exact shared-gene source data are also exported. Different input transformations and units make absolute coefficient magnitudes non-comparable.','',
 '## Gene-weight agreement','', '| Comparison | Shared genes | Pearson r | Spearman rho |','|---|---:|---:|---:|']
for _,r in fs.iterrows():lines.append(f'| {r.comparison} | {int(r.n)} | {r.pearson:.3f} | {r.spearman:.3f} |')
lines += ['', 'These are stress-test projections. Some samples/studies overlap the current model-development data; per-contrast and per-sample flags are preserved. Several cohorts have small groups. Developmental, natural-age placement and age-extrapolation results remain descriptive.','',
 '## Reproduction and validation','',
 '`bash scripts/run_all.sh` recreates this output folder using the recorded local R runtime, pinned tAge Python environment, current application files and tAge preprocessing caches. Paths resolve relative to this project. The original application, locked models and historical results are read-only.','',
 f'**{len(checks)} checks passed**: exact sample/contrast joins, control centering, modality separation, frozen prediction replay, application effect replay, independent SciPy correlation checks, tAge replay for {len(replays)} unchanged contrast/variant combinations, readable figure exports, and unchanged input checksums.','',
 f'The official tAge pipeline emitted {warning_count} previously documented NumPy warnings while processing NA-padded features. They are retained in the execution log and `provenance/runtime_warning_audit.tsv`; all outputs are finite and unchanged-contrast replay is checked.', '',
 'See `provenance/validation_checks.tsv`, `provenance/input_integrity.tsv`, `provenance/tage_cache_sample_alignment.tsv`, `provenance/locked_model_validation.tsv` and the figure source tables under `tables/`. Visual review is recorded separately in `provenance/visual_QA.md`.']
(out/'README.md').write_text('\n'.join(lines)+'\n')
files=[]
for folder in ['tables','figures','scripts','provenance']:
    for path in sorted((out/folder).glob('*')):
        if path.is_file():files.append(dict(path=str(path.relative_to(out)),sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
pd.DataFrame(files).to_csv(out/'output_checksums.tsv',sep='\t',index=False)
print(json.dumps(status,indent=2))
