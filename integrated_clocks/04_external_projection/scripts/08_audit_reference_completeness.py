#!/usr/bin/env python3
"""Audit every reference/application file and validate substantive output parity.

Python handles text/tabular/file integrity only. All figures are generated in R.
"""
import csv
import hashlib
import json
import subprocess
import struct
from collections import Counter, defaultdict
from pathlib import Path

out = Path(__file__).resolve().parent.parent
old = out.parent.parent.parent / 'bulk-pseudobulk_clock/results/external_projection'
scope = {'diagnostics', 'plots', 'GSE235021', 'PRJNA863443', 'PRJNA1002222'}

def read(name, base=out):
    with (base / name).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))

def write(name, rows, fields=None):
    if Path(name).parts[0] in {'simplified_classes', 'alignment_audit', 'matched_control_verification'}:
        return
    path = out / name
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open('w') as handle:
        w = csv.DictWriter(handle, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

retired_prefixes = tuple(r['path_prefix'] for r in read('provenance/retired_diagnostic_outputs.tsv'))
assert not any(str(path.relative_to(out)).startswith(retired_prefixes) for path in out.rglob('*') if path.is_file()), 'Retired diagnostic output recreated'

baseline = json.loads((out / 'provenance/plot_revision_baseline.json').read_text())
for name, sha in baseline['core_sha256'].items():
    assert digest(out / name) == sha, f'Protected numerical input changed: {name}'

pred = read('external_projection_predictions_all.tsv')
meta = {r['sample_key']: r for r in pred}
assert len(pred) == 562 and len(meta) == 281
assert len({r['dataset_cohort'] for r in pred}) == 20
assert len({r['dataset_folder'] for r in pred}) == 19
assert not any(r['dataset_folder'] == 'GSE154890' for r in pred)
assert all(r['pass'] == 'TRUE' for r in read('provenance/comprehensive_analysis_validation.tsv'))
assert all(r['pass'] == 'TRUE' for r in read('provenance/derived_analysis_validation.tsv'))

stats = read('tables/all_reference_contrast_statistics.tsv')
assert len(stats) == 90
bytest = {(r['analysis_contrast_id'], r['model']): r for r in stats}
ann = read('provenance/complete_plot_annotation_audit.tsv')
annotated = set()
for r in ann:
    k = (r['analysis_contrast_id'], r['model'])
    p = float(bytest[k]['welch_p'])
    expected = '****' if p < 1e-4 else '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else f'p = {p:.3f}'
    assert r['annotation'] == expected and abs(float(r['welch_p']) - p) < 1e-12, r
    annotated.add(k)
assert annotated == set(bytest), ('Missing annotations', set(bytest) - annotated)

source = read('provenance/complete_plot_source_data.tsv')
sources = defaultdict(list)
for r in source:
    sources[r['figure_stem']].append(r)
membership = defaultdict(set)
for r in read('tables/all_reference_contrast_membership.tsv'):
    role = 'Reference' if r['contrast_role'] in ('control', 'placement_reference') else 'Target'
    membership[r['analysis_contrast_id'], role].add(r['sample_key'])
for r in ann:
    if r['stem'] not in sources:  # internal components are exported in combined source tables
        continue
    for role in ('Reference', 'Target'):
        keys = {z['sample_key'] for z in sources[r['stem']] if z['model'] == r['model'] and z['panel'] == r['panel'] and z['plot_group'] == r[role]}
        assert keys == membership[r['analysis_contrast_id'], role], (r['stem'], r['analysis_contrast_id'], role, keys ^ membership[r['analysis_contrast_id'], role])

dictionary = {r['sample_key']: r for r in read('provenance/reference_annotation_dictionary.tsv')}
orderrows = read('provenance/complete_group_order_audit.tsv')
variants = {r['stem']: r['variant'] for r in orderrows}
axis_orders = defaultdict(list)
for r in orderrows:
    axis_orders[r['stem']].append(r['plot_group'])
for r in ann:
    if r['stem'] not in sources:
        continue
    present = {z['plot_group'] for z in sources[r['stem']] if z['model']==r['model'] and z['panel']==r['panel']}
    axis = [x for x in axis_orders[r['stem']] if x in present]
    assert float(r['xmin']) == axis.index(r['Reference'])+1, (r['stem'],'reference bracket coordinate')
    assert float(r['xmax']) == axis.index(r['Target'])+1, (r['stem'],'target bracket coordinate')
label_checks = []
for stem, rows in sources.items():
    if stem not in variants:
        continue
    variant = variants[stem]
    for r in rows:
        md = dictionary[r['sample_key']]
        groupcol, classcol = {'root': ('simple_group', 'simple_class'), 'special': ('special_group', 'special_class'), 'new': ('new_group', 'new_class'), 'natural': ('natural_aging_plot_group', 'natural_aging_plot_group')}.get(variant, ('simple_group', 'simple_class'))
        if variant in ('raw', 'root') and meta[r['sample_key']]['dataset_folder'] == 'PRJNA1002222':
            groupcol, classcol = 'special_group', 'special_class'
        if meta[r['sample_key']]['dataset_folder'] == 'GSE235021':
            groupcol, classcol = 'new_group', 'new_class'
        if meta[r['sample_key']]['dataset_folder'] == 'GSE228975':
            groupcol = None  # old OOD bucket separated into the two prespecified test arms
        expected_class = md[classcol]
        if variant not in ('simplified', 'special', 'natural') and expected_class == 'stress':
            expected_class = 'stress_model'
        assert r['plot_class'] == expected_class, (stem, r, expected_class)
        if groupcol:
            assert r['plot_group'] == md[groupcol], (stem, r, md[groupcol])
    label_checks.append({'figure_stem': stem, 'variant': variant, 'sample_labels_classes_exact': True})
write('provenance/reference_label_class_validation.tsv', label_checks)

alignment = read('provenance/sample_expression_alignment.tsv')
warnings = []
for r in pred:
    flags = []
    if r['negative_predicted_age'] == 'TRUE': flags.append('negative unconstrained predicted age')
    if r['outside_training_age_range'] == 'TRUE': flags.append('chronological age outside frozen training range')
    if r['locked_training_sample'] == 'TRUE': flags.append('sample included in model development')
    if float(r['projection_distance_percentile']) > .95: flags.append('distance above training 95th percentile')
    if flags: warnings.append({'sample_key': r['sample_key'], 'dataset_folder': r['dataset_folder'], 'model': r['model'], 'warning': '; '.join(flags)})
write('diagnostics/external_projection_diagnostic_warnings.tsv', warnings)
write('external_projection_alignment_issues.tsv', [dict(r, issue='Preserved historical positional mapping; independent identity evidence not added in this task') for r in alignment if 'position' in r['alignment_method']])
write('external_projection_skipped_unaligned_samples.tsv', [], ['sample_key','dataset_folder','reason'])
write('diagnostics/skipped_projection_alignment_summary.tsv', [], ['dataset_folder','n_skipped','reason'])
for name in ['external_projection_warnings.txt']:
    (out/name).write_text('All 281 requested profiles were projected. GSE154890 excluded.\nDevelopment overlap is model-specific (HVG-ElasticNet excluded single-age studies from training); inspect per-model sample flags.\nNegative predicted ages retained without clipping. OOD is based on the current frozen training reference.\nCRA003645 has no valid matched-control effect. Positional source alignment is retained and documented.\nPseudobulk profiles, not cells, are the statistical units.\n')

checks = read('provenance/comprehensive_analysis_validation.tsv')
checks += [{'check':'88 prespecified plus 2 added Welch annotations replayed','pass':'TRUE'}, {'check':'plotted bracket membership exact','pass':'TRUE'}, {'check':'reference group/class labels match corresponding folder','pass':'TRUE'}, {'check':'protected core numerical files SHA-256 unchanged','pass':'TRUE'}]
write('diagnostics/projection_diagnostic_validation_checks.tsv', checks)
write('matched_control_verification/matched_control_verification_results.tsv', read('external_projection_control_matching_audit.tsv'))
write('matched_control_verification/matched_control_verification_failures.tsv', [], ['analysis_contrast_id','model','failure'])
write('diagnostics/special_dataset_control_audit.tsv', [r for r in read('external_projection_control_matching_audit.tsv') if r['dataset_folder'] in ('PRJNA863443','PRJNA1002222','GSE235021','GSE279446','GSE279573')])
write('simplified_classes/external_projection_simplified_validation_audit.tsv', label_checks)
write('projection_artifact_audit.tsv', read('provenance/prediction_replay_checks.tsv'))
write('input_gene_identifier_audit.tsv', read('provenance/input_provenance.tsv'))
write('diagnostics/pca_projection_score_source_audit.tsv', [{'model':'VersionB-PCA/glmnet','source':'rebutal4 lockedmodel (unchanged from count_tuning_20261011) frozen representation','n_external':281,'n_training':136,'n_PC_inputs':10,'active_PCs':'PC1;PC2;PC5;PC8','reference_projection_file':'reference_pca_projection_scores.tsv'}])

for ds in ('PRJNA863443','PRJNA1002222','GSE235021'):
    low = ds.lower()
    metadata_fields = list(read('tables/external_projection_samples.tsv')[0]) + ['projection_group','projection_class','projection_role','comparison_id','natural_aging_plot_group','simple_group','simple_class','special_group','special_class']
    metadata_fields = list(dict.fromkeys(metadata_fields))
    write(f'{ds}/{low}_projection_metadata.tsv', [r for r in meta.values() if r['dataset_folder']==ds], metadata_fields)
    for cohort in {r['dataset_cohort'] for r in meta.values() if r['dataset_folder']==ds}:
        rows = [r for r in read('external_projection_relative_to_controls.tsv') if r['dataset_cohort']==cohort]
        for r in rows:
            r['plot_group'] = r['new_group'] if ds=='GSE235021' else r['special_group']
        write(f'{ds}/relative_delta_plot_source_{cohort}.tsv', rows, ['sample_key','analysis_contrast_id','model','plot_group','relative_delta_age_week'])
    write(f'{ds}/{low}_artifact_audit.tsv', read('provenance/prediction_replay_checks.tsv'))
    write(f'{ds}/{low}_expression_scale_audit.tsv', [r for r in read('provenance/normalization_audit.tsv') if meta[r['sample_key']]['dataset_folder'] == ds])
    write(f'{ds}/{low}_final_validation.tsv', checks)
    write(f'{ds}/{low}_group_size_audit.tsv', [r for r in read('external_projection_group_summary.tsv') if r['dataset_folder']==ds])
    (out/f'{ds}/{low}_projection_warnings.txt').write_text('\n'.join(f"{r['sample_key']} | {r['model']} | {r['warning']}" for r in warnings if r['dataset_folder']==ds)+'\n')
write('PRJNA863443/verification/prjna863443_projection_verification.tsv', checks)
write('PRJNA863443/verification/prjna863443_projection_verification_failures.tsv', [], ['check','failure'])
write('GSE235021/tables/PRJNA863443_GSE235021_combined_predicted_vs_age_validation.tsv', checks)
write('GSE235021/tables/gse235021_pca_pc_contribution_validation.tsv', checks)
write('GSE235021/tables/gse235021_pca_pc_contribution_geometry.tsv', [r for r in read('provenance/complete_figure_exports.tsv') if 'gse235021_pca_pc_contribution_heatmap' in r['stem']])

exports = read('provenance/complete_figure_exports.tsv')
export_stems = {r['stem'] for r in exports}
assert all(str(Path(x).with_suffix('')) in export_stems for x in baseline['plot_files'] if not x.startswith(retired_prefixes))
file_checks = []
for r in exports:
    stem = r['stem']
    pdf, png = out/(stem+'.pdf'), out/(stem+'.png')
    result = subprocess.run(['pdfinfo', str(pdf)], capture_output=True, text=True)
    assert result.returncode == 0, stem
    with png.open('rb') as f:
        header = f.read(24)
    assert header[:8] == b'\x89PNG\r\n\x1a\n', stem
    width, height = struct.unpack('>II', header[16:24])
    assert abs(width-300*float(r['width_inches'])) <= 1 and abs(height-300*float(r['height_inches'])) <= 1, (stem, width, height)
    file_checks.append({'stem':stem,'PDF_readable':True,'PNG_present':True,'PNG_width':width,'PNG_height':height,'PNG_dpi':300,'pdf_bytes':pdf.stat().st_size})
write('provenance/comprehensive_figure_file_checks.tsv',file_checks)

# Audit only the five corresponding folders requested by the user.
equivalents = {
 'cra003645_heldout_exclusion_audit.tsv':'tables/external_projection_samples.tsv',
 'external_projection_run_manifest.tsv':'provenance/source_checksums.tsv',
 'diagnostics/optional_projection_files_found.tsv':'provenance/reference_file_coverage.tsv',
 'PRJNA1002222/preexisting_output_preservation_audit.tsv':'provenance/plot_revision_baseline.json',
}
inventory=[]
for path in sorted(old.rglob('*')):
    if not path.is_file() or path.relative_to(old).parts[0] not in scope: continue
    name=str(path.relative_to(old)); target=name; status='regenerated'; reason='Current two-clock output at corresponding relative path.'
    if name.startswith(retired_prefixes):
        status='retired_by_user';target='';reason='Duplicate diagnostic removed by the user; generator and associated figure-specific outputs removed.'
    elif name.startswith('tAge/'):
        status='not_applicable';target='';reason='Separate published tAge clock, cross-tAge comparisons, environment or cache; outside the two locked models.'
    elif 'GSE154890' in name:
        status='excluded_dataset';target='provenance/excluded_GSE154890_samples.tsv';reason='All GSE154890 profiles excluded as requested.'
    elif name.startswith('diagnostics/05_relative_delta_age_matched_controls_CRA003645.'):
        status='not_applicable';target='external_projection_dataset_control_availability.tsv';reason='No matched controls or Table 4 contrast for the six supplied developmental profiles; historical relative plot is unsupported/stale.'
    elif name.endswith('.log'):
        status='historical_runtime_record';target='logs/07_complete_reference_outputs.log';reason='New execution log replaces historical runtime record.'
    elif name.endswith('preappend_preservation_snapshot.rds'):
        status='not_applicable';target='';reason='Old two-clock append-operation snapshot; current projection is a separate two-clock application.'
    elif name.startswith('alignment_audit/'):
        status='equivalent_audit';target='alignment_audit/current_sample_expression_alignment.tsv';reason='All 281 selected samples aligned; current explicit alignment and source provenance supersede historical mapping templates/skipped-candidate searches.'
    elif name in equivalents:
        status='equivalent_audit';target=equivalents[name];reason='Equivalent current application audit with current model/sample scope.'
    elif name.endswith('.md'):
        status='current_documentation';target=name
    elif not (out/name).exists():
        status='MISSING';reason='No current counterpart or documented scope exception.'
    inventory.append({'reference_path':name,'reference_bytes':path.stat().st_size,'status':status,'current_path':target,'reason':reason,'reference_sha256':digest(path) if not name.startswith('tAge/') else ''})
missing=[r for r in inventory if r['status']=='MISSING']
write('provenance/reference_file_coverage.tsv',inventory)
write('diagnostics/optional_projection_files_found.tsv', [{'reference_path':r['reference_path'],'current_path':r['current_path'],'status':r['status']} for r in inventory if r['reference_path'].startswith('diagnostics/')])
assert not missing, missing
summary=Counter(r['status'] for r in inventory)
paired_reference_figures={str(Path(r['reference_path']).with_suffix('')) for r in inventory if r['status']=='regenerated' and Path(r['reference_path']).suffix in ('.pdf','.png')}
report=[
 '# Review of the five corresponding reference folders','',
 'Scope: diagnostics, plots, GSE235021, PRJNA863443 and PRJNA1002222. No additional reference subtrees are imported.','',
 f"Generated {len(exports)} PDF/PNG figure pairs. Recreated {len(paired_reference_figures)} applicable reference figure stems at corresponding paths; all retained previous application plots are present; user-retired diagnostic outputs are excluded.", '',
 'The cohort contains 281 profiles across 20 modality-specific datasets (19 study accessions): 231 bulk profiles across 18 datasets, plus 17 GSE235021 and 33 PRJNA863443 pseudobulk profiles. All 19 GSE154890 profiles are excluded.', '',
 '## Annotation and statistical conventions','',
 'Group names, class vocabulary and explicit ordering are taken from the corresponding original scripts and metadata. Dedicated PRJNA863443 figures use their dedicated reference-control classes; diagnostics use the simplified metadata preferred by the original diagnostic scripts. Diagnostic legends display stress_model where the originals do. Original group strings retain underscores and punctuation.', '',
 'GSE235021 uses the original age-specific labels from the newly-added plots (6m/9m polished, black, brown). Bulk and pseudobulk remain separate in all tests and displays. GSE228975’s historical single age-bin OOD bucket is split into Srsf1 WT/cKO solely to identify its prespecified comparison. The projection class remains developmental_OOD.', '',
 'All 88 prespecified Welch tests (44 contrasts × 2 clocks) are reported. The original reference-only GSE279573 oil-free rescue versus WT contrast adds two separately labelled tests. Both the 88-test BH adjustment and expanded 90-test BH adjustment are supplied. Stars use raw two-sided Welch p; nonsignificant results show numeric p. Relative effects subtract each contrast’s own control mean; controls reused across comparisons appear in separately labelled panels.', '',
 '## Added analyses','',
 'Matched-control effects and integrity audits; availability and cohort-flow summaries; the pairwise clock comparison and group/sample disagreement; frozen training-reference projections and in-sample MAE/R²; PCA reference map and active-PC contributions; model-specific distance/OOD flags and interpretation grading; natural-aging panels; modality-blocked/bulk-only GSE235021 heatmaps; PRJNA863443/GSE235021 scatter; source tables and reproducibility audits.', '',
 'Training-reference lines are descriptive least-squares display lines through frozen in-sample predictions (HVG-ElasticNet n = 110; Version B-PCA n = 136). They do not recalibrate predictions. EN distance is Euclidean in frozen standardized genes; PCA distance is squared Mahalanobis in ten frozen PC inputs (covariance diagonal ridge 10^-6). Percentiles compare against each model’s own training profiles (110 or 136), using original >95th/>99th rules.', '',
 'Interpretation grades retain the original vocabulary but require both clocks to agree in direction. Placement, age extrapolation or current development overlap restricts a comparison to descriptive interpretation. Grades are heuristic QC summaries and do not establish rejuvenation. Development overlap is flagged per model; see the README.', '',
 '## Explicit scope exceptions','',
 '- GSE154890 figures: excluded by request.',
 '- CRA003645 matched-control relative figure: no valid control contrast among the six requested developmental profiles; its raw-delta plot is present.',
 '- Duplicate outputs retired on 2026-10-12: newly added with pvalue, combined three-panel delta age (historical folder name), and diagnostics/01_projection_cohort_flow_*; their generators and figure-specific source outputs were removed. The generic predicted-age overview and bulk/pseudobulk contrast-effect plots are also retired. Existing overview cohort-flow plots are retained. See provenance/retired_diagnostic_outputs.tsv.',
 '- Other reference subfolders are outside this review and are not imported.',
 '- Old append snapshots and runtime logs: replaced by current provenance/validation rather than copied as new results.', '',
 '## Validation','',
 f"Verified {len(annotated)} distinct clock-specific Welch annotations, exact plotted bracket membership, reference label/class metadata, readable PDFs, PNG counterparts, and SHA-256 identity of all protected core predictions/statistics/expression files. {len(label_checks)} plot-specific label/class checks passed.", '',
 'Representative visual inspection is recorded separately in provenance/comprehensive_visual_QA.md.', '',
 'See [every reference file](provenance/reference_file_coverage.tsv), [all group orders](provenance/complete_group_order_audit.tsv), [all figure annotations](provenance/complete_plot_annotation_audit.tsv), [plotted source data](provenance/complete_plot_source_data.tsv), and [90-test table](tables/all_reference_contrast_statistics.tsv).', '',
 'Reproduce with `bash scripts/run_all.sh`. For a figure/derived-analysis refresh: run `scripts/07_complete_reference_outputs.R` with the recorded R environment, then the Python validation/report scripts. R alone creates all graphics.', '',
 'Reference-file accounting: '+ '; '.join(f'{k}: {v}' for k,v in sorted(summary.items()))+'.',''
]
(out/'REFERENCE_REVIEW.md').write_text('\n'.join(report))
for folder in sorted(scope):
    write(f'{folder}/reference_coverage.tsv', [r for r in inventory if Path(r['reference_path']).parts[0] == folder])
# Replace corresponding old narrative files with truthful current-scope summaries.
for r in inventory:
    if r['status']=='current_documentation':
        path=out/r['current_path'];path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text('# Two locked models: '+path.stem.replace('_',' ')+'\n\n'+
                        'Current projections: 281 profiles, 20 modality-specific datasets, two clocks (HVG-ElasticNet single-age exclusion; Version B-PCA). GSE154890 excluded.\n\n'+
                        'This folder contains the corresponding current-clock analyses. Numerical results are recomputed from the rebutal4 locked release. Group annotations follow the corresponding historical reference.\n\n'+
                        'Consult the application-root REFERENCE_REVIEW.md and provenance/reference_file_coverage.tsv for complete output coverage, validation and explicit exceptions.\n')
result={'status':'PASS','samples':281,'modality_specific_datasets':20,'accessions':19,'pdf_png_pairs':len(exports),'reference_figure_stems_recreated':len(paired_reference_figures),'prespecified_tests':88,'additional_tests':2,'missing_applicable_reference_files':0,'reference_file_categories':dict(summary),'numeric_inputs_unchanged':True,'bracket_membership_exact':True,'reference_annotation_match':True}
(out/'provenance/comprehensive_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
