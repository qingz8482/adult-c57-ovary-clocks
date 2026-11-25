#!/usr/bin/env python3
"""Write evidence-grounded captions, audit coverage and verify final exports."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET

root = Path(__file__).resolve().parent.parent
dd = root / 'derived_data'

def read(name):
    with (dd / (name + '.tsv')).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows, fields):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

for r in read('input_files_used') + read('panel_f4_input_files_used'):
    assert hashlib.sha256(Path(r['input_file']).read_bytes()).hexdigest() == r['sha256'], r['input_file']
checks = read('numerical_validation') + read('independent_validation') + read('f4_validation')
assert all(r['passed'] == 'TRUE' for r in checks)
sample = read('panel_c_lopo_delta_mae_bootstrap_ci')
study = read('panel_c_lopo_study_bootstrap_ci')
slots = read('gsea_panel_slots')
back = max(float(r['absolute_discrepancy']) for r in read('gene_backprojection_validation'))
parts = max(float(r['absolute_discrepancy']) for r in read('sample_contribution_sum_validation'))
replay = max(float(r['max_prediction_error']) for r in read('conditional_oof_independent_replay'))
params = json.loads((root / '../../../locked_models/count_tuning_20261011/VersionB-PCA-glmnet/hyperparameters.json').read_text())
lam = float(params['lambda'])
cv = read('panel_a_cv_error_path')
chosen = min(cv, key=lambda r: abs(float(r['lambda']) - lam))
steroid = read('panel_f3_all_tested_steroid_terms')
dimensions = read('gsea_curve_dimensions')
curve_w = float(dimensions[0]['width_mm'])
curve_h = float(dimensions[0]['height_mm'])

caption = f'''**Supplementary Figure 2 | Model selection, component attribution and gene-level interpretation of the locked VersionB-PCA ovarian clock.**

**a,** Actual leave-one-study-out tuning curve for the winning 3,000-gene, 10-PC-input, alpha = 1 candidate. The line and ribbon show the mean of 14 held-out-study MAEs and ±1 study-level s.e.; the feature selector, gene scaling and PCA were fitted within each training fold in the original model-development run. The red line marks the deployed lambda ({lam:.4f}), chosen by the largest-lambda-within-one-s.e. rule for this candidate. The dashed line marks its minimum-MAE lambda. The upper axis reports nonzero coefficients in the all-136-sample deployment fit at each lambda, not fold-specific counts. These tuning errors informed model development and are not independent validation estimates.

**b,** Saved coefficient trajectories for all 10 PC inputs along the frozen model's regularization path. PC1, PC2, PC5 and PC8, which have nonzero deployed coefficients, are coloured by coefficient sign; the other PC paths are grey. The dashed line marks deployed lambda. Coefficients are in weeks per PC score unit. No age calibration is applied.

**c,** Conditional leave-one-active-PC-out diagnostic. The full four-active-PC baseline and each three-PC ablation were separately retuned and refitted using leave-one-study-out outer folds over 136 training profiles from 14 studies. Within each outer training set, alpha and lambda were selected by inner leave-one-study-out CV using the current six-alpha grid, study-balanced training weights, study-mean MAE and the one-s.e. lambda rule. Points show the sample-mean paired change in absolute error, MAE_without_PC minus MAE_full_refit. Bars are 95% percentile intervals from 5,000 paired sample bootstrap resamples (seed 20260630). The frozen deployment gene selection, scaling and PCA basis were held fixed, so this is a conditional component-ablation diagnostic, not fully nested prediction validation. Sample-bootstrap intervals do not account for study dependence; a separate panel and table report equally weighted study means with study bootstrap intervals. The refitted four-PC baseline is distinct from the original locked 10-input-PC fit.

**d,** Exact frozen contributions, coefficient × PC score, for all 136 training profiles. Profiles were ordered by chronological age, then sample key, and assigned deterministically to 20 nearly equal-count bins (six or seven profiles per bin). Each heatmap cell is a median PC contribution. Red increases predicted age; blue decreases it. At sample level, the frozen intercept plus the four PC contributions reproduces the deployed prediction (maximum discrepancy {parts:.2g} weeks). Medians of contributions need not sum to the median prediction. Grey tracks show independently standardized bin-median chronological and predicted ages; they are relative within-track scales, not comparable absolute ages.

**e,** Back-projected gene slopes in weeks per unit of input ln(1 + CPM): loading × PC coefficient divided by the frozen gene standard deviation. The Integrated column sums the four active-PC slopes. Rows are the deduplicated union of the strongest positive and negative per-PC and integrated weights, selected using the reference rule and capped at 40. Gene centring contributes to the intercept. Reconstruction uses actual training expression, not a truncated-PC reconstruction, and reproduces the PC-space predictor within {back:.2g} weeks.

**f,** Representative integrated-rank GSEA curves selected at within-rank BH FDR < 0.05: two older-aligned immune terms, one younger-aligned energy term and one younger-aligned lipid term. Candidates are ordered by FDR, decreasing absolute NES and term name; the two immune terms must have leading-edge Jaccard overlap <0.5. The GSEA ranking is the sum of loading × frozen PC coefficient on standardized genes, which differs from panel e's unstandardized-expression slopes. Curves use the current 3,000-gene ranking and saved mouse MSigDB 2026.1.Mm results for MH, M2:CP, M5:GO:BP and M7, restricted to 10–300 in-universe genes. The inherited weighted-KS method uses 500 gene-set permutations per set size; NES uses same-sign null normalization and empirical P uses all 500 permutations in its denominator. Insets give NES, within-rank BH FDR, in-universe set size (n) and leading-edge size (LE). Ticks mark gene-set members; shading marks the direction-appropriate leading-edge side. For negative ES, the saved peak rank identifies the hit immediately after the running-curve minimum. Pathway enrichment is conditional on the selected clock-gene universe and is not evidence of causal ageing effects.

**f2,** Two additional older- and two additional younger-aligned terms after excluding panel f, with FDR <0.05 and leading-edge Jaccard <0.5 against the previously displayed and newly selected terms.

**f3,** The requested steroid/cholesterol and cell-cycle categories, excluding terms in panels f/f2. No steroid/cholesterol term passes FDR <0.05 in the current model, so those two positions explicitly report the absence rather than substituting a nonsignificant term. The two displayed cell-cycle terms pass the threshold and have within-category leading-edge Jaccard <0.5.

**f4,** The top two integrated GSEA pathways within each of the nine functional groups displayed in diagnostic heatmap 13 (18 curves). Selection is restricted to the 70 displayed terms and follows the heatmap ranking rule: increasing q + P/1000, decreasing absolute NES, then pathway name. No exclusion of earlier panels, direction-balancing quota or additional redundancy filter is applied. All 18 chosen pathways have q <0.05 and nominal P <0.05. Functional labels are the heatmap's descriptive name-based assignments. Curves sharing genes are not independent biological observations. Group pairs are arranged in a 3-by-3 grid.

**Curve dimensions.** Each plotting area in panels f, f2, f3 and f4 is fixed to panel f's measured dimensions, {curve_w:.4f} × {curve_h:.4f} mm. Curve typography and line widths also match panel f. Panels f/f2/f3 have 78 × 108 mm canvases; f4 has a 234 × 185 mm canvas to accommodate 18 curves and group headings.

All expression-derived analyses use the same 136 model-development profiles. No external projection cohort or HVG comparison is included.
'''
(root / 'caption.md').write_text(caption)

lines = [
    '# Supplementary Figure 2: locked VersionB-PCA model interpretation', '',
    'Completed 2026-10-12 (Asia/Singapore). The analyses in the previous `supple2_PCA` figure package have been recomputed for the locked count-tuning release dated 2026-10-11: **3,000 genes, 136 profiles, 14 studies, 10 PC inputs and four active PCs (PC1, PC2, PC5, PC8)**.', '',
    '- [Assembled figure (PDF)](final/Supplementary_Figure_2_PCA_clock_model_interpretation.pdf)',
    '- [Assembled figure (PNG)](final/Supplementary_Figure_2_PCA_clock_model_interpretation.png)',
    '- [Full caption and methods](caption.md)',
    '- [Individual panels](panels/), including f2, f3, [f4 (nine functional groups)](panels/panel_f4.pdf) and the additional study-weighted ablation panel',
    '- [Source data and checks](derived_data/); [reference-output coverage](provenance/reference_output_coverage.tsv)', '',
    'The assembled figure and every panel are exported as PDF, SVG, 600-dpi PNG and TIFF. The assembled figure is 183 × 220 mm. PDF text is selectable; native R SVG exports retain vector graphics but can encode font glyphs as vector paths.', '',
    '## What was reproduced', '',
    '| Panel | Analysis | Current-model adaptation |',
    '|---|---|---|',
    '| a | Grouped-CV regularization error path | Uses the actual saved fold-local final-model tuning results for N=3,000, 10 PC inputs, alpha=1 |',
    '| b | Candidate-PC coefficient paths | Uses the saved weighted glmnet path; active PCs are 1, 2, 5 and 8 |',
    '| c | Retuned leave-one-PC-out ablation | 70 conditional outer refits, each with inner study-held-out tuning; no calibration |',
    '| d | Age-binned contribution heatmap | All 136 profiles, 20 age-ordered bins, four active-PC contributions |',
    '| e | Gene back-projection | Exact unstandardized-expression slopes and intercept verified against actual inputs |',
    '| f | Representative integrated GSEA curves | Current significant immune, energy and lipid categories; data-driven term labels |',
    '| f2 | Additional nonredundant GSEA curves | Excludes the actual current panel-f terms |',
    '| f3 | Steroid/cholesterol and cell-cycle curves | No significant steroid/cholesterol terms; two explicit absence panels and two qualifying cell-cycle curves |',
    '| f4 | Top two integrated pathways per heatmap function group | 18 curves in nine groups, selected from the displayed heatmap rows; original peaks, NES and FDR |', '',
    f'All f/f2/f3/f4 curve plotting areas measure {curve_w:.4f} × {curve_h:.4f} mm. Panels f/f2/f3 are 78 × 108 mm; f4 is 234 × 185 mm. The measured geometry is recorded in `derived_data/gsea_curve_dimensions.tsv`.', '',
    'The manuscript DOCX files accompanying the old folder were not analysis inputs and were not copied or edited. The earlier 117-sample/2,000-gene model, calibration and selected pathway claims were not reused. No HVG comparison was added.', '',
    '## Conditional ablation results', '',
    'Positive delta MAE means the refit predicted less accurately after omitting that PC. These results use the frozen deployment PCA basis; gene selection and PCA construction are not nested within the outer folds. Their bootstrap intervals describe the conditional paired error differences and do not include the uncertainty of re-developing the gene/PCA representation.', '',
    '| Omitted PC | Sample delta MAE, weeks (95% bootstrap interval) | Equal-study delta MAE, weeks (95% study-bootstrap interval) |',
    '|---|---:|---:|',
]
for row in sample:
    st = next(s for s in study if s['omitted_pc'] == row['omitted_pc'])
    def interval(r):
        return f"{float(r['delta_mae']):.2f} ({float(r['ci_low']):.2f}, {float(r['ci_high']):.2f})"
    lines.append(f"| {row['omitted_pc']} | {interval(row)} | {interval(st)} |")
lines += ['', 'PC2 and PC5 have positive intervals under both summaries. PC8 has a positive sample-bootstrap interval, but its study-bootstrap interval crosses zero. PC1 has a small negative point estimate with both intervals crossing zero. These are conditional diagnostics, not independent evidence of unique biological causation.', '',
          'The full ablation baseline retunes regression using the four active PCs; each omission retunes the remaining three. All 14 outer test studies are excluded from inner fitting, tuning and weight estimation. The original training normalization, gene scaling and rotation remain fixed for this conditional analysis. The six alpha candidates and 100 absolute lambda values match the current development recipe; unlike the earlier clock, no calibration is fitted.', '',
          '## Pathway selection', '',
          'Panel f4 selects the literal top two displayed terms within each heatmap function group, ranked by q + P/1000, decreasing absolute NES and name. Unlike f2/f3, f4 permits reuse of earlier pathways and does not apply a further redundancy filter. All selected terms have integrated-rank BH FDR <0.05. The 500-permutation GSEA results were reused from the current diagnostic analysis; enrichment curves, peaks and leading edges were independently reconstructed. The top categories are selected within the current clock-gene universe, with no reranking by raw gene slopes.', '',
          '| Panel | Category | Selected term or result |', '|---|---|---|']
for row in slots:
    lines.append(f"| {row['panel']} | {row['category']} | {row['gene_set_name'] or 'No terms pass FDR <0.05'} |")
lines += ['', f"All {len(steroid)} tested steroid/cholesterol terms are retained in `derived_data/panel_f3_all_tested_steroid_terms.tsv`, including nonsignificant results. The FDR threshold was not relaxed.", '',
          '## Validation and provenance', '',
          f"All **{len(checks)} numerical and independent checks** passed. All 70 saved outer predictions were independently replayed from expression-derived frozen scores and the saved tuning choices (maximum discrepancy {replay:.2g} weeks). All 70 alpha/lambda choices were recalculated from their inner tuning summaries. The 910 inner-fold memberships exclude the outer test study and keep inner train/test studies separate.", '',
          f"Prediction decomposition error: {parts:.2g} weeks. Gene back-projection error: {back:.2g} weeks. All displayed GSEA peaks and enrichment scores match the frozen tables, the running sums close at zero, and selection/redundancy thresholds pass. Input checksums cover the locked model, expression, tuning summary, fitting code and current GSEA tables.", '',
          'The sample-unit bootstrap is retained to match the reference panel. The additional study-unit bootstrap resamples the 14 paired study-mean differences, giving each study equal weight; both use 5,000 resamples and seed 20260630. Confidence intervals are not multiplicity-adjusted. Bin heatmaps use medians for display; exact sample-wise additivity is validated before aggregation.', '',
          '## Reproduction', '',
          'Run `bash scripts/run_all.sh` from this folder. The scripts write only inside this supplementary-analysis folder. To refresh figures without recomputing ablation fits, run `03_plot.R` using the existing Rscript environment recorded in `run_all.sh`. Use `03_plot.R --gsea-only` to export only f/f2/f3/f4. To recompute f4 selection and curves from a revised heatmap, first run `02_gsea_curves.R`, then the plot script, `python3 scripts/06_validate_f4.py` and `python3 scripts/05_report.py`. The reference scripts are preserved under `provenance/reference_scripts/` for audit, not executed against the current model.', '']
(root / 'README.md').write_text('\n'.join(lines))

# The reference package's figure/analysis files all have current counterparts.
ref = root.parents[5] / 'bulk-pseudobulk_clock/results/supple2_PCA/figures and legend/supplementary/pca_clock_model_interpretation'
coverage = []
for old in sorted(ref.rglob('*')):
    if not old.is_file():
        continue
    rel = old.relative_to(ref)
    if rel.parts[0] == 'scripts':
        new = Path('scripts/01_compute.R; scripts/02_gsea_curves.R; scripts/03_plot.R')
        status = 'adapted_into_current_pipeline'
        exists = True
    else:
        new = rel
        exists = (root / rel).is_file()
        status = 'recomputed_for_locked_VersionB_PCA'
    coverage.append(dict(reference=str(rel), output=str(new), exists=exists, status=status))
assert coverage and all(r['exists'] for r in coverage)
write(root / 'provenance/reference_output_coverage.tsv', coverage, ['reference', 'output', 'exists', 'status'])

pdf_checks = []
for path in sorted([*(root / 'panels').glob('*.pdf'), *(root / 'final').glob('*.pdf')]):
    tree = ET.fromstring(subprocess.check_output(['pdftotext', '-bbox', str(path), '-']))
    ns = {'x': 'http://www.w3.org/1999/xhtml'}
    pages = tree.findall('.//x:page', ns)
    count, outside = 0, 0
    for page in pages:
        width, height = float(page.get('width')), float(page.get('height'))
        for word in page.findall('.//x:word', ns):
            count += 1
            outside += (float(word.get('xMin')) < 0 or float(word.get('yMin')) < 0 or
                        float(word.get('xMax')) > width + 1 or float(word.get('yMax')) > height + 1)
    assert len(pages) == 1 and count > 0 and outside == 0, path
    pdf_checks.append(dict(file=str(path.relative_to(root)), pages=len(pages), selectable_words=count, words_outside_page=outside))
write(root / 'provenance/pdf_export_checks.tsv', pdf_checks, ['file', 'pages', 'selectable_words', 'words_outside_page'])
manifest = []
for path in sorted(root.rglob('*')):
    if path.is_file() and path.name != 'output_manifest.tsv':
        manifest.append(dict(path=str(path.relative_to(root)), bytes=path.stat().st_size, sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
write(root / 'provenance/output_manifest.tsv', manifest, ['path', 'bytes', 'sha256'])
print(json.dumps(dict(checks_passed=len(checks), reference_files_covered=len(coverage), pdfs_checked=len(pdf_checks), selected_GSEA_terms=sum(r['status']=='selected' for r in slots)), indent=2))
