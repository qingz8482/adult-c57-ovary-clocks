#!/usr/bin/env python3
"""Audit output coverage and write the local analysis index without plotting."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import xml.etree.ElementTree as ET

root = Path(__file__).resolve().parent.parent
development = root.parent.parent
clock_root = development.parents[2]
reference = clock_root / 'bulk-pseudobulk_clock/results/integrated_diagnostics'
folders = ['10_pca_diagnostics', '12_PCA deployment model GSEA']

def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows, fields):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        w.writeheader()
        w.writerows(rows)

for row in read(root / 'provenance/input_checksums.tsv'):
    assert hashlib.sha256(Path(row['path']).read_bytes()).hexdigest() == row['sha256'], row['path']

coverage = []
for folder in folders:
    for old in sorted((reference / folder).iterdir()):
        if not old.is_file():
            continue
        new = root / folder / old.name
        coverage.append(dict(reference=f'{folder}/{old.name}', output=f'{folder}/{old.name}',
                             exists=new.is_file(), bytes=new.stat().st_size if new.is_file() else 0,
                             status='recomputed_for_locked_VersionB_PCA'))
assert coverage and all(row['exists'] and row['bytes'] > 0 for row in coverage)
write(root / 'provenance/reference_output_coverage.tsv', coverage,
      ['reference', 'output', 'exists', 'bytes', 'status'])

checks = read(root / 'provenance/pca_validation.tsv') + read(root / 'provenance/gsea_validation.tsv')
assert all(row['passed'] == 'TRUE' for row in checks)
heatmap_checks = read(root / 'provenance/heatmap_revision_validation.tsv')
assert all(row['passed'].upper() == 'TRUE' for row in heatmap_checks)
pdf_checks = []
for folder in folders:
    for path in sorted((root / folder).glob('*.pdf')):
        tree = ET.fromstring(subprocess.check_output(['pdftotext', '-bbox', str(path), '-']))
        ns = {'x': 'http://www.w3.org/1999/xhtml'}
        pages = tree.findall('.//x:page', ns)
        count, clipped = 0, 0
        for page in pages:
            width, height = float(page.attrib['width']), float(page.attrib['height'])
            for word in page.findall('.//x:word', ns):
                count += 1
                a = word.attrib
                clipped += (float(a['xMin']) < 0 or float(a['yMin']) < 0 or
                            float(a['xMax']) > width + 1 or float(a['yMax']) > height + 1)
        assert len(pages) == 1 and count > 0 and clipped == 0, path
        pdf_checks.append(dict(file=str(path.relative_to(root)), pages=len(pages),
                               selectable_text_words=count, out_of_bounds_words=clipped))
write(root / 'provenance/pdf_export_checks.tsv', pdf_checks,
      ['file', 'pages', 'selectable_text_words', 'out_of_bounds_words'])
pc = read(root / folders[0] / 'final_pca_clock_all_input_pc_summary.tsv')
gsea = read(root / folders[1] / 'final_pca_pc_gsea_all_results.tsv')
integrated = read(root / folders[1] / 'final_pca_integrated_gsea_all_results.tsv')
active = [row['pc'] for row in pc if row['active_in_clock'] == 'TRUE']
summary = []
for name in active + ['integrated_PCA_clock']:
    rows = integrated if name == 'integrated_PCA_clock' else [r for r in gsea if r['pc'] == name]
    sig = [r for r in rows if float(r['empirical_p']) < .05 and float(r['fdr_within_pc']) < .1]
    summary.append(dict(rank=name, tested=len(rows), older=sum(r['direction'] == 'older_aligned' for r in sig),
                        younger=sum(r['direction'] == 'younger_aligned' for r in sig)))
write(root / 'provenance/result_counts.tsv', summary, ['rank', 'tested', 'older', 'younger'])

lines = [
    '# Locked VersionB-PCA diagnostics and deployment-model GSEA', '',
    'Completed 2026-10-12 (Asia/Singapore). Recomputed from the frozen count-tuning release dated 2026-10-11: **3,000 genes, 136 training profiles from 14 studies, 10 PC inputs, and four active PCs (PC1, PC2, PC5, PC8)**. The selected model was taken from the current deployment registry.', '',
    'PCA diagnostics, the complete deployment GSEA analysis, and the supplementary model-interpretation package are included. No HVG-versus-Version-B comparison was run. These are descriptive interpretations of the all-training deployment fit, not new validation estimates.', '',
    '## Outputs', '',
    '- [Supplementary Figure 2 analyses](supple2_PCA/README.md): grouped-CV and coefficient paths, conditional PC ablation, contribution and gene-weight heatmaps, and GSEA panels f/f2/f3/f4.',
    '- [PCA diagnostics](10_pca_diagnostics/): six PDF/PNG figure pairs covering variance, scores, coefficients, age associations and top 20 genes per active PC; complete loadings, selector scores, gene counts and prediction-replay tables.',
    '- [Deployment-model GSEA](12_PCA%20deployment%20model%20GSEA/): all 13 reference figure types, PC and integrated ranks, all tested results, older/younger terms, significant terms, leading-edge genes, collection summaries, permutation QC and the PC-plus-integrated heatmap source.',
    '- [GSEA result summary](12_PCA%20deployment%20model%20GSEA/final_pca_pc_gsea_summary.md).',
    '- [Functional GSEA review and revised heatmaps](12_PCA%20deployment%20model%20GSEA/functional_heatmap_review.md): top five terms in each direction of every PC, integrated functional representatives, and per-cell significance symbols without rank numbers. Exact ranks remain in the companion tables. Energy, lipid, other metabolism, RNA/protein processing, adhesion, neural and signalling groups are included when supported; full assignments and selection reasons are exported.',
    '- [Reference-file coverage](provenance/reference_output_coverage.tsv), [PCA validation](provenance/pca_validation.tsv), [GSEA validation](provenance/gsea_validation.tsv), [running-sum replay](provenance/gsea_running_sum_replay.tsv), scripts, logs and frozen mouse gene-set snapshots.', '',
    f'All **{len(coverage)} reference files** have corresponding newly computed outputs. The old benchmark/deployment rows are replaced by the current locked deployment scope; earlier model numerical results and biological interpretations are not reused.', '',
    '## Scope and ranking', '',
    'All 10 input PCs appear in the PCA diagnostics. PC3, PC4, PC6, PC7, PC9 and PC10 have exactly zero clock coefficients; coefficient-oriented GSEA is undefined for their all-zero rankings, so GSEA uses only the four active PCs. Primary/secondary displays show PC1/PC2 and PC5/PC8, respectively.', '',
    'PC gene rank = frozen loading × frozen PC coefficient. Integrated rank = sum of these contributions across active PCs, equal to the exported coefficient on standardized gene expression. It is not the coefficient per raw ln1p(CPM) unit. Positive ranks increase predicted age; positive/negative NES denotes enrichment toward the older/younger predictive side. No calibration is applied.', '',
    'Variance fractions use total variance across all 3,000 standardized input genes. Cumulative variance includes intervening PCs with zero regression coefficients. Per-PC gene counts are dense loadings, not independent selected-gene sets. Training age correlations are unweighted descriptive correlations.', '',
    '## GSEA method', '',
    'Mouse MSigDB 2026.1.Mm is loaded from saved native-mouse cache snapshots. Collections match the prior workflow: MH Hallmark, M2 canonical pathways, M5 GO biological processes and M7. Gene sets are intersected with the 3,000 selected clock genes; primary sizes are 10–300 genes. The max-500 sensitivity output reports eligible-set counts, as in the reference, rather than a separate GSEA run.', '',
    'The reference custom weighted-KS procedure is retained: weight exponent 1, seed 20260630, 500 random gene-set permutations per rank and set size. NES uses same-sign null normalization. Empirical P = (1 + null scores at least as extreme in the observed direction)/(501), using all permutations in the denominator. This preserves the reference calculation; it is not fgsea or phenotype-permutation GSEA. No refined permutation pass was run.', '',
    'BH q values are computed separately within each PC (across its tested collections), within each PC/collection, and for the integrated rank. The reported exploratory threshold is q < 0.10 and nominal P < 0.05. P values have coarse Monte Carlo resolution (minimum 1/501); q values do not control a joint family across all PCs. Top-ranked plots can include nonsignificant terms, marked separately. Pathway labels characterize the model within its selected-gene universe, not causal age effects or independent biological validation.', '',
    '| Rank | Tested pathways | Older-aligned at threshold | Younger-aligned at threshold |',
    '|---|---:|---:|---:|',
]
lines += [f"| {r['rank']} | {r['tested']} | {r['older']} | {r['younger']} |" for r in summary]
lines += ['', '## Verification and reproduction', '',
          f'The functional heatmap revision additionally passed {len(heatmap_checks)} checks for unchanged numerical results, complete PC top-term coverage, correct significance labels, representative-set redundancy and traceable assignments. Refresh only these figures with `Rscript scripts/06_refresh_gsea_figures.R --heatmaps-only` in the scRNA_r4 environment, then run `python3 scripts/07_validate_heatmap.py` and this report script.', '',
          f'The original PCA/GSEA suite passed all {len(checks)} numerical and scope checks; the supplementary package has its own validation record. Frozen predictions were replayed through PC scores and original-gene coefficients; integrated ranks reproduce the exported standardized coefficients. BH corrections were independently recalculated, and 500 pathway/rank examples were checked using a dense running-sum implementation, including leading-edge membership. Input checksums were verified before and after the analysis.', '',
          'All 19 PDFs were checked for selectable text and text extending beyond the page boundary. R-generated contact sheets in provenance support visual review. Run `bash scripts/run_all.sh` from this diagnostic folder to reproduce. The scripts use the existing R environment, saved model and prepared training inputs, and local gene-set snapshots. All generated outputs stay in this folder. PDF figures retain vector text; PNG previews use 300 dpi. These are report-sized diagnostic figures rather than assembled final manuscript panels.', '']
(root / 'README.md').write_text('\n'.join(lines))
manifest = []
for path in sorted(root.rglob('*')):
    if path.is_file() and path.name != 'output_manifest.tsv':
        manifest.append(dict(path=str(path.relative_to(root)), bytes=path.stat().st_size,
                             sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
write(root / 'provenance/output_manifest.tsv', manifest, ['path', 'bytes', 'sha256'])
print(json.dumps(dict(reference_outputs=len(coverage), checks_passed=len(checks), result_counts=summary), indent=2))
