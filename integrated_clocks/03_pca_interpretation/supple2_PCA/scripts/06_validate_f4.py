#!/usr/bin/env python3
"""Verify f4 selection, replay its running sums and check physical curve sizes."""
from pathlib import Path
from collections import defaultdict, Counter
import csv
import hashlib
import math

root = Path(__file__).resolve().parent.parent
dd = root / 'derived_data'
gd = root.parent / '12_PCA deployment model GSEA'


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


checks = []


def check(name, passed, detail):
    checks.append(dict(check=name, passed='TRUE' if passed else 'FALSE', detail=detail))
    assert passed, (name, detail)


source = read(gd / 'heatmap_term_selection_audit.tsv')
selected = read(dd / 'panel_f4_selected_terms.tsv')
source_by_group = defaultdict(list)
for row in source:
    source_by_group[row['biology_group']].append(row)
expected = []
for group, rows in source_by_group.items():
    rows.sort(key=lambda r: (float(r['integrated_q']) + float(r['integrated_p']) / 1000,
                            -abs(float(r['integrated_NES'])), r['gene_set_name']))
    expected.extend((group, r['gene_set_id']) for r in rows[:2])
check('two_pathways_per_displayed_function_group', len(selected) == 18 and len(source_by_group) == 9
      and Counter(r['category'] for r in selected) == Counter({g: 2 for g in source_by_group}),
      '18 curves, two per each of the nine heatmap groups.')
check('literal_top_two_integrated_terms', [(r['category'], r['gene_set_id']) for r in selected] == expected,
      'Candidates restricted to displayed heatmap rows; original q + P/1000, decreasing |NES|, term-name order. No exclusions.')
check('all_f4_terms_significant', all(float(r['fdr_within_pc']) < .05 and float(r['empirical_p']) < .05 for r in selected),
      'All selected integrated results meet q<0.05 and P<0.05.')
inputs = read(dd / 'input_files_used.tsv') + read(dd / 'panel_f4_input_files_used.tsv')
check('input_checksums_unchanged', all(hashlib.sha256(Path(r['input_file']).read_bytes()).hexdigest() == r['sha256'] for r in inputs),
      'Existing model/GSEA inputs and the f4 heatmap selection source verified.')
all_results = {r['gene_set_id']: r for r in read(gd / 'final_pca_integrated_gsea_all_results.tsv')}
check('f4_statistics_match_integrated_results', all(math.isclose(float(r[k]), float(all_results[r['gene_set_id']][k]),
      rel_tol=1e-12, abs_tol=1e-12) for r in selected for k in ['NES', 'ES', 'empirical_p', 'fdr_within_pc', 'peak_rank']),
      'Reported statistics copied from the original integrated GSEA, without new permutations.')
rank_rows = sorted(read(gd / 'final_pca_integrated_gsea_rank_vector.tsv'), key=lambda r: -float(r['rank_score']))
genes = [r['gene_id'] for r in rank_rows]
rank_values = [float(r['rank_score']) for r in rank_rows]
universe_sets = {r['gene_set_id']: set(r['genes'].split(';')) for r in read(gd / 'gene_set_universe_audit.tsv')}
curves = defaultdict(list)
for row in read(dd / 'panel_f4_enrichment_curves.tsv'):
    curves[row['gene_set_name']].append(row)
errors, es_errors, leading_matches = [], [], []
for term in selected:
    rows = sorted(curves[term['gene_set_name']], key=lambda r: int(r['rank']))
    assert len(rows) == 3000 and [r['gene_id'] for r in rows] == genes
    members = universe_sets[term['gene_set_id']]
    hits = [i for i, gene in enumerate(genes) if gene in members]
    hit_set = set(hits)
    assert len(hits) == int(term['gene_set_size_in_universe'])
    weight_sum = sum(abs(rank_values[i]) for i in hits)
    running, value = [], 0.0
    for i in range(len(genes)):
        value += abs(rank_values[i]) / weight_sum if i in hit_set else -1 / (len(genes) - len(hits))
        running.append(value)
    errors.append(max(abs(v - float(row['running_es'])) for v, row in zip(running, rows)))
    assert all((row['hit'] == 'TRUE') == (i in hit_set) for i, row in enumerate(rows))
    maximum, minimum = max(running), min(running)
    positive = abs(maximum) >= abs(minimum)
    es = maximum if positive else minimum
    peak = running.index(maximum) + 1 if positive else min(i + 1 for i in hits if i > running.index(minimum))
    assert peak == int(term['peak_rank'])
    es_errors.append(abs(es - float(term['ES'])))
    leading = {genes[i] for i in hits if (i + 1 <= peak if positive else i + 1 >= peak)}
    leading_matches.append(leading == set(term['leading_edge_genes'].split(';')))
    assert abs(running[-1]) < 1e-10
check('independent_running_sum_and_membership_replay', max(errors) < 1e-12,
      f'Maximum curve discrepancy across 54,000 rank positions: {max(errors):.3g}.')
check('independent_ES_peak_and_leading_edge_replay', max(es_errors) < 1e-12 and all(leading_matches),
      f'All 18 peak ranks and leading edges match; maximum ES discrepancy {max(es_errors):.3g}.')
dims = read(dd / 'gsea_curve_dimensions.tsv')
reference = [r for r in dims if r['panel'] == 'panel_f']
rw, rh = float(reference[0]['width_mm']), float(reference[0]['height_mm'])
check('all_curve_areas_match_panel_f', Counter(r['panel'] for r in dims) ==
      Counter(dict(panel_f=4, panel_f2=4, panel_f3=4, panel_f4=18, panel_f5=2)) and
      all(abs(float(r['width_mm']) - rw) < 1e-6 and abs(float(r['height_mm']) - rh) < 1e-6 for r in dims),
      f'All 32 plotting areas, including f3 empty slots, measure {rw:.6f} x {rh:.6f} mm on the PDF device.')
check('f2_f3_canvas_matches_panel_f', all(float(r['canvas_width_mm']) == 78 and float(r['canvas_height_mm']) == 108
      for r in dims if r['panel'] in ['panel_f', 'panel_f2', 'panel_f3']), 'f/f2/f3: 78 x 108 mm; f4: 234 x 185 mm.')
check('f5_canvas_one_row_of_panel_f', all(float(r['canvas_width_mm']) == 78 and float(r['canvas_height_mm']) == 54
      for r in dims if r['panel'] == 'panel_f5'), 'f5: 78 x 54 mm (one row of the panel-f grid).')
with (dd / 'f4_validation.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=['check', 'passed', 'detail'], delimiter='\t')
    writer.writeheader()
    writer.writerows(checks)
print(f'Passed {len(checks)} f4 and curve-size checks; 18 curves across nine groups.')
