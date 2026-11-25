#!/usr/bin/env python3
"""Independently check heatmap values/coverage and report functional evidence."""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import math

root = Path(__file__).resolve().parent.parent
out = root / '12_PCA deployment model GSEA'
pcs = ['PC1', 'PC2', 'PC5', 'PC8']


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def truth(value):
    return str(value).upper() == 'TRUE'


def number(row, key):
    return float(row[key])


def qpass(row, threshold):
    return number(row, 'fdr_within_pc') < threshold and number(row, 'empirical_p') < .05


checks = []


def check(name, condition, detail):
    checks.append({'check': name, 'passed': bool(condition), 'detail': detail})
    assert condition, name


snapshot = root / 'provenance/heatmap_revision_numeric_checksums.tsv'
check('saved_GSEA_inputs_unchanged', all(hashlib.sha256(Path(r['path']).read_bytes()).hexdigest() == r['sha256']
      for r in read(snapshot)), 'Both all-result tables and the tested gene-set universe retain their original hashes.')
raw = read(out / 'final_pca_pc_gsea_all_results.tsv') + read(out / 'final_pca_integrated_gsea_all_results.tsv')
lookup = {(r['pc'], r['gene_set_id']): r for r in raw}
selection = read(out / 'heatmap_term_selection_audit.tsv')
selected = {r['gene_set_id'] for r in selection}
source = read(out / 'final_pca_pc_plus_integrated_gsea_heatmap_source.tsv')
source_lookup = {(r['pc'], r['gene_set_id']): r for r in source}
groups = defaultdict(list)
for row in raw:
    groups[row['pc'], row['direction']].append(row)
ranks = {}
anchors = set()
for (pc, direction), rows in groups.items():
    rows.sort(key=lambda r: (number(r, 'fdr_within_pc') + number(r, 'empirical_p') / 1000,
                            -abs(number(r, 'NES')), r['gene_set_name']))
    for rank, row in enumerate(rows, 1):
        ranks[pc, row['gene_set_id']] = rank
        if pc in pcs and rank <= 5:
            anchors.add((pc, row['gene_set_id']))
check('all_PC_top_five_terms_preserved', len(anchors) == 40 and anchors <= set(source_lookup),
      'Five older- and five younger-aligned terms for each of four active PCs; no global truncation.')
check('complete_heatmap_matrix', len(source) == len(source_lookup) == 5 * len(selected)
      and set(source_lookup) == set(itertools.product(pcs + ['integrated_PCA_clock'], selected)),
      f'{len(selected)} rows by five ranks; every value is measured, with no missing cells or zero imputation.')
check('NES_P_and_q_match_saved_results', all(math.isclose(number(row, field), number(lookup[key], field),
      rel_tol=1e-12, abs_tol=1e-12) for key, row in source_lookup.items()
      for field in ['NES', 'empirical_p', 'fdr_within_pc']), 'Every displayed number agrees with the saved result.')
check('directional_ranks_match_full_results', all(int(row['directional_rank']) == ranks[key]
      for key, row in source_lookup.items()), 'Directional ranks use all tested pathways, before display selection.')


def expected_label(key, row):
    return '*' if qpass(row, .05) else '+' if qpass(row, .10) else ''


check('significance_only_labels_correct', all(row['cell_label'] == expected_label(key, row)
      and truth(row['pc_top_term']) == (key in anchors)
      and truth(row['significant_q05']) == qpass(row, .05)
      and truth(row['significant_q10']) == qpass(row, .10) for key, row in source_lookup.items()),
      '* means q<0.05; + means 0.05<=q<0.10; both require P<0.05. Tile labels contain no rank numbers.')
representatives = [r for r in selection if truth(r['integrated_function_representative'])]
check('integrated_representatives_significant', all(qpass(lookup['integrated_PCA_clock', r['gene_set_id']], .05)
      for r in representatives), f'All {len(representatives)} integrated representatives pass q<0.05 and P<0.05.')
check('every_selected_row_has_a_reason', selected == {key[1] for key in anchors} | {r['gene_set_id'] for r in representatives},
      'All rows are required PC top-five terms or significant integrated functional representatives.')
rep_groups = defaultdict(list)
for row in representatives:
    rep_groups[row['biology_group']].append(row['gene_set_id'])
check('maximum_four_representatives_per_group', all(len(v) <= 4 for v in rep_groups.values()), str(dict(Counter(r['biology_group'] for r in representatives))))
membership = {r['gene_set_id']: set(r['genes'].split(';')) for r in read(out / 'gene_set_universe_audit.tsv')}
overlaps = [len(membership[a] & membership[b]) / len(membership[a] | membership[b])
            for ids in rep_groups.values() for a, b in itertools.combinations(ids, 2)]
check('representatives_nonredundant_within_group', all(v < .5 for v in overlaps),
      f'Maximum pairwise Jaccard within a representative group: {max(overlaps):.4f}; PC anchors are exempt.')
check('requested_metabolic_groups_present', all(name in rep_groups for name in ['Mitochondrial / energy', 'Fatty acid / lipid']),
      'Both requested metabolic groups have significant integrated representatives.')
audit = read(out / 'heatmap_all_terms_function_audit.tsv')
check('all_tested_pathways_have_audited_assignments', len(audit) == len(raw) and {(r['pc'], r['gene_set_id']) for r in audit} == set(lookup),
      f'All {len(audit)} pathway/rank combinations have a function assignment, rank and significance flag.')
check('cytokinesis_and_respiratory_burst_classified_correctly',
      all(r['biology_group'] == 'Cell cycle / DNA' for r in audit if 'CYTOKINESIS' in r['gene_set_name'])
      and all(r['biology_group'] == 'Immune / inflammatory' for r in audit if 'RESPIRATORY_BURST' in r['gene_set_name']),
      'Cytokinesis is cell division; respiratory burst is not assigned to mitochondrial energy.')
top8 = read(out / 'top_terms_per_pc_by_function.tsv')
check('full_top_eight_per_PC_exported', len(top8) == 64 and all(ranks[r['pc'], r['gene_set_id']] == int(r['directional_rank']) <= 8 for r in top8),
      'All eight older- and eight younger-aligned terms per PC retain their groups and exact q values in the companion table.')
with (root / 'provenance/heatmap_revision_validation.tsv').open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=['check', 'passed', 'detail'], delimiter='\t')
    writer.writeheader()
    writer.writerows(checks)

order = list(dict.fromkeys(r['biology_group'] for r in selection))
integrated_audit = [r for r in audit if r['pc'] == 'integrated_PCA_clock']
summary = ['# Functional review of integrated and PC-specific GSEA', '',
    'Revised 2026-10-12. The saved GSEA results are unchanged. The previous union of top PC terms was truncated to 45 rows globally; this concealed metabolic functions and removed some leading PC5/PC8 terms.', '',
    '## Revised figure', '',
    '[Combined PC and integrated heatmap](13_gsea_nes_heatmap_top_terms_with_integrated_total.pdf) · [PC-only heatmap](05_gsea_nes_heatmap_top_terms.pdf)', '',
    f'The combined heatmap contains **{len(selected)} terms in {len(order)} functional groups**. All top five terms in each direction of each active PC are retained. Tiles display significance symbols only; rank numbers are omitted at the user\'s request. The full top eight per PC/direction remain in plots 03/04 and in [the grouped top-term table](top_terms_per_pc_by_function.tsv).', '',
    'Additional rows are up to four integrated representatives per named functional group at **q < 0.05 and P < 0.05**. Candidates are ordered by the existing FDR/P/absolute-NES rule. Up to two representatives per direction are selected first, then unused places are filled from the remaining candidates. Within a group, representatives must have pairwise Jaccard < 0.5 using the tested 3,000-gene universe. Required PC top terms are never removed by this redundancy filter. There is no global row cap.', '',
    '`*` marks q < 0.05; `+` marks 0.05 ≤ q < 0.10, both with P < 0.05. An unmarked tile does not meet q < 0.10 and P < 0.05. Tile colour always shows the saved NES, including nonsignificant values. q values are the original BH correction over all tested pathways within each rank, not recalculated on selected rows.', '',
    'The integrated NES is computed from the summed PC gene contributions, not an average of PC NES values. Positive and negative NES describe the older- and younger-predictive sides of the clock. They do not establish activation or repression with chronological age. In particular, UP/DN suffixes in immune-signature names refer to their source experiments.', '',
    'Functional groups are descriptive assignments based on pathway names, with explicit rule precedence; M7 immune signatures are assigned to the immune group. They are not statistically discovered clusters or independent biological modules. Overlapping genes can connect different groups. [Assignment rules](heatmap_function_assignment_rules.tsv), [all-term assignments](heatmap_all_terms_function_audit.tsv), and [row-selection reasons](heatmap_term_selection_audit.tsv) are provided.', '',
    '## Support in the full integrated GSEA', '',
    'Counts below use every tested integrated term, including terms absent from the heatmap. Groups are mutually exclusive display labels; pathway counts are not independent biological observations.', '',
    '| Functional group | Tested | Older, q < 0.05 | Younger, q < 0.05 | Heatmap rows |',
    '|---|---:|---:|---:|---:|']
for group in order + [g for g in sorted({r['biology_group'] for r in integrated_audit}) if g not in order]:
    rows = [r for r in integrated_audit if r['biology_group'] == group]
    counts = Counter(r['direction'] for r in rows if qpass(r, .05))
    summary.append(f"| {group} | {len(rows)} | {counts['older_aligned']} | {counts['younger_aligned']} | {sum(r['biology_group'] == group for r in selection)} |")
summary += ['', '## Selected metabolic evidence', '',
    '| Term | Integrated NES | Integrated q |', '|---|---:|---:|']
for name in ['REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT',
             'REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
             'HALLMARK_FATTY_ACID_METABOLISM', 'WP_FATTY_ACID_BETAOXIDATION',
             'GOBP_SPHINGOLIPID_BIOSYNTHETIC_PROCESS', 'WP_AMINO_ACID_METABOLISM', 'REACTOME_TRANSLATION']:
    row = next(r for r in integrated_audit if r['gene_set_name'] == name)
    summary.append(f"| {name} | {number(row, 'NES'):.3f} | {number(row, 'fdr_within_pc'):.4f} |")
summary += ['', 'Fatty-acid/energy terms are predominantly younger-aligned, while lipid-response and sphingolipid terms can be older-aligned. The evidence therefore supports separate metabolic functions with mixed directions, not a single direction for all lipid biology. The evidence table also includes significant pathways that were not selected as nonredundant heatmap representatives.', '',
    'No steroid/sterol group term passes the integrated q < 0.05 threshold; a significant steroid cluster was not forced into the plot. PC5 has no terms at q < 0.05, although its leading terms meet q < 0.10; their **+** symbols preserve this distinction.', '',
    '## Top five terms per PC and direction', '',
    '| PC | Direction | Rank | Term | Function | NES | q |', '|---|---|---:|---|---|---:|---:|']
for pc in pcs:
    for direction in ['older_aligned', 'younger_aligned']:
        for row in sorted((r for r in audit if r['pc'] == pc and r['direction'] == direction and int(r['directional_rank']) <= 5), key=lambda r: int(r['directional_rank'])):
            summary.append(f"| {pc} | {direction.removesuffix('_aligned')} | {row['directional_rank']} | {row['gene_set_name']} | {row['biology_group']} | {number(row, 'NES'):.3f} | {number(row, 'fdr_within_pc'):.4f} |")
summary += ['', '## Reproduction and verification', '',
    f'All {len(checks)} heatmap validation checks passed. The existing 500-permutation weighted-KS estimates and their Monte Carlo resolution are retained. See [validation](../provenance/heatmap_revision_validation.tsv) and the main diagnostic README for the original statistical method.', '',
    'Run `Rscript scripts/06_refresh_gsea_figures.R --heatmaps-only` using the existing scRNA_r4 environment, then `python3 scripts/07_validate_heatmap.py` and `python3 scripts/04_report.py` from the diagnostic directory. The full analysis runner also includes the new validation. Previous heatmaps and plotting source are archived in `provenance/heatmap_before_functional_revision_20261012/`.', '']
(out / 'functional_heatmap_review.md').write_text('\n'.join(summary))
print(f'Passed {len(checks)} checks: {len(selected)} heatmap terms, {len(order)} functional groups, 40 required PC/rank positions.')
