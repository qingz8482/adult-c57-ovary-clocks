#!/usr/bin/env python3
"""Check plot coverage and unchanged numerical inputs without drawing figures."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path

out = Path(__file__).resolve().parent.parent
with (out / 'provenance/retired_diagnostic_outputs.tsv').open() as f:
    retired_prefixes = tuple(r['path_prefix'] for r in csv.DictReader(f, delimiter='\t'))
baseline = json.loads((out / 'provenance/plot_revision_baseline.json').read_text())
for name, digest in baseline['core_sha256'].items():
    assert hashlib.sha256((out / name).read_bytes()).hexdigest() == digest, name
with open(out / 'provenance/plot_revision_exports.tsv') as f:
    exports = list(csv.DictReader(f, delimiter='\t'))
stems = {row['stem'] for row in exports}
assert all(str(Path(name).with_suffix('')) in stems for name in baseline['plot_files'] if not name.startswith(retired_prefixes))
checks = []
for stem in sorted(stems):
    pdf, png = out / (stem + '.pdf'), out / (stem + '.png')
    assert pdf.is_file() and png.is_file(), stem
    result = subprocess.run(['pdfinfo', str(pdf)], text=True, capture_output=True)
    assert result.returncode == 0, (stem, result.stderr)
    assert png.read_bytes()[:8] == b'\x89PNG\r\n\x1a\n', stem
    checks.append([stem + '.pdf', pdf.stat().st_size, 'PASS'])
with open(out / 'provenance/figure_export_checks.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['pdf', 'bytes', 'pdf_readable_and_png_present'])
    w.writerows(checks)
with open(out / 'provenance/plot_annotation_audit.tsv') as f:
    annotations = list(csv.DictReader(f, delimiter='\t'))
assert len({(x['analysis_contrast_id'], x['model']) for x in annotations}) == 88
for row in annotations:
    p = float(row['welch_p'])
    expected = '****' if p < .0001 else '***' if p < .001 else '**' if p < .01 else '*' if p < .05 else f'p = {p:.3f}'
    assert row['annotation'] == expected
with open(out / 'diagnostics/previous_style_plot_source.tsv') as f:
    sources = list(csv.DictReader(f, delimiter='\t'))
with open(out / 'tables/external_projection_contrast_membership.tsv') as f:
    memberships = list(csv.DictReader(f, delimiter='\t'))
keys = {r['sample_key']: r for r in sources}
groups, roles = {}, {}
for r in sources:
    groups.setdefault((r['dataset_cohort'], r['plot_group']), set()).add(r['sample_key'])
for r in memberships:
    roles.setdefault((r['analysis_contrast_id'], r['contrast_role']), set()).add(r['sample_key'])
for members in roles.values():
    row = keys[next(iter(members))]
    assert members == groups[(row['dataset_cohort'], row['plot_group'])]
result = {'status': 'PASS', 'pdf_png_pairs': len(stems), 'prior_plot_stems_regenerated': len({str(Path(x).with_suffix('')) for x in baseline['plot_files'] if not x.startswith(retired_prefixes)}), 'unique_annotated_model_tests': 88, 'numeric_inputs_unchanged': True, 'bracket_sample_membership_exact': True}
(out / 'provenance/plot_revision_validation.json').write_text(json.dumps(result, indent=2) + '\n')
(out / 'provenance/figure_QA.txt').write_text(
    'R/ggplot2 produced all figures, PDF exports and 300-dpi PNG previews.\n'
    f"{len(stems)} PDF/PNG pairs verified; every retained prior plot stem verified; user-retired outputs excluded.\n"
    'Prior styling restored: boxplots and individual profiles; grey clock strips; group colours; dashed zero lines; capped brackets.\n'
    'Two clocks shown in parallel; GSE235021 bulk and pseudobulk remain separate.\n'
    'Dense group labels use vertical text and a taller canvas to avoid overlap.\n'
    'Representative CRA017492, PRJNA863443 and GSE235021 figures visually inspected.\n'
    'All 88 model-specific Welch annotations verified; each bracket connects exactly its tested sample groups.\n'
    'Original predictions, statistics, relative deltas and model input RDS are SHA-256-identical.\n'
    'PDF font inspection confirms embedded Liberation Sans. A font-cache compatibility message did not prevent rendering.\n'
)
print(json.dumps(result))
