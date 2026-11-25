#!/usr/bin/env python3
"""Record SHA-256 of this run's core numerical outputs; later figure stages must leave them unchanged."""
import hashlib, json
from pathlib import Path
out = Path(__file__).resolve().parent.parent
core = ['external_projection_predictions.tsv', 'external_projection_contrast_statistics.tsv',
        'tables/external_projection_relative_to_controls.tsv', 'inputs/model_gene_expression.rds']
p = out / 'provenance/plot_revision_baseline.json'
b = json.loads(p.read_text())
b['core_sha256'] = {n: hashlib.sha256((out / n).read_bytes()).hexdigest() for n in core}
p.write_text(json.dumps(b, indent=1) + '\n')
print('Recorded numerical baseline for', len(core), 'files')
