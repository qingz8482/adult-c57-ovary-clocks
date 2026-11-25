#!/usr/bin/env python3
"""Inventory all application files and check text/table integrity without drawing."""
import ast
import csv
import json
from pathlib import Path

out = Path(__file__).resolve().parent.parent
rows = []
for path in sorted(out.rglob('*')):
    if not path.is_file() or path.name == 'application_file_review.tsv':
        continue
    rel = str(path.relative_to(out))
    kind, records = 'file present; included in SHA-256 manifest', ''
    if path.suffix in ('.tsv', '.csv'):
        with path.open(newline='') as f:
            reader = csv.reader(f, delimiter='\t' if path.suffix == '.tsv' else ',')
            header = next(reader, None)
            assert header, ('Missing table header', rel)
            assert len(set(header)) == len(header), ('Duplicate column names', rel)
            records = 0
            for line in reader:
                if not line:
                    continue
                assert len(line) == len(header), ('Ragged table row', rel, records+2)
                records += 1
        kind = 'tabular structure and unique columns verified'
    elif path.suffix == '.json':
        json.loads(path.read_text())
        kind = 'JSON parsed'
    elif path.suffix == '.py':
        ast.parse(path.read_text())
        kind = 'Python syntax parsed'
    elif path.suffix in ('.md', '.txt', '.R', '.sh', '.log'):
        path.read_text()
        kind = 'text readable'
    elif path.suffix in ('.pdf', '.png'):
        kind = 'figure; see comprehensive_figure_file_checks.tsv and visual QA'
    elif path.suffix == '.rds':
        kind = 'R input/cache; validated by projection/replay workflow; SHA-256 recorded'
    rows.append({'relative_path':rel,'bytes':path.stat().st_size,'records':records,'review':kind})
with (out/'provenance/application_file_review.tsv').open('w') as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
print(f'PASS: inventoried {len(rows)} application files; text, JSON and tabular structures checked.')
