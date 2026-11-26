"""Verify frozen-release artifacts without changing any files."""
import csv, hashlib
from pathlib import Path
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'checksums.tsv').open(), delimiter='\t'))
for row in rows:
    p=root/row['path']
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
        raise SystemExit('FAIL: '+row['path'])
print(f'PASS: {len(rows)} locked files match SHA-256 checksums')
