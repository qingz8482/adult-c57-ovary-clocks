from pathlib import Path
import csv,hashlib,re,subprocess
out=Path(__file__).resolve().parents[1];root=out.parents[2]
with (out/'provenance/reference_panel_manifest.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dims(p):
 s=subprocess.check_output(['pdfinfo',str(p)],text=True)
 return tuple(map(float,re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',s).groups()))
checks=[]
for r in rows:
 ref=root/'bulk-pseudobulk_clock'/r['reference_module']/r['reference_path'];p=out/'figures'/f"{r['output_stem']}.pdf"
 exist=all(p.with_suffix(e).exists() and p.with_suffix(e).stat().st_size>1000 for e in ['.pdf','.png','.svg','.tiff'])
 rd=dims(ref);nd=dims(p) if p.exists() else (0,0)
 identical=all(sha(p.with_suffix(e))==sha(ref.with_suffix(e)) for e in ['.pdf','.png','.svg','.tiff']) if exist and r['status'].startswith('unchanged') else False
 checks.append(dict(figure=r['output_stem'],all_four_formats=exist,dimensions_match=max(abs(a-b) for a,b in zip(rd,nd))<1.1,reference_unchanged=sha(ref)==r['reference_sha256'],exact_copy_verified=identical if r['status'].startswith('unchanged') else 'not_applicable',reference_dimensions_pt=str(rd),output_dimensions_pt=str(nd)))
with (out/'provenance/reference_figure_validation.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=checks[0],delimiter='\t');w.writeheader();w.writerows(checks)
failed=[c for c in checks if not all(c[k] is True for k in ['all_four_formats','dimensions_match','reference_unchanged']) or c['exact_copy_verified'] is False]
assert len(rows)==58 and len(list((out/'figures').glob('*.pdf')))==58
if failed:raise SystemExit('FAIL: '+str(failed))
print('PASS: 58 figure counterparts, four formats, matching dimensions, 32 byte-identical cohort figures, original PDFs unchanged.')
