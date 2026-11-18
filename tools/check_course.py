"""Check complete course sources, meaningful-visual quotas and built links.

An image count is a mechanical proxy. Human review still has to check chemical
correctness, whether a visual is substantive, and projected readability.
"""
import csv
import json
import re
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit
from zipfile import ZipFile
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
units=json.loads((ROOT/'plan/units.json').read_text())
errors=[]
inventory=[]

def visible_source(text):
    text=re.sub(r'^---\n.*?\n---\n','',text,flags=re.S)
    text=re.sub(r'::: notes\n.*?\n:::', '',text,flags=re.S)
    text=re.sub(r'```.*?```','',text,flags=re.S)
    return text

for unit in units:
    slug=unit['slug']
    pages=[ROOT/slug/(name+'.qmd') for name in ['lec','exercise','quiz']]
    for f in pages:
        if not f.exists(): errors.append(f'Missing {f.relative_to(ROOT)}')
    if not all(f.exists() for f in pages): continue
    lecture=pages[0].read_text()
    body=visible_source(lecture)
    sections=re.split(r'^## .+$',body,flags=re.M)[1:]
    # One automatic title slide is present in these Reveal decks.
    total=1+len(sections)
    visual=sum(bool(re.search(r'!\[[^\]]+\]\([^)]+\.svg\)',section)) for section in sections)
    if total<15: errors.append(f'{slug}: only {total} slides')
    if visual<=total/2: errors.append(f'{slug}: visual quota fails ({visual}/{total})')
    quiz=pages[2].read_text()
    if len(re.findall(r'^## ',quiz,re.M))<5: errors.append(f'{slug}: fewer than five quiz questions')
    if len(re.findall(r'collapse\s*=\s*["\']true',quiz))<5:
        errors.append(f'{slug}: fewer than five collapsible answer explanations')
    for f in pages:
        text=f.read_text()
        for alt,ref in re.findall(r'!\[([^\]]*)\]\(([^)]+)\)',text):
            if not alt.strip(): errors.append(f'{f.relative_to(ROOT)}: empty image alt text')
            if '://' not in ref and not (f.parent/ref.split('#')[0]).exists():
                errors.append(f'{f.relative_to(ROOT)}: missing image {ref}')
    inventory.append({'unit':slug,'slides':total,'visual_slides':visual,'visual_percent':round(100*visual/total,1)})

for f in ROOT.glob('*/figures/*.svg'):
    try:
        xml=ET.parse(f)
        if xml.find('{http://www.w3.org/2000/svg}title') is None:
            errors.append(f'{f.relative_to(ROOT)}: missing SVG title')
    except ET.ParseError as e: errors.append(f'{f.relative_to(ROOT)}: invalid SVG: {e}')

class Links(HTMLParser):
    def __init__(self): super().__init__(); self.refs=[]; self.ids=set(); self.missing_alt=[]
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'id' in a: self.ids.add(a['id'])
        field={'a':'href','img':'src','script':'src','link':'href','iframe':'src','object':'data'}.get(tag)
        if field and a.get(field): self.refs.append(a[field])
        if tag=='img':
            if a.get('data-src'): self.refs.append(a['data-src'])
            if not a.get('alt','').strip(): self.missing_alt.append(a.get('src') or a.get('data-src','unnamed image'))

docs=ROOT/'docs'
pages={}
for f in docs.rglob('*.html'):
    parser=Links(); parser.feed(f.read_text(errors='replace')); pages[f.resolve()]=parser
if not (docs/'index.html').exists(): errors.append('Built docs/index.html is missing: run make render')
for f,parser in pages.items():
    for image in parser.missing_alt: errors.append(f'{f.relative_to(ROOT)}: missing built image alt: {image}')
    for ref in parser.refs:
        u=urlsplit(ref)
        if u.scheme or u.netloc or not u.path: continue
        target=(docs/unquote(u.path.lstrip('/')) if u.path.startswith('/') else f.parent/unquote(u.path)).resolve()
        if not target.exists(): errors.append(f'{f.relative_to(ROOT)}: missing built target {ref}')
        elif target.is_dir() and not (target/'index.html').exists():
            errors.append(f'{f.relative_to(ROOT)}: directory without index: {ref}')
        elif target.suffix=='.qmd': errors.append(f'{f.relative_to(ROOT)}: unresolved source link {ref}')

for forbidden in ['instructor','plan','tools','contentplan.html','README.html']:
    if (docs/forbidden).exists(): errors.append(f'Private/planning output present: docs/{forbidden}')
archive=ROOT/'downloads/organic-chemistry-student.zip'
if archive.exists():
    with ZipFile(archive) as z:
        for required in ['report.Rmd','organic-chemistry.Rproj','data/amounts.csv','data/README.md']:
            if required not in z.namelist(): errors.append(f'Archive missing {required}')
        for name in z.namelist():
            if name.startswith(('instructor/','plan/','tools/')): errors.append(f'Private archive input: {name}')
else: errors.append('Student ZIP missing')

with (ROOT/'plan/slide-inventory.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['unit','slides','visual_slides','visual_percent'])
    writer.writeheader();writer.writerows(inventory)
print(f'Units: {len(inventory)}/{len(units)}; slides: {sum(x["slides"] for x in inventory)}; '
      f'visual slides: {sum(x["visual_slides"] for x in inventory)}; built HTML pages: {len(pages)}')
for row in inventory: print(f'  {row["unit"]}: {row["visual_slides"]}/{row["slides"]} visual slides')
if errors:
    for error in errors: print('ERROR:',error)
    raise SystemExit(1)
print('Source structure, SVGs, per-deck quotas, built local links and download boundaries pass.')
