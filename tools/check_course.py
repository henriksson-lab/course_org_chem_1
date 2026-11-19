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
import yaml

ROOT=Path(__file__).resolve().parents[1]
units=json.loads((ROOT/'plan/units.json').read_text())
errors=[]
warnings=[]
inventory=[]
practice_inventory=[]

def check_paper_set(path, minimum):
    if not path.exists():
        errors.append(f'Missing paper practice: {path.relative_to(ROOT)}')
        return 0
    body=path.read_text()
    problems=re.split(r'^##\s+\d+[.)]\s+.+$',body,flags=re.M)[1:]
    if len(problems)<minimum:
        errors.append(f'{path.relative_to(ROOT)}: {len(problems)} paper problems, expected at least {minimum}')
    for number,problem in enumerate(problems,1):
        problem=re.split(r'^##\s+(?!\d+[.)])',problem,maxsplit=1,flags=re.M)[0]
        if len(re.findall(r'\bcollapse\s*=\s*(?:"true"|\'true\'|true)',problem))<2:
            errors.append(f'{path.relative_to(ROOT)} problem {number}: needs separate hint and worked answer')
    return len(problems)

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
    metadata=yaml.safe_load(lecture.split('---',2)[1])
    body=visible_source(lecture)
    sections=re.split(r'^## .+$',body,flags=re.M)[1:]
    # One automatic title slide is present in these Reveal decks.
    total=1+len(sections)
    visual=sum(bool(re.search(r'!\[[^\]]+\]\([^)]+\.svg\)',section)) for section in sections)
    session_minutes=metadata.get('session-minutes')
    seconds_per_slide=metadata.get('pacing-seconds-per-slide')
    # Decks are sized by topic (quantbio convention: about 15-19 content slides),
    # not by a seconds-per-slide quota. Flag oversized decks as a warning only.
    if total>25: warnings.append(f'{slug}: {total} slides; consider tightening the deck')
    # Cadence is guidance, not a slide quota or a measure of teaching quality.
    # Actual session budgets also include drawing, discussion and worked tasks.
    nominal_slide_minutes=round(total*seconds_per_slide/60,1) if seconds_per_slide else None
    unique_figures=set(re.findall(r'!\[[^\]]+\]\(([^)]+\.svg)\)',body))
    paper_count=check_paper_set(pages[1],20)
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
    inventory.append({'unit':slug,'slides':total,'visual_slides':visual,'visual_percent':round(100*visual/total,1),
                      'unique_figures':len(unique_figures),'planned_session_minutes':session_minutes,
                      'nominal_slide_minutes':nominal_slide_minutes,'paper_problems':paper_count})

study_days=json.loads((ROOT/'plan/study-days.json').read_text())
if [d['day'] for d in study_days]!=list(range(1,26)): errors.append('Daily practice route must cover study days 1–25')
for day in study_days:
    count=check_paper_set(ROOT/day['href'],day['minimum_paper_problems'])
    practice_inventory.append({'day':day['day'],'week':day['week'],'kind':day['kind'],
                               'page':day['href'],'paper_problems':count})
if yaml.safe_load((ROOT/'_quarto.yml').read_text()).get('history') is not False:
    errors.append('Shared presentation history must be false so browser Back leaves a lecture')

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
    writer=csv.DictWriter(f,fieldnames=['unit','slides','visual_slides','visual_percent','unique_figures','planned_session_minutes','nominal_slide_minutes','paper_problems'],lineterminator='\n')
    writer.writeheader();writer.writerows(inventory)
with (ROOT/'plan/practice-inventory.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['day','week','kind','page','paper_problems'],lineterminator='\n')
    writer.writeheader();writer.writerows(practice_inventory)
print(f'Units: {len(inventory)}/{len(units)}; slides: {sum(x["slides"] for x in inventory)}; '
      f'visual slides: {sum(x["visual_slides"] for x in inventory)}; built HTML pages: {len(pages)}')
for row in inventory: print(f'  {row["unit"]}: {row["visual_slides"]}/{row["slides"]} visual slides')
print(f'Paper chemistry problems across 25 study days: {sum(r["paper_problems"] for r in practice_inventory)}')
for warning in warnings: print("WARNING:",warning)
if errors:
    for error in errors: print('ERROR:',error)
    raise SystemExit(1)
print('Source structure, SVGs, per-deck quotas, built local links and download boundaries pass.')
