"""Require the learner-facing Week → Day → materials scaffold."""
import json
import re
from pathlib import Path
import yaml

ROOT=Path(__file__).resolve().parents[1]
cfg=yaml.safe_load((ROOT/'_quarto.yml').read_text())
days=json.loads((ROOT/'plan/study-days.json').read_text())
demos=json.loads((ROOT/'plan/demos.json').read_text())
weeks=[x for x in cfg['website']['sidebar']['contents'] if x.get('section','').startswith('Week ')]
errors=[]
if len(weeks)!=5:errors.append('Expected five week sections')
for number,week in enumerate(weeks,1):
    groups=week.get('contents',[])
    if len(groups)!=5:errors.append(f'Week {number}: expected five day groups')
    for position,group in enumerate(groups):
        n=(number-1)*5+position+1
        if not group.get('section','').startswith(f'Day {n}:'):
            errors.append(f'Week {number}: missing Day {n} section');continue
        day=days[n-1];demo=demos[n-1]
        leaves=group.get('contents',[])
        exercise=[x for x in leaves if x.get('href')==day['href']]
        if len(exercise)!=1 or exercise[0].get('text')!='Exercises':
            errors.append(f'Day {n}: exercise link must be labelled Exercises')
        if not any(x.get('href')==day['href']+'#explore' and x.get('text','').startswith('Demo:') for x in leaves):
            errors.append(f'Day {n}: missing embedded demo link')
        if not any(x.get('href','').endswith('/lec.qmd') for x in leaves):
            errors.append(f'Day {n}: missing lecture/review link')
        if day['kind']=='unit' and not any(x.get('href')==day['unit']+'/quiz.qmd' for x in leaves):
            errors.append(f'Day {n}: missing self-check')
        for leaf in leaves:
            if 'contents' in leaf:errors.append(f'Day {n}: unexpected nested material group')
            if not (ROOT/leaf.get('href','').split('#')[0]).is_file():
                errors.append(f'Day {n}: missing target {leaf.get("href")}')
        if not (ROOT/demo['include_href']).is_file():errors.append(f'Day {n}: missing demo content include')
for page in [ROOT/d['href'] for d in days]:
    text=page.read_text()
    if re.search(r'standalone (?:paper/model workshop|(?:practice|problem) bank)',text,re.I):
        errors.append(str(page.relative_to(ROOT))+': workshop/bank boilerplate remains')
if 'Practice and worked examples' in (ROOT/'_quarto.yml').read_text():
    errors.append('Old exercise label remains in sidebar')
if errors:
    for error in errors:print('ERROR:',error)
    raise SystemExit(1)
print('Five weeks / 25 day groups each contain lecture, Exercises and the embedded demo.')
