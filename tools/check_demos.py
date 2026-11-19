"""Check the daily demo pages, local scripts and direct sidebar destinations."""
import json
import re
import subprocess
from pathlib import Path
import yaml
from source_includes import expand_includes

ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'plan/demos.json').read_text())
config=yaml.safe_load((ROOT/'_quarto.yml').read_text())
errors=[]
def links(items):
    result=[]
    for item in items:
        if isinstance(item,dict):
            if 'href' in item: result.append(item['href'])
            result+=links(item.get('contents',[]))
    return result
sidebar=links(config['website']['sidebar']['contents'])
if [d['day'] for d in manifest]!=list(range(1,26)):
    errors.append('Demo route must map exactly study days 1–25')
for demo in manifest:
    source=ROOT/demo['href']
    if not source.exists(): errors.append('Missing '+demo['href']);continue
    body=expand_includes(ROOT/demo['exercise_href'])
    if not re.search(r'id=["\']'+demo['root_id']+r'["\']',body):
        errors.append(demo['href']+': missing matching mount root')
    if '../interactives/core.js' not in body or '../'+demo['module'] not in body:
        errors.append(demo['href']+': missing local runtime scripts')
    if demo['menu_href'] not in sidebar: errors.append(demo['menu_href']+': missing direct sidebar link')
    if '{#explore}' not in (ROOT/demo['exercise_href']).read_text():
        errors.append(demo['exercise_href']+': missing Explore anchor')
    target=ROOT/'docs'/Path(demo['exercise_href']).with_suffix('.html')
    if not target.exists(): errors.append('Missing built '+str(target.relative_to(ROOT)))
    if not (ROOT/'docs'/demo['module']).exists(): errors.append('Missing deployed runtime '+demo['module'])
    if 'window.location.replace' not in source.read_text(): errors.append(demo['href']+': old demo URL must redirect to its exercise')
for script in (ROOT/'interactives').glob('*.js'):
    result=subprocess.run(['node','--check',str(script)],capture_output=True,text=True)
    if result.returncode: errors.append(result.stderr)
    if re.search(r'\bfetch\s*\(|^\s*(?:import|export)\s',script.read_text(),re.M):
        errors.append(str(script.relative_to(ROOT))+': demos should use local classic scripts')
if errors:
    for error in errors:print('ERROR:',error)
    raise SystemExit(1)
print('All 25 daily demos have local runtimes and direct sidebar links.')
