"""Inventory original course figures and all their source-slide uses."""
import csv
import re
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
uses={}
for page in ROOT.glob('*/*.qmd'):
    if page.parts[-2] in ['docs','instructor']: continue
    for alt,target in re.findall(r'!\[([^\]]*)\]\(([^)]+)\)',page.read_text()):
        if '://' in target: continue
        figure=(page.parent/target).resolve()
        uses.setdefault(figure,[]).append(str(page.relative_to(ROOT)))
rows=[]
for figure in sorted(list(ROOT.glob('*/figures/*.svg'))+list((ROOT/'assets').glob('*.png'))):
    rows.append({'file':str(figure.relative_to(ROOT)),
                 'used_in':'; '.join(sorted(set(uses.get(figure.resolve(),[])))),
                 'source':'Original course artwork' if figure.suffix=='.svg' else 'Original R plot of constructed teaching data',
                 'licence':'Course owner to specify reuse terms',
                 'modifications':'Original drawing; no external image copied',
                 'review_status':'Author checked; cross-unit scientific spot checks; classroom review pending'})
with (ROOT/'assets/credits.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['file','used_in','source','licence','modifications','review_status'])
    writer.writeheader();writer.writerows(rows)
print(f'Inventoried {len(rows)} original visual assets.')
