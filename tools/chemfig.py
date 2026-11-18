"""Small, dependency-free SVG drawing helpers for original course diagrams.

Coordinates are in a 960 × 400 canvas. Chemical correctness is the author's
responsibility; these helpers do not infer bonds, valence or stereochemistry.
"""
from html import escape
from pathlib import Path

INK = '#153747'
BLUE = '#176b91'
RED = '#af3b45'
GREEN = '#23765d'

class Figure:
    def __init__(self, title, description='', width=960, height=400):
        self.width, self.height = width, height
        self.parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img"><title>{escape(title)}</title><desc>{escape(description or title)}</desc>',
          '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="#176b91"/></marker></defs>',
          f'<rect width="{width}" height="{height}" rx="12" fill="#f7fafc"/>']
    def text(self, x, y, label, size=26, color=INK, anchor='middle', weight='normal'):
        self.parts.append(f'<text x="{x}" y="{y}" font-family="DejaVu Sans, sans-serif" font-size="{size}" fill="{color}" text-anchor="{anchor}" font-weight="{weight}">{escape(str(label))}</text>')
        return self
    def line(self, x1, y1, x2, y2, color=INK, width=3, dash=None, arrow=False, both=False):
        attrs = (' stroke-dasharray="'+str(dash)+'"') if dash else ''
        attrs += ' marker-end="url(#arrow)"' if arrow else ''
        attrs += ' marker-start="url(#arrow)"' if both else ''
        self.parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"{attrs}/>')
        return self
    def path(self, d, color=BLUE, width=3, fill='none', arrow=False):
        attrs = ' marker-end="url(#arrow)"' if arrow else ''
        self.parts.append(f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="{fill}"{attrs}/>')
        return self
    def circle(self, x, y, radius=10, color=INK, fill='white', width=3):
        self.parts.append(f'<circle cx="{x}" cy="{y}" r="{radius}" stroke="{color}" fill="{fill}" stroke-width="{width}"/>')
        return self
    def rect(self, x, y, width, height, fill='#e4f0f4', color='none'):
        self.parts.append(f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="8" fill="{fill}" stroke="{color}"/>')
        return self
    def wedge(self, x1, y1, x2, y2, spread=12, hashed=False):
        import math
        dx, dy = x2-x1, y2-y1
        length = math.hypot(dx,dy)
        nx, ny = -dy/length*spread, dx/length*spread
        if hashed:
            for t in [0.2,0.35,0.5,0.65,0.8,0.95]:
                self.line(x1+t*dx-t*nx,y1+t*dy-t*ny,x1+t*dx+t*nx,y1+t*dy+t*ny,width=2)
        else:
            self.parts.append(f'<polygon points="{x1},{y1} {x2+nx},{y2+ny} {x2-nx},{y2-ny}" fill="{INK}"/>')
        return self
    def save(self, path):
        path=Path(path)
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text('\n'.join(self.parts+['</svg>']))

def reaction(title, left, right, condition='', note='', path=None):
    f=Figure(title, note)
    f.text(230,190,left,32).text(730,190,right,32)
    f.line(400,180,560,180,arrow=True,color=BLUE)
    f.text(480,140,condition,21).text(480,305,note,23)
    if path: f.save(path)
    return f
