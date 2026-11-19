"""Figures for the Week 5 spectroscopy lecture (RDKit skeletal structures via molfig).

Spectra are schematic, simulated teaching traces drawn with chemfig. Shifts and
splittings follow standard literature values for the compounds shown; the line
spacing inside a multiplet is exaggerated so the pattern is visible.
"""
import math
from pathlib import Path
from rdkit import Chem
from chemfig import Figure, BLUE, RED, GREEN, INK
import molfig
from molfig import place

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'spectroscopy'/'figures'
GREY = '#8aa0ab'
# Pale shades used to highlight one H environment in the colour of its signal.
PALE = {BLUE: '#cfe3ec', RED: '#f3d5d8', GREEN: '#d2eadf', INK: '#dde3e7'}

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

def with_h(smiles, carbons):
    """Indices of carbons plus their hydrogens when the molecule is drawn with explicit_h."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    hs = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'H'
          and a.GetNeighbors()[0].GetIdx() in carbons]
    return list(carbons) + hs

def mol(f, smiles, x, y, colour=None, **kw):
    """place() with a per-call highlight colour (molfig has one global HIGHLIGHT)."""
    old = molfig.HIGHLIGHT
    if colour: molfig.HIGHLIGHT = PALE.get(colour, colour)
    try:
        return place(f, smiles, x, y, **kw)
    finally:
        molfig.HIGHLIGHT = old

class Axis:
    """Horizontal shift axis; values decrease to the right, as on a real spectrum."""
    def __init__(self, f, x0, x1, y, hi, lo, ticks, label=None, size=18):
        self.f, self.x0, self.x1, self.y, self.hi, self.lo = f, x0, x1, y, hi, lo
        f.line(x0, y, x1, y, color=INK, width=2)
        for t in ticks:
            f.line(self.x(t), y, self.x(t), y + 7, color=INK, width=2)
            f.text(self.x(t), y + 27, t, size, color=INK)
        if label: f.text((x0 + x1)/2, y + 56, label, 20)
    def x(self, v):
        return self.x0 + (self.hi - v)/(self.hi - self.lo)*(self.x1 - self.x0)

def multiplet(ax, ppm, lines, n_h, colour, k=26, gap=8, width=3):
    """Stick multiplet with binomial line heights; total height is proportional to n_h."""
    weights = [math.comb(lines - 1, i) for i in range(lines)]
    total = sum(weights)
    top = ax.y
    for i, w in enumerate(weights):
        xx = ax.x(ppm) + (i - (lines - 1)/2)*gap
        h = k*n_h*w/total
        ax.f.line(xx, ax.y, xx, ax.y - h, color=colour, width=width)
        top = min(top, ax.y - h)
    return top

def broad(ax, ppm, height, half, colour):
    """A broad exchangeable O-H hump."""
    x, y = ax.x(ppm), ax.y
    pts = ' '.join(f'L {x + d:.1f} {y - height*math.exp(-(d/half)**2*3):.1f}'
                   for d in range(-int(2*half), int(2*half) + 1, 2))
    ax.f.path('M' + pts[1:], color=colour, width=3)

def ir_trace(f, x0, x1, y_top, y_bot, bands, colour=BLUE, wmax=4000, wmin=600):
    """Transmittance trace: bands are (wavenumber, depth 0-1, half-width)."""
    pts = []
    for i in range(301):
        w = wmax - (wmax - wmin)*i/300
        t = 0.95 - sum(d*math.exp(-((w - c)/s)**2) for c, d, s in bands)
        pts.append(f'{x0 + (x1 - x0)*i/300:.1f} {y_bot - (y_bot - y_top)*t:.1f}')
    f.path('M ' + ' L '.join(pts), color=colour, width=3)
    return lambda w: x0 + (wmax - w)/(wmax - wmin)*(x1 - x0)

ETAC = 'CC(=O)OCC'       # ethyl ethanoate: 0 acetyl CH3, 1 C=O, 2 =O, 3 O, 4 CH2, 5 CH3
BUAC = 'CCCC(=O)O'       # butanoic acid: 0 CH3, 1 CH2, 2 CH2, 3 C=O, 4 =O, 5 OH
ACETOIN = 'CC(=O)C(O)C'  # 3-hydroxybutan-2-one: 0 CH3, 1 C=O, 2 =O, 3 CH, 4 OH, 5 CH3
DIOXANE = 'C1COCCO1'

# 1. Four C4H8O2 candidates, each with one degree of unsaturation.
f = Figure('Four structures with formula C4H8O2', 'Ethyl ethanoate, butanoic acid, 3-hydroxybutan-2-one and 1,4-dioxane, each with its single degree of unsaturation shaded', 960, 400)
for i, (smi, name, sub, hl) in enumerate([
        (ETAC, 'ethyl ethanoate', 'C=O (ester)', [1, 2]),
        (BUAC, 'butanoic acid', 'C=O (acid)', [3, 4]),
        (ACETOIN, '3-hydroxybutan-2-one', 'C=O (ketone)', [1, 2]),
        (DIOXANE, '1,4-dioxane', 'one ring', [0, 1, 2, 3, 4, 5])]):
    x = 120 + 240*i
    mol(f, smi, x, 150, highlight_atoms=hl, flip=(smi == ETAC), bond_px=42, font_px=24)
    f.text(x, 275, name, 21, weight='bold').text(x, 308, sub, 20, color=BLUE)
f.text(480, 370, 'All C₄H₈O₂ · DBE = (2×4 + 2 − 8)/2 = 1', 24)
save('candidates', f)

# 2. Three IR traces with the responsible group shaded.
f = Figure('IR traces of an ester, a carboxylic acid and an alcohol', 'Simulated transmittance traces: ethyl ethanoate shows C=O near 1740 and no O-H; butanoic acid shows a very broad O-H from 2500 to 3300 plus C=O near 1710; ethanol shows a broad O-H near 3350 and no C=O', 960, 490)
rows = [
    (ETAC, 'ester', [1, 2], [(2980, .12, 60), (1740, .62, 25), (1240, .45, 35), (1050, .3, 25)], RED, [(1740, 'C=O')]),
    (BUAC, 'acid', [3, 4, 5], [(2950, .45, 380), (1710, .6, 28), (1250, .3, 35)], GREEN, [(2950, 'very broad O–H'), (1710, 'C=O')]),
    ('CCO', 'alcohol', [2], [(3350, .5, 170), (2970, .15, 45), (1050, .5, 25)], BLUE, [(3350, 'broad O–H')]),
]
for r, (smi, name, hl, bands, col, labels) in enumerate(rows):
    yt, yb = 50 + 125*r, 145 + 125*r
    mol(f, smi, 120, (yt + yb)/2 - 4, colour=col, highlight_atoms=hl, flip=(smi == ETAC), bond_px=40, font_px=21)
    f.text(120, yb + 8, name, 19, weight='bold', color=col)
    wx = ir_trace(f, 260, 930, yt, yb, bands, colour=col)
    for w, lab in labels:
        f.text(wx(w), yb + 10, lab, 17, color=col, weight='bold')
for w in (4000, 3000, 2000, 1000):
    x = 260 + (4000 - w)/3400*670
    f.line(x, 420, x, 427, width=2).text(x, 446, w, 17)
f.line(260, 420, 930, 420, width=2)
f.text(595, 474, 'wavenumber / cm⁻¹  ·  transmittance: dips mean absorption', 18)
save('ir-groups', f)

# 3. U1 proton spectrum with each signal tied to its H environment.
f = Figure('Proton NMR of U1 mapped onto ethyl ethanoate', 'Schematic 1H NMR: quartet at 4.12 ppm (2H, O-CH2), singlet at 2.05 ppm (3H, acetyl CH3), triplet at 1.25 ppm (3H, ethyl CH3), each shown under a copy of ethyl ethanoate with that environment shaded in the same colour', 960, 490)
ax = Axis(f, 60, 900, 200, 5, 0, [5, 4, 3, 2, 1, 0], 'chemical shift δ / ppm', size=17)
signals = [(4.12, 4, 2, RED, [4], '4.12 · q · 2H', 175), (2.05, 1, 3, GREEN, [0], '2.05 · s · 3H', 480),
           (1.25, 3, 3, BLUE, [5], '1.25 · t · 3H', 760)]
for ppm, lines, nh, col, carbons, lab, mx in signals:
    top = multiplet(ax, ppm, lines, nh, col, k=42, gap=10, width=4)
    f.text(ax.x(ppm), top - 14, lab, 20, color=col, weight='bold')
    mol(f, ETAC, mx, 395, colour=col, highlight_atoms=with_h(ETAC, carbons), explicit_h=True, flip=True, bond_px=40, font_px=19)
save('u1-proton', f)

# 4. Counting H environments: letters mark equivalent sets.
f = Figure('Counting proton environments', 'Ethyl ethanoate has three H environments, butanoic acid four and 1,4-dioxane one; equivalent hydrogens share a letter', 960, 380)
for x, smi, notes, name, n in [
        (150, ETAC, {0: 'a', 4: 'b', 5: 'c'}, 'ethyl ethanoate', '3 signals'),
        (500, BUAC, {0: 'a', 1: 'b', 2: 'c', 5: 'd'}, 'butanoic acid', '4 signals'),
        (800, DIOXANE, {0: 'a', 1: 'a', 3: 'a', 4: 'a'}, '1,4-dioxane', '1 signal (8H)')]:
    mol(f, smi, x, 150, numbers=notes, flip=(smi == ETAC), bond_px=54, font_px=32)
    f.text(x, 290, name, 22, weight='bold').text(x, 325, n, 22, color=BLUE)
save('environments', f)

# 5. Chemical-shift map.
f = Figure('Where protons appear on the shift scale', 'Typical 1H shift ranges: alkyl C-H 0.9 to 1.7, H-C-C=O 2.0 to 2.6, H-C-O 3.3 to 4.3, alkene 4.6 to 6.5, aromatic 6.5 to 8, aldehyde 9.5 to 10, carboxylic acid 10 to 12 ppm', 960, 360)
ax = Axis(f, 40, 920, 250, 12, 0, list(range(12, -1, -2)), 'chemical shift δ / ppm', size=17)
f.text(40, 330, '← downfield: H near O or C=O', 18, anchor='start', color=RED).text(920, 330, 'upfield: plain alkyl H →', 18, anchor='end', color=INK)
bands = [(0.9, 1.7, 'alkyl C–H', INK, 0), (2.0, 2.6, 'H–C–C=O', GREEN, 1), (3.3, 4.3, 'H–C–O', RED, 0),
         (4.6, 6.5, 'C=C–H', GREY, 1), (6.5, 8.0, 'Ar–H', GREY, 0), (9.5, 10.0, 'CHO', BLUE, 1), (10.0, 12.0, 'COOH', BLUE, 0)]
for lo, hi, lab, col, row in bands:
    x0, x1 = ax.x(hi), ax.x(lo)
    y = 150 - 70*row
    f.rect(x0, y, x1 - x0, 46, fill=PALE.get(col, '#e7ecef'))
    f.text((x0 + x1)/2, y - 10, lab, 19, color=col if col != GREY else INK, weight='bold')
for ppm, col in [(4.12, RED), (2.05, GREEN), (1.25, BLUE)]:
    f.line(ax.x(ppm), 210, ax.x(ppm), 248, color=col, width=4)
f.text(ax.x(3.1), 233, 'U1', 18, color=INK, weight='bold')
save('shift-scale', f)

# 6. Integration: area, not height.
f = Figure('Integrals measure area, not line height', 'Schematic U1 spectrum with integral step curves; relative step heights 4, 6 and 6 reduce to 2:3:3, totalling 8 H', 960, 400)
ax = Axis(f, 60, 900, 300, 5, 0, [5, 4, 3, 2, 1, 0], 'chemical shift δ / ppm', size=17)
level = 250
steps = []
for ppm, lines, nh, col in [(4.12, 4, 2, RED), (2.05, 1, 3, GREEN), (1.25, 3, 3, BLUE)]:
    multiplet(ax, ppm, lines, nh, col, k=40, gap=10, width=4)
    steps.append((ppm, nh, col))
y = 260
pts = [f'M {ax.x(5):.1f} {y}']
for ppm, nh, col in steps:
    xa, xb = ax.x(ppm) - 32, ax.x(ppm) + 32
    pts.append(f'L {xa:.1f} {y}')
    y2 = y - 22*nh
    pts.append(f'C {xa + 32:.1f} {y} {xb - 32:.1f} {y2} {xb:.1f} {y2}')
    f.text(xb + 22, (y + y2)/2 + 8, f'{2*nh}', 24, color=col, weight='bold', anchor='start')
    y = y2
pts.append(f'L {ax.x(0):.1f} {y}')
f.path(' '.join(pts), color=GREY, width=2.5)
f.text(480, 50, 'step heights 4 : 6 : 6  →  2 : 3 : 3  →  2 + 3 + 3 = 8 H = C₄H₈O₂', 22)
save('integrals', f)

# 7. Singlet versus doublet: what sits next door.
f = Figure('A methyl next to C=O is a singlet; next to CH it is a doublet', 'Left: acetyl CH3 of ethyl ethanoate beside a carbonyl carbon with no H gives a 3H singlet. Right: the CH3 of 3-hydroxybutan-2-one beside a CH with one H gives a 3H doublet', 960, 480)
for x, smi, hl, col, lines, lab, sub in [
        (240, ETAC, [0], GREEN, 1, 'next door: C=O carbon, 0 H', 'n = 0 → singlet'),
        (720, ACETOIN, [5], BLUE, 2, 'next door: CH, 1 H', 'n = 1 → doublet')]:
    mol(f, smi, x, 125, colour=col, highlight_atoms=with_h(smi, hl), explicit_h=True,
        flip=(smi == ETAC), bond_px=48, font_px=22)
    f.text(x, 260, lab, 22, weight='bold')
    base = 380
    f.line(x - 110, base, x + 110, base, width=2)
    for k in range(lines):
        xx = x + (k - (lines - 1)/2)*26
        f.line(xx, base, xx, base - 90/lines, color=col, width=5)
    f.text(x, 418, sub, 23, color=col, weight='bold')
f.text(480, 466, 'n H on the neighbouring carbons → n + 1 lines', 20)
save('singlet-doublet', f)

# 8. Carbon-13 spectrum with every carbon labelled by its shift.
f = Figure('Carbon-13 NMR of U1 assigned to ethyl ethanoate', 'Schematic decoupled 13C NMR with four lines at 171, 60, 21 and 14 ppm; ethyl ethanoate drawn with each carbon labelled by its shift', 960, 440)
ax = Axis(f, 60, 900, 360, 200, 0, [200, 150, 100, 50, 0], 'chemical shift δ / ppm', size=17)
for ppm, col, h, dx in [(171, RED, 60, 0), (60, BLUE, 100, 0), (21, GREEN, 110, -16), (14, INK, 105, 16)]:
    f.line(ax.x(ppm), 360, ax.x(ppm), 360 - h, color=col, width=4)
    f.text(ax.x(ppm) + dx, 360 - h - 10, ppm, 20, color=col, weight='bold')
mol(f, ETAC, 340, 130, numbers={0: '21', 1: '171', 4: '60', 5: '14'}, highlight_atoms=[1, 2], flip=True, bond_px=70, font_px=32)
f.text(760, 90, 'C=O: 160–220 ppm', 20, color=RED).text(760, 122, 'C–O: 50–80 ppm', 20, color=BLUE).text(760, 154, 'alkyl C: 10–40 ppm', 20, color=INK)
save('carbon', f)

# 9. Symmetry: equivalent carbons share one signal.
f = Figure('Symmetry merges carbon signals', '1,4-dioxane has 4 carbons and 1 signal, pentan-3-one 5 carbons and 3 signals, isopropyl ethanoate 5 carbons and 4 signals; equivalent carbons share a letter', 960, 380)
for x, smi, notes, name, n in [
        (160, DIOXANE, {0: 'a', 1: 'a', 3: 'a', 4: 'a'}, '1,4-dioxane', '4 C → 1 signal'),
        (480, 'CCC(=O)CC', {0: 'a', 1: 'b', 2: 'c', 4: 'b', 5: 'a'}, 'pentan-3-one', '5 C → 3 signals'),
        (800, 'CC(=O)OC(C)C', {0: 'a', 1: 'b', 4: 'c', 5: 'd', 6: 'd'}, 'isopropyl ethanoate', '5 C → 4 signals')]:
    mol(f, smi, x, 150, numbers=notes, bond_px=55, font_px=32)
    f.text(x, 290, name, 22, weight='bold').text(x, 325, n, 22, color=BLUE)
save('symmetry', f)

# 10. Full assignment map for U1 (flip=True: acetyl CH3 top left, ethyl on the right).
f = Figure('Every observation assigned to ethyl ethanoate', 'Ethyl ethanoate with callouts: acetyl CH3 2.05 s 3H and 21 ppm; C=O IR 1740 and 171 ppm; O-CH2 4.12 q 2H and 60 ppm; CH3 1.25 t 3H and 14 ppm', 960, 440)
c = mol(f, ETAC, 480, 190, bond_px=80, font_px=32, highlight_atoms=[1, 2], flip=True)
callouts = [(0, (135, 60), (200, 85), GREEN, ['¹H 2.05 · s · 3H', '¹³C 21']),
            (2, (250, 380), (330, 330), RED, ['IR 1740 cm⁻¹ (C=O)', '¹³C 171']),
            (4, (640, 330), (600, 290), BLUE, ['¹H 4.12 · q · 2H', '¹³C 60']),
            (5, (850, 70), (790, 95), INK, ['¹H 1.25 · t · 3H', '¹³C 14'])]
for atom, (tx, ty), (lx, ly), col, lines in callouts:
    f.line(c[atom][0], c[atom][1], lx, ly, color=col, width=2, dash=5)
    for j, t in enumerate(lines):
        f.text(tx, ty + 28*j, t, 21, color=col, weight='bold')
f.text(850, 400, 'IR: no broad O–H', 20, color=INK)
save('assignment', f)

def compare(title, desc, name, obs, pred, smi, hi, ticks, extra=None):
    """Observed U1 spectrum above a candidate's predicted spectrum, same axis."""
    f = Figure(title, desc, 960, 440)
    for row, (label, sigs, y) in enumerate([('U1 observed', obs, 160), (name + ': predicted', pred, 350)]):
        ax = Axis(f, 40, 640, y, hi, 0, ticks if row else [], None, size=16)
        f.text(40, y - 125, label, 20, anchor='start', weight='bold', color=INK if row == 0 else RED)
        for ppm, lines, nh, col, lab, dy in sigs:
            if lines == 0:
                broad(ax, ppm, 22*nh, 12, col)
                top = ax.y - 22*nh
            else:
                top = multiplet(ax, ppm, lines, nh, col, k=24, gap=5, width=3)
            f.text(ax.x(ppm), top - 10 - dy, lab, 16, color=col, weight='bold')
    f.text(340, 432, "δ / ppm", 18)
    mol(f, smi, 795, 290, **(extra or {}))
    f.text(795, 420, name, 20, weight='bold')
    return f

OBS = [(4.12, 4, 2, INK, 'q 2H', 0), (2.05, 1, 3, INK, 's 3H', 0), (1.25, 3, 3, INK, 't 3H', 0)]
# 11. Butanoic acid predicts a different spectrum.
save('acid-test', compare(
    'U1 against the butanoic acid prediction',
    'Observed U1 signals at 4.12 q 2H, 2.05 s 3H and 1.25 t 3H compared with butanoic acid predicted at 11.5 broad 1H, 2.35 t 2H, 1.68 sextet 2H and 0.98 t 3H',
    'butanoic acid', OBS,
    [(11.5, 0, 1, RED, 'br 1H', 0), (2.35, 3, 2, RED, 't 2H', 0), (1.68, 6, 2, RED, '6 lines 2H', 36),
     (0.98, 3, 3, RED, 't 3H', 62)],
    BUAC, 12, [12, 10, 8, 6, 4, 2, 0], dict(colour=RED, highlight_atoms=[5], bond_px=44, font_px=24)))

# 12. 3-Hydroxybutan-2-one predicts a CH quartet of 1H and a methyl doublet.
save('acetoin-test', compare(
    'U1 against the 3-hydroxybutan-2-one prediction',
    'Observed U1 signals compared with 3-hydroxybutan-2-one predicted at 4.25 q 1H, about 3.6 broad OH 1H, 2.20 s 3H and 1.40 d 3H',
    '3-hydroxybutan-2-one', OBS,
    [(4.25, 4, 1, RED, 'q 1H', 0), (3.6, 0, 1, RED, 'OH', 0), (2.20, 1, 3, RED, 's 3H', 0), (1.40, 2, 3, RED, 'd 3H', 0)],
    ACETOIN, 5, [5, 4, 3, 2, 1, 0], dict(colour=RED, highlight_atoms=with_h(ACETOIN, [3]), explicit_h=True, bond_px=38, font_px=19)))

print('Built 12 spectroscopy figures.')
