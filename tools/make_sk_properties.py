"""Figures for the Day 2 properties and IR lecture (RDKit via molfig, plots via chemfig)."""
import math
from contextlib import contextmanager
from pathlib import Path

from rdkit import Chem
import molfig
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'properties'/'figures'
GREY = '#6b7f8a'
PALE = '#a2c8d8'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

# molfig parses SMILES with default settings, which strips explicit [H] atoms.
# Hydrogen-bond pictures need the O-H / N-H hydrogen as a real atom, so we
# temporarily hand molfig a parser that keeps the [H] written in the SMILES.
class _KeepH:
    def __getattr__(self, name):
        return getattr(Chem, name)
    @staticmethod
    def MolFromSmiles(smiles):
        p = Chem.SmilesParserParams()
        p.removeHs = False
        return Chem.MolFromSmiles(smiles, p)

@contextmanager
def keep_h():
    old = molfig.Chem
    molfig.Chem = _KeepH()
    try:
        yield
    finally:
        molfig.Chem = old

def placeh(fig, smiles, x, y, **kw):
    with keep_h():
        return place(fig, smiles, x, y, **kw)

def angle(a, b):
    """Direction from point a to point b, degrees, 0 = right, 90 = up."""
    return math.degrees(math.atan2(-(b[1] - a[1]), b[0] - a[0]))

def pairs(fig, c, atom, neighbours, n, distance=26, spread=None):
    """Draw n lone pairs on `atom`, spread evenly through the widest gap between its bonds."""
    at = c[atom]
    angs = sorted(angle(at, c[j]) % 360 for j in neighbours)
    if not angs:
        start, gap = 90 - 180, 360
    else:
        gaps = [((angs[(k + 1) % len(angs)] - angs[k]) % 360 or 360, angs[k]) for k in range(len(angs))]
        gap, start = max(gaps)
    if spread is not None and n == 2:
        mid = start + gap/2
        dirs = [mid - spread, mid + spread]
    else:
        dirs = [start + gap*(k + 1)/(n + 1) for k in range(n)]
    for a in dirs:
        lone_pair(fig, at, a, distance=distance)

def dash(fig, a, b, shorten=(18, 22), color=BLUE):
    """Dashed hydrogen-bond contact from a to b, pulled back from both labels."""
    d = math.hypot(b[0] - a[0], b[1] - a[1])
    ux, uy = (b[0] - a[0])/d, (b[1] - a[1])/d
    fig.line(a[0] + ux*shorten[0], a[1] + uy*shorten[0], b[0] - ux*shorten[1], b[1] - uy*shorten[1],
             color=color, width=3, dash='7 6')

def arrow(fig, a, b, color=BLUE, width=4):
    fig.line(a[0], a[1], b[0], b[1], color=color, width=width, arrow=True)

# 1. Polar bonds: electronegativity differences in ethanol.
f = Figure('Bond polarity in ethanol', 'Ethanol with partial charges on C, O and H, and a Pauling electronegativity scale', 960, 400)
c = placeh(f, 'CCO[H]', 260, 170, bond_px=95)
f.text(c[1][0], c[1][1] - 34, 'δ+', 28, color=RED)
f.text(c[2][0] + 4, c[2][1] + 52, 'δ−', 28, color=BLUE)
f.text(c[3][0] + 6, c[3][1] - 30, 'δ+', 28, color=RED)
f.text(250, 320, 'C–O and O–H are polar;', 22).text(250, 352, 'C–H is nearly nonpolar', 22)
x0, x1 = 580, 900
f.line(x0, 210, x1, 210, width=3)
for sym, en, col in [('H', 2.20, INK), ('C', 2.55, INK), ('N', 3.04, BLUE), ('O', 3.44, RED)]:
    x = x0 + (en - 2.0)/(3.6 - 2.0)*(x1 - x0)
    f.circle(x, 210, 9, color=col, fill=col)
    f.text(x, 180, sym, 30, color=col, weight='bold').text(x, 250, f'{en:.2f}', 21)
f.text(740, 110, 'electronegativity (Pauling)', 23, weight='bold')
f.line(600, 300, 880, 300, color=BLUE, arrow=True)
f.text(740, 335, 'pulls bonding electrons harder', 21, color=BLUE)
save('polar-bonds', f)

# 2. Partial charge versus formal charge.
f = Figure('Partial charge and formal charge are different ideas', 'Ethanol oxygen: delta minus, formal charge 0; ethoxide oxygen: formal charge minus one', 960, 400)
c = placeh(f, 'CCO[H]', 250, 150, bond_px=85)
pairs(f, c, 2, [1, 3], 2)
f.text(c[2][0] - 40, c[2][1] + 44, 'δ−', 26, color=BLUE)
c = place(f, 'CC[O-]', 700, 150, bond_px=85)
pairs(f, c, 2, [1], 3, distance=28)
f.text(250, 285, 'ethanol', 25, weight='bold').text(250, 322, 'O: 2 bonds, 2 lone pairs', 22, color=BLUE)
f.text(250, 356, 'formal charge 6 − 4 − 2 = 0', 22, color=BLUE)
f.text(700, 285, 'ethoxide', 25, weight='bold').text(700, 322, 'O: 1 bond, 3 lone pairs', 22, color=BLUE)
f.text(700, 356, 'formal charge 6 − 6 − 1 = −1', 22, color=BLUE)
save('partial-formal', f)

# 3. Bond dipoles add as vectors.
f = Figure('Bond dipoles add as vectors', 'Linear carbon dioxide: two bond dipoles cancel. Bent water: they add to a net dipole', 960, 400)
c = place(f, 'O=C=O', 240, 190, bond_px=85)
xl, xc, xr = sorted(p[0] for p in c.values())
arrow(f, (xc - 12, 130), (xl + 8, 130))
arrow(f, (xc + 12, 130), (xr - 8, 130))
f.text(xl, 100, 'δ−', 24, color=BLUE).text(xc, 100, 'δ+', 24, color=RED).text(xr, 100, 'δ−', 24, color=BLUE)
f.text(240, 290, 'CO₂ (linear)', 25, weight='bold').text(240, 330, 'bond dipoles cancel · net 0', 22, color=BLUE)
c = place(f, 'O', 680, 160, bond_px=95, explicit_h=True, rotate=180)
o = c[0]
for h in (c[1], c[2]):
    # bond dipole drawn parallel to O-H, offset to the outside of the angle, pointing at O
    dx, dy = o[0] - h[0], o[1] - h[1]
    d = math.hypot(dx, dy)
    nx, ny = -dy/d, dx/d
    if nx*(h[0] - o[0]) < 0:
        nx, ny = -nx, -ny
    a = (h[0] + dx*0.2 + nx*24, h[1] + dy*0.2 + ny*24)
    b = (h[0] + dx*0.75 + nx*24, h[1] + dy*0.75 + ny*24)
    arrow(f, a, b)
hy = c[1][1]
arrow(f, (o[0], hy + 40), (o[0], o[1] + 22), color=RED, width=5)
f.text(o[0] + 14, hy + 34, 'net dipole', 21, color=RED, anchor='start')
f.text(o[0] + 30, o[1] - 14, 'δ−', 24, color=BLUE, anchor='start')
f.text(680, 300, 'H₂O (bent)', 25, weight='bold').text(680, 340, 'bond dipoles add · polar', 22, color=BLUE)
save('dipoles', f)

# 4. Three molecules of the same size, three boiling points.
f = Figure('Same size, different forces', 'Butane, propanone and propan-1-ol: similar molar mass, boiling points minus 1, 56 and 97 degrees Celsius', 960, 440)
items = [('CCCC', 'butane', '58 g/mol · bp −1 °C', 'dispersion only', {}),
         ('CC(C)=O', 'propanone', '58 g/mol · bp 56 °C', '+ dipole–dipole', dict(highlight_atoms=[1, 3])),
         ('CCCO', 'propan-1-ol', '60 g/mol · bp 97 °C', '+ hydrogen bonds', dict(highlight_atoms=[3]))]
for i, (smi, name, sub, force, kw) in enumerate(items):
    x = 160 + 320*i
    place(f, smi, x, 140, bond_px=70, **kw)
    f.text(x, 270, name, 25, weight='bold').text(x, 307, sub, 22, color=BLUE)
    f.text(x, 345, force, 22, color=[GREY, INK, RED][i])
f.line(90, 395, 870, 395, color=BLUE, arrow=True)
f.text(480, 428, 'stronger attraction between molecules → higher boiling point', 21, color=BLUE)
save('same-size', f)

# 5. Hydrogen bonds: donor H, acceptor lone pair.
f = Figure('A hydrogen bond needs a donor and an acceptor', 'Ethanol O-H donating to water oxygen; water O-H donating to dimethyl ether oxygen', 960, 400)
c = placeh(f, 'CCO[H]', 130, 150, bond_px=70)
w = place(f, 'O', 380, 175, bond_px=65, explicit_h=True, rotate=-90)
pairs(f, w, 0, [1, 2], 2, spread=38)
dash(f, c[3], w[0])
f.text(240, 290, 'ethanol donates · water accepts', 22, weight='bold')
w = place(f, 'O', 600, 150, bond_px=65, explicit_h=True, rotate=-90)
e = place(f, 'COC', 830, 190, bond_px=70, rotate=0)
pairs(f, e, 1, [0, 2], 2, spread=38)
h = max((w[1], w[2]), key=lambda p: p[1])
dash(f, h, e[1])
f.text(720, 290, 'water donates · ether accepts', 22, weight='bold')
f.text(480, 360, 'donor: H on O or N    acceptor: a lone pair on O or N', 23, color=BLUE)
save('hbond', f)

# 6. Donor and acceptor inventory of oxygen compounds.
f = Figure('Donor and acceptor inventory', 'Ethanol, dimethyl ether, propanone and methyl acetate with acceptor oxygens shaded and donor counts', 960, 400)
items = [('CCO[H]', 'ethanol', 1, 1, [2]),
         ('COC', 'dimethyl ether', 0, 1, [1]),
         ('CC(C)=O', 'propanone', 0, 1, [3]),
         ('CC(=O)OC', 'methyl acetate', 0, 2, [2, 3])]
for i, (smi, name, don, acc, hl) in enumerate(items):
    x = 120 + 240*i
    placeh(f, smi, x, 150, bond_px=58, font_px=26, highlight_atoms=hl, highlight_bonds=[])
    f.text(x, 270, name, 23, weight='bold')
    f.text(x, 302, f'donor O–H: {don}', 21, color=RED if don else GREY)
    f.text(x, 330, f'acceptor O: {acc}', 21, color=BLUE)
f.text(480, 378, 'Every O can accept. Only an O–H can donate.', 23)
save('inventory', f)

# 7. Amines: count N-H bonds and the lone pair.
f = Figure('Amine donors and acceptors', 'Methylamine, dimethylamine, trimethylamine and methylammonium with N-H bonds and lone pairs', 960, 400)
items = [('CN([H])[H]', 'methylamine', '2 N–H · lone pair', 1),
         ('CN([H])C', 'dimethylamine', '1 N–H · lone pair', 1),
         ('CN(C)C', 'trimethylamine', '0 N–H · lone pair', 1),
         ('C[N+]([H])([H])[H]', 'methylammonium', '3 N–H · no lone pair', 0)]
for i, (smi, name, sub, lp) in enumerate(items):
    x = 120 + 240*i
    c = placeh(f, smi, x, 160, bond_px=58, font_px=26)
    if lp:
        nb = [a.GetIdx() for a in _KeepH.MolFromSmiles(smi).GetAtomWithIdx(1).GetNeighbors()]
        pairs(f, c, 1, nb, 1, distance=28)
    f.text(x, 285, name, 22, weight='bold').text(x, 322, sub, 20, color=RED if not lp else BLUE)
f.text(480, 370, 'Protonation turns the lone pair into an N–H bond.', 23)
save('amines', f)

# 8. Back to the hook: ethanol vs dimethyl ether in the liquid.
f = Figure('Why ethanol boils 100 °C higher than dimethyl ether', 'Two ethanol molecules joined by an O-H to O hydrogen bond; two dimethyl ether molecules with no donor', 960, 420)
a = placeh(f, 'CCO[H]', 130, 120, bond_px=62)
b = placeh(f, '[H]OCC', 340, 235, bond_px=62, rotate=0)
pairs(f, b, 1, [0, 2], 2, spread=45, distance=24)
dash(f, a[3], b[1])
f.text(240, 330, 'ethanol · bp 78 °C', 25, weight='bold')
f.text(240, 368, 'O–H···O hydrogen bonds hold the liquid', 21, color=BLUE)
place(f, 'COC', 640, 120, bond_px=62)
place(f, 'COC', 820, 230, bond_px=62)
f.text(730, 330, 'dimethyl ether · bp −24 °C', 25, weight='bold')
f.text(730, 368, 'acceptor O but no O–H donor', 21, color=BLUE)
f.line(480, 60, 480, 380, color='#d6e2e8', width=2)
save('boiling', f)

# 9. Size and shape: dispersion.
f = Figure('Chain length and branching change the boiling point', 'Boiling points of straight-chain alcohols from methanol to hexan-1-ol, and pentane versus 2,2-dimethylpropane', 960, 420)
data = [(1, 65), (2, 78), (3, 97), (4, 118), (5, 138), (6, 157)]
X0, Y0, XW, YH = 90, 330, 380, 250
f.line(X0, Y0, X0 + XW, Y0, width=2).line(X0, Y0, X0, Y0 - YH - 10, width=2)
for t in (50, 100, 150):
    y = Y0 - t/170*YH
    f.line(X0 - 6, y, X0, y, width=2).text(X0 - 12, y + 7, str(t), 18, anchor='end')
pts = []
for n, bp in data:
    x = X0 + 30 + (n - 1)*(XW - 60)/5
    y = Y0 - bp/170*YH
    pts.append((x, y))
    f.text(x, Y0 + 28, str(n), 19)
f.path('M ' + ' L '.join(f'{x:.0f},{y:.0f}' for x, y in pts), color=BLUE, width=3)
for x, y in pts:
    f.circle(x, y, 7, color=BLUE, fill=BLUE)
f.text(X0 + XW/2, Y0 + 62, 'carbons in CₙH₂ₙ₊₁OH', 21)
f.text(X0 + 10, 50, 'bp (°C) of alcohols', 21, anchor='start', weight='bold')
place(f, 'CCCCC', 610, 150, bond_px=50)
place(f, 'CC(C)(C)C', 825, 150, bond_px=50)
f.text(610, 255, 'pentane', 22, weight='bold').text(610, 290, 'bp 36 °C', 22, color=BLUE)
f.text(825, 255, '2,2-dimethylpropane', 20, weight='bold').text(825, 290, 'bp 10 °C', 22, color=BLUE)
f.text(740, 345, 'both C₅H₁₂: the compact ball', 21).text(740, 375, 'has less contact area', 21)
save('size-shape', f)

# 10. Solubility: polar head against hydrocarbon tail.
f = Figure('A polar head and a hydrocarbon tail compete in water', 'Methanol, butan-1-ol and hexan-1-ol with the OH shaded and their water solubility', 960, 400)
items = [('CO', 'methanol', 'mixes completely'),
         ('CCCCO', 'butan-1-ol', '≈ 7 g per 100 mL'),
         ('CCCCCCO', 'hexan-1-ol', '≈ 0.6 g per 100 mL')]
for i, (smi, name, sub) in enumerate(items):
    x = [100, 375, 740][i]
    n = len(smi)
    place(f, smi, x, 160, bond_px=50, highlight_atoms=[n - 1], highlight_bonds=[])
    f.text(x, 280, name, 25, weight='bold').text(x, 318, sub, 22, color=BLUE)
f.text(480, 370, 'One OH can keep about three or four carbons in water, not six.', 22)
save('solubility', f)

# 11. Charge changes solubility: benzoic acid and benzoate.
f = Figure('A charge makes a molecule far more water soluble', 'Benzoic acid, poorly soluble, and sodium benzoate, very soluble, with the carboxyl group shaded', 960, 400)
place(f, 'OC(=O)c1ccccc1', 250, 160, bond_px=60, highlight_atoms=[0, 1, 2], rotate=0)
place(f, '[O-]C(=O)c1ccccc1', 690, 160, bond_px=60, highlight_atoms=[0, 1, 2], rotate=0)
f.text(880, 215, 'Na⁺', 28, color=INK)
f.text(250, 290, 'benzoic acid', 25, weight='bold').text(250, 328, '≈ 0.3 g per 100 mL water', 22, color=BLUE)
f.text(690, 290, 'sodium benzoate', 25, weight='bold').text(690, 328, '> 50 g per 100 mL water', 22, color=BLUE)
f.line(400, 160, 530, 160, color=BLUE, arrow=True).text(465, 140, 'base', 20, color=BLUE)
f.text(480, 375, 'Water surrounds an ion strongly: ion–dipole hydration.', 22)
save('ionised', f)

# 12. IR absorption needs a changing dipole: CO2 stretches.
f = Figure('A vibration absorbs IR only if it changes the dipole', 'Carbon dioxide symmetric stretch: no dipole change, no IR band. Asymmetric stretch: dipole changes, absorbs near 2349 per centimetre', 960, 420)
for x, mode, sub, col, arrows in [
        (250, 'symmetric stretch', 'dipole stays zero → no IR band', GREY, [(-1, -1), (0, 0), (1, 1)]),
        (710, 'asymmetric stretch', 'dipole changes → absorbs ≈ 2349 cm⁻¹', RED, [(-1, 1), (0, -1), (1, 1)])]:
    c = place(f, 'O=C=O', x, 170, bond_px=85)
    for (k, d), (ax, ay) in zip(arrows, sorted(c.values())):
        if d == 0:
            continue
        arrow(f, (ax - 22*d, ay - 50), (ax + 22*d, ay - 50), color=col, width=4)
    f.text(x, 280, mode, 25, weight='bold').text(x, 320, sub, 22, color=col)
f.text(480, 385, 'Bonds behave like springs: stiffer bonds and lighter atoms vibrate faster.', 22)
save('ir-co2', f)

# Simulated spectra shared by the next figures.
def trace(fig, x0, y0, w, h, bands, color=BLUE, absorbance=False, lo=500, hi=4000):
    """bands: (centre, depth 0-1, width) in cm-1; draws a schematic trace."""
    pts = []
    for i in range(401):
        wn = hi - i*(hi - lo)/400
        a = sum(d*math.exp(-((wn - c0)/wd)**2) for c0, d, wd in bands)
        a = min(a, 0.95)
        v = a if absorbance else 1 - a
        pts.append((x0 + i*w/400, y0 + h - v*h*0.9 - (0 if absorbance else 0.05*h)))
    fig.path('M ' + ' L '.join(f'{x:.1f},{y:.1f}' for x, y in pts), color=color, width=3)

def wn_x(x0, w, wn, lo=500, hi=4000):
    return x0 + (hi - wn)/(hi - lo)*w

def axis(fig, x0, y0, w, h, ylabel, ticks=True):
    fig.line(x0, y0 + h, x0 + w, y0 + h, width=2).line(x0, y0 + h, x0, y0 - 5, width=2)
    if ticks:
        for wn in (4000, 3000, 2000, 1000):
            x = wn_x(x0, w, wn)
            fig.line(x, y0 + h, x, y0 + h + 6, width=2).text(x, y0 + h + 26, str(wn), 17)
    fig.text(x0 - 12, y0 + h/2, ylabel, 18, anchor='end')

KETONE = [(2950, 0.25, 50), (1715, 0.85, 30), (1365, 0.3, 25), (1220, 0.35, 30)]

# 13. Transmittance and absorbance display the same absorption.
f = Figure('Transmittance and absorbance show the same band', 'The same carbonyl absorption near 1715 per centimetre as a dip on a transmittance plot and a peak on an absorbance plot. Schematic, not measured data', 960, 420)
for x0, absb, title in [(70, False, 'transmittance: absorption points down'), (540, True, 'absorbance: absorption points up')]:
    w, h = 360, 220
    f.line(x0, 70 + h, x0 + w, 70 + h, width=2).line(x0, 70 + h, x0, 60, width=2)
    trace(f, x0, 70, w, h, KETONE, color=GREEN if absb else BLUE, absorbance=absb)
    for wn in (4000, 3000, 2000, 1000):
        x = wn_x(x0, w, wn)
        f.line(x, 70 + h, x, 76 + h, width=2).text(x, 70 + h + 26, str(wn), 17)
    xb = wn_x(x0, w, 1715)
    f.text(xb + 10, 70 + h - 20 if not absb else 95, 'C=O', 21, color=RED, anchor='start')
    f.text(x0 + w/2, 40, title, 22, weight='bold')
    f.text(x0 + w/2, 70 + h + 60, 'wavenumber (cm⁻¹), decreasing →', 19)
f.text(480, 410, 'schematic', 16, color=GREY)
save('ir-axes', f)

# 14. Where the useful bands are.
f = Figure('IR regions worth knowing', 'Approximate IR ranges: alcohol O-H 3200 to 3600, acid O-H 2500 to 3300, N-H 3300 to 3500, C-H 2850 to 3000, C=O 1650 to 1800 per centimetre, fingerprint below 1500', 960, 400)
X0, W = 250, 650
rows = [('alcohol O–H', 3200, 3600, 'broad', PALE),
        ('acid O–H', 2500, 3300, 'very broad', PALE),
        ('N–H', 3300, 3500, 'medium', '#b8d8c9'),
        ('C–H', 2850, 3000, 'almost always', '#d6e2e8'),
        ('C=O', 1650, 1800, 'strong, sharp', '#efc1c5'),
        ('fingerprint', 500, 1500, 'compare, don’t assign', '#e9edf0')]
for i, (name, lo, hi, note, col) in enumerate(rows):
    y = 40 + 46*i
    f.text(X0 - 15, y + 26, name, 21, anchor='end', weight='bold')
    xa, xb = wn_x(X0, W, hi), wn_x(X0, W, lo)
    f.rect(xa, y + 6, xb - xa, 30, fill=col)
    tx = xb + 10 if xb < X0 + W - 200 else xa - 10
    f.text(tx, y + 28, note, 18, color=INK, anchor='start' if xb < X0 + W - 200 else 'end')
yax = 40 + 46*len(rows) + 10
f.line(X0, yax, X0 + W, yax, width=2)
for wn in (4000, 3500, 3000, 2500, 2000, 1500, 1000, 500):
    x = wn_x(X0, W, wn)
    f.line(x, yax, x, yax + 6, width=2).text(x, yax + 28, str(wn), 17)
f.text(X0 + W/2, yax + 60, 'wavenumber (cm⁻¹)', 19)
save('ir-regions', f)

# 15. Three simulated spectra: alcohol, ketone, acid.
f = Figure('Ethanol, propanone and acetic acid in the IR', 'Schematic transmittance spectra: ethanol with a broad O-H band and no C=O; propanone with a strong C=O and no O-H; acetic acid with a very broad O-H and a C=O', 960, 490)
specs = [('CCO', 'ethanol', [(3350, 0.6, 180), (2950, 0.35, 40), (1050, 0.6, 40), (1380, 0.15, 30)], [(3350, 'broad O–H')]),
         ('CC(C)=O', 'propanone', KETONE, [(1715, 'C=O')]),
         ('CC(=O)O', 'acetic acid', [(3000, 0.55, 380), (1710, 0.85, 30), (1290, 0.45, 30), (940, 0.2, 40)], [(3000, 'very broad O–H'), (1710, 'C=O')])]
X0, W, H = 300, 620, 105
for i, (smi, name, bands, labels) in enumerate(specs):
    y0 = 25 + 140*i
    place(f, smi, 110, y0 + 50, bond_px=42, font_px=22)
    f.text(110, y0 + 118, name, 21, weight='bold')
    f.line(X0, y0 + H, X0 + W, y0 + H, width=1.5, color='#c4d3db')
    trace(f, X0, y0, W, H, bands)
    for wn, lab in labels:
        f.text(wn_x(X0, W, wn), y0 + H + 2, lab, 18, color=RED)
for wn in (4000, 3000, 2000, 1000):
    x = wn_x(X0, W, wn)
    f.text(x, 25 + 140*2 + H + 34, str(wn), 17)
f.text(X0 + W/2, 478, 'wavenumber (cm⁻¹) · schematic transmittance traces', 18, color=GREY)
save('ir-spectra', f)

# 16. Choosing among candidates.
f = Figure('Which candidates fit: strong C=O, no O-H', 'Ethanol, dimethyl ether, propanone and propanal with predicted O-H and C=O bands; only the two carbonyl compounds fit', 960, 420)
items = [('CCO', 'ethanol', 'O–H yes · C=O no', False, [2]),
         ('COC', 'dimethyl ether', 'O–H no · C=O no', False, []),
         ('CC(C)=O', 'propanone', 'O–H no · C=O yes', True, [1, 3]),
         ('CCC=O', 'propanal', 'O–H no · C=O yes', True, [2, 3])]
for i, (smi, name, pred, ok, hl) in enumerate(items):
    x = 120 + 240*i
    if ok:
        f.rect(x - 110, 30, 220, 330, fill='#e4f0f4')
    place(f, smi, x, 140, bond_px=50, font_px=24, highlight_atoms=hl)
    f.text(x, 265, name, 23, weight='bold').text(x, 302, pred, 20, color=BLUE)
    f.text(x, 342, 'fits' if ok else 'ruled out', 22, color=GREEN if ok else RED, weight='bold')
f.text(480, 400, 'Observation: strong band near 1715 cm⁻¹, no broad O–H. Two candidates survive.', 21)
save('ir-candidates', f)

print('Built 16 property and IR figures.')
