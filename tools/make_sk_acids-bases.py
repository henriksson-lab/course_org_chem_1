"""Skeletal-structure figures for the Day 3 acids-and-bases lecture (RDKit via molfig)."""
import math
from pathlib import Path
from rdkit import Chem
import molfig
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'acids-bases'/'figures'

# molfig parses SMILES with RDKit defaults, which strip explicit [H] atoms. Proton-transfer
# drawings need the transferred H as a real atom, so parse with removeHs=False here.
class _KeepH:
    def __getattr__(self, name):
        return getattr(Chem, name)
    @staticmethod
    def MolFromSmiles(smiles):
        p = Chem.SmilesParserParams()
        p.removeHs = False
        return Chem.MolFromSmiles(smiles, p)
molfig.Chem = _KeepH()

GREY = '#5b6b75'
PALE = '#e4f0f4'
PALE_GREEN = '#dcefe6'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

def eq_arrow(f, x1, x2, y, color=BLUE):
    """Equilibrium half-arrows (⇌) drawn as two offset lines with barbs."""
    f.line(x1, y - 6, x2, y - 6, color=color, width=3)
    f.path(f'M {x2} {y - 6} L {x2 - 14} {y - 16}', color=color, width=3)
    f.line(x1, y + 6, x2, y + 6, color=color, width=3)
    f.path(f'M {x1} {y + 6} L {x1 + 14} {y + 16}', color=color, width=3)

def mid(a, b, t=0.5):
    return (a[0] + (b[0] - a[0])*t, a[1] + (b[1] - a[1])*t)

# 1. Aspirin: the same drug in stomach and blood.
f = Figure('Aspirin in the stomach and in the blood', '', 960, 430)
place(f, 'CC(=O)Oc1ccccc1C(=O)O', 250, 170, bond_px=50, font_px=26, highlight_atoms=[10, 11, 12])
place(f, 'CC(=O)Oc1ccccc1C(=O)[O-]', 710, 170, bond_px=50, font_px=26, highlight_atoms=[10, 11, 12])
f.text(250, 330, 'stomach, pH ≈ 2', 25, weight='bold').text(250, 368, 'about 97 % neutral HA', 22, color=BLUE)
f.text(710, 330, 'blood, pH 7.4', 25, weight='bold').text(710, 368, 'more than 99.9 % A⁻', 22, color=BLUE)
f.text(480, 412, 'aspirin, pKa ≈ 3.5', 22, color=GREY)
save('aspirin', f)

# 2. Proton transfer from ethanoic acid to hydroxide, with both curved arrows.
f = Figure('Ethanoic acid and hydroxide: two curved arrows', '', 960, 400)
a = place(f, 'CC(=O)O[H]', 130, 160, bond_px=56, rotate=-30)
O_a, H_a = a[3], a[4]
b = place(f, '[H][O-]', 410, H_a[1], bond_px=56)
O_b = b[1]
lone_pair(f, O_b, 90); lone_pair(f, O_b, 270); lone_pair(f, O_b, 180, distance=26)
curved_arrow(f, (O_b[0] - 32, O_b[1] - 4), (H_a[0] + 6, H_a[1] - 10), bend=0.45)
curved_arrow(f, mid(O_a, H_a), (O_a[0] + 4, O_a[1] + 30), bend=-0.7, shorten=4)
f.text(300, H_a[1] + 10, '+', 32)
f.line(490, 190, 565, 190, color=BLUE, arrow=True)
place(f, 'CC(=O)[O-]', 650, 160, bond_px=52)
f.text(758, 200, '+', 32)
place(f, '[H]O[H]', 870, 195, bond_px=46)
f.text(250, 320, 'acid + base', 22, color=GREY)
f.text(760, 320, 'conjugate base + conjugate acid', 22, color=GREY)
f.text(250, 372, 'arrow 1: lone pair → H', 23, color=RED).text(700, 372, 'arrow 2: O–H bond pair → O', 23, color=RED)
save('proton-transfer', f)

# 3. Protonating methylamine with a generic acid H–A.
f = Figure('Methylamine accepts a proton from H–A', '', 960, 360)
m = place(f, 'C[NH2]', 130, 170, bond_px=64)
N = m[1]
lone_pair(f, N, 0, distance=62)
hx, hy = 310, N[1]
f.text(hx, hy + 11, 'H', 32).line(hx + 16, hy, hx + 84, hy, width=3).text(hx + 104, hy + 11, 'A', 32)
curved_arrow(f, (N[0] + 70, N[1] - 6), (hx - 2, hy - 14), bend=0.6, shorten=6)
curved_arrow(f, (hx + 50, hy + 4), (hx + 104, hy + 18), bend=-0.8, shorten=4)
f.line(480, 175, 560, 175, color=BLUE, arrow=True)
place(f, 'C[NH3+]', 680, 170, bond_px=64)
f.text(800, 185, '+', 32).text(865, 185, 'A⁻', 32)
f.text(160, 270, 'N: 3 bonds + 1 lone pair', 22, color=GREY).text(160, 300, 'neutral', 22, color=BLUE)
f.text(680, 270, 'N: 4 bonds, no lone pair', 22, color=GREY).text(680, 300, 'charge +1', 22, color=BLUE)
f.text(480, 340, 'total charge: 0 before, (+1) + (−1) = 0 after', 22, color=RED)
save('amine-protonation', f)

# 4. Conjugate pairs: remove one H⁺, lower the charge by one.
f = Figure('Four acids and their conjugate bases', '', 960, 470)
pairs = [
    ('CC(=O)O', 'CC(=O)[O-]', 'ethanoic acid', 'ethanoate', [3]),
    ('CC[OH]', 'CC[O-]', 'ethanol', 'ethoxide', [2]),
    ('C[NH3+]', 'C[NH2]', 'methylammonium', 'methylamine', [1]),
    ('[NH4+]', 'N', 'ammonium', 'ammonia', [0]),
]
for i, (acid, base, an, bn, hl) in enumerate(pairs):
    x = 120 + 240*i
    place(f, acid, x, 85, bond_px=44, font_px=26, highlight_atoms=hl)
    f.text(x, 190, an, 21, weight='bold')
    f.line(x, 205, x, 255, color=BLUE, arrow=True).text(x + 12, 237, '− H⁺', 21, color=BLUE, anchor='start')
    place(f, base, x, 315, bond_px=44, font_px=26, highlight_atoms=hl)
    f.text(x, 415, bn, 21, weight='bold')
f.text(480, 455, 'one H fewer, charge one lower: 0 → −1 or +1 → 0', 22, color=RED)
save('conjugate-pairs', f)

# 5. Water as base and as acid.
f = Figure('Water accepts a proton from HCl and donates one to ammonia', '', 960, 360)
y1, y2 = 105, 255
for y, left, right1, right2, role, col in (
        (y1, 'HCl', '[H][OH+][H]', 'Cl⁻', 'water = base', BLUE),
        (y2, 'NH₃', '[H][O-]', 'NH₄⁺', 'water = acid', GREEN)):
    f.text(70, y + 10, left, 30).text(140, y + 10, '+', 30)
    place(f, '[H]O[H]', 225, y, bond_px=46)
    eq_arrow(f, 310, 390, y)
    if y == y1:
        place(f, right1, 480, y - 10, bond_px=46)
        f.text(575, y + 10, '+', 30).text(635, y + 10, right2, 30)
    else:
        f.text(450, y + 10, right2, 30).text(525, y + 10, '+', 30)
        place(f, right1, 600, y, bond_px=46)
    f.text(840, y + 8, role, 24, weight='bold', color=col)
f.text(480, 340, 'H₃O⁺/H₂O and H₂O/OH⁻ are both conjugate pairs', 22, color=GREY)
save('water-both', f)

# 6. pKa scale with structures.
f = Figure('Approximate aqueous pKa values on a log scale', '', 960, 470)
x0, x1, ay = 100, 900, 240
def px(pka): return x0 + (x1 - x0)*pka/17
f.line(x0 - 20, ay, x1 + 20, ay, width=3)
for t in range(0, 18, 2):
    f.line(px(t), ay - 8, px(t), ay + 8, width=2).text(px(t), ay + 30, t, 18, color=GREY)
f.text(x1 + 22, ay + 6, 'pKa', 19, color=GREY, anchor='start')
acids = [
    ('ClCC(=O)O', 2.9, 'chloroethanoic acid', 'up'),
    ('CC(=O)O', 4.8, 'ethanoic acid', 'down'),
    ('[NH4+]', 9.2, 'ammonium', 'up'),
    ('C[NH3+]', 10.6, 'methylammonium', 'down'),
    ('CCO', 16.0, 'ethanol', 'up'),
]
for smi, pka, name, side in acids:
    x = px(pka)
    f.circle(x, ay, 7, color=RED, fill=RED, width=1)
    if side == 'up':
        place(f, smi, x, 110, bond_px=44, font_px=25)
        f.line(x, 165, x, ay - 12, color=GREY, width=1.5, dash='4 4')
        f.text(x, 32, f'{name}  {pka:g}', 19, weight='bold')
    else:
        place(f, smi, x, 345, bond_px=44, font_px=25)
        f.text(x, 428, f'{name}  {pka:g}', 19, weight='bold')
f.line(260, 455, 90, 455, color=BLUE, arrow=True).text(270, 462, 'stronger acid', 19, color=RED, anchor='start')
f.line(700, 455, 870, 455, color=BLUE, arrow=True).text(690, 462, 'weaker acid', 19, color=BLUE, anchor='end')
save('pka-scale', f)

# 7. Ethoxide versus ethanoate: where does the negative charge sit?
f = Figure('The negative charge in ethoxide and ethanoate', '', 960, 400)
place(f, 'CC[O-]', 170, 160, bond_px=64, highlight_atoms=[2])
f.text(170, 280, 'ethoxide', 25, weight='bold').text(170, 318, 'charge on one O', 22, color=BLUE)
f.text(170, 360, 'ethanol pKa ≈ 16', 22, color=GREY)
c1 = place(f, 'CC(=O)[O-]', 500, 160, bond_px=64, highlight_atoms=[2, 3])
f.text(625, 175, '↔', 40, color=BLUE)
c2 = place(f, 'CC([O-])=O', 760, 160, bond_px=64, highlight_atoms=[2, 3])
f.text(630, 280, 'ethanoate', 25, weight='bold').text(630, 318, 'charge shared by two equivalent O', 22, color=BLUE)
f.text(630, 360, 'ethanoic acid pKa ≈ 4.8', 22, color=GREY)
f.line(320, 40, 320, 360, color='#c9d6dd', width=2)
save('ethoxide-ethanoate', f)

# 8. Resonance in ethanoate: curved arrows relate the two contributors.
f = Figure('Two resonance contributors of ethanoate', '', 960, 380)
r1 = place(f, 'CC(=O)[O-]', 230, 170, bond_px=72)
C, Od, Os = r1[1], r1[2], r1[3]
lone_pair(f, Os, -90, distance=30)
curved_arrow(f, (Os[0] - 4, Os[1] + 34), mid(C, Os), bend=-0.7, shorten=6)
curved_arrow(f, (C[0] + 10, (C[1] + Od[1])/2 + 6), (Od[0] + 20, Od[1] + 4), bend=0.7, shorten=2)
f.text(480, 190, '↔', 52, color=BLUE)
r2 = place(f, 'CC([O-])=O', 730, 170, bond_px=72)
f.text(480, 305, 'only electrons move: every nucleus stays put', 24, weight='bold')
f.text(480, 345, 'real ion: both C–O bonds the same length (about 126 pm)', 22, color=BLUE)
save('resonance', f)

# 9. Inductive effect of chlorine on ethanoic acid.
f = Figure('Chlorine substituents lower the pKa of ethanoic acid', '', 960, 380)
items = [
    ('CC(=O)O', 'ethanoic acid', '4.8', []),
    ('ClCC(=O)O', 'chloroethanoic acid', '2.9', [0]),
    ('OC(=O)C(Cl)(Cl)Cl', 'trichloroethanoic acid', '0.7', [4, 5, 6]),
]
for i, (smi, name, pka, hl) in enumerate(items):
    x = 160 + 320*i
    place(f, smi, x, 140, bond_px=58, highlight_atoms=hl)
    f.text(x, 270, name, 23, weight='bold').text(x, 308, f'pKa {pka}', 24, color=RED)
f.text(480, 360, 'Cl pulls electron density through σ bonds and stabilises the carboxylate', 21, color=GREY)
save('induction', f)

# 10. Predicting the position of an acid-base equilibrium.
f = Figure('Ethanoic acid and methylamine: equilibrium lies to the right', '', 960, 400)
y = 130
place(f, 'CC(=O)O', 120, y, bond_px=55, highlight_atoms=[3])
f.text(225, y + 12, '+', 30)
place(f, 'C[NH2]', 305, y + 15, bond_px=55)
eq_arrow(f, 390, 480, y + 15)
place(f, 'CC(=O)[O-]', 590, y, bond_px=55)
f.text(690, y + 12, '+', 30)
place(f, 'C[NH3+]', 790, y + 15, bond_px=55, highlight_atoms=[1])
f.text(120, 245, 'acid, pKa 4.8', 23, color=RED, weight='bold')
f.text(305, 245, 'base', 23, color=GREY)
f.text(590, 245, 'base', 23, color=GREY)
f.text(790, 245, 'acid, pKa 10.6', 23, color=BLUE, weight='bold')
f.parts.append('<text x="480" y="315" font-family="DejaVu Sans, sans-serif" font-size="28" fill="%s" '
               'text-anchor="middle" font-weight="bold">K ≈ 10<tspan dy="-12" font-size="19">10.6 − 4.8</tspan>'
               '<tspan dy="12"> = 10</tspan><tspan dy="-12" font-size="19">5.8</tspan><tspan dy="12"> ≈ 6.3 × 10</tspan>'
               '<tspan dy="-12" font-size="19">5</tspan></text>' % INK)
f.text(480, 360, 'the side with the weaker acid (higher pKa) is favoured', 22, color=GREY)
save('equilibrium', f)

# 11. Ionisation curves for pKa 4.8 and 6.0.
f = Figure('Fraction present as A⁻ against pH for pKa 4.8 and 6.0', '', 960, 450)
gx0, gx1, gy0, gy1 = 110, 900, 380, 50   # pH 1..9 → x, fraction 0..1 → y
def gx(ph): return gx0 + (gx1 - gx0)*(ph - 1)/8
def gy(fr): return gy0 + (gy1 - gy0)*fr
f.line(gx0, gy0, gx1, gy0, width=2).line(gx0, gy0, gx0, gy1 - 10, width=2)
for ph in range(1, 10):
    f.line(gx(ph), gy0, gx(ph), gy0 + 7, width=2).text(gx(ph), gy0 + 30, ph, 19, color=GREY)
for fr in (0, 0.5, 1):
    f.line(gx0 - 7, gy(fr), gx0, gy(fr), width=2).text(gx0 - 14, gy(fr) + 7, f'{fr:g}', 19, color=GREY, anchor='end')
    if fr: f.line(gx0, gy(fr), gx1, gy(fr), color='#d5e0e6', width=1, dash='5 5')
f.text(505, 440, 'pH', 21).text(gx0 - 14, gy1 - 18, 'fraction A⁻', 19, anchor='start')
for pka, col, lab in ((4.8, BLUE, 'pKa 4.8'), (6.0, GREEN, 'pKa 6.0')):
    pts = ' L '.join(f'{gx(p/20):.1f} {gy(1/(1 + 10**(pka - p/20))):.1f}' for p in range(20, 181))
    f.path('M ' + pts, color=col, width=4)
    f.circle(gx(pka), gy(0.5), 7, color=col, fill=col, width=1)
f.text(gx(4.8) - 18, gy(0.5) - 8, 'pKa 4.8', 21, color=BLUE, anchor='end', weight='bold')
f.text(gx(6.0) + 18, gy(0.5) + 26, 'pKa 6.0', 21, color=GREEN, anchor='start', weight='bold')
f.line(gx(4.8), gy0, gx(4.8), gy(0.5), color=GREY, width=1.5, dash='4 4')
fr6 = 1/(1 + 10**1.2)
f.circle(gx(4.8), gy(fr6), 6, color=GREEN, fill='white', width=3)
f.text(gx(4.8) + 40, gy0 - 8, '5.9 % at pH 4.8', 20, color=GREEN, anchor='start', weight='bold')
save('curve', f)

# 12. Populations: a ratio of definite species, not half-molecules.
f = Figure('Illustrative populations one unit below, at and one unit above pKa', '', 960, 420)
panels = [('pH = pKa − 1', 1, 10, 'A⁻ : HA = 1 : 10', '1/11 ≈ 9.1 % A⁻'),
          ('pH = pKa', 5, 5, 'A⁻ : HA = 1 : 1', '5/10 = 50 % A⁻'),
          ('pH = pKa + 1', 10, 1, 'A⁻ : HA = 10 : 1', '10/11 ≈ 90.9 % A⁻')]
for i, (title, na, nh, ratio, frac) in enumerate(panels):
    cx = 160 + 320*i
    f.rect(cx - 140, 60, 280, 210, fill='#eef4f7')
    f.text(cx, 40, title, 23, weight='bold')
    species = ['A'] * na + ['H'] * nh
    for k, s in enumerate(species):
        r, c = divmod(k, 4)
        x = cx - 96 + 64*c; yy = 100 + 58*r
        if s == 'A':
            f.circle(x, yy, 24, color=GREEN, fill=PALE_GREEN).text(x, yy + 7, 'A⁻', 19, color=GREEN, weight='bold')
        else:
            f.circle(x, yy, 24, color=BLUE, fill=PALE).text(x, yy + 7, 'HA', 19, color=BLUE, weight='bold')
    f.text(cx, 310, ratio, 22).text(cx, 345, frac, 22, color=RED, weight='bold')
f.text(480, 400, 'each particle is either HA or A⁻; the fraction describes the population', 21, color=GREY)
save('population', f)

# 13. Glycine at pH 7.4: each group follows its own pKa.
f = Figure('Glycine as the neutral form and as the zwitterion at pH 7.4', '', 960, 400)
place(f, 'NCC(=O)O', 220, 150, bond_px=64, highlight_atoms=[0])
f.text(220, 285, 'drawn neutral form', 23, weight='bold').text(220, 322, 'almost absent in water', 21, color=GREY)
f.line(400, 160, 520, 160, color=BLUE, arrow=True).text(460, 140, 'pH 7.4', 21, color=BLUE)
g = place(f, '[NH3+]CC(=O)[O-]', 730, 150, bond_px=64, highlight_atoms=[0, 4])
f.text(730, 285, 'zwitterion', 23, weight='bold')
f.text(730, 318, 'COOH, pKa 2.3: pH far above → COO⁻', 20, color=GREY)
f.text(730, 344, 'NH₃⁺, pKa 9.6: pH below → stays NH₃⁺', 20, color=GREY)
f.text(220, 360, 'compare pH with each pKa', 21, color=RED)
save('glycine', f)

# 14. Extraction: aqueous speciation is only one part of the model.
f = Figure('Aqueous ionisation and partition in a two-phase extraction', '', 960, 420)
f.rect(100, 40, 520, 140, fill='#f3ead6')
f.rect(100, 190, 520, 170, fill=PALE)
f.text(120, 72, 'organic layer', 22, anchor='start', weight='bold')
f.text(120, 222, 'water at a set pH', 22, anchor='start', weight='bold')
f.circle(400, 115, 26, color=BLUE, fill=PALE).text(400, 122, 'HA', 20, color=BLUE, weight='bold')
f.circle(400, 300, 26, color=BLUE, fill=PALE).text(400, 307, 'HA', 20, color=BLUE, weight='bold')
eq_arrow(f, 440, 520, 300)
f.circle(560, 300, 26, color=GREEN, fill=PALE_GREEN).text(560, 307, 'A⁻', 20, color=GREEN, weight='bold')
f.line(385, 268, 385, 148, color=BLUE, arrow=True).line(415, 148, 415, 268, color=BLUE, arrow=True)
f.text(650, 105, 'HA between layers:', 21, color=RED, anchor='start', weight='bold')
f.text(650, 135, 'needs partition coefficient', 21, color=RED, anchor='start')
f.text(650, 163, 'and phase volumes', 21, color=RED, anchor='start')
f.text(650, 285, 'HA ⇌ A⁻ in water:', 21, color=GREEN, anchor='start', weight='bold')
f.text(650, 315, 'pKa and pH (today)', 21, color=GREEN, anchor='start')
f.text(480, 400, 'A⁻ stays in the water; the ionisation curve alone cannot predict recovery', 21, color=GREY)
save('extraction', f)

print('Built 14 acids-bases skeletal figures.')
