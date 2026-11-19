"""Skeletal-structure and mechanism figures for the Day 16 carbonyl-addition lecture.

Molecules are drawn with RDKit through molfig; curved arrows and lone pairs are
attached to the atom coordinates that `place` returns.

Local workarounds (molfig itself is shared and left untouched):
- SMILES are parsed with removeHs=False, so an explicit [H] (the proton that
  moves in a proton transfer) is drawn as its own atom with a bond.
- `orient` searches molfig's rotate/flip options so that chosen bonds point in
  chosen directions; reactant and product then share one layout.
- `ring_template` temporarily seeds RDKit's 2D layout with ring-shaped
  coordinates, so an open-chain molecule can be drawn curled up, ready to close.
- Generic groups use CXSMILES atom labels, e.g. '* |$Nu$|'.
"""
import math
from contextlib import contextmanager
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point2D

import molfig
from chemfig import Figure, BLUE, RED, INK, GREEN
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'carbonyl-addition'/'figures'

# --- workaround: keep explicit [H] atoms -----------------------------------
_parse = Chem.MolFromSmiles
def parse(smiles, *args, **kwargs):
    if args or kwargs:
        return _parse(smiles, *args, **kwargs)
    p = Chem.SmilesParserParams()
    p.removeHs = False
    return _parse(smiles, p)
molfig.Chem.MolFromSmiles = parse

# --- workaround: seed 2D coordinates from a template ------------------------
_compute = rdDepictor.Compute2DCoords
@contextmanager
def ring_template(coords):
    """coords: {atom_index: (x, y)} in RDKit units (bond length 1.5)."""
    def seeded(mol, **kw):
        return _compute(mol, coordMap={i: Point2D(*xy) for i, xy in coords.items()})
    molfig.rdDepictor.Compute2DCoords = seeded
    try:
        yield
    finally:
        molfig.rdDepictor.Compute2DCoords = _compute

def polygon(n=6, r=1.5, start=90):
    return [(r*math.cos(math.radians(start - 360*k/n)), r*math.sin(math.radians(start - 360*k/n))) for k in range(n)]

# --- geometry helpers --------------------------------------------------------
def ang(a, b):
    """Direction a -> b in degrees, 0 = right, 90 = up (screen y points down)."""
    return math.degrees(math.atan2(a[1] - b[1], b[0] - a[0]))

def angdiff(a, b):
    return abs((a - b + 180) % 360 - 180)

def orient(smiles, want, **kw):
    """Find rotate/flip so that each bond (i, j) in `want` points at its angle."""
    best = None
    for flip in (False, True):
        for rot in range(0, 360, 5):
            molfig.molecule_svg(smiles, rotate=rot, flip=flip, **kw)
            c = molfig.molecule_svg.last_coords
            err = sum(angdiff(ang(c[i], c[j]), a) for (i, j), a in want.items())
            if best is None or err < best[0]:
                best = (err, rot, flip)
    return dict(rotate=best[1], flip=best[2])

def oplace(f, smiles, x, y, want, **kw):
    return place(f, smiles, x, y, **orient(smiles, want, **kw), **kw)

def neighbours(smiles, i):
    return [n.GetIdx() for n in parse(smiles).GetAtomWithIdx(i).GetNeighbors()]

def free_angle(p, smiles, i):
    """Direction pointing away from the bonds of atom i."""
    ns = neighbours(smiles, i)
    if not ns:
        return 90
    vx = sum(math.cos(math.radians(ang(p[i], p[n]))) for n in ns)
    vy = sum(math.sin(math.radians(ang(p[i], p[n]))) for n in ns)
    if math.hypot(vx, vy) < 1e-6:
        return ang(p[i], p[ns[0]]) + 60
    return math.degrees(math.atan2(-vy, -vx))

def lps(f, p, smiles, i, angles, distance=26, color=BLUE):
    """Lone pairs at angles relative to the free direction of atom i."""
    base = free_angle(p, smiles, i)
    for a in angles:
        lone_pair(f, p[i], base + a, distance=distance, color=color)

def at_angle(at, angle, distance):
    a = math.radians(angle)
    return (at[0] + distance*math.cos(a), at[1] - distance*math.sin(a))

def tag(f, p, smiles, i, label, avoid=(), distance=40, color=GREEN, size=23):
    """Bold atom label (A, B, C) put in the emptiest direction around atom i.

    Bond directions, the side where RDKit writes an attached H, the charge
    corner and any angles in `avoid` (lone pairs, arrows) are kept clear.
    """
    mol = parse(smiles)
    atom = mol.GetAtomWithIdx(i)
    taken = [ang(p[i], p[n]) for n in neighbours(smiles, i)] + list(avoid)
    if atom.GetSymbol() != 'C':
        dx = sum(p[n][0] - p[i][0] for n in neighbours(smiles, i))
        if atom.GetTotalNumHs():
            taken += [180, 165, 195] if dx > 1 else [0, 15, -15]
        if atom.GetFormalCharge():
            taken += [35]
    # Emptiest direction, with a mild preference for the upper half (captions sit below).
    best = max(range(0, 360, 15), key=lambda a: (min(angdiff(a, t) for t in taken) if taken else 0)
               + 25*math.sin(math.radians(a)))
    d = distance + 10*abs(math.cos(math.radians(best)))
    x, y = at_angle(p[i], best, d)
    f.text(round(x), round(y + size*0.35), label, size, color=color, weight='bold')

def bond_arrow(f, a, b, side=1, bend=0.5, reach=26):
    """Curved arrow moving the a-b bond pair onto atom b, bowing out on one side."""
    perp = ang(a, b) + 90*side
    start = at_angle(mid(a, b), perp, 7)
    end = at_angle(at_angle(b, ang(b, a), 4), perp, reach)
    curved_arrow(f, start, end, bend=bend*side, shorten=2)

def mid(a, b, t=0.5):
    return (a[0] + t*(b[0] - a[0]), a[1] + t*(b[1] - a[1]))

def react_arrow(f, x1, x2, y, label='', below=''):
    f.line(x1, y, x2, y, color=BLUE, arrow=True)
    if label: f.text((x1 + x2)/2, y - 16, label, 21, color=BLUE)
    if below: f.text((x1 + x2)/2, y + 32, below, 21, color=BLUE)

def eq_arrow(f, x1, x2, y, label='', below=''):
    f.line(x1, y - 6, x2, y - 6, color=BLUE, arrow=True)
    f.line(x2, y + 6, x1, y + 6, color=BLUE, arrow=True)
    if label: f.text((x1 + x2)/2, y - 20, label, 21, color=BLUE)
    if below: f.text((x1 + x2)/2, y + 38, below, 21, color=BLUE)

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

UP, DL, DR, LEFT, RIGHT, UL, UR = 90, 210, 330, 180, 0, 150, 30
ETHANAL = 'CC=O'                       # 0 CH3, 1 C, 2 O
ETHANAL_UP = {(1, 2): UP, (1, 0): DR}  # O up, CH3 down-right
KW = dict(bond_px=70, font_px=30)


# 1. The carbonyl group: polar and planar ---------------------------------------
f = Figure('The carbonyl group is polar and flat', '', 960, 400)
p = oplace(f, ETHANAL, 250, 160, ETHANAL_UP, bond_px=84, font_px=34)
C, O = p[1], p[2]
f.text(C[0] - 40, C[1] + 8, 'δ+', 30, color=BLUE, weight='bold')
f.text(O[0] + 78, O[1] + 10, 'δ−', 30, color=RED, weight='bold')
lone_pair(f, O, 25, distance=31, color=RED); lone_pair(f, O, 155, distance=31, color=RED)
f.text(250, 300, 'ethanal (acetaldehyde)', 25, weight='bold')
f.text(250, 340, 'O pulls electrons: C is δ+', 22, color=BLUE)
p = oplace(f, 'CC(C)=O', 720, 160, {(1, 3): UP, (1, 0): DL}, bond_px=84, font_px=34)
C = p[1]
f.text(C[0], C[1] + 50, '120°', 23, color=BLUE)
f.text(C[0] + 64, C[1] - 34, '120°', 23, color=BLUE)
f.text(720, 300, 'propanone (acetone)', 25, weight='bold')
f.text(720, 340, 'flat: open to attack from either face', 22, color=BLUE)
save('polar', f)

# 2. The general pattern: Nu adds to C, the pi pair moves onto O ------------------
f = Figure('Nucleophile adds to carbon, pi pair moves to oxygen', '', 960, 400)
ACETONE = 'CC(C)=O'
p = oplace(f, ACETONE, 340, 215, {(1, 3): UP, (1, 0): DL}, **KW)
C, O = p[1], p[3]
nu = (120, 205)
f.text(nu[0], nu[1] + 11, 'Nu', 32, color=BLUE, weight='bold').text(nu[0] + 30, nu[1] - 12, '−', 28, color=BLUE)
lone_pair(f, (nu[0] + 20, nu[1]), 0, distance=26)
curved_arrow(f, (nu[0] + 52, nu[1]), C, bend=0.3, shorten=8)
bond_arrow(f, C, O, side=-1)
react_arrow(f, 480, 580, 205)
PROD = 'CC(C)([O-])* |$;;;;Nu$|'
oplace(f, PROD, 760, 215, {(1, 3): UP, (1, 4): LEFT, (1, 0): RIGHT, (1, 2): 270}, **KW)
f.text(340, 350, 'planar carbonyl', 24, weight='bold').text(760, 350, 'tetrahedral alkoxide', 24, weight='bold')
f.text(480, 50, 'one σ bond in, one π bond out: carbon stays at four bonds', 22, color=BLUE)
save('pattern', f)

# 3. Hydroxide adds to ethanal (oxygens labelled A and B) ------------------------
f = Figure('Hydroxide adds to ethanal', '', 960, 400)
h = place(f, '[OH-]', 110, 205, **KW)
Ob = h[0]
f.text(Ob[0] - 10, Ob[1] - 38, 'B', 23, color=GREEN, weight='bold')
p = oplace(f, ETHANAL, 330, 215, ETHANAL_UP, **KW)
C, Oa = p[1], p[2]
tag(f, p, ETHANAL, 2, 'A', avoid=[0, 20, 340])
lone_pair(f, (Ob[0] + 22, Ob[1]), 0, distance=24)
curved_arrow(f, (Ob[0] + 54, Ob[1]), C, bend=0.3, shorten=8)
bond_arrow(f, C, Oa, side=-1)
react_arrow(f, 470, 570, 205)
INT = 'CC(O)[O-]'                      # 0 CH3, 1 C, 2 O(B)H, 3 O(A)-
INT_W = {(1, 3): UP, (1, 2): LEFT, (1, 0): DR}
p = oplace(f, INT, 760, 215, INT_W, **KW)
tag(f, p, INT, 3, 'A'); tag(f, p, INT, 2, 'B')
f.text(220, 350, 'hydroxide + ethanal', 24, weight='bold')
f.text(740, 350, 'alkoxide, net charge −1', 24, weight='bold')
save('hydroxide-add', f)

# 4. Water protonates the alkoxide: hydrate + regenerated hydroxide ---------------
f = Figure('Water protonates the alkoxide to give the hydrate', '', 960, 400)
p = oplace(f, INT, 160, 220, INT_W, **KW)
Oa = p[3]
tag(f, p, INT, 3, 'A', avoid=[0, 20]); tag(f, p, INT, 2, 'B')
W = '[H]O[H]'
w = oplace(f, W, 345, 150, {(1, 0): LEFT, (1, 2): DR}, bond_px=60, font_px=30)
H, Ow = w[0], w[1]
f.text(Ow[0] + 36, Ow[1] + 2, 'C', 23, color=GREEN, weight='bold')
lone_pair(f, Oa, 0, distance=26)
curved_arrow(f, at_angle(Oa, 0, 30), H, bend=0.4, shorten=10)
bond_arrow(f, H, Ow, side=1)
react_arrow(f, 450, 540, 210)
HYD = 'CC(O)O'                         # 0 CH3, 1 C, 2 O(B)H, 3 O(A)H
p = oplace(f, HYD, 660, 220, {(1, 3): UP, (1, 2): LEFT, (1, 0): DR}, **KW)
tag(f, p, HYD, 3, 'A'); tag(f, p, HYD, 2, 'B')
h = place(f, '[OH-]', 880, 210, **KW)
f.text(h[0][0] - 10, h[0][1] - 38, 'C', 23, color=GREEN, weight='bold')
f.text(220, 345, 'alkoxide + water', 24, weight='bold')
f.text(660, 345, 'hydrate (gem-diol)', 24, weight='bold')
f.text(885, 345, 'OH⁻ back', 24, weight='bold')
f.text(480, 388, 'net: CH₃CHO + H₂O ⇌ CH₃CH(OH)₂, with OH⁻ as catalyst', 21, color=BLUE)
save('hydrate-finish', f)

# 5. How much hydrate? It depends on the carbonyl ---------------------------------
f = Figure('Fraction of carbonyl present as hydrate in water', '', 960, 400)
items = [('C=O', 'methanal', '> 99.9 %', dict(explicit_h=True)), ('CC=O', 'ethanal', '≈ 50 %', {}),
         ('CC(C)=O', 'propanone', '≈ 0.1 %', {}), ('O=CC(Cl)(Cl)Cl', 'chloral', '> 99.9 %', {})]
for i, (smi, name, pct, kw) in enumerate(items):
    x = 115 + 235*i
    place(f, smi, x, 135, bond_px=60 if i == 3 else 64, font_px=30, **kw)
    f.text(x, 265, name, 25, weight='bold').text(x, 310, pct, 28, color=RED, weight='bold')
f.text(480, 365, 'share present as the hydrate at equilibrium in water, 25 °C (approximate)', 21, color=BLUE)
save('hydrate-extent', f)

# 6. Acid activation: the carbonyl oxygen takes a proton from hydronium ----------
f = Figure('Hydronium protonates the carbonyl oxygen', '', 960, 400)
p = oplace(f, ETHANAL, 140, 225, ETHANAL_UP, **KW)
C, O = p[1], p[2]
HY = '[H][OH2+]'
hy = oplace(f, HY, 300, 130, {(0, 1): RIGHT}, bond_px=60, font_px=30)
H, Oh = hy[0], hy[1]
lone_pair(f, O, 30, distance=26)
curved_arrow(f, at_angle(O, 30, 30), H, bend=0.3, shorten=10)
bond_arrow(f, H, Oh, side=1)
react_arrow(f, 430, 510, 215)
PC = 'CC=[OH+]'
oplace(f, PC, 630, 225, ETHANAL_UP, **KW)
f.text(740, 230, '↔', 38)
oplace(f, 'C[CH+]O', 850, 225, {(1, 2): UP, (1, 0): DR}, **KW)
f.text(230, 345, 'ethanal + H₃O⁺', 24, weight='bold')
f.text(740, 345, 'protonated carbonyl (+ H₂O)', 24, weight='bold')
f.text(740, 383, 'right: the resonance form with C⁺', 21, color=BLUE)
save('acid-activation', f)

# 7. Methanol attacks the protonated carbonyl -----------------------------------
f = Figure('Methanol adds to the protonated carbonyl', '', 960, 400)
m = oplace(f, 'CO', 110, 215, {(1, 0): DL}, **KW)
Ob = m[1]
f.text(Ob[0] - 4, Ob[1] - 38, 'B', 23, color=GREEN, weight='bold')
lone_pair(f, Ob, 300, distance=26)
p = oplace(f, PC, 340, 225, ETHANAL_UP, **KW)
C, Oa = p[1], p[2]
f.text(Oa[0] - 40, Oa[1] + 8, 'A', 23, color=GREEN, weight='bold')
curved_arrow(f, at_angle(Ob, 300, 30), C, bend=0.3, shorten=8)
bond_arrow(f, C, Oa, side=-1)
react_arrow(f, 480, 570, 215)
OXH = 'CC(O)[O+]([H])C'                # 0 CH3, 1 C, 2 O(A)H, 3 O(B)+, 4 H, 5 CH3
OXH_W = {(1, 3): LEFT, (1, 2): 60, (1, 0): 300, (3, 4): 120, (3, 5): 240}
p = oplace(f, OXH, 760, 215, OXH_W, **KW)
tag(f, p, OXH, 2, 'A'); tag(f, p, OXH, 3, 'B')
f.text(230, 350, 'methanol + protonated ethanal', 24, weight='bold')
f.text(740, 350, 'B now has three bonds: +1', 24, weight='bold')
save('methanol-attack', f)

# 8. Water removes the proton from oxygen B: the hemiacetal -----------------------
f = Figure('Water deprotonates the oxonium to give the hemiacetal', '', 960, 400)
p = oplace(f, OXH, 230, 220, OXH_W, bond_px=66, font_px=30)
Ob, H = p[3], p[4]
tag(f, p, OXH, 2, 'A'); tag(f, p, OXH, 3, 'B')
w = place(f, 'O', H[0] + 70, H[1] - 70, bond_px=60, font_px=30)
Ow = w[0]
lone_pair(f, Ow, 285, distance=28)
curved_arrow(f, at_angle(Ow, 285, 32), H, bend=0.3, shorten=8)
bond_arrow(f, H, Ob, side=-1, reach=24)
react_arrow(f, 430, 520, 215)
HEMI = 'CC(O)OC'                       # 0 CH3, 1 C, 2 O(A)H, 3 O(B), 4 CH3
HEMI_W = {(1, 3): LEFT, (1, 2): 60, (1, 0): 300, (3, 4): 240}
p = oplace(f, HEMI, 670, 220, HEMI_W, bond_px=66, font_px=30)
tag(f, p, HEMI, 2, 'A'); tag(f, p, HEMI, 3, 'B')
place(f, '[OH3+]', 885, 215, bond_px=60, font_px=30)
f.text(230, 345, 'oxonium + water', 24, weight='bold')
f.text(670, 345, 'hemiacetal', 24, weight='bold')
f.text(885, 345, 'H₃O⁺ back', 24, weight='bold')
f.text(480, 388, 'A (the old C=O oxygen) ends as OH; B (from methanol) ends in OCH₃', 21, color=BLUE)
save('hemiacetal', f)

# 9. Four classes, read at one carbon ------------------------------------------
f = Figure('Aldehyde, hydrate, hemiacetal and acetal', '', 960, 400)
ACETAL_W = {(1, 4): UP, (1, 2): LEFT, (1, 0): DR, (2, 3): DL, (4, 5): UR}
items = [('CC=O', 'aldehyde', 'C=O', [1, 2], ETHANAL_UP),
         ('CC(O)O', 'hydrate', 'OH + OH', [1, 2, 3], {(1, 3): UP, (1, 2): LEFT, (1, 0): DR}),
         ('CC(O)OC', 'hemiacetal', 'OH + OR', [1, 2, 3], HEMI_W),
         ('CC(OC)OC', 'acetal', 'OR + OR', [1, 2, 4], ACETAL_W)]
for i, (smi, name, sub, hl, want) in enumerate(items):
    x = 125 + 237*i
    oplace(f, smi, x, 150, want, bond_px=60, font_px=28, highlight_atoms=hl)
    f.text(x, 285, name, 25, weight='bold').text(x, 325, sub, 23, color=BLUE)
f.text(480, 375, 'read the two groups on the shaded carbon', 21, color=BLUE)
save('classes', f)

# 10. Hemiacetal to acetal: protonate OH, lose water, add a second methanol --------
f = Figure('From hemiacetal to acetal under acid catalysis', '', 960, 430)
PH = 'CC(OC)[OH2+]'                    # 0 CH3, 1 C, 2 OMe O, 3 CH3, 4 OH2+
p = oplace(f, PH, 165, 185, {(1, 4): UP, (1, 2): LEFT, (1, 0): DR, (2, 3): UL}, bond_px=70, font_px=28)
C, Ome, Ow = p[1], p[2], p[4]
lone_pair(f, Ome, 290, distance=26)
curved_arrow(f, at_angle(Ome, 290, 30), at_angle(mid(C, Ome, 0.6), 270, 6), bend=0.5, shorten=4)
bond_arrow(f, C, Ow, side=1)
f.text(160, 315, '① OH protonated', 22, weight='bold').text(160, 347, 'water can now leave', 20, color=BLUE)
react_arrow(f, 290, 370, 190, '−H₂O')
oplace(f, 'CC=[O+]C', 485, 190, {(1, 2): RIGHT, (1, 0): DL, (2, 3): UR}, bond_px=62, font_px=28)
f.text(485, 315, '② oxocarbenium ion', 22, weight='bold').text(485, 347, 'very electrophilic C', 20, color=BLUE)
react_arrow(f, 610, 690, 190, '+CH₃OH', '−H⁺')
oplace(f, 'CC(OC)OC', 820, 190, ACETAL_W, bond_px=62, font_px=28, highlight_atoms=[1, 2, 4])
f.text(820, 315, '③ acetal', 22, weight='bold').text(820, 347, 'two OCH₃ on one C', 20, color=BLUE)
f.text(480, 405, 'net: CH₃CHO + 2 CH₃OH ⇌ CH₃CH(OCH₃)₂ + H₂O', 23)
save('acetal-steps', f)

# 11. Water decides the position of the acetal equilibrium -----------------------
f = Figure('Removing or adding water shifts the acetal equilibrium', '', 960, 400)
oplace(f, ETHANAL, 110, 160, ETHANAL_UP, bond_px=62, font_px=30)
f.text(250, 172, '+ 2 CH₃OH', 27)
eq_arrow(f, 345, 495, 160, 'H⁺ catalyst')
oplace(f, 'CC(OC)OC', 640, 160, ACETAL_W, bond_px=62, font_px=30)
f.text(810, 172, '+ H₂O', 27)
f.rect(110, 270, 350, 92, fill='#e4f0f4')
f.text(285, 306, 'remove water  →', 24, color=GREEN, weight='bold').text(285, 342, 'acetal builds up', 22)
f.rect(500, 270, 350, 92, fill='#f6e6e7')
f.text(675, 306, '←  excess water + acid', 24, color=RED, weight='bold').text(675, 342, 'acetal hydrolyses back', 22)
save('acetal-equilibrium', f)

# 12. Stereocentre or not? ------------------------------------------------------
f = Figure('A new stereocentre forms only if all four groups differ', '', 960, 400)
oplace(f, 'CC(O)OC', 260, 150, HEMI_W, bond_px=74, font_px=30, highlight_atoms=[1])
f.text(260, 290, 'hemiacetal from ethanal', 25, weight='bold').text(260, 330, 'CH₃, H, OH, OCH₃: stereocentre', 22, color=GREEN)
oplace(f, 'CC(C)(O)OC', 700, 150, {(1, 3): UP, (1, 4): LEFT},
       bond_px=74, font_px=30, highlight_atoms=[1])
f.text(700, 290, 'hemiketal from propanone', 25, weight='bold').text(700, 330, 'two identical CH₃: no stereocentre', 22, color=RED)
f.text(480, 380, 'either face of the flat aldehyde can be attacked: a 50:50 (racemic) mixture forms', 20, color=BLUE)
save('stereocentre', f)

# 13. 5-Hydroxypentanal closes to a six-membered cyclic hemiacetal ---------------
f = Figure('5-Hydroxypentanal forms a cyclic hemiacetal', '', 960, 420)
hx = polygon(6, 1.5, 150)              # hexagon vertices, clockwise from upper left
OPEN = 'OCCCCC=O'                      # 0 O(B)H, 1..4 CH2, 5 CHO carbon, 6 O(A)
with ring_template({0: hx[0], 1: hx[5], 2: hx[4], 3: hx[3], 4: hx[2], 5: hx[1]}):
    p = place(f, OPEN, 240, 190, bond_px=66, font_px=30, numbers={5: 'C1', 1: 'C5'})
Ob, C, Oa = p[0], p[5], p[6]
f.text(Ob[0] - 46, Ob[1] + 8, 'B', 23, color=GREEN, weight='bold')
tag(f, p, OPEN, 6, 'A', avoid=[0, 20, 340])
lone_pair(f, Ob, 0, distance=26)
curved_arrow(f, at_angle(Ob, 0, 30), C, bend=-0.35, shorten=10)
bond_arrow(f, C, Oa, side=-1)
react_arrow(f, 420, 530, 190, 'H⁺ transfers')
RING = 'OC1CCCCO1'                     # 0 O(A)H, 1 C1, 2..5 CH2, 6 ring O(B)
with ring_template({6: hx[0], 1: hx[1], 2: hx[2], 3: hx[3], 4: hx[4], 5: hx[5]}):
    p = place(f, RING, 720, 190, bond_px=66, font_px=30, numbers={1: 'C1', 5: 'C5'}, highlight_atoms=[6])
f.text(p[6][0] - 46, p[6][1] + 8, 'B', 23, color=GREEN, weight='bold')
tag(f, p, RING, 0, 'A')
f.text(240, 362, '5-hydroxypentanal, C₅H₁₀O₂', 23, weight='bold')
f.text(720, 362, 'oxan-2-ol, C₅H₁₀O₂', 23, weight='bold')
f.text(480, 402, 'ring = oxygen B + five carbons; same formula, nothing lost', 21, color=BLUE)
save('cyclisation', f)

# 14. Glucose: open chain and beta-pyranose ring ---------------------------------
f = Figure('D-Glucose exists mainly as a cyclic hemiacetal', '', 960, 440)
OPEN_GLC = 'O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO'   # (2R,3S,4R,5R): D-glucose
place(f, OPEN_GLC, 245, 190, bond_px=52, font_px=24, numbers={1: 'C1', 8: 'C5'}, highlight_atoms=[0, 1, 8, 9])
eq_arrow(f, 455, 540, 200)
BETA = 'OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O'   # beta-D-glucopyranose
place(f, BETA, 730, 190, bond_px=58, font_px=24, numbers={4: 'C1', 2: 'C5'}, highlight_atoms=[2, 3, 4, 5])
f.text(245, 360, 'open-chain D-glucose (a trace)', 24, weight='bold').text(245, 398, 'aldehyde at C1, OH on C5', 21, color=BLUE)
f.text(730, 360, 'β-D-glucopyranose', 24, weight='bold').text(730, 398, 'C5 oxygen in the ring, new OH on C1', 21, color=BLUE)
save('glucose', f)

# 15. Other nucleophiles: hydride (Week 3) and cyanide ----------------------------
f = Figure('Hydride and cyanide follow the same pattern', '', 960, 460)
SMALL = dict(bond_px=46, font_px=26)
rows = [('CC(C)[O-]', 'CC(C)O', 'propan-2-ol', 'H⁻', 'from NaBH₄', {(1, 3): UP, (1, 0): DL}),
        ('CC(C)([O-])C#N', 'CC(C)(O)C#N', 'a cyanohydrin', '⁻C≡N', 'from NaCN',
         {(1, 4): RIGHT, (4, 5): RIGHT})]
for row, (inter, prod, name, nu, src, want) in enumerate(rows):
    y = 100 + 210*row
    oplace(f, 'CC(C)=O', 100, y, {(1, 3): UP, (1, 0): DL}, **SMALL)
    react_arrow(f, 185, 285, y + 10, nu, src)
    oplace(f, inter, 425, y, want, **SMALL)
    react_arrow(f, 580, 655, y + 10, 'H⁺')
    oplace(f, prod, 805, y, want, **SMALL)
    f.text(805, y + (82 if row == 0 else 100), name, 22, weight='bold')
f.text(480, 445, 'alkoxide first, then a proton: the same two steps as hydration', 20, color=BLUE)
save('other-nucleophiles', f)

# 16. Aldehyde versus ester tetrahedral intermediate -----------------------------
f = Figure('Aldehyde and ester tetrahedral intermediates', '', 960, 420)
oplace(f, INT, 240, 180, INT_W, **KW)
f.text(240, 315, 'from ethanal + OH⁻', 24, weight='bold')
f.text(240, 350, 'H and CH₃ can\'t leave: it stays added', 20, color=BLUE)
EST = 'CC(O)(OC)[O-]'                  # 0 CH3, 1 C, 2 OH, 3 OMe O, 4 CH3, 5 O-
p = oplace(f, EST, 700, 180, {(1, 5): UP, (1, 2): LEFT, (1, 3): RIGHT, (1, 0): 270, (3, 4): UR}, **KW)
C, Ome, Om = p[1], p[3], p[5]
lone_pair(f, Om, 180, distance=26, color=RED)
curved_arrow(f, at_angle(Om, 180, 30), at_angle(mid(C, Om), 180, 8), bend=-0.7, shorten=3)
bond_arrow(f, C, Ome, side=-1)
f.text(700, 315, 'from methyl ethanoate + OH⁻', 24, weight='bold')
f.text(700, 350, 'CH₃O⁻ can leave: C=O re-forms (Day 17)', 20, color=RED)
f.text(480, 400, 'ask what else is attached to the carbonyl carbon', 21)
save('aldehyde-vs-ester', f)

print('Built 16 carbonyl-addition figures.')
