"""Curved-arrow and energy-profile figures for the Day 8 mechanisms lecture.

Molecules are RDKit skeletal drawings placed with molfig; every curved arrow is
anchored on the atom coordinates that `place` returns, so tails sit on drawn lone
pairs or bond midpoints and heads land on the accepting atom or new-bond position.
"""
import math
import sys
from pathlib import Path

from rdkit import Chem as _Chem

import molfig
from chemfig import Figure, BLUE, RED, INK, GREEN
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'mechanisms'/'figures'


class _KeepH:
    """molfig parses SMILES with RDKit's default, which deletes explicit [H] atoms.
    Mechanisms need the transferring H drawn as an atom, so keep bracket H here."""
    def __getattr__(self, name):
        return getattr(_Chem, name)

    @staticmethod
    def MolFromSmiles(smiles):
        p = _Chem.SmilesParserParams()
        p.removeHs = False
        return _Chem.MolFromSmiles(smiles, p)


molfig.Chem = _KeepH()

from rdkit.Chem import rdDepictor as _rdDepictor


class _Layout:
    """molfig always computes its own 2D layout. RDKit draws CH3–X with explicit H
    as a cross, which puts an H exactly where an SN2 nucleophile must approach.
    `place_xy` passes hand-set coordinates (RDKit units, bond = 1.5, y up)."""
    pending = None

    def __getattr__(self, name):
        return getattr(_rdDepictor, name)

    def Compute2DCoords(self, mol, **kw):
        _rdDepictor.Compute2DCoords(mol, **kw)
        if _Layout.pending:
            conf = mol.GetConformer()
            for i, (x, y) in _Layout.pending.items():
                conf.SetAtomPosition(i, (x, y, 0))


molfig.rdDepictor = _Layout()


def polar(angles, r=1.5):
    return [(r*math.cos(math.radians(a)), r*math.sin(math.radians(a))) for a in angles]


def place_xy(f, smiles, x, y, coords, **kw):
    _Layout.pending = dict(enumerate(coords))
    try:
        return place(f, smiles, x, y, **kw)
    finally:
        _Layout.pending = None


def methyl_x(f, x_symbol, x, y, side='right', **kw):
    """CH3–X with X on `side` and the three H folded back, leaving the far side open.
    Atom order: 0 = C, 1 = X, 2–4 = H."""
    s = 1 if side == 'right' else -1
    hs = [112, 232, 295] if s == 1 else [68, 308, 245]
    pts = [(0, 0), (1.5*s, 0)] + polar(hs)
    return place_xy(f, f'C({x_symbol})([H])([H])[H]', x, y, pts, **kw)


def save(name, f):
    f.save(OUT/f'sk-{name}.svg')


def lp_at(at, angle, distance=24):
    a = math.radians(angle)
    return at[0] + distance*math.cos(a), at[1] - distance*math.sin(a)


def lp(f, at, *angles, distance=24):
    for a in angles:
        lone_pair(f, at, a, distance=distance)


def mid(a, b):
    return (a[0] + b[0])/2, (a[1] + b[1])/2


def push(f, start, end, bend=0.4, tail=6, head=16, color=RED):
    """Curved arrow with separate trimming at the tail (lone pair/bond) and head (atom)."""
    (x1, y1), (x2, y2) = start, end
    d = math.hypot(x2 - x1, y2 - y1)
    ux, uy = (x2 - x1)/d, (y2 - y1)/d
    curved_arrow(f, (x1 + ux*tail, y1 + uy*tail), (x2 - ux*head, y2 - uy*head), bend=bend, shorten=0, color=color)


def react(f, x, y, length=90, eq=False, label=None):
    if eq:
        f.line(x - length/2, y - 5, x + length/2, y - 5, color=INK, width=2.5, arrow=True)
        f.line(x + length/2, y + 5, x - length/2, y + 5, color=INK, width=2.5, arrow=True)
    else:
        f.line(x - length/2, y, x + length/2, y, color=INK, width=3, arrow=True)
    if label:
        f.text(x, y - 22, label, 19, color=BLUE)


def plus(f, x, y):
    f.text(x, y + 12, '+', 36)


BP = 56          # bond_px; RDKit draws bonds at 1.5 x this (84 px)


def tag(f, at, label, dx, dy, color=GREEN, size=22):
    """Small bold atom label (e.g. oxygen A/B) offset from an atom position."""
    f.text(at[0] + dx, at[1] + dy + size*0.35, label, size, color=color, weight='bold')


def caption(f, y, *lines, size=22, color=BLUE, gap=32):
    for i, line in enumerate(lines):
        f.text(480, y + i*gap, line, size, color=color)


# 1. Proton transfer: two arrows, one lone pair in and one bond pair out.
def proton_transfer():
    f = Figure('Proton transfer from methanol to hydroxide',
               'Hydroxide oxygen B uses a lone pair to take the methanol H; the O–H bond pair stays on methanol oxygen A.', 960, 400)
    y = 170
    hyd = place(f, '[H][O-]', 95, y, flip=True, bond_px=BP)
    lp(f, hyd[1], 100, -100)
    lp(f, hyd[1], -10, distance=30)
    tag(f, hyd[1], 'B', 0, 48)
    meoh = place(f, 'CO[H]', 345, y, flip=True, bond_px=BP)
    lp(f, meoh[1], 50, 130)
    tag(f, meoh[1], 'A', 0, -48)
    plus(f, 200, y + 30)
    push(f, lp_at(hyd[1], -10, 30), meoh[2], bend=0.55)
    push(f, mid(meoh[1], meoh[2]), lp_at(meoh[1], -75, 24), bend=-0.9, tail=2, head=0)
    react(f, 515, y, eq=True)
    w = place(f, '[H]O[H]', 665, y - 15, bond_px=BP)
    lp(f, w[1], -50, -130)
    tag(f, w[1], 'B', 48, 0)
    plus(f, 775, y)
    meo = place(f, 'C[O-]', 865, y, bond_px=BP)
    lp(f, meo[1], 100, -100)
    lp(f, meo[1], -10, distance=30)
    tag(f, meo[1], 'A', 0, 48)
    f.text(110, 290, 'tail on a lone pair', 21, color=RED)
    f.text(345, 290, 'tail on a bond', 21, color=RED)
    caption(f, 350, 'Oxygen A stays with carbon; oxygen B becomes water. Total charge −1 on both sides.')
    save('proton-transfer', f)
    return hyd, meoh


def bromide(f, x, y, angles=(90, 0, -90, 180)):
    at = place(f, '[Br-]', x, y, bond_px=BP)[0]
    lp(f, at, *angles, distance=26)
    return at


def chloride(f, x, y, angles=(90, 0, -90, 180)):
    at = place(f, '[Cl-]', x, y, bond_px=BP)[0]
    lp(f, at, *angles, distance=26)
    return at


def hydroxide(f, x, y, toward=0, flip=True, rotate=0):
    """H–O(−) with three lone pairs; one points along `toward` (degrees) for an arrow tail."""
    a = place(f, '[H][O-]', x, y, flip=flip, rotate=rotate, bond_px=BP)
    o = a[1]
    lp(f, o, toward + 100, toward - 100)
    lp(f, o, toward, distance=30)
    return o, lp_at(o, toward, 30)


# 2. Ammonia takes a proton from HCl.
def ammonia():
    f = Figure('Ammonia takes a proton from hydrogen chloride',
               'The nitrogen lone pair forms a new N–H bond while the H–Cl bond pair leaves with chlorine; NH4+ and Cl− result.', 960, 400)
    y = 175
    n = place(f, '[H]N([H])[H]', 110, y, rotate=-90, bond_px=BP)
    lp(f, n[1], 0, distance=26)
    hcl = place(f, '[H]Cl', 330, y, flip=True, bond_px=BP)
    lp(f, hcl[1], 90, -90, 0, distance=26)
    push(f, lp_at(n[1], 0, 26), hcl[0], bend=0.55)
    push(f, mid(hcl[0], hcl[1]), lp_at(hcl[1], -90, 26), bend=-0.7, tail=2, head=4)
    react(f, 490, y)
    nh4 = place(f, '[H][N+]([H])([H])[H]', 650, y, bond_px=BP)
    plus(f, 775, y)
    chloride(f, 860, y)
    f.text(110, 300, 'N: 3 bonds + 1 pair', 21, color=BLUE)
    f.text(650, 300, 'N: 4 bonds, no pair → +1', 21, color=BLUE)
    caption(f, 360, 'The tail starts on nitrogen’s pair, never on the proton. Net charge 0 → 0.')
    save('ammonia', f)


# 3. Electron-pair sources and acceptors.
def source_sink():
    f = Figure('Electron-pair sources and electron-pair acceptors',
               'Top: hydroxide, ammonia and water with lone pairs shaded as sources. Bottom: acetaldehyde carbonyl carbon, the carbon of bromomethane and the H of HCl shaded as acceptors with partial charges.', 960, 505)
    f.text(30, 40, 'sources: a lone pair (nucleophile or base)', 23, color=BLUE, anchor='start', weight='bold')
    f.text(30, 275, 'acceptors: an electron-poor atom (electrophile or acid)', 23, color=RED, anchor='start', weight='bold')
    y1, y2 = 130, 390
    o = place(f, '[H][O-]', 160, y1, flip=True, bond_px=BP, highlight_atoms=[1])[1]
    lp(f, o, 100, -100); lp(f, o, -10, distance=30)
    n = place(f, '[H]N([H])[H]', 480, y1, bond_px=BP, highlight_atoms=[1])[1]
    lp(f, n, 90, distance=26)
    w = place(f, '[H]O[H]', 800, y1 + 15, rotate=180, bond_px=BP, highlight_atoms=[1])[1]
    lp(f, w, 50, 130)
    c = place(f, '[H]C(C)=O', 160, y2, bond_px=BP, highlight_atoms=[1])
    f.text(c[1][0] + 28, c[1][1] + 8, 'δ+', 21, color=RED, anchor='start')
    f.text(c[3][0] + 22, c[3][1] + 4, 'δ−', 21, color=BLUE, anchor='start')
    m = methyl_x(f, 'Br', 450, y2, bond_px=BP, highlight_atoms=[0])
    f.text(m[0][0] + 6, m[0][1] - 18, 'δ+', 21, color=RED, anchor='start')
    f.text(m[1][0] + 4, m[1][1] - 26, 'δ−', 21, color=BLUE, anchor='start')
    h = place(f, '[H]Cl', 800, y2, flip=True, bond_px=BP, highlight_atoms=[0])
    f.text(h[0][0], h[0][1] - 30, 'δ+', 21, color=RED)
    f.text(h[1][0], h[1][1] - 30, 'δ−', 21, color=BLUE)
    for x, lab in [(160, 'hydroxide'), (480, 'ammonia'), (800, 'water')]:
        f.text(x, 240, lab, 21)
    for x, lab in [(160, 'carbonyl C'), (450, 'C bonded to Br'), (800, 'H bonded to Cl')]:
        f.text(x, 492, lab, 21)
    save('source-sink', f)


# 4. Same hydroxide, two roles.
def base_vs_nucleophile():
    f = Figure('Hydroxide as a nucleophile and as a base',
               'Top: hydroxide lone pair attacks the carbon of bromomethane, so it acts as a nucleophile. Bottom: hydroxide takes a proton from hydronium, so it acts as a base.', 960, 470)
    for row, y in enumerate((125, 345)):
        o, tail = hydroxide(f, 95, y)
        plus(f, 190, y + 30)
        if row == 0:
            m = methyl_x(f, 'Br', 300, y, bond_px=BP)
            lp(f, m[1], 90, 0, -90, distance=26)
            push(f, tail, m[0], bend=0.55, head=12)
            push(f, mid(m[0], m[1]), lp_at(m[1], -90, 26), bend=-0.7, tail=2, head=4)
            prod, role, how = 'CH₃OH + Br⁻', 'nucleophile', 'pair goes to C'
        else:
            h3o = place(f, '[H][O+]([H])[H]', 330, y, rotate=-90, bond_px=BP)
            lp(f, h3o[1], 0, distance=26)
            hl = min((h3o[i] for i in (0, 2, 3)), key=lambda q: q[0])
            push(f, tail, hl, bend=0.55)
            push(f, mid(hl, h3o[1]), lp_at(h3o[1], 240, 24), bend=0.8, tail=2, head=2)
            prod, role, how = 'H₂O + H₂O', 'base', 'pair goes to H'
        react(f, 520, y)
        f.text(600, y + 10, prod, 28, anchor='start')
        f.text(600, y + 55, f'{role}: {how}', 22, color=RED if row == 0 else BLUE, anchor='start', weight='bold')
    f.line(40, 235, 920, 235, color='#c9d6dd', width=2)
    save('base-vs-nucleophile', f)


# 5. Concerted SN2 at a methyl carbon.
def sn2():
    f = Figure('One-step substitution of bromomethane by hydroxide',
               'The hydroxide lone pair forms the C–O bond at the back of the carbon while the C–Br bond pair leaves with bromide, all in one step; the CH3 umbrella flips.', 960, 400)
    y = 165
    o, tail = hydroxide(f, 90, y)
    m = methyl_x(f, 'Br', 330, y, bond_px=BP)
    lp(f, m[1], 90, 0, -90, distance=26)
    push(f, tail, m[0], bend=0.5, head=12)
    push(f, mid(m[0], m[1]), lp_at(m[1], -90, 26), bend=-0.7, tail=2, head=4)
    react(f, 520, y)
    pts = [(0, 0), (-1.5, 0), (-2.4, 0.9)] + polar([80, 325, 260])
    p = place_xy(f, 'C(O[H])([H])([H])[H]', 690, y, pts, bond_px=BP)
    lp(f, p[1], 205, 280)
    plus(f, 805, y + 10)
    bromide(f, 880, y)
    f.text(210, 290, 'two arrows, one step', 22, color=RED, weight='bold')
    f.text(690, 290, 'H’s fold the other way', 21, color=BLUE)
    caption(f, 350, 'OH⁻ + CH₃Br → CH₃OH + Br⁻     net charge −1 → −1')
    save('sn2', f)


# 6. Ammonia as a neutral nucleophile, then a second ammonia as the base.
def amine_two_steps():
    f = Figure('Ammonia substitution, then deprotonation',
               'Step 1: ammonia attacks bromomethane to give methylammonium and bromide. Step 2: a second ammonia removes an N–H proton to give methylamine and ammonium.', 960, 470)
    y = 120
    n = place(f, '[H]N([H])[H]', 85, y, rotate=-90, bond_px=BP)
    lp(f, n[1], 0, distance=26)
    plus(f, 175, y + 30)
    m = methyl_x(f, 'Br', 300, y, bond_px=BP)
    lp(f, m[1], 90, 0, -90, distance=26)
    push(f, lp_at(n[1], 0, 26), m[0], bend=0.5, head=12)
    push(f, mid(m[0], m[1]), lp_at(m[1], -90, 26), bend=-0.7, tail=2, head=4)
    react(f, 490, y, label='step 1')
    pts = [(0, 0), (1.5, 0), (1.5, 1.5), (1.5, -1.5), (3.0, 0)] + polar([100, 215, 280])
    p = place_xy(f, 'C([N+]([H])([H])[H])([H])([H])[H]', 640, y, pts, bond_px=40)
    plus(f, 815, y + 10)
    bromide(f, 885, y)
    f.line(40, 240, 920, 240, color='#c9d6dd', width=2)
    y = 355
    n2 = place(f, '[H]N([H])[H]', 85, y, rotate=-90, bond_px=BP)
    lp(f, n2[1], 0, distance=26)
    plus(f, 175, y + 30)
    pts = [(1.5, 0), (0, 0), (-0.75, 1.3), (-0.75, -1.3), (-1.5, 0)]
    q = place_xy(f, 'C[N+]([H])([H])[H]', 330, y, pts, bond_px=BP)
    push(f, lp_at(n2[1], 0, 26), q[4], bend=0.5)
    push(f, mid(q[4], q[1]), lp_at(q[1], -120, 24), bend=0.8, tail=2, head=2)
    react(f, 490, y, label='step 2')
    f.text(580, y + 10, 'CH₃NH₂ + NH₄⁺', 28, anchor='start')
    f.text(580, y + 55, 'neutral methylamine needs a base', 21, color=BLUE, anchor='start')
    f.text(640, 215, 'N: 4 bonds → +1', 21, color=BLUE)
    save('amine-two-steps', f)


# 7. Heterolysis bookkeeping: methyl versus tert-butyl.
def heterolysis():
    f = Figure('Heterolysis of a C–Br bond: methyl versus tert-butyl',
               'Both C–Br bond pairs move to bromine on paper. The methyl cation is far too unstable to form, while the tert-butyl cation does form in polar solvents, which is SN1 in Week 3.', 960, 470)
    y = 120
    m = methyl_x(f, 'Br', 150, y, bond_px=BP)
    lp(f, m[1], 90, 0, -90, distance=26)
    push(f, mid(m[0], m[1]), lp_at(m[1], -90, 26), bend=-0.7, tail=2, head=4)
    react(f, 370, y)
    place(f, '[H][C+]([H])[H]', 520, y, bond_px=BP)
    plus(f, 600, y + 10)
    bromide(f, 670, y)
    f.text(735, y - 5, 'methyl cation:', 20, color=RED, weight='bold', anchor='start')
    f.text(735, y + 25, 'never forms', 20, color=RED, weight='bold', anchor='start')
    f.line(40, 240, 920, 240, color='#c9d6dd', width=2)
    y = 360
    t = place(f, 'CC(C)(C)Br', 150, y, rotate=-120, bond_px=BP)
    br = t[4]
    a = math.degrees(math.atan2(-(br[1] - t[1][1]), br[0] - t[1][0]))
    lp(f, br, a, a + 90, a - 90, distance=26)
    push(f, mid(t[1], br), lp_at(br, a - 90, 26), bend=-0.7, tail=2, head=4)
    react(f, 370, y)
    place(f, 'C[C+](C)C', 520, y, bond_px=BP)
    plus(f, 600, y + 10)
    bromide(f, 670, y)
    f.text(735, y - 5, 'tert-butyl cation:', 20, color=GREEN, weight='bold', anchor='start')
    f.text(735, y + 25, 'real (SN1, Week 3)', 20, color=GREEN, weight='bold', anchor='start')
    save('heterolysis', f)


def acetaldehyde_addition(f, y, x0=0, nuc='oh'):
    """Hydroxide below-left of acetaldehyde; returns the carbonyl atom coordinates."""
    c = place(f, '[H]C(C)=O', x0 + 330, y - 10, bond_px=BP)
    lp(f, c[3], 35, 145, distance=26)
    o, tail = hydroxide(f, x0 + 110, y + 60, toward=20)
    push(f, tail, c[1], bend=0.35, head=12)
    push(f, mid(c[1], c[3]), lp_at(c[3], -15, 34), bend=-1.0, tail=2, head=2)
    return c


# 8. Hydroxide adds to a carbonyl.
def carbonyl_addition():
    f = Figure('Hydroxide adds to acetaldehyde',
               'The hydroxide lone pair forms a C–O bond to the carbonyl carbon while the C=O pi pair moves onto oxygen, giving a tetrahedral alkoxide with charge −1.', 960, 340)
    y = 140
    c = acetaldehyde_addition(f, y)
    f.text(c[1][0] + 30, c[1][1] + 8, 'δ+', 21, color=RED, anchor='start')
    react(f, 530, y)
    p = place(f, 'CC([H])(O)[O-]', 720, y, bond_px=BP)
    lp(f, p[4], 105, 175, 240, distance=26)
    caption(f, 305, 'C keeps four bonds; the π pair becomes a lone pair on O. Charge −1 → −1.')
    save('carbonyl-addition', f)


def orient(smiles, i, j, angle, **kw):
    """Rotation that makes the bond j→i point along `angle` (degrees, 0 = right, 90 = up)."""
    d = place(Figure('', '', 10, 10), smiles, 0, 0, **kw)
    a0 = math.degrees(math.atan2(-(d[i][1] - d[j][1]), d[i][0] - d[j][0]))
    return angle - a0


def lp_away(f, at, nbr, count=3, spread=None):
    """Lone pairs spread on the side of `at` facing away from its one neighbour."""
    b = math.degrees(math.atan2(-(nbr[1] - at[1]), nbr[0] - at[0]))
    spread = spread or (90 if count == 3 else 60)
    offs = {1: [0], 2: [-spread, spread], 3: [-spread, 0, spread]}[count]
    for o in offs:
        lone_pair(f, at, b + 180 + o, distance=26)
    return [lp_at(at, b + 180 + o, 26) for o in offs]


# 9. Hydride addition to acetone, then protonation by water.
def hydride_two_steps():
    f = Figure('Hydride addition to acetone, then protonation by water',
               'Step 1: a hydride from NaBH4, drawn as H minus, adds to the carbonyl carbon of acetone and the pi pair moves to oxygen. Step 2: the alkoxide takes a proton from water, giving propan-2-ol and hydroxide.', 960, 470)
    y = 115
    c = place(f, 'CC(C)=O', 330, y - 5, bond_px=BP)
    lp(f, c[3], 35, 145, distance=26)
    h = place(f, '[H-]', 120, y + 65, bond_px=BP)[0]
    lone_pair(f, h, 20, distance=26)
    push(f, lp_at(h, 20, 26), c[1], bend=0.35, head=12)
    push(f, mid(c[1], c[3]), lp_at(c[3], -15, 34), bend=-1.0, tail=2, head=2)
    f.text(120, y + 97, 'from NaBH₄', 19, color=BLUE)
    react(f, 520, y, label='step 1')
    p = place(f, 'CC(C)([H])[O-]', 720, y, bond_px=BP)
    lp_away(f, p[4], p[1])
    f.line(40, 235, 920, 235, color='#c9d6dd', width=2)
    y = 355
    r = orient('CC(C)([H])[O-]', 4, 1, 0)
    a = place(f, 'CC(C)([H])[O-]', 120, y, rotate=r, bond_px=BP)
    tails = lp_away(f, a[4], a[1])
    plus(f, 245, y + 35)
    w = place(f, '[H]O[H]', 380, y, rotate=orient('[H]O[H]', 0, 1, 180), bond_px=BP)
    lp(f, w[1], 40, 105, distance=26)
    push(f, tails[1], w[0], bend=0.5)
    push(f, mid(w[0], w[1]), lp_at(w[1], 235, 26), bend=-0.9, tail=2, head=2)
    react(f, 540, y, label='step 2')
    f.text(610, y + 10, '(CH₃)₂CHOH + HO⁻', 28, anchor='start')
    f.text(610, y + 55, 'water supplies the H', 21, color=BLUE, anchor='start')
    save('hydride-two-steps', f)


# 10. Resonance versus reaction.
def resonance():
    f = Figure('Acetate resonance contributors compared with acid dissociation',
               'Top: two acetate contributors linked by a double-headed resonance arrow; only electrons are redrawn. Bottom: acetic acid giving acetate and hydronium, where a proton moves to another molecule.', 960, 470)
    y = 140
    a = place(f, 'CC(=O)[O-]', 230, y, bond_px=66)
    lp(f, a[2], 35, 145, distance=26)
    t3 = lp_away(f, a[3], a[1])
    top = min(t3, key=lambda q: q[1])
    push(f, top, mid(a[1], a[3]), bend=0.9, tail=4, head=2)
    push(f, mid(a[1], a[2]), lp_at(a[2], -15, 34), bend=-1.0, tail=2, head=2)
    f.line(425, y, 535, y, color=INK, width=3, arrow=True, both=True)
    pts = [(-1.3, -0.75), (0, 0), (0, 1.5), (1.3, -0.75)]
    b = place_xy(f, 'CC([O-])=O', 730, y, pts, bond_px=66)
    lp_away(f, b[2], b[1])
    lp_away(f, b[3], b[1], count=2)
    f.text(480, y - 25, 'resonance', 20, color=BLUE)
    caption(f, 270, 'Same nuclei, electrons redrawn: one real ion, both C–O bonds equal.', size=21)
    f.line(40, 305, 920, 305, color='#c9d6dd', width=2)
    f.text(480, 375, 'CH₃COOH + H₂O  ⇌  CH₃COO⁻ + H₃O⁺', 30)
    caption(f, 430, 'A reaction: a proton moves to a different molecule.', size=21, color=RED)
    save('resonance', f)


# ---- energy profiles (chemfig only) ----
def axes(f, x0, y0, x1, ytop, xlabel='reaction progress'):
    f.line(x0, y0, x0, ytop, color=INK, width=2.5, arrow=True)
    f.line(x0, y0, x1, y0, color=INK, width=2.5, arrow=True)
    f.text(x0 - 16, (y0 + ytop)/2, 'free energy G', 20)
    f.parts[-1] = f.parts[-1].replace('<text ', f'<text transform="rotate(-90 {x0 - 16} {(y0 + ytop)/2})" ', 1)
    f.text((x0 + x1)/2, y0 + 30, xlabel, 20)


def profile(f, pts, gy, color=BLUE, width=4, dash=None):
    """pts: [(x, G), ...] alternating minima and maxima; smooth cosine joins. gy maps G to y."""
    d = []
    for (xa, ga), (xb, gb) in zip(pts, pts[1:]):
        for k in range(0 if not d else 1, 41):
            t = k/40
            x = xa + (xb - xa)*t
            g = ga + (gb - ga)*(1 - math.cos(math.pi*t))/2
            d.append(f'{x:.1f},{gy(g):.1f}')
    extra = f' stroke-dasharray="{dash}"' if dash else ''
    f.parts.append(f'<polyline points="{" ".join(d)}" fill="none" stroke="{color}" stroke-width="{width}"{extra}/>')


def span(f, x, ya, yb, label, color, side='right', size=21):
    f.line(x, ya, x, yb, color=color, width=2.5)
    for yy in (ya, yb):
        f.line(x - 7, yy, x + 7, yy, color=color, width=2.5)
    f.text(x + (10 if side == 'right' else -10), (ya + yb)/2 + 7, label, size, color=color,
           anchor='start' if side == 'right' else 'end', weight='bold')


def level(f, x0, x1, y):
    f.line(x0, y, x1, y, color='#9fb3bd', width=1.5, dash='6 5')


# 11. One step: barrier, reverse barrier and overall ΔG.
def profile_one_step():
    f = Figure('One-step reaction profile with forward and reverse barriers',
               'Reactants at 0, transition state at +45 and products at −20 kJ/mol. Forward barrier 45, reverse barrier 65, overall ΔG −20 kJ/mol.', 960, 420)
    gy = lambda g: 250 - 3.6*g
    axes(f, 70, 370, 920, 30)
    pts = [(110, 0), (250, 0), (450, 45), (650, -20), (820, -20)]
    profile(f, pts, gy)
    level(f, 250, 560, gy(0)); level(f, 450, 860, gy(45)); level(f, 650, 860, gy(-20))
    f.text(170, gy(0) - 14, 'reactants 0', 20)
    f.text(740, gy(-20) + 30, 'products −20', 20)
    f.text(450, gy(45) - 16, 'transition state ‡  +45', 20, weight='bold')
    span(f, 540, gy(0), gy(45), 'forward ΔG‡ = 45', RED)
    span(f, 840, gy(-20), gy(45), 'reverse 65', GREEN, side='left')
    span(f, 200, gy(0), gy(-20), 'ΔG = −20', BLUE, side='right')
    level(f, 110, 650, gy(-20))
    f.text(900, 400, 'kJ/mol, schematic', 17, color=BLUE, anchor='end')
    save('profile-one-step', f)


# 12. Two steps: an intermediate sits in a minimum between two transition states.
def profile_two_step():
    f = Figure('Two-step reaction profile with an intermediate',
               'Two maxima (TS1 and TS2) separated by one minimum, the alkoxide intermediate, so two elementary steps.', 960, 420)
    gy = lambda g: 250 - 3.4*g
    axes(f, 70, 370, 920, 30)
    pts = [(110, 0), (200, 0), (340, 45), (480, 15), (620, 35), (760, -15), (880, -15)]
    profile(f, pts, gy)
    f.text(150, gy(0) - 14, 'reactants', 20)
    f.text(340, gy(45) - 14, 'TS1 ‡', 21, color=RED, weight='bold')
    f.text(620, gy(35) - 14, 'TS2 ‡', 21, color=RED, weight='bold')
    f.text(480, gy(15) + 34, 'intermediate', 21, color=GREEN, weight='bold')
    f.text(480, gy(15) + 60, 'CH₃CH(OH)O⁻', 20, color=GREEN)
    f.text(820, gy(-15) - 14, 'products', 20)
    for xa, xb, lab in [(200, 480, 'step 1'), (480, 760, 'step 2')]:
        f.line(xa + 8, 335, xb - 8, 335, color=BLUE, width=2, arrow=True, both=True)
        f.text((xa + xb)/2, 325, lab, 19, color=BLUE)
    f.text(900, 60, 'maximum = TS (not isolable)', 19, anchor='end', color=RED)
    f.text(900, 88, 'minimum = intermediate', 19, anchor='end', color=GREEN)
    save('profile-two-step', f)


# 13. Favourable is not the same as fast.
def favourable_vs_fast():
    f = Figure('Favourable but slow versus fast but unfavourable',
               'Left: products far below reactants but a high barrier, like sugar in air. Right: a low barrier but products above reactants, so the reaction is fast but little product is present at equilibrium.', 960, 420)
    for x0, pts, title, colour, notes in [
        (40, [(80, 0), (150, 0), (240, 75), (340, -60), (420, -60)], 'favourable but slow', RED,
         [('high barrier → slow', 0), ('big −ΔG → K ≫ 1', 1)]),
        (500, [(540, 0), (610, 0), (700, 22), (790, 18), (880, 18)], 'fast but unfavourable', BLUE,
         [('low barrier → fast', 0), ('ΔG > 0 → K < 1', 1)]),
    ]:
        gy = lambda g: 210 - 1.8*g
        f.line(x0 + 20, 340, x0 + 20, 50, color=INK, width=2.5, arrow=True)
        f.line(x0 + 20, 340, x0 + 410, 340, color=INK, width=2.5, arrow=True)
        f.text(x0 + 215, 40, title, 23, color=colour, weight='bold')
        profile(f, pts, gy, color=colour)
        level(f, pts[1][0], pts[-1][0], gy(0))
        for line, k in notes:
            f.text(x0 + 215, 378 + 30*k, line, 20)
    save('favourable-vs-fast', f)


# 14. Catalysed versus uncatalysed path between the same endpoints.
def catalysis_profile():
    f = Figure('Catalysed and uncatalysed paths with the same endpoints',
               'The uncatalysed path has one barrier of 45 kJ/mol. The catalysed path has two smaller humps, the highest at +25. Both start at 0 and end at −20 kJ/mol, so K is unchanged.', 960, 440)
    gy = lambda g: 255 - 3.6*g
    axes(f, 70, 395, 920, 30)
    profile(f, [(110, 0), (230, 0), (450, 45), (670, -20), (860, -20)], gy, color=RED)
    profile(f, [(110, 0), (230, 0), (340, 25), (450, 8), (560, 20), (670, -20), (860, -20)], gy, color=GREEN)
    f.text(900, 45, 'uncatalysed: one barrier, 45 kJ/mol', 21, color=RED, weight='bold', anchor='end')
    f.text(900, 74, 'catalysed: highest barrier 25 kJ/mol', 21, color=GREEN, weight='bold', anchor='end')
    level(f, 110, 900, gy(0)); level(f, 670, 900, gy(-20))
    span(f, 880, gy(0), gy(-20), 'ΔG = −20', BLUE, side='left')
    f.text(480, 368, 'same start, same end → same K; lower barrier → faster both ways', 20, color=BLUE)
    save('catalysis-profile', f)


# 15. Base-catalysed hydration: hydroxide is used, then given back.
def catalytic_hydration():
    f = Figure('Base-catalysed hydration of acetaldehyde',
               'Step 1: hydroxide adds to acetaldehyde to give an alkoxide intermediate. Step 2: the alkoxide takes a proton from water, giving the hydrate and regenerating hydroxide.', 960, 470)
    y = 110
    c = acetaldehyde_addition(f, y)
    react(f, 530, y, label='step 1')
    p = place(f, 'CC([H])(O)[O-]', 720, y, bond_px=BP)
    lp(f, p[4], 105, 175, 240, distance=26)
    f.text(850, y + 10, 'intermediate', 20, color=GREEN, weight='bold')
    f.line(40, 235, 920, 235, color='#c9d6dd', width=2)
    y = 345
    r = orient('CC([H])(O)[O-]', 4, 1, 0)
    a = place(f, 'CC([H])(O)[O-]', 130, y, rotate=r, bond_px=BP)
    tails = lp_away(f, a[4], a[1])
    plus(f, 255, y + 40)
    w = place(f, '[H]O[H]', 390, y, rotate=orient('[H]O[H]', 0, 1, 180), bond_px=BP)
    lp(f, w[1], 40, 105, distance=26)
    push(f, tails[1], w[0], bend=0.5)
    push(f, mid(w[0], w[1]), lp_at(w[1], 235, 26), bend=-0.9, tail=2, head=2)
    react(f, 545, y, label='step 2')
    f.text(610, y + 10, 'CH₃CH(OH)₂ + HO⁻', 28, anchor='start')
    f.text(610, y + 55, 'hydroxide is given back', 21, color=GREEN, anchor='start', weight='bold')
    save('catalytic-hydration', f)


FIGURES = [proton_transfer, ammonia, source_sink, base_vs_nucleophile, sn2, amine_two_steps, heterolysis, carbonyl_addition,
           hydride_two_steps, resonance, profile_one_step, profile_two_step, favourable_vs_fast,
           catalysis_profile, catalytic_hydration]

if __name__ == '__main__':
    import sys
    only = sys.argv[1:]
    for fn in FIGURES:
        if not only or fn.__name__ in only:
            fn()
    print(f'Built {len(FIGURES)} mechanism figures.')


