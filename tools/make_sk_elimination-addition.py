"""Skeletal-structure and mechanism figures for the Day 12 elimination and addition lecture.

Molecules are drawn with RDKit through molfig; curved arrows and lone pairs are
attached to the atom coordinates that `place` returns.

Two workarounds for molfig live here rather than in the shared module:
- `keep_h` keeps hydrogens written as [H] in a SMILES, so a mechanism can show
  the one β-H that reacts without drawing every other H;
- `put` orients a molecule so that a chosen bond points in a chosen direction
  and pins a chosen atom to a chosen point on the canvas.
"""
import contextlib
import math
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor

from chemfig import Figure, BLUE, RED, INK, GREEN
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'elimination-addition'/'figures'
GREY = '#8aa1ac'


def save(name, f):
    f.save(OUT/f'sk-{name}.svg')


@contextlib.contextmanager
def keep_h():
    """molfig parses SMILES with RDKit defaults, which drop explicit [H] atoms."""
    original = Chem.MolFromSmiles
    params = Chem.SmilesParserParams()
    params.removeHs = False
    Chem.MolFromSmiles = lambda s, *a, **k: original(s, params)
    try:
        yield
    finally:
        Chem.MolFromSmiles = original


def _rotation(smiles, a, b, angle, flip, keep):
    """Rotation (degrees, counter-clockwise) that makes bond a→b point at `angle` (0 = right, 90 = up)."""
    params = Chem.SmilesParserParams()
    params.removeHs = not keep
    mol = Chem.MolFromSmiles(smiles, params)
    rdDepictor.Compute2DCoords(mol, canonOrient=True)
    conf = mol.GetConformer()
    pa, pb = conf.GetAtomPosition(a), conf.GetAtomPosition(b)
    dx = (pb.x - pa.x)*(-1 if flip else 1)
    return angle - math.degrees(math.atan2(pb.y - pa.y, dx))


def put(fig, smiles, at=None, xy=(480, 200), bond=None, angle=0, flip=False, keep=False, **kw):
    """Place a molecule. bond=(a, b) with `angle` orients it; at=atom index pins that atom to xy."""
    if bond is not None:
        kw['rotate'] = _rotation(smiles, bond[0], bond[1], angle, flip, keep)
        kw['flip'] = flip
    ctx = keep_h() if keep else contextlib.nullcontext()
    with ctx:
        x, y = xy
        if at is not None:
            probe = Figure('probe')
            pos = place(probe, smiles, 0, 0, **kw)
            x, y = xy[0] - pos[at][0], xy[1] - pos[at][1]
        return place(fig, smiles, x, y, **kw)


def mid(p, q, t=0.5):
    return (p[0] + (q[0] - p[0])*t, p[1] + (q[1] - p[1])*t)


def toward(p, q, dist):
    """Point `dist` px from p in the direction of q (negative: away from q)."""
    d = math.hypot(q[0] - p[0], q[1] - p[1]) or 1
    return (p[0] + (q[0] - p[0])*dist/d, p[1] + (q[1] - p[1])*dist/d)


def deg(p, q):
    """Screen direction from p to q in degrees (0 = right, 90 = up)."""
    return math.degrees(math.atan2(-(q[1] - p[1]), q[0] - p[0]))


def lg_arrow(f, c, x, side=1):
    """Leaving-group arrow: from the middle of bond c–x to just beside atom x."""
    d = math.hypot(x[0] - c[0], x[1] - c[1]) or 1
    ux, uy = (x[0] - c[0])/d, (x[1] - c[1])/d
    end = (x[0] - side*uy*30 + ux*8, x[1] + side*ux*30 + uy*8)
    curved_arrow(f, mid(c, x), end, bend=0.5*side, shorten=3)


def rxn_arrow(f, x1, x2, y, above='', below='', color=BLUE):
    f.line(x1, y, x2, y, arrow=True, color=color)
    if above: f.text((x1 + x2)/2, y - 14, above, 19)
    if below: f.text((x1 + x2)/2, y + 30, below, 19)


def add_h(f, at, angle, length):
    """Draw an explicit C–H bond and H label from atom position `at` (RDKit drops [H] atoms
    from the layout we want, so the one reacting H is drawn by hand)."""
    a = math.radians(angle)
    ux, uy = math.cos(a), -math.sin(a)
    f.line(at[0], at[1], at[0] + ux*length*0.72, at[1] + uy*length*0.72, width=2.5)
    hx, hy = at[0] + ux*length, at[1] + uy*length
    f.text(hx, hy + 10, 'H', 28)
    return (hx, hy)


def with_beta_h(f, smiles, xy, bond_px, leaving=None):
    """Molecule whose atom 0 is Cβ and atom 1 is Cα, Cβ→Cα horizontal, plus a hand-drawn β-H
    pointing up-left. If `leaving` is given, flip so that atom points down (anti to the H)."""
    for flip in (False, True):
        n = len(f.parts)
        s = put(f, smiles, at=0, xy=xy, bond=(0, 1), angle=0, flip=flip, bond_px=bond_px)
        if leaving is None or s[leaving][1] > s[1][1]:
            break
        del f.parts[n:]
    H = add_h(f, s[0], 120, bond_px)
    return H, s


def bromopropane_with_h(f, xy, bond_px):
    """2-bromopropane with its reacting β-H anti to Br. Returns H, Cβ, Cα, Br positions."""
    H, s = with_beta_h(f, 'CC(Br)C', xy, bond_px, leaving=2)
    return H, s[0], s[1], s[2]


# 1. Elimination and addition: the same bonds, opposite directions. -------------
f = Figure('Elimination makes a C=C bond; addition consumes it',
           'Bromoethane loses HBr to give ethene; ethene adds HBr to give bromoethane', 960, 360)
place(f, 'CCBr', 190, 150, bond_px=75, highlight_atoms=[0, 1])
place(f, 'C=C', 690, 150, bond_px=75, highlight_atoms=[0, 1])
f.text(850, 162, '+  HBr', 32)
f.line(370, 125, 560, 125, arrow=True, color=BLUE)
f.line(560, 185, 370, 185, arrow=True, color=GREEN)
f.text(465, 108, 'elimination', 22, color=BLUE, weight='bold')
f.text(465, 222, 'addition', 22, color=GREEN, weight='bold')
f.text(190, 285, 'bromoethane, C₂H₅Br', 23, weight='bold')
f.text(760, 285, 'ethene + HBr', 23, weight='bold')
f.text(480, 335, 'Elimination removes H and Br from neighbouring carbons and makes a π bond.', 19, color=BLUE)
save('overview', f)

# 2. Alpha and beta carbons. ----------------------------------------------------
f = Figure('Alpha and beta carbons in 2-bromobutane',
           '2-Bromobutane with C2 marked alpha and C1 and C3 marked beta; bromomethane has no beta carbon', 960, 400)
p = put(f, 'CC(Br)CC', xy=(300, 150), bond=(0, 1), angle=30, bond_px=85,
        numbers={0: 'C1', 1: 'C2', 3: 'C3', 4: 'C4'}, highlight_atoms=[0, 1, 3])
br_up = p[2][1] < p[1][1]
for i, lab, col in [(1, 'α', RED), (0, 'β', BLUE), (3, 'β', BLUE)]:
    x, y = p[i]
    f.text(x, y + (62 if br_up else -40), lab, 36, color=col, weight='bold')
f.text(300, 335, '2-bromobutane', 23, weight='bold')
f.text(300, 368, 'α carries the Br; β carbons are its neighbours', 20, color=BLUE)
place(f, 'CBr', 760, 165, bond_px=85, explicit_h=True, font_px=26)
f.text(760, 335, 'bromomethane', 23, weight='bold')
f.text(760, 368, 'no β carbon, so no β-elimination', 20, color=BLUE)
f.line(560, 50, 560, 380, color='#c9d8df', width=2)
save('alpha-beta', f)

# 3. E2: three arrows in one step. ------------------------------------------------
f = Figure('E2: three electron pairs move in one step',
           'Ethoxide removes a beta hydrogen from 2-bromopropane while the C–H pair forms the C=C bond and bromide leaves', 960, 330)
H, Cb, Ca, Br = bromopropane_with_h(f, (370, 185), 76)
b = put(f, 'CC[O-]', at=2, xy=(H[0] - 150, H[1] + 10), bond=(1, 2), angle=0, bond_px=60)
O = b[2]
lone_pair(f, O, 15, distance=24)
lone_pair(f, O, -75, distance=24)
curved_arrow(f, toward(O, H, 30), H, bend=0.35, shorten=12)
curved_arrow(f, mid(Cb, H), mid(Cb, Ca), bend=0.7, shorten=2)
lg_arrow(f, Ca, Br, side=-1)
rxn_arrow(f, 590, 670, 200, above='one step')
place(f, 'C=CC', 750, 195, bond_px=58)
f.text(890, 170, '+  EtOH', 24).text(890, 230, '+  Br⁻', 24)
f.text(750, 275, 'propene', 22, weight='bold')
save('e2', f)

# 4. Several beta carbons give several alkenes. ---------------------------------
f = Figure('2-Bromobutane can give three alkenes',
           'Loss of H from C1 gives but-1-ene; loss of H from C3 gives E- or Z-but-2-ene', 960, 400)
put(f, 'CC(Br)CC', xy=(125, 185), bond=(0, 1), angle=30, bond_px=55, numbers={0: '1', 3: '3'}, font_px=24)
rxn_arrow(f, 225, 305, 185, above='EtO⁻', below='− HBr')
items = [('C=CCC', (0, 1), 'but-1-ene', 'H lost from C1', 'minor', [0, 1]),
         ('C/C=C/C', (1, 2), '(E)-but-2-ene', 'H lost from C3', 'major', [1, 2]),
         ('C/C=C\\C', (1, 2), '(Z)-but-2-ene', 'H lost from C3', 'less', [1, 2])]
for i, (smi, bnd, name, src, amt, hl) in enumerate(items):
    x = 430 + 205*i
    put(f, smi, xy=(x, 175), bond=bnd, angle=0 if i else 30, bond_px=50, highlight_atoms=hl, font_px=24)
    f.text(x, 270, name, 22, weight='bold')
    f.text(x, 302, src, 19, color=BLUE)
    f.text(x, 332, amt, 19, color=RED)
f.text(480, 382, 'A small base favours the more substituted alkene; a bulky base (t-BuO⁻) gives more but-1-ene.', 18)
save('butenes', f)

# 5. The E/Z test. -------------------------------------------------------------
f = Figure('E/Z needs two different groups on each alkene carbon',
           'But-1-ene has two hydrogens on one alkene carbon, so no E/Z; but-2-ene has H and CH3 on each end', 960, 400)
items = [('[H]C([H])=C([H])CC', (1, 3), 'but-1-ene', 'C1 carries H and H', 'no E/Z possible', [1, 3]),
         ('C/C([H])=C([H])/C', (1, 3), '(E)-but-2-ene', 'each end: H and CH₃', 'CH₃ groups opposite', [1, 3]),
         ('C/C([H])=C([H])\\C', (1, 3), '(Z)-but-2-ene', 'each end: H and CH₃', 'CH₃ groups same side', [1, 3])]
for i, (smi, bnd, name, sub, note, hl) in enumerate(items):
    x = 165 + 315*i
    put(f, smi, keep=True, xy=(x, 140), bond=bnd, angle=0, bond_px=54, font_px=26, highlight_atoms=hl)
    f.text(x, 292, name, 23, weight='bold')
    f.text(x, 324, sub, 20, color=BLUE)
    f.text(x, 356, note, 20, color=RED if i == 0 else INK)
save('ez-test', f)

# 6. Substitution or elimination: where does the pair go? ------------------------
f = Figure('Attack at carbon or at a beta hydrogen',
           'Hydroxide attacking the alpha carbon of 2-bromopropane from the back gives substitution; removing the beta H gives elimination', 960, 420)
for x0, kind in [(240, 'sub'), (720, 'elim')]:
    if kind == 'sub':
        s = put(f, 'CC(Br)C', at=1, xy=(x0 + 20, 175), bond=(1, 2), angle=-60, bond_px=70)
        Ca, Br = s[1], s[2]
        tgt = toward(Ca, Br, -120)          # back side of the C–Br bond
        o = put(f, '[OH-]', at=0, xy=tgt, bond_px=60)
        O = o[0]
        lone_pair(f, O, deg(O, Ca), distance=26)
        curved_arrow(f, toward(O, Ca, 30), Ca, bend=0.25, shorten=10)
        lg_arrow(f, Ca, Br, side=-1)
        f.text(x0, 350, 'pair goes to the α carbon', 22, weight='bold')
        f.text(x0, 383, 'substitution (SN2) → propan-2-ol', 20, color=BLUE)
    else:
        H, Cb, Ca, Br = bromopropane_with_h(f, (x0 - 20, 190), 70)
        o = put(f, '[OH-]', at=0, xy=(H[0] - 130, H[1] + 5), bond_px=60)
        O = o[0]
        lone_pair(f, O, 15, distance=26)
        curved_arrow(f, toward(O, H, 30), H, bend=0.35, shorten=12)
        curved_arrow(f, mid(Cb, H), mid(Cb, Ca), bend=0.7, shorten=2)
        lg_arrow(f, Ca, Br, side=-1)
        f.text(x0, 350, 'pair goes to a β hydrogen', 22, weight='bold')
        f.text(x0, 383, 'elimination (E2) → propene', 20, color=BLUE)
f.line(480, 30, 480, 395, color='#c9d8df', width=2)
save('sub-vs-elim', f)

# 7. E1: ionise first, then lose a beta H. ----------------------------------------
f = Figure('E1: ionisation, then loss of a beta hydrogen',
           'tert-Butyl bromide ionises to a tertiary carbocation and bromide; water then removes a beta H to give 2-methylpropene and hydronium', 960, 455)
f.text(30, 35, 'Step 1 (slow): C–Br breaks', 21, color=RED, anchor='start', weight='bold')
s = put(f, 'CC(C)(C)Br', at=1, xy=(170, 140), bond=(1, 4), angle=0, bond_px=58)
lg_arrow(f, s[1], s[4], side=-1)
rxn_arrow(f, 330, 410, 145)
put(f, 'C[C+](C)C', at=1, xy=(500, 145), bond=(1, 0), angle=180, bond_px=58)
f.text(640, 155, '+  Br⁻', 26)
f.text(30, 250, 'Step 2 (fast): water takes a β-H', 21, color=RED, anchor='start', weight='bold')
H, c = with_beta_h(f, 'C[C+](C)C', (280, 368), 62)
Cb, Cp = c[0], c[1]
for flip in (False, True):
    n = len(f.parts)
    w = put(f, '[H]O[H]', keep=True, at=1, xy=(H[0] - 140, H[1] + 15), bond=(1, 0), angle=120, flip=flip,
            bond_px=40, font_px=24)
    if w[2][1] > w[1][1] and w[2][0] < w[1][0]:
        break
    del f.parts[n:]
O = w[1]
lone_pair(f, O, 0, distance=24)
curved_arrow(f, toward(O, H, 30), H, bend=0.35, shorten=12)
curved_arrow(f, mid(Cb, H), mid(Cb, Cp), bend=0.7, shorten=2)
rxn_arrow(f, 480, 560, 368)
put(f, 'C=C(C)C', at=1, xy=(660, 368), bond=(1, 0), angle=180, bond_px=58)
f.text(830, 378, '+  H₃O⁺', 26)
f.text(830, 428, '2-methylpropene', 20, weight='bold')
save('e1', f)

# 8. One carbocation, two fates. ----------------------------------------------
f = Figure('A tertiary carbocation can be captured or lose a beta H',
           'The tert-butyl cation in water gives tert-butanol by capture (SN1) or 2-methylpropene by beta-H loss (E1)', 960, 420)
put(f, 'C[C+](C)C', at=1, xy=(480, 100), bond=(1, 0), angle=-90, bond_px=55)
f.text(480, 30, 'tert-butyl cation in water', 21, weight='bold')
f.line(410, 170, 300, 215, arrow=True, color=BLUE)
f.line(550, 170, 660, 215, arrow=True, color=BLUE)
f.text(330, 168, 'H₂O bonds to C⁺', 19, color=BLUE, anchor='end')
f.text(630, 168, 'H₂O removes a β-H', 19, color=BLUE, anchor='start')
put(f, 'CC(C)(C)O', at=1, xy=(230, 290), bond=(1, 4), angle=90, bond_px=44, highlight_atoms=[4])
f.text(230, 385, 'tert-butanol + H₃O⁺', 21, weight='bold')
f.text(230, 406, 'SN1', 19, color=RED, weight='bold')
put(f, 'C=C(C)C', at=1, xy=(730, 300), bond=(1, 0), angle=90, bond_px=44, highlight_atoms=[0, 1])
f.text(730, 385, '2-methylpropene + H₃O⁺', 21, weight='bold')
f.text(730, 406, 'E1', 19, color=RED, weight='bold')
save('e1-sn1', f)

# 9. Energy profiles: E2 one hump, E1 two humps. --------------------------------
f = Figure('Energy profiles of E2 and E1',
           'E2 has a single transition state; E1 has two transition states with a carbocation intermediate between them', 960, 400)
for x0, label, two in [(70, 'E2: one step', False), (530, 'E1: two steps', True)]:
    f.line(x0, 330, x0, 70, width=2, arrow=True, color=INK)
    f.line(x0, 330, x0 + 380, 330, width=2, arrow=True, color=INK)
    f.text(x0 + 190, 362, 'reaction progress', 18)
    f.text(x0 + 6, 62, 'energy', 17, anchor='start')
    f.text(x0 + 190, 34, label, 23, weight='bold')
    if not two:
        d = f'M {x0+20} 260 C {x0+120} 260 {x0+150} 110 {x0+195} 110 C {x0+240} 110 {x0+270} 290 {x0+370} 290'
        f.path(d, color=RED, width=4)
        f.text(x0 + 195, 96, 'one transition state', 17, color=RED)
        f.text(x0 + 20, 285, 'RBr + B⁻', 17, anchor='start').text(x0 + 370, 315, 'alkene', 17, anchor='end')
    else:
        d = (f'M {x0+20} 260 C {x0+80} 260 {x0+100} 100 {x0+140} 100 C {x0+175} 100 {x0+185} 200 {x0+215} 200 '
             f'C {x0+240} 200 {x0+250} 160 {x0+275} 160 C {x0+300} 160 {x0+320} 290 {x0+370} 290')
        f.path(d, color=RED, width=4)
        f.text(x0 + 140, 86, 'slow: C–Br breaks', 17, color=RED)
        f.text(x0 + 215, 228, 'carbocation', 17, color=BLUE, weight='bold')
        f.text(x0 + 20, 285, 'RBr', 17, anchor='start').text(x0 + 370, 315, 'alkene', 17, anchor='end')
save('energy', f)

# 10. Electrophilic addition of HBr to propene, both steps. ----------------------
f = Figure('HBr adds to propene in two steps',
           'The propene pi bond takes H from HBr to give the secondary carbocation and bromide; bromide then bonds to the cation to give 2-bromopropane', 960, 450)
f.text(30, 42, 'Step 1 (slow): the π bond takes H⁺', 21, color=RED, anchor='start', weight='bold')
pr = put(f, 'C=CC', at=0, xy=(110, 195), bond=(0, 1), angle=30, bond_px=72, numbers={0: '1', 1: '2'})
hb = put(f, '[H]Br', keep=True, at=0, xy=(pr[1][0] + 60, pr[1][1] - 60), bond=(0, 1), angle=0, bond_px=70)
H, Br = hb[0], hb[1]
curved_arrow(f, mid(pr[0], pr[1]), H, bend=0.35, shorten=8)
lg_arrow(f, H, Br, side=-1)
rxn_arrow(f, 420, 500, 150)
c = put(f, 'C[CH+]C', at=1, xy=(640, 135), bond=(1, 0), angle=-150, bond_px=72)
f.text(820, 160, '+  Br⁻', 26)
f.text(640, 215, 'secondary carbocation', 20, weight='bold')
f.line(30, 250, 930, 250, color='#c9d8df', width=2)
f.text(30, 290, 'Step 2 (fast): Br⁻ bonds to C⁺', 21, color=RED, anchor='start', weight='bold')
c2 = put(f, 'C[CH+]C', at=1, xy=(200, 360), bond=(1, 0), angle=-150, bond_px=72)
b2 = put(f, '[Br-]', at=0, xy=(c2[1][0] + 150, c2[1][1] - 20), bond_px=60)
lone_pair(f, b2[0], 180, distance=30)
curved_arrow(f, toward(b2[0], c2[1], 34), c2[1], bend=-0.4, shorten=20)
rxn_arrow(f, 420, 500, 375)
put(f, 'CC(C)Br', at=1, xy=(640, 385), bond=(1, 3), angle=90, bond_px=56)
f.text(820, 385, '2-bromopropane', 21, weight='bold')
save('hbr', f)

# 11. Markovnikov: which carbocation? -------------------------------------------
f = Figure('H adds to the end that gives the more stable carbocation',
           'Protonating C1 of propene gives a secondary cation and 2-bromopropane; protonating C2 would give a primary cation', 960, 420)
put(f, 'C=CC', at=1, xy=(480, 75), bond=(0, 1), angle=30, bond_px=60, numbers={0: '1', 1: '2'})
f.line(420, 120, 290, 170, arrow=True, color=BLUE)
f.line(540, 120, 670, 170, arrow=True, color=GREY, dash='8 6')
f.text(350, 128, 'H⁺ to C1', 20, color=BLUE, anchor='end')
f.text(610, 128, 'H⁺ to C2', 20, color=GREY, anchor='start')
put(f, 'C[CH+]C', at=1, xy=(220, 215), bond=(1, 0), angle=-150, bond_px=62, highlight_atoms=[1])
f.text(220, 300, 'secondary: two C neighbours', 21, weight='bold')
f.text(220, 332, 'formed → 2-bromopropane', 21, color=BLUE)
put(f, 'CC[CH2+]', at=1, xy=(740, 205), bond=(1, 0), angle=-150, bond_px=62, highlight_atoms=[2])
f.text(740, 300, 'primary: one C neighbour', 21, weight='bold')
f.text(740, 332, 'too unstable → no 1-bromopropane', 21, color=RED)
f.text(480, 395, 'Alkyl groups stabilise C⁺: tertiary > secondary > primary.', 20)
save('markovnikov', f)

# 12. A planar cation is attacked from both faces. ---------------------------------
f = Figure('A planar carbocation gives both enantiomers',
           'But-1-ene plus HBr goes through a planar secondary cation; bromide attacks either face to give R- and S-2-bromobutane 50:50', 960, 400)
put(f, 'C=CCC', xy=(115, 175), bond=(0, 1), angle=30, bond_px=46)
f.text(110, 265, 'but-1-ene', 21, weight='bold')
rxn_arrow(f, 205, 275, 175, above='H⁺')
put(f, 'C[CH+]CC', at=1, xy=(390, 160), bond=(1, 0), angle=-150, bond_px=60, highlight_atoms=[1])
f.text(390, 265, 'planar C⁺', 21, weight='bold')
f.text(390, 297, 'Br⁻ attacks either face', 19, color=BLUE)
f.line(505, 150, 590, 100, arrow=True, color=BLUE)
f.line(505, 200, 590, 250, arrow=True, color=BLUE)
put(f, 'C[C@@H](Br)CC', at=1, xy=(740, 105), bond=(1, 0), angle=-150, bond_px=52)
put(f, 'C[C@H](Br)CC', at=1, xy=(740, 285), bond=(1, 0), angle=-150, bond_px=52)
f.text(890, 100, '(R)', 24, color=RED, weight='bold')
f.text(890, 280, '(S)', 24, color=RED, weight='bold')
f.text(480, 380, '50 : 50, a racemic mixture of 2-bromobutane', 20)
save('faces', f)

# 13. Br2 adds through a bromonium ion. ------------------------------------------
f = Figure('Br2 adds through a bridged bromonium ion',
           'Cyclohexene pi electrons attack Br2 while the far Br leaves as bromide, giving a bridged bromonium ion', 960, 380)
cy = put(f, 'C1=CCCCC1', xy=(150, 190), bond=(1, 0), angle=0, bond_px=62)
top = mid(cy[0], cy[1])
bb = put(f, 'BrBr', at=0, xy=(top[0] + 110, top[1] - 60), bond=(0, 1), angle=0, bond_px=75)
B1, B2 = bb[0], bb[1]
curved_arrow(f, top, B1, bend=0.4, shorten=10)
lg_arrow(f, B1, B2, side=-1)
rxn_arrow(f, 410, 490, 190)
for flip in (False, True):
    n = len(f.parts)
    bro = put(f, 'C1CCC2[Br+]C2C1', xy=(650, 175), bond=(3, 5), angle=0, flip=flip, bond_px=62)
    if bro[4][1] < bro[3][1]:
        break
    del f.parts[n:]
f.text(650, 355, 'bromonium ion', 21, weight='bold')
f.text(830, 200, '+  Br⁻', 26)
f.text(260, 340, 'π pair to one Br; Br–Br pair leaves', 19, color=RED)
save('bromonium', f)

# 14. Anti addition: back-side opening, trans product. ---------------------------
f = Figure('Bromide opens the bromonium ion from the back face',
           'Bromide attacks one carbon of the bromonium ion from the face opposite the bridge, giving trans-1,2-dibromocyclohexane as a racemic pair and no cis isomer', 960, 420)
# Bromonium drawn as ring + wedged bridge (the bridge sits above the ring plane).
ring = put(f, 'C1CCCCC1', xy=(220, 215), bond=(1, 0), angle=0, bond_px=62)
c1, c2 = ring[1], ring[0]
brx, bry = (c1[0] + c2[0])/2, c1[1] - 70
f.wedge(c1[0], c1[1], brx - 8, bry + 14, spread=7)
f.wedge(c2[0], c2[1], brx + 8, bry + 14, spread=7)
f.text(brx, bry + 6, 'Br⁺', 26, color='#7a3b2e')
f.text(c1[0] - 95, c1[1] + 95, 'Br⁻', 26, color='#7a3b2e')
f.text(c1[0] - 95, c1[1] + 125, 'from below', 17, color=BLUE)
curved_arrow(f, (c1[0] - 80, c1[1] + 70), c1, bend=-0.3, shorten=8)
curved_arrow(f, mid(c1, (brx - 8, bry + 14)), (brx - 34, bry - 4), bend=-0.6, shorten=2)
f.text(200, 385, 'bridge on top (wedges)', 19, color=BLUE)
rxn_arrow(f, 335, 400, 215)
put(f, 'Br[C@@H]1CCCC[C@H]1Br', xy=(530, 205), bond_px=60)
f.text(530, 340, 'trans: anti addition', 21, weight='bold')
f.text(530, 372, 'formed, as a racemic pair', 19, color=BLUE)
f.line(680, 50, 680, 390, color='#c9d8df', width=2)
put(f, 'Br[C@@H]1CCCC[C@@H]1Br', xy=(820, 205), bond_px=60)
f.text(820, 340, 'cis: syn addition', 21, weight='bold')
f.text(820, 372, 'not formed', 19, color=RED)
save('anti', f)

# 15. Hydration and dehydration: conditions choose the direction. ------------------
f = Figure('Cyclohexanol and cyclohexene interconvert',
           'Concentrated acid and heat with distillation dehydrate cyclohexanol to cyclohexene; dilute aqueous acid hydrates cyclohexene back', 960, 380)
put(f, 'OC1CCCCC1', xy=(165, 175), bond=(1, 0), angle=0, bond_px=60, highlight_atoms=[0, 1])
f.text(165, 300, 'cyclohexanol', 22, weight='bold')
f.line(340, 150, 590, 150, arrow=True, color=BLUE)
f.line(590, 205, 340, 205, arrow=True, color=GREEN)
f.text(465, 133, 'elimination (dehydration)', 19, color=BLUE, weight='bold')
f.text(465, 178, 'conc. H₃PO₄, heat, distil', 18, color=BLUE)
f.text(465, 238, 'addition (hydration)', 19, color=GREEN, weight='bold')
f.text(465, 262, 'dilute H₃O⁺, excess water', 18, color=GREEN)
put(f, 'C1=CCCCC1', xy=(720, 175), bond=(1, 0), angle=0, bond_px=60, highlight_atoms=[0, 1])
f.text(870, 185, '+  H₂O', 24)
f.text(740, 300, 'cyclohexene (bp 83 °C)', 22, weight='bold')
f.text(480, 352, 'Same bonds either way. Removing a product (Le Chatelier) decides which way it runs.', 19)
save('hydration', f)

# 16. Aconitase: an elimination and an addition in the citric acid cycle. ----------
f = Figure('Aconitase eliminates water and adds it back',
           'Citrate loses water to give cis-aconitate, then water adds the other way round to give isocitrate; drawn as the free acids', 960, 380)
items = [('OC(=O)CC(O)(CC(=O)O)C(=O)O', 'citrate', 'OH on the central C', [3, 4, 5]),
         ('OC(=O)/C=C(/CC(=O)O)C(=O)O', 'cis-aconitate', 'C=C, water removed', [3, 4]),
         ('OC(=O)C(O)C(CC(=O)O)C(=O)O', 'isocitrate', 'OH now one C over', [3, 4, 5])]
for i, (smi, name, sub, hl) in enumerate(items):
    x = 160 + 320*i
    place(f, smi, x, 170, bond_px=42, font_px=20, highlight_atoms=hl)
    f.text(x, 330, name, 23, weight='bold')
    f.text(x, 362, sub, 19, color=BLUE)
for x, lab, col in [(300, '− H₂O', BLUE), (620, '+ H₂O', GREEN)]:
    f.line(x - 5, 325, x + 35, 325, arrow=True, color=col)
    f.text(x + 15, 305, lab, 18, color=col)
save('aconitase', f)

print('Built 16 elimination-addition figures.')
