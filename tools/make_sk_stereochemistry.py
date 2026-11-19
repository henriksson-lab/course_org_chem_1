"""Skeletal and projection figures for the Day 7 stereochemistry lecture.

Molecules are drawn with RDKit through molfig; wedges and hashes come from the
@/@@ tags in each SMILES. Every R/S and E/Z label written on a figure is
checked against RDKit's CIP labeller before the figure is saved, and every
Fischer projection is converted to a 3D model and labelled the same way.
"""
import math
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdCIPLabeler

import molfig
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'stereochemistry'/'figures'
GREY = '#5b6f7a'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

# ---------------------------------------------------------------- verification

def _parse(smiles):
    p = Chem.SmilesParserParams()
    p.removeHs = False
    return Chem.MolFromSmiles(smiles, p)

def cip(smiles):
    """{atom_index: 'R'/'S'} and {bond_index: 'E'/'Z'} from RDKit's new CIP labeller."""
    m = _parse(smiles)
    rdCIPLabeler.AssignCIPLabels(m)
    atoms = {a.GetIdx(): a.GetProp('_CIPCode') for a in m.GetAtoms() if a.HasProp('_CIPCode')}
    bonds = {b.GetIdx(): b.GetProp('_CIPCode') for b in m.GetBonds() if b.HasProp('_CIPCode')}
    # Cross-check atoms with FindMolChiralCenters (non-legacy implementation).
    m2 = Chem.MolFromSmiles(smiles)
    found = dict(Chem.FindMolChiralCenters(m2, includeUnassigned=True, useLegacyImplementation=False))
    for i, lab in atoms.items():
        assert found.get(i) == lab, (smiles, i, lab, found)
    return atoms, bonds

def expect(smiles, atoms=None, bonds=None):
    a, b = cip(smiles)
    if atoms is not None: assert a == atoms, (smiles, a, atoms)
    if bonds is not None: assert sorted(b.values()) == sorted(bonds), (smiles, b, bonds)
    return smiles

def sense(coords, i1, i2, i3):
    """Rotation 1->2->3 as seen on the page (screen y points down)."""
    (x1, y1), (x2, y2), (x3, y3) = coords[i1], coords[i2], coords[i3]
    cross = (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1)
    return 'clockwise' if cross > 0 else 'counterclockwise'

# Calibrate SMILES @/@@ against 3D geometry with an embedded model, then use it
# to turn a Fischer projection into a SMILES whose CIP label RDKit computes.
def _tag_sign():
    m = Chem.AddHs(Chem.MolFromSmiles('[C@](F)(Cl)(Br)I'))
    AllChem.EmbedMolecule(m, randomSeed=7)
    p = m.GetConformer().GetPositions()
    v = [p[i] - p[0] for i in (1, 2, 3, 4)]
    return np.sign(np.linalg.det(np.array([v[1] - v[0], v[2] - v[0], v[3] - v[0]])))
AT_SIGN = _tag_sign()

FISCHER_VECTORS = {'top': (0, 1, -1), 'right': (1, 0, 1), 'bottom': (0, -1, -1), 'left': (-1, 0, 1)}

def fischer_label(top, right, bottom, left):
    """CIP label of a one-centre Fischer projection; groups are SMILES fragments.

    Fischer convention: horizontal bonds toward the viewer (+z), vertical bonds away.
    """
    groups = dict(top=top, right=right, bottom=bottom, left=left)
    order = ['top', 'right', 'bottom', 'left']
    v = [np.array(FISCHER_VECTORS[k], float) for k in order]
    s = np.sign(np.linalg.det(np.array([v[1] - v[0], v[2] - v[0], v[3] - v[0]])))
    tag = '@' if s == AT_SIGN else '@@'
    smi = f'[C{tag}]({groups["top"]})({groups["right"]})({groups["bottom"]}){groups["left"]}'
    a, _ = cip(smi)
    return a[0]

FRAG = {'COOH': 'C(=O)O', 'CH₃': 'C', 'OH': 'O', 'H': '[H]'}

def lactic_fischer(top, right, bottom, left):
    return fischer_label(*(FRAG[g] for g in (top, right, bottom, left)))

# ------------------------------------------------------------------- drawing

@contextmanager
def keep_explicit_h():
    """molfig parses SMILES with default settings, which deletes [H] atoms.
    Swap in a parser that keeps them, so the H on a stereocentre is drawn with
    its own wedge or hash (a workaround; molfig itself is shared and unchanged)."""
    class Shim:
        def __getattr__(self, name): return getattr(Chem, name)
        @staticmethod
        def MolFromSmiles(s): return _parse(s)
    old = molfig.Chem
    molfig.Chem = Shim()
    try:
        yield
    finally:
        molfig.Chem = old

def place_h(fig, smiles, x, y, **kw):
    with keep_explicit_h():
        return place(fig, smiles, x, y, **kw)

def badge(f, xy, n, color=BLUE, r=16):
    f.circle(xy[0], xy[1], r, color=color, fill='white', width=2.5)
    f.text(xy[0], xy[1] + 7, n, 20, color=color, weight='bold')

def outward(coords, centre, atom, dist):
    cx, cy = coords[centre]; ax, ay = coords[atom]
    d = math.hypot(ax - cx, ay - cy)
    return ax + (ax - cx)/d*dist, ay + (ay - cy)/d*dist

def priorities(f, coords, centre, ranks, dist=42, extra=None):
    """Circled CIP priority next to each substituent atom. extra: {atom: px} for wide labels."""
    for atom, n in ranks.items():
        badge(f, outward(coords, centre, atom, dist + (extra or {}).get(atom, 0)), n)

def arc_arrow(f, centre, angles, r, clockwise, color=RED):
    """Arc around `centre` through the given screen angles (degrees, y down), arrow at the end."""
    cx, cy = centre
    mid = 'arc' + color.lstrip('#')
    if not any(f'id="{mid}"' in p for p in f.parts):
        f.parts.insert(1, f'<defs><marker id="{mid}" viewBox="0 0 10 10" refX="7" refY="5" markerWidth="5" '
                          f'markerHeight="5" orient="auto"><path d="M 0 0 L 10 5 L 0 10 z" fill="{color}"/></marker></defs>')
    pt = lambda a: (cx + r*math.cos(math.radians(a)), cy + r*math.sin(math.radians(a)))
    x0, y0 = pt(angles[0])
    d = f'M {x0:.1f} {y0:.1f}'
    for a1, a2 in zip(angles, angles[1:]):
        sweep = (a2 - a1) % 360 if clockwise else (a1 - a2) % 360
        x, y = pt(a2)
        d += f' A {r} {r} 0 {1 if sweep > 180 else 0} {1 if clockwise else 0} {x:.1f} {y:.1f}'
    f.parts.append(f'<path d="{d}" stroke="{color}" stroke-width="3" fill="none" marker-end="url(#{mid})"/>')

def angle(coords, centre, atom):
    cx, cy = coords[centre]; ax, ay = coords[atom]
    return math.degrees(math.atan2(ay - cy, ax - cx))

def fischer(f, x, y, top, right, bottom, left, arm=78, size=27, ranks=None):
    f.line(x - arm, y, x + arm, y, width=3).line(x, y - arm, x, y + arm, width=3)
    f.text(x, y - arm - 14, top, size).text(x, y + arm + 34, bottom, size)
    f.text(x - arm - 12, y + 9, left, size, anchor='end').text(x + arm + 12, y + 9, right, size, anchor='start')
    if ranks:
        pos = {'top': (x + 58, y - arm - 22), 'bottom': (x + 52, y + arm + 25),
               'left': (x - arm + 10, y - 26), 'right': (x + arm - 10, y - 26)}
        for k, n in ranks.items():
            badge(f, pos[k], n, r=14)

# ------------------------------------------------------------------ figures

LACTIC_S = expect('C[C@]([H])(O)C(=O)O', atoms={1: 'S'})   # atoms: 0 CH3, 1 C*, 2 H, 3 O, 4 C(OOH)
LACTIC_R = expect('C[C@@]([H])(O)C(=O)O', atoms={1: 'R'})

# 1. Three kinds of isomer pair.
f = Figure('Constitutional isomers, enantiomers and diastereomers',
           'Ethanol and dimethyl ether; (R)- and (S)-lactic acid; (Z)- and (E)-but-2-ene', 960, 460)
cols = [
    ('CCO', 'COC', 'constitutional isomers', 'different connectivity', {}, {}),
    (LACTIC_S, LACTIC_R, 'enantiomers', 'mirror images', {}, dict(flip=True)),
    (expect('C/C=C\\C', bonds=['Z']), expect('C/C=C/C', bonds=['E']), 'diastereomers',
     'stereoisomers, not mirrors', {}, {}),
]
for i, (a, b, name, sub, ka, kb) in enumerate(cols):
    x = 160 + 320*i
    if i: f.line(x - 160, 30, x - 160, 430, color='#c9d6dc', width=2)
    pl = place_h if '[H]' in a else place
    pl(f, a, x, 95, bond_px=46, font_px=22, **ka)
    pl(f, b, x, 245, bond_px=46, font_px=22, **kb)
    f.text(x, 375, name, 24, weight='bold').text(x, 410, sub, 21, color=BLUE)
save('relationships', f)

# 2. A stereocentre needs four different groups.
f = Figure('Propan-2-ol has no stereocentre; lactic acid has one',
           'Propan-2-ol carries two identical methyl groups; lactic acid carries OH, H, CH3 and COOH', 960, 400)
place(f, 'CC(O)C', 250, 150, bond_px=72, highlight_atoms=[1])
place(f, 'CC(O)C(=O)O', 700, 150, bond_px=72, highlight_atoms=[1])
f.text(250, 290, 'propan-2-ol', 25, weight='bold').text(250, 328, 'H, OH, CH₃, CH₃', 22, color=BLUE)
f.text(250, 364, 'two groups the same → not a stereocentre', 21, color=GREY)
f.text(700, 290, 'lactic acid', 25, weight='bold').text(700, 328, 'H, OH, CH₃, COOH', 22, color=BLUE)
f.text(700, 364, 'four different groups → stereocentre', 21, color=GREEN)
save('stereocentre', f)

# 3. Spot the stereocentres (drawn without stereo on purpose: it's a question).
f = Figure('Which of these molecules has a stereocentre?',
           'Glycine, alanine, glycerol and ibuprofen drawn as plain skeletal structures', 960, 360)
for x, smi, name, kw in [(115, 'NCC(=O)O', 'glycine', {}), (345, 'CC(N)C(=O)O', 'alanine', {}),
                         (570, 'OCC(O)CO', 'glycerol', {}),
                         (795, 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'ibuprofen', dict(bond_px=30, font_px=20))]:
    place(f, smi, x, 150, **{'bond_px': 48, 'font_px': 23, **kw})
    f.text(x, 315, name, 25, weight='bold')
save('find-centres', f)

# 4. Enantiomers: lactic acid and its mirror image.
f = Figure('(S)- and (R)-lactic acid are mirror images',
           'Wedge-and-hash drawings of (S)-lactic acid and its mirror image (R)-lactic acid either side of a dashed mirror line',
           960, 400)
place_h(f, LACTIC_S, 260, 165, bond_px=78)
place_h(f, LACTIC_R, 700, 165, bond_px=78, flip=True)
f.line(480, 40, 480, 300, color=BLUE, width=2, dash='8 8').text(480, 330, 'mirror', 21, color=BLUE)
f.text(260, 330, '(S)-lactic acid', 25, weight='bold').text(260, 366, 'the one your muscles make', 21, color=GREY)
f.text(700, 330, '(R)-lactic acid', 25, weight='bold').text(700, 366, 'same atoms, same bonds', 21, color=GREY)
save('enantiomers', f)

# 5. CIP rule 1: atomic number of the attached atom.
BRCL = expect('O[C@@]([H])(Br)Cl', atoms={1: 'R'})
BRPROP = expect('C[C@@]([H])(O)CBr', atoms={1: 'R'})
f = Figure('Priority starts with atomic number of the attached atom',
           'Bromochloromethanol ranked Br > Cl > O > H; 1-bromopropan-2-ol ranked OH > CH2Br > CH3 > H',
           960, 440)
c = place_h(f, BRCL, 250, 175, bond_px=66)
priorities(f, c, 1, {3: 1, 4: 2, 0: 3, 2: 4}, extra={0: 18, 3: -8, 4: -8})
f.text(250, 385, 'Br > Cl > O > H', 25, weight='bold').text(250, 420, 'Z = 35 > 17 > 8 > 1', 22, color=BLUE)
c = place_h(f, BRPROP, 700, 180, bond_px=70)
priorities(f, c, 1, {3: 1, 4: 2, 0: 3, 2: 4}, extra={3: 18})
f.text(700, 385, 'OH > CH₂Br > CH₃ > H', 25, weight='bold')
f.text(700, 420, 'CH₂Br is attached through C, so O wins', 22, color=BLUE)
save('cip-rule1', f)

# 6. Breaking a tie: lactic acid, COOH versus CH3.
f = Figure('Breaking the carbon tie in lactic acid',
           'Lactic acid with priorities OH 1, COOH 2, CH3 3, H 4; the COOH carbon carries (O, O, O) and the CH3 carbon (H, H, H)',
           960, 420)
c = place_h(f, LACTIC_S, 250, 190, bond_px=80, highlight_atoms=[0, 4])
priorities(f, c, 1, {3: 1, 4: 2, 0: 3, 2: 4}, extra={3: 18})
f.text(700, 120, 'Both attach through C: a tie.', 24, weight='bold')
f.text(700, 175, 'COOH carbon carries  (O, O, O)', 24, color=RED)
f.text(700, 210, 'C=O counts O twice', 20, color=GREY)
f.text(700, 265, 'CH₃ carbon carries  (H, H, H)', 24, color=INK)
f.text(700, 330, 'First difference: O beats H', 24, color=BLUE, weight='bold')
f.text(700, 366, 'so COOH ranks above CH₃', 22, color=BLUE)
save('cip-ties', f)

# 7. Assign R/S with the lowest group pointing away (hashed H).
f = Figure('Assigning (S)-lactic acid with H pointing away',
           'With H hashed away from the viewer, OH to COOH to CH3 runs counterclockwise, so the centre is S', 960, 420)
c = place_h(f, LACTIC_S, 270, 200, bond_px=88)
assert sense(c, 3, 4, 0) == 'counterclockwise'      # H (atom 2) is hashed in this layout
priorities(f, c, 1, {3: 1, 4: 2, 0: 3, 2: 4}, extra={3: 18})
arc_arrow(f, c[1], [angle(c, 1, 3), angle(c, 1, 4), angle(c, 1, 0)], 58, clockwise=False)
for i, (t, col) in enumerate([('1  Rank the four groups', INK), ('2  Put 4 (H) at the back', INK),
                              ('3  Trace 1 → 2 → 3', INK), ('counterclockwise  →  S', BLUE)]):
    f.text(560, 130 + 62*i, t, 25 if i < 3 else 28, color=col, anchor='start', weight='bold' if i == 3 else 'normal')
f.text(560, 375, 'clockwise would be R', 21, color=GREY, anchor='start')
save('rs-assign', f)

# 8. Same molecule, H toward you: reverse the apparent answer.
f = Figure('When H points toward you, reverse the answer',
           'The same (S)-lactic acid drawn twice: H hashed reads counterclockwise S directly; H wedged appears clockwise and is reversed to S',
           960, 420)
c = place_h(f, LACTIC_S, 240, 175, bond_px=78)
assert sense(c, 3, 4, 0) == 'counterclockwise'
arc_arrow(f, c[1], [angle(c, 1, 3), angle(c, 1, 4), angle(c, 1, 0)], 52, clockwise=False)
priorities(f, c, 1, {3: 1, 4: 2, 0: 3}, extra={3: 18})
f.text(240, 340, 'H away: read it directly', 23, weight='bold').text(240, 378, 'counterclockwise → S', 23, color=BLUE)
c = place_h(f, LACTIC_S, 720, 175, bond_px=78, flip=True)   # mirrored layout, same S molecule: H now wedged
assert sense(c, 3, 4, 0) == 'clockwise'
arc_arrow(f, c[1], [angle(c, 1, 3), angle(c, 1, 4), angle(c, 1, 0)], 52, clockwise=True)
priorities(f, c, 1, {3: 1, 4: 2, 0: 3}, extra={3: 18})
f.text(720, 340, 'H toward you: reverse', 23, weight='bold').text(720, 378, 'looks clockwise → really S', 23, color=BLUE)
f.line(480, 50, 480, 300, color='#c9d6dc', width=2)
save('h-forward', f)

# 9. Fischer projection and its 3D meaning.
assert lactic_fischer('COOH', 'H', 'CH₃', 'OH') == 'S'
f = Figure('A Fischer projection and the 3D shape it stands for',
           'Fischer projection of (S)-lactic acid next to the same molecule with horizontal bonds as wedges and vertical bonds as hashes',
           960, 420)
fischer(f, 220, 190, 'COOH', 'H', 'CH₃', 'OH')
f.text(480, 200, '=', 40, color=BLUE)
x, y = 720, 190
f.wedge(x, y, x - 92, y, spread=11).wedge(x, y, x + 92, y, spread=11)
f.wedge(x, y, x, y - 82, spread=10, hashed=True).wedge(x, y, x, y + 82, spread=10, hashed=True)
f.text(x, y - 96, 'COOH', 27).text(x, y + 116, 'CH₃', 27)
f.text(x - 104, y + 9, 'OH', 27, anchor='end').text(x + 104, y + 9, 'H', 27, anchor='start')
f.text(220, 378, 'Fischer projection', 23, weight='bold')
f.text(720, 378, 'horizontal = toward you · vertical = away', 22, color=BLUE)
save('fischer', f)

# 10. R/S straight from a Fischer projection.
f = Figure('Assigning R/S on a Fischer projection',
           'Fischer projection of lactic acid with OH 1 on the left, COOH 2 at the top, CH3 3 at the bottom and H 4 on the right; the trace looks clockwise but H is horizontal, so the centre is S',
           960, 420)
fischer(f, 250, 190, 'COOH', 'H', 'CH₃', 'OH', arm=86, ranks={'left': 1, 'top': 2, 'bottom': 3, 'right': 4})
arc_arrow(f, (250, 190), [180, 270, 90], 46, clockwise=True)
for i, (t, col) in enumerate([('1 → 2 → 3 looks clockwise', INK), ('but H (4) is horizontal,', INK),
                              ('so it points toward you', INK), ('reverse:  S', BLUE)]):
    f.text(560, 135 + 58*i, t, 25 if i < 3 else 30, color=col, anchor='start', weight='bold' if i == 3 else 'normal')
f.text(560, 375, 'this is L-lactic acid, the muscle form', 21, color=GREY, anchor='start')
save('fischer-rs', f)

# 11. Legal and illegal Fischer moves.
moves = [('starting point', ('COOH', 'H', 'CH₃', 'OH')),
         ('rotate 180°', ('CH₃', 'OH', 'COOH', 'H')),
         ('rotate 90°', ('OH', 'COOH', 'H', 'CH₃')),
         ('swap OH and H', ('COOH', 'OH', 'CH₃', 'H'))]
labels = [lactic_fischer(*g) for _, g in moves]
assert labels == ['S', 'S', 'R', 'R'], labels
f = Figure('What you may and may not do to a Fischer projection',
           'Starting from (S)-lactic acid: a 180 degree rotation keeps S; a 90 degree rotation or one swap gives R', 960, 400)
for i, ((name, g), lab) in enumerate(zip(moves, labels)):
    x = 120 + 240*i
    fischer(f, x, 160, *g, arm=46, size=22)
    f.text(x, 315, name, 23, weight='bold')
    same = lab == 'S'
    f.text(x, 355, f'({lab}): ' + ('same molecule' if same else 'enantiomer'), 22, color=GREEN if same else RED)
save('fischer-moves', f)

# 12. Two stereocentres: the four threonines (C2 = atom 3, C3 = atom 1).
THR = [('C[C@@H](O)[C@H](N)C(=O)O', '(2S,3R)', 'L-threonine, in proteins', {}),
       ('C[C@H](O)[C@@H](N)C(=O)O', '(2R,3S)', 'its enantiomer', dict(flip=True)),
       ('C[C@H](O)[C@H](N)C(=O)O', '(2S,3S)', 'L-allothreonine', {}),
       ('C[C@@H](O)[C@@H](N)C(=O)O', '(2R,3R)', 'its enantiomer', dict(flip=True))]
for smi, lab, _, _ in THR:
    a, _ = cip(smi)
    assert f'(2{a[3]},3{a[1]})' == lab, (smi, a, lab)
f = Figure('Threonine has two stereocentres and four stereoisomers',
           'The four stereoisomers of threonine: (2S,3R), (2R,3S), (2S,3S) and (2R,3R), arranged as two pairs of enantiomers',
           960, 485)
for i, (smi, lab, sub, kw) in enumerate(THR):
    x = 210 + 540*(i % 2); y = 95 + 240*(i//2)
    place(f, smi, x, y, bond_px=42, font_px=22, **kw)
    f.text(x, y + 117, lab, 23, weight='bold').text(x, y + 140, sub, 20, color=GREY)
for y in (95, 335):
    f.line(375, y, 585, y, color=BLUE, width=2.5, arrow=True, both=True)
    f.text(480, y - 12, 'enantiomers', 20, color=BLUE)
f.line(480, 135, 480, 290, color=RED, width=2.5, arrow=True, both=True)
f.text(492, 215, 'diastereomers', 20, color=RED, anchor='start')
save('threonine', f)

# 13. Tartaric acid: (R,R), (S,S) and meso, as Fischer projections.
def tartaric_labels(upper_left, lower_left):
    """Each centre separately; the other half is treated as a substituent."""
    up = fischer_label('C(=O)O', 'O' if upper_left == 'H' else '[H]', 'C(O)C(=O)O',
                       'O' if upper_left == 'OH' else '[H]')
    lo = fischer_label('C(O)C(=O)O', 'O' if lower_left == 'H' else '[H]', 'C(=O)O',
                       'O' if lower_left == 'OH' else '[H]')
    return up, lo
TART = [('H', 'OH'), ('OH', 'H'), ('H', 'H')]
tl = [tartaric_labels(*t) for t in TART]
assert tl == [('R', 'R'), ('S', 'S'), ('R', 'S')], tl
expect('O[C@@H](C(=O)O)[C@@H](O)C(=O)O', atoms={1: 'R', 5: 'R'})
expect('O[C@H](C(=O)O)[C@@H](O)C(=O)O', atoms={1: 'S', 5: 'R'})
f = Figure('Tartaric acid: two chiral forms and one meso form',
           'Fischer projections of (R,R)-, (S,S)- and meso-tartaric acid; the meso form has a mirror plane between its two halves',
           960, 440)
names = [('(2R,3R)', 'chiral, from grapes'), ('(2S,3S)', 'its enantiomer'), ('meso (2R,3S)', 'mirror plane → achiral')]
for i, ((ul, ll), (lab, sub)) in enumerate(zip(TART, names)):
    x, y = 170 + 310*i, 185
    other = lambda s: 'H' if s == 'OH' else 'OH'
    f.line(x, y - 110, x, y + 110, width=3)
    f.line(x - 62, y - 50, x + 62, y - 50, width=3).line(x - 62, y + 50, x + 62, y + 50, width=3)
    f.text(x, y - 124, 'COOH', 25).text(x, y + 142, 'COOH', 25)
    f.text(x - 72, y - 41, ul, 25, anchor='end').text(x + 72, y - 41, other(ul), 25, anchor='start')
    f.text(x - 72, y + 59, ll, 25, anchor='end').text(x + 72, y + 59, other(ll), 25, anchor='start')
    if i == 2:
        f.line(x - 125, y, x + 125, y, color=BLUE, width=2.5, dash='9 7')
    f.text(x, 375, lab, 24, weight='bold').text(x, 410, sub, 21, color=GREEN if i == 2 else GREY)
save('tartaric', f)

# 14. E/Z in but-2-ene.
Z = expect('C/C([H])=C([H])\\C', bonds=['Z']); E = expect('C/C([H])=C([H])/C', bonds=['E'])
f = Figure('(Z)- and (E)-but-2-ene', 'Z-but-2-ene has both methyl groups on the same side; E-but-2-ene has them on opposite sides',
           960, 400)
place_h(f, Z, 250, 150, bond_px=85, highlight_atoms=[0, 5])
place_h(f, E, 710, 150, bond_px=85, highlight_atoms=[0, 5])
f.text(250, 300, '(Z)-but-2-ene', 25, weight='bold').text(250, 338, 'CH₃ > H at each end · same side', 21, color=BLUE)
f.text(710, 300, '(E)-but-2-ene', 25, weight='bold').text(710, 338, 'CH₃ > H at each end · opposite sides', 21, color=BLUE)
f.text(480, 382, 'Z = zusammen (together)  ·  E = entgegen (opposite)', 21, color=GREY)
save('but-2-ene', f)

# 15. When does E/Z apply?
EPENT = expect('C/C([H])=C(\\C)CC', bonds=['E'])
for smi in ('[H]C([H])=C([H])C', 'CC([H])=C(C)C'):
    assert not cip(smi)[1], smi
f = Figure('Does this alkene have E and Z forms?',
           'Propene and 2-methylbut-2-ene have two identical groups on one alkene carbon; (E)-3-methylpent-2-ene has two different groups on both',
           960, 400)
for x, smi, name, sub, col, hl in [
        (165, '[H]C([H])=C([H])C', 'propene', 'H, H on one end → no', RED, [0, 2]),
        (480, 'CC([H])=C(C)C', '2-methylbut-2-ene', 'CH₃, CH₃ on one end → no', RED, [4, 5]),
        (790, EPENT, '(E)-3-methylpent-2-ene', 'ethyl > methyl → E', GREEN, [0, 5, 6])]:
    place_h(f, smi, x, 150, bond_px=62, highlight_atoms=hl, font_px=24)
    f.text(x, 300, name, 23, weight='bold').text(x, 338, sub, 21, color=col)
save('ez-test', f)

# 16. Carvone: the two smells.
RCARV = expect('CC1=CC[C@H](CC1=O)C(C)=C', atoms={4: 'R'})
SCARV = expect('CC1=CC[C@@H](CC1=O)C(C)=C', atoms={4: 'S'})
f = Figure('(R)- and (S)-carvone', '(R)-(−)-carvone smells of spearmint; its enantiomer (S)-(+)-carvone smells of caraway',
           960, 420)
place(f, RCARV, 250, 160, bond_px=62, highlight_atoms=[4])
place(f, SCARV, 710, 160, bond_px=62, highlight_atoms=[4], flip=True)
f.line(480, 40, 480, 300, color=BLUE, width=2, dash='8 8')
f.text(250, 330, '(R)-(−)-carvone', 25, weight='bold').text(250, 368, 'spearmint', 23, color=GREEN)
f.text(710, 330, '(S)-(+)-carvone', 25, weight='bold').text(710, 368, 'caraway (rye bread)', 23, color=GREEN)
save('carvone', f)

# 17. Enantiomeric excess as bars.
f = Figure('Enantiomeric excess as a bar',
           'Bars for 50:50, 70:30 and 95:5 R:S mixtures; the matched R/S pairs form a racemic part and the leftover R is the excess',
           960, 420)
X0, W = 230, 560
for i, (r, s) in enumerate([(50, 50), (70, 30), (95, 5)]):
    y = 70 + 105*i
    ee = r - s
    f.text(X0 - 20, y + 30, f'{r} : {s}', 26, anchor='end', weight='bold')
    if ee: f.rect(X0, y, W*ee/100, 44, fill=BLUE)
    f.rect(X0 + W*ee/100, y, W*(2*s)/100, 44, fill='#d4dee3')
    if ee > 30: f.text(X0 + W*ee/200, y + 30, f'excess R {ee}%', 21, color='white', weight='bold')
    f.text(X0 + W*ee/100 + W*s/100, y + 30, f'racemic {2*s}%' if s >= 15 else '', 21, color=INK)
    f.text(X0 + W + 14, y + 30, f'ee {ee}%', 23, color=BLUE, anchor='start', weight='bold')
f.text(480, 395, 'Pair every S with one R. Whatever R is left over is the excess.', 22, color=GREY)
save('ee-bars', f)

print('Built 17 stereochemistry figures; all R/S and E/Z labels checked with RDKit.')
