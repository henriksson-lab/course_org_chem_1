"""Skeletal-structure and mechanism figures for the Day 11 substitution lecture (RDKit via molfig).

Curved arrows and lone pairs use molfig's curved_arrow/lone_pair at the atom
coordinates returned by place(). Stereodescriptors in the labels are checked
against RDKit's CIP labeller when the script runs (check_cip), so a wrong @/@@
fails the build.
"""
import math
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler, rdDepictor
from chemfig import Figure, BLUE, RED, INK, GREEN
import molfig
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'substitution'/'figures'


class _KeepH:
    """molfig parses SMILES with RDKit defaults, which strip explicit [H] atoms.
    Mechanism drawings need particular H atoms shown (a beta-H, the H on a
    stereocentre, the H that moves in a proton transfer), so parse with
    removeHs=False for this script only."""
    def __getattr__(self, name):
        return getattr(Chem, name)

    @staticmethod
    def MolFromSmiles(smiles):
        p = Chem.SmilesParserParams()
        p.removeHs = False
        return Chem.MolFromSmiles(smiles, p)


class _Layout:
    """molfig always computes its own 2D layout. For SN2/SN1 panels we need the
    textbook geometry (leaving group straight out to one side, the other three
    groups splayed away from it), so this script can hand molfig fixed
    coordinates (RDKit units, bond length 1.5). RDKit still chooses the wedges
    from these coordinates and the @/@@ tags."""
    coords = None

    def __getattr__(self, name):
        return getattr(rdDepictor, name)

    def Compute2DCoords(self, mol, canonOrient=True):
        if self.coords is None:
            return rdDepictor.Compute2DCoords(mol, canonOrient=canonOrient)
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, (x, y) in self.coords.items():
            conf.SetAtomPosition(i, (x, y, 0.0))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)


molfig.Chem = _KeepH()
_layout = _Layout()
molfig.rdDepictor = _layout


def place_at(fig, smiles, x, y, coords, **kw):
    """place() with a fixed layout: coords {atom: (x, y)} at RDKit scale, y up."""
    _layout.coords = coords
    try:
        return place(fig, smiles, x, y, **kw)
    finally:
        _layout.coords = None


def polar(origin, angle_deg, length=1.5):
    a = math.radians(angle_deg)
    return (origin[0] + length*math.cos(a), origin[1] + length*math.sin(a))


def check_cip(smiles, atom, expected):
    p = Chem.SmilesParserParams(); p.removeHs = False
    m = Chem.MolFromSmiles(smiles, p)
    rdCIPLabeler.AssignCIPLabels(m)
    got = m.GetAtomWithIdx(atom).GetProp('_CIPCode')
    assert got == expected, f'{smiles}: atom {atom} is {got}, label says {expected}'


def save(name, f):
    f.save(OUT/f'sk-{name}.svg')


def mid(a, b):
    return ((a[0] + b[0])/2, (a[1] + b[1])/2)


def lp(fig, at, angle, distance=24):
    """Draw a lone pair and return its centre, the tail of an arrow."""
    lone_pair(fig, at, angle, distance)
    a = math.radians(angle)
    return (at[0] + distance*math.cos(a), at[1] - distance*math.sin(a))


def push(fig, start, end, bend=0.35, shorten=10):
    curved_arrow(fig, start, end, bend=bend, shorten=shorten)


def react_arrow(fig, x1, x2, y, label=''):
    fig.line(x1, y, x2, y, color=BLUE, width=3, arrow=True)
    if label:
        fig.text((x1 + x2)/2, y - 14, label, 20, color=BLUE)


def plus(fig, x, y):
    fig.text(x, y + 10, '+', 32)


B = 46   # bond_px: molfig draws bonds about 1.5 × this long
FS = 26  # atom label size
M = dict(bond_px=B, font_px=FS)
BR = '#7a3b2e'
RULE = '#c9d8df'
C0 = (0, 0)


def to_atom(fig, a, b, side=1, reach=34):
    """Arrow for a bond pair moving onto atom b: from the a–b bond midpoint,
    bulging to one side, ending just past b."""
    d = math.hypot(b[0] - a[0], b[1] - a[1])
    vx, vy = (b[0] - a[0])/d, (b[1] - a[1])/d
    nx, ny = -vy*side, vx*side
    m = mid(a, b)
    start = (m[0] + nx*5, m[1] + ny*5)
    end = (b[0] + vx*4 + nx*reach, b[1] + vy*4 + ny*reach)
    ex, ey = end[0] - start[0], end[1] - start[1]
    e = math.hypot(ex, ey)
    s = 1 if (ey/e)*nx - (ex/e)*ny > 0 else -1
    curved_arrow(fig, start, end, bend=0.45*s, shorten=2)


def to_bond(fig, a, b, c, side=1):
    """Arrow for the a–b bond pair moving to become a b–c pi bond."""
    s, e = mid(a, b), mid(b, c)
    curved_arrow(fig, s, e, bend=0.6*side, shorten=6)


def umbrella(centre, groups, lengths=None):
    lengths = lengths or {}
    return {atom: polar(centre, ang, lengths.get(atom, 1.5)) for atom, ang in groups.items()}


def energy_axes(f, ox, oy, top, right):
    f.line(ox, oy, ox, top, arrow=True).line(ox, oy, right, oy, arrow=True)
    f.parts.append(f'<text x="{ox - 22}" y="{(oy + top)/2}" font-family="DejaVu Sans, sans-serif" font-size="22" '
                   f'fill="{INK}" text-anchor="middle" transform="rotate(-90 {ox - 22} {(oy + top)/2})">energy</text>')
    f.text((ox + right)/2, oy + 32, 'reaction progress', 22)


def step_label(f, y, label):
    f.text(30, y, label, 23, color=BLUE, anchor='start', weight='bold')


# ---------------------------------------------------------------------------
# 1. The three parts of a substitution, in an ether synthesis.
f = Figure('Ethoxide replaces bromide on bromoethane',
           'Ethoxide is the nucleophile, the CH2 carbon of bromoethane is the electrophilic carbon and bromide is the leaving group.',
           960, 340)
y = 150
nu = place(f, 'CC[O-]', 105, y, **M)
plus(f, 210, y)
sub = place(f, 'CCBr', 330, y, highlight_atoms=[1], **M)
react_arrow(f, 440, 530, y)
place(f, 'CCOCC', 680, y, highlight_atoms=[2], **M)
plus(f, 828, y)
place(f, '[Br-]', 893, y, **M)
f.text(sub[1][0], sub[1][1] - 34, 'δ+', 22, color=RED)
f.text(sub[2][0] + 6, sub[2][1] - 36, 'δ−', 22, color=BLUE)
f.text(330, 62, 'leaving group: Br', 22, color=GREEN, weight='bold')
for x, label, sub_, col in [(105, 'nucleophile', 'lone-pair donor', BLUE), (320, 'electrophilic C', 'bonded to Br', RED),
                            (680, 'diethyl ether', 'new C–O bond', INK), (893, 'bromide', 'kept both e⁻', INK)]:
    f.text(x, 255, label, 24, color=col, weight='bold').text(x, 287, sub_, 20, color=BLUE)
save('parts', f)

# ---------------------------------------------------------------------------
# 2. Classify the carbon that carries the leaving group.
f = Figure('Methyl, primary, secondary and tertiary bromides',
           'Bromomethane, bromoethane, 2-bromopropane and 2-bromo-2-methylpropane with the carbon bearing Br shaded.',
           960, 340)
items = [('CBr', 0, 'methyl', '0 C neighbours'),
         ('CCBr', 1, 'primary (1°)', '1 C neighbour'),
         ('CC(C)Br', 1, 'secondary (2°)', '2 C neighbours'),
         ('CC(C)(C)Br', 1, 'tertiary (3°)', '3 C neighbours')]
for i, (smi, c, label, sub_) in enumerate(items):
    x = 120 + 240*i
    place(f, smi, x, 130, highlight_atoms=[c], bond_px=50, font_px=28)
    f.text(x, 262, label, 25, weight='bold').text(x, 298, sub_, 21, color=BLUE)
save('substrates', f)

# ---------------------------------------------------------------------------
# 3. SN2: two arrows, one step (methyl bromide, all H drawn).
METHYL_BR = {1: C0, 4: polar(C0, 0), **umbrella(C0, {0: 115, 2: 245, 3: 200}, {3: 1.15})}
METHYL_FLIP = {1: C0, 4: polar(C0, 180), **umbrella(C0, {0: 65, 2: 295, 3: 340}, {3: 1.15})}
f = Figure('SN2: hydroxide displaces bromide from bromomethane in one step',
           'A hydroxide lone pair attacks carbon from the side opposite bromine while the C-Br bonding pair moves onto bromine, giving methanol and bromide.',
           960, 360)
y = 165
oh = place(f, '[OH-]', 80, y, **M)
mb = place_at(f, '[H]C([H])([H])Br', 320, y, METHYL_BR, **M)
tail = lp(f, oh[0], 0, 34)
push(f, tail, mb[1], bend=0.3, shorten=12)
to_atom(f, mb[1], mb[4], side=-1)
react_arrow(f, 470, 560, y)
place_at(f, '[H]C([H])([H])O', 690, y, METHYL_FLIP, **M)
plus(f, 815, y)
place(f, '[Br-]', 885, y, **M)
f.text(240, 290, 'HO⁻ attacks the back of C', 23, color=RED)
f.text(240, 324, 'while Br⁻ leaves from the front', 23, color=RED)
f.text(760, 290, 'one step · no intermediate', 23, weight='bold')
f.text(760, 324, 'charge −1 on both sides', 22, color=BLUE)
save('sn2-arrows', f)

# ---------------------------------------------------------------------------
# 4. SN2 energy profile with the transition state.
f = Figure('SN2 energy profile: one step, one transition state',
           'Energy against reaction progress rises to a single maximum, the transition state with partial O-C and C-Br bonds, then falls to products.',
           960, 400)
energy_axes(f, 90, 340, 40, 600)
f.path('M 105 230 C 200 230 245 95 330 95 C 415 95 470 270 590 270', color=BLUE, width=4)
f.text(105, 262, 'HO⁻ + CH₃Br', 21, anchor='start')
f.text(590, 302, 'CH₃OH + Br⁻', 21, anchor='end')
f.text(330, 78, 'transition state ‡', 22, color=RED, weight='bold')
f.line(185, 230, 185, 95, color=INK, width=2, dash='6 5').line(185, 95, 320, 95, color=INK, width=2, dash='6 5')
f.text(175, 172, 'Eₐ', 24, anchor='end')
cx, cy = 790, 205
f.text(cx, 70, 'the structure at the top', 22, weight='bold')
for dx, dy in [(0, -62), (-48, 40), (48, 40)]:
    f.line(cx + dx*0.3, cy + dy*0.3, cx + dx*0.78, cy + dy*0.78, width=3)
    f.text(cx + dx, cy + dy + 9, 'H', 26)
f.text(cx, cy + 9, 'C', 28)
f.line(cx - 22, cy, cx - 88, cy, width=3, dash='7 6').line(cx + 22, cy, cx + 88, cy, width=3, dash='7 6')
f.text(cx - 118, cy + 10, 'HO', 28, color=RED).text(cx - 118, cy - 26, 'δ−', 20, color=BLUE)
f.text(cx + 115, cy + 10, 'Br', 28, color=BR).text(cx + 115, cy - 26, 'δ−', 20, color=BLUE)
f.text(cx - 165, cy + 20, '[', 64).text(cx + 162, cy + 20, ']', 64).text(cx + 184, cy - 26, '‡', 26)
f.text(cx, 310, 'C–O half made, C–Br half broken', 21, color=BLUE)
f.text(cx, 340, 'three H flat, like an umbrella', 21, color=BLUE)
f.text(cx, 370, 'caught turning inside out', 21, color=BLUE)
save('sn2-profile', f)

# ---------------------------------------------------------------------------
# 5. Inversion: (S)-2-bromobutane to (R)-butan-2-ol.
S_BROMO = 'C[C@]([H])(Br)CC'      # 0 CH3, 1 C*, 2 H, 3 Br, 4 CH2, 5 CH3
R_OL = 'C[C@@]([H])(O)CC'         # same atom order, opposite tag after backside attack
check_cip(S_BROMO, 1, 'S')
check_cip(R_OL, 1, 'R')
f = Figure('Backside attack inverts (S)-2-bromobutane to (R)-butan-2-ol',
           'Hydroxide attacks C2 of (S)-2-bromobutane opposite bromine. The three other groups flip like an umbrella, giving (R)-butan-2-ol. CIP priority numbers are shown at each group.',
           960, 400)
y = 170
co = {1: C0, 3: polar(C0, 0), 0: polar(C0, 110), 2: polar(C0, 200, 1.2), 4: polar(C0, 245)}
co[5] = polar(co[4], 185)
oh = place(f, '[OH-]', 75, y - 30, **M)
r = place_at(f, S_BROMO, 300, y, co, numbers={3: 1, 4: 2, 0: 3, 2: 4}, **M)
tail = lp(f, oh[0], 0, 34)
push(f, tail, r[1], bend=0.25, shorten=12)
to_atom(f, r[1], r[3], side=-1)
react_arrow(f, 470, 560, y)
co2 = {1: C0, 3: polar(C0, 180), 0: polar(C0, 70), 2: polar(C0, 340, 1.2), 4: polar(C0, 295)}
co2[5] = polar(co2[4], 355)
place_at(f, R_OL, 710, y, co2, numbers={3: 1, 4: 2, 0: 3, 2: 4}, **M)
plus(f, 850, y)
place(f, '[Br-]', 905, y, **M)
f.text(300, 330, '(S)-2-bromobutane', 25, weight='bold')
f.text(710, 330, '(R)-butan-2-ol', 25, weight='bold')
f.text(480, 375, 'blue numbers: CIP priorities · Br and OH are both priority 1 here', 21, color=BLUE)
save('inversion', f)

# ---------------------------------------------------------------------------
# 6. Crowding at the back of the carbon.
f = Figure('Groups behind the carbon block backside attack',
           'Bromomethane, bromoethane, 2-bromopropane and 2-bromo-2-methylpropane drawn with Br to the right and the three other groups on the attack side. SN2 slows from methyl to secondary and stops at tertiary.',
           960, 340)
crowd = [('[H]C([H])([H])Br', 'methyl', 'fast', BLUE),
         ('CC([H])([H])Br', 'primary', 'slower', BLUE),
         ('CC(C)([H])Br', 'secondary', 'slow', BLUE),
         ('CC(C)(C)Br', 'tertiary', 'no SN2', RED)]
for i, (smi, label, rate, col) in enumerate(crowd):
    x = 135 + 235*i
    a = place_at(f, smi, x + 15, 135, {1: C0, 4: polar(C0, 0), **umbrella(C0, {0: 120, 2: 240, 3: 180}, {3: 1.2})},
                 highlight_atoms=[1], bond_px=40, font_px=24)
    f.line(a[1][0] - 128, a[1][1] - 56, a[1][0] - 62, a[1][1] - 28, color=col, width=3, dash='8 6', arrow=(col == BLUE))
    if col == RED:
        f.text(a[1][0] - 70, a[1][1] - 14, '✕', 34, color=RED)
    f.text(x, 250, label, 25, weight='bold').text(x, 286, rate, 23, color=col)
f.text(480, 326, 'dashed arrow: the path the nucleophile needs, straight at the back of C', 20, color=BLUE)
save('steric', f)

# ---------------------------------------------------------------------------
# 7. Activating an alcohol: protonate, then let bromide do SN2.
f = Figure('Protonating ethanol turns OH into a water leaving group',
           'Step 1: ethanol takes a proton from hydronium to give the ethyloxonium ion. Step 2: bromide attacks carbon in an SN2 step and neutral water leaves, giving bromoethane.',
           960, 420)
y1, y2 = 115, 320
step_label(f, 34, '1  protonation: OH becomes OH₂⁺')
et = place(f, 'CCO', 120, y1, **M)
h3o = place_at(f, '[H][OH2+]', 315, y1 + 18, {0: (-1.3, 0), 1: C0}, **M)
t = lp(f, et[2], 55, 28)
push(f, t, h3o[0], bend=0.35, shorten=10)
to_atom(f, h3o[0], h3o[1], side=1)
react_arrow(f, 430, 520, y1)
place(f, 'CC[OH2+]', 650, y1, **M)
plus(f, 790, y1)
place(f, 'O', 860, y1, **M)
f.text(650, y1 + 75, 'ethyloxonium ion', 21, color=BLUE)
f.line(40, 212, 920, 212, color=RULE, width=2)
step_label(f, 246, '2  SN2 by bromide: water leaves')
br = place(f, '[Br-]', 90, y2 + 25, **M)
ox = place(f, 'CC[OH2+]', 270, y2, **M)
t = lp(f, br[0], 20, 30)
push(f, t, ox[1], bend=0.3, shorten=12)
to_atom(f, ox[1], ox[2], side=-1)
react_arrow(f, 430, 520, y2)
place(f, 'CCBr', 650, y2, **M)
plus(f, 790, y2)
place(f, 'O', 860, y2, **M)
f.text(650, y2 + 72, 'bromoethane', 21, color=BLUE)
f.text(860, y2 + 72, 'water', 21, color=BLUE)
save('activation', f)

# ---------------------------------------------------------------------------
# 8. A neutral nucleophile: ammonia, then a separate proton transfer.
f = Figure('Ammonia gives methylammonium first, then a second NH3 removes H+',
           'Step 1: ammonia attacks bromomethane to give methylammonium and bromide. Step 2: another ammonia takes a proton from methylammonium, giving methylamine and ammonium.',
           960, 420)
y1, y2 = 120, 320
step_label(f, 30, '1  SN2: N gains a fourth bond and a + charge')
n1 = place(f, 'N', 80, y1, **M)
mb = place_at(f, '[H]C([H])([H])Br', 300, y1, METHYL_BR, **M)
t = lp(f, n1[0], 0, 48)
push(f, t, mb[1], bend=0.3, shorten=12)
to_atom(f, mb[1], mb[4], side=-1)
react_arrow(f, 440, 530, y1)
place_at(f, '[H]C([H])([H])[NH3+]', 660, y1, METHYL_FLIP, **M)
plus(f, 795, y1)
place(f, '[Br-]', 865, y1, **M)
f.text(555, y1 + 72, 'methylammonium', 21, color=BLUE)
f.line(40, 212, 920, 212, color=RULE, width=2)
step_label(f, 246, '2  proton transfer to a second NH₃')
# 0 C, 1 N+, 2 H that moves, 3-4 the other N-H
ma = place_at(f, 'C[N+]([H])([H])[H]', 150, y2, {0: polar(C0, 180), 1: C0, 2: polar(C0, 0, 1.2),
                                                 3: polar(C0, 90, 1.1), 4: polar(C0, 270, 1.1)}, **M)
n2 = place(f, 'N', 385, y2, **M)
t = lp(f, n2[0], 180, 30)
push(f, t, ma[2], bend=0.4, shorten=12)
to_atom(f, ma[2], ma[1], side=1, reach=20)
react_arrow(f, 440, 530, y2)
place(f, 'CN', 650, y2, **M)
plus(f, 780, y2)
place(f, '[NH4+]', 860, y2, **M)
f.text(650, y2 + 72, 'methylamine', 21, color=BLUE)
f.text(860, y2 + 72, 'ammonium', 21, color=BLUE)
save('ammonia', f)

# ---------------------------------------------------------------------------
# 9. SN1 step 1: ionisation, and why tertiary cations are the ones that form.
TBU_BR = {1: C0, 4: polar(C0, 0), **umbrella(C0, {0: 120, 2: 240, 3: 180})}
TBU_CAT = {1: C0, **umbrella(C0, {0: 90, 2: 210, 3: 330})}
f = Figure('SN1 starts with ionisation to a carbocation',
           'Top: the C-Br bond of 2-bromo-2-methylpropane breaks with both electrons going to bromine, giving the tert-butyl cation and bromide. Bottom: carbocation stability rises from methyl to tertiary.',
           960, 420)
y1, y2 = 110, 305
tb = place_at(f, 'CC(C)(C)Br', 220, y1, TBU_BR, **M)
to_atom(f, tb[1], tb[4], side=-1)
react_arrow(f, 380, 500, y1, 'slow')
place_at(f, 'C[C+](C)C', 640, y1, TBU_CAT, **M)
plus(f, 770, y1)
place(f, '[Br-]', 840, y1, **M)
f.text(640, y1 + 92, 'tert-butyl cation: 3 bonds, flat', 21, color=BLUE)
f.line(40, 222, 920, 222, color=RULE, width=2)
cats = [('[CH3+]', 'methyl'), ('C[CH2+]', 'primary'), ('C[CH+]C', 'secondary'), ('C[C+](C)C', 'tertiary')]
for i, (smi, label) in enumerate(cats):
    x = 130 + 233*i
    place(f, smi, x, y2 - 15, bond_px=38, font_px=24)
    f.text(x, y2 + 50, label, 22, weight='bold')
f.line(110, 405, 850, 405, color=BLUE, width=3, arrow=True)
f.text(480, 395, 'more alkyl groups → more stable cation → SN1 possible', 20, color=BLUE)
save('sn1-ionise', f)

# ---------------------------------------------------------------------------
# 10. SN1 steps 2 and 3: water captures the cation, a second water removes H+.
f = Figure('Water captures the carbocation, then a second water removes a proton',
           'Step 2: a lone pair on water bonds to the cationic carbon, giving a protonated alcohol with a three-bonded, positive oxygen. Step 3: another water takes a proton, giving tert-butanol and hydronium.',
           960, 420)
y1, y2 = 115, 320
step_label(f, 34, '2  capture: O now has three bonds and a + charge')
w = place(f, 'O', 75, y1, **M)
cat = place_at(f, 'C[C+](C)C', 270, y1, TBU_CAT, **M)
t = lp(f, w[0], 0, 44)
push(f, t, cat[1], bend=0.3, shorten=14)
react_arrow(f, 420, 510, y1)
place_at(f, 'CC(C)(C)[OH2+]', 660, y1, {1: C0, 4: polar(C0, 180), **umbrella(C0, {0: 60, 2: 300, 3: 0})}, **M)
f.text(870, y1 - 6, 'oxonium ion', 21, color=BLUE).text(870, y1 + 22, '(protonated', 19, color=BLUE).text(870, y1 + 46, 'alcohol)', 19, color=BLUE)
f.line(40, 212, 920, 212, color=RULE, width=2)
step_label(f, 246, '3  deprotonation by a second water')
# 0-3 tert-butyl, 4 O+, 5 H that moves, 6 H that stays
O_POS = polar(C0, 0)
ox = place_at(f, 'CC(C)(C)[O+]([H])[H]', 150, y2, {1: C0, 4: O_POS, 5: polar(O_POS, 40, 1.2), 6: polar(O_POS, 320, 1.2),
                                                   **umbrella(C0, {0: 120, 2: 240, 3: 180})}, **M)
w2 = place(f, 'O', 380, y2 - 30, **M)
t = lp(f, w2[0], 180, 44)
push(f, t, ox[5], bend=0.3, shorten=12)
to_atom(f, ox[5], ox[4], side=1)
react_arrow(f, 470, 560, y2)
place_at(f, 'CC(C)(C)O', 690, y2, {1: C0, 4: polar(C0, 0), **umbrella(C0, {0: 120, 2: 240, 3: 180})}, **M)
plus(f, 805, y2)
place(f, '[OH3+]', 875, y2, **M)
f.text(690, y2 + 88, 'tert-butanol', 21, color=BLUE)
f.text(875, y2 + 88, 'hydronium', 21, color=BLUE)
save('sn1-capture', f)

# ---------------------------------------------------------------------------
# 11. SN1 energy profile: two intermediates, the first barrier is the highest.
f = Figure('SN1 energy profile: three steps, two intermediates',
           'Energy profile for tert-butyl bromide in water. The first barrier, ionisation, is highest. Two dips are the carbocation and the protonated alcohol, both intermediates.',
           960, 400)
energy_axes(f, 90, 340, 40, 640)
f.path('M 105 240 C 150 240 170 70 215 70 C 260 70 270 190 305 190 C 335 190 345 150 370 150 '
       'C 395 150 405 215 440 215 C 470 215 480 180 500 180 C 525 180 540 270 625 270', color=BLUE, width=4)
f.text(105, 270, 'R–Br', 21, anchor='start')
f.text(215, 55, 'ionisation (slow)', 21, color=RED, weight='bold')
f.text(305, 222, 'R⁺', 22, weight='bold')
f.text(440, 248, 'R–OH₂⁺', 22, weight='bold')
f.text(625, 302, 'R–OH', 21, anchor='end')
for x, y in [(305, 190), (440, 215)]:
    f.circle(x, y, 6, color=BLUE, fill=BLUE, width=1)
tx = 800
f.text(tx, 90, 'rate = k[R–Br]', 28, color=RED, weight='bold')
f.text(tx, 135, 'only the slow first step', 21)
f.text(tx, 163, 'sets the rate', 21)
f.text(tx, 215, 'dips = intermediates', 21, color=BLUE)
f.text(tx, 243, '(real, short-lived species)', 21, color=BLUE)
f.text(tx, 295, 'peaks = transition states', 21, color=BLUE)
f.text(tx, 323, '(never isolable)', 21, color=BLUE)
save('sn1-profile', f)

# ---------------------------------------------------------------------------
# 12. A flat cation has two faces: one enantiomer in, both out.
R_SUB = 'CCC[C@](C)(Br)CC'     # 0-2 propyl, 3 C*, 4 CH3, 5 Br, 6-7 ethyl
CAT = 'CCC[C+](C)CC'           # 0-2 propyl, 3 C+, 4 CH3, 5-6 ethyl
OL_A = 'CCC[C@](C)(O)CC'
OL_B = 'CCC[C@@](C)(O)CC'
check_cip(R_SUB, 3, 'R')
check_cip(OL_A, 3, 'R')
check_cip(OL_B, 3, 'S')


def skeleton(c, me, pr, et, extra=None):
    """3-methylhexane skeleton around carbon c, given the atom indices."""
    co = {c: C0, me: polar(C0, 90), pr[0]: polar(C0, 210), et[0]: polar(C0, 330)}
    co[pr[1]] = polar(co[pr[0]], 150)
    co[pr[2]] = polar(co[pr[1]], 210)
    co[et[1]] = polar(co[et[0]], 30)
    co.update(extra or {})
    return co


f = Figure('A planar carbocation gives both enantiomers',
           '(R)-3-bromo-3-methylhexane ionises to a flat carbocation. Water can attack either face, so both (R)- and (S)-3-methylhexan-3-ol form.',
           960, 460)
yc = 220
s = place_at(f, R_SUB, 140, yc, skeleton(3, 4, (2, 1, 0), (6, 7), {5: polar(C0, 270, 1.2)}), bond_px=38, font_px=24)
f.text(140, 400, '(R)-3-bromo-', 21, weight='bold').text(140, 426, '3-methylhexane', 21, weight='bold')
react_arrow(f, 280, 350, yc, '−Br⁻')
place_at(f, CAT, 470, yc, skeleton(3, 4, (2, 1, 0), (5, 6)), bond_px=38, font_px=24)
f.text(470, 400, 'flat cation', 21, weight='bold').text(470, 426, 'two open faces', 21, color=BLUE)
f.line(590, yc - 20, 660, yc - 95, color=BLUE, width=3, arrow=True)
f.line(590, yc + 20, 660, yc + 95, color=BLUE, width=3, arrow=True)
f.text(600, yc - 85, 'water from', 20, color=BLUE, anchor='end').text(600, yc - 60, 'the back', 20, color=BLUE, anchor='end')
f.text(600, yc + 80, 'water from', 20, color=BLUE, anchor='end').text(600, yc + 105, 'the front', 20, color=BLUE, anchor='end')
OL_CO = skeleton(3, 4, (2, 1, 0), (6, 7), {5: polar(C0, 270, 1.2)})
place_at(f, OL_A, 805, 85, OL_CO, bond_px=34, font_px=22)
place_at(f, OL_B, 805, 305, OL_CO, bond_px=34, font_px=22)
f.text(805, 170, '(R)-3-methylhexan-3-ol', 21, weight='bold')
f.text(805, 390, '(S)-3-methylhexan-3-ol', 21, weight='bold')
f.text(805, 430, 'mirror images, about 1 : 1', 21, color=BLUE)
save('sn1-stereo', f)

# ---------------------------------------------------------------------------
# 13. The same cation, two jobs for water: SN1 versus E1.
f = Figure('Water can attack the cation or remove a beta-H from it',
           'Left: water acts as a nucleophile at the cationic carbon, the SN1 path to tert-butanol. Right: water acts as a base on a beta C-H, giving 2-methylpropene and hydronium, the E1 path.',
           960, 430)
f.line(480, 30, 480, 400, color=RULE, width=2)
f.text(240, 45, 'water as nucleophile (SN1)', 23, color=BLUE, weight='bold')
f.text(720, 45, 'water as base (E1)', 23, color=RED, weight='bold')
w = place(f, 'O', 80, 150, **M)
cat = place_at(f, 'C[C+](C)C', 270, 150, TBU_CAT, **M)
t = lp(f, w[0], 0, 44)
push(f, t, cat[1], bend=0.3, shorten=14)
f.line(240, 245, 240, 300, color=BLUE, width=3, arrow=True)
place_at(f, 'CC(C)(C)[OH2+]', 240, 350, {1: C0, 4: polar(C0, 180), **umbrella(C0, {0: 60, 2: 300, 3: 0})}, bond_px=40, font_px=24)
f.text(240, 415, '→ tert-butanol after H⁺ loss', 20, color=BLUE)
beta = place_at(f, '[H]C[C+](C)C', 640, 160, {2: C0, 1: polar(C0, 330), 0: polar(polar(C0, 330), 30, 1.2),
                                                3: polar(C0, 90), 4: polar(C0, 210)}, **M)
w2 = place(f, 'O', 885, 105, **M)
t = lp(f, w2[0], 200, 44)
push(f, t, beta[0], bend=0.3, shorten=12)
to_bond(f, beta[0], beta[1], beta[2], side=1)
f.text(beta[1][0] + 4, beta[1][1] + 40, 'β-C', 20, color=BLUE)
f.line(720, 245, 720, 300, color=BLUE, width=3, arrow=True)
place(f, 'C=C(C)C', 650, 350, bond_px=40, font_px=24)
plus(f, 755, 350)
place(f, '[OH3+]', 830, 350, bond_px=40, font_px=24)
f.text(720, 415, '2-methylpropene: Day 12', 20, color=RED)
save('competition', f)

print('Built 13 substitution figures.')
