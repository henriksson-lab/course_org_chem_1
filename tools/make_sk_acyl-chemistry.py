"""Skeletal-structure figures for the Day 17 acyl chemistry lecture (RDKit via molfig).

Run from the repository root: python3 tools/make_sk_acyl-chemistry.py
Writes acyl-chemistry/figures/sk-*.svg.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rdkit import Chem
import molfig
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'acyl-chemistry'/'figures'


class _KeepExplicitH:
    """molfig parses SMILES with RDKit defaults, which drop explicit [H] atoms.

    Mechanisms need the transferred proton drawn as an atom, so keep any [H]
    written in the SMILES (local workaround; molfig itself is unchanged).
    """
    def __getattr__(self, name):
        return getattr(Chem, name)

    def MolFromSmiles(self, smiles):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        return Chem.MolFromSmiles(smiles, params)


molfig.Chem = _KeepExplicitH()


def save(name, f):
    f.save(OUT/f'sk-{name}.svg')


def mol(f, smiles, x, y, bond=None, deg=90, flip=False, anchor=None, **kw):
    """place() with an orientation: rotate so the bond a->b points at deg (0 = right, 90 = up).

    anchor: put this atom (rather than the drawing's centre) at (x, y).
    """
    rot = 0
    if bond:
        a, b = bond
        molfig.molecule_svg(smiles, flip=flip, **kw)
        c = molfig.molecule_svg.last_coords
        now = math.degrees(math.atan2(-(c[b][1] - c[a][1]), c[b][0] - c[a][0]))
        rot = deg - now
        # molfig rotates in RDKit coordinates (y up); check the sign empirically.
        molfig.molecule_svg(smiles, flip=flip, rotate=rot, **kw)
        c = molfig.molecule_svg.last_coords
        got = math.degrees(math.atan2(-(c[b][1] - c[a][1]), c[b][0] - c[a][0]))
        if abs(((got - deg + 180) % 360) - 180) > 1:
            rot = -rot
    if anchor is not None:
        _, w, h = molfig.molecule_svg(smiles, flip=flip, rotate=rot, **kw)
        ax, ay = molfig.molecule_svg.last_coords[anchor]
        x, y = x - (ax - w/2), y - (ay - h/2)
    return place(f, smiles, x, y, rotate=rot, flip=flip, **kw)


def arrow(f, x1, x2, y, label=None, lsize=21):
    f.line(x1, y, x2, y, color=BLUE, arrow=True)
    if label:
        f.text((x1 + x2)/2, y - 14, label, lsize, color=BLUE)


def equilibrium(f, x1, x2, y, label=None, lsize=21):
    f.line(x1, y - 6, x2, y - 6, color=BLUE, arrow=True)
    f.line(x2, y + 6, x1, y + 6, color=BLUE, arrow=True)
    if label:
        f.text((x1 + x2)/2, y - 20, label, lsize, color=BLUE)


def resonance(f, x1, x2, y):
    f.line(x1, y, x2, y, color=BLUE, arrow=True, both=True)


def plus(f, x, y):
    f.text(x, y + 11, '+', 34)


def mid(p, q, t=0.5):
    return (p[0] + t*(q[0] - p[0]), p[1] + t*(q[1] - p[1]))


def lp_at(at, deg, d=24):
    a = math.radians(deg)
    return (at[0] + d*math.cos(a), at[1] - d*math.sin(a))


def push_from_lp(f, at, deg, end, bend=0.35, d=24):
    """Lone pair on an atom plus a curved arrow starting from it."""
    lone_pair(f, at, deg, distance=d)
    curved_arrow(f, lp_at(at, deg, d + 4), end, bend=bend, shorten=6)


# 1. The acyl family, with an ether for contrast.
f = Figure('Acid, ester, amide and ether', 'Acetic acid, ethyl acetate and acetamide share an acyl carbon; diethyl ether has none', 960, 400)
items = [
    ('CC(=O)O', 'carboxylic acid', 'Y = OH', [1, 2, 3], (1, 2)),
    ('CC(=O)OCC', 'ester', 'Y = OR', [1, 2, 3], (1, 2)),
    ('CC(N)=O', 'amide', 'Y = NH₂', [1, 2, 3], (1, 3)),
    ('CCOCC', 'ether', 'no C=O at all', [2], None),
]
for i, (smi, fam, sub, hl, bond) in enumerate(items):
    x = 125 + 237*i
    mol(f, smi, x, 170, bond=bond, highlight_atoms=hl, bond_px=44 if bond is None else 50, font_px=26)
    f.text(x, 300, fam, 25, weight='bold').text(x, 338, sub, 22, color=BLUE)
f.text(480, 382, 'acetic acid · ethyl acetate · acetamide · diethyl ether', 20)
save('acyl-family', f)

# Tetrahedral intermediates are drawn as a cross: old carbonyl O up, CH3 down,
# the incoming group left and the group that can leave right (flip=True).

# 2. Step 1: hydroxide adds to ethyl acetate.
f = Figure('Hydroxide adds to the ester carbonyl', 'Two curved arrows: hydroxide lone pair to the acyl carbon, C=O pi pair to oxygen', 960, 380)
c = mol(f, 'CC(=O)OCC', 230, 150, bond=(1, 2), bond_px=66, font_px=28)
oh = mol(f, '[OH-]', c[1][0] - 70, c[1][1] + 120, anchor=0, bond_px=66, font_px=28)
push_from_lp(f, oh[0], 60, c[1], bend=-0.2)
curved_arrow(f, mid(c[1], c[2]), c[2], bend=0.9, shorten=4)
arrow(f, 430, 530, 185)
t = mol(f, 'CC([O-])(O)OCC', 730, 175, bond=(1, 2), flip=True, bond_px=66, font_px=28)
f.text(230, 360, 'ethyl acetate + hydroxide', 24, weight='bold')
f.text(730, 360, 'tetrahedral intermediate', 24, weight='bold')
save('addition', f)

# 3. The tetrahedral intermediate, with every oxygen labelled.
f = Figure('Tracking the three oxygens', 'Tetrahedral intermediate with oxygens a (old carbonyl, now O minus), b (still bonded to ethyl) and c (from hydroxide)', 960, 380)
t = mol(f, 'CC([O-])(O)OCC', 270, 185, bond=(1, 2), flip=True, bond_px=80, font_px=32, highlight_atoms=[1])
f.text(t[2][0] + 48, t[2][1] + 10, 'a', 30, color=GREEN, weight='bold')
f.text(t[3][0] - 62, t[3][1] + 10, 'c', 30, color=GREEN, weight='bold')
f.text(t[4][0], t[4][1] + 50, 'b', 30, color=GREEN, weight='bold')
rows = [('a', 'old C=O oxygen, now O⁻', '1 bond + 3 lone pairs'),
        ('b', 'old ester oxygen', 'still bonded to ethyl'),
        ('c', 'came from hydroxide', 'neutral OH, keeps its H')]
for k, (lab, txt, sub) in enumerate(rows):
    y = 95 + 95*k
    f.text(540, y, lab, 30, color=GREEN, weight='bold', anchor='start')
    f.text(580, y, txt, 24, anchor='start')
    f.text(580, y + 34, sub, 21, color=BLUE, anchor='start')
f.text(270, 362, 'shaded C: four single bonds', 22, color=BLUE)
save('intermediate', f)

# 4. Step 2: the intermediate collapses and ethoxide leaves.
f = Figure('The intermediate collapses', 'O minus lone pair reforms C=O while the C-OEt bond pair leaves with ethoxide', 960, 360)
t = mol(f, 'CC([O-])(O)OCC', 220, 170, bond=(1, 2), flip=True, bond_px=70, font_px=28)
push_from_lp(f, t[2], 180, mid(t[1], t[2]), bend=-0.5, d=26)
curved_arrow(f, mid(t[1], t[4]), t[4], bend=-0.9, shorten=4)
arrow(f, 420, 510, 180)
mol(f, 'CC(=O)O', 620, 155, bond=(1, 2), bond_px=66, font_px=28)
plus(f, 725, 180)
mol(f, '[O-]CC', 865, 185, bond=(0, 1), deg=-30, bond_px=66, font_px=28)
f.text(220, 340, 'tetrahedral intermediate', 24, weight='bold')
f.text(620, 340, 'acetic acid', 24, weight='bold')
f.text(865, 340, 'ethoxide', 24, weight='bold')
save('elimination', f)

# 5. Aldehyde versus ester: only one has a group that can leave.
f = Figure('Same addition, different fate', 'Tetrahedral intermediates from acetaldehyde and from methyl acetate; only the ester one has a leaving group', 960, 370)
mol(f, 'CC([O-])(O)[H]', 240, 150, bond=(1, 2), bond_px=66, font_px=28, highlight_atoms=[4])
mol(f, 'CC([O-])(O)OC', 720, 150, bond=(1, 2), flip=True, bond_px=66, font_px=28, highlight_atoms=[4, 5])
f.text(240, 305, 'from an aldehyde', 25, weight='bold').text(240, 342, 'H⁻ would have to leave: no', 22, color=RED)
f.text(720, 305, 'from an ester', 25, weight='bold').text(720, 342, 'CH₃O⁻ can leave: yes', 22, color=GREEN)
f.line(480, 40, 480, 345, color='#c9d6dd', width=2)
save('aldehyde-vs-ester', f)

# 6. Step 3: proton transfer from the acid to ethoxide.
f = Figure('The proton transfer that drives hydrolysis', 'Ethoxide removes the acid proton, giving acetate and ethanol', 960, 420)
a = mol(f, 'CC(=O)O[H]', 170, 165, bond=(1, 2), bond_px=64, font_px=28)
e = mol(f, '[O-]CC', a[4][0] + 60, a[4][1] + 125, bond=(0, 1), deg=-30, anchor=0, bond_px=64, font_px=28)
push_from_lp(f, e[0], 100, a[4], bend=-0.25)
curved_arrow(f, mid(a[3], a[4]), a[3], bend=-0.9, shorten=4)
arrow(f, 470, 560, 200, 'K ≈ 10¹¹', lsize=26)
mol(f, 'CC(=O)[O-]', 650, 165, bond=(1, 2), bond_px=64, font_px=28)
plus(f, 738, 200)
mol(f, 'OCC', 860, 205, bond=(0, 1), deg=-30, bond_px=58, font_px=28)
f.text(150, 395, 'acetic acid, pKa 4.8', 23, weight='bold')
f.text(650, 395, 'acetate', 23, weight='bold')
f.text(850, 395, 'ethanol, pKa 16', 23, weight='bold')
save('proton-transfer', f)

# 7. Net base-promoted hydrolysis, then acid work-up.
f = Figure('Saponification and work-up', 'Ethyl acetate plus hydroxide gives acetate and ethanol; adding acid later gives acetic acid', 960, 470)
f.text(20, 40, 'in base', 24, color=BLUE, weight='bold', anchor='start')
mol(f, 'CC(=O)OCC', 150, 120, bond=(1, 2), bond_px=54, font_px=26)
plus(f, 285, 135)
f.text(345, 147, 'OH⁻', 30)
arrow(f, 400, 520, 140)
mol(f, 'CC(=O)[O-]', 620, 120, bond=(1, 2), bond_px=54, font_px=26)
plus(f, 705, 135)
mol(f, 'OCC', 835, 145, bond=(0, 1), deg=-30, bond_px=54, font_px=26)
f.line(30, 235, 930, 235, color='#c9d6dd', width=2)
f.text(20, 275, 'work-up', 24, color=BLUE, weight='bold', anchor='start')
mol(f, 'CC(=O)[O-]', 160, 350, bond=(1, 2), bond_px=54, font_px=26)
plus(f, 275, 365)
f.text(345, 377, 'H₃O⁺', 30)
arrow(f, 400, 520, 370)
mol(f, 'CC(=O)O', 620, 350, bond=(1, 2), bond_px=54, font_px=26)
plus(f, 720, 365)
f.text(800, 377, 'H₂O', 30)
f.text(480, 455, 'OH⁻ is used up, one per ester. Work-up only moves a proton.', 21, color=BLUE)
save('saponification', f)

# 8. Fischer esterification, step 1: protonate the carbonyl oxygen.
f = Figure('Acid activates the carbonyl', 'Hydronium protonates the carbonyl oxygen of acetic acid; the positive charge is shared onto carbon', 960, 360)
a = mol(f, 'CC(=O)O', 150, 195, bond=(1, 2), bond_px=64, font_px=28)
h = mol(f, '[H][OH2+]', a[2][0] + 125, a[2][1] - 10, bond=(0, 1), deg=0, anchor=0, bond_px=64, font_px=28)
push_from_lp(f, a[2], 30, h[0], bend=-0.4, d=26)
curved_arrow(f, mid(h[0], h[1]), h[1], bend=0.9, shorten=4)
arrow(f, 400, 480, 200)
mol(f, 'CC(=[OH+])O', 590, 195, bond=(1, 2), bond_px=64, font_px=28)
resonance(f, 690, 760, 200)
mol(f, 'C[C+](O)O', 860, 195, bond=(1, 2), bond_px=64, font_px=28)
f.text(200, 342, 'acetic acid + H₃O⁺', 24, weight='bold')
f.text(725, 342, 'protonated acid: C more δ+', 24, weight='bold')
save('fischer-activate', f)

# 9. Fischer step 2: neutral ethanol adds.
f = Figure('Ethanol adds to the activated carbonyl', 'Ethanol oxygen lone pair attacks carbon while the C=O pi pair moves to the protonated oxygen', 960, 420)
p = mol(f, 'CC(=[OH+])O', 190, 165, bond=(1, 2), bond_px=64, font_px=28)
e = mol(f, 'OCC', p[1][0] - 60, p[1][1] + 125, bond=(0, 1), deg=-30, anchor=0, bond_px=64, font_px=28)
push_from_lp(f, e[0], 70, p[1], bend=-0.25)
curved_arrow(f, mid(p[1], p[2]), p[2], bend=0.9, shorten=4)
arrow(f, 420, 510, 200)
mol(f, 'CC(O)(O)[OH+]CC', 720, 185, bond=(1, 2), flip=True, bond_px=64, font_px=28, highlight_atoms=[4])
f.text(200, 405, 'protonated acid + ethanol', 24, weight='bold')
f.text(720, 405, 'ethanol O now has 3 bonds: O⁺', 24, weight='bold')
save('fischer-attack', f)

# 10. Fischer step 3: proton shuffle turns OH into a leaving group.
f = Figure('A proton shuffle makes a leaving group', 'Water removes the proton from the ethanol oxygen; hydronium then protonates one OH to give OH2 plus', 960, 400)
mol(f, 'CC(O)(O)[OH+]CC', 140, 175, bond=(1, 2), flip=True, bond_px=42, font_px=26, highlight_atoms=[4])
arrow(f, 275, 345, 185, '− H⁺')
mol(f, 'CC(O)(O)OCC', 480, 175, bond=(1, 2), flip=True, bond_px=42, font_px=26)
arrow(f, 615, 685, 185, '+ H⁺')
mol(f, 'CC(O)([OH2+])OCC', 820, 175, bond=(1, 2), flip=True, bond_px=42, font_px=26, highlight_atoms=[3])
f.text(145, 320, 'O⁺ on ethanol O', 22, weight='bold')
f.text(480, 320, 'all neutral', 22, weight='bold')
f.text(820, 320, 'OH₂⁺: water can leave', 22, weight='bold')
f.text(480, 372, 'H₂O takes the first proton, H₃O⁺ gives the second', 21, color=BLUE)
save('fischer-shuffle', f)

# 11. Fischer step 4: water leaves, then deprotonation gives the ester.
f = Figure('Water leaves and the ester forms', 'An OH lone pair reforms C=O and pushes out water; water then removes the last proton', 960, 360)
t = mol(f, 'CC(O)([OH2+])OCC', 175, 175, bond=(1, 2), flip=True, bond_px=68, font_px=28)
push_from_lp(f, t[2], 165, mid(t[1], t[2]), bend=-0.5, d=26)
curved_arrow(f, mid(t[1], t[3]), t[3], bend=0.9, shorten=4)
arrow(f, 350, 415, 185, '− H₂O')
mol(f, 'CC(=[OH+])OCC', 525, 165, bond=(1, 2), bond_px=48, font_px=26)
arrow(f, 640, 705, 185, '− H⁺')
mol(f, 'CC(=O)OCC', 820, 165, bond=(1, 2), bond_px=48, font_px=26)
f.text(175, 340, 'OH₂⁺ leaves as H₂O', 23, weight='bold')
f.text(525, 340, 'protonated ester', 23, weight='bold')
f.text(820, 340, 'ethyl acetate', 23, weight='bold')
save('fischer-collapse', f)

# 12. The overall Fischer equilibrium with the ethanol oxygen labelled.
f = Figure('Fischer esterification is an equilibrium', 'Acetic acid plus ethanol gives ethyl acetate plus water; the ester single-bond oxygen comes from ethanol', 960, 360)
mol(f, 'CC(=O)O', 100, 150, bond=(1, 2), bond_px=60, font_px=28)
plus(f, 200, 160)
mol(f, '[18OH]CC', 290, 175, bond=(0, 1), deg=-30, bond_px=60, font_px=28, highlight_atoms=[0])
equilibrium(f, 390, 510, 175, 'H₃O⁺ (cat.)')
mol(f, 'CC(=O)[18O]CC', 650, 150, bond=(1, 2), bond_px=60, font_px=28, highlight_atoms=[3])
plus(f, 790, 160)
f.text(860, 187, 'H₂O', 32)
f.text(480, 300, '¹⁸O from ethanol ends up in the ester, not in the water', 23, color=BLUE)
f.text(480, 338, 'K ≈ 4: use excess alcohol or remove water', 23, color=BLUE)
save('fischer-overall', f)

# 13. Amide resonance.
f = Figure('Amide resonance', 'Acetamide and its charge-separated contributor with a C=N double bond', 960, 350)
a = mol(f, 'CC(N)=O', 250, 160, bond=(1, 3), bond_px=74, font_px=30)
push_from_lp(f, a[2], -120, mid(a[1], a[2]), bend=-0.7, d=28)
curved_arrow(f, mid(a[1], a[3]), a[3], bend=0.9, shorten=4)
resonance(f, 420, 540, 175)
mol(f, 'CC([O-])=[NH2+]', 700, 160, bond=(1, 2), bond_px=74, font_px=30)
f.text(250, 325, 'neutral contributor', 24, weight='bold')
f.text(700, 325, 'charge-separated: C=N, O⁻, N⁺', 24, weight='bold')
save('amide-resonance', f)

# 14. Amide versus amine nitrogen.
f = Figure('Amide N and amine N', 'In ethylamine the N lone pair is free; in acetamide it is shared with the carbonyl', 960, 330)
m1 = mol(f, 'CCN', 240, 140, bond=(1, 2), deg=-30, bond_px=74, font_px=30, highlight_atoms=[2])
lone_pair(f, m1[2], 90, distance=28)
m2 = mol(f, 'CC(N)=O', 720, 130, bond=(1, 3), bond_px=74, font_px=30, highlight_atoms=[2])
lone_pair(f, m2[2], -100, distance=30)
f.text(240, 265, 'ethylamine: lone pair free', 24, weight='bold').text(240, 302, 'a base (pKaH 10.7)', 22, color=BLUE)
f.text(720, 265, 'acetamide: lone pair shared', 24, weight='bold').text(720, 302, 'barely basic (pKaH ≈ 0)', 22, color=BLUE)
f.line(480, 30, 480, 310, color='#c9d6dd', width=2)
save('amide-amine', f)

# 15. Acid plus amine: proton transfer, not amide.
f = Figure('Acid plus amine gives a salt', 'Acetic acid and methylamine transfer a proton to give methylammonium acetate, not an amide', 960, 400)
mol(f, 'CC(=O)O', 110, 150, bond=(1, 2), bond_px=56, font_px=26)
plus(f, 205, 160)
mol(f, 'CN', 270, 175, bond=(0, 1), deg=-30, bond_px=56, font_px=26)
arrow(f, 345, 445, 175, 'fast')
mol(f, 'CC(=O)[O-]', 545, 150, bond=(1, 2), bond_px=56, font_px=26)
plus(f, 640, 160)
mol(f, 'C[NH3+]', 715, 175, bond=(0, 1), deg=-30, bond_px=56, font_px=26)
f.text(850, 185, 'a salt', 26, weight='bold', color=GREEN)
f.text(110, 270, 'pKa 4.8', 22, color=BLUE).text(715, 270, 'pKa 10.6', 22, color=BLUE)
f.text(480, 325, 'K ≈ 10⁵·⁸  (10.6 − 4.8 = 5.8)', 24)
f.text(480, 368, 'No amide yet: CH₃CONHCH₃ needs activation or strong heating', 22, color=RED)
save('acid-amine', f)

# 16. Ester and amide links in biology.
f = Figure('Ester and amide links in biology', 'Acetylcholine, an ester, and the dipeptide glycylalanine with one amide (peptide) link', 960, 340)
mol(f, 'CC(=O)OCC[N+](C)(C)C', 250, 145, bond=(1, 2), bond_px=56, font_px=26, highlight_atoms=[1, 2, 3])
mol(f, 'NCC(=O)N[C@@H](C)C(=O)O', 720, 145, bond=(2, 3), bond_px=56, font_px=26, highlight_atoms=[2, 3, 4])
f.text(250, 285, 'acetylcholine', 24, weight='bold').text(250, 320, 'ester: an enzyme cuts >10⁴ per second', 21, color=BLUE)
f.text(720, 285, 'Gly-Ala dipeptide', 24, weight='bold').text(720, 320, 'amide (peptide) link: very slow in water', 21, color=BLUE)
f.line(480, 30, 480, 330, color='#c9d6dd', width=2)
save('bio-links', f)

# 17. Aspirin hydrolysis: why an old bottle smells of vinegar.
f = Figure('Aspirin hydrolysis', 'Aspirin plus water gives salicylic acid and acetic acid', 960, 330)
mol(f, 'CC(=O)Oc1ccccc1C(=O)O', 150, 145, bond=(1, 2), bond_px=44, font_px=22, highlight_atoms=[1, 2, 3])
plus(f, 318, 155)
f.text(370, 167, 'H₂O', 30)
arrow(f, 420, 495, 160)
mol(f, 'Oc1ccccc1C(=O)O', 610, 145, bond=(1, 2), deg=90, bond_px=44, font_px=22, highlight_atoms=[0])
plus(f, 740, 155)
mol(f, 'CC(=O)O', 860, 145, bond=(1, 2), bond_px=44, font_px=22, highlight_atoms=[1, 2, 3])
f.text(150, 285, 'aspirin', 24, weight='bold').text(150, 318, 'ester shaded', 21, color=BLUE)
f.text(610, 285, 'salicylic acid', 24, weight='bold')
f.text(860, 285, 'acetic acid', 24, weight='bold').text(860, 318, 'smells of vinegar', 21, color=BLUE)
save('aspirin', f)

print('Built 17 skeletal structure figures.')
