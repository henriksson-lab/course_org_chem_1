"""Skeletal-structure figures for the Day 18 amines and aromatics lecture (RDKit via molfig)."""
import math
from pathlib import Path
from rdkit import Chem
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'amines-aromatics'/'figures'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

def neighbours(smiles, idx, explicit_h=False):
    mol = Chem.MolFromSmiles(smiles)
    if explicit_h:
        mol = Chem.AddHs(mol)
    return [n.GetIdx() for n in mol.GetAtomWithIdx(idx).GetNeighbors()]

def pair_out(fig, smiles, xy, idx, distance=26, turn=0, explicit_h=False):
    """Lone pair on atom idx, pointing away from its bonds (turn rotates it, degrees)."""
    ax, ay = xy[idx]
    sx = sum(xy[n][0] - ax for n in neighbours(smiles, idx, explicit_h))
    sy = sum(xy[n][1] - ay for n in neighbours(smiles, idx, explicit_h))
    angle = math.degrees(math.atan2(sy, -sx)) + turn
    lone_pair(fig, (ax, ay), angle, distance=distance)
    return angle

def mid(a, b):
    return ((a[0] + b[0])/2, (a[1] + b[1])/2)

def double_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE]

def ring_centre(xy, atoms):
    return (sum(xy[i][0] for i in atoms)/len(atoms), sum(xy[i][1] for i in atoms)/len(atoms))

def count_label(fig, xy, ring, a, b, label, push=34, size=24, color=GREEN):
    """Put a label just inside the ring next to bond a-b."""
    cx, cy = ring_centre(xy, ring)
    mx, my = mid(xy[a], xy[b])
    d = math.hypot(mx - cx, my - cy) or 1
    fig.text(mx - (mx - cx)/d*push, my - (my - cy)/d*push + size*0.35, label, size, color=color, weight='bold')

# 1. Amine classes: count carbons on N.
f = Figure('Primary, secondary and tertiary amines', 'Methylamine, dimethylamine and trimethylamine with nitrogen shaded and its lone pair drawn', 960, 400)
for x, smi, name, sub in [(160, 'CN', 'methylamine', 'primary · 1 C on N'),
                          (480, 'CNC', 'dimethylamine', 'secondary · 2 C on N'),
                          (800, 'CN(C)C', 'trimethylamine', 'tertiary · 3 C on N')]:
    n = 1
    xy = place(f, smi, x, 165, bond_px=70, highlight_atoms=[n])
    if smi == 'CN':
        lone_pair(f, xy[n], 90, distance=36)
    else:
        lone_pair(f, xy[n], -90, distance=30)
    f.text(x, 300, name, 25, weight='bold').text(x, 340, sub, 22, color=BLUE)
save('amine-classes', f)

# 2. Protonation of methylamine by hydronium, with both curved arrows.
f = Figure('Methylamine takes a proton from hydronium', 'Curved arrows: nitrogen lone pair to H, and the O-H bond pair back to oxygen, giving methylammonium and water', 960, 360)
xa = place(f, 'CN', 110, 170, bond_px=70)
lp = pair_out(f, 'CN', xa, 1, distance=30, turn=60)
f.text(225, 182, '+', 34)
xh = place(f, '[OH3+]', 350, 170, bond_px=80, explicit_h=True, font_px=26)
# the H closest to the amine receives the arrow
hs = [i for i in xh if i != 0]
h = min(hs, key=lambda i: xh[i][0])
lpx = xa[1][0] + 30*math.cos(math.radians(lp)); lpy = xa[1][1] - 30*math.sin(math.radians(lp))
curved_arrow(f, (lpx, lpy), xh[h], bend=-0.45, shorten=10)
m = mid(xh[h], xh[0])
curved_arrow(f, (m[0], m[1] + 6), (xh[0][0] - 4, xh[0][1] + 16), bend=0.7, shorten=2)
f.line(470, 170, 570, 170, arrow=True, color=BLUE)
place(f, 'C[NH3+]', 680, 170, bond_px=70)
f.text(785, 182, '+', 34)
f.text(860, 182, 'H₂O', 32, color=RED)
f.text(140, 300, 'base: 3 bonds + 1 pair', 22, color=BLUE)
f.text(345, 300, 'acid: H₃O⁺', 22, color=BLUE)
f.text(680, 300, 'N: 4 bonds, charge +1', 22, color=BLUE)
f.text(860, 300, 'water', 22, color=BLUE)
save('protonation', f)

# 3. Tertiary amine, its conjugate acid, quaternary ammonium.
save_items = [('CN(C)C', 'trimethylamine', '3 C · 0 N–H · neutral'),
              ('C[NH+](C)C', 'trimethylammonium', '3 C · 1 N–H · +1'),
              ('C[N+](C)(C)C', 'tetramethylammonium', '4 C · 0 N–H · +1')]
f = Figure('Tertiary amine, its conjugate acid and a quaternary ammonium ion', 'Trimethylamine, trimethylammonium and tetramethylammonium with carbon count, N-H count and charge', 960, 400)
for i, (smi, name, sub) in enumerate(save_items):
    x = 160 + 320*i
    xy = place(f, smi, x, 165, bond_px=66, highlight_atoms=[1])
    if i == 0:
        pair_out(f, smi, xy, 1, distance=30)
    f.text(x, 300, name, 23, weight='bold').text(x, 340, sub, 22, color=BLUE)
save('tertiary-quaternary', f)

# 4. Conjugate-acid pKa bars with the neutral bases drawn underneath.
f = Figure('Conjugate-acid pKa values of five nitrogen bases', 'Bar chart: ethylammonium 10.7, benzylammonium 9.3, ammonium 9.2, pyridinium 5.2, anilinium 4.6, with each base drawn below its bar', 960, 470)
bars = [('CCN', 'ethylamine', 10.7), ('NCc1ccccc1', 'benzylamine', 9.3), ('N', 'ammonia', 9.2),
        ('c1ccncc1', 'pyridine', 5.2), ('Nc1ccccc1', 'aniline', 4.6)]
base_y, scale = 260, 17
f.line(60, base_y, 920, base_y, width=2)

for i, (smi, name, pka) in enumerate(bars):
    x = 110 + 185*i
    h = pka*scale
    f.rect(x - 45, base_y - h, 90, h, fill=BLUE if pka > 7 else GREEN)
    f.text(x, base_y - h - 12, f'{pka}', 24, weight='bold')
    place(f, smi, x, 335, bond_px=34, font_px=20, highlight_atoms=[a for a, s in enumerate(Chem.MolFromSmiles(smi).GetAtoms()) if s.GetSymbol() == 'N'])
    f.text(x, 445, name, 22, weight='bold')
f.text(920, 40, 'conjugate-acid pKa · higher bar = stronger base', 22, anchor='end', color=BLUE)
save('pka-bars', f)

# 5. Proton transfer between two nitrogen bases.
f = Figure('Anilinium gives its proton to ethylamine', 'Anilinium plus ethylamine in equilibrium with aniline plus ethylammonium; pKa 4.6 and 10.7 give K about 10^6', 960, 440)
place(f, '[NH3+]c1ccccc1', 150, 85, bond_px=44, font_px=24, highlight_atoms=[0])
f.text(300, 95, '+', 32)
place(f, 'CCN', 410, 85, bond_px=50, font_px=24, highlight_atoms=[2])
f.text(150, 172, 'acid, pKa 4.6', 21, color=RED).text(410, 172, 'base', 21, color=BLUE)
f.line(270, 185, 270, 240, arrow=True, color=BLUE, width=3)
f.line(295, 240, 295, 185, arrow=True, color=BLUE, width=1.5)
place(f, 'Nc1ccccc1', 150, 300, bond_px=44, font_px=24, highlight_atoms=[0])
f.text(300, 310, '+', 32)
place(f, 'CC[NH3+]', 410, 300, bond_px=50, font_px=24, highlight_atoms=[2])
f.text(150, 387, 'base', 21, color=BLUE).text(410, 387, 'acid, pKa 10.7', 21, color=RED)
f.text(590, 150, 'K ≈ 10^(10.7 − 4.6)', 27, anchor='start', weight='bold')
f.text(590, 195, '= 10^6.1 ≈ 1 × 10⁶', 27, anchor='start', weight='bold')
f.text(590, 260, 'products win: the proton', 23, anchor='start', color=BLUE)
f.text(590, 295, 'ends up on the stronger base', 23, anchor='start', color=BLUE)
save('proton-transfer', f)

# 6. Amide resonance: the N pair is shared with C=O.
f = Figure('The amide nitrogen pair is shared with the carbonyl', 'Ethanamide drawn neutral and as the charge-separated contributor, with curved arrows from the N lone pair and the C=O pi bond', 960, 400)
smi = 'CC(N)=O'
xy = place(f, smi, 260, 170, bond_px=78, highlight_atoms=[2])
ang = pair_out(f, smi, xy, 2, distance=30, turn=-40)
lp = (xy[2][0] + 30*math.cos(math.radians(ang)), xy[2][1] - 30*math.sin(math.radians(ang)))
curved_arrow(f, lp, mid(xy[1], xy[2]), bend=0.6, shorten=6)
m = mid(xy[1], xy[3])
curved_arrow(f, (m[0] + 8, m[1] + 10), (xy[3][0] + 18, xy[3][1] + 10), bend=-0.6, shorten=2)
f.line(430, 170, 530, 170, color=BLUE, arrow=True, both=True)
place(f, 'CC(=[NH2+])[O-]', 700, 170, bond_px=78, highlight_atoms=[2])
f.text(260, 315, 'neutral ethanamide', 24, weight='bold').text(260, 352, 'N: 3 bonds + 1 pair', 22, color=BLUE)
f.text(700, 315, 'charge-separated contributor', 24, weight='bold').text(700, 352, 'N +1 · O −1 · net 0', 22, color=BLUE)
save('amide-resonance', f)

# 7. A peptide-like fragment has two different nitrogens.
f = Figure('Two kinds of nitrogen in one molecule', 'H2N-CH2-C(=O)-NH-CH3 with the amine nitrogen and the amide nitrogen labelled', 960, 360)
smi = 'NCC(=O)NC'
xy = place(f, smi, 480, 160, bond_px=72, highlight_atoms=[0, 4])
f.text(xy[0][0] + 55, xy[0][1] - 5, 'amine N', 24, weight='bold', color=BLUE, anchor='start')
f.text(xy[0][0] + 55, xy[0][1] + 25, 'pair free · basic', 21, color=BLUE, anchor='start')
f.text(xy[4][0] - 55, xy[4][1] - 5, 'amide N', 24, weight='bold', color=RED, anchor='end')
f.text(xy[4][0] - 55, xy[4][1] + 25, 'pair shared with C=O', 21, color=RED, anchor='end')
f.text(480, 335, 'the peptide bond in every protein is this amide', 22)
save('peptide', f)

# 8. Aniline versus benzylamine.
f = Figure('Aniline and benzylamine', 'Aniline with N directly on the ring, and benzylamine where a CH2 separates N from the ring; conjugate-acid pKa 4.6 and 9.3', 960, 400)
xy = place(f, 'Nc1ccccc1', 260, 165, bond_px=66, highlight_atoms=[0])
f.text(260, 315, 'aniline', 25, weight='bold').text(260, 352, 'N on the ring: pair shared · pKa(BH⁺) 4.6', 21, color=BLUE)
xy = place(f, 'NCc1ccccc1', 700, 165, bond_px=66, highlight_atoms=[1])
f.text(700, 315, 'benzylamine', 25, weight='bold').text(700, 352, 'CH₂ breaks the link · pKa(BH⁺) 9.3', 21, color=BLUE)
save('aniline-benzylamine', f)

# 9. Benzene: count pi electrons.
f = Figure('Benzene has six pi electrons', 'Benzene in Kekule form with 2 pi electrons counted at each double bond, and the four aromaticity checks', 960, 400)
smi = 'c1ccccc1'
xy = place(f, smi, 240, 190, bond_px=95)
for a, b in double_bonds(smi):
    count_label(f, xy, range(6), a, b, '2', push=-34)
checks = ['① ring', '② every ring atom has a p orbital', '③ flat, so the p orbitals overlap', '④ 4n + 2 π electrons: 2, 6, 10, …']
for i, c in enumerate(checks):
    f.text(470, 115 + 52*i, c, 25, anchor='start')
f.text(470, 345, 'benzene: 2 + 2 + 2 = 6 = 4(1) + 2', 25, anchor='start', weight='bold', color=GREEN)
save('benzene', f)

# 10. A ring is not enough.
f = Figure('Aromatic, nonaromatic and antiaromatic rings', 'Benzene aromatic with 6 pi electrons; cyclohexa-1,3-diene and cyclohexane nonaromatic; cyclobutadiene antiaromatic with 4', 960, 400)
items = [('c1ccccc1', 'benzene', 'aromatic', '6 π, closed loop', GREEN, []),
         ('C1=CC=CCC1', 'cyclohexa-1,3-diene', 'nonaromatic', 'two CH₂ break the loop', INK, [4, 5]),
         ('C1CCCCC1', 'cyclohexane', 'nonaromatic', 'no π electrons', INK, []),
         ('C1=CC=C1', 'cyclobutadiene', 'antiaromatic', '4 π: wrong count', RED, [])]
for i, (smi, name, kind, sub, col, hl) in enumerate(items):
    x = 120 + 240*i
    place(f, smi, x, 150, bond_px=62, highlight_atoms=hl)
    f.text(x, 275, name, 22, weight='bold').text(x, 312, kind, 23, weight='bold', color=col).text(x, 346, sub, 19, color=BLUE)
save('ring-types', f)

# 11. Pyridine and pyridinium.
f = Figure('Pyridine keeps its lone pair outside the aromatic sextet', 'Pyridine with its in-plane N lone pair taking a proton to give aromatic pyridinium, pKa 5.2', 960, 400)
smi = 'c1ccncc1'
xy = place(f, smi, 230, 165, bond_px=68, highlight_atoms=[3])
pair_out(f, smi, xy, 3, distance=32)
for a, b in double_bonds(smi):
    count_label(f, xy, range(6), a, b, '2', push=-30, size=22)
f.line(400, 165, 540, 165, arrow=True, color=BLUE)
f.text(470, 145, '+ H⁺', 24)
xy2 = place(f, 'c1cc[nH+]cc1', 720, 165, bond_px=68, highlight_atoms=[3])
f.text(230, 340, 'pyridine: 6 π + a separate pair', 23, weight='bold')
f.text(720, 340, 'pyridinium: still 6 π, aromatic', 23, weight='bold')
f.text(470, 380, 'the pair points out of the ring, in the plane: it is not part of the sextet', 21, color=BLUE)
save('pyridine', f)

# 12. Pyrrole: the N pair is part of the sextet.
f = Figure('Pyrrole uses its nitrogen pair in the sextet', 'Pyrrole with 2 pi electrons from each C=C and 2 from the N lone pair, total 6', 960, 400)
smi = 'c1cc[nH]c1'
xy = place(f, smi, 270, 175, bond_px=90, highlight_atoms=[3])
for a, b in double_bonds(smi):
    count_label(f, xy, range(5), a, b, '2', push=-34)
# N lone pair drawn inside the ring, where it joins the pi system
cx, cy = ring_centre(xy, range(5))
nx, ny = xy[3]
ang = math.degrees(math.atan2(-(cy - ny), cx - nx))
lone_pair(f, (nx, ny), ang, distance=30, color=GREEN)
f.text(nx + (cx - nx)*0.75, ny + (cy - ny)*0.75 + 8, '2', 24, color=GREEN, weight='bold')
f.text(500, 120, '2 + 2 + 2 = 6 π electrons', 27, anchor='start', weight='bold', color=GREEN)
f.text(500, 175, 'the N pair is in the p orbital', 23, anchor='start')
f.text(500, 215, 'protonating N breaks the sextet', 23, anchor='start')
f.text(500, 280, 'pyrrole N is barely basic', 24, anchor='start', weight='bold', color=RED)
f.text(500, 318, '(pyridine N: pKa(BH⁺) 5.2)', 22, anchor='start', color=BLUE)
save('pyrrole', f)

# 13. Imidazole has one nitrogen of each kind.
f = Figure('Imidazole has one pyridine-like and one pyrrole-like nitrogen', 'Imidazole with labelled nitrogens, protonated at the pyridine-like nitrogen to give imidazolium, pKa about 7', 960, 400)
smi = 'c1cnc[nH]1'
xy = place(f, smi, 300, 170, bond_px=78, highlight_atoms=[2, 4])
pair_out(f, smi, xy, 2, distance=30)
f.text(xy[2][0] - 45, xy[2][1] + 8, 'pyridine-like', 21, color=BLUE, weight='bold', anchor='end')
f.text(xy[4][0] - 35, xy[4][1] + 8, 'pyrrole-like', 21, color=GREEN, weight='bold', anchor='end')
f.line(470, 170, 580, 170, arrow=True, color=BLUE)
f.text(525, 150, '+ H⁺', 24)
place(f, 'c1c[nH+]c[nH]1', 740, 170, bond_px=78, highlight_atoms=[2, 4])
f.text(300, 330, 'imidazole', 25, weight='bold')
f.text(740, 330, 'imidazolium, pKa ≈ 7', 25, weight='bold')
f.text(480, 375, 'histidine side chain: half protonated near pH 6–7', 22, color=BLUE)
save('imidazole', f)

# 14. Glycine across the pH scale.
f = Figure('Glycine at low, middle and high pH', 'Three glycine forms on a pH axis with pKa 2.3 and 9.6: cation net +1, zwitterion net 0, anion net -1', 960, 430)
forms = [('[NH3+]CC(=O)O', 'pH 1', 'net +1'), ('[NH3+]CC(=O)[O-]', 'pH 6 · zwitterion', 'net 0'), ('NCC(=O)[O-]', 'pH 12', 'net −1')]
for i, (smi, ph, net) in enumerate(forms):
    x = 160 + 320*i
    place(f, smi, x, 125, bond_px=52, font_px=24)
    f.text(x, 245, ph, 23, weight='bold').text(x, 280, net, 23, color=BLUE)
y = 340
f.line(60, y, 900, y, color=INK, width=3, arrow=True)
f.text(915, y + 8, 'pH', 22, anchor='start')
for ph in range(0, 15, 2):
    px = 60 + ph*58
    f.line(px, y - 6, px, y + 6, width=2).text(px, y + 30, ph, 18)
for pka, label in [(2.3, 'pKa 2.3 · COOH'), (9.6, 'pKa 9.6 · NH₃⁺')]:
    px = 60 + pka*58
    f.line(px, y - 22, px, y + 22, color=RED, width=4).text(px, y + 62, label, 20, color=RED)
save('glycine-ph', f)

# 15. Protonated fraction of a simple amine against pH.
f = Figure('Fraction of an amine protonated, pKa 10', 'Curve of fraction BH+ against pH for pKa 10, falling from 0.91 at pH 9 to 0.5 at pH 10 to 0.09 at pH 11', 960, 400)
x0, x1, y0, y1 = 120, 900, 330, 50
X = lambda ph: x0 + (ph - 7)/(13 - 7)*(x1 - x0)
Y = lambda fr: y0 - fr*(y0 - y1)
f.line(x0, y0, x1, y0, width=2).line(x0, y0, x0, y1, width=2)
for ph in range(7, 14):
    f.line(X(ph), y0, X(ph), y0 + 6, width=2).text(X(ph), y0 + 30, ph, 18)
for fr in (0, 0.5, 1):
    f.text(x0 - 12, Y(fr) + 7, fr, 18, anchor='end')
f.text(510, 392, 'pH', 22).text(x0 - 12, y1 - 18, 'fraction BH⁺', 20, anchor='start', color=BLUE)
pts = ' '.join(f'{"M" if k == 0 else "L"} {X(7 + k*0.05):.1f} {Y(1/(1 + 10**(7 + k*0.05 - 10))):.1f}' for k in range(121))
f.path(pts, color=BLUE, width=4)
for ph, fr, lab in [(9, 0.909, '0.91'), (10, 0.5, '0.50'), (11, 0.0909, '0.09')]:
    f.circle(X(ph), Y(fr), 8, color=RED, fill=RED, width=1)
    f.text(X(ph) + 16, Y(fr) - 12, f'pH {ph}: {lab}', 21, anchor='start', color=RED)
f.text(880, 90, '[B]/[BH⁺] = 10^(pH − pKa)', 23, anchor='end', weight='bold')
save('fraction-curve', f)

print('Built 15 skeletal structure figures for amines-aromatics.')
