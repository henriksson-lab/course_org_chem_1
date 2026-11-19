"""Skeletal-structure figures for the Week 5 integration lecture (RDKit via molfig)."""
from contextlib import contextmanager
from pathlib import Path
from rdkit import Chem
from chemfig import Figure, BLUE, RED, INK, GREEN
import molfig
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'integration'/'figures'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

@contextmanager
def keep_h():
    """molfig parses SMILES with default settings, which strips explicit [H].
    Here we want one chosen H drawn (stereocentre H, acidic O-H, N-H), so parse
    with removeHs=False for the duration of the block."""
    orig = Chem.MolFromSmiles
    def parse(s, *a):
        p = Chem.SmilesParserParams(); p.removeHs = False
        return orig(s, p)
    molfig.Chem.MolFromSmiles = parse
    try:
        yield
    finally:
        molfig.Chem.MolFromSmiles = orig

def arrow(f, x1, x2, y, label='', sub=''):
    f.line(x1, y, x2, y, color=INK, width=3, arrow=True)
    if label: f.text((x1 + x2)/2, y - 18, label, 21, color=BLUE)
    if sub: f.text((x1 + x2)/2, y + 34, sub, 21, color=BLUE)

# Writing the H explicitly moves it in the neighbour order, so the @/@@ tags differ
# from C[C@H](O)C(=O)O; the CIP labels are asserted below.
LACTIC_H = 'C[C@]([H])(O)C(=O)O'       # (S)-lactic acid with the stereocentre H kept
LACTIC_H_R = 'C[C@@]([H])(O)C(=O)O'    # (R)-lactic acid

def cip(smiles):
    p = Chem.SmilesParserParams(); p.removeHs = False
    m = Chem.MolFromSmiles(smiles, p)
    return [lab for _, lab in Chem.FindMolChiralCenters(m, useLegacyImplementation=False)]

for smi, want in [(LACTIC_H, 'S'), (LACTIC_H_R, 'R'), ('C[C@H](O)C(=O)O', 'S'),
                  ('C[C@H]([NH3+])C(=O)[O-]', 'S'), ('C[C@H]([N+]([H])([H])[H])C(=O)[O-]', 'S')]:
    assert cip(smi) == [want], (smi, cip(smi), want)

# 1. Lactic acid read locally: groups and the stereocentre.
f = Figure('Lactic acid read group by group', '', 960, 400)
c = place(f, 'C[C@H](O)C(=O)O', 330, 170, bond_px=80, highlight_atoms=[1], font_px=30)
f.text(c[1][0], c[1][1] + 70, '*', 40, color=RED, weight='bold')
f.text(720, 90, 'secondary alcohol', 25, weight='bold', color=RED)
f.text(720, 122, 'OH carbon has two C neighbours', 21, color=BLUE)
f.text(720, 185, 'carboxylic acid', 25, weight='bold', color=RED)
f.text(720, 217, 'OH on a C=O carbon', 21, color=BLUE)
f.text(720, 280, 'stereocentre (*)', 25, weight='bold', color=RED)
f.text(720, 312, 'OH, COOH, CH₃ and H', 21, color=BLUE)
f.text(330, 365, 'lactic acid · C₃H₆O₃ · 90.08 g/mol', 24)
save('lactic-features', f)

# 2. CIP priorities at the stereocentre, with the COOH vs CH3 tie-break.
with keep_h():
    f = Figure('CIP priorities at the lactic acid stereocentre', '', 960, 400)
    c = place(f, LACTIC_H, 260, 180, bond_px=80, font_px=30,
              numbers={3: '1', 4: '2', 0: '3', 2: '4'})
    f.text(260, 360, '(S)-lactic acid with priorities 1–4', 23, weight='bold')
    f.rect(500, 60, 200, 210, fill='#e4f0f4').rect(730, 60, 200, 210, fill='#f3f3f3')
    f.text(600, 105, 'COOH carbon', 24, weight='bold').text(830, 105, 'CH₃ carbon', 24, weight='bold')
    f.text(600, 160, '(O, O, O)', 30, color=RED).text(830, 160, '(H, H, H)', 30)
    f.text(600, 205, 'C=O counts O twice', 19, color=BLUE).text(830, 205, 'three H', 19, color=BLUE)
    f.text(715, 245, 'first atoms tie (C vs C)', 19, color=BLUE)
    f.text(715, 320, 'OH  >  COOH  >  CH₃  >  H', 28, weight='bold')
    save('lactic-cip', f)

# 2b. Reading R/S: trace 1-2-3, and reverse when H points toward you.
def rs_arc(f, centre, r=46, color=BLUE):
    """Anticlockwise (on screen) arc from the left of the centre round to the top."""
    cx, cy = centre
    f.path(f'M {cx - r:.1f} {cy + 8:.1f} A {r} {r} 0 1 0 {cx - 8:.1f} {cy - r:.1f}', color=color, width=3, arrow=True)

with keep_h():
    f = Figure('Tracing 1 to 2 to 3 at the lactic acid stereocentre', '', 960, 420)
    c = place(f, LACTIC_H, 240, 170, bond_px=82, font_px=30, numbers={3: '1', 4: '2', 0: '3', 2: '4'})
    rs_arc(f, c[1])
    f.text(240, 335, 'H away (hashed)', 24, weight='bold')
    f.text(240, 372, 'anticlockwise → S', 24, color=BLUE)
    c = place(f, LACTIC_H_R, 720, 170, bond_px=82, font_px=30, numbers={3: '1', 4: '2', 0: '3', 2: '4'})
    rs_arc(f, c[1])
    f.line(480, 40, 480, 390, color='#c9d6dc', width=2)
    f.text(720, 335, 'H toward you (wedge)', 24, weight='bold')
    f.text(720, 372, 'looks anticlockwise → reverse → R', 24, color=RED)
    save('lactic-rs', f)

# 3. The two enantiomers as mirror images.
with keep_h():
    f = Figure('(S)- and (R)-lactic acid are mirror images', '', 960, 420)
    place(f, LACTIC_H, 240, 170, bond_px=78, font_px=30)
    place(f, LACTIC_H_R, 720, 170, bond_px=78, font_px=30, flip=True)
    f.line(480, 40, 480, 300, color=BLUE, width=2, dash='8 8')
    f.text(240, 320, '(S)-lactic acid', 25, weight='bold').text(720, 320, '(R)-lactic acid', 25, weight='bold')
    f.text(240, 358, 'your muscles · (+) rotation', 21, color=BLUE)
    f.text(720, 358, 'some bacteria · (−) rotation', 21, color=BLUE)
    f.text(480, 400, 'same connectivity, opposite arrangement: not superimposable', 21)
    save('lactic-enantiomers', f)

# 4. Deprotonation by hydroxide, with both curved arrows.
with keep_h():
    f = Figure('Hydroxide removes the carboxyl proton', '', 960, 400)
    c = place(f, 'C[C@H](O)C(=O)O[H]', 205, 140, bond_px=70, font_px=28, flip=True, highlight_atoms=[6])
    hx, hy = c[6]
    o = place(f, '[OH-]', hx + 80, hy + 95, font_px=30)
    ox, oy = o[0][0] + 4, o[0][1]
    lone_pair(f, (ox, oy), 90, distance=24)
    curved_arrow(f, (ox, oy - 32), (hx + 6, hy + 14), bend=-0.35, shorten=4)
    mx, my = (c[5][0] + hx)/2, (c[5][1] + hy)/2
    curved_arrow(f, (mx + 2, my + 8), (c[5][0] + 4, c[5][1] + 20), bend=-0.6, shorten=0)
    arrow(f, 500, 580, 170)
    place(f, 'C[C@H](O)C(=O)[O-]', 720, 140, bond_px=70, font_px=28, flip=True, highlight_atoms=[5])
    f.text(895, 185, '+  H₂O', 28)
    f.text(240, 345, 'lactic acid + OH⁻', 23, weight='bold')
    f.text(730, 345, 'lactate + water', 23, weight='bold')
    f.text(480, 385, 'only the O–H bond breaks: the stereocentre is not touched', 21, color=BLUE)
    save('deprotonation', f)

# 5. Oxidation of the secondary alcohol to the ketone of pyruvic acid.
f = Figure('Lactic acid is oxidised to pyruvic acid', '', 960, 380)
place(f, 'C[C@H](O)C(=O)O', 220, 160, bond_px=74, font_px=28, highlight_atoms=[1, 2])
arrow(f, 410, 550, 165, 'oxidation', '−2 H')
place(f, 'CC(=O)C(=O)O', 740, 160, bond_px=74, font_px=28, highlight_atoms=[1, 2])
f.text(220, 300, 'lactic acid', 25, weight='bold').text(740, 300, 'pyruvic acid', 25, weight='bold')
f.text(220, 340, 'C2 tetrahedral · stereocentre', 21, color=BLUE)
f.text(740, 340, 'C2 trigonal planar · no R/S', 21, color=BLUE)
save('pyruvate', f)

# 6. Esterification with the ethanol oxygen tracked.
f = Figure('Ethanol oxygen ends up between the acyl carbon and ethyl', '', 960, 400)
place(f, 'C[C@H](O)C(=O)O', 105, 150, bond_px=52, font_px=26, highlight_atoms=[3, 5])
f.text(222, 165, '+', 34)
e = place(f, 'CCO', 322, 150, bond_px=52, font_px=26, highlight_atoms=[2])
f.text(e[2][0], e[2][1] + 44, 'O*', 24, color=RED, weight='bold')
arrow(f, 400, 490, 160, 'H⁺, heat')
p = place(f, 'C[C@H](O)C(=O)OCC', 680, 150, bond_px=52, font_px=26, highlight_atoms=[5])
f.text(p[5][0], p[5][1] + 44, 'O*', 24, color=RED, weight='bold')
f.text(900, 165, '+ H₂O', 26)
f.text(205, 300, 'lactic acid + ethanol', 23, weight='bold')
f.text(680, 300, 'ethyl lactate', 23, weight='bold')
f.text(480, 360, 'the acid loses its OH; the alcohol keeps its O* and loses only its H', 21, color=BLUE)
save('ester-label', f)

# 7. Addition (aldehyde) versus acyl substitution (ester).
f = Figure('Aldehyde addition versus ester substitution', '', 960, 470)
f.text(20, 30, 'aldehyde + cyanide: addition', 22, anchor='start', weight='bold')
place(f, 'CC=O', 100, 105, bond_px=52, font_px=24, highlight_atoms=[1, 2])
arrow(f, 180, 260, 110, 'CN⁻')
place(f, 'CC([O-])C#N', 380, 105, bond_px=52, font_px=24, highlight_atoms=[1, 2])
arrow(f, 510, 590, 110, 'H⁺')
place(f, 'CC(O)C#N', 720, 105, bond_px=52, font_px=24, highlight_atoms=[1, 2])
f.text(720, 205, 'stays tetrahedral: H⁻ can’t leave', 20, color=RED)
f.text(20, 255, 'ester + hydroxide: substitution', 22, anchor='start', weight='bold')
place(f, 'CC(=O)OC', 100, 335, bond_px=52, font_px=24, highlight_atoms=[1, 2])
arrow(f, 180, 260, 340, 'OH⁻')
place(f, 'CC([O-])(O)OC', 380, 335, bond_px=52, font_px=24, highlight_atoms=[1, 2])
arrow(f, 510, 590, 340, '− CH₃O⁻')
place(f, 'CC(=O)O', 720, 335, bond_px=52, font_px=24, highlight_atoms=[1, 2])
f.text(720, 430, 'C=O comes back', 20, color=GREEN)
f.text(380, 450, 'both: Nu adds, C=O π pair moves onto O', 20, color=BLUE)
save('add-vs-sub', f)

# 8. Alanine zwitterion: two charged sites, net zero.
f = Figure('L-alanine at pH 7 is a zwitterion', '', 960, 360)
c = place(f, 'C[C@H]([NH3+])C(=O)[O-]', 330, 160, bond_px=80, font_px=30, highlight_atoms=[2, 5])
f.text(720, 110, 'N: 4 bonds, 0 lone pairs  →  +1', 24, color=BLUE)
f.text(720, 160, 'O: 1 bond, 3 lone pairs  →  −1', 24, color=RED)
f.text(720, 225, 'net charge 0', 28, weight='bold')
f.text(720, 262, 'but two charged sites', 22)
f.text(330, 330, 'L-alanine, C₃H₇NO₂', 24, weight='bold')
save('alanine', f)

# 9. Raising the pH removes an ammonium proton.
with keep_h():
    f = Figure('At high pH the ammonium group loses H+', '', 960, 360)
    c = place(f, 'C[C@H]([N+]([H])([H])[H])C(=O)[O-]', 200, 160, bond_px=62, font_px=26, highlight_atoms=[2, 3])
    arrow(f, 390, 560, 165, '+ OH⁻', '− H₂O')
    place(f, 'C[C@H](N)C(=O)[O-]', 740, 160, bond_px=62, font_px=26, highlight_atoms=[2])
    f.text(200, 300, 'zwitterion · net 0', 24, weight='bold').text(740, 300, 'anion · net −1', 24, weight='bold')
    f.text(200, 338, 'N⁺: four bonds', 21, color=BLUE).text(740, 338, 'N: three bonds + lone pair', 21, color=BLUE)
    save('alanine-ph', f)

# 10. Unknown U1: ethyl acetate with its 1H NMR stick spectrum.
f = Figure('Ethyl acetate assigned against its 1H NMR spectrum', '', 960, 420)
c = place(f, 'CC(=O)OCC', 220, 150, bond_px=74, font_px=28, highlight_atoms=[1, 2])
for idx, lab, col in [(0, 'a', GREEN), (4, 'b', RED), (5, 'c', BLUE)]:
    f.text(c[idx][0], c[idx][1] + (42 if idx != 4 else -26), lab, 26, color=col, weight='bold')
f.text(220, 300, 'CH₃COOCH₂CH₃ · C₄H₈O₂', 23, weight='bold')
f.text(220, 338, 'DBE = (2·4 + 2 − 8)/2 = 1: the C=O', 21, color=BLUE)
x0, x1, base = 470, 930, 300
ppm = lambda d: x1 - (d/5.0)*(x1 - x0)
f.line(x0, base, x1, base, color=INK, width=2)
for d in range(0, 6):
    f.line(ppm(d), base, ppm(d), base + 8, color=INK, width=2).text(ppm(d), base + 30, str(d), 18)
f.text((x0 + x1)/2, base + 58, 'δ / ppm', 20)
def multiplet(d, heights, col, lab):
    n = len(heights); sp = 7
    for i, h in enumerate(heights):
        x = ppm(d) + (i - (n - 1)/2)*sp
        f.line(x, base, x, base - h, color=col, width=2.5)
    f.text(ppm(d), base - max(heights) - 14, lab, 19, color=col, weight='bold')
multiplet(4.12, [25, 75, 75, 25], RED, 'b · 4.12 q 2H')
multiplet(2.05, [150], GREEN, 'a · 2.05 s 3H')
multiplet(1.25, [55, 110, 55], BLUE, 'c · 1.25 t 3H')
save('u1-ester', f)

# 11. The hydroxyketone alternative predicts different signals.
f = Figure('Two C4H8O2 candidates predict different spectra', '', 960, 420)
place(f, 'CC(=O)OCC', 240, 140, bond_px=70, font_px=28, highlight_atoms=[4, 5])
place(f, 'CC(=O)C(O)C', 720, 140, bond_px=70, font_px=28, highlight_atoms=[3, 4, 5])
f.line(480, 40, 480, 400, color='#c9d6dc', width=2)
f.text(240, 270, 'ethyl acetate', 25, weight='bold').text(720, 270, '3-hydroxybutan-2-one', 25, weight='bold')
for i, (l, r) in enumerate([('O–CH₂: 2H quartet ✓', 'O–CH: 1H quartet ✗'),
                             ('CH₃ next to CH₂: triplet ✓', 'CH₃ next to CH: doublet ✗'),
                             ('no O–H ✓', 'O–H: broad IR band ✗')]):
    f.text(240, 310 + 34*i, l, 21, color=GREEN).text(720, 310 + 34*i, r, 21, color=RED)
save('alternative', f)

print('Built 12 skeletal structure figures.')
