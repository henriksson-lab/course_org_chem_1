"""Skeletal-structure figures for the Day 13 alcohols and redox lecture (RDKit via molfig)."""
from pathlib import Path
from contextlib import contextmanager
from rdkit import Chem
import molfig
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place, curved_arrow, lone_pair

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'alcohols-redox'/'figures'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

# molfig parses SMILES with RDKit's default, which strips written [H] atoms.
# To show only the H atoms that matter (on the carbinol carbon, on water),
# parse with removeHs=False for that one call. Works around molfig without editing it.
class _KeepH:
    def __getattr__(self, name):
        return getattr(Chem, name)
    @staticmethod
    def MolFromSmiles(smiles):
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        return Chem.MolFromSmiles(smiles, ps)

def place_h(fig, smiles, x, y, **kw):
    """place() that keeps H atoms written as [H] in the SMILES."""
    saved = molfig.Chem
    molfig.Chem = _KeepH()
    try:
        return place(fig, smiles, x, y, **kw)
    finally:
        molfig.Chem = saved

def rotation_to_point(smiles, atom, anchor, target_deg=0, keep_h=True):
    """Find the rotate= angle that points the anchor->atom bond along target_deg (0 = right)."""
    import math
    probe = Figure('probe', '', 960, 400)
    c = (place_h if keep_h else place)(probe, smiles, 480, 200)
    (x1, y1), (x2, y2) = c[anchor], c[atom]
    current = math.degrees(math.atan2(-(y2 - y1), x2 - x1))
    return target_deg - current

def arrow(f, x1, x2, y, above='', below='', color=BLUE, size=21):
    f.line(x1, y, x2, y, color=color, width=3, arrow=True)
    if above: f.text((x1 + x2)/2, y - 14, above, size, color=color)
    if below: f.text((x1 + x2)/2, y + 30, below, size, color=color)

ROT_O_RIGHT = rotation_to_point('CC([H])(C)[O-]', 4, 1, 0)

# 1. Hook: how the liver oxidises ethanol.
f = Figure('Ethanol metabolism: ethanol to acetaldehyde to acetate',
           'Alcohol dehydrogenase oxidises ethanol to ethanal; aldehyde dehydrogenase oxidises ethanal to acetate', 960, 380)
for x, smi, hl, name, sub in [(115, 'CCO', [1], 'ethanol', 'what you drink'),
                              (480, 'CC=O', [1], 'ethanal (acetaldehyde)', 'toxic: flushing, nausea'),
                              (850, 'CC(=O)[O-]', [1], 'acetate', 'used for energy')]:
    place(f, smi, x, 150, highlight_atoms=hl, bond_px=66)
    f.text(x, 275, name, 24, weight='bold').text(x, 310, sub, 21, color=BLUE)
arrow(f, 225, 375, 150, 'ADH', 'NAD⁺ → NADH', size=19)
arrow(f, 590, 740, 150, 'ALDH', 'NAD⁺ → NADH', size=19)
f.text(480, 360, 'disulfiram (Antabuse) blocks ALDH, so ethanal builds up', 21, color=RED)
save('metabolism', f)

# 2. Classify the carbinol carbon, with its H atoms drawn.
f = Figure('Primary, secondary and tertiary alcohols',
           'Ethanol, propan-2-ol and 2-methylpropan-2-ol with the OH-bearing carbon shaded and its H atoms drawn', 960, 420)
for x, smi, hl, name, cls, sub in [
        (160, 'CC([H])([H])O', [1], 'ethanol', 'primary (1°)', '1 C neighbour · 2 H'),
        (480, 'CC([H])(C)O', [1], 'propan-2-ol', 'secondary (2°)', '2 C neighbours · 1 H'),
        (800, 'CC(C)(C)O', [1], '2-methylpropan-2-ol', 'tertiary (3°)', '3 C neighbours · no H')]:
    place_h(f, smi, x, 150, highlight_atoms=hl, bond_px=66)
    f.text(x, 290, cls, 25, weight='bold').text(x, 325, name, 21).text(x, 360, sub, 21, color=BLUE)
save('classes', f)

# 3. Oxidation states of the changing carbon along ethanol -> ethanal -> ethanoic acid.
f = Figure('Oxidation state of the changing carbon',
           'Ethanol, ethanal and ethanoic acid with the changing carbon shaded: oxidation states minus 1, plus 1 and plus 3', 960, 440)
for x, smi, name, calc, val in [
        (150, 'CC([H])([H])O', 'ethanol', '+1 − 2', '−1'),
        (480, 'CC([H])=O', 'ethanal', '+2 − 1', '+1'),
        (810, 'CC(=O)O', 'ethanoic acid', '+2 + 1', '+3')]:
    place_h(f, smi, x, 140, highlight_atoms=[1], bond_px=66)
    f.text(x, 270, name, 24, weight='bold').text(x, 310, f'C–O and C–H: {calc}', 21)
    f.text(x, 360, val, 40, color=RED, weight='bold')
arrow(f, 270, 360, 140, '[O]')
arrow(f, 600, 690, 140, '[O]')
f.text(480, 420, 'each C–O bond +1 · each C–H −1 · each C–C 0', 21, color=BLUE)
save('ox-states', f)

# 4. Primary alcohol in the lab: stop at the aldehyde or go on to the acid.
f = Figure('Propan-1-ol oxidation in the lab',
           'PCC stops at propanal; aqueous dichromate in acid continues to propanoic acid', 960, 400)
place(f, 'CCCO', 120, 200, highlight_atoms=[2], bond_px=54)
place(f, 'CCC=O', 520, 110, highlight_atoms=[2], bond_px=58)
place(f, 'CCC(=O)O', 520, 300, highlight_atoms=[2], bond_px=58)
f.line(255, 175, 420, 115, color=BLUE, arrow=True).line(255, 225, 420, 290, color=BLUE, arrow=True)
f.text(325, 125, 'PCC, CH₂Cl₂', 21, color=BLUE)
f.text(345, 310, 'Na₂Cr₂O₇, H₂SO₄, H₂O', 21, color=BLUE)
f.text(120, 110, 'propan-1-ol (1°)', 23, weight='bold')
f.text(650, 115, 'propanal', 23, weight='bold', anchor='start')
f.text(650, 150, 'no water: stops here', 21, color=GREEN, anchor='start')
f.text(650, 305, 'propanoic acid', 23, weight='bold', anchor='start')
f.text(650, 340, 'in water: oxidised again', 21, color=RED, anchor='start')
save('primary', f)

# 5. Secondary alcohol: the carbonyl stays where the OH was.
f = Figure('Butan-2-ol oxidises to butan-2-one, not butanal',
           'Butan-2-ol with numbered carbons gives butan-2-one; butanal is crossed out because the C=O moved to C1', 960, 380)
place(f, 'CC(O)CC', 150, 140, numbers={0: 1, 1: 2, 3: 3, 4: 4}, highlight_atoms=[1], bond_px=66)
arrow(f, 290, 410, 150, 'oxidation')
place(f, 'CC(=O)CC', 560, 140, numbers={0: 1, 1: 2, 3: 3, 4: 4}, highlight_atoms=[1, 2], bond_px=66)
place(f, 'CCCC=O', 830, 140, bond_px=46, font_px=24)
f.line(760, 95, 900, 190, color=RED, width=4).line(900, 95, 760, 190, color=RED, width=4)
f.text(150, 270, 'butan-2-ol (2°)', 24, weight='bold')
f.text(560, 270, 'butan-2-one', 24, weight='bold').text(560, 305, 'C=O stays on C2', 21, color=GREEN)
f.text(830, 270, 'butanal', 24, weight='bold').text(830, 305, 'C=O moved to C1', 21, color=RED)
f.text(480, 360, 'the OH carbon keeps both of its carbon neighbours', 21, color=BLUE)
save('secondary', f)

# 6. Tertiary alcohol: no H on the carbinol carbon.
f = Figure('Why tertiary alcohols resist simple oxidation',
           'Propan-2-ol has one H on its OH carbon and oxidises to propanone; 2-methylpropan-2-ol has none', 960, 400)
place_h(f, 'CC([H])(C)O', 150, 150, highlight_atoms=[1, 2], bond_px=62)
arrow(f, 255, 345, 160, 'oxidation')
place(f, 'CC(C)=O', 440, 150, highlight_atoms=[1], bond_px=62)
f.text(150, 285, 'propan-2-ol', 24, weight='bold').text(150, 320, 'C–H to lose ✓', 21, color=GREEN)
f.text(440, 285, 'propanone', 24, weight='bold')
f.line(570, 40, 570, 340, color='#c9d6dd', width=2)
place(f, 'CC(C)(C)O', 790, 150, highlight_atoms=[1], bond_px=62)
f.text(790, 285, '2-methylpropan-2-ol', 24, weight='bold').text(790, 320, 'no C–H on this carbon ✗', 21, color=RED)
f.text(480, 380, 'C=O here would need a fifth bond to carbon', 21, color=BLUE)
save('tertiary', f)

# 7. Reduction runs the map backwards.
f = Figure('NaBH4 reduces aldehydes and ketones to alcohols',
           'Propanal gives propan-1-ol and propanone gives propan-2-ol after NaBH4 then water', 960, 440)
for y, smi1, smi2, n1, n2, c in [
        (95, 'CCC=O', 'CCCO', 'propanal', 'propan-1-ol (1°)', 'aldehyde → primary'),
        (290, 'CC(C)=O', 'CC(C)O', 'propanone', 'propan-2-ol (2°)', 'ketone → secondary')]:
    place(f, smi1, 170, y, highlight_atoms=[2] if 'CCC=O' == smi1 else [1], bond_px=58)
    arrow(f, 330, 560, y, '1. NaBH₄  2. H₂O')
    place(f, smi2, 720, y, highlight_atoms=[2] if smi1 == 'CCC=O' else [1], bond_px=58)
    f.text(170, y + 95, n1, 22, weight='bold').text(720, y + 95, n2, 22, weight='bold')
    f.text(445, y + 32, c, 20, color=GREEN)
save('reduction-map', f)

# 8. Step 1: hydride delivery to the carbonyl carbon.
f = Figure('Hydride delivery to propanone gives an alkoxide',
           'Curved arrows: hydride pair to the carbonyl carbon, C=O pi pair onto oxygen, giving an alkoxide', 960, 400)
c = place(f, 'CC(C)=O', 330, 200, bond_px=72, font_px=30)
hx, hy = c[1][0] - 190, c[1][1] + 20
f.text(hx, hy + 10, 'H:⁻', 34, color=BLUE, weight='bold')
curved_arrow(f, (hx + 20, hy - 10), c[1], bend=-0.35, color=RED)
mid = ((c[1][0] + c[3][0])/2, (c[1][1] + c[3][1])/2)
curved_arrow(f, mid, (c[3][0] + 18, c[3][1]), bend=-0.9, color=RED, shorten=6)
lone_pair(f, c[3], 160); lone_pair(f, c[3], 20)
arrow(f, 450, 560, 200)
p = place_h(f, 'CC([H])(C)[O-]', 740, 200, bond_px=72, font_px=30, highlight_atoms=[1, 2])
pass
f.text(330, 330, 'propanone', 24, weight='bold')
f.text(740, 330, 'alkoxide', 24, weight='bold').text(740, 365, 'new C–H · O carries −1', 21, color=BLUE)
f.text(170, 365, 'H:⁻ = delivered by BH₄⁻', 21, color=BLUE)
save('hydride', f)

# 9. Step 2: proton transfer from water. Water is drawn by hand so its H points at the alkoxide.
f = Figure('Water protonates the alkoxide',
           'The alkoxide oxygen takes a proton from water, giving propan-2-ol and hydroxide', 960, 380)
a = place_h(f, 'CC([H])(C)[O-]', 150, 170, bond_px=62, font_px=28, rotate=ROT_O_RIGHT)
ox, oy = a[4]
for ang in (95, 265, 330): lone_pair(f, a[4], ang)
hx = ox + 95
f.text(hx, oy + 10, 'H', 28).line(hx + 13, oy, hx + 47, oy, width=2.5).text(hx + 72, oy + 10, 'OH', 28, color=RED)
curved_arrow(f, (ox + 24, oy + 14), (hx + 2, oy + 12), bend=-0.6, color=RED, shorten=6)
curved_arrow(f, (hx + 30, oy - 4), (hx + 66, oy - 22), bend=0.7, color=RED, shorten=2)
arrow(f, 510, 590, oy)
place_h(f, 'CC([H])(C)O', 700, 170, bond_px=62, font_px=28, rotate=ROT_O_RIGHT)
f.text(885, oy + 10, '+ OH⁻', 30)
f.text(150, 300, 'alkoxide', 23, weight='bold').text(hx + 40, 300, 'water', 23, weight='bold')
f.text(700, 300, 'propan-2-ol', 23, weight='bold')
f.text(480, 355, 'charge: −1 on the left, −1 on the right', 21, color=BLUE)
save('workup', f)

# 10. A planar ketone attacked from either face gives both enantiomers.
f = Figure('Reducing butan-2-one gives racemic butan-2-ol',
           'Hydride adds to either face of planar butan-2-one, giving (R)- and (S)-butan-2-ol in equal amounts', 960, 400)
place(f, 'CC(=O)CC', 170, 170, highlight_atoms=[1], bond_px=64)
f.text(170, 300, 'butan-2-one', 24, weight='bold').text(170, 335, 'planar C=O, two faces', 21, color=BLUE)
f.line(320, 170, 440, 110, color=BLUE, arrow=True).line(320, 170, 440, 250, color=BLUE, arrow=True)
f.text(400, 177, 'NaBH₄', 20, color=BLUE)
place(f, 'C[C@@H](O)CC', 580, 110, highlight_atoms=[1], bond_px=58)
place(f, 'C[C@H](O)CC', 580, 270, highlight_atoms=[1], bond_px=58)
f.text(720, 120, '(R)-butan-2-ol', 23, weight='bold', anchor='start')
f.text(720, 280, '(S)-butan-2-ol', 23, weight='bold', anchor='start')
f.text(820, 200, '50 : 50', 26, color=RED)
f.text(480, 380, 'new stereocentre · achiral reagent · racemic product', 21, color=BLUE)
save('racemic', f)

# 11. NaBH4 scope: aldehydes and ketones yes, esters and amides no.
f = Figure('What NaBH4 reduces',
           'NaBH4 reduces a ketone such as cyclohexanone, but not ethyl acetate or acetamide under ordinary conditions', 960, 380)
for x, smi, hl, name, verdict, col in [
        (160, 'O=C1CCCCC1', [0, 1], 'cyclohexanone', 'reduced ✓', GREEN),
        (480, 'CC(=O)OCC', [1, 2, 3], 'ethyl acetate (ester)', 'no reaction ✗', RED),
        (800, 'CC(N)=O', [1, 2, 3], 'acetamide (amide)', 'no reaction ✗', RED)]:
    place(f, smi, x, 150, highlight_atoms=hl, bond_px=58)
    f.text(x, 275, name, 23, weight='bold').text(x, 315, verdict, 23, color=col, weight='bold')
f.text(480, 360, 'same C=O, different neighbours, different reactivity', 21, color=BLUE)
save('scope', f)

# 12. Pick the direction from the target.
f = Figure('Ethanol and ethanal interconvert',
           'Oxidation with PCC takes ethanol to ethanal; NaBH4 then water takes ethanal back to ethanol', 960, 380)
place(f, 'CCO', 190, 180, highlight_atoms=[1], bond_px=76)
place(f, 'CC=O', 770, 180, highlight_atoms=[1], bond_px=76)
f.line(330, 140, 630, 140, color=RED, arrow=True).line(630, 220, 330, 220, color=GREEN, arrow=True)
f.text(480, 122, 'oxidation (PCC): C–H lost, C–O gained', 21, color=RED)
f.text(480, 258, 'reduction (1. NaBH₄  2. H₂O): C–H gained', 21, color=GREEN)
f.text(190, 300, 'ethanol', 24, weight='bold').text(190, 340, 'C: −1', 24, color=BLUE)
f.text(770, 300, 'ethanal', 24, weight='bold').text(770, 340, 'C: +1', 24, color=BLUE)
save('route', f)

# 13. Activating OH: protonate it so water can leave. Hydronium is drawn by hand.
f = Figure('Protonation turns OH into a good leaving group',
           'Ethanol takes a proton from hydronium to give an oxonium ion; water can then leave', 960, 400)
e = place(f, 'CCO', 140, 200, bond_px=70, font_px=30)
ox, oy = e[2]
lone_pair(f, e[2], 90); lone_pair(f, e[2], 270)
hx, hy = ox + 110, oy - 85
f.text(hx, hy + 10, 'H', 28).line(hx + 13, hy, hx + 47, hy, width=2.5).text(hx + 82, hy + 10, 'OH₂⁺', 28, color=RED)
curved_arrow(f, (ox + 2, oy - 26), (hx - 6, hy + 4), bend=0.35, color=RED, shorten=6)
curved_arrow(f, (hx + 30, hy - 4), (hx + 70, hy - 22), bend=0.7, color=RED, shorten=2)
arrow(f, 470, 560, 190)
place_h(f, 'CC[OH2+]', 720, 190, bond_px=70, font_px=30, highlight_atoms=[2])
f.text(140, 300, 'ethanol', 23, weight='bold').text(hx + 45, hy + 55, 'H₃O⁺', 23, weight='bold')
f.text(720, 300, 'ethyloxonium ion', 23, weight='bold')
f.text(240, 360, 'leaving OH⁻: pKa(H₂O) ≈ 15.7 ✗', 21, color=RED)
f.text(720, 360, 'leaving H₂O: pKa(H₃O⁺) ≈ −1.7 ✓', 21, color=GREEN)
save('activation', f)

# 14. Aspirin synthesis amounts: n = m/M, then n / coefficient.
f = Figure('Amounts in the aspirin synthesis',
           'Salicylic acid 1.3812 g is 10.00 mmol and acetic anhydride 1.53135 g is 15.00 mmol, so salicylic acid limits aspirin to 10.00 mmol', 960, 460)
cols = (220, 490, 790)
place(f, 'OC(=O)c1ccccc1O', cols[0], 125, bond_px=42, font_px=22)
f.text(355, 135, '+', 34)
place(f, 'CC(=O)OC(C)=O', cols[1], 125, bond_px=42, font_px=22)
arrow(f, 610, 665, 130, '1 : 1')
place(f, 'CC(=O)Oc1ccccc1C(=O)O', cols[2], 125, bond_px=42, font_px=22)
rows = [('', 'salicylic acid', 'acetic anhydride', 'aspirin'),
        ('m', '1.3812 g', '1.53135 g', '?'),
        ('M', '138.12 g/mol', '102.09 g/mol', '180.16 g/mol'),
        ('n = m/M', '10.00 mmol', '15.00 mmol', 'at most 10.00 mmol')]
for r, row in enumerate(rows):
    y = 270 + 42*r
    f.text(30, y, row[0], 22, weight='bold', anchor='start')
    for x, cell in zip(cols, row[1:]):
        col = RED if (r == 3 and x == cols[0]) else INK
        f.text(x, y, cell, 22, color=col, weight='bold' if r == 0 else 'normal')
f.text(480, 445, 'compare n ÷ coefficient: 10.00 < 15.00, so salicylic acid limits', 21, color=BLUE)
save('aspirin-amounts', f)

# 15. Yield bars: theoretical, isolated, wet.
f = Figure('Isolated yield against the theoretical maximum',
           'Bars: theoretical maximum 1.8016 g, isolated 1.20 g (66.6 percent), wet sample 1.95 g (108 percent, impossible as pure product)', 960, 400)
base, scale = 320, 130
for x, val, label, sub, col in [(200, 1.8016, 'theoretical', '1.8016 g · 100%', BLUE),
                                (480, 1.20, 'isolated, dry', '1.20 g · 66.6%', GREEN),
                                (760, 1.95, 'isolated, wet', '1.95 g · "108%"', RED)]:
    h = val*scale
    f.rect(x - 70, base - h, 140, h, fill=col)
    f.text(x, base + 35, label, 23, weight='bold').text(x, base + 68, sub, 21, color=col)
f.line(100, base - 1.8016*scale, 860, base - 1.8016*scale, color=INK, width=2, dash='8 6')
f.text(880, base - 1.8016*scale + 7, 'max', 20, anchor='start')
f.line(100, base, 860, base, color=INK, width=2)
save('yield', f)

print('Built 15 alcohols-redox figures.')
