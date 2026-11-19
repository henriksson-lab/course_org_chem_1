"""Skeletal-structure figures for the Day 1 structures lecture (RDKit via molfig)."""
from pathlib import Path
from chemfig import Figure, BLUE, RED, INK
from molfig import place

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'structures'/'figures'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

def panels(title, items, height=400, mol_y=170, label_y=None, sub_y=None, **common):
    """items: (smiles, label, sublabel, kwargs) laid out evenly across the canvas."""
    f = Figure(title, title, 960, height)
    n = len(items)
    label_y = label_y or height - 95
    sub_y = sub_y or height - 55
    for i, (smi, label, sub, kw) in enumerate(items):
        x = 960*(i + 0.5)/n
        place(f, smi, x, mol_y, **{**common, **kw})
        if label: f.text(x, label_y, label, 25, weight='bold')
        if sub: f.text(x, sub_y, sub, 22, color=BLUE)
    return f

# 1. Constitutional isomers: same formula, different oxygen neighbours.
f = panels('Ethanol and dimethyl ether share C2H6O', [
    ('CCO', 'ethanol', 'O bonded to C and H · bp 78 °C', dict(highlight_atoms=[2])),
    ('COC', 'dimethyl ether', 'O bonded to C and C · bp −24 °C', dict(highlight_atoms=[1])),
], mol_y=150, label_y=290, sub_y=330, bond_px=80)
f.text(480, 378, 'Both C₂H₆O', 25)
save('isomers', f)

# 2. Three representations of ethanol.
f = Figure('Ethanol as formula, condensed and skeletal structure', '', 960, 400)
f.text(160, 185, 'C₂H₆O', 40).text(480, 185, 'CH₃CH₂OH', 40)
place(f, 'CCO', 800, 170, bond_px=70)
for x, label, sub in [(160, 'molecular formula', 'which atoms'), (480, 'condensed', 'connectivity, as text'),
                      (800, 'skeletal', 'connectivity, as a drawing')]:
    f.text(x, 305, label, 25, weight='bold').text(x, 345, sub, 22, color=BLUE)
save('representations', f)

# 3. Pentane: every end and corner is a carbon.
f = Figure('Pentane: count ends and corners', '', 960, 400)
place(f, 'CCCCC', 300, 170, numbers={i: i + 1 for i in range(5)}, bond_px=70)
place(f, 'CCCCC', 720, 170, explicit_h=True, bond_px=48, font_px=20)
f.text(300, 320, 'skeletal: 5 corners and ends', 24, weight='bold')
f.text(720, 320, 'all atoms drawn', 24, weight='bold')
f.text(480, 370, '3 + 2 + 2 + 2 + 3 = 12 H  →  C₅H₁₂', 25, color=BLUE)
save('pentane', f)

# 4. Branches and rings.
save('branches', panels('Butane, 2-methylpropane and cyclobutane', [
    ('CCCC', 'butane', 'C₄H₁₀', {}),
    ('CC(C)C', '2-methylpropane', 'C₄H₁₀', dict(highlight_atoms=[1])),
    ('C1CCC1', 'cyclobutane', 'C₄H₈', {}),
], mol_y=160, bond_px=70))

# 5. Bond order and hydrogen count, all atoms shown.
save('bond-order', panels('Ethane, ethene and ethyne with every H drawn', [
    ('CC', 'ethane', 'C₂H₆', dict(explicit_h=True)),
    ('C=C', 'ethene', 'C₂H₄', dict(explicit_h=True)),
    ('C#C', 'ethyne', 'C₂H₂', dict(explicit_h=True)),
], mol_y=165, font_px=24))

# 6. Functional groups, each in a real molecule with the group shaded.
f = Figure('Six functional groups in real molecules', '', 960, 470)
groups = [
    ('CCO', 'alcohol', 'ethanol', [2]),
    ('CCN', 'amine', 'ethylamine', [2]),
    ('CCOCC', 'ether', 'diethyl ether', [2]),
    ('CC=O', 'aldehyde', 'ethanal', [1, 2]),
    ('CC(=O)O', 'carboxylic acid', 'acetic acid', [1, 2, 3]),
    ('CC(=O)OCC', 'ester', 'ethyl acetate', [1, 2, 3]),
]
for i, (smi, fam, name, hl) in enumerate(groups):
    x = 160 + 320*(i % 3); y = 80 + 230*(i//3)
    place(f, smi, x, y, highlight_atoms=hl, bond_px=46, font_px=24)
    f.text(x, y + 95, fam, 24, weight='bold').text(x, y + 125, name, 21, color=BLUE)
save('functional-groups', f)

# 7. What sits next to C=O decides the family.
save('carbonyl-family', panels('The neighbour of C=O decides the family', [
    ('CC=O', 'aldehyde', 'H', dict(highlight_atoms=[1, 2])),
    ('CC(C)=O', 'ketone', 'C', dict(highlight_atoms=[1, 3])),
    ('CC(=O)O', 'acid', 'OH', dict(highlight_atoms=[1, 2, 3])),
    ('CC(=O)OC', 'ester', 'OR', dict(highlight_atoms=[1, 2, 3])),
    ('CC(N)=O', 'amide', 'N', dict(highlight_atoms=[1, 2, 3])),
], mol_y=160, bond_px=52, font_px=26))

# 8. Naming: number from the end nearest the group.
save('naming', panels('Numbering puts the OH on the lowest carbon', [
    ('OCCCC', 'butan-1-ol', 'OH on C1', dict(numbers={1: 1, 2: 2, 3: 3, 4: 4}, highlight_atoms=[0])),
    ('CC(O)C', 'propan-2-ol', 'OH on C2', dict(numbers={0: 1, 1: 2, 3: 3}, highlight_atoms=[2])),
], height=330, mol_y=115, label_y=255, sub_y=292, bond_px=70))

# 9. Two molecules from biology.
save('biomolecules', panels('Glycine and lactic acid', [
    ('NCC(=O)O', 'glycine', 'amine + carboxylic acid', dict(highlight_atoms=[0, 2, 3, 4])),
    ('CC(O)C(=O)O', 'lactic acid', 'alcohol + carboxylic acid', dict(highlight_atoms=[2, 3, 4, 5])),
], mol_y=160, bond_px=70))

print('Built 9 skeletal structure figures.')
