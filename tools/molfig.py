"""Skeletal structure drawings from SMILES, rendered with RDKit.

`molecule_svg` returns a self-contained <svg> fragment sized to the molecule.
`place` embeds it into a chemfig.Figure canvas centred on (x, y), so panel
labels, arrows and captions keep using the existing Figure helpers.

Requires RDKit (`pip install rdkit`).
"""
import re
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from chemfig import INK, BLUE, RED, GREEN

rdDepictor.SetPreferCoordGen(False)

def _rgb(hexcolour):
    h = hexcolour.lstrip('#')
    return tuple(int(h[i:i+2], 16)/255 for i in (0, 2, 4))

HIGHLIGHT = '#cfe3ec'
ELEMENT_COLOURS = {6: INK, 1: INK, 8: RED, 7: BLUE, 16: '#8a6d00', 17: GREEN, 35: '#7a3b2e'}

def molecule_svg(smiles, numbers=None, highlight_atoms=(), highlight_bonds=None,
                 bond_px=60, font_px=28, rotate=0, flip=False, explicit_h=False):
    """Draw one molecule.

    numbers: {atom_index: label} drawn next to atoms (e.g. chain locants).
    highlight_atoms: atom indices to shade; bonds between two highlighted atoms
    are shaded too unless highlight_bonds is given explicitly.
    rotate: degrees, to orient the drawing; flip mirrors left-right.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles}')
    if explicit_h:
        mol = Chem.AddHs(mol)
    rdDepictor.Compute2DCoords(mol, canonOrient=True)
    if rotate or flip:
        import math
        conf = mol.GetConformer()
        a = math.radians(rotate)
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            x = -p.x if flip else p.x
            conf.SetAtomPosition(i, (x*math.cos(a) - p.y*math.sin(a), x*math.sin(a) + p.y*math.cos(a), 0))
    for idx, label in (numbers or {}).items():
        mol.GetAtomWithIdx(idx).SetProp('atomNote', str(label))
    highlight_atoms = list(highlight_atoms)
    if highlight_bonds is None:
        hs = set(highlight_atoms)
        highlight_bonds = [b.GetIdx() for b in mol.GetBonds()
                           if b.GetBeginAtomIdx() in hs and b.GetEndAtomIdx() in hs]
    # Fixed bond length on a canvas sized to the molecule keeps every panel at one scale.
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    scale = bond_px/1.5
    margin = 2.2*font_px
    w = int((max(xs) - min(xs))*scale + 2*margin)
    h = int((max(ys) - min(ys))*scale + 2*margin)
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    o = drawer.drawOptions()
    o.fixedBondLength = bond_px
    o.clearBackground = False
    o.bondLineWidth = 2.5
    o.minFontSize = o.maxFontSize = font_px
    o.annotationFontScale = 0.65
    o.padding = 0.0
    o.updateAtomPalette({k: _rgb(v) for k, v in ELEMENT_COLOURS.items()})
    o.setSymbolColour(_rgb(INK))
    o.highlightBondWidthMultiplier = 12
    o.highlightRadius = 0.35
    colour = _rgb(HIGHLIGHT)
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, highlightBonds=list(highlight_bonds),
                        highlightAtomColors={i: colour for i in highlight_atoms},
                        highlightBondColors={i: colour for i in highlight_bonds})
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Atom notes (locants) in course blue; RDKit ignores the annotation colour option.
    svg = re.sub(r"(<path class='note' d='[^']*' fill=')#[0-9A-Fa-f]{6}", r"\g<1>" + BLUE, svg)
    # Highlighted atoms lose their element colour in RDKit; restore it.
    for i in highlight_atoms:
        col = ELEMENT_COLOURS.get(mol.GetAtomWithIdx(i).GetAtomicNum(), INK)
        svg = re.sub(r"(<path class='atom-%d' d='[^']*' fill=')#[0-9A-Fa-f]{6}" % i, r"\g<1>" + col, svg)
    body = svg[svg.index('<!-- END OF HEADER -->') + len('<!-- END OF HEADER -->'):svg.rindex('</svg>')]
    return body, w, h

def place(fig, smiles, x, y, **kwargs):
    """Embed a molecule into a chemfig.Figure centred on (x, y). Returns (w, h)."""
    body, w, h = molecule_svg(smiles, **kwargs)
    fig.parts.append(f'<svg x="{x - w/2:.1f}" y="{y - h/2:.1f}" width="{w:.1f}" height="{h:.1f}" '
                     f'viewBox="0 0 {w:.1f} {h:.1f}" overflow="visible">{body}</svg>')
    return w, h

def atom_indices(smiles):
    """Helper for authors: list atom index, symbol and neighbours of a SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    return [(a.GetIdx(), a.GetSymbol(), [n.GetIdx() for n in a.GetNeighbors()]) for a in mol.GetAtoms()]
