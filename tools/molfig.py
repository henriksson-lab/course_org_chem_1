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
    coords = {i: (drawer.GetDrawCoords(i).x, drawer.GetDrawCoords(i).y) for i in range(mol.GetNumAtoms())}
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # Atom notes (locants) in course blue; RDKit ignores the annotation colour option.
    svg = re.sub(r"(<path class='note' d='[^']*' fill=')#[0-9A-Fa-f]{6}", r"\g<1>" + BLUE, svg)
    # Highlighted atoms lose their element colour in RDKit; restore it.
    for i in highlight_atoms:
        col = ELEMENT_COLOURS.get(mol.GetAtomWithIdx(i).GetAtomicNum(), INK)
        svg = re.sub(r"(<path class='atom-%d' d='[^']*' fill=')#[0-9A-Fa-f]{6}" % i, r"\g<1>" + col, svg)
    body = svg[svg.index('<!-- END OF HEADER -->') + len('<!-- END OF HEADER -->'):svg.rindex('</svg>')]
    molecule_svg.last_coords = coords
    return body, w, h

def place(fig, smiles, x, y, **kwargs):
    """Embed a molecule into a chemfig.Figure centred on (x, y).

    Returns a dict {atom_index: (x, y)} of atom positions in figure coordinates,
    for attaching arrows, lone pairs or labels. Atom indices follow SMILES order
    (see atom_indices); with explicit_h=True the added H atoms come last.
    """
    body, w, h = molecule_svg(smiles, **kwargs)
    x0, y0 = x - w/2, y - h/2
    fig.parts.append(f'<svg x="{x0:.1f}" y="{y0:.1f}" width="{w:.1f}" height="{h:.1f}" '
                     f'viewBox="0 0 {w:.1f} {h:.1f}" overflow="visible">{body}</svg>')
    return {i: (x0 + cx, y0 + cy) for i, (cx, cy) in molecule_svg.last_coords.items()}

def curved_arrow(fig, start, end, bend=0.35, color=RED, width=2.5, shorten=14):
    """Electron-pushing arrow from start to end (figure coordinates).

    bend > 0 curves to the left of the start→end direction, < 0 to the right.
    shorten pulls both ends back so the arrow does not cover atom labels.
    """
    import math
    (x1, y1), (x2, y2) = start, end
    d = math.hypot(x2 - x1, y2 - y1) or 1
    ux, uy = (x2 - x1)/d, (y2 - y1)/d
    x1, y1, x2, y2 = x1 + ux*shorten, y1 + uy*shorten, x2 - ux*shorten, y2 - uy*shorten
    mx, my = (x1 + x2)/2, (y1 + y2)/2
    cx, cy = mx + uy*bend*d, my - ux*bend*d
    mid = 'head' + color.lstrip('#')
    if not any(f'id="{mid}"' in p for p in fig.parts):
        fig.parts.insert(1, f'<defs><marker id="{mid}" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" '
                            f'markerHeight="6" orient="auto-start-reverse"><path d="M 0 0 L 10 5 L 0 10 z" fill="{color}"/></marker></defs>')
    fig.parts.append(f'<path d="M {x1:.1f} {y1:.1f} Q {cx:.1f} {cy:.1f} {x2:.1f} {y2:.1f}" stroke="{color}" '
                     f'stroke-width="{width}" fill="none" marker-end="url(#{mid})"/>')

def lone_pair(fig, at, angle_deg, distance=24, color=BLUE):
    """Two dots at `distance` px from an atom position, in direction angle_deg (0 = right, 90 = up)."""
    import math
    a = math.radians(angle_deg)
    px, py = at[0] + distance*math.cos(a), at[1] - distance*math.sin(a)
    for side in (-1, 1):
        fig.circle(px + side*5*math.sin(a), py + side*5*math.cos(a), 3, color=color, fill=color, width=1)

def atom_indices(smiles):
    """Helper for authors: list atom index, symbol and neighbours of a SMILES."""
    mol = Chem.MolFromSmiles(smiles)
    return [(a.GetIdx(), a.GetSymbol(), [n.GetIdx() for n in a.GetNeighbors()]) for a in mol.GetAtoms()]
