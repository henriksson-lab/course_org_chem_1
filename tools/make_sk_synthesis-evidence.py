"""Figures for the Week 5 synthesis-and-evidence lecture (RDKit via molfig, plots via chemfig)."""
import math
from pathlib import Path
from chemfig import Figure, BLUE, RED, INK, GREEN
from molfig import place, molecule_svg

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'synthesis-evidence'/'figures'
GREY = '#8aa0ab'
PALE = '#dcebf1'
PALE_RED = '#f6e1e3'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

def extent(smiles, **kw):
    """Width and height of the atom positions (no label padding)."""
    molecule_svg(smiles, **kw)
    xs = [c[0] for c in molecule_svg.last_coords.values()]
    ys = [c[1] for c in molecule_svg.last_coords.values()]
    return max(xs) - min(xs), max(ys) - min(ys)

def centred(f, smiles, x, y, **kw):
    """place() centres the padded canvas; this centres the atoms themselves on (x, y)."""
    _, mw, mh = molecule_svg(smiles, **kw)
    xs = [c[0] for c in molecule_svg.last_coords.values()]
    ys = [c[1] for c in molecule_svg.last_coords.values()]
    return place(f, smiles, x + mw/2 - (min(xs) + max(xs))/2, y + mh/2 - (min(ys) + max(ys))/2, **kw)

def row(f, items, y, gap=22, pad=18):
    """Lay out a reaction left to right, centred on the 960-px canvas.

    items: ('mol', smiles, kwargs), ('text', str) or ('arrow', length, above[, colour]).
    Returns the x centre of each item.
    """
    widths = [extent(it[1], **it[2])[0] + 2*pad if it[0] == 'mol' else 26 if it[0] == 'text' else it[1]
              for it in items]
    x = (960 - sum(widths) - gap*(len(items) - 1))/2
    assert x > 0, f'row too wide: {960 - 2*x:.0f} px'
    centres = []
    for it, w in zip(items, widths):
        cx = x + w/2
        if it[0] == 'mol':
            centred(f, it[1], cx, y, **it[2])
        elif it[0] == 'text':
            f.text(cx, y + 12, it[1], 34)
        else:
            col = it[3] if len(it) > 3 else BLUE
            f.line(x, y, x + w, y, color=col, arrow=True)
            if it[2]: f.text(cx, y - 16, it[2], 21, color=col)
        centres.append(cx)
        x += w + gap
    return centres

# 1. Route planning: ethanol -> ethanoic acid -> ethyl ethanoate; the second ethanol portion sits on the arrow.
f = Figure('Two steps from ethanol to ethyl ethanoate',
           'Ethanol is oxidised to ethanoic acid, which is esterified with a second portion of ethanol', 960, 400)
k = dict(bond_px=44)
cx = row(f, [('mol', 'CCO', k), ('arrow', 110, ''), ('mol', 'CC(=O)O', k), ('arrow', 160, ''),
             ('mol', 'CC(=O)OCC', dict(k, highlight_atoms=[3, 4, 5]))], 200, gap=18)
f.text(cx[1], 184, '[O]', 22, color=BLUE).text(cx[1], 232, 'oxidation', 20, color=BLUE)
centred(f, 'OCC', cx[3], 130, bond_px=36, highlight_atoms=[0, 1, 2])
f.text(cx[3], 70, '+ ethanol', 20, color=BLUE)
f.text(cx[3], 232, 'H⁺, esterification', 20, color=BLUE)
for i, label in [(0, 'ethanol'), (2, 'ethanoic acid'), (4, 'ethyl ethanoate')]:
    f.text(cx[i], 310, label, 23, weight='bold')
f.text(480, 370, 'shaded: the O–CH₂CH₃ part comes from a second ethanol portion', 22, color=BLUE)
save('route', f)

# 2. Oxidation endpoint: ethanol -> ethanal -> ethanoic acid.
f = Figure('Oxidation of ethanol can stop at ethanal or go on to ethanoic acid',
           'Ethanol, ethanal and ethanoic acid in a row with two oxidation arrows; the second arrow is red', 960, 380)
k = dict(bond_px=56)
cx = row(f, [('mol', 'CCO', k), ('arrow', 130, '[O]'), ('mol', 'CC=O', dict(k, highlight_atoms=[1, 2])),
             ('arrow', 130, '[O]', RED), ('mol', 'CC(=O)O', dict(k, highlight_atoms=[1, 2, 3]))], 160, gap=24)
f.text(cx[3], 200, 'easy to', 19, color=RED).text(cx[3], 222, 'overshoot', 19, color=RED)
for i, label, sub in [(0, 'ethanol', 'alcohol'), (2, 'ethanal', 'aldehyde'), (4, 'ethanoic acid', 'carboxylic acid')]:
    f.text(cx[i], 290, label, 24, weight='bold').text(cx[i], 325, sub, 21, color=BLUE)
save('oxidation', f)

# 3. The aspirin reaction: acetic anhydride drawn over the arrow, molar masses underneath.
f = Figure('Salicylic acid plus acetic anhydride gives aspirin plus acetic acid',
           'Skeletal equation for aspirin synthesis: acetic anhydride over the arrow transfers its shaded acetyl group to the shaded phenolic oxygen; acetic acid is the by-product', 960, 460)
k = dict(bond_px=50)
cx = row(f, [('mol', 'OC(=O)c1ccccc1O', dict(k, rotate=267, highlight_atoms=[9])), ('arrow', 250, ''),
             ('mol', 'OC(=O)c1ccccc1OC(C)=O', dict(k, rotate=168, highlight_atoms=[9, 10, 11, 12]))], 240, gap=16)
centred(f, 'CC(=O)OC(C)=O', cx[1], 165, bond_px=42, font_px=24, highlight_atoms=[0, 1, 2])
f.text(cx[1], 280, '+ acetic anhydride', 23, color=BLUE).text(cx[1], 306, 'C₄H₆O₃ · 102.09 g/mol', 21, color=BLUE)
f.text(cx[1], 348, 'by-product: acetic acid', 23, color=GREY).text(cx[1], 374, 'C₂H₄O₂ · 60.05 g/mol', 21, color=GREY)
for i, name, formula in [(0, 'salicylic acid', 'C₇H₆O₃ · 138.12 g/mol'), (2, 'aspirin', 'C₉H₈O₄ · 180.16 g/mol')]:
    f.text(cx[i], 404, name, 25, weight='bold').text(cx[i], 436, formula, 22, color=BLUE)
save('aspirin', f)

# Bar helper for reaction-extent charts.
def bars(f, x0, base, items, scale, width=90, gap=150):
    """items: (label, value, valuetext, colour, outline). Bars grow upward from base."""
    for i, (label, value, vt, colour, outline) in enumerate(items):
        x = x0 + i*gap
        h = value*scale
        f.rect(x - width/2, base - h, width, h, fill=colour, color=outline)
        f.text(x, base - h - 12, vt, 22, weight='bold')
        f.text(x, base + 30, label, 21)

# 4. Limiting reagent: as planned and with less anhydride.
f = Figure('The smaller n divided by coefficient limits the product',
           'Bar charts of reaction extent: 10.0 versus 15.0 mmol, then 10.0 versus 5.00 mmol after reducing the anhydride', 960, 430)
base, s = 320, 14
for x0, title, a, b, lim, result in [
        (170, 'as planned', 10.0, 15.0, 0, 'max aspirin 10.0 mmol = 1.80 g'),
        (650, 'only 0.51045 g anhydride', 10.0, 5.0, 1, 'max aspirin 5.00 mmol = 0.901 g')]:
    f.text(x0 + 75, 50, title, 24, weight='bold')
    items = [('salicylic acid', a, f'{a:.1f} / 1', BLUE if lim == 0 else PALE, RED if lim == 0 else 'none'),
             ('anhydride', b, f'{b:.1f} / 1', GREEN if lim == 1 else PALE, RED if lim == 1 else 'none')]
    bars(f, x0, base, items, s)
    lx = x0 + 150*lim
    f.text(lx, base - (a if lim == 0 else b)*s/2 + 8, 'limits', 21, color='white', weight='bold')
    f.text(x0 + 75, 395, result, 22, color=RED)
f.line(480, 40, 480, 400, color=GREY, width=2, dash='6 6')
save('limiting', f)

# 5. Coefficients: 2A + B -> P.
f = Figure('With 2 A + B → P, divide each amount by its coefficient',
           'Paired bars: A 10 mmol raw becomes 5 mmol extent; B 8 mmol raw stays 8 mmol extent; A limits', 960, 440)
f.text(480, 48, '2 A + B → P', 30, weight='bold')
base, s = 320, 20
for x0, label, raw, coef in [(300, 'A', 10, 2), (660, 'B', 8, 1)]:
    ext = raw/coef
    f.rect(x0 - 110, base - raw*s, 90, raw*s, fill=PALE)
    f.text(x0 - 65, base - raw*s - 12, f'{raw} mmol', 22)
    col = RED if label == 'A' else BLUE
    f.rect(x0 + 20, base - ext*s, 90, ext*s, fill=col)
    f.text(x0 + 65, base - ext*s - 12, f'{raw}/{coef} = {ext:g}', 22, weight='bold')
    f.text(x0 - 65, base + 30, 'amount n', 20, color=GREY).text(x0 + 65, base + 30, 'n / coefficient', 20, color=BLUE)
    f.text(x0, base + 66, label, 28, weight='bold')
f.line(120, base, 840, base, color=INK, width=2)
f.text(480, 425, 'A limits: extent 5 mmol, so at most 5 mmol P', 23, color=RED)
save('coefficients', f)

# 6. Apparent isolated yield against isolated mass.
f = Figure('Apparent yield rises in proportion to the weighed mass',
           'Line plot of apparent isolated yield versus isolated mass with theoretical mass 1.8016 g; points at 1.10, 1.20 and 1.95 g; region above 100 percent shaded red', 960, 450)
X0, Y0, W, H = 130, 380, 740, 320          # axes origin and size
xmax, ymax = 2.2, 125
px = lambda m: X0 + m/xmax*W
py = lambda y: Y0 - y/ymax*H
f.rect(X0, py(ymax), W, py(100) - py(ymax), fill=PALE_RED)
f.text(X0 + 14, py(108) + 7, 'above 100 %: not all product', 20, color=RED, anchor='start')
f.line(X0, Y0, X0 + W + 10, Y0, arrow=True).line(X0, Y0, X0, py(ymax) - 10, arrow=True)
for m in [0, 0.5, 1.0, 1.5, 2.0]:
    f.line(px(m), Y0, px(m), Y0 + 8, width=2).text(px(m), Y0 + 32, f'{m:g}', 19)
for y in [0, 50, 100]:
    f.line(X0 - 8, py(y), X0, py(y), width=2).text(X0 - 14, py(y) + 7, f'{y}', 19, anchor='end')
f.text(X0 + W/2, 440, 'isolated mass / g', 21)
f.text(X0 - 14, py(ymax) - 18, 'yield / %', 20, anchor='start')
f.line(px(0), py(0), px(xmax), py(100*xmax/1.8016), color=BLUE, width=3)
f.line(px(1.8016), Y0, px(1.8016), py(100), color=GREY, width=2, dash='6 5')
f.text(px(1.8016) + 6, Y0 - 12, 'theoretical 1.8016 g', 18, color=GREY, anchor='start')
for m, col, dx, dy, anchor in [(1.10, BLUE, -14, -8, 'end'), (1.20, BLUE, 14, 30, 'start'), (1.95, RED, -16, 7, 'end')]:
    y = 100*m/1.8016
    f.circle(px(m), py(y), 8, color=col, fill=col)
    f.text(px(m) + dx, py(y) + dy, f'{m:.2f} g → {y:.1f} %', 20, color=col, anchor=anchor, weight='bold')
save('yield', f)

# TLC plate helper.
def plate(f, x, top, w, h, lanes, front_cm=5.0, show_front=True, labels=True):
    """lanes: list of (label, [spot_cm, ...], style). Returns (origin_y, cm->px)."""
    origin = top + h - 40
    cm = (origin - (top + 30))/front_cm
    f.rect(x, top, w, h, fill='white', color=GREY)
    f.line(x + 12, origin, x + w - 12, origin, color=INK, width=2, dash='4 5')
    if show_front:
        f.line(x + 12, origin - front_cm*cm, x + w - 12, origin - front_cm*cm, color=BLUE, width=2, dash='8 5')
    n = len(lanes)
    for i, (label, spots, style) in enumerate(lanes):
        lx = x + w*(i + 0.5)/n
        for s in spots:
            if style == 'faint':
                f.parts.append(f'<ellipse cx="{lx}" cy="{origin - s*cm}" rx="15" ry="10" fill="none" stroke="{GREY}" stroke-width="2" stroke-dasharray="4 3"/>')
            else:
                f.parts.append(f'<ellipse cx="{lx}" cy="{origin - s*cm}" rx="15" ry="10" fill="{style}"/>')
        if labels: f.text(lx, top + h + 32, label, 24, weight='bold')
    return origin, cm

# 7. TLC: start, product reference and crude.
f = Figure('A TLC plate with starting material, product reference and crude lanes',
           'Three lanes: S with a spot at 1.2 cm, P with a spot at 3.2 cm, and crude with spots at both heights', 960, 430)
origin, cm = plate(f, 300, 30, 360, 330, [('S', [1.2], INK), ('P', [3.2], BLUE), ('crude', [1.2, 3.2], RED)])
f.text(290, origin + 7, 'origin', 23, anchor='end')
f.text(290, origin - 5*cm + 7, 'solvent front', 23, color=BLUE, anchor='end')
f.text(680, origin - 3.2*cm + 7, 'matches P', 23, color=BLUE, anchor='start')
f.text(680, origin - 1.2*cm + 7, 'matches S', 23, color=INK, anchor='start')
f.text(150, 140, 'silica (polar)', 23, color=GREY)
f.text(150, 170, 'stays put', 23, color=GREY)
f.text(810, 60, 'crude: at least', 23, color=RED)
f.text(810, 88, 'two components', 23, color=RED)
save('tlc', f)

# 8. Rf from distances, and the missing front.
f = Figure('Rf is spot distance over front distance, both from the origin',
           'Left: one lane with arrows of 3.2 cm to the spot and 5.0 cm to the front, Rf 0.64. Right: a plate without a marked front, Rf cannot be calculated', 960, 400)
origin, cm = plate(f, 110, 25, 170, 310, [('', [3.2], BLUE)], labels=False)
f.line(315, origin, 315, origin - 3.2*cm, color=BLUE, arrow=True, width=4)
f.line(435, origin, 435, origin - 5.0*cm, color=INK, arrow=True, width=4)
f.text(325, origin - 1.6*cm, '3.2 cm', 26, color=BLUE, anchor='start')
f.text(445, origin - 2.6*cm, '5.0 cm', 26, anchor='start')
f.text(290, 380, 'Rf = 3.2 / 5.0 = 0.64', 28, weight='bold')
f.line(570, 25, 570, 390, color=GREY, width=2, dash='6 6')
o2, cm2 = plate(f, 650, 25, 170, 310, [('', [3.0], BLUE)], show_front=False, labels=False)
f.line(855, o2, 855, o2 - 3.0*cm2, color=BLUE, arrow=True, width=4)
f.text(865, o2 - 1.5*cm2, '3.0 cm', 26, color=BLUE, anchor='start')
f.text(735, 72, 'front?', 28, color=RED, weight='bold')
f.text(765, 380, 'no front: no Rf', 28, color=RED, weight='bold')
save('rf', f)

# 9. Co-migration and invisible components.
f = Figure('One spot can hide a mixture',
           'Lanes A, B and A plus B all show a single spot at Rf 0.64; a fourth lane shows a dashed outline for a component that the stain does not reveal', 960, 340)
origin, cm = plate(f, 230, 20, 500, 260, [('A', [3.2], BLUE), ('B', [3.2], GREEN), ('A + B', [3.2], INK), ('C', [1.8], 'faint')])
f.text(218, origin - 3.2*cm + 9, 'Rf 0.64', 26, anchor='end')
f.text(840, origin - 1.8*cm - 6, 'C: not seen', 24, color=GREY)
f.text(840, origin - 1.8*cm + 24, 'under UV', 24, color=GREY)
f.text(120, 185, 'same Rf,', 24, color=RED).text(120, 213, 'different', 24, color=RED).text(120, 241, 'compounds', 24, color=RED)
save('comigration', f)

# Spectrum helpers.
def ir_trace(f, x0, y0, w, h, bands, label, colour):
    """bands: (wavenumber, depth 0-1, half-width). x axis 4000 -> 500 cm-1, transmittance down."""
    pts = []
    for i in range(401):
        nu = 4000 - i*3500/400
        t = 1.0
        for c, d, s in bands:
            t -= d*math.exp(-((nu - c)/s)**2)
        pts.append((x0 + i*w/400, y0 + (1 - max(t, 0.02))*h))
    f.path('M ' + ' L '.join(f'{x:.1f} {y:.1f}' for x, y in pts), color=colour, width=3)
    f.text(x0 + w + 10, y0 + 20, label, 22, color=colour, anchor='start', weight='bold')
    return lambda nu: x0 + (4000 - nu)/3500*w

# 10. Schematic IR: salicylic acid vs aspirin.
f = Figure('Schematic IR: what changes from salicylic acid to aspirin',
           'Two schematic IR traces. Both keep a broad acid O-H band. Salicylic acid has a phenol O-H band; aspirin adds an ester C=O band near 1750 next to the acid C=O near 1690', 960, 460)
X0, W = 60, 720
common = [(2950, 0.30, 330), (1610, 0.30, 15), (1480, 0.25, 15), (1300, 0.35, 25), (760, 0.35, 20)]
xs = ir_trace(f, X0, 60, W, 150, common + [(3230, 0.45, 70), (1660, 0.80, 22)], 'salicylic acid', INK)
ir_trace(f, X0, 250, W, 150, common + [(1750, 0.80, 18), (1690, 0.75, 18), (1190, 0.55, 22)], 'aspirin', BLUE)
f.line(X0, 420, X0 + W, 420, width=2)
for nu in [4000, 3000, 2000, 1500, 1000]:
    f.line(xs(nu), 420, xs(nu), 428, width=2).text(xs(nu), 450, str(nu), 18)
f.text(X0 + W + 10, 450, 'cm⁻¹', 18, anchor='start')
f.text(xs(3230), 52, 'phenol O–H', 19, color=RED)
f.text(xs(2700), 245, 'broad acid O–H (both)', 19, color=INK)
f.text(xs(1760) - 8, 412, 'ester C=O ~1750', 19, color=RED, anchor='end')
f.text(xs(1690) + 6, 412, 'acid C=O ~1690', 19, color=BLUE, anchor='start')
f.text(xs(1640), 212, 'acid C=O ~1660', 19, color=INK, anchor='start')
save('ir', f)

# 11. Schematic 1H NMR of an aspirin sample with one unexplained signal.
f = Figure('Schematic proton NMR of an isolated aspirin sample with one extra signal',
           'Stick spectrum from 12 to 0 ppm: broad 1H acid signal near 11, four 1H aromatic signals between 7.1 and 8.1, a 3H singlet at 2.3, and an extra red singlet near 2.1', 960, 400)
X0, W, base = 60, 840, 300
xp = lambda d: X0 + (12 - d)/12*W
f.line(X0, base, X0 + W, base, width=2)
for d in range(12, -1, -2):
    f.line(xp(d), base, xp(d), base + 8, width=2).text(xp(d), base + 32, str(d), 19)
f.text(X0 + W/2, base + 66, 'chemical shift δ / ppm', 21)
unit = 55
f.path(f'M {xp(11.6):.1f} {base} Q {xp(11):.1f} {base - 1.4*unit:.1f} {xp(10.4):.1f} {base}', color=BLUE, width=3)
f.text(xp(11), base - 1.0*unit - 14, '1H, broad', 19, color=BLUE).text(xp(11), base - 1.0*unit - 40, 'COOH', 19, weight='bold')
for d in [8.12, 7.62, 7.34, 7.14]:
    f.line(xp(d), base, xp(d), base - unit, color=BLUE, width=3)
f.text(xp(7.6), base - unit - 14, '4 × 1H', 19, color=BLUE).text(xp(7.6), base - unit - 40, 'aromatic', 19, weight='bold')
f.line(xp(2.33), base, xp(2.33), base - 3*unit, color=BLUE, width=3)
f.text(xp(2.33) + 10, base - 3*unit + 10, '3H, s', 19, color=BLUE, anchor='start')
f.text(xp(2.33) + 10, base - 3*unit + 36, 'acetyl CH₃', 19, anchor='start', weight='bold')
f.line(xp(2.05), base, xp(2.05), base - 1.6*unit, color=RED, width=3)
f.text(xp(2.05) + 10, base - 1.6*unit + 30, 'extra: ?', 20, color=RED, anchor='start', weight='bold')
f.text(480, 385, 'schematic, not a measured spectrum', 18, color=GREY)
save('nmr', f)

print('Built 11 synthesis-evidence figures.')
