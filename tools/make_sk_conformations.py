"""Figures for the Week 2 conformations lecture.

Skeletal structures come from RDKit (molfig). Newman projections, chairs and
energy profiles are drawn with chemfig: the chairs are projected from ideal 3D
chair coordinates, so axial and equatorial bonds point the right way.
"""
import math
import numpy as np
from pathlib import Path
from chemfig import Figure, BLUE, RED, GREEN, INK
from molfig import place

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'conformations'/'figures'
GREY = '#6b7f8a'
LIGHT = '#cfe3ec'

def save(name, f):
    f.save(OUT/f'sk-{name}.svg')

# ---------------------------------------------------------------- Newman
def newman(f, cx, cy, front, rear, dihedral, r=46, arm=112, size=24, front_col=BLUE, rear_col=INK,
           highlight=None):
    """Newman projection. front/rear: three labels, first one at 12 o'clock (front)
    and at `dihedral` degrees clockwise from it (rear). Eclipsed rear bonds are
    nudged 14° so both sets stay visible."""
    shown = dihedral + (14 if dihedral % 120 == 0 else 0)
    f.circle(cx, cy, r, color=INK, fill='white', width=3)
    for j, lab in enumerate(rear):
        a = math.radians(-90 + shown + j*120)
        f.line(cx + r*math.cos(a), cy + r*math.sin(a), cx + arm*math.cos(a), cy + arm*math.sin(a),
               color=rear_col, width=3)
        d = arm + (18 + 6*len(lab))
        f.text(cx + d*math.cos(a), cy + d*math.sin(a) + size*0.35, lab, size, color=rear_col)
    for j, lab in enumerate(front):
        a = math.radians(-90 + j*120)
        f.line(cx, cy, cx + (arm - 22)*math.cos(a), cy + (arm - 22)*math.sin(a), color=front_col, width=3.5)
        d = arm - 22 + (16 + 6*len(lab))
        if dihedral % 120 == 0:
            d -= 4
        f.text(cx + d*math.cos(a), cy + d*math.sin(a) + size*0.35, lab, size, color=front_col, weight='bold')
    f.circle(cx, cy, 5, color=INK, fill=INK)

# ---------------------------------------------------------------- chairs
CC, ANG = 1.54, math.radians(109.5)
R_RING = 2*CC*math.sin(ANG/2)/math.sqrt(3)
H_RING = math.sqrt(CC**2 - R_RING**2)/2

def _unit(v):
    n = math.sqrt(sum(c*c for c in v)); return [c/n for c in v]

def _cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

class Chair:
    """Ideal chair. Atom j sits at angle -60*j (0 = right end, 1 = front right, 2 = front left,
    3 = left end, 4 = back left, 5 = back right). Even atoms are low, odd atoms high,
    unless flipped. Screen projection views the ring from slightly above."""
    def __init__(self, cx, cy, scale=95, flipped=False, elev=13, azim=7):
        self.cx, self.cy, self.s = cx, cy, scale
        self.e = math.radians(elev)
        self.az = math.radians(azim)
        sign = -1 if flipped else 1
        self.xyz = []
        for j in range(6):
            t = math.radians(-60*j)
            y = (-H_RING if j % 2 == 0 else H_RING)*sign
            self.xyz.append([R_RING*math.cos(t), y, R_RING*math.sin(t)])  # z > 0 is the back
        self.flipped = flipped

    def p(self, v):
        x0, y, z0 = v
        x = x0*math.cos(self.az) - z0*math.sin(self.az)
        z = x0*math.sin(self.az) + z0*math.cos(self.az)
        return (self.cx + self.s*x, self.cy - self.s*(y*math.cos(self.e) + z*math.sin(self.e)))

    def atom(self, j):
        return self.p(self.xyz[j])

    def draw(self, f, color=INK, width=3.5, skip=()):
        for j in range(6):
            if j in skip or (j + 1) % 6 in skip and False:
                continue
            a, b = self.atom(j), self.atom((j + 1) % 6)
            f.line(round(a[0], 1), round(a[1], 1), round(b[0], 1), round(b[1], 1), color=color, width=width)
        return self

    def direction(self, j, kind):
        """Unit 3D vector of the axial or equatorial bond at atom j."""
        a = self.xyz[j]
        u1 = _unit([n - c for n, c in zip(self.xyz[(j - 1) % 6], a)])
        u2 = _unit([n - c for n, c in zip(self.xyz[(j + 1) % 6], a)])
        b = _unit([-(p + q) for p, q in zip(u1, u2)])
        n = _unit(_cross(u1, u2))
        al = ANG/2
        d1 = [bb*math.cos(al) + nn*math.sin(al) for bb, nn in zip(b, n)]
        d2 = [bb*math.cos(al) - nn*math.sin(al) for bb, nn in zip(b, n)]
        ax, eq = (d1, d2) if abs(d1[1]) > abs(d2[1]) else (d2, d1)
        return ax if kind == 'ax' else eq

    def is_up(self, j, kind):
        return self.direction(j, kind)[1] > 0

    def sub(self, f, j, kind, label='', length=0.95, color=INK, size=24, width=3, lab_dist=0.42, dx=0, dy=0):
        """Draw a substituent bond; returns its end point on screen."""
        a = self.xyz[j]
        d = self.direction(j, kind)
        end = [c + length*e for c, e in zip(a, d)]
        x1, y1 = self.atom(j); x2, y2 = self.p(end)
        f.line(round(x1, 1), round(y1, 1), round(x2, 1), round(y2, 1), color=color, width=width)
        if label:
            lx, ly = self.p([c + (length + lab_dist)*e for c, e in zip(a, d)])
            f.text(round(lx + dx, 1), round(ly + size*0.35 + dy, 1), label, size, color=color)
        return x2, y2

# ---------------------------------------------------------------- plotting helpers
def axes(f, x0, y0, x1, y1, xlabel, ylabel):
    f.line(x0, y1, x1, y1, width=2.5, arrow=True)
    f.line(x0, y1, x0, y0, width=2.5, arrow=True)
    f.text((x0 + x1)/2, y1 + 62, xlabel, 22)
    f.parts.append(f'<text x="{x0 - 50}" y="{(y0 + y1)/2}" font-family="DejaVu Sans, sans-serif" font-size="22" '
                   f'fill="{INK}" text-anchor="middle" transform="rotate(-90 {x0 - 50} {(y0 + y1)/2})">{ylabel}</text>')

def curve(f, pts, color=BLUE, width=4, dash=None):
    d = 'M ' + ' L '.join(f'{x:.1f},{y:.1f}' for x, y in pts)
    extra = f' stroke-dasharray="{dash}"' if dash else ''
    f.parts.append(f'<path d="{d}" stroke="{color}" stroke-width="{width}" fill="none"{extra}/>')

RT = 8.314462618*298.15  # J/mol

# ======================================================================
# 1. Conformers vs constitutional isomers.
f = Figure('Conformers versus constitutional isomers',
           'Anti and gauche butane are conformers; butane and 2-methylpropane are constitutional isomers', 960, 440)
f.rect(20, 20, 450, 400, fill='#eef5f8')
f.rect(490, 20, 450, 400, fill='#f3f0ea')
f.text(245, 60, 'conformers: rotate a bond', 25, weight='bold', color=BLUE)
newman(f, 135, 200, ['CH₃', 'H', 'H'], ['CH₃', 'H', 'H'], 180, r=34, arm=78, size=21)
newman(f, 355, 200, ['CH₃', 'H', 'H'], ['CH₃', 'H', 'H'], 60, r=34, arm=78, size=21)
f.text(135, 362, 'anti', 23, weight='bold').text(355, 362, 'gauche', 23, weight='bold')
f.text(245, 402, 'same bonds, different shape', 22, color=BLUE)
f.text(715, 60, 'isomers: break and remake bonds', 25, weight='bold', color=RED)
place(f, 'CCCC', 610, 190, bond_px=60)
place(f, 'CC(C)C', 830, 190, bond_px=60, highlight_atoms=[1])
f.text(610, 300, 'butane', 23, weight='bold').text(830, 300, '2-methylpropane', 23, weight='bold')
f.text(715, 345, 'both C₄H₁₀', 22).text(715, 402, 'different bonded neighbours', 22, color=RED)
save('conformer-isomer', f)

# 2. Choosing the bond and the view: butane C2→C3.
f = Figure('Looking along the C2–C3 bond of butane',
           'Numbered butane with an eye looking from C2 toward C3, and the resulting Newman projection', 960, 420)
at = place(f, 'CCCC', 250, 190, numbers={0: 'C1', 1: 'C2', 2: 'C3', 3: 'C4'}, bond_px=90, highlight_atoms=[1, 2])
(x1, y1), (x2, y2) = at[1], at[2]
ux, uy = (x2 - x1), (y2 - y1); n = math.hypot(ux, uy); ux, uy = ux/n, uy/n
ex, ey = x1 - 120*ux, y1 - 120*uy
f.line(round(ex, 1), round(ey, 1), round(x1 - 30*ux, 1), round(y1 - 30*uy, 1), color=BLUE, width=4, arrow=True)
f.parts.append(f'<ellipse cx="{ex - 22*ux:.1f}" cy="{ey - 22*uy:.1f}" rx="20" ry="12" fill="white" stroke="{BLUE}" stroke-width="3"/>'
               f'<circle cx="{ex - 22*ux:.1f}" cy="{ey - 22*uy:.1f}" r="6" fill="{BLUE}"/>')
f.text(250, 345, 'eye at C2, looking toward C3', 24)
f.line(470, 200, 560, 200, color=GREY, width=3, arrow=True)
newman(f, 735, 195, ['CH₃', 'H', 'H'], ['CH₃', 'H', 'H'], 180)
f.text(735, 395, 'front C2 = dot (blue) · rear C3 = circle', 22)
save('butane-view', f)

# 3. Ethane: staggered and eclipsed.
f = Figure('Staggered and eclipsed ethane', 'Newman projections of staggered (60 degree offset) and eclipsed (0 degree offset) ethane', 960, 400)
newman(f, 250, 185, ['H', 'H', 'H'], ['H', 'H', 'H'], 60)
newman(f, 710, 185, ['H', 'H', 'H'], ['H', 'H', 'H'], 0)
f.text(250, 355, 'staggered · rear bonds 60° off', 24, weight='bold')
f.text(710, 355, 'eclipsed · rear bonds lined up', 24, weight='bold')
f.text(710, 388, '(drawn slightly offset so you can see them)', 19, color=GREY)
save('ethane-newman', f)

# 4. Ethane energy profile, 12 kJ/mol barrier.
f = Figure('Ethane torsional energy', 'Energy versus H–C–C–H dihedral angle: three equal minima (staggered) and 12 kJ/mol maxima (eclipsed)', 960, 420)
X0, X1, Y0, Y1 = 130, 900, 50, 330
sx = lambda a: X0 + (X1 - X0 - 30)*a/360
sy = lambda e: Y1 - 40 - 15*e
axes(f, X0, Y0, X1, Y1, 'H–C–C–H dihedral angle (°)', 'energy (kJ/mol)')
curve(f, [(sx(a), sy(6*(1 + math.cos(math.radians(3*a))))) for a in range(0, 361, 2)])
for a in range(0, 361, 60):
    f.text(sx(a), Y1 + 28, str(a), 19)
for e in (0, 4, 8, 12):
    f.line(X0 - 6, sy(e), X0, sy(e), width=2).text(X0 - 14, sy(e) + 7, str(e), 18, anchor='end')
f.line(sx(120), sy(12), sx(120), sy(0), color=RED, width=2, dash='6 5')
f.text(sx(120) + 10, sy(6), '12 kJ/mol', 21, color=RED, anchor='start')
f.text(sx(0) + 8, sy(12) - 14, 'eclipsed', 20, anchor='start')
for a in (60, 180, 300):
    f.text(sx(a), sy(0) + 28, 'staggered', 20, color=GREEN)
save('ethane-energy', f)

# 5. Butane through a full turn.
f = Figure('Butane through one full turn', 'Six Newman projections of butane along C2–C3 at methyl dihedrals 0, 60, 120, 180, 240 and 300 degrees', 960, 585)
names = [(0, 'syn', 'eclipsed', 'CH₃/CH₃'), (60, 'gauche', 'staggered', ''), (120, '', 'eclipsed', 'CH₃/H'),
         (180, 'anti', 'staggered', ''), (240, '', 'eclipsed', 'CH₃/H'), (300, 'gauche', 'staggered', '')]
for i, (a, nm, kind, note) in enumerate(names):
    x, y = 160 + 320*(i % 3), 115 + 285*(i//3)
    newman(f, x, y, ['CH₃', 'H', 'H'], ['CH₃', 'H', 'H'], a, r=26, arm=60, size=22)
    col = RED if kind == 'eclipsed' else GREEN
    f.text(x, y + 128, f'{a}° · {kind}', 27, weight='bold', color=col)
    f.text(x, y + 158, nm or f'{note} eclipsed', 25, color=BLUE if nm else RED)
save('butane-turn', f)

# 6. Butane energy profile with values.
E = {0: 19.0, 60: 3.8, 120: 16.0, 180: 0.0}
# Fit E(φ) = a0 + a1 cos φ + a2 cos 2φ + a3 cos 3φ through the four reference points.
A = np.array([[1, math.cos(math.radians(p)), math.cos(math.radians(2*p)), math.cos(math.radians(3*p))] for p in E])
coef = np.linalg.solve(A, np.array(list(E.values())))
Ef = lambda p: float(coef @ [1, math.cos(math.radians(p)), math.cos(math.radians(2*p)), math.cos(math.radians(3*p))])
f = Figure('Butane torsional energy', 'Energy versus methyl dihedral: anti 0, gauche 3.8, CH3/H eclipsed 16, CH3/CH3 eclipsed 19 kJ/mol', 960, 440)
X0, X1, Y0, Y1 = 130, 900, 40, 350
sy = lambda e: Y1 - 40 - 12*e
axes(f, X0, Y0, X1, Y1, 'CH₃–C2–C3–CH₃ dihedral angle (°)', 'energy (kJ/mol)')
curve(f, [(sx(a), sy(Ef(a))) for a in range(0, 361, 2)])
for a in range(0, 361, 60):
    f.text(sx(a), Y1 + 28, str(a), 19)
for e in (0, 5, 10, 15, 20):
    f.line(X0 - 6, sy(e), X0, sy(e), width=2).text(X0 - 14, sy(e) + 7, str(e), 18, anchor='end')
f.text(sx(0) + 10, sy(19) - 6, '19 · CH₃/CH₃ eclipsed', 20, color=RED, anchor='start')
f.text(sx(120), sy(16) - 14, '16 · CH₃/H eclipsed', 20, color=RED)
f.text(sx(240), sy(16) - 14, '16', 20, color=RED)
f.text(sx(60), sy(3.8) + 30, 'gauche 3.8', 20, color=GREEN)
f.text(sx(300), sy(3.8) + 30, 'gauche 3.8', 20, color=GREEN)
f.text(sx(180), sy(0) + 28, 'anti 0', 21, color=GREEN, weight='bold')
save('butane-energy', f)

# 7. Boltzmann populations for three free-energy gaps.
f = Figure('Populations from free-energy gaps', 'Two-state populations at 298 K for gaps of 0, 3500 and 7000 J/mol', 960, 420)
X0, X1, Y0, Y1 = 120, 920, 50, 320
axes(f, X0, Y0, X1, Y1, '', 'fraction')
for v in (0, 0.5, 1.0):
    f.line(X0 - 6, Y1 - 250*v, X0, Y1 - 250*v, width=2).text(X0 - 14, Y1 - 250*v + 7, f'{v:g}', 18, anchor='end')
for i, g in enumerate((0, 3500, 7000)):
    w = math.exp(-g/RT); pa, pb = 1/(1 + w), w/(1 + w)
    cx = 250 + 270*i
    for k, (p, nm, col) in enumerate(((pa, 'A', BLUE), (pb, 'B', '#7fb3c9'))):
        x = cx - 75 + 80*k
        f.rect(x, Y1 - 250*p, 70, 250*p, fill=col)
        f.text(x + 35, Y1 - 250*p - 10, f'{p:.3f}', 20)
        f.text(x + 35, Y1 + 26, nm, 20, weight='bold')
    f.text(cx - 5, Y1 + 62, f'ΔG(B) = {g} J/mol', 21, color=BLUE)
f.text(520, 40, 'A has ΔG = 0 in every case · T = 298.15 K', 20, color=GREY)
save('populations', f)

# 8. Equilibrium versus rate.
f = Figure('Equilibrium versus rate', 'Two energy profiles with the same free-energy gap between A and B but different barriers', 960, 420)
X0, X1, Y0, Y1 = 110, 900, 40, 330
axes(f, X0, Y0, X1, Y1, 'rotation coordinate', 'energy')
def hump(barrier, gap):
    pts = []
    for i in range(101):
        t = i/100
        base = gap*(1 - math.cos(math.pi*t))/2
        pts.append((X0 + 60 + (X1 - X0 - 140)*t, Y1 - 30 - 9*(base + barrier*math.sin(math.pi*t)**2)))
    return pts
GAP = 6  # schematic units
curve(f, hump(12, GAP), color=BLUE)
curve(f, hump(26, GAP), color=RED, dash='10 7')
f.text(X0 + 60, Y1 - 2, 'A', 24, weight='bold').text(X1 - 80, Y1 - 30 - 9*GAP - 12, 'B', 24, weight='bold')
f.line(X0 + 60, Y1 - 30, X1 - 30, Y1 - 30, color=GREEN, width=2, dash='5 5')
f.line(X1 - 40, Y1 - 30, X1 - 40, Y1 - 30 - 9*GAP, color=GREEN, width=3, both=True)
f.text(X1 - 50, Y1 - 4, 'ΔG between A and B: sets the ratio', 20, color=GREEN, anchor='end')
f.text(380, Y1 - 30 - 9*22, 'high barrier: slow', 21, color=RED, anchor='end')
f.text(640, Y1 - 30 - 9*2.2, 'low barrier: fast', 21, color=BLUE, anchor='start')
save('equilibrium-rate', f)

# 9. Rings: a flat drawing versus the real chair.
f = Figure('Flat ring drawing versus the chair', 'Skeletal cyclopropane and cyclohexane, and the chair conformation of cyclohexane', 960, 400)
place(f, 'C1CC1', 150, 165, bond_px=80)
place(f, 'C1CCCCC1', 460, 165, bond_px=65)
c = Chair(790, 175, scale=85); c.draw(f)
f.text(150, 295, 'cyclopropane', 23, weight='bold').text(150, 330, '60° angles: strained', 21, color=RED)
f.text(460, 295, 'cyclohexane, as drawn', 23, weight='bold').text(460, 330, 'if flat: 120° angles', 21, color=RED)
f.text(790, 295, 'real shape: chair', 23, weight='bold').text(790, 330, '~109.5°, almost no strain', 21, color=GREEN)
f.line(590, 170, 650, 170, color=GREY, width=3, arrow=True)
save('rings', f)

# 10. Chair positions: six axial (alternating up/down) and six equatorial.
f = Figure('Axial and equatorial positions in a chair', 'Chair cyclohexane with six axial bonds alternating up and down in blue and six equatorial bonds around the rim in green', 960, 440)
c = Chair(330, 215, scale=110); c.draw(f)
for j in range(6):
    c.sub(f, j, 'ax', color=BLUE, length=0.95)
    c.sub(f, j, 'eq', color=GREEN, length=0.95)
c.sub(f, 3, 'ax', 'up', color=BLUE, size=20, lab_dist=0.25)
c.sub(f, 0, 'ax', 'down', color=BLUE, size=20, lab_dist=0.3)
f.text(740, 120, 'axial', 28, color=BLUE, weight='bold')
f.text(740, 158, 'vertical, alternating up/down', 21, color=BLUE)
f.text(740, 248, 'equatorial', 28, color=GREEN, weight='bold')
f.text(740, 286, 'around the rim, tilted slightly', 21, color=GREEN)
f.text(740, 316, 'opposite to that carbon’s axial', 21, color=GREEN)
f.text(740, 360, 'every carbon: one of each', 21)
save('chair-positions', f)

# 11. Ring flip: axial-up becomes equatorial-up.
f = Figure('A ring flip swaps axial and equatorial', 'Methylcyclohexane: axial-up methyl on the left becomes equatorial-up methyl after the ring flip', 960, 420)
a = Chair(200, 215, scale=82); a.draw(f)
b = Chair(770, 215, scale=82, flipped=True); b.draw(f)
j = 3  # left-end carbon: high in chair a, low in chair b
a.sub(f, j, 'ax', 'CH₃', color=RED, width=4, lab_dist=0.35)
a.sub(f, j, 'eq', 'H', color=GREY, size=21, length=0.6, lab_dist=0.3)
b.sub(f, j, 'eq', 'CH₃', color=RED, width=4, length=0.75, lab_dist=0.5, dy=-4)
b.sub(f, j, 'ax', 'H', color=GREY, size=21, length=0.6, lab_dist=0.3)
f.line(365, 205, 465, 205, color=BLUE, width=3, arrow=True).line(465, 235, 365, 235, color=BLUE, width=3, arrow=True)
f.text(415, 185, 'ring flip', 21, color=BLUE)
f.text(200, 345, 'axial · up', 25, weight='bold').text(770, 345, 'equatorial · up', 25, weight='bold')
f.text(480, 395, 'axial ⇄ equatorial changes; up stays up', 23, color=BLUE)
save('ring-flip', f)

# 12. Methylcyclohexane: 1,3-diaxial interactions.
f = Figure('Why methyl prefers equatorial', 'Axial methyl crowds the axial hydrogens on C3 and C5 (1,3-diaxial interactions); equatorial methyl points away', 960, 430)
a = Chair(240, 240, scale=95)
a.draw(f)
mx, my = a.sub(f, 3, 'ax', 'CH₃', color=RED, width=4, lab_dist=0.35)
hpts = [a.sub(f, k, 'ax', 'H', color=BLUE, size=22, lab_dist=0.3) for k in (1, 5)]
for hx, hy in hpts:
    f.line(round(mx, 1), round(my, 1), round(hx, 1), round(hy, 1), color=RED, width=2.5, dash='5 5')
for k, lab in ((3, '1'), (1, '3'), (5, '5')):
    x, y = a.atom(k); f.text(round(x, 1), round(y + 30, 1), lab, 19, color=GREY)
b = Chair(720, 240, scale=95, flipped=True); b.draw(f)
b.sub(f, 3, 'eq', 'CH₃', color=RED, width=4, lab_dist=0.5, dy=-4)
f.text(240, 375, 'axial: two 1,3-diaxial clashes', 23, weight='bold', color=RED)
f.text(720, 375, 'equatorial: points away', 23, weight='bold', color=GREEN)
f.text(480, 415, 'ΔG ≈ 7.3 kJ/mol → about 95 : 5 equatorial : axial at 25 °C', 22, color=BLUE)
save('methylcyclohexane', f)

# 13. cis/trans-1,2-dimethylcyclohexane: skeletal with wedges and as chairs.
f = Figure('cis and trans in a ring', 'cis- and trans-1,2-dimethylcyclohexane as wedge drawings and as chairs: cis is axial/equatorial, trans is diequatorial', 960, 450)
place(f, 'C[C@H]1CCCC[C@H]1C', 125, 160, bond_px=46)
place(f, 'C[C@H]1CCCC[C@@H]1C', 605, 160, bond_px=46)
c1 = Chair(315, 175, scale=54); c1.draw(f)
c1.sub(f, 1, 'ax', 'CH₃', color=RED, width=4, size=20, length=0.8, lab_dist=0.4)
c1.sub(f, 0, 'eq', 'CH₃', color=RED, width=4, size=20, length=0.7, lab_dist=0.6, dy=-4)
c2 = Chair(795, 175, scale=54); c2.draw(f)
c2.sub(f, 1, 'eq', 'CH₃', color=RED, width=4, size=20, length=0.7, lab_dist=0.5, dy=6)
c2.sub(f, 0, 'eq', 'CH₃', color=RED, width=4, size=20, length=0.7, lab_dist=0.6, dy=-4)
f.text(235, 310, 'cis: same side', 25, weight='bold').text(235, 345, 'chair: one axial, one equatorial', 21, color=BLUE)
f.text(715, 310, 'trans: opposite sides', 25, weight='bold').text(715, 345, 'chair: both equatorial (best)', 21, color=BLUE)
f.line(480, 30, 480, 420, color=GREY, width=1.5, dash='6 6')
f.text(480, 425, 'a ring flip changes axial/equatorial, never cis/trans', 22, color=RED)
save('cis-trans', f)

# 14. beta-D-glucose: every substituent equatorial.
f = Figure('β-D-glucose in its chair', 'Chair of beta-D-glucopyranose with all four OH groups and the CH2OH group equatorial and every axial position holding H', 960, 420)
g = Chair(360, 215, scale=125)
# atom index: 0 = C1 (right end), 1 = C2, 2 = C3, 3 = C4 (left end), 4 = C5, 5 = ring O
g.draw(f)
ox, oy = g.atom(5)
f.parts.append(f'<circle cx="{ox:.1f}" cy="{oy:.1f}" r="18" fill="#f7fafc"/>')
f.text(round(ox, 1), round(oy + 9, 1), 'O', 28, color=RED, weight='bold')
for j, lab, dx, dy in ((0, 'OH', 8, 0), (1, 'OH', 4, 6), (2, 'HO', -6, 8), (3, 'HO', -10, 0), (4, 'CH₂OH', -30, -4)):
    g.sub(f, j, 'eq', lab, color=RED if lab != 'CH₂OH' else INK, width=4, length=0.75, lab_dist=0.45, dx=dx, dy=dy, size=25)
x, y = g.atom(0); f.text(round(x + 4, 1), round(y + 36, 1), 'C1', 20, color=BLUE, weight='bold')
x, y = g.atom(4); f.text(round(x + 22, 1), round(y + 26, 1), 'C5', 20, color=BLUE, weight='bold')
f.text(360, 385, 'axial positions all carry H (not drawn)', 21, color=GREY)
f.text(800, 150, 'all five big groups', 23, weight='bold')
f.text(800, 185, 'equatorial', 23, weight='bold', color=GREEN)
f.text(800, 260, 'no 1,3-diaxial clashes', 21, color=BLUE)
save('glucose', f)

print('Built 14 conformation figures.')
