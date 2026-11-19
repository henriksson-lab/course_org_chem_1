"""Original SVG chemistry illustrations for the four foundation units."""
from pathlib import Path
import math
from chemfig import Figure, BLUE, RED, GREEN, INK
ROOT=Path(__file__).resolve().parents[1]

def save(unit,num,name,f):
    f.save(ROOT/unit/'figures'/f'{num:02d}-{name}.svg')

def bond(f,a,b,order=1,color=INK):
    x1,y1=a;x2,y2=b
    if order==1: f.line(x1,y1,x2,y2,color=color)
    else:
        d=math.hypot(x2-x1,y2-y1);dx=-(y2-y1)/d*5;dy=(x2-x1)/d*5
        f.line(x1+dx,y1+dy,x2+dx,y2+dy,color=color).line(x1-dx,y1-dy,x2-dx,y2-dy,color=color)
    return f

def zigzag(f,points,labels=False):
    for a,b in zip(points,points[1:]):bond(f,a,b)
    if labels:
        for i,(x,y) in enumerate(points):f.text(x,y+42,str(i+1),20,color=BLUE)
    return f

def axes(f,xlabel,ylabel):
    f.line(105,300,870,300,arrow=True).line(105,300,105,65,arrow=True)
    f.text(480,380,xlabel,23).text(180,45,ylabel,22)

def newman(f,cx,cy,front,back,angle=60,title=''):
    r=44; arm=110
    f.circle(cx,cy,r)
    for j,label in enumerate(back):
        a=math.radians(-90+angle+j*120)
        f.line(cx+r*math.cos(a),cy+r*math.sin(a),cx+arm*math.cos(a),cy+arm*math.sin(a))
        f.text(cx+(arm+42 if angle==0 else arm+30)*math.cos(a),cy+(arm+42 if angle==0 else arm+30)*math.sin(a)+8,label,23)
    f.circle(cx,cy,5,fill=INK)
    for j,label in enumerate(front):
        a=math.radians(-90+j*120)
        f.line(cx,cy,cx+(arm-10)*math.cos(a),cy+(arm-10)*math.sin(a),color=BLUE)
        f.text(cx+(arm-2 if angle==0 else arm+23)*math.cos(a),cy+(arm-2 if angle==0 else arm+23)*math.sin(a)+8,label,23,color=BLUE)
    if title:f.text(cx,365,title,24)

def fischer(f,x,y,left='OH',right='H',top='COOH',bottom='CH₃',label=''):
    f.line(x-70,y,x+70,y).line(x,y-65,x,y+65)
    f.text(x-105,y+8,left,25).text(x+105,y+8,right,25)
    f.text(x,y-85,top,26).text(x,y+105,bottom,26)
    if label:f.text(x,y+155,label,23)

# Structures
u='structures'
f=Figure('Ethanol in three complementary representations')
f.text(175,130,'C₂H₆O',40).text(175,235,'molecular formula',24)
f.text(480,130,'CH₃CH₂OH',35).text(480,235,'condensed structure',24)
zigzag(f,[(690,160),(750,110)]);f.line(750,110,800,140);f.text(825,160,'OH',30)
f.text(775,235,'skeletal structure',24).text(480,335,'Connectivity is visible in the last two representations.',25)
save(u,1,'representations',f)
f=Figure('Counting the hidden carbons and hydrogens in pentane')
zigzag(f,[(210,180),(340,110),(470,180),(600,110),(730,180)],True)
for x,y,h in [(210,180,'CH₃'),(340,110,'CH₂'),(470,180,'CH₂'),(600,110,'CH₂'),(730,180,'CH₃')]:f.text(x,y-35,h,24,color=BLUE)
f.text(480,330,'5 carbons; 3 + 2 + 2 + 2 + 3 = 12 hydrogens → C₅H₁₂',26)
save(u,2,'implicit-atoms',f)
f=Figure('Bond order controls the number of attached hydrogens')
for x,name,h,order in [(210,'ethane','H₃C     CH₃',1),(480,'ethene','H₂C     CH₂',2),(750,'ethyne','HC       CH',3)]:
 f.text(x,105,name,27).text(x-70,195,['H₃C','H₂C','HC'][order-1],28).text(x+70,195,['CH₃','CH₂','CH'][order-1],28)
 for dy in ([0] if order==1 else [-5,5] if order==2 else [-9,0,9]):f.line(x-25,185+dy,x+25,185+dy)
 f.text(x,285,['C₂H₆','C₂H₄','C₂H₂'][order-1],29)
f.text(480,350,'Neutral carbon has a total bond order of four.',25)
save(u,3,'bond-order',f)
f=Figure('Constitutional isomers share a formula but not connectivity')
f.text(245,170,'CH₃—CH₂—OH',35).text(720,170,'CH₃—O—CH₃',35)
f.text(245,245,'ethanol: alcohol',27).text(720,245,'dimethyl ether: ether',27)
f.text(480,330,'Both C₂H₆O; the oxygen has different neighbours.',27)
save(u,4,'constitutional-isomers',f)
f=Figure('Functional groups identify local patterns')
for x,y,s,name in [(180,120,'R—OH','alcohol'),(480,120,'R—NH₂','amine'),(780,120,'R—O—R′','ether'),(180,270,'R—C(=O)H','aldehyde'),(480,270,'R—C(=O)OH','carboxylic acid'),(780,270,'R—C(=O)OR′','ester')]:
 f.text(x,y,s,25).text(x,y+50,name,23,color=BLUE)
save(u,5,'functional-groups',f)
f=Figure('An oxygen charge audit counts bonds and lone pairs')
for x,name,bonds,pairs,charge in [(190,'water',2,2,'0'),(480,'hydroxide',1,3,'−1'),(770,'hydronium',3,1,'+1')]:
 f.text(x,42,name,26).text(x,180,'O',38)
 for j in range(bonds):
  a=math.radians(25+j*110);f.line(x+20*math.cos(a),166+20*math.sin(a),x+72*math.cos(a),166+72*math.sin(a));f.text(x+93*math.cos(a),174+93*math.sin(a),'H',25)
 for angle in {2:[225,315],1:[115,205,295],3:[315]}[bonds]:
  a=math.radians(angle)
  px,py=x+35*math.cos(a),166+35*math.sin(a)
  for side in [-1,1]:
   f.circle(px-side*4*math.sin(a),py+side*4*math.cos(a),2.5,color=BLUE,fill=BLUE,width=1)
 if charge!='0': f.text(x+55,160,'−' if charge=='−1' else '+',26,color=RED)
 f.text(x,275,f'{bonds} bond{"s"*(bonds!=1)}; {pairs} lone pair{"s"*(pairs!=1)}',22).text(x,330,f'formal charge {charge}',24)
save(u,6,'oxygen-charge',f)
f=Figure('Geometry around carbon depends on its bonding')
for cx,geom,angle in [(200,'tetrahedral','≈109.5°'),(490,'trigonal planar','≈120°'),(780,'linear','180°')]:
 f.text(cx,60,geom,24).text(cx,200,'C',29)
 if cx==200:
  f.line(cx-17,185,cx-75,115).line(cx+17,185,cx+75,115).wedge(cx,205,cx-60,280).wedge(cx+10,205,cx+70,280,hashed=True)
 elif cx==490:
  for a in [-90,30,150]:
   rad=math.radians(a);f.line(cx+18*math.cos(rad),190+18*math.sin(rad),cx+100*math.cos(rad),190+100*math.sin(rad))
 else:f.line(cx-110,190,cx-20,190).line(cx+20,190,cx+110,190)
 f.text(cx,350,angle,26)
save(u,7,'geometry',f)
f=Figure('Name a chain by following connectivity')
zigzag(f,[(160,210),(310,120),(460,210),(610,120)])
f.line(610,120,665,153);f.text(695,170,'OH',30)
for x,y,n in [(160,210,'4'),(310,120,'3'),(460,210,'2'),(610,120,'1')]:f.text(x,y+45,n,23,color=BLUE)
f.text(480,325,'Four-carbon chain; number from the OH end → butan-1-ol',26)
save(u,8,'simple-naming',f)
f=Figure('A carbonyl carbon has three atom neighbours and four bond orders')
f.text(460,190,'C',36).text(460,70,'O',36).text(245,190,'CH₃',34).text(675,190,'OH',34)
bond(f,(450,150),(450,92),2);f.line(420,180,295,180).line(495,180,620,180)
f.text(480,285,'CH₃C(=O)OH: two oxygens, one acidic O—H',27)
f.text(480,340,'Three attached atoms; bond-order sum = 1 + 2 + 1 = 4',23)
save(u,9,'carbonyl-reading',f)
f=Figure('From measured mass to chemical amount')
f.rect(70,125,245,110).rect(360,125,245,110).rect(650,125,245,110)
f.text(192,175,'1.3812 g',30).text(192,210,'salicylic acid',23)
f.text(482,175,'÷ 138.12 g/mol',28).text(482,210,'molar mass',23)
f.text(772,175,'0.01000 mol',28).text(772,210,'10.00 mmol',23)
f.line(315,180,355,180,arrow=True,color=BLUE).line(605,180,645,180,arrow=True,color=BLUE)
f.text(480,330,'g ÷ (g/mol) = mol; the formula supplies the conversion factor.',24)
save(u,10,'mass-to-amount',f)

# Properties
u='properties'
f=Figure('A polar bond has unequal electron sharing')
f.text(220,175,'C',44).text(550,175,'O',44).line(255,160,510,160)
f.text(220,110,'δ+',30,color=RED).text(550,110,'δ−',30,color=BLUE)
f.line(290,230,495,230,arrow=True,color=BLUE).text(480,320,'Bond dipole points toward the more electronegative atom.',25)
save(u,1,'bond-polarity',f)
f=Figure('Bond dipoles combine as vectors')
f.text(240,175,'C',40).text(90,175,'O',38).text(390,175,'O',40)
bond(f,(120,160),(205,160),2);bond(f,(275,160),(355,160),2)
f.line(190,100,125,100,arrow=True,color=BLUE).line(290,100,355,100,arrow=True,color=BLUE)
f.text(240,290,'CO₂: opposite dipoles cancel',25)
f.text(725,110,'O',40).text(625,220,'H',32).text(825,220,'H',32)
f.line(700,130,650,185).line(750,130,800,185)
f.line(655,160,697,115,arrow=True,color=BLUE).line(795,160,753,115,arrow=True,color=BLUE)
f.text(725,290,'H₂O: bent; dipoles add',25)
save(u,2,'molecular-dipoles',f)
f=Figure('A hydrogen bond joins donor H and acceptor lone pair')
f.text(135,180,'R',34).text(225,180,'O',34).text(320,180,'H',34)
f.line(155,167,199,167).line(250,167,294,167)
f.text(600,180,'O',34).text(755,180,'R′',34).text(600,270,'H',32)
f.line(625,167,715,167).line(600,193,600,233)
f.text(563,175,':',28).line(345,165,545,165,color=BLUE,dash='8 7')
f.text(225,310,'donor O—H',26).text(650,310,'acceptor O lone pair',26)
f.text(480,355,'Dashed contact is an intermolecular interaction, not a covalent O—H bond.',22)
save(u,3,'hydrogen-bond',f)
f=Figure('An acceptor need not be a donor')
f.text(235,150,'CH₃CH₂—O—H',34).text(725,150,'CH₃—O—CH₃',34)
f.text(235,235,'O accepts; O—H donates',25).text(725,235,'O accepts; no O—H donor',25)
f.text(480,330,'Compare the bonds actually present, not only the presence of oxygen.',24)
save(u,4,'donor-acceptor',f)
f=Figure('Stronger liquid-phase attractions raise the boiling temperature','Illustrative normal boiling points in degrees Celsius; comparison of constitutional isomers.')
axes(f,'Compound','Normal boiling point (°C)')
for x,v,name in [(310,-24,'dimethyl ether'),(680,78,'ethanol')]:
 baseline=245;top=245-v*1.4
 if v<0: f.rect(x-65,baseline,130,abs(v)*1.4,fill='#a2c8d8')
 else:f.rect(x-65,top,130,v*1.4,fill='#23765d')
 f.text(x,top-15 if v>0 else 230,str(v),25).text(x,340,name,23)
f.line(140,245,835,245,color=BLUE,dash='5 5')
save(u,5,'boiling-comparison',f)
f=Figure('A polar head and a hydrocarbon chain compete in water')
for y,n in [(135,2),(260,6)]:
 points=[(170+i*75,y+(i%2)*35) for i in range(n)]
 zigzag(f,points);last=points[-1];f.line(last[0],last[1],last[0]+45,last[1]-20);f.text(last[0]+75,last[1]-18,'OH',28)
 f.text(775,y+20,'shorter chain' if n==2 else 'longer chain',25)
f.text(480,355,'For related alcohols, a larger nonpolar region usually reduces water solubility.',23)
save(u,6,'solubility-balance',f)
f=Figure('Vibrations change bond length; the molecule does not break')
for x,length in [(240,95),(700,155)]:
 f.circle(x-length/2,175,28,fill='#d8eee3').circle(x+length/2,175,35,fill='#e4f0f4')
 f.line(x-length/2+28,175,x+length/2-35,175,color=BLUE)
 f.text(x-length/2,185,'H',25).text(x+length/2,185,'O',25)
 f.text(x,280,'shorter separation' if length==95 else 'longer separation',25)
f.line(400,175,540,175,arrow=True,color=BLUE).line(540,205,400,205,arrow=True,color=BLUE)
f.text(480,355,'IR absorption excites molecular vibrations with a changing dipole.',24)
save(u,7,'bond-vibration',f)
f=Figure('An IR spectrum is evidence about functional groups','Simulated transmission trace illustrating broad alcohol O-H and carbonyl absorptions, not an experimental spectrum.')
axes(f,'Wavenumber (cm⁻¹), decreasing →','Transmittance (%, schematic)')
pts=[]
for i in range(401):
 wn=4000-i*3500/400;v=92-48*math.exp(-((wn-3350)/260)**2)-62*math.exp(-((wn-1720)/65)**2)-15*math.exp(-((wn-1150)/90)**2)
 pts.append((110+i*745/400,300-v*2.2))
f.path('M '+' L '.join(f'{x:.1f},{y:.1f}' for x,y in pts),color=BLUE)
for wn in [4000,3000,2000,1000,500]:f.text(110+(4000-wn)/3500*745,330,str(wn),19)
f.text(245,210,'broad O—H',23).text(600,255,'C=O',24).text(775,60,'SIMULATED',18,color=RED)
save(u,8,'ir-trace',f)
f=Figure('Two spectral clues distinguish an alcohol and a ketone')
for x,s,name in [(250,'CH₃CH₂OH','ethanol'),(710,'CH₃C(=O)CH₃','propanone')]:
 f.text(x,110,s,31).text(x,165,name,25)
 if x==250:
  f.path('M 100 230 Q 160 230 190 285 Q 230 335 275 250 Q 295 230 390 230',color=BLUE)
  f.text(x,355,'broad O—H; no C=O',24)
 else:
  f.path('M 560 230 L 680 230 L 710 310 L 740 230 L 860 230',color=GREEN)
  f.text(x,355,'strong C=O; no O—H',24)
save(u,9,'ir-candidate-comparison',f)
f=Figure('Ionisation changes interaction with water')
f.text(240,150,'CH₃COOH',36).text(720,150,'CH₃COO⁻',36)
for cx in [240,720]:
 for dx,dy in [(-110,70),(110,70),(0,-65)]:
  f.text(cx+dx,150+dy,'H₂O',24,color=BLUE).line(cx+dx*.75,150+dy*.75,cx+dx*.35,150+dy*.35,color=BLUE,dash='5 5')
f.text(240,305,'neutral acid',25).text(720,305,'ion–dipole hydration',25)
f.text(480,365,'Charge can change solubility; pH alone is not a complete solubility model.',23)
save(u,10,'hydration',f)

# Conformations
u='conformations'
f=Figure('Rotation changes conformation without changing connectivity')
f.text(245,145,'CH₃—CH₃',38).path('M 225 185 C 150 300 340 300 285 185',arrow=True)
f.text(690,145,'C—C sigma bond',30).text(690,220,'same connected atoms',27)
f.text(480,345,'A conformer is a different arrangement reached by bond rotation.',26)
save(u,1,'rotation',f)
f=Figure('A Newman projection looks along a carbon–carbon bond')
f.text(225,90,'view along C—C',27).circle(170,190,30,fill='#e4f0f4').circle(320,190,30,fill='#d8eee3').line(200,190,290,190)
f.line(65,190,125,190,arrow=True,color=BLUE).text(170,200,'C',28).text(320,200,'C',28)
newman(f,700,190,['H','H','H'],['H','H','H'],60,'front dot; rear circle')
f.text(230,320,'nearest carbon = front',25)
save(u,2,'newman-view',f)
f=Figure('Ethane has staggered and eclipsed conformations')
newman(f,250,185,['H','H','H'],['H','H','H'],60,'staggered: 60° offset')
# small rear offset visually separates labels while retaining eclipsed alignment
newman(f,710,185,['H','H','H'],['H','H','H'],0,'eclipsed: 0° offset')
# Hide coincident-label ambiguity by label front/rear colours in legend
f.text(480,38,'Blue: front bonds; dark: rear bonds',22)
save(u,3,'ethane-conformers',f)
f=Figure('The ethane torsional energy has three equivalent minima','Qualitative model only; barrier height is not measured in this drawing.')
axes(f,'Dihedral angle (degrees)','Relative energy (schematic)')
pts=[(110+t*745/360,260-130*(1+math.cos(math.radians(3*t)))/2) for t in range(361)]
f.path('M '+' L '.join(f'{x:.1f},{y:.1f}' for x,y in pts),color=BLUE,width=4)
for t in [0,60,120,180,240,300,360]:f.text(110+t*745/360,333,str(t),20)
f.text(500,80,'Eclipsed maxima; staggered minima',25)
save(u,4,'ethane-energy',f)
f=Figure('Butane: locate the two methyl groups')
newman(f,250,185,['CH₃','H','H'],['CH₃','H','H'],180,'anti: methyl groups 180° apart')
newman(f,710,185,['CH₃','H','H'],['CH₃','H','H'],60,'gauche: methyl groups 60° apart')
save(u,5,'butane-anti-gauche',f)
f=Figure('Butane conformers differ in energy','Qualitative torsional profile; anti is lowest and methyl-methyl eclipsing highest.')
axes(f,'C1–C2–C3–C4 dihedral angle (°)','Relative energy (schematic)')
vals=[(0,220),(60,65),(120,175),(180,0),(240,175),(300,65),(360,220)]
pts=[(110+a*745/360,285-e*.8) for a,e in vals]
f.path('M '+' L '.join(f'{x:.1f},{y:.1f}' for x,y in pts),color=BLUE,width=4)
for a,e in vals:f.text(110+a*745/360,330,str(a),20)
f.text(480,255,'anti',24).text(240,180,'gauche',22)
save(u,6,'butane-energy',f)
f=Figure('Rings restrict rotation and can store strain')
for cx,n,name in [(235,3,'cyclopropane'),(710,6,'cyclohexane (flat sketch)')]:
 pts=[(cx+105*math.cos(2*math.pi*j/n-math.pi/2),185+100*math.sin(2*math.pi*j/n-math.pi/2)) for j in range(n)]
 for a,b in zip(pts,pts[1:]+pts[:1]):bond(f,a,b)
 f.text(cx,335,name,25)
f.text(480,45,'A flat ring drawing does not establish the actual 3D shape.',24)
save(u,7,'ring-strain',f)
def chair(f,cx,cy,scale=1,flipped=False):
    # Projection of a D3d six-carbon chair: alternating z = +/-sqrt(1/32).
    pts=[]
    for j in range(6):
        a=j*math.pi/3
        z=(1 if j%2==0 else -1)*math.sqrt(1/32)*(-1 if flipped else 1)
        pts.append((cx+160*scale*math.cos(a),cy+160*scale*(.4*math.sin(a)-.92*z)))
    for a,b in zip(pts,pts[1:]+pts[:1]):bond(f,a,b)
    return pts
f=Figure('A chair is a nonplanar cyclohexane conformation')
pts=chair(f,340,205)
x,y=pts[0]
f.line(x,y,x,y-100,color=BLUE).text(x,y-120,'axial up',24,color=BLUE)
f.line(x,y,x+90,y+30,color=GREEN).text(x+140,y+65,'equatorial down',23,color=GREEN)
f.text(480,325,'At one carbon: one axial and one equatorial position.',24)
f.text(480,375,'Axial positions alternate up/down around the ring.',24)
save(u,8,'chair-directions',f)
f=Figure('A ring flip exchanges axial and equatorial, preserving up/down')
a=chair(f,235,220,.72,False)
b=chair(f,700,220,.72,True)
x,y=a[0];f.line(x,y,x,y-85,color=BLUE).text(x,y-105,'CH₃ (up)',24)
x,y=b[0];f.line(x,y,x+65,y-25,color=BLUE).text(x+60,y-50,'CH₃ (up)',24)
f.text(235,325,'axial-up methyl',25).text(700,325,'equatorial-up methyl',25)
f.line(415,200,505,200,arrow=True,color=BLUE).line(505,225,415,225,arrow=True,color=BLUE)
f.text(480,375,'Both ring geometry and local bond direction change; up/down is retained.',22)
save(u,9,'ring-flip',f)
f=Figure('A free-energy gap predicts unequal conformer populations','Two-state teaching model at 298.15 kelvin; supplied free energies include degeneracy.')
axes(f,'Conformer','Predicted population')
w=math.exp(-3500/(8.314462618*298.15));p=[1/(1+w),w/(1+w)]
for x,v,name in [(310,p[0],'A: ΔG = 0'),(680,p[1],'B: ΔG = 3500 J/mol')]:
 f.rect(x-70,300-v*225,140,v*225,fill='#176b91').text(x,280-v*225,f'{100*v:.1f}%',27).text(x,335,name,22)
save(u,10,'conformer-populations',f)

# Stereochemistry
u='stereochemistry'
f=Figure('Classify the difference before assigning a descriptor')
f.text(230,65,'CH₃CH₂OH',32).text(710,65,'CH₃OCH₃',32)
f.line(365,55,565,55,arrow=True,color=BLUE).text(480,25,'different connectivity',20)
f.parts.append('<g transform="scale(.75)">')
newman(f,333.33,300,['CH₃','H','H'],['CH₃','H','H'],180,'')
newman(f,953.33,300,['CH₃','H','H'],['CH₃','H','H'],60,'')
f.parts.append('</g>')
f.text(480,390,'Same connected atoms below: anti and gauche butane.',23)
save(u,1,'classification',f)
f=Figure('A tetrahedral stereocentre has four different substituents')
f.text(480,200,'C',35).line(465,170,400,95).line(500,175,580,105)
f.text(380,80,'COOH',29).text(610,95,'CH₃',29)
f.wedge(465,210,365,290).wedge(495,210,590,290,hashed=True)
f.text(330,325,'OH',30).text(620,325,'H',30)
f.text(480,50,'lactic acid: CH₃CH(OH)COOH',28)
save(u,2,'stereocentre',f)
f=Figure('Mirror images need not be superimposable')
fischer(f,245,175,left='OH',right='H',label='one lactic-acid configuration')
fischer(f,715,175,left='H',right='OH',label='its mirror image')
f.line(480,55,480,330,color=BLUE,dash='8 7').text(480,370,'A mirror swaps left/right; it does not change connectivity.',24)
save(u,3,'enantiomers',f)
f=Figure('CIP priority begins at the directly attached atom')
f.text(480,160,'OH  >  COOH  >  CH₃  >  H',36)
for x,n,reason in [(175,'1','O beats C'),(400,'2','C attached to O,O,O'),(625,'3','C attached to H,H,H'),(835,'4','H is lowest')]:
 f.text(x,95,n,33,color=BLUE).text(x,255,reason,21)
f.text(480,340,'A double bond is counted using duplicate atoms for the comparison.',23)
save(u,4,'cip-priority',f)
f=Figure('Fischer projections make an explicit depth convention')
fischer(f,250,175,left='OH',right='H')
f.text(700,120,'Horizontal bonds: toward you',27).text(700,200,'Vertical bonds: away from you',27)
f.text(700,280,'The central crossing is carbon.',25)
f.text(480,375,'Rotating a Fischer projection 90° changes its implied configuration.',23)
save(u,5,'fischer-depth',f)
f=Figure('Assign S to this Fischer projection of lactic acid')
fischer(f,245,175,left='OH (1)',right='H (4)',top='COOH (2)',bottom='CH₃ (3)')
f.text(705,115,'1 → 2 → 3 appears clockwise',25).text(705,180,'4 points toward you: invert',25).text(705,260,'Configuration = S',34,color=BLUE)
f.text(480,375,'With 4 pointing away, clockwise would mean R.',25)
save(u,6,'rs-worked',f)
f=Figure('A double bond can lock two configurations')
for cx,mode in [(245,'Z'),(715,'E')]:
 bond(f,(cx-55,180),(cx+55,180),2)
 f.line(cx-55,180,cx-125,110).line(cx-55,180,cx-125,250)
 f.line(cx+55,180,cx+125,110).line(cx+55,180,cx+125,250)
 f.text(cx-150,95,'CH₃',25).text(cx-150,285,'H',25)
 f.text(cx+150,95,'CH₃' if mode=='Z' else 'H',25).text(cx+150,285,'H' if mode=='Z' else 'CH₃',25)
 f.text(cx,345,mode+'-but-2-ene',27)
f.text(480,45,'On each C: CH₃ has higher priority than H.',25)
save(u,7,'ez-alkenes',f)
f=Figure('A meso structure is achiral despite two stereocentres')
x=290
f.line(x,85,x,295).line(x-70,145,x+70,145).line(x-70,235,x+70,235)
f.text(x,65,'COOH',27).text(x,330,'COOH',27)
for y in [145,235]:f.text(x-105,y+8,'H',26).text(x+105,y+8,'OH',26)
f.line(100,190,485,190,color=BLUE,dash='7 6')
f.text(715,130,'meso-tartaric acid',28).text(715,200,'two stereocentres',26).text(715,265,'internal symmetry → achiral',25)
f.text(480,375,'Number of stereocentres alone does not establish chirality.',24)
save(u,8,'meso',f)
f=Figure('A mixture can have an enantiomeric excess')
for i in range(10):
 x=125+i*78;f.circle(x,165,24,color=BLUE if i<7 else GREEN,fill='#e4f0f4' if i<7 else '#d8eee3').text(x,173,'R' if i<7 else 'S',23)
f.text(480,260,'70% R + 30% S → 40% ee in favour of R',31)
f.text(480,330,'Major percentage = (100 + ee)/2 for a two-enantiomer mixture.',24)
save(u,9,'enantiomeric-excess',f)
f=Figure('Configuration and optical rotation answer different questions')
f.text(250,110,'R / S',43).text(250,185,'arrangement of substituents',25)
f.line(510,70,510,310,color=BLUE,dash='8 8')
f.text(725,110,'(+) / (−)',43).text(725,185,'measured optical rotation',25)
f.line(630,260,820,260,arrow=True,color=BLUE).path('M 660 290 Q 725 220 790 290',color=GREEN)
f.text(480,360,'R does not mean (+); S does not mean (−).',29)
save(u,10,'rotation-vs-configuration',f)
print('Built 40 foundation SVG diagrams.')
