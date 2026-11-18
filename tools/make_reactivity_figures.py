"""Original molecular/mechanistic diagrams for the five reactivity units."""
from pathlib import Path
import math
from chemfig import Figure, BLUE, RED, GREEN, INK

ROOT=Path(__file__).resolve().parents[1]

def save(unit,number,name,f):
    # Long captions are two readable lines rather than text outside the canvas.
    import re, textwrap
    from html import unescape, escape
    for i,part in enumerate(f.parts):
        match=re.search(r'<text x="480" y="([\d.]+)".*?font-size="([\d.]+)".*?>(.*?)</text>',part)
        if match and float(match[1]) >= 330 and len(unescape(match[3])) > 78:
            lines=textwrap.wrap(unescape(match[3]),width=76)
            f.parts[i]='\n'.join(re.sub(r'y="[\d.]+"',f'y="{345+25*j}"',part,count=1).replace(match[3],escape(line)) for j,line in enumerate(lines))
    f.save(ROOT/unit/'figures'/f'{number:02d}-{name}.svg')

def carbonyl(f,x,y,left='CH₃',right='H',charge=False):
    f.text(x,y+9,'C',31)
    f.line(x-7,y-20,x-7,y-78).line(x+7,y-20,x+7,y-78)
    f.text(x,y-95,'O',32,color=RED)
    f.line(x-22,y,x-95,y+27).text(x-135,y+42,left,28)
    f.line(x+22,y,x+112,y+27).text(x+135,y+42,right,28)
    if charge:
        f.text(x+43,y-25,'δ+',23,color=BLUE).text(x+40,y-100,'δ−',23,color=RED)

def tetra(f,x,y,left='CH₃',right='H',top='O⁻',bottom='OH'):
    f.text(x,y+9,'C',30)
    for dx,dy in [(0,-65),(-105,20),(120,20),(0,85)]: f.line(x+dx*.2,y+dy*.2,x+dx,y+dy)
    f.text(x,y-85,top,29,color=RED)
    f.text(x-140,y+35,left,27).text(x+140,y+35,right,27)
    f.text(x,y+115,bottom,28,color=BLUE)

def scheme(title,left,right,condition='',note='',size=29):
    f=Figure(title)
    f.text(230,190,left,size).text(730,190,right,size)
    f.line(410,180,550,180,color=BLUE,arrow=True)
    f.text(480,125,condition,23).text(480,330,note,22)
    return f

def axes(f,ylab='Free energy',xlab='Reaction coordinate'):
    f.line(100,315,870,315,arrow=True).line(100,315,100,60,arrow=True)
    f.text(200,45,ylab,23).text(490,375,xlab,24)

# Mechanisms: electron flow is shown only where the drawing specifies a pair.
u='mechanisms'
f=Figure('Lone-pair donation and bond cleavage in proton transfer')
f.text(200,180,':OH⁻',34).text(510,180,'H',34).text(745,180,'OCH₃',34)
f.line(535,170,655,170)
f.path('M 170 153 Q 330 50 496 153',arrow=True)
f.path('M 595 170 Q 650 230 700 190',arrow=True)
f.text(480,305,'CH₃OH + OH⁻ ⇌ CH₃O⁻ + H₂O',29)
f.text(480,355,'Both arrows carry a pair; total charge remains −1.',23)
save(u,1,'pair-flow',f)

f=Figure('Heterolysis gives both bonding electrons to one atom')
f.text(235,165,'CH₃',32).text(410,165,'Br',32)
f.line(280,155,370,155)
f.path('M 325 155 Q 375 70 410 120',arrow=True)
f.text(700,165,'CH₃⁺ + Br⁻',31)
f.line(470,155,555,155,arrow=True,color=BLUE)
f.text(480,285,'A bookkeeping example of heterolysis',28)
f.text(480,335,'This drawing does not claim methyl bromide readily forms CH₃⁺.',21)
save(u,2,'heterolysis',f)

f=Figure('A nucleophile supplies a pair; an electrophile accepts it')
carbonyl(f,690,210,charge=True)
f.text(200,230,':OH⁻',34,color=BLUE)
f.path('M 160 195 Q 420 65 665 190',arrow=True)
f.text(200,315,'electron-pair donor',24,color=BLUE)
f.text(690,350,'electron-poor carbon',24)
save(u,3,'source-sink',f)

f=Figure('Carbonyl addition requires two linked arrows')
carbonyl(f,380,205,charge=True)
f.text(100,220,':OH⁻',30,color=BLUE)
f.path('M 72 185 Q 200 105 355 188',arrow=True)
f.path('M 388 155 Q 445 110 410 92',arrow=True)
tetra(f,745,190)
f.text(480,375,'Form C—O while moving the C=O π pair onto oxygen.',24)
save(u,4,'carbonyl-pairs',f)

f=Figure('Substitution forms one bond while another breaks')
f.text(140,195,':OH⁻',31,color=BLUE).text(425,195,'CH₃',32).text(690,195,'Br',32)
f.line(475,185,650,185)
f.path('M 105 165 Q 245 65 390 166',arrow=True)
f.path('M 565 185 Q 660 90 689 145',arrow=True)
f.text(480,315,'HO—CH₃ + Br⁻',34)
f.text(480,365,'A concerted SN2 step: two arrows, no carbocation intermediate.',23)
save(u,5,'concerted-substitution',f)

f=Figure('A proton transfer can neutralise an alkoxide')
f.text(160,170,':O⁻—CH₂CH₃',29).text(540,170,'H',31).text(740,170,'O—H',31)
f.line(560,160,690,160)
f.path('M 75 140 Q 280 25 525 140',arrow=True)
f.path('M 625 160 Q 690 230 715 188',arrow=True)
f.text(480,300,'CH₃CH₂O⁻ + H₂O ⇌ CH₃CH₂OH + OH⁻',28)
f.text(480,355,'Write a real proton donor; do not add an unexplained H⁺.',23)
save(u,6,'proton-donor',f)

f=Figure('A one-step energy profile has one transition state')
axes(f)
f.path('M 125 255 C 250 255 285 95 430 95 C 575 95 640 270 845 270',width=4)
f.text(430,80,'‡',35).text(210,235,'reactants',22).text(745,245,'products',22)
f.line(125,255,430,255,dash='6 5',width=2)
f.line(430,255,430,95,arrow=True,color=RED)
f.text(480,185,'ΔG‡',25,color=RED)
save(u,7,'one-step-energy',f)

f=Figure('A multistep path contains an intermediate minimum')
axes(f)
f.path('M 125 270 C 200 270 230 85 315 85 C 370 85 380 205 455 205 C 515 205 530 135 600 135 C 700 135 715 280 845 280',width=4)
f.text(315,70,'TS₁',23).text(600,115,'TS₂',23)
f.text(455,245,'intermediate',23)
f.text(770,65,'Schematic; not to scale',20)
save(u,8,'intermediate-energy',f)

f=Figure('A catalyst changes the path, not the endpoint difference')
axes(f)
f.path('M 125 265 C 260 265 300 85 460 85 C 620 85 650 285 845 285',color=RED,width=4)
f.path('M 125 265 C 260 265 300 175 460 175 C 620 175 650 285 845 285',color=GREEN,width=4)
f.text(690,80,'uncatalysed',23,color=RED).text(700,165,'catalysed',23,color=GREEN)
save(u,9,'catalyst-profile',f)

f=Figure('Check atom identity and total charge after a step')
f.text(210,130,'HO⁻ + CH₃Br',32).text(735,130,'CH₃OH + Br⁻',32)
f.line(420,120,540,120,arrow=True,color=BLUE)
for x in [240,720]:
 f.text(x,230,'C₁  H₄  O₁  Br₁',29).text(x,300,'total charge −1',27)
f.text(480,365,'Electrons move; atoms are neither created nor silently removed.',23)
save(u,10,'bookkeeping',f)

# Substitution.
u='substitution'
f=Figure('Substrate classification counts carbon neighbours of C—Br')
for x,formula,kind in [(170,'CH₃CH₂—Br','primary'),(480,'(CH₃)₂CH—Br','secondary'),(790,'(CH₃)₃C—Br','tertiary')]:
 f.text(x,160,formula,28).circle(x+70,150,62,color=BLUE,fill='none')
 f.text(x,285,kind,28)
f.text(480,355,'Count carbon groups bonded directly to the carbon bearing Br.',22)
save(u,1,'substrates',f)

f=Figure('A methyl SN2 example follows a concerted pathway')
f.text(135,190,':OH⁻',32).text(420,190,'CH₃',32).text(690,190,'Br',32)
f.line(470,180,650,180)
f.path('M 100 160 Q 240 65 385 165',arrow=True)
f.path('M 560 180 Q 645 85 688 140',arrow=True)
f.text(480,305,'CH₃OH + Br⁻',35).text(480,360,'Backside approach opposite the leaving group.',24)
save(u,2,'sn2-step',f)

f=Figure('SN2 inverts the arrangement at the reacting carbon')
for x,invert in [(260,False),(710,True)]:
 f.text(x,195,'C',31)
 f.line(x-20,180,x-95,140).text(x-135,130,'CH₃',27)
 f.line(x-20,200,x-95,250).text(x-135,278,'H',27)
 f.wedge(x+12,180,x+80,105,hashed=invert).text(x+105,85,'R',29)
 f.wedge(x+12,205,x+80,280,hashed=not invert).text(x+125,310,'Br' if not invert else 'OH',28)
f.line(420,185,535,185,arrow=True,color=BLUE)
f.text(480,365,'Follow the three retained groups; R/S labels also depend on priorities.',22)
save(u,3,'inversion',f)

f=Figure('A simple SN2 rate model depends on both reactants')
axes(f,'Relative initial rate','Relative [OH⁻] at fixed [CH₃Br]')
f.line(125,300,800,85,color=BLUE,width=4)
for j in range(1,5):
 x=125+j*150;y=300-j*48
 f.circle(x,y,7,fill=BLUE,color=BLUE)
f.text(615,85,'rate = k[CH₃Br][OH⁻]',27)
save(u,4,'sn2-rate',f)

f=Figure('Steric crowding makes backside access harder')
for cx,groups,label in [(220,1,'less crowded'),(710,3,'more crowded')]:
 f.circle(cx,200,22,fill='white').text(cx,209,'C',26)
 f.line(cx+25,200,cx+100,200).text(cx+130,209,'Br',28)
 for angle in ([135] if groups==1 else [100,180,260]):
  a=math.radians(angle);x=cx+80*math.cos(a);y=200+80*math.sin(a)
  f.line(cx+25*math.cos(a),200+25*math.sin(a),x,y).circle(x,y,38,fill='#dbe8ed')
  f.text(x,y+8,'CH₃',22)
 f.text(cx,340,label,26)
f.line(50,200,155,200,arrow=True,color=BLUE)
f.line(470,200,595,200,arrow=True,color=BLUE)
save(u,5,'steric-access',f)

f=Figure('SN1 starts with ionisation to a carbocation and leaving group')
f.text(235,175,'(CH₃)₃C',32).text(490,175,'Br',32)
f.line(330,165,450,165)
f.path('M 390 165 Q 450 75 490 130',arrow=True)
f.text(480,300,'(CH₃)₃C⁺ + Br⁻',34)
f.text(480,355,'Suitable polar solvent stabilises separated charge; ionisation costs energy.',21)
save(u,6,'sn1-ionisation',f)

f=Figure('Water captures the carbocation through an oxygen lone pair')
f.text(235,190,'H₂O:',34,color=BLUE).text(640,190,'C⁺(CH₃)₃',34)
f.path('M 275 175 Q 430 65 580 160',arrow=True)
f.text(480,305,'(CH₃)₃C—OH₂⁺',34)
f.text(480,355,'The attacking water oxygen now has three bonds and charge +1.',23)
save(u,7,'sn1-capture',f)

f=Figure('A second water molecule removes a proton from oxonium')
f.text(140,175,'H₂O:',31).text(470,175,'H',31).text(700,175,'O⁺H—C(CH₃)₃',29)
f.line(490,165,575,165)
f.path('M 180 155 Q 300 55 455 150',arrow=True)
f.path('M 530 165 Q 570 240 605 195',arrow=True)
f.text(480,300,'(CH₃)₃COH + H₃O⁺',33)
f.text(480,355,'Proton transfer restores a neutral alcohol.',25)
save(u,8,'sn1-deprotonation',f)

f=Figure('A planar carbocation permits approach from either face')
f.path('M 240 170 L 610 110 L 760 250 L 390 310 Z',color='#9fb9c6',fill='#e4eef2',width=2)
f.text(480,215,'C⁺',34)
for ex,ey,tx,ty,label in [(370,245,335,263,'R¹'),(595,230,635,247,'R²'),(440,155,420,140,'R³')]:
 f.line(480+(ex-480)*.2,205+(ey-205)*.2,ex,ey)
 f.text(tx,ty,label,27)
f.line(480,45,480,177,arrow=True,color=BLUE,dash='9 6')
f.line(480,350,480,238,arrow=True,color=BLUE,dash='9 6')
f.text(535,65,'Nu: above face',24,anchor='start').text(620,335,'Nu: below face',24)
f.text(170,155,'plane of',22).text(170,190,'three C bonds',22)
f.text(480,390,'Oblique schematic: approach is from either side of the plane, not along it.',21)
save(u,9,'planar-cation',f)

f=Figure('The same substrate may have competing reaction pathways')
f.text(480,80,'Secondary R—Br',32)
f.line(440,110,220,220,arrow=True,color=BLUE).line(520,110,740,220,arrow=True,color=BLUE)
f.text(200,270,'Substitution: R—Nu',29).text(750,270,'Elimination: C=C',29)
f.text(200,325,'nucleophilic attack at C',22).text(750,325,'base removes β-H',22)
f.text(480,380,'Substrate, nucleophile/base, solvent and temperature all matter.',22)
save(u,10,'competition',f)

# Elimination and addition.
u='elimination-addition'
f=Figure('Elimination loses groups from neighbouring carbons')
f.text(250,175,'CH₃—CH₂—Br',33)
f.text(720,175,'CH₂=CH₂',34)
f.line(445,165,540,165,arrow=True,color=BLUE)
f.text(480,120,'base',26)
f.text(480,280,'CH₃CH₂Br + OH⁻ → CH₂=CH₂ + H₂O + Br⁻',28)
f.text(480,345,'A formal E2 example; primary substrates can strongly favour SN2.',22)
save(u,1,'elimination-balance',f)

f=Figure('An E2 step moves three electron pairs together')
f.text(90,225,':B⁻',31).text(300,100,'H',31)
f.text(400,230,'C',31).text(650,230,'C',31).text(840,135,'Br',31)
f.line(315,115,385,200).line(425,220,620,220).line(670,200,810,140)
f.text(400,300,'β-carbon',23).text(650,300,'α-carbon',23)
f.path('M 65 190 Q 160 55 280 90',arrow=True)
f.path('M 350 155 Q 445 110 510 205',arrow=True)
f.path('M 745 170 Q 855 220 860 160',arrow=True)
f.text(480,375,'B—H forms; C—H → C=C; C—Br → Br⁻. Other C substituents omitted.',21)
save(u,2,'e2-arrows',f)

f=Figure('Anti-periplanar means the reacting C—H and C—Br bonds oppose')
f.circle(480,200,65,color=INK,fill='none').circle(480,200,8,fill=INK)
for ang,label in [(-90,'H'),(30,'R'),(150,'H')]:
 a=math.radians(ang);f.line(480,200,480+125*math.cos(a),200+125*math.sin(a));f.text(480+155*math.cos(a),208+155*math.sin(a),label,28)
for ang,label in [(90,'Br'),(210,'R'),(330,'H')]:
 a=math.radians(ang);f.line(480+65*math.cos(a),200+65*math.sin(a),480+125*math.cos(a),200+125*math.sin(a));f.text(480+155*math.cos(a),208+155*math.sin(a),label,28)
f.text(165,190,'Look down',25).text(165,230,'the C—C bond',25)
f.text(790,190,'H / Br dihedral',23).text(790,235,'180°',33)
save(u,3,'anti-newman',f)

f=Figure('Rotation changes which beta hydrogen is aligned')
for cx,angle,lab in [(250,90,'anti: 180°'),(710,-90,'eclipsed: 0°')]:
 f.circle(cx,190,55,fill='none').circle(cx,190,7,fill=INK)
 f.line(cx,190,cx,75)
 a=math.radians(angle);f.line(cx+55*math.cos(a),190+55*math.sin(a),cx+115*math.cos(a),190+115*math.sin(a))
 if angle==90:
  f.text(cx,60,'H',29).text(cx,340,'Br',29)
 else:
  f.text(cx,60,'H / Br',29)
 f.text(cx,385,lab,25)
f.text(480,35,'Partial Newman views: other substituents omitted',22)
save(u,4,'rotation-geometry',f)

f=Figure('A secondary bromide can give more than one alkene')
f.text(480,65,'CH₃—CH(Br)—CH₂—CH₃',32)
f.line(420,100,220,210,arrow=True,color=BLUE).line(540,100,740,210,arrow=True,color=BLUE)
f.text(220,270,'CH₂=CHCH₂CH₃',29).text(740,270,'CH₃CH=CHCH₃',29)
f.text(220,325,'but-1-ene',25).text(740,325,'but-2-ene (E and Z)',25)
f.text(480,380,'Product proportions depend on base, geometry and conditions.',23)
save(u,5,'alkene-products',f)

f=Figure('HBr addition begins with protonation of the alkene')
f.text(270,225,'CH₃—CH=CH₂',31).text(680,115,'H',31).text(850,115,'Br',31)
f.line(705,105,810,105)
f.path('M 330 210 Q 485 65 660 95',arrow=True)
f.path('M 760 105 Q 815 45 850 75',arrow=True)
f.text(480,320,'CH₃—CH⁺—CH₃ + Br⁻',32)
f.text(480,375,'Representative ionic pathway without peroxides.',23)
save(u,6,'hbr-protonation',f)

f=Figure('Bromide captures the carbocation')
f.text(210,180,':Br⁻',34,color=BLUE).text(650,180,'CH₃—CH⁺—CH₃',32)
f.path('M 175 148 Q 390 50 650 145',arrow=True)
f.text(480,305,'CH₃—CH(Br)—CH₃',34)
f.text(480,365,'2-bromopropane: explain the route using relative carbocation stability.',22)
save(u,7,'hbr-capture',f)

f=Figure('A bromonium bridge directs anti addition')
f.text(420,200,'C',31).text(620,200,'C',31).line(445,190,595,190)
f.text(520,80,'Br⁺',31).line(430,165,500,95).line(610,165,540,95)
f.text(520,310,':Br⁻',31,color=BLUE)
f.path('M 485 280 Q 335 250 405 208',arrow=True)
f.path('M 465 132 Q 390 60 480 55',arrow=True)
f.text(480,375,'Backside opening gives anti addition; alkene substituents omitted.',22)
save(u,8,'bromonium',f)

f=Figure('E1 shares a carbocation with a possible SN1 route')
f.text(480,80,'(CH₃)₃C⁺',34)
f.line(440,115,220,230,arrow=True,color=BLUE).line(520,115,740,230,arrow=True,color=BLUE)
f.text(205,275,'(CH₃)₃COH',31).text(750,275,'(CH₃)₂C=CH₂',31)
f.text(205,325,'water capture → SN1',23).text(750,325,'β-H loss → E1',23)
f.text(480,380,'Ionisation does not by itself identify the final product.',23)
save(u,9,'e1-competition',f)

f=Figure('Addition and elimination change the local bond inventory')
f.text(210,170,'C—C',42).text(740,170,'C=C',42)
f.line(360,130,590,130,arrow=True,color=BLUE).text(480,90,'elimination',25)
f.line(590,220,360,220,arrow=True,color=GREEN).text(480,270,'addition',25)
f.text(480,350,'Track which σ bonds break/form and which π bond appears/disappears.',22)
save(u,10,'bond-inventory',f)

# Alcohol redox.
u='alcohols-redox'
f=Figure('A primary alcohol can be oxidised in two stages')
f.text(170,180,'CH₃CH₂OH',30).text(480,180,'CH₃CHO',30).text(790,180,'CH₃COOH',30)
f.line(290,170,365,170,arrow=True,color=BLUE).line(590,170,665,170,arrow=True,color=BLUE)
f.text(330,120,'[O]',25).text(630,120,'[O]',25)
f.text(170,280,'ethanol',25).text(480,280,'ethanal',25).text(790,280,'ethanoic acid',25)
f.text(480,360,'[O] summarises oxidation; actual reagents and conditions must be specified.',21)
save(u,1,'primary-map',f)

f=Figure('A secondary alcohol gives a ketone')
f.text(240,170,'CH₃CH(OH)CH₃',32).text(730,170,'CH₃C(=O)CH₃',32)
f.line(425,160,555,160,arrow=True,color=BLUE).text(480,110,'oxidation',25)
f.text(240,280,'propan-2-ol',27).text(730,280,'propanone',27)
f.text(480,350,'The carbon skeleton is retained in this selected transformation.',23)
save(u,2,'secondary-map',f)

f=Figure('The carbon bearing OH needs a hydrogen for simple alcohol oxidation')
for x,label,groups in [(240,'secondary',('CH₃','CH₃','H')),(720,'tertiary',('CH₃','CH₃','CH₃'))]:
 tetra(f,x,170,left=groups[0],right=groups[1],top='OH',bottom=groups[2])
 f.text(x,340,label,26)
save(u,3,'tertiary-boundary',f)

f=Figure('Oxidation state counts C—O and C—H bonds differently')
for x,formula,value in [(180,'CH₃CH₂OH','−1'),(480,'CH₃CHO','+1'),(780,'CH₃COOH','+3')]:
 f.text(x,145,formula,30).text(x,235,value,42,color=BLUE)
f.text(480,320,'For the changing carbon: each C—O bond order +1; each C—H −1.',22)
f.text(480,365,'C—C bonds contribute zero; oxidation raises the assigned value.',22)
save(u,4,'oxidation-count',f)

f=Figure('A hydride equivalent adds to carbonyl carbon')
carbonyl(f,330,200,left='CH₃',right='CH₃')
f.text(90,205,'H:⁻',31,color=BLUE)
f.path('M 115 175 Q 220 80 305 178',arrow=True)
f.path('M 338 145 Q 405 90 360 83',arrow=True)
tetra(f,745,190,left='CH₃',right='CH₃',bottom='H')
f.text(480,375,'H:⁻ represents delivery from a reagent such as BH₄⁻, not free hydride in water.',20)
save(u,5,'hydride-addition',f)

f=scheme('A proton donor converts the alkoxide to alcohol','(CH₃)₂CHO⁻ + H₂O','(CH₃)₂CHOH + OH⁻','proton transfer','A separate proton-transfer step follows carbonyl reduction.',25)
save(u,6,'workup',f)

f=Figure('Protonation converts a poor leaving group into water')
f.text(230,170,'CH₃CH₂—OH',32).text(730,170,'CH₃CH₂—OH₂⁺',32)
f.line(430,160,540,160,arrow=True,color=BLUE).text(480,110,'H⁺ transfer',25)
f.text(230,280,'OH⁻: poor leaving group',24).text(730,280,'H₂O: better leaving group',24)
f.text(480,350,'Protonation enables selected reactions; it does not guarantee one pathway.',21)
save(u,7,'activate-oh',f)

f=Figure('A synthetic IR comparison tests a functional-group change','Simulated teaching traces, not experimental spectra.')
f.line(100,310,870,310,arrow=True).line(100,310,100,55,arrow=True)
for shift,col,carbonylpeak in [(0,BLUE,False),(80,RED,True)]:
 pts=[]
 for j in range(400):
  wn=4000-j*3500/399;x=120+j*700/399
  val=30*math.exp(-((wn-3350)/250)**2) if not carbonylpeak else 0
  val+=60*math.exp(-((wn-1720)/45)**2) if carbonylpeak else 0
  pts.append((x,105+shift+val))
 f.path('M '+' L '.join(f'{x:.1f},{y:.1f}' for x,y in pts),color=col,width=3)
f.text(230,70,'alcohol',23,color=BLUE).text(710,185,'carbonyl product',23,color=RED)
for x,t in [(120,'4000'),(470,'2250'),(820,'500')]:f.text(x,340,t,20)
f.text(480,380,'Wavenumber / cm⁻¹ · Simulated, vertically offset traces',22)
f.text(210,35,'Transmittance (arbitrary)',20)
save(u,8,'ir-redox',f)

f=Figure('Choose the transformation before choosing a reagent')
f.text(190,130,'CH₃CH₂OH',30).text(750,130,'CH₃CHO',30)
f.line(340,120,600,120,arrow=True,color=BLUE).text(480,80,'controlled oxidation',24)
f.text(190,270,'CH₃CHO',30).text(750,270,'CH₃CH₂OH',30)
f.line(340,260,600,260,arrow=True,color=GREEN).text(480,220,'NaBH₄, then proton donor',23)
f.text(480,365,'Preserving an aldehyde requires conditions that limit further oxidation.',22)
save(u,9,'route-choice',f)

f=Figure('Theoretical yield is capped by stoichiometric availability')
for x,height,lab in [(220,140,'10.0 mmol acid'),(480,210,'15.0 mmol anhydride'),(740,140,'10.0 mmol product max')]:
 f.rect(x-60,300-height,120,height,fill='#2c7894')
 f.text(x,340,lab,21)
f.line(100,300,870,300)
f.text(480,65,'Supplied balanced 1:1 teaching reaction',27)
f.text(480,390,'Amount divided by stoichiometric coefficient identifies the limiting reagent.',21)
save(u,10,'limiting-reagent',f)

# Carbonyl addition.
u='carbonyl-addition'
f=Figure('The carbonyl group is polar and approximately planar')
carbonyl(f,480,215,left='CH₃',right='H',charge=True)
f.path('M 375 242 Q 480 325 585 242',color=GREEN)
f.text(480,345,'Three bond directions around carbon: approximately 120°',23)
save(u,1,'carbonyl-geometry',f)

f=Figure('Hydroxide adds through a lone pair; the pi pair moves to oxygen')
carbonyl(f,345,205)
f.text(85,220,':OH⁻',30,color=BLUE)
f.path('M 55 185 Q 195 85 320 182',arrow=True)
f.path('M 353 150 Q 415 85 375 87',arrow=True)
tetra(f,745,190)
f.text(480,375,'Ethanal + OH⁻ ⇌ CH₃CH(OH)O⁻',28)
save(u,2,'hydrate-addition',f)

f=Figure('Water protonates the alkoxide and regenerates hydroxide')
f.text(170,175,'CH₃CH(OH)O:⁻',28).text(515,175,'H',30).text(740,175,'OH',30)
f.line(540,165,705,165)
f.path('M 300 150 Q 390 65 500 145',arrow=True)
f.path('M 620 165 Q 710 235 730 195',arrow=True)
f.text(480,300,'CH₃CH(OH)₂ + OH⁻',33)
f.text(480,355,'Hydration is reversible; the equilibrium depends on the carbonyl structure.',22)
save(u,3,'hydrate-protonation',f)

f=Figure('Acid catalysis begins by protonating carbonyl oxygen')
carbonyl(f,295,220)
f.text(295,100,':',28,color=RED).text(575,105,'H',31).text(780,105,'OH₂⁺',31)
f.line(595,95,720,95)
f.path('M 290 90 Q 420 20 560 80',arrow=True)
f.path('M 650 95 Q 710 150 742 127',arrow=True)
f.text(670,260,'CH₃CH=OH⁺ + H₂O',30)
f.text(480,365,'Oxygen accepts H⁺; the carbonyl carbon becomes more electrophilic.',22)
save(u,4,'acid-activation',f)

f=Figure('Methanol attacks an activated carbonyl')
f.text(110,225,'CH₃OH:',29,color=BLUE)
carbonyl(f,425,210)
f.rect(375,80,100,44,fill='#f7fafc').text(425,115,'O⁺H',30,color=RED)
f.path('M 180 195 Q 290 75 400 187',arrow=True)
f.path('M 433 155 Q 490 105 454 95',arrow=True)
tetra(f,745,190,top='OH',bottom='O⁺HCH₃')
f.text(480,375,'The incoming oxygen has three bonds before deprotonation.',23)
save(u,5,'methanol-attack',f)

f=Figure('Deprotonation produces a neutral hemiacetal')
f.text(225,160,'CH₃CH(OH)O⁺HCH₃',27)
f.line(425,150,550,150,arrow=True,color=BLUE)
f.text(730,160,'CH₃CH(OH)OCH₃',29)
f.text(480,95,'water accepts the proton',24)
f.text(480,280,'Oxonium + H₂O → hemiacetal + H₃O⁺',27)
f.text(480,345,'The acid catalyst is regenerated; proton transfer is a real step.',22)
save(u,6,'hemiacetal-neutral',f)

f=Figure('Hemiacetals and acetals have different groups at the same carbon')
tetra(f,240,190,top='OH',bottom='OCH₃')
tetra(f,720,190,top='OCH₃',bottom='OCH₃')
f.text(240,355,'hemiacetal: OH + OR',25).text(720,355,'acetal: OR + OR',25)
save(u,7,'hemiacetal-acetal',f)

f=Figure('Water removal can shift an acid-catalysed acetal equilibrium')
f.text(480,115,'CH₃CHO + 2 CH₃OH ⇌ CH₃CH(OCH₃)₂ + H₂O',28)
f.line(240,210,720,210,arrow=True,color=BLUE)
f.text(480,175,'remove water: favour acetal',25)
f.line(720,285,240,285,arrow=True,color=GREEN)
f.text(480,330,'add excess water: favour hydrolysis',25)
f.text(480,380,'Acid catalyses both directions; it does not set the equilibrium constant.',21)
save(u,8,'water-equilibrium',f)

f=Figure('An internal alcohol can form a cyclic hemiacetal')
f.text(230,185,'HO(CH₂)₄CHO',31)
f.line(405,175,535,175,arrow=True,color=BLUE)
verts=[(670,130),(755,80),(840,130),(840,230),(755,280),(670,230)]
for a,b in zip(verts,verts[1:]+verts[:1]): f.line(*a,*b)
f.rect(652,111,36,40,fill='#f7fafc').text(670,141,'O',30,color=RED)
f.line(755,80,755,35).text(755,25,'OH',25,color=RED)
f.text(230,280,'5-hydroxypentanal',25)
f.text(480,365,'A six-membered ring forms without changing the total molecular formula.',22)
save(u,9,'cyclic-hemiacetal',f)

f=Figure('The product depends on the substituent already at the carbonyl carbon')
carbonyl(f,240,190,left='CH₃',right='H')
carbonyl(f,720,190,left='CH₃',right='OR')
f.text(240,315,'aldehyde: addition',27)
f.text(720,315,'ester: possible addition–elimination',23)
f.text(480,375,'A tetrahedral intermediate can expel OR only through a suitable pathway.',21)
save(u,10,'addition-vs-substitution',f)

print('Built 50 original reactivity diagrams across five units.')
