"""Original vector diagrams for the five application units; no external artwork."""
from pathlib import Path
import math
from chemfig import Figure, BLUE, RED, GREEN
ROOT=Path(__file__).resolve().parents[1]
def save(f,unit,name): f.save(ROOT/unit/'figures'/(name+'.svg'))
def acyl(f,x,y,left='CH₃',right='OCH₂CH₃',carbonyl=True,top='O',bottom=None):
    f.text(x,y+10,'C',31).text(x,y-101,top,31)
    f.line(x-5,y-20,x-5,y-75)
    if carbonyl: f.line(x+5,y-20,x+5,y-75)
    f.line(x-20,y,x-95,y+25).text(x-105,y+35,left,28,anchor='end')
    f.line(x+20,y,x+95,y+25).text(x+105,y+35,right,28,anchor='start')
    if bottom:
        f.line(x,y+20,x,y+80).text(x,y+120,bottom,28)
    return f

def ring(f,cx,cy,r=70,n=6,aromatic=True,nitrogen=None):
    pts=[(cx+r*math.cos(-math.pi/2+2*math.pi*i/n),cy+r*math.sin(-math.pi/2+2*math.pi*i/n)) for i in range(n)]
    for i,((x,y),(xx,yy)) in enumerate(zip(pts,pts[1:]+pts[:1])):
        if nitrogen is not None and i==nitrogen:
            # Bond ends stop short of the letter at this vertex.
            t=.18; x,y=x+(xx-x)*t,y+(yy-y)*t
        if nitrogen is not None and (i+1)%n==nitrogen:
            t=.18;xx,yy=xx+(x-xx)*t,yy+(y-yy)*t
        f.line(x,y,xx,yy)
    if aromatic: f.circle(cx,cy,r*.52,color=BLUE,fill='none',width=2)
    if nitrogen is not None:
        x,y=pts[nitrogen];f.text(x,y+10,'N' if n==6 else 'NH',28)
    return pts

def axes(f,xlabel,ylabel):
    f.line(100,305,880,305,arrow=True).line(100,305,100,55,arrow=True)
    f.text(490,375,xlabel,23).text(175,40,ylabel,21)

def ir(name,peaks,title,note):
    f=Figure(title);axes(f,'Wavenumber / cm⁻¹ (decreases to the right)','Transmittance / %')
    pts=[]
    for j in range(500):
        wn=4000-j*3400/499
        val=94-sum(depth*math.exp(-.5*((wn-pos)/width)**2) for pos,width,depth in peaks)
        pts.append((100+j*760/499,305-val*2.3))
    f.path('M '+' L '.join(f'{x:.1f},{y:.1f}' for x,y in pts),width=3)
    for wn in [4000,3000,2000,1000]: f.text(100+(4000-wn)/3400*760,335,wn,18)
    f.text(570,65,note,19,color=RED)
    for v in [0,50,100]: f.text(75,310-v*2.3,v,17)
    return f

def nmr(title,kind='proton'):
    f=Figure(title);axes(f,'Chemical shift / ppm (decreases to the right)','Schematic intensity')
    if kind=='proton':
        for center,n,base in [(1.25,3,100),(2.05,1,150),(4.12,4,95)]:
            for k in range(n):
                shift=center+(k-(n-1)/2)*.07;x=100+(5-shift)/5*760
                amp=base*([1,2,1][k]/2 if n==3 else [1,3,3,1][k]/3 if n==4 else 1)
                f.line(x,304,x,304-amp,color=BLUE,width=3)
        for shift in [0,1,2,3,4,5]:f.text(100+(5-shift)/5*760,335,shift,18)
        f.text(220,125,'4.12 q, 2H',22).text(550,110,'2.05 s, 3H',22).text(740,165,'1.25 t, 3H',22)
    else:
        for shift in [14,21,60,171]:
            x=100+(200-shift)/200*760;f.line(x,304,x,130,color=BLUE,width=3)
            if shift==14: f.line(x,130,850,120,color=BLUE,width=1).text(860,115,shift,21)
            elif shift==21: f.line(x,130,750,100,color=BLUE,width=1).text(740,90,shift,21)
            else: f.text(x,110,shift,21)
        for shift in [0,50,100,150,200]:f.text(100+(200-shift)/200*760,335,shift,18)
    f.text(590,65,'Simulated teaching spectrum; not measured',19,color=RED)
    return f

def tlc(title,lanes=None):
    f=Figure(title);f.rect(230,45,500,300,fill='#fff',color=BLUE)
    f.line(250,310,710,310,dash='6 5').line(250,90,710,90,dash='6 5')
    lanes=lanes or [('S',[.24]),('P',[.64]),('Crude',[.24,.64])]
    for j,(name,rs) in enumerate(lanes):
        x=320+j*160;f.text(x,377,name,23)
        for rf in rs:f.circle(x,310-rf*220,10,fill=BLUE,color=BLUE)
    f.text(160,95,'front',20).text(155,315,'origin',20)
    return f

# Acyl chemistry: exact bond diagrams and elementary electron-pair movement.
u='acyl-chemistry'
f=Figure('Carboxylic acid derivatives share an acyl carbon')
for x,r in [(175,'OH'),(480,'OR′'),(785,'NH₂')]: acyl(f,x,185,left='R',right=r);f.text(x,340,{'OH':'acid','OR′':'ester','NH₂':'amide'}[r],25)
save(f,u,'01-acyl-family')
f=Figure('Addition starts by moving the carbonyl pi pair to oxygen');acyl(f,480,170,left='R',right='Y');f.text(480,320,':Nu⁻',30)
f.path('M 460 287 Q 375 240 465 180',arrow=True).path('M 486 115 Q 550 100 500 60',arrow=True)
f.text(480,380,'Nu lone pair → carbon; C=O π pair → oxygen',23);save(f,u,'02-addition-arrows')
f=Figure('The tetrahedral intermediate has four single bonds');acyl(f,480,170,left='R',right='Y',carbonyl=False,top='O⁻',bottom='Nu')
f.text(480,380,'Carbon has four bonds; total charge remains −1.',24);save(f,u,'03-tetrahedral')
f=Figure('Collapse reforms the carbonyl and releases a leaving group');acyl(f,480,170,left='R',right='Y',carbonyl=False,top=':O⁻',bottom='Nu')
f.path('M 459 70 Q 410 110 480 130',arrow=True).path('M 542 186 Q 565 125 587 179',arrow=True)
f.text(480,380,'O lone pair → C–O π bond; C–Y bond pair → Y',22);save(f,u,'04-collapse-arrows')
f=Figure('Acyl substitution replaces Y while retaining a carbonyl')
f.text(250,155,'R–C(=O)–Y + Nu⁻',30).text(715,245,'R–C(=O)–Nu + Y⁻',30)
f.line(370,195,575,195,arrow=True,color=BLUE);f.text(480,340,'Addition followed by elimination; no free acyl carbocation.',24);save(f,u,'05-net-substitution')
f=Figure('Acid-catalysed esterification is an equilibrium')
f.text(480,115,'CH₃COOH + HOCH₂CH₃',32).text(480,265,'CH₃COOCH₂CH₃ + H₂O',32)
f.line(390,155,390,215,arrow=True,color=BLUE).line(570,215,570,155,arrow=True,color=BLUE)
f.text(480,195,'H⁺ catalyst',23).text(480,355,'Proton transfers activate the carbonyl and make water a leaving group.',21);save(f,u,'06-ester-equilibrium')
f=Figure('Base-promoted ester hydrolysis traps the carboxylate')
f.text(480,100,'CH₃COOCH₂CH₃ + OH⁻',33).line(480,130,480,205,arrow=True,color=BLUE)
f.text(480,255,'CH₃COO⁻ + HOCH₂CH₃',33).text(480,345,'Final acid–base step favours carboxylate; acid work-up gives acid.',22);save(f,u,'07-hydrolysis')
f=Figure('Amide nitrogen donates into the carbonyl')
f.text(245,170,'R–C(=O)–NH₂',32).text(730,170,'R–C(–O⁻)=NH₂⁺',32)
f.line(425,160,535,160,arrow=True,both=True,color=BLUE)
f.text(480,285,'Partial C–N double-bond character limits rotation.',26)
f.text(480,340,'Contributors describe one amide; N lone pair is less available.',22);save(f,u,'08-amide-resonance')
f=Figure('Acid plus amine first gives an ammonium carboxylate salt')
f.text(480,110,'RCOOH + R′NH₂',34).line(480,140,480,215,arrow=True,color=BLUE)
f.text(480,265,'RCOO⁻ + R′NH₃⁺',34).text(480,350,'Amide synthesis requires suitable activation/conditions.',25);save(f,u,'09-salt-not-amide')
f=Figure('Biological ester and amide links share an acyl motif')
acyl(f,245,185,left='R',right='OR′');acyl(f,715,185,left='R',right='NHR′')
f.text(245,340,'ester linkage',27).text(715,340,'peptide (amide) linkage',27);save(f,u,'10-bio-links')

# Nitrogen, aromaticity and amino-acid charge.
u='amines-aromatics'
f=Figure('Count carbon substituents on nitrogen, not total N bonds')
for x,label,name in [(170,'CH₃–NH₂','primary'),(480,'(CH₃)₂NH','secondary'),(790,'(CH₃)₃N','tertiary')]:
    f.text(x,175,label,31).text(x,270,name,27)
f.text(480,345,'Each neutral example has one nitrogen lone pair.',24);save(f,u,'01-amine-family')
f=Figure('Amine protonation creates four bonds and positive charge')
f.text(250,155,'CH₃–NH₂: + H⁺',32).text(720,245,'CH₃–NH₃⁺',34)
f.line(360,195,600,195,arrow=True,color=BLUE);f.text(480,345,'The original lone pair becomes the new N–H bonding pair.',24);save(f,u,'02-amine-protonation')
f=Figure('Amine and amide nitrogen have different electron availability')
f.text(240,160,'CH₃CH₂–NH₂:',32);acyl(f,700,180,left='CH₃',right='NH₂')
f.text(240,285,'Localised lone pair',25).text(700,335,'Lone pair delocalised into C=O',23);save(f,u,'03-amine-amide')
f=Figure('Aniline nitrogen can donate into the aromatic ring')
pts=ring(f,310,200,85);f.line(310,115,310,65).text(310,45,':NH₂',27)
f.line(286,64,270,130,color=BLUE,dash='5 5')
f.text(700,170,'Lone pair overlaps with ring π system',23).text(700,245,'Aniline is less basic than ethylamine',23)
f.text(480,355,'Compare conjugate-acid pKa values measured in the same solvent.',21);save(f,u,'04-aniline')
f=Figure('Benzene has six pi electrons in a cyclic conjugated system');ring(f,300,195,105)
for i in range(6):
    a=-math.pi/2+i*math.pi/3;x=300+105*math.cos(a);y=195+105*math.sin(a)
    f.line(x,y-32,x,y+32,color=GREEN,width=8)
f.text(700,155,'Six p orbitals overlap',27).text(700,225,'Six π electrons = 4n + 2, n = 1',25)
f.text(480,355,'Aromaticity also requires a planar, continuously conjugated ring.',22);save(f,u,'05-aromaticity')
f=Figure('A ring alone does not make a compound aromatic');ring(f,270,200,90,aromatic=True);ring(f,700,200,90,aromatic=False)
f.text(270,335,'Benzene: continuous p-orbital ring',22).text(700,335,'Cyclohexane: no π network',22);save(f,u,'06-ring-comparison')
f=Figure('Pyridine has a lone pair outside its aromatic pi system');ring(f,300,210,90,nitrogen=0)
f.text(300,83,'••',30,color=BLUE);f.text(700,145,'Ring N contributes one π electron',24).text(700,220,'Lone pair can bind H⁺',25)
f.text(480,355,'Protonating this lone pair does not remove the aromatic sextet.',22);save(f,u,'07-pyridine')
f=Figure('Pyrrole uses its nitrogen lone pair in the aromatic sextet');ring(f,300,215,95,n=5,nitrogen=0)
f.text(300,75,'••',29,color=GREEN);f.line(320,82,355,167,color=GREEN,dash='5 5')
f.text(700,145,'Four π electrons from two C=C bonds',23).text(700,215,'+ two from N lone pair = six',24)
f.text(480,355,'Pyrrole nitrogen is much less basic than pyridine nitrogen.',23);save(f,u,'08-pyrrole')
f=Figure('Glycine charge depends on pH')
for x,formula,charge in [(175,'H₃N⁺–CH₂–COOH','net +1'),(480,'H₃N⁺–CH₂–COO⁻','net 0'),(785,'H₂N–CH₂–COO⁻','net −1')]:
    f.text(x,165,formula,23).text(x,245,charge,28)
f.line(130,310,830,310,arrow=True,color=BLUE).text(480,365,'Increasing pH; two distinct deprotonation steps',24);save(f,u,'09-glycine')
f=Figure('Basicity is an equilibrium comparison, not a reaction-speed ranking')
f.text(480,95,'B + H⁺ ⇌ BH⁺',34)
for x,name,pka in [(230,'ethylammonium','≈10.7'),(730,'anilinium','≈4.6')]:
    f.text(x,200,name,29).text(x,265,'aqueous pKa '+pka,26)
f.text(480,350,'Larger conjugate-acid pKa → stronger base in this comparison.',23);save(f,u,'10-basicity-comparison')

# Spectral inference. All synthetic.
u='spectroscopy'
f=Figure('A formula constrains but does not specify connectivity')
f.text(480,75,'C₄H₈O₂: degree of unsaturation = (2×4 + 2 − 8)/2 = 1',27)
for x,formula,name in [(230,'CH₃COOCH₂CH₃','ethyl ethanoate'),(730,'CH₃CH₂CH₂COOH','butanoic acid')]:
    f.text(x,190,formula,31).text(x,275,name,27)
f.text(480,350,'Both formulas fit; distinguish them using independent observations.',22);save(f,u,'01-formula-candidates')
save(ir('ester',[(1740,22,62),(1240,45,28),(1050,32,20)],'Simulated ester IR spectrum','Teaching trace: strong C=O; no broad O–H'),u,'02-ester-ir')
save(ir('acid',[(1710,28,55),(2950,280,45),(1230,40,20)],'Simulated carboxylic acid IR spectrum','Teaching trace: C=O plus very broad acid O–H'),u,'03-acid-ir')
save(ir('alcohol',[(3330,170,50),(1050,40,35)],'Simulated alcohol IR spectrum','Teaching trace: broad O–H; no C=O'),u,'04-alcohol-ir')
save(nmr('Simulated proton NMR: unknown U1'),u,'05-proton-nmr')
f=Figure('An ethyl group gives a linked quartet and triplet')
f.text(480,90,'CH₃–CH₂–O–',43)
f.text(245,145,'CH₃ sees two neighbouring H',24).text(715,145,'CH₂ sees three neighbouring H',24)
for x,h in [(170,70),(215,140),(260,70)]:f.line(x,310,x,310-h,color=BLUE,width=5)
for x,h in [(650,45),(690,135),(730,135),(770,45)]:f.line(x,310,x,310-h,color=GREEN,width=5)
f.text(220,355,'triplet, integral 3H',23).text(710,355,'quartet, integral 2H',23);save(f,u,'06-ethyl-coupling')
f=Figure('A carbonyl interrupts simple vicinal proton coupling')
f.text(480,155,'CH₃—C(=O)—O—CH₂—CH₃',40)
f.line(175,195,290,195,color=RED,width=5).text(240,260,'3H singlet',28,color=RED)
f.text(480,345,'No adjacent carbon-bound H across C=O for the acetyl methyl.',23);save(f,u,'07-acetyl-singlet')
save(nmr('Simulated carbon NMR: unknown U1','carbon'),u,'08-carbon-nmr')
f=Figure('Attach each spectral clue to the proposed structure')
for x,label in [(190,'CH₃'),(310,'C(=O)'),(450,'O'),(590,'CH₂'),(730,'CH₃')]: f.text(x,160,label,34)
for x1,x2 in [(230,260),(365,420),(480,550),(630,690)]: f.line(x1,150,x2,150)
f.line(190,180,120,260,color=BLUE).text(120,300,'2.05 ppm, 3H s',20)
f.line(310,180,330,240,color=BLUE).text(330,275,'IR C=O',20).text(330,310,'171 ppm ¹³C',20)
f.line(590,180,590,270,color=BLUE).text(590,310,'4.12 ppm, 2H q',20)
f.line(730,180,800,240,color=BLUE).text(800,280,'1.25 ppm, 3H t',20)
save(f,u,'09-assignment-map')
f=Figure('Predictions can reject a plausible-looking candidate')
f.text(480,110,'Butanoic acid predicts: acid O–H + terminal ethyl + CH₂–COOH',25)
f.text(480,195,'U1 shows: no broad O–H + O–CH₂ quartet + acetyl singlet',25)
f.line(160,265,800,265,color=RED,width=4)
f.text(480,330,'Reject the acid candidate using several independent observations.',23);save(f,u,'10-candidate-test')

# Synthesis and evaluation.
u='synthesis-evidence'
f=Figure('Two simple transformations connect an alcohol and an ester')
f.text(180,155,'CH₃CH₂OH',29).text(480,155,'CH₃COOH',29).text(790,155,'CH₃COOCH₂CH₃',25)
f.line(280,145,375,145,arrow=True,color=BLUE).line(580,145,650,145,arrow=True,color=BLUE)
f.text(325,220,'oxidation',23).text(625,220,'esterification',23)
f.text(480,330,'Example conditions must be chosen to favour the intended transformation.',22);save(f,u,'01-short-route')
f=Figure('An amount table follows a balanced equation')
f.text(480,85,'salicylic acid + acetic anhydride → aspirin + ethanoic acid',25)
for x,name,amount in [(230,'salicylic acid','10.0 mmol'),(730,'acetic anhydride','15.0 mmol')]:
    f.text(x,175,name,28).text(x,245,amount,33)
f.text(480,345,'One mole of each reagent gives one mole of aspirin (1:1:1).',23);save(f,u,'02-stoichiometry')
f=Figure('Compare amount divided by stoichiometric coefficient')
for x,amount,label in [(270,10,'salicylic acid'),(690,15,'acetic anhydride')]:
    f.rect(x-65,300-amount*12,130,amount*12,fill=BLUE if amount==10 else GREEN)
    f.text(x,330,label,24).text(x,285-amount*12,f'{amount} mmol / 1',25)
f.text(480,65,'The smaller amount/coefficient limits theoretical product.',26);save(f,u,'03-limiting-reagent')
f=Figure('Theoretical mass is a stoichiometric prediction')
f.text(190,160,'0.0100 mol',31).line(300,155,395,155,arrow=True,color=BLUE)
f.text(490,160,'× 180.16 g/mol',30).line(620,155,700,155,arrow=True,color=BLUE)
f.text(805,160,'1.8016 g',31).text(480,275,'Theoretical aspirin mass',30)
f.text(480,345,'Use consistent units; retain precision until final reporting.',23);save(f,u,'04-theoretical-mass')
f=Figure('Isolated yield compares recovered mass with theoretical mass')
f.rect(210,115,540,70,fill='#e4f0f4').rect(210,230,360,70,fill='#d8eee3')
f.text(480,160,'Theoretical 1.8016 g',30).text(390,275,'Isolated 1.20 g',30)
f.text(480,365,'100 × 1.20 / 1.8016 = 66.6% isolated yield',28);save(f,u,'05-yield')
f=tlc('Simulated TLC: crude contains two visible components');save(f,u,'06-tlc')
f=Figure('Rf measures distance relative to the solvent front')
f.line(240,310,240,70,color=BLUE,width=4).line(230,310,660,310,dash='6 5').line(230,90,660,90,dash='6 5')
f.circle(450,169,12,fill=BLUE);f.line(500,310,500,169,arrow=True,color=BLUE)
f.text(700,185,'spot: 3.2 cm',27).text(700,255,'front: 5.0 cm',27)
f.text(480,370,'Rf = 3.2 / 5.0 = 0.64 (dimensionless)',27);save(f,u,'07-rf')
save(ir('synthesis',[(1750,25,45),(1695,30,30),(2900,260,32),(1200,50,23)],'Simulated aspirin-like IR teaching trace','Two carbonyl regions; acid O–H retained'),u,'08-product-ir')
f=Figure('One spot is evidence under one separation condition')
f=tlc('Two compounds can co-migrate on a TLC plate',[('A',[.64]),('B',[.64]),('A+B',[.64])])
f.text(480,32,'Schematic counterexample: different compounds, same Rf',21,color=RED);save(f,u,'09-comigration')
f=Figure('Trace a reported result back to its observations')
f.rect(35,70,270,100).text(170,110,'Starting-reagent',24).text(170,145,'masses and molar masses',19)
f.rect(365,70,270,100).text(500,110,'Theoretical product',23).text(500,145,'1.8016 g',28)
f.rect(35,235,270,100).text(170,275,'Isolated mass',25).text(170,315,'1.20 g',28)
f.rect(660,235,270,100).text(795,275,'Isolated yield',25).text(795,315,'66.6%',28)
f.line(310,120,360,120,arrow=True,color=BLUE)
f.line(640,120,795,225,arrow=True,color=BLUE)
f.line(310,285,650,285,arrow=True,color=GREEN)
f.text(480,385,'Yield uses both the theoretical amount and the independently weighed isolate.',20)
save(f,u,'10-reproducible-chain')

# Integration: related representations and independent evidence.
u='integration'
f=Figure('Lactic acid combines several familiar chemical features')
f.text(480,160,'CH₃—CH(OH)—C(=O)OH',43)
for x,label in [(240,'alkyl group'),(480,'secondary alcohol'),(770,'carboxylic acid')]:
    f.line(x,205,x,255,color=BLUE).text(x,290,label,22)
f.text(480,330,'stereocentre',22)
save(f,u,'01-lactic-features')
f=Figure('Keep connectivity fixed when identifying stereoisomers')
for x,s in [(270,1),(690,-1)]:
    f.text(x,185,'C',30).line(x-s*20,170,x-s*95,115).text(x-s*130,105,'COOH',26)
    f.line(x+s*20,170,x+s*95,115).text(x+s*130,105,'CH₃',26)
    f.wedge(x-s*10,195,x-s*80,270,hashed=False);f.text(x-s*110,310,'OH',27)
    f.wedge(x+s*10,195,x+s*80,270,hashed=True);f.text(x+s*110,310,'H',27)
f.line(480,65,480,335,dash='6 6',color=BLUE);f.text(480,375,'Mirror-image pair; same connectivity, opposite tetrahedral arrangement.',22);save(f,u,'02-lactic-mirrors')
f=Figure('Acid–base change preserves the carbon skeleton')
f.text(480,100,'CH₃CH(OH)COOH + OH⁻',33).line(480,130,480,210,arrow=True,color=BLUE)
f.text(480,260,'CH₃CH(OH)COO⁻ + H₂O',33)
f.text(480,345,'Removing the carboxyl proton does not invert the stereocentre.',24);save(f,u,'03-lactic-ionisation')
f=Figure('Oxidation changes the alcohol carbon into a carbonyl carbon')
f.text(245,170,'CH₃CH(OH)COOH',31).text(730,170,'CH₃C(=O)COOH',31)
f.line(405,155,570,155,arrow=True,color=BLUE).text(480,110,'oxidation',23)
f.text(480,280,'Lactic acid → pyruvic acid',29).text(480,345,'The former tetrahedral stereocentre becomes trigonal planar.',24);save(f,u,'04-lactic-redox')
f=Figure('Esterification replaces the acidic OH as a net transformation')
f.text(480,100,'CH₃CH(OH)COOH + CH₃CH₂OH',30).line(480,130,480,210,arrow=True,color=BLUE)
f.text(480,255,'CH₃CH(OH)COOCH₂CH₃ + H₂O',30)
f.text(480,345,'Acid-catalysed equilibrium; secondary alcohol remains in this target.',23);save(f,u,'05-lactic-ester')
f=Figure('Choose a mechanism from the site and the structural change')
for x,draw,cap in [(240,'R₂C=O → R₂C(OH)Nu','addition: carbonyl lost'),(720,'RC(=O)Y → RC(=O)Nu','substitution: carbonyl retained')]:
    f.text(x,160,draw,27).text(x,265,cap,25)
f.text(480,350,'Both begin with nucleophilic attack; the available leaving group matters.',22);save(f,u,'06-mechanism-choice')
f=Figure('The same molecule can carry acid and basic sites')
f.text(480,150,'H₃N⁺—CH(CH₃)—COO⁻',43)
f.text(245,260,'protonated amine',27).text(730,260,'carboxylate',27)
f.line(265,210,300,180,color=BLUE).line(725,210,710,180,color=BLUE)
f.text(480,345,'Net zero charge does not mean absence of formal charges.',25);save(f,u,'07-amino-acid')
f=Figure('Require agreement between different measurements')
f.text(480,95,'C₄H₈O₂: one degree of unsaturation',31)
f.text(240,195,'IR: carbonyl; no broad O–H',25).text(730,195,'¹H NMR: 3H s + 2H q + 3H t',25)
f.line(240,230,420,290,arrow=True,color=BLUE).line(730,230,550,290,arrow=True,color=BLUE)
f.text(480,330,'CH₃COOCH₂CH₃ fits both observations',28);save(f,u,'08-evidence-combination')
f=Figure('High isolated yield cannot rescue contradictory identity evidence')
f.rect(100,105,300,45,fill='#dce5ea').rect(100,180,283,45,fill='#23765d')
f.text(250,90,'Theoretical: 1.8016 g',23).text(250,265,'Isolated: 1.70 g (94.4%)',23)
f.line(520,260,865,260)
for x,h,col in [(560,100,BLUE),(650,150,BLUE),(710,80,RED),(815,125,BLUE)]: f.line(x,260,x,260-h,color=col,width=4)
f.text(700,80,'NMR: unexpected extra signal',23).text(710,305,'Extra set conflicts with pure product',19,color=RED)
f.text(480,365,'Investigate identity, impurities and residual solvent before claiming success.',21);save(f,u,'09-yield-identity')
f=Figure('A scientific explanation links claim, observation and chemical reason')
f.text(480,75,'Claim: the ester candidate best fits this supplied evidence',26)
f.text(300,175,'CH₃–C(=O)–O–CH₂–CH₃',32)
f.line(220,210,220,275,color=BLUE).text(220,315,'C=O region in IR',23)
f.line(415,210,650,245,color=GREEN).text(710,215,'O–CH₂: 2H quartet',24).text(710,280,'No broad acid O–H',24)
f.text(480,375,'Limit: this supports the candidate among the structures tested.',22);save(f,u,'10-claim-evidence')
print('Built 50 original application-unit SVG figures.')
