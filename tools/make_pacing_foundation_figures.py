"""Additional original diagrams for short, staged foundation lectures."""
import math
from pathlib import Path
from chemfig import Figure, BLUE, RED, GREEN
ROOT=Path(__file__).resolve().parents[1]
def save(u,n,f):f.save(ROOT/u/'figures'/f'{n:02d}-pace.svg')
def pair(title,a,b,la='',lb=''):
 f=Figure(title);f.text(245,175,a,33).text(715,175,b,33)
 f.text(245,270,la,24).text(715,270,lb,24);return f

def newman(f,cx,cy,front,back,angle=60):
 f.circle(cx,cy,38)
 for groups,off,start,col in [(back,angle,38,'#153747'),(front,0,0,BLUE)]:
  for j,label in enumerate(groups):
   a=math.radians(-90+off+j*120);r=103
   f.line(cx+start*math.cos(a),cy+start*math.sin(a),cx+r*math.cos(a),cy+r*math.sin(a),color=col)
   rr=139 if off else 125
   f.text(cx+rr*math.cos(a),cy+rr*math.sin(a)+7,label,22,color=col)
 f.circle(cx,cy,4,fill=BLUE)
def fischer(f,x,y,left='OH',right='H',top='COOH',bottom='CH₃'):
 f.line(x-55,y,x+55,y).line(x,y-65,x,y+65)
 f.text(x-95,y+8,left,25).text(x+95,y+8,right,25).text(x,y-90,top,25).text(x,y+105,bottom,25)
def curve(title,fn,xmin,xmax,ylab='fraction',xticks=None):
 f=Figure(title);f.line(100,300,870,300,arrow=True).line(100,300,100,65,arrow=True)
 pts=[(100+i*750/200,300-220*fn(xmin+i*(xmax-xmin)/200)) for i in range(201)]
 f.path('M '+' L '.join(f'{x:.2f},{y:.2f}' for x,y in pts),width=4)
 for v in (xticks or [xmin,(xmin+xmax)/2,xmax]):f.text(100+(v-xmin)/(xmax-xmin)*750,335,f'{v:g}',22)
 f.text(190,45,ylab,24)
 for v in [0,.5,1]:f.text(70,307-v*220,v,20)
 return f

def zig(f,pts):
 for a,b in zip(pts,pts[1:]):f.line(*a,*b)

u='structures'
rows=[
('Methanol retains an explicit oxygen','CH₃OH','H₃C—O—H','one carbon, one oxygen','four total H'),
('Changing view does not change neighbours','CH₃CH₂OH','HOCH₂CH₃','O attached to C and H','same connected molecule'),
('Propane hydrogen audit','CH₃—CH₂—CH₃','C₃H₈','terminal C: three H each','middle C: two H'),
('Branching changes a connected framework','CH₃CH₂CH₂CH₃','(CH₃)₃CH','butane: unbranched','2-methylpropane: branched'),
('The alkene carbons each have one H','CH₃—CH=CH—CH₃','C₄H₈','each central C: 1 + 2 + 1 = 4','double bond counts twice'),
('A terminal alkyne has one terminal H','CH₃—C≡C—H','C₃H₄','middle C: single + triple','terminal C: triple + H'),
('Nitrogen changes charge when it gains H','CH₃NH₂:','CH₃NH₃⁺','three N bonds, one lone pair','four N bonds, no lone pair'),
('Proton loss changes oxygen charge','CH₃CH₂OH','CH₃CH₂O⁻','two O bonds, two lone pairs','one O bond, three lone pairs'),
('Locate the carbonyl carbon neighbours','CH₃—C(=O)—H','CH₃—C(=O)—CH₃','aldehyde: one attached H','ketone: two carbon groups'),
('Acid and ester differ beyond oxygen count','CH₃—C(=O)—OH','CH₃—C(=O)—OCH₃','carboxylic acid','ester'),
('An ether lacks the adjacent carbonyl','CH₃CH₂—O—CH₃','CH₃C(=O)—O—CH₃','ether','ester'),
('Count the atoms inside parentheses','CH₃CH(OH)COOH','C₃H₆O₃','one H on central carbon','H on each OH is explicit'),
('A nitrogen group is not an extra carbon','CH₃CH₂NH₂','C₂H₇N','N symbol is explicit','two N–H bonds'),
('Number from the principal group end','CH₃CH₂CH₂OH','HOCH₂CH₂CH₃','propan-1-ol','same name after redrawing'),
('Place the OH locant explicitly','CH₃CH(OH)CH₃','propan-2-ol','OH on carbon 2','three-carbon parent chain'),
('Two familiar groups in one molecule','H₂N—CH₂—COOH','CH₃CH(OH)COOH','amine + carboxylic acid','alcohol + carboxylic acid'),
('Use units to decide the calculation','1.3812 g ÷ 138.12 g/mol','0.01000 mol','g cancels','10.00 mmol'),
('A molar mass is a conversion factor','0.00500 mol × 180.16 g/mol','0.9008 g','mol cancels','predicted mass'),
('A missing vertex changes formula','CH₃CH₂CH₂OH','CH₃CH₂OH','three carbons','two carbons'),
('Two oxygens do not uniquely define a group','HOCH₂CH₂OH','CH₃C(=O)OH','diol: no carbonyl','acid: carbonyl + hydroxyl')]
for n,row in enumerate(rows,11):save(u,n,pair(*row))
# One explicit branched bond diagram replaces the condensed-only branching view.
f=Figure('A branch adds a carbon vertex, not an invisible label')
zig(f,[(210,240),(400,140),(590,240)]);f.line(400,140,400,60)
for x,y,t in [(210,285,'CH₃'),(400,115,'CH'),(590,285,'CH₃'),(400,40,'CH₃')]:f.text(x,y,t,28)
f.text(480,355,'(CH₃)₃CH has four carbon atoms and ten H atoms.',25);save(u,14,f)

u='properties'
rows=[
('An O–H bond is polar even when neutral','R—Oδ−—Hδ+','R—O⁻','partial charge','formal charge'),
('A nonpolar region can dominate a large molecule','CH₃OH','CH₃(CH₂)₅OH','one OH / one carbon','one OH / six carbons'),
('A ketone oxygen accepts from water','CH₃C(=O)CH₃','H—O—H','carbonyl O: acceptor','water O–H: donor'),
('Alcohol and amine both offer donor H','CH₃CH₂OH','CH₃CH₂NH₂','O–H donor; O accepts','N–H donor; N accepts'),
('Tertiary amine has no N–H donor','(CH₃)₃N:','(CH₃)₂NH','acceptor only at N','donor and acceptor at N'),
('Compare ester with ether','CH₃COOCH₃','CH₃OCH₃','two O acceptor sites; no O–H','one O acceptor; no O–H'),
('Neutrality does not mean nonpolarity','CH₃CH₂OH','CH₃CH₂O⁻','neutral, polar molecule','charged conjugate base'),
('An acid OH has a distinctive local environment','CH₃C(=O)OH','CH₃CH₂OH','acid OH: very broad IR region','alcohol OH: broad IR region'),
('A carbonyl signal tests a drawn bond','CH₃C(=O)CH₃','CH₃CH₂CH₂OH','C=O expected','no C=O expected'),
('Do not turn a property into unique identity','CH₃C(=O)CH₃','CH₃CH₂CHO','ketone: C=O band','aldehyde: C=O band')]
for n,row in enumerate(rows,11):save(u,n,pair(*row))
# Additional IR plots use deliberately schematic Gaussian absorptions.
for n,center,width,kind in [(21,3350,210,'alcohol O–H'),(22,1715,35,'carbonyl'),(23,2950,330,'acid O–H')]:
 f=curve('Schematic '+kind+' transmittance',lambda q,c=center,w=width:.9-.65*math.exp(-((q-c)/w)**2),4000,500,'transmittance (fraction)',[4000,3000,2000,1000,500])
 f.text(480,380,'Wavenumber / cm⁻¹; simulated teaching feature',22);save(u,n,f)
for n,title,ys,labels in [(24,'Illustrative normal boiling-point comparison',[-24,78],['dimethyl ether','ethanol']),(25,'A larger chain increases nonpolar surface',[1,6],['methanol C count','hexan-1-ol C count'])]:
 f=Figure(title)
 for x,v,label in zip([260,700],ys,labels):
  base=260;h=(v+30)*1.6 if n==24 else v*28
  f.rect(x-65,base-h,130,h,fill=BLUE).text(x,base-h-15,str(v),27).text(x,325,label,24)
 f.text(480,375,'°C at about one atmosphere' if n==24 else 'Same one OH group; different nonpolar region',23);save(u,n,f)
for n,absorbance in [(26,False),(27,True)]:
 fn=lambda q:.7*math.exp(-((q-1715)/65)**2)
 f=curve('Same band, different ordinate',fn if absorbance else lambda q:1-fn(q),4000,500,'Absorbance (schematic)' if absorbance else 'Transmittance (fraction)',[4000,2000,500])
 f.text(480,380,'Wavenumber / cm⁻¹; simulated feature',23);save(u,n,f)
for n,forms,labels in [(28,['H₂O','CH₃OCH₃'],['donor + acceptor','acceptor, no O–H']),(29,['CH₃COOH','CH₃COO⁻'],['neutral acid','carboxylate; ion–dipole hydration']),(30,['C=O present','C=O absent'],['supports a carbonyl candidate','tests a different structural feature'])]:save(u,n,pair('Choose evidence appropriate to the claim',*forms,*labels))

u='acids-bases'
rows=[
('Conjugates differ by one H plus','CH₃COOH ⇌ CH₃COO⁻ + H⁺','NH₄⁺ ⇌ NH₃ + H⁺','acid → conjugate base','charge decreases by one'),
('An alcohol proton loss gives alkoxide','CH₃CH₂OH','CH₃CH₂O⁻ + H⁺','two O bonds','one O bond; negative charge'),
('A neutral amine accepts one proton','CH₃NH₂: + H–A','CH₃NH₃⁺ + A⁻','N pair → acidic H','H–A pair → A'),
('Water can donate or accept a proton','H₃O⁺ / H₂O','H₂O / OH⁻','water as conjugate base','water as acid'),
('Resonance is not proton transfer','CH₃C(=O)O⁻ ↔ CH₃C(O⁻)=O','CH₃COOH ⇌ CH₃COO⁻ + H⁺','same atoms and nuclei','different molecular species'),
('Stabilise a base to strengthen its acid','CH₃CH₂O⁻','CH₃COO⁻','charge largely localised','equivalent O delocalisation'),
('The pKa difference sets an equilibrium ratio','pKa(HA) = 4.8','pKa(BH⁺) = 10.6','reactant acid','product acid is weaker'),
('A weaker acid can be the product','HA + B ⇌ A⁻ + BH⁺','K ≈ 10^(10.6 − 4.8)','compare the two acids','about 6.3 × 10⁵'),
('Reverse the reaction, invert the constant','A⁻ + BH⁺ ⇌ HA + B','K(reverse) ≈ 1.6 × 10⁻⁶','same equilibrium','reciprocal constant'),
('A substituent acts through sigma bonds','CH₃COOH: pKa ≈ 4.8','ClCH₂COOH: pKa ≈ 2.9','same acid family','electron withdrawal stabilises base')]
for n,row in enumerate(rows,11):save(u,n,pair(*row))
for n,pka in [(21,4.8),(22,6.0)]:
 f=curve('Predicted monoprotic-acid ionisation',lambda ph,p=pka:1/(1+10**(p-ph)),1,9,'fraction A⁻',[1,3,4.8,6,9])
 f.text(480,380,f'pH; supplied pKa = {pka}; dilute aqueous model',22);save(u,n,f)
for n,ph,ratio in [(23,3.8,.1),(24,4.8,1),(25,5.8,10),(26,6.8,100)]:
 f=Figure('Ratio and fraction answer different questions')
 f.text(480,75,f'pH {ph}; pKa 4.8',30)
 f.text(245,170,f'[A⁻]/[HA] = {ratio:g}',31).text(715,170,f'fraction A⁻ = {ratio/(1+ratio):.4f}',29)
 f.rect(170,255,620,40,fill='#e4f0f4').rect(170,255,620*ratio/(1+ratio),40,fill=GREEN)
 f.text(480,355,'Fraction = ratio / (1 + ratio)',26);save(u,n,f)
for n,count in [(27,1),(28,5),(29,9)]:
 f=Figure('Species fraction describes a molecular population')
 for i in range(10):
  x=115+i*80;f.circle(x,165,26,color=GREEN if i<count else BLUE,fill='#d8eee3' if i<count else '#e4f0f4');f.text(x,173,'A⁻' if i<count else 'HA',21)
 f.text(480,285,f'{count} A⁻ and {10-count} HA in an illustrative ten-particle sketch',26)
 f.text(480,345,'Actual solution species interconvert; each molecule has a definite form.',22);save(u,n,f)
save(u,30,pair('pH is not the same quantity as titrant volume','input: specified pH','input: added titrant volume','speciation curve computes fraction','titration model must compute pH'))

u='conformations'
for n,front,back,angle,caption in [
(11,['H','H','H'],['H','H','H'],60,'ethane: staggered'),(12,['H','H','H'],['H','H','H'],0,'ethane: eclipsed'),
(13,['CH₃','H','H'],['CH₃','H','H'],180,'butane: anti'),(14,['CH₃','H','H'],['CH₃','H','H'],60,'butane: +60° gauche'),
(15,['CH₃','H','H'],['CH₃','H','H'],300,'butane: −60° gauche'),(16,['CH₃','H','H'],['CH₃','H','H'],0,'butane: methyl/methyl eclipse'),
(17,['CH₃','H','H'],['CH₃','H','H'],120,'butane: methyl/H eclipse'),(18,['CH₃','H','H'],['CH₃','H','H'],240,'butane: another methyl/H eclipse')]:
 f=Figure(caption);newman(f,480,190,front,back,angle);f.text(480,375,caption,25);save(u,n,f)
for n,title,a,b,la,lb in [
(19,'Label the chosen bond before projecting','C1—C2—C3—C4','view C2 → C3','C2 front','C3 rear'),
(20,'Connectivity does not change during rotation','CH₃CH₂CH₂CH₃','CH₃CH₂CH₂CH₃','anti view','gauche view'),
(21,'Chair labels describe different relationships','axial / equatorial','up / down','local conformation','side of ring'),
(22,'A ring flip follows a substituent','axial-up CH₃','equatorial-up CH₃','initial chair','flipped chair'),
(23,'Follow a down substituent too','axial-down H','equatorial-down H','initial chair','flipped chair'),
(24,'Cis means on the same side of the ring','up + up','down + down','cis relationship','also cis'),
(25,'Trans means opposite sides of the ring','up + down','down + up','trans relationship','also trans'),
(26,'Numerical model inputs have units','ΔG = 3500 J/mol','RT = 2478.96 J/mol','relative state free energy','at T = 298.15 K')]:save(u,n,pair(title,a,b,la,lb))
for n,dg in [(27,3500),(28,0),(29,7000)]:
 w=math.exp(-dg/(8.314462618*298.15));p=1/(1+w);f=Figure('A two-state equilibrium population')
 for x,v,label in [(260,p,'A: ΔG = 0'),(700,1-p,f'B: ΔG = {dg} J/mol')]:
  f.rect(x-65,300-v*220,130,v*220,fill=BLUE).text(x,280-v*220,f'{v:.3f}',28).text(x,345,label,24)
 f.text(480,45,'298.15 K; free energies include degeneracy',24);save(u,n,f)
save(u,30,pair('Free-energy difference and barrier are distinct','G(B) − G(A)','G‡ − G(A)','equilibrium population ratio','interconversion rate'))

u='stereochemistry'
for n,left,right,caption in [(11,'OH','H','S-lactic acid'),(12,'H','OH','R-lactic acid')]:
 f=Figure(caption);fischer(f,480,185,left,right);f.text(480,375,caption+'; horizontal bonds toward you',24);save(u,n,f)
save(u,13,pair('Two identical groups prevent a simple stereocentre','CH₃CH(OH)CH₃','CH₃CH(OH)COOH','two CH₃ substituents','four different substituents'))
save(u,14,pair('Atomic number resolves the first comparison','Br > Cl > O > N > C > H','OH > COOH > CH₃ > H','compare attached atoms first','then resolve carbon tie'))
save(u,15,pair('Resolve the carbon tie one shell out','COOH carbon: O, O, O','CH₃ carbon: H, H, H','duplicate atom for C=O','first difference favours acid group'))
f=Figure('Rotating a Fischer projection 180 degrees preserves configuration');fischer(f,245,185,'OH','H');fischer(f,715,185,'H','OH','CH₃','COOH');f.text(480,375,'Both drawings represent the same S configuration.',25);save(u,16,f)
f=Figure('One interchange inverts a configuration');fischer(f,245,185,'OH','H');fischer(f,715,185,'H','OH');f.text(480,375,'Swap OH and H only: opposite configurations.',25);save(u,17,f)
for n,mode in [(18,'Z'),(19,'E')]:
 f=Figure(mode+'-but-2-ene');f.line(390,185,570,185).line(390,195,570,195)
 for x,top,bottom,direction in [(390,'CH₃','H',-1),(570,'CH₃' if mode=='Z' else 'H','H' if mode=='Z' else 'CH₃',1)]:
  f.line(x,190,x+direction*90,110).line(x,190,x+direction*90,270)
  f.text(x+direction*120,95,top,29).text(x+direction*120,310,bottom,29)
 f.text(480,375,'Higher-priority groups '+('together' if mode=='Z' else 'opposite'),26);save(u,n,f)
save(u,20,pair('Check whether E/Z is defined','CH₂=CHCH₃','(CH₃)₂C=CHCH₃','one end has H and H','one end has CH₃ and CH₃'))
save(u,21,pair('Enantiomers invert every stereocentre','(R,R) ↔ (S,S)','(R,S) ↔ (S,R)','generic unsymmetrical two-centre framework','test symmetry before counting distinct species'))
save(u,22,pair('Diastereomers are not mirror images','(R,R) versus (R,S)','E versus Z','one centre differs in this framework','alkene configurations also qualify'))
f=Figure('Meso symmetry is a property of the whole structure');x=480;f.line(x,85,x,300)
f.text(x,55,'COOH',28).text(x,340,'COOH',28)
for y in [145,235]:f.line(x-65,y,x+65,y).text(x-100,y+8,'H',25).text(x+100,y+8,'OH',25)
f.line(220,190,740,190,dash='7 6',color=BLUE);f.text(480,385,'Equivalent halves: meso-tartaric acid is achiral.',24);save(u,23,f)
for n,r in [(24,50),(25,70),(26,95),(27,90)]:
 f=Figure('Two-enantiomer composition and excess')
 f.rect(120,140,720,70,fill=GREEN).rect(120,140,720*r/100,70,fill=BLUE)
 f.text(245,105,f'{r}% R',30).text(715,105,f'{100-r}% S',30)
 f.text(480,290,f'ee = {abs(2*r-100)}%'+(' (racemic)' if r==50 else ' in favour of R'),30)
 f.text(480,355,'R% + S% = 100; ee = |R% − S%|',26);save(u,n,f)
save(u,28,pair('Configuration and optical rotation are independent labels','R / S','(+) / (−)','assigned from geometry and priorities','measured under stated conditions'))
save(u,29,pair('A chiral environment can distinguish enantiomers','same achiral bulk properties','different chiral interactions possible','same solvent, pressure, temperature','binding differences need measurements'))
save(u,30,pair('A stereochemical claim needs a complete representation','four different substituents','explicit depth or projection convention','find a candidate centre','then assign and compare'))
print('Built 100 additional foundation pacing diagrams.')
