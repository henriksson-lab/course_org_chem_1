"""Rebuild the original vector figures for the acids/bases pilot unit."""
from pathlib import Path
from math import log10
from chemfig import Figure, reaction, BLUE, RED, GREEN

out=Path(__file__).resolve().parents[1]/'acids-bases'/'figures'
def save(f,name): f.save(out/(name+'.svg'))

f=Figure('A proton transfer changes two bonds','Hydroxide accepts the acidic proton from ethanoic acid; water and ethanoate form.')
f.text(100,160,'CH₃',32).text(220,160,'C',32).text(220,65,'O',32)
f.line(150,150,195,150).line(213,120,213,80).line(225,120,225,80)
f.line(245,150,305,150).text(330,160,'O',32)
f.line(350,150,405,150).text(430,160,'H',32)
f.text(510,160,'+',32).text(600,160,'O',32).text(638,160,'H',32).text(670,140,'−',26)
f.circle(570,142,2,fill=BLUE,width=1).circle(570,154,2,fill=BLUE,width=1)
f.path('M 570 142 Q 500 40 435 140',arrow=True)
f.path('M 377 150 Q 375 220 333 178',arrow=True)
f.text(480,275,'CH₃COO⁻  +  H₂O',34)
f.text(480,345,'Lone pair → H; O—H bond pair → acid oxygen',24)
save(f,'01-proton-transfer')

f=Figure('Conjugate pairs differ by one proton')
for x,acid,base in [(245,'CH₃COOH','CH₃COO⁻'),(715,'CH₃NH₃⁺','CH₃NH₂')]:
 f.text(x,100,acid,34).line(x,130,x,235,arrow=True,color=BLUE).text(x,285,base,34).text(x+90,195,'− H⁺',25)
f.text(480,355,'The conjugate base has one less H and one less positive charge.',23)
save(f,'02-conjugate-pairs')

f=Figure('Resonance distributes charge without moving nuclei')
f.text(235,165,'CH₃—C(=O)—O⁻',32).text(725,165,'CH₃—C(—O⁻)=O',32)
f.line(425,165,530,165,arrow=True,both=True,color=BLUE)
f.text(480,280,'One ion, two equivalent contributors',28)
f.text(480,330,'Actual C—O bonds are equivalent; the ion does not switch structures.',22)
save(f,'03-resonance')

f=Figure('A lower pKa means a stronger acid','Aqueous approximate pKa comparison; values depend on conditions.')
f.line(100,220,860,220,arrow=True,color=BLUE)
for x,name,value in [(190,'ethanoic acid','4.8'),(430,'phenol','10'),(740,'ethanol','16')]:
 f.line(x,205,x,235).text(x,165,name,26).text(x,275,value,29)
f.text(480,90,'Approximate aqueous pKa',28).text(480,345,'Stronger acid ←                         → Weaker acid',25)
save(f,'04-pka-scale')

f=Figure('Use the conjugate acids to predict equilibrium')
f.text(480,110,'CH₃COOH + CH₃NH₂ ⇌ CH₃COO⁻ + CH₃NH₃⁺',31)
f.text(240,190,'Acid pKa ≈ 4.8',27).text(730,190,'Product acid pKa ≈ 10.6',27)
f.line(220,240,760,240,arrow=True,color=GREEN,width=5)
f.text(480,300,'K ≈ 10^(10.6 − 4.8) ≈ 6 × 10⁵',29)
f.text(480,355,'Equilibrium favours the weaker acid/base pair.',24)
save(f,'05-equilibrium')

f=Figure('Induction stabilises a nearby negative charge')
f.text(240,130,'CH₃—COO⁻',35).text(715,130,'ClCH₂—COO⁻',35)
f.text(240,220,'CH₃COOH: pKa ≈ 4.8',26).text(715,220,'ClCH₂COOH: pKa ≈ 2.9',26)
f.line(735,165,625,165,arrow=True,color=BLUE)
f.text(480,315,'Electron withdrawal stabilises the conjugate base.',26)
f.text(480,360,'Compare the same solvent and similar measurement conditions.',21)
save(f,'06-induction')

f=Figure('Ionisation is a continuous change','Predicted fraction of A minus in dilute aqueous solution for a monoprotic acid of pKa 4.8.')
f.line(115,310,865,310,arrow=True).line(115,310,115,70,arrow=True)
pts=[]
for j in range(181):
 ph=1+j*8/180; frac=1/(1+10**(4.8-ph)); pts.append((115+(ph-1)/8*720,310-frac*220))
f.path('M '+' L '.join(f'{x:.1f},{y:.1f}' for x,y in pts),color=BLUE,width=4)
for ph in [1,3,4.8,7,9]:
 x=115+(ph-1)/8*720; f.text(x,340,ph,20)
f.text(480,380,'pH (pKa = 4.8)',23).text(160,55,'fraction A⁻',22)
for v in [0,0.5,1]: f.text(85,315-v*220,v,20)
save(f,'07-ionisation-curve')

f=Figure('At pH = pKa the two forms have equal concentrations')
for j in range(10):
 x=140+(j%5)*80;y=125+(j//5)*100
 f.circle(x,y,26,fill='#e4f0f4' if j<5 else '#d8eee3',color=BLUE if j<5 else GREEN)
 f.text(x,y+8,'HA' if j<5 else 'A⁻',22)
f.text(735,145,'[A⁻] / [HA] = 1',30).text(735,220,'50% / 50%',32)
f.text(480,350,'A ratio of concentrations, not a molecule frozen halfway through losing H.',22)
save(f,'08-half-ionised')

f=Figure('A one-unit pH change gives a tenfold ratio change')
for x,lab,frac,ratio in [(210,'pKa − 1',1/11,'1 : 10'),(480,'pKa',0.5,'1 : 1'),(750,'pKa + 1',10/11,'10 : 1')]:
 f.rect(x-70,90,140,200,fill='#e4f0f4').rect(x-70,290-200*frac,140,200*frac,fill='#23765d')
 f.text(x,65,lab,26).text(x,325,f'{100*frac:.1f}% A⁻',24).text(x,370,ratio,24)
save(f,'09-ratios')

f=Figure('Aqueous ionisation does not determine partition by itself')
f.rect(180,80,600,105,fill='#e8edf4').rect(180,185,600,150,fill='#d8eee3')
f.text(480,135,'Organic phase: partitioning needs more information',25)
f.text(480,240,'Aqueous phase: HA ⇌ H⁺ + A⁻',31)
f.text(480,290,'pH and pKa describe this acid–base ratio',24)
f.text(480,375,'Extraction also depends on partition coefficients and phase volumes.',22)
save(f,'10-model-boundary')
print('Built 10 original acid/base diagrams.')
