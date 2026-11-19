(function (global) {
'use strict';
const R = 8.314462618;
const carbon = (x,y) => ({element:'C',x,y});
const oxygen = (x,y) => ({element:'O',x,y});
const nitrogen = (x,y) => ({element:'N',x,y});
// Neutral, ordinary-valence examples only; every edge stores its bond order.
const molecules = {
 ethanol:{name:'Ethanol',condensed:'CH₃CH₂OH',atoms:[carbon(80,175),carbon(220,100),oxygen(360,175)],bonds:[[0,1,1],[1,2,1]],groups:{alcohol:[2]},donors:[2],acceptors:[2],boiling:78},
 ether:{name:'Dimethyl ether',condensed:'CH₃OCH₃',atoms:[carbon(80,175),oxygen(220,100),carbon(360,175)],bonds:[[0,1,1],[1,2,1]],groups:{ether:[1]},donors:[],acceptors:[1],boiling:-24},
 propane:{name:'Propane',condensed:'CH₃CH₂CH₃',atoms:[carbon(80,175),carbon(220,100),carbon(360,175)],bonds:[[0,1,1],[1,2,1]],groups:{},donors:[],acceptors:[]},
 propanone:{name:'Propanone',condensed:'CH₃C(=O)CH₃',atoms:[carbon(80,175),carbon(220,175),carbon(360,175),oxygen(220,65)],bonds:[[0,1,1],[1,2,1],[1,3,2]],groups:{ketone:[1,3]},donors:[],acceptors:[3],boiling:56},
 acid:{name:'Ethanoic acid',condensed:'CH₃C(=O)OH',atoms:[carbon(80,175),carbon(220,175),oxygen(220,65),oxygen(360,175)],bonds:[[0,1,1],[1,2,2],[1,3,1]],groups:{'carboxylic acid':[1,2,3]},donors:[3],acceptors:[2]},
 lactic:{name:'Lactic acid',condensed:'CH₃CH(OH)COOH',atoms:[carbon(50,185),carbon(165,125),oxygen(165,40),carbon(285,185),oxygen(285,65),oxygen(410,185)],bonds:[[0,1,1],[1,2,1],[1,3,1],[3,4,2],[3,5,1]],groups:{alcohol:[2],'carboxylic acid':[3,4,5]},donors:[2,5],acceptors:[2,4]},
 glycine:{name:'Glycine (neutral drawn form)',condensed:'H₂NCH₂COOH',atoms:[nitrogen(55,175),carbon(170,120),carbon(285,175),oxygen(285,65),oxygen(410,175)],bonds:[[0,1,1],[1,2,1],[2,3,2],[2,4,1]],groups:{amine:[0],'carboxylic acid':[2,3,4]},donors:[0,4],acceptors:[0,3]}
};
function atomHydrogens(molecule) {
 const m=typeof molecule==='string'?molecules[molecule]:molecule;
 if(!m) throw new Error('Unknown molecule');
 const valence={C:4,O:2,N:3};
 return m.atoms.map((a,i)=>valence[a.element]-m.bonds.reduce((sum,b)=>sum+(b[0]===i||b[1]===i?b[2]:0),0));
}
function molecularFormula(key) {
 const m=molecules[key],counts={C:0,H:0,N:0,O:0};
 if(!m) throw new Error('Unknown molecule');
 m.atoms.forEach(a=>counts[a.element]++);counts.H=atomHydrogens(m).reduce((a,b)=>a+b,0);
 return ['C','H','N','O'].filter(k=>counts[k]).map(k=>k+(counts[k]>1?counts[k]:'')).join('');
}
function siteInventory(key) {
 const m=molecules[key];if(!m)throw new Error('Unknown molecule');
 return {donorSites:m.donors.length,acceptorSites:m.acceptors.length,boilingC:m.boiling??null};
}
function groupMatch(key,group) {
 const m=molecules[key];if(!m)throw new Error('Unknown molecule');
 return {present:Object.prototype.hasOwnProperty.call(m.groups,group),atomIndices:(m.groups[group]||[]).slice(),groups:Object.keys(m.groups)};
}
function acidFraction(pH,pKa) {
 if(!Number.isFinite(pH)||!Number.isFinite(pKa))throw new Error('pH and pKa must be finite');
 return 1/(1+10**(pKa-pH));
}
function extractionModel(pH,pKa,P,Vaq,Vorg) {
 if(![pH,pKa,P,Vaq,Vorg].every(Number.isFinite)||P<=0||Vaq<=0||Vorg<=0)throw new Error('Positive partition coefficient and volumes required');
 const ratio=10**(pH-pKa),D=P/(1+ratio),organicFraction=D*Vorg/(Vaq+D*Vorg);
 const aqueousNeutralFraction=(1-organicFraction)/(1+ratio);
 return {ratio,D,organicFraction,aqueousNeutralFraction,aqueousIonisedFraction:1-organicFraction-aqueousNeutralFraction};
}
function butaneModel(angle) {
 if(!Number.isFinite(angle))throw new Error('Finite angle required');
 const phi=((angle%360)+360)%360,anchors=[1.9,.6,1.3,0,1.3,.6,1.9],i=Math.floor(phi/60),t=(phi-i*60)/60;
 const energy=anchors[i]+(anchors[i+1]-anchors[i])*(1-Math.cos(Math.PI*t))/2;
 let conformation='Between canonical arrangements';
 if(Math.abs(phi-180)<1e-8)conformation='Anti; staggered';
 else if(Math.abs(phi-60)<1e-8||Math.abs(phi-300)<1e-8)conformation='Gauche; staggered';
 else if(Math.abs(phi)<1e-8)conformation='Methyl/methyl eclipsed';
 else if(Math.abs(phi-120)<1e-8||Math.abs(phi-240)<1e-8)conformation='Methyl/H eclipsed';
 return {angle:phi,energy,conformation};
}
const initialFischer=['COOH','H','CH₃','OH']; // top, right, bottom, left
function fischerConfiguration(groups) {
 if(groups.length!==4||new Set(groups).size!==4||!initialFischer.every(g=>groups.includes(g)))throw new Error('Four lactic-acid substituents required');
 const positions=[[0,-1],[1,0],[0,1],[-1,0]],p=['OH','COOH','CH₃'].map(g=>positions[groups.indexOf(g)]);
 const determinant=(p[1][0]-p[0][0])*(p[2][1]-p[0][1])-(p[1][1]-p[0][1])*(p[2][0]-p[0][0]);
 const clockwise=determinant>0,lowestToward=[1,3].includes(groups.indexOf('H'));
 return (lowestToward?!clockwise:clockwise)?'R':'S';
}
function rotateFischer(groups,quarterTurns) {
 const k=((quarterTurns%4)+4)%4;
 return groups.map((_,j)=>groups[(j-k+4)%4]);
}
function swapFischer(groups,a,b) {const out=groups.slice();[out[a],out[b]]=[out[b],out[a]];return out;}
const mechanisms={
 acid:{name:'Hydroxide deprotonates ethanoic acid',netCharge:-1,product:'CH₃COO⁻ + H₂O',first:['Hydroxide lone pair → acidic H','Acidic H → hydroxide O','Hydroxide lone pair → methyl C'],second:['Acid O–H bond pair → acid O','Acid O–H bond pair → H','Acid O lone pair → base O'],reason:'The base pair makes a new bond to H; the original O–H pair stays on acid O. Products contain carboxylate and water; total charge stays −1.'},
 carbonyl:{name:'Hydroxide adds to ethanal',netCharge:-1,product:'CH₃CH(OH)O⁻',first:['Hydroxide lone pair → carbonyl C','Carbonyl C → hydroxide O','Hydroxide lone pair → methyl C'],second:['Carbonyl π pair → carbonyl O','Keep C=O unchanged','Carbonyl π pair → attached H'],reason:'The new C–O bond forms while C=O opens to O⁻. The original carbonyl carbon now has four single bonds. Keeping C=O would give five bond orders at carbon.'},
 amine:{name:'Methylamine accepts a proton from hydronium',netCharge:1,product:'CH₃NH₃⁺ + H₂O',first:['Nitrogen lone pair → acidic H','Acidic H → nitrogen','Nitrogen lone pair → donor O'],second:['Hydronium O–H bond pair → O','Hydronium O–H bond pair → H','Remove a nitrogen electron'],reason:'The N pair becomes a new N–H bond. Nitrogen has four bonds and charge +1; hydronium becomes neutral water. Net +1 is conserved.'}
};
function mechanismModel(key,first,second,charge) {
 const m=mechanisms[key];if(!m)throw new Error('Unknown elementary step');
 return {validFirst:Number(first)===0,validSecond:Number(second)===0,validCharge:Number(charge)===m.netCharge,correct:Number(first)===0&&Number(second)===0&&Number(charge)===m.netCharge,expectedCharge:m.netCharge,carbonBondOrder:key==='carbonyl'&&Number(first)===0?(Number(second)===0?4:Number(second)===1?5:null):null};
}
function conformerModel(deltaG,temperature,reference=0) {
 if(![deltaG,temperature,reference].every(Number.isFinite)||temperature<=0)throw new Error('Finite energy and positive kelvin temperature required');
 const relativeWeight=Math.exp(-deltaG/(R*temperature)),fractionA=1/(1+relativeWeight);
 return {deltaG,temperature,reference,weightA:Math.exp(-reference/(R*temperature)),weightB:Math.exp(-(reference+deltaG)/(R*temperature)),fractionA,fractionB:1-fractionA};
}
function checks() {
 const near=(a,b)=>Math.abs(a-b)<1e-9;
 const ex=extractionModel(4.8,4.8,10,50,50),s=initialFischer;
 return [
  {name:'Ethanol formula and implicit hydrogen count',pass:molecularFormula('ethanol')==='C2H6O'&&atomHydrogens('ethanol').join(',')==='3,2,1'},
  {name:'Carbonyl bond order changes H count',pass:molecularFormula('propanone')==='C3H6O'&&molecularFormula('acid')==='C2H4O2'},
  {name:'Ether accepts but has no donor',pass:siteInventory('ether').donorSites===0&&siteInventory('ether').acceptorSites===1},
  {name:'Functional-group distinction',pass:groupMatch('acid','carboxylic acid').present&&!groupMatch('acid','alcohol').present},
  {name:'Acid midpoint and one-unit change',pass:near(acidFraction(4.8,4.8),.5)&&near(acidFraction(5.8,4.8),10/11)},
  {name:'Extraction midpoint D and mass balance',pass:near(ex.D,5)&&near(ex.organicFraction,5/6)&&near(ex.organicFraction+ex.aqueousNeutralFraction+ex.aqueousIonisedFraction,1)},
  {name:'Higher pH reduces neutral-only organic recovery',pass:extractionModel(2,4.8,10,50,50).organicFraction>extractionModel(8,4.8,10,50,50).organicFraction},
  {name:'Butane anti is below gauche and eclipsed',pass:butaneModel(180).energy<butaneModel(60).energy&&butaneModel(60).energy<butaneModel(0).energy},
  {name:'Butane periodicity',pass:near(butaneModel(360).energy,butaneModel(0).energy)},
  {name:'Lactic Fischer baseline and mirror',pass:fischerConfiguration(s)==='S'&&fischerConfiguration(swapFischer(s,1,3))==='R'},
  {name:'Fischer 180 preserves;90 changes implied layout',pass:fischerConfiguration(rotateFischer(s,2))==='S'&&fischerConfiguration(rotateFischer(s,1))==='R'},
  {name:'Any one substituent interchange inverts',pass:[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]].every(p=>fischerConfiguration(swapFischer(s,...p))==='R')},
  {name:'Correct elementary steps conserve supplied charges',pass:['acid','carbonyl','amine'].every(k=>mechanismModel(k,0,0,mechanisms[k].netCharge).correct)},
  {name:'Retained pi bond gives five carbon bond orders',pass:mechanismModel('carbonyl',0,1,-1).carbonBondOrder===5},
  {name:'Conformer fractions normalise',pass:near(conformerModel(3500,298.15).fractionA+conformerModel(3500,298.15).fractionB,1)},
  {name:'Reference energy does not change fractions',pass:near(conformerModel(3500,298.15,0).fractionA,conformerModel(3500,298.15,2000).fractionA)},
  {name:'Equal free energies give equal populations',pass:near(conformerModel(0,298.15).fractionA,.5)}
 ];
}
const models={molecules,atomHydrogens,molecularFormula,siteInventory,groupMatch,acidFraction,extractionModel,butaneModel,initialFischer,fischerConfiguration,rotateFischer,swapFischer,mechanismModel,conformerModel,checks};
if(typeof module!=='undefined'&&module.exports)module.exports=models;
if(typeof document==='undefined')return;
const {initDemo,selectControl,rangeControl,numberControl,button,state,svg,text,line,circle,fmt,plot,esc}=global.ChemDemo;
const C={ink:'#153747',blue:'#176b91',green:'#23765d',red:'#af3b45'};
const sub=s=>String(s).replace(/[0-9]/g,d=>'₀₁₂₃₄₅₆₇₈₉'[Number(d)]);
function chemicalDrawing(key,{mode='skeletal',reveal=false,highlight=[],donors=false,acceptors=false,title=null}={}) {
 const m=molecules[key],hs=atomHydrogens(m);let b='';
 const labelled=i=>mode!=='skeletal'||m.atoms[i].element!=='C';
 for(const [i,j,order] of m.bonds){
  const a=m.atoms[i],z=m.atoms[j],dx=z.x-a.x,dy=z.y-a.y,len=Math.hypot(dx,dy),ux=dx/len,uy=dy/len;
  const start=labelled(i)?25:0,end=labelled(j)?25:0;
  const col=highlight.includes(i)&&highlight.includes(j)?C.red:C.ink;
  const offsets=order===2?[-4,4]:[0];
  offsets.forEach(o=>{b+=line(a.x+ux*start-uy*o,a.y+uy*start+ux*o,z.x-ux*end-uy*o,z.y-uy*end+ux*o,col,3);});
 }
 m.atoms.forEach((a,i)=>{
  if(highlight.includes(i))b+=circle(a.x,a.y,30,'#fbe3e1',C.red);
  if(donors&&m.donors.includes(i))b+=circle(a.x,a.y,36,'none',C.green);
  if(acceptors&&m.acceptors.includes(i))b+=`<rect x="${a.x-31}" y="${a.y-31}" width="62" height="62" fill="none" stroke="${C.blue}" stroke-width="3" stroke-dasharray="6 4"/>`;
  let label=a.element;
  if(a.element==='C'&&mode==='skeletal')label='';
  else if(a.element!=='C'||reveal)label+=hs[i]?('H'+(hs[i]>1?sub(hs[i]):'')):'';
  if(label)b+=text(a.x,a.y+9,label,29);
  if(reveal&&a.element==='C'&&mode==='skeletal')b+=text(a.x,a.y+40,hs[i]+' H',24,C.blue);
 });
 b+=text(240,260,title||m.name,27);
 return svg(m.name+' molecular structure',b,480,295);
}
function smallPlot(title,series,xd,yd,xlabel,ylabel,marker=null) {
 const L=67,Rp=455,T=45,B=275,x=v=>L+(v-xd[0])/(xd[1]-xd[0])*(Rp-L),y=v=>B-(v-yd[0])/(yd[1]-yd[0])*(B-T);
 let b=line(L,B,Rp,B)+line(L,B,L,T)+text(65,26,ylabel,24,C.ink,'start');
 [0,.5,1].forEach(t=>{const xv=xd[0]+t*(xd[1]-xd[0]),yv=yd[0]+t*(yd[1]-yd[0]);b+=text(x(xv),B+34,fmt(xv),25)+text(L-12,y(yv)+8,fmt(yv),25,C.ink,'end');});
 series.forEach(s=>{b+=`<polyline points="${s.points.map(p=>x(p[0])+','+y(p[1])).join(' ')}" fill="none" stroke="${s.color||C.blue}" stroke-width="4"/>`;});
 if(marker)b+=circle(x(marker[0]),y(marker[1]),7,C.red);
 b+=text(260,355,xlabel,26);return svg(title,b,480,380);
}
function setRange(input,value){input.value=value;input.dispatchEvent(new Event('input',{bubbles:true}));}
const moleculeOptions=keys=>keys.map(k=>[k,molecules[k].name]);
function mount01(){const ui=initDemo('day01');if(!ui)return;let key='ethanol',mode='skeletal',reveal=false;
 const molecule=selectControl(ui.controls,'Molecule',moleculeOptions(['ethanol','propane','propanone','acid','ether']),key,v=>{key=v;update();});
 const representation=selectControl(ui.controls,'Representation',[['skeletal','Skeletal bonds'],['atoms','Label carbon atoms'],['condensed','Condensed formula']],mode,v=>{mode=v;update();});
 const hidden=selectControl(ui.controls,'Carbon-bound hydrogen counts',[['hide','Hide counts: predict first'],['reveal','Reveal counts']], 'hide',v=>{reveal=v==='reveal';update();});
 button(ui.controls,'Reset',()=>{key='ethanol';mode='skeletal';reveal=false;molecule.value=key;representation.value=mode;hidden.value='hide';update();});
 function update(){const hs=atomHydrogens(key),m=molecules[key];ui.visual.innerHTML=mode==='condensed'?svg(m.name+' condensed formula',text(240,115,m.condensed,37)+text(240,200,'Connectivity is explicit.',27),480,260):chemicalDrawing(key,{mode,reveal});
  ui.result.textContent=reveal?'Formula: '+sub(molecularFormula(key))+'. Carbon H counts: '+m.atoms.map((a,i)=>a.element==='C'?hs[i]:null).filter(v=>v!==null).join(', ')+'.':'Count carbon endpoints/corners, then predict the hidden H. Oxygen-bound H is shown explicitly.';
  ui.detail.textContent='Ordinary neutral C has bond-order sum four; O has two. A double bond uses two bond orders at each end. Changing the representation does not change connectivity.';
  state(ui,{molecule:key,representation:mode,reveal,formula:molecularFormula(key),hydrogens:hs});}update();}
function mount02(){const ui=initDemo('day02');if(!ui)return;let a='ethanol',b='ether',show='both';
 const keys=['ethanol','ether','propanone','acid'];
 const ca=selectControl(ui.controls,'First molecule',moleculeOptions(keys),a,v=>{a=v;update();});
 const cb=selectControl(ui.controls,'Second molecule',moleculeOptions(keys),b,v=>{b=v;update();});
 const sites=selectControl(ui.controls,'Highlight interaction sites',[['both','Donors and acceptors'],['donors','Donor sites only'],['acceptors','Acceptor sites only'],['none','No highlights: predict first']],show,v=>{show=v;update();});
 button(ui.controls,'Compare alcohol and ether',()=>{a='ethanol';b='ether';ca.value=a;cb.value=b;update();});
 button(ui.controls,'Reset',()=>{a='ethanol';b='ether';show='both';ca.value=a;cb.value=b;sites.value=show;update();});
 function update(){const x=siteInventory(a),y=siteInventory(b),opts={donors:show==='both'||show==='donors',acceptors:show==='both'||show==='acceptors'};
  ui.visual.innerHTML=chemicalDrawing(a,opts)+chemicalDrawing(b,opts);
  const boiling=[a,b].map(k=>molecules[k].boiling===undefined?molecules[k].name+': no boiling value supplied':molecules[k].name+': supplied approximate normal boiling point '+molecules[k].boiling+' °C').join('; ');
  ui.result.textContent=show==='none'?'Predict donor and acceptor sites for each molecule before revealing highlights.':molecules[a].name+': '+x.donorSites+' donor site(s), '+x.acceptorSites+' acceptor site(s). '+molecules[b].name+': '+y.donorSites+' donor site(s), '+y.acceptorSites+' acceptor site(s).';
  ui.detail.textContent='Solid green circles = donor sites; dashed blue squares = acceptor sites. Count sites, not donor H atoms. '+boiling+'. These supplied property examples are not calculated from site count. Carboxylic-acid OH is a donor; its protonated hydroxyl O is not counted as a useful acceptor here.';
  state(ui,{first:a,second:b,highlight:show,firstSites:x,secondSites:y});}update();}
function mount03(){const ui=initDemo('day03');if(!ui)return;let pH=4.8,pKa=4.8;
 const ph=rangeControl(ui.controls,'Aqueous pH',0,14,.1,pH,'',v=>{pH=v;update();});
 const pk=rangeControl(ui.controls,'Supplied acid pKa',2,10,.1,pKa,'',v=>{pKa=v;update();});
 button(ui.controls,'Set pH = pKa',()=>setRange(ph,pKa));
 button(ui.controls,'One pH unit above pKa',()=>setRange(ph,Math.min(14,pKa+1)));
 button(ui.controls,'Reset',()=>{setRange(pk,4.8);setRange(ph,4.8);});
 function update(){const f=acidFraction(pH,pKa),points=Array.from({length:141},(_,i)=>[i/10,acidFraction(i/10,pKa)]);
  ui.visual.innerHTML=smallPlot('Predicted aqueous ionisation',[{points}],[0,14],[0,1],'pH','Fraction A⁻',[pH,f])+svg('Acid and conjugate base',text(240,65,'HA ⇌ H⁺ + A⁻',32)+text(240,120,'At the marker: '+fmt(f*100)+'% A⁻',27),480,150);
  ui.result.textContent='[A⁻]/[HA] = '+fmt(10**(pH-pKa))+'. Fraction A⁻ = '+fmt(f,4)+'; fraction HA = '+fmt(1-f,4)+'.';
  ui.detail.textContent='The curve is calculated, not measured. A monoprotic acid in a dilute aqueous approximation is assumed. At pH = pKa the two concentrations are equal. Fractions describe a molecular population, not a fraction of a proton on every molecule.';
  state(ui,{pH,pKa,fractionAminus:f,fractionHA:1-f});}update();}
function mount04(){const ui=initDemo('day04');if(!ui)return;let key='lactic',group='alcohol',reveal=false;
 const mol=selectControl(ui.controls,'Molecule to inspect',moleculeOptions(['lactic','glycine','acid','ethanol','ether','propanone']),key,v=>{key=v;reveal=false;update();});
 const g=selectControl(ui.controls,'Functional group to locate',['alcohol','ether','ketone','carboxylic acid','amine'],group,v=>{group=v;reveal=false;update();});
 button(ui.controls,'Highlight the selected pattern',()=>{reveal=true;update();});
 button(ui.controls,'Hide highlight',()=>{reveal=false;update();});
 button(ui.controls,'Reset',()=>{key='lactic';group='alcohol';reveal=false;mol.value=key;g.value=group;update();});
 const patterns={alcohol:'An O–H attached to a saturated carbon is the alcohol site here; acid OH is part of a different motif.',ether:'A simple ether has C–O–C without that O being attached to an acyl carbon.',ketone:'A ketone carbonyl carbon is attached to two carbon groups.', 'carboxylic acid':'The complete acid motif is C(=O)OH, including both oxygens.',amine:'This neutral amine site has an available N pair and no directly attached carbonyl.'};
 function update(){const m=groupMatch(key,group);ui.visual.innerHTML=chemicalDrawing(key,{highlight:reveal?m.atomIndices:[]});
  ui.result.textContent=reveal?(m.present?'The selected '+group+' pattern is highlighted.':'This drawn molecule has no '+group+' pattern.'):'Locate the complete pattern before revealing. One oxygen-containing group does not classify every O atom the same way.';
  ui.detail.textContent=patterns[group]+(reveal?' Other groups present: '+(m.groups.join(', ')||'none from this set')+'.':'');
  state(ui,{molecule:key,selectedGroup:group,revealed:reveal,present:m.present,highlightedAtoms:reveal?m.atomIndices:[]});}update();}
function mount05(){const ui=initDemo('day05');if(!ui)return;let pH=4.8,pKa=4.8,P=10,Vaq=50,Vorg=50;
 const inputs={};
 inputs.ph=rangeControl(ui.controls,'Fixed aqueous pH',0,12,.1,pH,'',v=>{pH=v;update();});
 inputs.pka=rangeControl(ui.controls,'Supplied aqueous pKa',2,9,.1,pKa,'',v=>{pKa=v;update();});
 inputs.p=rangeControl(ui.controls,'Neutral-species partition coefficient P',.1,20,.1,P,'',v=>{P=v;update();});
 inputs.aq=rangeControl(ui.controls,'Aqueous phase volume',5,100,5,Vaq,'mL',v=>{Vaq=v;update();});
 inputs.org=rangeControl(ui.controls,'Organic phase volume',5,100,5,Vorg,'mL',v=>{Vorg=v;update();});
 button(ui.controls,'Mostly neutral acid (pH 2)',()=>setRange(inputs.ph,2));
 button(ui.controls,'Mostly ionised acid (pH 8)',()=>setRange(inputs.ph,8));
 button(ui.controls,'Reset',()=>{setRange(inputs.pka,4.8);setRange(inputs.p,10);setRange(inputs.aq,50);setRange(inputs.org,50);setRange(inputs.ph,4.8);});
 function update(){const m=extractionModel(pH,pKa,P,Vaq,Vorg),oh=100,top=65,sep=top+oh;
  let b=`<rect x="75" y="65" width="330" height="200" rx="9" fill="#dbeff5" stroke="${C.ink}" stroke-width="3"/><rect x="77" y="67" width="326" height="${Math.max(1,oh-2)}" fill="#e6e6f4"/>`;
  b+=line(75,sep,405,sep,C.blue,2)+text(240,Math.max(93,top+oh/2),'Organic: HA only',25)+text(240,Math.min(246,sep+(200-oh)/2),'Aqueous: HA + A⁻',25);
  b+=text(240,310,'HA ⇌ H⁺ + A⁻ in water',28)+text(240,355,'Organic solute: '+fmt(m.organicFraction*100)+'%',28);
  ui.visual.innerHTML=svg('Two-phase neutral-only partition model',b,480,390);
  ui.result.textContent='D = '+fmt(m.D,4)+'. Total solute fractions: organic '+fmt(m.organicFraction,4)+', aqueous HA '+fmt(m.aqueousNeutralFraction,4)+', aqueous A⁻ '+fmt(m.aqueousIonisedFraction,4)+'.';
  ui.detail.textContent='P = [HA]organic/[HA]aqueous. Only neutral HA partitions into the organic phase; A⁻ stays in water. The aqueous pH is externally fixed, the phases are immiscible, volumes stay constant and equilibrium is reached. D = P/(1 + 10^(pH − pKa)); organic fraction = D Vorg/(Vaq + D Vorg). Layer ordering and equal layer heights are schematic, not density or phase-volume predictions.';
  state(ui,{pH,pKa,P,Vaq,Vorg,...m});}update();}
function mount06(){const ui=initDemo('day06');if(!ui)return;let angle=180;
 const slider=rangeControl(ui.controls,'C1–C2–C3–C4 methyl dihedral',0,360,5,angle,'degrees',v=>{angle=v;update();});
 button(ui.controls,'Anti (180°)',()=>setRange(slider,180));button(ui.controls,'Gauche (60°)',()=>setRange(slider,60));button(ui.controls,'Eclipsed (0°)',()=>setRange(slider,0));button(ui.controls,'Reset',()=>setRange(slider,180));
 function update(){const m=butaneModel(angle);let b=circle(240,220,38,'white',C.ink);
  for(const [offset,rear,col] of [[angle,true,C.ink],[0,false,C.blue]])for(let j=0;j<3;j++){
   const a=(-90+offset+j*120)*Math.PI/180,start=rear?38:0,end=rear?140:105,rr=rear?175:130;
   b+=line(240+start*Math.cos(a),220+start*Math.sin(a),240+end*Math.cos(a),220+end*Math.sin(a),col,3);
   b+=text(240+rr*Math.cos(a),220+rr*Math.sin(a)+8,j===0?'CH₃':'H',26,col);
  }
  b+=circle(240,220,4,C.blue)+text(240,445,m.conformation,25);
  ui.visual.innerHTML=svg('Butane Newman projection viewed C2 toward C3',b,480,480)+smallPlot('Qualitative butane torsional energy',[{points:Array.from({length:181},(_,i)=>[i*2,butaneModel(i*2).energy])}],[0,360],[0,2],'Methyl dihedral / degrees','Schematic relative energy',[m.angle,m.energy]);
  ui.result.textContent=m.conformation+'. Relative illustrative energy '+fmt(m.energy)+' (arbitrary units). Front carbon = blue bonds from central dot; rear carbon = dark bonds from circle.';
  ui.detail.textContent='The viewing direction is C2 → C3. The two CH₃ groups are C1 and C4. Turning the rear set preserves connectivity. The energy curve is a smooth qualitative interpolation through canonical arrangements, not measured free energies, a population calculation or a kinetic model.';
  state(ui,{inputAngle:angle,...m});}update();}
function mount07(){const ui=initDemo('day07');if(!ui)return;let groups=initialFischer.slice(),operation='Initial lactic-acid Fischer layout.';
 const pairs=[[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]],names=['top','right','bottom','left'];let chosen=0;
 const sel=selectControl(ui.controls,'Groups to interchange',pairs.map((p,i)=>[String(i),names[p[0]]+' and '+names[p[1]]]),'0',v=>{chosen=Number(v);});
 button(ui.controls,'Swap selected groups',()=>{const p=pairs[chosen];groups=swapFischer(groups,...p);operation='One pair was interchanged. This inverts the implied configuration.';update();});
 button(ui.controls,'Mirror left and right',()=>{groups=swapFischer(groups,1,3);operation='Left/right reflection swaps one pair, giving the enantiomer.';update();});
 button(ui.controls,'Rotate Fischer layout 180°',()=>{groups=rotateFischer(groups,2);operation='A 180° in-plane Fischer rotation preserves the implied configuration.';update();});
 button(ui.controls,'Rotate Fischer layout 90°',()=>{groups=rotateFischer(groups,1);operation='A 90° page-layout rotation changes the implied Fischer configuration; it is not a rigid 3D rotation of the molecule.';update();});
 button(ui.controls,'Reset',()=>{groups=initialFischer.slice();chosen=0;sel.value='0';operation='Initial lactic-acid Fischer layout.';update();});
 function update(){const config=fischerConfiguration(groups),priority={OH:1,COOH:2,'CH₃':3,H:4},pos=[[240,63],[390,178],[240,293],[90,178]];
  let b=line(145,170,335,170)+line(240,90,240,260);
  groups.forEach((g,i)=>{b+=text(pos[i][0],pos[i][1],g+' ('+priority[g]+')',28);});
  b+=text(240,355,'Current configuration: '+config,30);
  ui.visual.innerHTML=svg('Interactive lactic-acid Fischer projection',b,480,390);
  ui.result.textContent=operation+' Current assignment: '+config+'.';
  ui.detail.textContent='Horizontal bonds point toward you; vertical bonds point away. Priorities: OH > COOH > CH₃ > H. For the carbon tie, COOH compares as O,O,O and CH₃ as H,H,H. Clockwise 1→2→3 gives R only when group 4 points away; invert when it points toward. No optical-rotation sign is predicted.';
  state(ui,{groups:groups.slice(),configuration:config,operation});}update();}
function curved(d,color){return `<path d="${d}" fill="none" stroke="${color}" stroke-width="3" marker-end="url(#chem-arrow)"/>`;}
function mechanismDrawing(key,a,b,checked){
 const col=i=>checked&&i!==0?C.red:C.blue;let s='';
 if(key==='acid'){
  s=text(120,165,'CH₃C(=O)O',28)+line(225,155,277,155)+text(300,165,'H',29)+text(400,165,'OH⁻',29)+circle(356,142,2,C.blue)+circle(356,153,2,C.blue);
  s+=curved(a===0?'M356 142 Q330 55 301 140':a===1?'M300 140 Q330 55 375 140':'M356 142 Q230 35 55 140',col(a));
  s+=curved(b===0?'M250 155 Q250 230 207 180':b===1?'M250 155 Q275 210 300 180':'M205 140 Q295 30 380 140',col(b));
 }else if(key==='carbonyl'){
  s=text(180,165,'C',30)+text(180,65,'O',30)+text(65,165,'CH₃',28)+text(280,165,'H',29)+line(95,155,153,155)+line(205,155,255,155)+line(174,125,174,85)+line(186,125,186,85)+text(400,165,'OH⁻',29)+circle(356,142,2,C.blue)+circle(356,153,2,C.blue);
  s+=curved(a===0?'M356 142 Q290 65 195 145':a===1?'M195 145 Q290 65 375 140':'M356 142 Q240 30 70 137',col(a));
  if(b!==1)s+=curved(b===0?'M186 108 Q240 80 190 68':'M186 108 Q280 75 281 138',col(b));
  else s+=text(240,220,'π bond incorrectly retained',23,C.red);
 }else{
  s=text(110,165,'CH₃NH₂',29)+circle(187,142,2,C.blue)+circle(187,153,2,C.blue)+text(290,165,'H',30)+line(315,155,365,155)+text(410,165,'OH₂⁺',29);
  s+=curved(a===0?'M187 142 Q240 60 287 140':a===1?'M287 140 Q240 60 185 140':'M187 142 Q290 30 390 140',col(a));
  s+=curved(b===0?'M340 155 Q365 230 395 180':b===1?'M340 155 Q315 220 290 180':'M185 155 Q150 235 115 195',col(b));
 }
 s+=text(240,275,checked&&a===0&&b===0?mechanisms[key].product:'Proposed pair movements',29)+text(240,345,'Blue = proposed; red = a checked error',23);
 return svg('Supplied elementary step with proposed electron-pair arrows',s,480,375);
}
function mount08(){const ui=initDemo('day08');if(!ui)return;let key='acid',a=1,b=1,charge=0,checked=false;
 const mol=selectControl(ui.controls,'Elementary step',Object.entries(mechanisms).map(([k,v])=>[k,v.name]),key,v=>{key=v;a=1;b=1;charge=0;checked=false;rebuild();update();});
 const extra=document.createElement('div');ui.controls.append(extra);let first,second,q;
 function rebuild(){extra.innerHTML='';first=selectControl(extra,'First proposed electron-pair arrow',mechanisms[key].first.map((v,i)=>[String(i),v]),String(a),v=>{a=Number(v);checked=false;update();});second=selectControl(extra,'Second proposed electron-pair arrow',mechanisms[key].second.map((v,i)=>[String(i),v]),String(b),v=>{b=Number(v);checked=false;update();});q=selectControl(extra,'Proposed total product charge',[['-1','−1'],['0','0'],['1','+1']],String(charge),v=>{charge=Number(v);checked=false;update();});}
 button(ui.controls,'Check arrows and charge',()=>{checked=true;update();});
 button(ui.controls,'Show a valid elementary step',()=>{a=0;b=0;charge=mechanisms[key].netCharge;checked=true;rebuild();update();});
 button(ui.controls,'Clear proposed arrows',()=>{a=1;b=1;charge=0;checked=false;rebuild();update();});
 button(ui.controls,'Reset',()=>{key='acid';mol.value=key;a=1;b=1;charge=0;checked=false;rebuild();update();});
 function update(){const m=mechanismModel(key,a,b,charge);ui.visual.innerHTML=mechanismDrawing(key,a,b,checked);
  ui.result.textContent=checked?(m.correct?'The two pair movements and net charge are consistent.':(!m.validFirst?'Repair the first arrow: electron pairs start on a bond or lone pair. ':'')+(!m.validSecond?(key==='carbonyl'?(b===1?'Keeping C=O as the new bond forms gives five bond orders at carbon. ':'The carbonyl π pair must go onto O, not the attached H. '):'The donor bond pair must return to the atom losing H⁺. '):'')+(!m.validCharge?'The total product charge must be '+m.expectedCharge+'.':'')):'Predict the products, choose pair movements, and check the total product charge.';
  ui.detail.textContent=checked?mechanisms[key].reason:'Full-headed curved arrows describe electron-pair redistribution, not atom trajectories. The options include deliberately invalid proposals. An invalid proposal is not a claim that a stable product with impossible valence exists.';
  state(ui,{mechanism:key,firstArrow:a,secondArrow:b,proposedCharge:charge,checked,...m});}rebuild();update();}
function mount09(){const ui=initDemo('day09');if(!ui)return;let deltaG=3500,temperature=298.15,reference=0;
 const dg=rangeControl(ui.controls,'Free energy B minus A',-5000,10000,100,deltaG,'J/mol',v=>{deltaG=v;update();});
 const temp=rangeControl(ui.controls,'Temperature',250,350,.05,temperature,'K',v=>{temperature=v;update();});
 const ref=rangeControl(ui.controls,'Common free-energy reference offset',0,3000,100,reference,'J/mol',v=>{reference=v;update();});
 button(ui.controls,'Equal-state free energies',()=>setRange(dg,0));button(ui.controls,'Reset',()=>{setRange(dg,3500);setRange(temp,298.15);setRange(ref,0);});
 function update(){const m=conformerModel(deltaG,temperature,reference);let b='';
  [[130,m.fractionA,'A'],[350,m.fractionB,'B']].forEach(([x,p,label])=>{b+=`<rect x="${x-55}" y="${280-p*210}" width="110" height="${p*210}" fill="${label==='A'?C.blue:C.green}"/>`+text(x,315,label,29)+text(x,Math.max(45,260-p*210),fmt(p*100)+'%',26);});
  b+=text(240,355,'Predicted equilibrium population',26);
  ui.visual.innerHTML=svg('Two-state conformer free-energy populations',b,480,390)+svg('Supplied state free energies',text(240,55,'G(A) = '+fmt(reference)+' J/mol',27)+text(240,115,'G(B) = '+fmt(reference+deltaG)+' J/mol',27)+text(240,175,'ΔG(B − A) = '+fmt(deltaG)+' J/mol',27),480,210);
  ui.result.textContent='Fraction A = '+fmt(m.fractionA,4)+'; fraction B = '+fmt(m.fractionB,4)+'. Relative B/A weight = '+fmt(Math.exp(-deltaG/(R*temperature)),4)+'.';
  ui.detail.textContent='Supplied effective state free energies already include degeneracy. Use w = exp(−G/RT), then divide each weight by their sum. R = 8.314462618 J/(mol K). A common reference offset changes raw weights but not normalised fractions. This is an equilibrium model; no interconversion barrier or rate is supplied.';
  state(ui,m);}update();}
[mount01,mount02,mount03,mount04,mount05,mount06,mount07,mount08,mount09].forEach(mount=>mount());
})(typeof globalThis!=='undefined'?globalThis:window);
