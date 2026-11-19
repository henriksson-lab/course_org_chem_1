(function (global) {
'use strict';
const AMINES = {
 methylamine:{name:'Methylamine',neutral:'CH₃NH₂',acid:'CH₃NH₃⁺',pKa:10.6},
 aniline:{name:'Aniline',neutral:'PhNH₂',acid:'PhNH₃⁺',pKa:4.6},
 pyridine:{name:'Pyridine',neutral:'C₅H₅N',acid:'C₅H₅NH⁺',pKa:5.2}
};
const ACIDS = {
 ethanoic:{name:'Ethanoic acid',HA:'CH₃COOH',A:'CH₃COO⁻',pKa:4.8,charge:0},
 methylammonium:{name:'Methylammonium',HA:'CH₃NH₃⁺',A:'CH₃NH₂',pKa:10.6,charge:1},
 phenol:{name:'Phenol',HA:'C₆H₅OH',A:'C₆H₅O⁻',pKa:10.0,charge:0}
};
const BASES = {
 methylamine:{name:'Methylamine',B:'CH₃NH₂',BH:'CH₃NH₃⁺',pKa:10.6,charge:0},
 pyridine:{name:'Pyridine',B:'C₅H₅N',BH:'C₅H₅NH⁺',pKa:5.2,charge:0},
 phenoxide:{name:'Phenoxide',B:'C₆H₅O⁻',BH:'C₆H₅OH',pKa:10.0,charge:-1}
};
const near = (a,b,tol=1e-10)=>Math.abs(a-b)<tol;
function amineModel(molecule,pH) {
 const a=AMINES[molecule];
 if(!a||!Number.isFinite(pH)||pH<0||pH>14)return {valid:false,message:'Choose a supplied molecule and pH from 0 to 14.'};
 const protonated=1/(1+10**(pH-a.pKa));
 return {valid:true,molecule,pH,pKa:a.pKa,protonated,neutral:1-protonated,
  dominant:protonated>0.5?'protonated':protonated<0.5?'neutral':'equal',
  neutralFormula:a.neutral,protonatedFormula:a.acid};
}
function protonTransferModel(acidKey,baseKey) {
 const acid=ACIDS[acidKey],base=BASES[baseKey];
 if(!acid||!base)return {valid:false,message:'Select one supplied acid pair and base pair.'};
 const delta=base.pKa-acid.pKa,K=10**delta;
 const direction=delta>1?'products favoured':delta< -1?'reactants favoured':'neither side strongly favoured by K';
 return {valid:true,acidKey,baseKey,acid,base,delta,K,direction,
  totalCharge:acid.charge+base.charge,
  samePair:acid.HA===base.BH&&acid.A===base.B};
}
function electronModel(molecule,donor,acceptor,paired) {
 const isAmide=molecule==='amide',name=isAmide?'ethanamide':'ethylamine';
 if(!['amide','amine'].includes(molecule))return {valid:false,message:'Choose amine or amide.'};
 const base={molecule,name,donor,acceptor,paired:acceptor==='CN_pi'&&!!paired};
 if(donor==='atom_label')return {...base,valid:false,kind:'invalid',message:'An atom label or positive charge is not an electron pair. Start at an actual lone pair or bond.'};
 if(donor==='N_pair'&&acceptor==='CN_pi') {
  if(!isAmide)return {...base,valid:false,kind:'invalid',message:'Ethylamine has no neighbouring carbonyl to accept the coupled move. Making C=N while retaining both C–H bonds would overbond this carbon.'};
  if(!paired)return {...base,valid:false,kind:'invalid',message:'The N pair cannot form C=N while the C=O pi bond stays intact: carbon would have five bond orders. Move the C=O pi pair to O too.'};
  return {...base,valid:true,kind:'amide_resonance',charges:{N:1,O:-1,net:0},message:'Valid coupled contributor: N pair forms C=N and the C=O pi pair goes to O. N becomes +1, O −1; both drawings describe one delocalised amide.'};
 }
 if(donor==='N_pair'&&acceptor==='H_atom')return {...base,valid:true,kind:'protonation',charges:{N:1,net:1},message:isAmide?'This is a formal N-protonation reaction, not resonance. Amide N is much less basic than a simple amine; ordinary strong-acid protonation of an amide often favours O. The drawing alone does not predict the predominant site.':'N uses its lone pair to accept H⁺, giving four N bonds and +1 charge. This changes protonation state: it is a reaction, not a resonance contributor.'};
 if(donor==='C_O_pi'&&acceptor==='O_atom') {
  if(!isAmide)return {...base,valid:false,kind:'invalid',message:'Ethylamine has no C=O pi pair. Choose an electron source that exists in this molecule.'};
  return {...base,valid:true,kind:'polar_carbonyl',charges:{C:1,O:-1,N:0,net:0},message:'Moving C=O pi electrons to O gives a polar carbonyl contributor with C⁺ and O⁻. It preserves net charge but, without N donation, does not yet show the contributor explaining partial C–N double-bond character.'};
 }
 return {...base,valid:false,kind:'invalid',message:'That source/target combination is not one of the valid pair moves drawn here. For amide donation use N pair → C–N pi, coupled to C=O pi → O; for protonation use N pair → H⁺.'};
}
const SAMPLES={
 U1:{name:'Unknown U1',formula:'C₄H₈O₂',carbonyl:true,OH:false,
   nmr:'1.25 ppm t, 3H; 2.05 ppm s, 3H; 4.12 ppm q, 2H',
   peaks:[{ppm:1.25,n:3,H:3},{ppm:2.05,n:1,H:3},{ppm:4.12,n:4,H:2}]},
 U2:{name:'Unknown U2',formula:'C₂H₆O',carbonyl:false,OH:true,
   nmr:'1.18 ppm t, 3H; 3.65 ppm q, 2H; exchangeable broad OH, 1H',
   peaks:[{ppm:1.18,n:3,H:3},{ppm:3.65,n:4,H:2},{ppm:2.4,n:1,H:1,exchange:true}]}
};
const CANDIDATES={
 ethyl_ethanoate:{name:'Ethyl ethanoate',formula:'C₄H₈O₂',structure:'CH₃C(=O)OCH₂CH₃',carbonyl:true,OH:false,fits:'U1'},
 butanoic_acid:{name:'Butanoic acid',formula:'C₄H₈O₂',structure:'CH₃CH₂CH₂C(=O)OH',carbonyl:true,OH:true,fits:null},
 ethanol:{name:'Ethanol',formula:'C₂H₆O',structure:'CH₃CH₂OH',carbonyl:false,OH:true,fits:'U2'},
 dimethyl_ether:{name:'Dimethyl ether',formula:'C₂H₆O',structure:'CH₃OCH₃',carbonyl:false,OH:false,fits:null}
};
function spectralCandidateModel(sampleKey,candidateKey,evidence) {
 const sample=SAMPLES[sampleKey],candidate=CANDIDATES[candidateKey];
 if(!sample||!candidate)return {valid:false,message:'Select one supplied unknown and candidate.'};
 const level=['formula','IR','all'].includes(evidence)?evidence:'all';
 const tests=[{label:'Molecular formula',pass:sample.formula===candidate.formula}];
 if(level!=='formula') {
  tests.push({label:'Carbonyl IR pattern',pass:sample.carbonyl===candidate.carbonyl});
  tests.push({label:'O–H IR region',pass:sample.OH===candidate.OH});
 }
 if(level==='all')tests.push({label:'Full supplied proton-NMR pattern',pass:candidate.fits===sampleKey});
 return {valid:true,sampleKey,candidateKey,evidence:level,sample,candidate,tests,
  fits:tests.every(t=>t.pass),conflicts:tests.filter(t=>!t.pass).map(t=>t.label)};
}
function tlcModel(front,spot,recorded,composition) {
 const geometric=Number.isFinite(front)&&Number.isFinite(spot)&&front>0&&spot>=0&&spot<=front;
 return {valid:geometric,front,spot,recorded:!!recorded,composition,
  Rf:geometric&&recorded?spot/front:null,
  componentCount:composition==='co_migrate'?2:1,visibleFeatures:1,
  message:!geometric?'Distances must obey 0 ≤ spot ≤ front, with front > 0.':!recorded?'Rf is unavailable: the actual solvent front was not recorded.':'Rf compares two distances under these conditions; it is not a unique identity code.'};
}
const NODES={
 ethanol:{name:'Ethanol',formula:'CH₃CH₂OH',group:'primary alcohol'},
 ethanal:{name:'Ethanal',formula:'CH₃CHO',group:'aldehyde'},
 acid:{name:'Ethanoic acid',formula:'CH₃COOH',group:'carboxylic acid'},
 ester:{name:'Ethyl ethanoate',formula:'CH₃COOCH₂CH₃',group:'ester'},
 carboxylate:{name:'Ethanoate',formula:'CH₃COO⁻',group:'carboxylate'}
};
const ACTIONS={
 oxidize_aldehyde:'Controlled oxidation to aldehyde',oxidize_acid:'Oxidation to acid',
 reduce:'Selected carbonyl reduction, then proton donor',esterify:'Esterification with a separate ethanol portion',
 hydrolyse:'Base-promoted ester hydrolysis',workup:'Acid work-up'};
function routeStep(node,action) {
 const key=node+':'+action;
 const edges={
 'ethanol:oxidize_aldehyde':['ethanal','Alcohol → aldehyde; conditions must limit further oxidation.','Loss of alcohol O–H and appearance of carbonyl would support this functional-group change.'],
 'ethanol:oxidize_acid':['acid','Primary alcohol → acid in the selected overall oxidation map.','Acid-group evidence plus full structure checks would test the target.'],
 'ethanal:oxidize_acid':['acid','Aldehyde → carboxylic acid, retaining the carbon skeleton.','Compare complete acid/carbonyl evidence, not only the starting formula.'],
 'ethanal:reduce':['ethanol','Carbonyl reduction gives alkoxide, then a proton donor gives alcohol.','Appearance of alcohol O–H and loss of C=O are consistent; identity still needs other evidence.'],
 'acid:esterify':['ester','Acid + a separate ethanol portion ⇌ ester + water under suitable acid catalysis.','An ethyl ester pattern should be tested; equilibrium and isolation losses limit the amount.'],
 'ester:hydrolyse':['carboxylate','Ester + OH⁻ → carboxylate + ethanol; hydroxide is consumed overall.','A basic mixture is represented as carboxylate, not neutral acid.'],
 'carboxylate:workup':['acid','Carboxylate + H₃O⁺ → acid + water. This is proton transfer, not a new carbon skeleton.','Recorded protonation conditions distinguish the acid work-up from the preceding basic mixture.']
 };
 const edge=edges[key];
 const equations={
 'ethanol:oxidize_aldehyde':['CH₃CH₂OH + [O]','CH₃CHO + H₂O'],
 'ethanol:oxidize_acid':['CH₃CH₂OH + 2[O]','CH₃COOH + H₂O'],
 'ethanal:oxidize_acid':['CH₃CHO + [O]','CH₃COOH'],
 'ethanal:reduce':['CH₃CHO + H:⁻ + H₂O','CH₃CH₂OH + OH⁻'],
 'acid:esterify':['CH₃COOH + CH₃CH₂OH','CH₃COOCH₂CH₃ + H₂O'],
 'ester:hydrolyse':['CH₃COOCH₂CH₃ + OH⁻','CH₃COO⁻ + CH₃CH₂OH'],
 'carboxylate:workup':['CH₃COO⁻ + H₃O⁺','CH₃COOH + H₂O']
 };
 return edge?{valid:true,node,action,next:edge[0],reason:edge[1],evidence:edge[2],equation:equations[key]}:
 {valid:false,node,action,next:node,reason:'The supplied small toolkit does not support that transformation at this starting functional group. Choose a change matching the current structure.',evidence:'An unsupported arrow is not repaired by a plausible target name.'};
}
function routeModel(nodes,target,last) {
 return {nodes:nodes.slice(),target,current:nodes[nodes.length-1],achieved:nodes[nodes.length-1]===target,
  steps:nodes.length-1,last:last||null};
}
function evidenceStrengthModel(IR,NMR,claim) {
 const fullIR=IR==='full',fullNMR=NMR==='full',integrals=fullNMR||NMR==='relative';
 const outcomes={
 ester_best:{supported:fullNMR,reason:fullNMR?'The full formula-consistent proton pattern favours the ester among the supplied candidates; IR can add an independent check.':'Formula alone or partial clues do not distinguish all supplied candidates. Complete shifts, integrations and splitting provide the discriminating pattern.'},
 acid_challenged:{supported:fullIR||fullNMR,reason:fullIR?'The complete, adequately sensitive IR covers the acid O–H region and lacks the expected broad feature. This challenges the acid candidate.':fullNMR?'The full proton pattern gives O–CH₂ plus acetyl methyl rather than the acid candidate’s pattern.':'A carbonyl-only IR window cannot establish that O–H is absent; fragmentary NMR does not complete the candidate test.'},
 OH_absent:{supported:fullIR,reason:fullIR?'The relevant O–H region was measured with stated adequate sensitivity, so its absence can be used cautiously.':'Unmeasured or unavailable O–H data are missing evidence, not an observed absence.'},
 integrals_eight:{supported:integrals,reason:integrals?'Reliable relative areas 3:3:2 normalise to the eight H of the supplied formula. Areas do not alone locate every group.':'Peak line counts or heights are not reliable integrated H areas. The formula supplies eight total H, but these incomplete NMR data do not account for them by integration.'},
 OCH2_assigned:{supported:fullNMR,reason:fullNMR?'A 2H quartet near 4.12 ppm with a 3H triplet supports the O–CH₂–CH₃ assignment in this candidate comparison.':'A ratio without chemical shifts or a splitting pattern without integration is insufficient for this complete assignment.'},
 purity:{supported:false,reason:'These selected qualitative clues do not establish sample purity. A mixture can preserve familiar features; quantitative composition requires suitable independent evidence.'},
 unique:{supported:false,reason:'The supplied candidate comparison can favour a structure without excluding every imaginable alternative. “Proves the only possible structure” is too strong.'}
 };
 const outcome=outcomes[claim]||outcomes.unique;
 return {IR,NMR,claim,supported:outcome.supported,reason:outcome.reason,fullIR,fullNMR,integrals};
}
function integratedModel(mass,source,pattern,spots,conversionData,claim) {
 const valid=Number.isFinite(mass)&&mass>=0,theoretical=0.0100*180.16;
 const apparentYield=valid?100*mass/theoretical:null;
 const sameSample=source==='same',identitySupport=sameSample&&pattern==='consistent';
 const conversion=conversionData==='residual2'?80:null;
 const inconsistent=conversion!==null&&apparentYield>conversion+0.2;
 const outcomes={
 yield:{supported:valid,reason:'The mass calculation gives an apparent isolated yield under the product-mass assumption. It does not by itself establish identity, purity or conversion.'},
 conversion:{supported:conversion!==null,reason:conversion!==null?'The supplied quantitative same-batch residual assay reports 2.00 mmol starting material from 10.00 mmol initially: conversion = 80.0%.':'Isolated mass is not a starting-material assay. Without a supplied residual-amount measurement, conversion is unknown.'},
 identity:{supported:identitySupport,reason:identitySupport?'A consistent spectrum explicitly linked to this isolate supports the target assignment. Treat it as supporting evidence, not a universal proof of purity.':source==='unrelated'?'An unrelated simulated spectrum is not a measurement of this isolate and cannot identify it.':source==='none'?'No sample-linked structural spectrum is supplied. A mass or TLC position alone cannot identify the isolate.':'Contradictory major signals challenge the assumed product identity or composition; investigate before claiming success.'},
 purity:{supported:false,reason:spots===1?'One visible TLC feature can hide co-migration or undetected components. Neither this spot count nor mass proves purity.':'Two detected TLC features give evidence of multiple components under the stated conditions; a pure-product claim is unsupported.'},
 success:{supported:false,reason:'Separate the claims: calculate amount, test sample-linked identity, assess composition and measure conversion independently. The combined data do not establish pure product with known conversion merely from isolated yield.'}
 };
 const out=outcomes[claim]||outcomes.success;
 return {valid,mass,theoretical,apparentYield,source,pattern,spots,conversionData,conversion,
  claim,supported:out.supported,reason:out.reason,sameSample,identitySupport,inconsistent};
}
function checks() {
 const out=[];const check=(name,pass)=>out.push({name,pass:!!pass});
 check('Amine at pH = pKa is half protonated',near(amineModel('methylamine',10.6).protonated,.5));
 check('Higher pH decreases ammonium fraction',amineModel('aniline',8).protonated<amineModel('aniline',3).protonated);
 check('One pH unit below pKa gives 10/11 protonated',near(amineModel('pyridine',4.2).protonated,10/11));
 check('Invalid pH is rejected',!amineModel('pyridine',20).valid);
 const transfer=protonTransferModel('ethanoic','methylamine');
 check('pKa difference gives K near 10^5.8',near(transfer.delta,5.8)&&near(transfer.K/10**5.8,1));
 check('Same acid/base pair has K = 1',near(protonTransferModel('methylammonium','methylamine').K,1));
 check('Weak acid versus weak base can favour reactants',protonTransferModel('phenol','pyridine').K<1);
 check('Amide contributor preserves net zero',electronModel('amide','N_pair','CN_pi',true).charges.net===0);
 check('Unpaired N donation would overbond carbon',!electronModel('amide','N_pair','CN_pi',false).valid);
 check('Amine has no carbonyl pi pair',!electronModel('amine','C_O_pi','O_atom',false).valid);
 check('Protonation is not resonance',electronModel('amine','N_pair','H_atom',false).kind==='protonation');
 check('U1 ester fits full supplied data',spectralCandidateModel('U1','ethyl_ethanoate','all').fits);
 check('Formula alone cannot reject the U1 acid',spectralCandidateModel('U1','butanoic_acid','formula').fits);
 check('Complete U1 data reject acid',!spectralCandidateModel('U1','butanoic_acid','all').fits);
 check('U2 ethanol fits; ether OH conflicts',spectralCandidateModel('U2','ethanol','all').fits&&!spectralCandidateModel('U2','dimethyl_ether','IR').fits);
 check('TLC Rf = 3.2 / 5.0',near(tlcModel(5,3.2,true,'single').Rf,.64));
 check('Missing TLC front does not yield an Rf',tlcModel(5,3.2,false,'single').Rf===null);
 check('Out-of-range spot rejected',!tlcModel(3,4,true,'single').valid);
 check('Co-migrating components can give one feature',tlcModel(5,3.2,true,'co_migrate').componentCount===2&&tlcModel(5,3.2,true,'co_migrate').visibleFeatures===1);
 check('Ethanol can lead to acid then ester',routeStep('ethanol','oxidize_acid').next==='acid'&&routeStep('acid','esterify').next==='ester');
 check('Base hydrolysis gives carboxylate',routeStep('ester','hydrolyse').next==='carboxylate');
 check('Esterification supplies ethanol and water explicitly',routeStep('acid','esterify').equation[0].includes('CH₃CH₂OH')&&routeStep('acid','esterify').equation[1].includes('H₂O'));
 check('Unsupported route does not change current structure',!routeStep('ethanol','workup').valid&&routeStep('ethanol','workup').next==='ethanol');
 check('Missing O-H window is not an observed absence',!evidenceStrengthModel('carbonyl','none','OH_absent').supported);
 check('Relative areas can account for eight H',evidenceStrengthModel('none','relative','integrals_eight').supported);
 check('Qualitative spectra do not prove purity',!evidenceStrengthModel('full','full','purity').supported);
 const report=integratedModel(1.20,'same','consistent',1,'none','yield');
 check('Teaching mass gives 66.6 percent apparent yield',Math.abs(report.apparentYield-66.60746)<.0001);
 check('Isolated mass alone does not supply conversion',integratedModel(1.20,'same','consistent',1,'none','conversion').conversion===null);
 check('Residual two of ten mmol gives 80 percent conversion',integratedModel(1.20,'same','consistent',1,'residual2','conversion').conversion===80);
 check('Unrelated simulated spectrum cannot identify isolate',!integratedModel(1.20,'unrelated','consistent',1,'none','identity').supported);
 check('Contradictory spectrum challenges identity',!integratedModel(1.20,'same','contradictory',1,'none','identity').supported);
 check('Apparent yield above measured conversion triggers audit',integratedModel(1.70,'same','consistent',1,'residual2','yield').inconsistent);
 check('Negative isolated mass rejected',!integratedModel(-1,'same','consistent',1,'none','yield').valid);
 check('Blank or nonfinite mass rejected',!integratedModel(NaN,'same','consistent',1,'none','yield').valid);
 return out;
}
const api={AMINES,ACIDS,BASES,SAMPLES,CANDIDATES,NODES,ACTIONS,amineModel,protonTransferModel,electronModel,
 spectralCandidateModel,tlcModel,routeStep,routeModel,evidenceStrengthModel,integratedModel,checks};
if(typeof module!=='undefined'&&module.exports)module.exports=api;
global.ChemApplications=api;
const C=global.ChemDemo;
if(typeof document==='undefined'||!C)return;
const T=C.text,L=C.line,O=C.circle;
function arrowPath(d,color='#176b91') {return `<path d="${d}" fill="none" stroke="${color}" stroke-width="3" marker-end="url(#chem-arrow)"/>`;}
function resetRange(input,value) {input.value=value;input.dispatchEvent(new Event('input',{bubbles:true}));}
function ringBody(x,y,{nitrogen=false,protonated=false,r=55}={}) {
 const points=Array.from({length:6},(_,i)=>[x+r*Math.cos(-Math.PI/2+i*Math.PI/3),y+r*Math.sin(-Math.PI/2+i*Math.PI/3)]);
 let b='';for(let i=0;i<6;i++) {
  let a=points[i].slice(),z=points[(i+1)%6].slice();
  if(nitrogen&&i===0)a=[a[0]+(z[0]-a[0])*.2,a[1]+(z[1]-a[1])*.2];
  if(nitrogen&&i===5)z=[z[0]+(a[0]-z[0])*.2,z[1]+(a[1]-z[1])*.2];
  b+=L(...a,...z);
 }
 b+=O(x,y,r*.5,'none','#176b91');
 if(nitrogen)b+=T(x,y-r+8,protonated?'NH⁺':'N',24);
 return b;
}
function amineStructure(key,x,y,protonated) {
 if(key==='methylamine')return T(x-65,y,'CH₃',29)+L(x-30,y-9,x+20,y-9)+T(x+65,y,protonated?'NH₃⁺':'NH₂',29)+(!protonated?T(x+65,y-40,'••',22):'');
 if(key==='aniline')return ringBody(x,y+15)+L(x,y-40,x,y-70)+T(x,y-85,protonated?'NH₃⁺':'NH₂',28)+(!protonated?T(x,y-118,'••',22):'');
 return ringBody(x,y+20,{nitrogen:true,protonated})+(!protonated?T(x,y-63,'••',22):'');
}
function irPlot(sampleKey,coverage='full') {
 const s=SAMPLES[sampleKey],points=[];
 for(let wn=4000;wn>=600;wn-=10) {
  let y=94;
  if(s.carbonyl)y-=60*Math.exp(-.5*((wn-1740)/25)**2)+25*Math.exp(-.5*((wn-1240)/50)**2);
  if(s.OH)y-=48*Math.exp(-.5*((wn-3330)/165)**2)+30*Math.exp(-.5*((wn-1050)/40)**2);
  points.push([wn,y]);
 }
 let result=C.plot({title:'Simulated infrared spectrum; not a measured sample',series:[{points}],xDomain:[4000,600],yDomain:[0,100],xLabel:'Wavenumber / cm⁻¹',yLabel:'Transmittance / %'});
 if(coverage==='carbonyl') {
  const left=85,right=680,x=wn=>left+(4000-wn)/3400*(right-left);
  const mask=`<rect x="${left}" y="45" width="${x(1850)-left}" height="265" fill="#e5ebef"/>`+
    `<rect x="${x(1600)}" y="45" width="${right-x(1600)}" height="265" fill="#e5ebef"/>`+
    T(260,140,'Unmeasured region',24)+T(260,185,'Do not infer O–H absence',21);
  result=result.replace('</svg>',mask+'</svg>');
 }
 return result;
}
function nmrPlot(sampleKey,{mode='full',extra=false}={}) {
 const s=SAMPLES[sampleKey];let b=L(75,270,675,270);
 if(mode==='relative') {
  [3,3,2].forEach((n,j)=>{b+=`<rect x="${150+j*180}" y="${260-n*45}" width="95" height="${n*45}" fill="#176b91"/>`+T(197+j*180,300,'Group '+String.fromCharCode(65+j),21)+T(197+j*180,95,n+' area units',21);});
  b+=T(375,350,'Relative integrated areas only; no shifts or splitting supplied',20);
 } else {
  const x=p=>75+(5-p)/5*600;
  s.peaks.forEach((p,j)=>{
   const h=p.n===4?[1,3,3,1]:p.n===3?[1,2,1]:[1];
   if(p.exchange) {
    const centre=mode==='multiplicity'?180+j*180:x(p.ppm),points=[];
    for(let k=0;k<=80;k++) {const dx=-55+k*110/80;points.push(`${centre+dx},${270-65*Math.exp(-.5*(dx/19)**2)}`);}
    b+=`<polyline points="${points.join(' ')}" fill="none" stroke="#176b91" stroke-width="3"/>`;
   } else h.forEach((v,k)=>{const xx=mode==='multiplicity'?180+j*180+(k-(p.n-1)/2)*12:x(p.ppm+(k-(p.n-1)/2)*.06);b+=L(xx,270,xx,270-v/Math.max(...h)*130,'#176b91',3);});
   if(mode==='full')b+=T(x(p.ppm),p.exchange?35:j%2?70:105,`${p.ppm} ppm; ${p.H}H${p.exchange?' (exchangeable)':''}`,19);
  });
  if(extra)b+=L(390,270,390,110,'#af3b45',4)+T(420,90,'Extra major group',21,'#af3b45');
  if(mode==='full')for(let p=0;p<=5;p++)b+=T(x(p),305,p,20);
  b+=T(375,350,mode==='multiplicity'?'Splitting groups only; shifts and areas unavailable':'Chemical shift / ppm; supplied areas, not peak heights',20);
 }
 return C.svg('Simulated proton-NMR teaching evidence; not a measured isolate',b,720,380);
}
function mount18() {
 const ui=C.initDemo('day18');if(!ui)return;
 let molecule='methylamine',pH=7;
 const mol=C.selectControl(ui.controls,'Nitrogen base',Object.keys(AMINES).map(k=>[k,AMINES[k].name]),molecule,v=>{molecule=v;update();});
 const ph=C.rangeControl(ui.controls,'Aqueous pH',0,14,.1,pH,'',v=>{pH=v;update();});
 C.button(ui.controls,'At supplied pKa',()=>resetRange(ph,AMINES[molecule].pKa));
 C.button(ui.controls,'One pH unit below pKa',()=>resetRange(ph,AMINES[molecule].pKa-1));
 C.button(ui.controls,'Reset',()=>{molecule='methylamine';mol.value=molecule;resetRange(ph,7);});
 function update() {
  const m=amineModel(molecule,pH);C.state(ui,m);
  let body=amineStructure(molecule,180,170,false)+amineStructure(molecule,540,170,true)+T(180,290,'B: neutral',25)+T(540,290,'BH⁺: protonated',25);
  body+=`<rect x="90" y="325" width="${540*m.neutral}" height="28" fill="#176b91"/><rect x="${90+540*m.neutral}" y="325" width="${540*m.protonated}" height="28" fill="#23765d"/>`;
  ui.visual.innerHTML=C.svg('Neutral and protonated structures with their calculated population fractions',body,720,380)+
   C.plot({title:'Calculated dilute aqueous amine protonation model',series:[{points:Array.from({length:141},(_,j)=>[j/10,amineModel(molecule,j/10).protonated])}],xDomain:[0,14],yDomain:[0,1],xLabel:'pH',yLabel:'Fraction BH⁺',markers:[{x:pH,y:m.protonated}]});
  ui.result.textContent=`${C.fmt(100*m.protonated,4)}% BH⁺; ${C.fmt(100*m.neutral,4)}% neutral B. Supplied pKa(BH⁺) = ${m.pKa}.`;
  ui.detail.textContent='Fraction BH⁺ = 1 / (1 + 10^(pH − pKa)). Approximate compatible aqueous values and a dilute one-site model are supplied. Fractions describe populations, not partly protonated individual molecules. This predicts neither an extraction recovery nor a nucleophilic reaction rate.';
 }
 update();
}
function mount19() {
 const ui=C.initDemo('day19');if(!ui)return;
 let acid='ethanoic',base='methylamine';
 const a=C.selectControl(ui.controls,'Proton donor HA',Object.keys(ACIDS).map(k=>[k,ACIDS[k].name]),acid,v=>{acid=v;update();});
 const b=C.selectControl(ui.controls,'Proton acceptor B',Object.keys(BASES).map(k=>[k,BASES[k].name]),base,v=>{base=v;update();});
 C.button(ui.controls,'Equal pair (K = 1)',()=>{acid='methylammonium';base='methylamine';a.value=acid;b.value=base;update();});
 C.button(ui.controls,'Reactant-favoured case',()=>{acid='phenol';base='pyridine';a.value=acid;b.value=base;update();});
 C.button(ui.controls,'Reset',()=>{acid='ethanoic';base='methylamine';a.value=acid;b.value=base;update();});
 function update() {
  const m=protonTransferModel(acid,base);C.state(ui,m);
  let body=T(360,65,'Proton transfer: conjugate pairs',28)+T(180,140,m.acid.HA,32)+T(540,140,m.base.B,32)+T(360,140,'+',30);
  body+=L(300,180,420,180,'#176b91',3,true)+L(420,205,300,205,'#176b91',3,true);
  body+=T(180,260,m.acid.A,32)+T(540,260,m.base.BH,32)+T(360,260,'+',30);
  body+=T(180,330,`Donor acid pKa ${m.acid.pKa}`,23)+T(540,330,`Product acid pKa ${m.base.pKa}`,23);
  const delta=Math.max(-12,Math.min(12,m.delta)),x=360+delta*22;
  body+=L(95,390,625,390)+L(360,380,360,400)+O(x,390,9,'#af3b45')+T(180,435,'Reactant-favoured',22)+T(540,435,'Product-favoured',22);
  ui.visual.innerHTML=C.svg('All species in a compatible aqueous proton-transfer comparison',body,720,465);
  ui.result.textContent=`ΔpKa = ${C.fmt(m.delta)}; K ≈ ${C.fmt(m.K)}: ${m.direction}.`;
  ui.detail.textContent=(m.samePair?'This is exchange between the same conjugate pair, so the formal equilibrium constant is one. ':'')+
   `K ≈ 10^(pKa of BH − pKa of HA), using the supplied compatible approximate aqueous values. Net reacting charge is ${m.totalCharge} on both sides. K describes equilibrium, not reaction speed or the exact final composition without starting amounts and activities.`;
 }
 update();
}
function mount20() {
 const ui=C.initDemo('day20');if(!ui)return;
 let molecule='amide',donor='N_pair',acceptor='CN_pi',paired=true;
 const mol=C.selectControl(ui.controls,'Local nitrogen environment',[['amide','Ethanamide: CH₃C(=O)NH₂'],['amine','Ethylamine: CH₃CH₂NH₂']],molecule,v=>{molecule=v;update();});
 const don=C.selectControl(ui.controls,'Arrow tail: electron source',[['N_pair','Nitrogen lone pair'],['C_O_pi','Carbonyl pi pair'],['atom_label','Atom label / positive charge']],donor,v=>{donor=v;update();});
 const acc=C.selectControl(ui.controls,'Arrow head: destination',[['CN_pi','Form a C–N pi bond'],['O_atom','Carbonyl oxygen'],['H_atom','H⁺ (proton transfer)']],acceptor,v=>{acceptor=v;update();});
 const pair=C.selectControl(ui.controls,'Coupled carbonyl pi pair → O',[['yes','Include the second pair move'],['no','Leave the C=O pi bond unchanged']],'yes',v=>{paired=v==='yes';update();});
 C.button(ui.controls,'Valid amide contributor',()=>{molecule='amide';donor='N_pair';acceptor='CN_pi';paired=true;sync();});
 C.button(ui.controls,'Compare formal amine protonation',()=>{molecule='amine';donor='N_pair';acceptor='H_atom';paired=false;sync();});
 C.button(ui.controls,'Reset',()=>{molecule='amide';donor='N_pair';acceptor='CN_pi';paired=true;sync();});
 function sync(){mol.value=molecule;don.value=donor;acc.value=acceptor;pair.value=paired?'yes':'no';update();}
 function update() {
  const m=electronModel(molecule,donor,acceptor,paired);C.state(ui,m);pair.disabled=acceptor!=='CN_pi';
  let body='';
  if(molecule==='amide') {
   for(const x of [180,540]) {
    const charged=x===540;
    body+=T(x,180,'C',28)+T(x,85,charged?'O⁻':'O',28)+T(x-105,210,'CH₃',27)+T(x+100,210,charged?'NH₂⁺':'NH₂',27);
    body+=L(x-5,155,x-5,108)+(!charged?L(x+5,155,x+5,108):'');
    body+=L(x-20,180,x-70,198)+L(x+20,180,x+65,198);
    if(charged)body+=L(x+17,190,x+62,208);else body+=T(x+100,245,'••',22);
   }
   body+=L(340,175,380,175,'#176b91',3,true)+L(380,175,340,175,'#176b91',3,true);
   body+=T(360,275,'Contributor drawings of ONE delocalised amide',23);
   if(m.kind==='amide_resonance')body+=arrowPath('M 280 230 Q 320 145 235 198')+arrowPath('M 185 130 Q 240 100 202 80');
   if(m.kind==='polar_carbonyl')body+=T(360,320,'Minor polar contributor: CH₃–C⁺(–O⁻)–NH₂',24);
  }else body+=T(360,140,'CH₃–CH₂–NH₂:',35)+T(360,220,'No adjacent carbonyl pi system',25);
  if(m.kind==='protonation')body+=T(360,310,molecule==='amide'?'Formal N protonation: CH₃C(=O)NH₃⁺':'CH₃CH₂NH₂ + H⁺ → CH₃CH₂NH₃⁺',27);
  if(!m.valid)body+=T(360,310,'Repair the proposed arrow before drawing a product.',24,'#af3b45');
  ui.visual.innerHTML=C.svg('Nitrogen electron availability and charge audit; resonance contributors are not separate species',body,720,370);
  ui.result.textContent=(m.valid?'Bookkeeping check: ':'Repair needed: ')+m.message;
  ui.detail.textContent='The two amide contributor drawings coexist as representations of one molecule: changing the selection highlights a proposed electron move, not an equilibrium switching between structures. For the N→C–N move, include the coupled C=O→O pair. Protonation adds a nucleus and changes the species, so it is a separate reaction. H⁺ abbreviates a proton donor; a full donor mechanism must also show its bond cleavage.';
 }
 update();
}
function mount21() {
 const ui=C.initDemo('day21');if(!ui)return;
 let sample='U1',candidate='ethyl_ethanoate',evidence='all';
 const sam=C.selectControl(ui.controls,'Supplied simulated unknown',Object.keys(SAMPLES).map(k=>[k,SAMPLES[k].name]),sample,v=>{sample=v;update();});
 const can=C.selectControl(ui.controls,'Test a candidate',Object.keys(CANDIDATES).map(k=>[k,CANDIDATES[k].name]),candidate,v=>{candidate=v;update();});
 const ev=C.selectControl(ui.controls,'Evidence used',[['formula','Formula only'],['IR','Formula and complete IR'],['all','Formula, IR and proton NMR']],evidence,v=>{evidence=v;update();});
 C.button(ui.controls,'Compare formula-only acid',()=>{sample='U1';candidate='butanoic_acid';evidence='formula';sync();});
 C.button(ui.controls,'Second unknown',()=>{sample='U2';candidate='ethanol';evidence='all';sync();});
 C.button(ui.controls,'Reset',()=>{sample='U1';candidate='ethyl_ethanoate';evidence='all';sync();});
 function sync(){sam.value=sample;can.value=candidate;ev.value=evidence;update();}
 function update() {
  const m=spectralCandidateModel(sample,candidate,evidence);C.state(ui,m);
  let body=T(360,60,`${m.sample.name}: ${m.sample.formula}`,30)+T(360,125,m.candidate.structure,33)+T(360,175,'Candidate formula: '+m.candidate.formula,24);
  m.tests.forEach((t,j)=>{body+=O(90,220+j*40,7,t.pass?'#23765d':'#af3b45')+T(120,228+j*40,`${t.pass?'Agrees':'Conflicts'}: ${t.label}`,24,'#153747','start');});
  ui.visual.innerHTML=C.svg('Candidate structure tested against visible simulated evidence',body,720,390)+
   (evidence!=='formula'?irPlot(sample):'')+(evidence==='all'?nmrPlot(sample):'');
  ui.result.textContent=m.fits?'This candidate fits the currently selected evidence; consider whether alternatives also fit.':'This candidate conflicts with: '+m.conflicts.join('; ')+'.';
  ui.detail.textContent=(evidence==='all'?'Supplied simulated proton NMR: '+m.sample.nmr+'. ':'')+
   'All traces are fixed teaching simulations. Interpret transmittance dips as absorption; integrated H areas are supplied separately from line heights. Use the simple first-order n+1 pattern only for the stated examples. Formula alone can leave several candidates. Agreement in this finite comparison supports a candidate; it does not establish purity or identify an unrelated real sample.';
 }
 update();
}
function mount22() {
 const ui=C.initDemo('day22');if(!ui)return;
 let front=5,spot=3.2,recorded=true,composition='single';
 const fr=C.rangeControl(ui.controls,'Solvent front from origin',1,10,.1,front,'cm',v=>{front=v;sp.max=front;if(spot>front){spot=front;resetRange(sp,spot);}update();});
 const sp=C.rangeControl(ui.controls,'Spot centre from origin',0,front,.05,spot,'cm',v=>{spot=v;update();});
 const rec=C.selectControl(ui.controls,'Was the front recorded?',[['yes','Yes: use the measured front'],['no','No: denominator unavailable']],'yes',v=>{recorded=v==='yes';update();});
 const mix=C.selectControl(ui.controls,'Known teaching composition',[['single','One supplied component'],['co_migrate','Two components with the same retention']],composition,v=>{composition=v;update();});
 C.button(ui.controls,'Spot halfway to front',()=>resetRange(sp,front/2));
 C.button(ui.controls,'Reset',()=>{front=5;spot=3.2;recorded=true;composition='single';rec.value='yes';mix.value='single';resetRange(fr,5);resetRange(sp,3.2);});
 function update() {
  const m=tlcModel(front,spot,recorded,composition);C.state(ui,m);
  const y=v=>320-v*23;let body=`<rect x="100" y="65" width="380" height="280" fill="#fff" stroke="#176b91" stroke-width="3"/>`;
  body+=L(120,320,460,320,'#153747',2,false,'5 4')+T(560,330,'Origin',22);
  if(recorded)body+=L(120,y(front),460,y(front),'#23765d',2,false,'6 4')+T(565,y(front)+7,front.toFixed(1)+' cm front',22);
  else body+=T(290,100,'Front unrecorded',25,'#af3b45');
  body+=O(290,y(spot),12,'#176b91')+L(500,320,500,y(spot),'#176b91',2,true)+T(600,y(spot)+7,spot.toFixed(2)+' cm spot',22);
  body+=T(360,390,m.componentCount===2?'Two known components can produce this one visible feature.':'One visible feature does not by itself prove a single component.',22);
  ui.visual.innerHTML=C.svg('TLC geometry with a common origin and optionally missing solvent front',body,720,430);
  ui.result.textContent=m.Rf===null?'Rf unavailable without a recorded valid front.':`Rf = ${spot.toFixed(2)} / ${front.toFixed(1)} = ${C.fmt(m.Rf)} (dimensionless).`;
  ui.detail.textContent=m.message+' The sliders are constrained to 0 ≤ spot ≤ front. Both distances use the same origin and units. The composition selector supplies a counterexample: co-migrating compounds can share one feature. Rf also depends on eluent, stationary phase and conditions; spot count alone does not establish purity or unique identity.';
 }
 update();
}
function mount23() {
 const ui=C.initDemo('day23');if(!ui)return;
 let nodes=['ethanol'],target='ester',action='oxidize_acid',last=null;
 const tar=C.selectControl(ui.controls,'Target functional group',[['ester','Ethyl ethanoate'],['acid','Ethanoic acid']],target,v=>{target=v;update();});
 const act=C.selectControl(ui.controls,'Proposed next transformation',Object.keys(ACTIONS).map(k=>[k,ACTIONS[k]]),action,v=>{action=v;update();});
 C.button(ui.controls,'Apply the proposed step',()=>{last=routeStep(nodes[nodes.length-1],action);if(last.valid)nodes.push(last.next);update();});
 C.button(ui.controls,'Undo last accepted step',()=>{if(nodes.length>1)nodes.pop();last=null;update();});
 C.button(ui.controls,'Reset',()=>{nodes=['ethanol'];target='ester';action='oxidize_acid';last=null;tar.value=target;act.value=action;update();});
 function update() {
  const m=routeModel(nodes,target,last);C.state(ui,{...m,action});
  const current=NODES[m.current];let body=T(360,45,'Current: '+current.name+' · '+current.group,26);
  if(last&&last.valid) {
   body+=T(360,110,last.equation[0],31)+L(360,140,360,185,'#176b91',3,true)+T(360,235,last.equation[1],31);
   if(last.action==='esterify')body+=T(360,275,'Separate ethanol portion supplied; water is a product.',22);
  } else body+=T(360,135,current.formula,38)+T(360,220,'Target: '+NODES[target].formula,31);
  const path=nodes.map(k=>NODES[k].name).join(' → ');
  if(path.length<=65)body+=T(360,340,path,19);
  else body+=T(360,330,`${nodes.length-1} accepted steps; use Undo to simplify.`,22);
  ui.visual.innerHTML=C.svg('Functional-group route with explicit extra ethanol and water in esterification',body,720,385);
  ui.result.textContent=m.achieved?'The current structure matches the selected target. Explain the accepted bond changes and needed evidence.':last?last.reason:'Choose a transformation that matches the current functional group, then apply it.';
  ui.detail.textContent=(last?last.evidence+' ':'')+'This is a short route-reasoning model, not a laboratory protocol or a guarantee of selectivity. [O] is a net oxidation equivalent, not a specified reagent. H:⁻ abbreviates delivery from a suitable reducing reagent, not free hydride in water; a real mechanism includes reagent chemistry. Oxidation and work-up conditions are specified only at the functional-group level. Esterification uses a separate ethanol portion and is an equilibrium. Base hydrolysis ends at carboxylate until an explicit acid work-up is selected.';
 }
 update();
}
function mount24() {
 const ui=C.initDemo('day24');if(!ui)return;
 let IR='carbonyl',NMR='relative',claim='OH_absent';
 const ir=C.selectControl(ui.controls,'Available IR coverage',[['full','Complete trace; adequate O–H sensitivity'],['carbonyl','Carbonyl window only'],['none','IR unavailable']],IR,v=>{IR=v;update();});
 const nm=C.selectControl(ui.controls,'Available proton-NMR information',[['full','Shifts, multiplicities and reliable relative integrals'],['relative','Relative integrated areas only (3:3:2)'],['multiplicity','Splitting groups only; no areas or shifts'],['none','NMR unavailable']],NMR,v=>{NMR=v;update();});
 const cl=C.selectControl(ui.controls,'Claim to audit',[
  ['OH_absent','The O–H absence is observed evidence'],['ester_best','The ester best fits this candidate comparison'],
  ['acid_challenged','The acid candidate is challenged'],['integrals_eight','The integrals account for eight H'],
  ['OCH2_assigned','An O–CH₂–CH₃ group is assigned'],['purity','The sample is proven pure'],['unique','This proves the only possible molecular structure']
 ],claim,v=>{claim=v;update();});
 C.button(ui.controls,'Supply the full evidence',()=>{IR='full';NMR='full';ir.value=IR;nm.value=NMR;update();});
 C.button(ui.controls,'Reset missing-region challenge',()=>{IR='carbonyl';NMR='relative';claim='OH_absent';ir.value=IR;nm.value=NMR;cl.value=claim;update();});
 function update() {
  const m=evidenceStrengthModel(IR,NMR,claim);C.state(ui,m);
  let body=T(360,60,'Fixed simulated sample: C₄H₈O₂',29)+T(360,120,'Candidate set: ester, butanoic acid, hydroxyketone',23)+T(360,185,m.supported?'Claim supported within the stated comparison':'Claim not supported by this evidence',25,m.supported?'#23765d':'#af3b45');
  ui.visual.innerHTML=C.svg('Audit the strength of a claim when evidence is incomplete',body,720,235)+
   (IR==='none'?C.svg('Infrared evidence unavailable',T(360,110,'IR not supplied: no absence inference is possible',25),720,180):irPlot('U1',IR))+
   (NMR==='none'?C.svg('Proton NMR evidence unavailable',T(360,110,'NMR not supplied',26),720,180):nmrPlot('U1',{mode:NMR==='relative'?'relative':NMR==='multiplicity'?'multiplicity':'full'}));
  ui.result.textContent=m.reason;
  ui.detail.textContent='Changing coverage hides observations rather than changing the molecule. Missing O–H coverage must not look like a flat measured baseline. Relative areas can be normalised using the formula but do not supply chemical shifts. Only the selected observations should be used to justify the current claim; changing coverage does not invent a different molecular sample. None of these qualitative candidate tests proves purity or every possible alternative false.';
 }
 update();
}
function claimSpectrum(pattern,source) {
 let b=T(360,35,'Hypothetical target-reference pattern',25)+L(100,155,650,155);
 const peaks=[[190,60],[330,95],[550,70]];
 peaks.forEach(([x,h])=>{b+=L(x,155,x,155-h,'#176b91',3);});
 b+=T(360,215,source==='same'?'Pattern explicitly linked to this isolate':'Unrelated teaching pattern',24,source==='same'?'#153747':'#af3b45')+L(100,340,650,340);
 peaks.forEach(([x,h])=>{b+=L(x,340,x,340-h,'#176b91',3);});
 if(pattern==='contradictory')b+=L(430,340,430,250,'#af3b45',4)+T(500,245,'Extra group',21,'#af3b45');
 b+=T(360,390,'Schematic comparison; shifts not assigned',20);
 return C.svg('Hypothetical reference versus observation; provenance and extra signals change the claim',b,720,420);
}
function mount25() {
 const ui=C.initDemo('day25');if(!ui)return;
 let mass=1.2,source='same',pattern='consistent',spots=1,conversionData='none',claim='yield';
 const ma=C.rangeControl(ui.controls,'Isolated material mass',.8,2.2,.01,mass,'g',v=>{mass=v;update();});
 const so=C.selectControl(ui.controls,'Structural spectrum provenance',[['same','Explicitly linked to this isolate'],['unrelated','Unrelated simulated teaching sample'],['none','No structural spectrum supplied']],source,v=>{source=v;update();});
 const pa=C.selectControl(ui.controls,'Supplied spectral pattern',[['consistent','Consistent target pattern'],['contradictory','Additional major unexpected group']],pattern,v=>{pattern=v;update();});
 const sp=C.selectControl(ui.controls,'Detected TLC features under one condition',[['1','One visible spot'],['2','Two resolved visible spots']],'1',v=>{spots=Number(v);update();});
 const co=C.selectControl(ui.controls,'Quantitative starting-material assay',[['none','Not measured'],['residual2','Same batch: 2.00 mmol remains from 10.00 mmol']],conversionData,v=>{conversionData=v;update();});
 const cl=C.selectControl(ui.controls,'Claim to evaluate',[['yield','Calculate apparent isolated yield'],['conversion','State starting-material conversion'],['identity','Support target identity for this isolate'],['purity','Claim the isolate is proven pure'],['success','Claim pure product and known conversion from yield alone']],claim,v=>{claim=v;update();});
 C.button(ui.controls,'Unrelated-spectrum trap',()=>{source='unrelated';claim='identity';so.value=source;cl.value=claim;update();});
 C.button(ui.controls,'High-mass contradiction',()=>{conversionData='residual2';co.value=conversionData;resetRange(ma,1.9);});
 C.button(ui.controls,'Reset',()=>{mass=1.2;source='same';pattern='consistent';spots=1;conversionData='none';claim='yield';so.value=source;pa.value=pattern;sp.value='1';co.value=conversionData;cl.value=claim;resetRange(ma,1.2);});
 function update() {
  const m=integratedModel(mass,source,pattern,spots,conversionData,claim);C.state(ui,m);pa.disabled=source==='none';
  let body=T(360,45,'Fixed theoretical product: 10.00 mmol × 180.16 g/mol',24);
  body+=`<rect x="90" y="90" width="540" height="38" fill="#dde7eb"/><rect x="90" y="150" width="${Math.min(600,540*mass/m.theoretical)}" height="38" fill="#176b91"/>`;
  body+=T(360,118,'Theoretical 1.8016 g',24)+T(360,180,'Isolated '+mass.toFixed(2)+' g',24,'#fff');
  body+=T(360,235,`Apparent yield ${C.fmt(m.apparentYield,4)}%`,30)+T(360,290,m.conversion===null?'Conversion: not measured':'Conversion from residual assay: 80.0%',26);
  body+=T(360,345,source==='unrelated'?'Spectrum belongs to a different teaching sample':source==='none'?'No sample-linked structural spectrum':'Structural evidence is explicitly linked to this isolate',22,source==='unrelated'?'#af3b45':'#153747');
  ui.visual.innerHTML=C.svg('Independent mass, conversion and structural evidence streams',body,720,380)+
   (source!=='none'?claimSpectrum(pattern,source):'');
  ui.result.textContent=(m.supported?'Supported conditional claim: ':'Insufficient or conflicting evidence: ')+m.reason;
  ui.detail.textContent=(m.inconsistent?'Audit needed: apparent product-mass yield exceeds the supplied 80% conversion. Check purity, losses, quantitative assumptions and whether records belong to the same batch. ':'')+
   (m.apparentYield>100?'An apparent mass percentage above 100% conflicts with the pure-product maximum; investigate retained solvent, impurities or input errors. ':'')+
   'All observations here are fixed hypothetical teaching scenarios. The unassigned pattern comparison is schematic, not an aspirin spectrum or the U1 spectrum from another demo; the arithmetic product is intentionally unnamed. “Same isolate” stipulates a pattern comparison explicitly linked to that hypothetical sample. A real report needs actual linked spectra and appropriate composition checks. One TLC spot does not prove purity. The optional conversion assay is stipulated quantitative and same-batch; isolated mass alone cannot supply it.';
 }
 update();
}
[mount18,mount19,mount20,mount21,mount22,mount23,mount24,mount25].forEach(mount=>mount());
})(typeof globalThis!=='undefined'?globalThis:window);
