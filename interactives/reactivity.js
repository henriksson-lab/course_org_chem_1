(function (global) {
"use strict";

// These deliberately small models express the stated teaching assumptions.
// They do not select real laboratory conditions or infer a mechanism from a name.
const finiteNonnegative = value => Number.isFinite(value) && value >= 0;
const wrapAngle = angle => ((angle % 360) + 360) % 360;

function enantiomerModel(rPercent) {
  if (!Number.isFinite(rPercent) || rPercent < 0 || rPercent > 100)
    return {valid:false, error:"R composition must be between 0 and 100%."};
  const sPercent=100-rPercent;
  return {valid:true,rPercent,sPercent,ee:Math.abs(rPercent-sPercent),
    major:rPercent===sPercent?"neither (racemic)":rPercent>sPercent?"R":"S"};
}

function rateModel(substrate,nucleophile) {
  if (![substrate,nucleophile].every(finiteNonnegative))
    return {valid:false,error:"Concentrations must be finite and nonnegative."};
  const k1=0.010, k2=0.100;
  if(!Number.isFinite(k1*substrate)||!Number.isFinite(k2*substrate*nucleophile))
    return {valid:false,error:"These inputs overflow the finite teaching calculation."};
  return {valid:true,substrate,nucleophile,k1,k2,sn1:k1*substrate,
    sn2:k2*substrate*nucleophile, equalAtNu:k1/k2};
}

function e2Model(rotation,hydrogen) {
  if (!Number.isFinite(rotation) || !Number.isInteger(hydrogen) || hydrogen<0 || hydrogen>2)
    return {valid:false,error:"Choose one of the three front-carbon hydrogens and a finite rotation."};
  const frontAngle=hydrogen*120,brAngle=wrapAngle(rotation);
  const dihedral=wrapAngle(brAngle-frontAngle);
  const deviation=Math.abs(dihedral-180);
  return {valid:true,rotation:brAngle,hydrogen,frontAngle,dihedral,deviation,
    aligned:deviation<=10+1e-9,exactAnti:deviation<1e-9};
}

const oxidationCases={
  ethanol:{name:"Ethanol",left:"CH₃",right:"H",top:"OH",bottom:"H",co:1,ch:2,next:"ethanal"},
  ethanal:{name:"Ethanal",left:"CH₃",right:"H",top:"O",bottom:null,co:2,ch:1,next:"acid"},
  acid:{name:"Ethanoic acid",left:"CH₃",right:"OH",top:"O",bottom:null,co:3,ch:0,next:null},
  secondary:{name:"Propan-2-ol",left:"CH₃",right:"CH₃",top:"OH",bottom:"H",co:1,ch:1,next:"ketone"},
  ketone:{name:"Propanone",left:"CH₃",right:"CH₃",top:"O",bottom:null,co:2,ch:0,next:null},
  tertiary:{name:"2-Methylpropan-2-ol",left:"CH₃",right:"CH₃",top:"OH",bottom:"CH₃",co:1,ch:0,next:null}
};
function oxidationModel(key) {
  const item=oxidationCases[key];
  return item?{valid:true,key,...item,value:item.co-item.ch,
    boundary:key==="tertiary"?"No H is attached to the OH-bearing carbon. A simple carbon-skeleton-retaining C=O product would give carbon five bond orders.":"Use only the selected transformations supplied here."}
    :{valid:false,error:"Unknown teaching structure."};
}

const sn1States=[
  {name:"Tertiary bromide before ionisation",group:"Br",centralCharge:0,oxygenCharge:0,waters:2,br:false,hydronium:false},
  {name:"Carbocation after ionisation",group:null,centralCharge:1,oxygenCharge:0,waters:2,br:true,hydronium:false},
  {name:"Oxonium after water capture",group:"OH₂⁺",centralCharge:0,oxygenCharge:1,waters:1,br:true,hydronium:false},
  {name:"Neutral alcohol after proton transfer",group:"OH",centralCharge:0,oxygenCharge:0,waters:0,br:true,hydronium:true}
];
const sn1Correct=["ionise","waterAttack","deprotonate"];
const sn1Reasons={
  ionise:"Move the C–Br bond pair onto Br. Carbon becomes +1; Br becomes −1. Both waters remain available.",
  waterAttack:"A water oxygen lone pair forms C–O. The incoming oxygen retains two H atoms, so it has three bonds and charge +1.",
  deprotonate:"A second water accepts H from oxonium; the original O–H pair returns to that oxygen. Alcohol, H₃O⁺ and Br⁻ result.",
  eraseH:"Deleting H without a recipient loses an atom and changes charge. Draw a real proton acceptor and the returning O–H pair.",
  neutralCapture:"Water cannot gain a third bond and remain the ordinary neutral alcohol while retaining its two H atoms. Form oxonium first.",
  methylCation:"The specified substrate is tertiary. Replacing its three methyl substituents with a methyl cation would change the atoms and mechanism."
};
function sn1Model(stage,action) {
  if (!Number.isInteger(stage) || stage<0 || stage>3) return {valid:false,error:"Unknown mechanism state."};
  const expected=sn1Correct[stage],accepted=stage<3 && action===expected;
  const item=sn1States[stage];
  let rejected=sn1Reasons[action]||"Choose the electron-pair operation appropriate to the current species.";
  if(action==="ionise"&&stage>0)rejected="The C–Br bond has already broken: Br⁻ is now a spectator. Find the electron-pair operation appropriate to the current carbon/oxygen state.";
  if(action==="waterAttack"&&stage!==1)rejected=stage===0?"There is no carbocation yet in this selected SN1 route. Ionise the tertiary C–Br bond first.":"Carbon already has its C–O bond and four bond orders. Adding another water directly would overbond it.";
  if(action==="deprotonate"&&stage!==2)rejected="The oxonium species that loses this proton has not formed. Identify the current species before using the oxonium deprotonation step.";
  return {valid:true,stage,...item,totalCharge:item.centralCharge+item.oxygenCharge+(item.hydronium?1:0)-(item.br?1:0),expected,
    accepted,nextStage:accepted?stage+1:stage,
    feedback:stage===3?"The selected substitution branch is complete. Bromide remains a spectator; hydronium carries the balancing positive charge.":
      accepted?sn1Reasons[action]:rejected};
}

function yieldModel(acidMass,anhydrideMass,isolatedMass) {
  if (![acidMass,anhydrideMass,isolatedMass].every(finiteNonnegative))
    return {valid:false,error:"All masses must be present, finite and nonnegative. A blank field is not a zero measurement."};
  const acidMol=acidMass/138.12,anhydrideMol=anhydrideMass/102.09;
  const extent=Math.min(acidMol,anhydrideMol),theoreticalMass=extent*180.16;
  const limiting=Math.abs(acidMol-anhydrideMol)<1e-12?"both (equal stoichiometric amounts)":acidMol<anhydrideMol?"salicylic acid":"acetic anhydride";
  const yieldPercent=theoreticalMass>0?100*isolatedMass/theoreticalMass:null;
  if(!Number.isFinite(theoreticalMass)||(yieldPercent!==null&&!Number.isFinite(yieldPercent)))
    return {valid:false,error:"These masses are too extreme for a finite yield calculation. Use laboratory-scale teaching amounts."};
  return {valid:true,acidMol,anhydrideMol,extent,theoreticalMass,isolatedMass,limiting,
    yieldPercent,aboveMaximum:yieldPercent!==null&&yieldPercent>100};
}

const additionStates={
  hydrate:[
    {name:"Before hydroxide addition",form:"carbonyl",a:0,b:-1,c:0,partner:"none",side:"OHᴮ⁻ + H₂Oᶜ",explain:"A is carbonyl oxygen; B is incoming hydroxide oxygen; C belongs to the water proton donor."},
    {name:"After hydroxide addition",form:"tetra",a:-1,b:0,c:0,partner:"none",side:"H₂Oᶜ",explain:"Oᴮ forms the new C–O bond and remains OH. The C=Oᴬ π pair returns to A, giving Oᴬ⁻."},
    {name:"After water proton donation",form:"tetra",a:0,b:0,c:-1,partner:"waterDonates",side:"OHᶜ⁻",explain:"Water C donates H to A and becomes hydroxide C. Hydroxide species is regenerated, but its labelled oxygen need not be B."}
  ],
  hemiacetal:[
    {name:"Before acid activation",form:"carbonyl",a:0,b:0,c:1,partner:"none",side:"CH₃OᴮH + H₃Oᶜ⁺",explain:"A is carbonyl O; B belongs to methanol; C belongs to the hydronium/water proton-transfer partner."},
    {name:"After carbonyl O protonation",form:"activated",a:1,b:0,c:0,partner:"hydroniumDonates",side:"CH₃OᴮH + H₂Oᶜ",explain:"Hydronium C supplies H to A and becomes water. Protonated carbonyl Oᴬ has formal charge +1."},
    {name:"After methanol attack",form:"tetra",a:0,b:1,c:0,partner:"none",side:"H₂Oᶜ",explain:"Methanol B attacks carbon; the carbonyl π pair returns to A. A is now neutral OH, while B has three bonds and charge +1."},
    {name:"After water accepts the extra proton",form:"tetra",a:0,b:0,c:1,partner:"waterAccepts",side:"H₃Oᶜ⁺",explain:"Water C accepts H from B; the Oᴮ–H pair returns to B. Neutral hemiacetal and hydronium are formed."}
  ]
};
function additionModel(route,stage) {
  const items=additionStates[route];
  if (!items || !Number.isInteger(stage) || !items[stage]) return {valid:false,error:"Unknown route or step."};
  const s=items[stage];
  return {valid:true,route,stage,...s,totalCharge:s.a+s.b+s.c,last:items.length-1};
}

const acylStates=[
  {name:"Ester plus hydroxide",form:"ester",totalCharge:-1,side:"OH⁻",explain:"Methyl ethanoate has an electrophilic acyl carbon. Begin by addition, not direct displacement at a saturated carbon."},
  {name:"Tetrahedral addition intermediate",form:"tetra",totalCharge:-1,side:"",explain:"The original carbonyl oxygen is O⁻. Carbon also carries OH and OCH₃. Collapse can restore C=O; reversal to reactants is also possible."},
  {name:"Transient acid–methoxide pair",form:"acid",totalCharge:-1,side:"CH₃O⁻",explain:"Collapse has expelled methoxide and formed an acid. In this basic model the acid–base pair does not represent the final dominant mixture."},
  {name:"Carboxylate plus methanol: proton-transfer trap",form:"carboxylate",totalCharge:-1,side:"CH₃OH",explain:"Methoxide removes the acid OH proton. Carboxylate and methanol result. Acid work-up is not included in this model."}
];
const acylReasons={
  add:"A hydroxide pair forms C–O while the C=O π pair returns to the original oxygen, giving the tetrahedral intermediate.",
  collapse:"The O⁻ pair reforms C=O as the acyl C–OCH₃ bond pair goes to methoxy oxygen. Acid plus methoxide is the bookkeeping product pair.",
  reverse:"Reforming C=O while expelling OH⁻ returns to the ester reactants. This is a plausible reverse addition; product trapping has not happened yet.",
  trap:"Methoxide accepts the acid proton. The acid O–H pair stays with carboxylate oxygen, giving carboxylate plus methanol.",
  direct:"For the selected basic ester-hydrolysis route, attack first produces a tetrahedral intermediate. Do not replace carbonyl addition with an unsupported direct saturated-carbon SN2 drawing.",
  hydride:"No H is bonded to this ester acyl carbon, and hydride is not an appropriate ordinary leaving group here. Check actual groups and leaving ability.",
  methyl:"Expelling CH₃⁻ would break the retained C–C skeleton and generate a very poor leaving group. The oxygen-based group is the appropriate leaving option in the selected collapse.",
  stop:"Neutral acid plus methoxide has a strongly favoured proton transfer. Stopping here misses the carboxylate product under the stated basic conditions."
};
function acylModel(stage,action) {
  if (!Number.isInteger(stage)||stage<0||stage>3)return {valid:false,error:"Unknown acyl state."};
  let accepted=false,nextStage=stage;
  if(stage===0&&action==="add"){accepted=true;nextStage=1;}
  if(stage===1&&action==="collapse"){accepted=true;nextStage=2;}
  if(stage===1&&action==="reverse"){accepted=true;nextStage=0;}
  if(stage===2&&action==="trap"){accepted=true;nextStage=3;}
  let feedback=acylReasons[action]||"Choose an electron-pair operation appropriate to the current structure.";
  if(!accepted){
    if(stage===0&&["collapse","reverse"].includes(action))feedback="There is no tetrahedral O⁻ intermediate to collapse yet. Begin by nucleophilic addition to the carbonyl.";
    if(stage===1&&action==="add")feedback="This tetrahedral carbon already has four single bonds. A further direct addition would overbond carbon. Consider collapse or reverse addition.";
    if(stage!==2&&action==="trap")feedback="The neutral acid–methoxide pair has not formed in this state. Do not transfer a proton between species that are not yet present.";
    if(stage===2&&["collapse","reverse","add"].includes(action))feedback="The carbonyl has already been restored. The strongly favoured next operation in this selected basic model is proton transfer from the acid to methoxide.";
  }
  return {valid:true,stage,...acylStates[stage],accepted,nextStage,
    trapK:10**(15.5-4.8),feedback:stage===3?acylStates[3].explain:feedback};
}

function checks() {
  const close=(a,b,tol=1e-8)=>Math.abs(a-b)<tol;
  return [
    {name:"70:30 gives 40% R ee",pass:enantiomerModel(70).ee===40&&enantiomerModel(70).major==="R"},
    {name:"racemate has no major enantiomer",pass:enantiomerModel(50).ee===0},
    {name:"invalid enantiomer percentage rejected",pass:!enantiomerModel(101).valid},
    {name:"supplied rate laws cross at 0.1 M Nu",pass:close(rateModel(.05,.1).sn1,rateModel(.05,.1).sn2)},
    {name:"SN1 model independent of Nu",pass:rateModel(.1,0).sn1===rateModel(.1,.2).sn1},
    {name:"SN2 model zero without Nu",pass:rateModel(.1,0).sn2===0},
    {name:"negative concentration rejected",pass:!rateModel(-1,.1).valid},
    {name:"first hydrogen anti at180deg",pass:e2Model(180,0).exactAnti},
    {name:"second hydrogen anti at300deg",pass:e2Model(300,1).exactAnti},
    {name:"third hydrogen anti at60deg",pass:e2Model(60,2).exactAnti},
    {name:"eclipsed is not anti",pass:!e2Model(0,0).aligned},
    {name:"ethanol changing carbon is minus one",pass:oxidationModel("ethanol").value===-1},
    {name:"acid carbon is plus three",pass:oxidationModel("acid").value===3},
    {name:"tertiary centre has no simple oxidation next",pass:oxidationModel("tertiary").next===null},
    {name:"SN1 ionisation accepted first",pass:sn1Model(0,"ionise").nextStage===1},
    {name:"neutral water capture shortcut rejected",pass:!sn1Model(1,"neutralCapture").accepted},
    {name:"all SN1 states include zero total charge",pass:sn1States.every((s,i)=>sn1Model(i).totalCharge===0)},
    {name:"default theoretical mass is1.8016g",pass:close(yieldModel(1.3812,1.53135,1.2).theoreticalMass,1.8016)},
    {name:"default yield is about66.607%",pass:close(yieldModel(1.3812,1.53135,1.2).yieldPercent,66.6074600355,1e-7)},
    {name:"zero theoretical amount has undefined yield",pass:yieldModel(0,1,0).yieldPercent===null},
    {name:"negative and blank masses rejected",pass:!yieldModel(-1,1,1).valid&&!yieldModel(NaN,1,1).valid},
    {name:"above maximum is flagged rather than clamped",pass:yieldModel(1.3812,1.53135,1.95).aboveMaximum},
    {name:"hydrate original oxygen initially becomes minus",pass:additionModel("hydrate",1).a===-1},
    {name:"hydrate preserves net minus one",pass:additionStates.hydrate.every((s,i)=>additionModel("hydrate",i).totalCharge===-1)},
    {name:"hemiacetal incoming O initially gains plus",pass:additionModel("hemiacetal",2).b===1},
    {name:"acid route preserves net plus one",pass:additionStates.hemiacetal.every((s,i)=>additionModel("hemiacetal",i).totalCharge===1)},
    {name:"ester addition precedes collapse",pass:acylModel(0,"add").nextStage===1&&!acylModel(0,"collapse").accepted},
    {name:"tetrahedral reversal is a valid branch",pass:acylModel(1,"reverse").nextStage===0&&acylModel(1,"reverse").accepted},
    {name:"methoxide acid-base trap favoured",pass:acylModel(2,"trap").nextStage===3&&acylModel(2).trapK>1e10}
  ];
}

const api={enantiomerModel,rateModel,e2Model,oxidationModel,sn1Model,yieldModel,additionModel,acylModel,checks};
if(typeof module!=="undefined"&&module.exports)module.exports=api;
if(typeof document==="undefined"||!global.ChemDemo)return;
const {initDemo,selectControl,rangeControl,numberControl,button,state,svg,text,line,circle,fmt,plot}=global.ChemDemo;
const numerical=input=>input.value.trim()===""?NaN:Number(input.value);
const change=(input,value)=>{input.value=value;input.dispatchEvent(new Event("input",{bubbles:true}));};
const say=(ui,message,detail)=>{ui.result.textContent=message;ui.detail.textContent=detail;};
const note=(label,body)=>text(360,325,label,23)+text(360,355,body,21);

// A drawn carbonyl/tetrahedral site with retained carbon and H/groups visible.
function carbonSite({form="carbonyl",top="O",right="H",bottom=null,left="CH₃",cx=230,cy=175}) {
  let b=text(cx,cy+9,"C",31)+text(cx-125,cy+9,left,27)+line(cx-25,cy,cx-91,cy);
  b+=text(cx+130,cy+9,right,27)+line(cx+23,cy,cx+93,cy);
  b+=text(cx,cy-90,top,28,"#af3b45");
  if(form==="carbonyl"||form==="activated")b+=line(cx-5,cy-22,cx-5,cy-65)+line(cx+5,cy-22,cx+5,cy-65);
  else b+=line(cx,cy-22,cx,cy-65);
  if(bottom)b+=line(cx,cy+20,cx,cy+73)+text(cx,cy+108,bottom,27,"#176b91");
  return b;
}

function mount10() {
  const ui=initDemo("day10");if(!ui)return;
  const r=rangeControl(ui.controls,"R among R+S",0,100,1,70,"%",update);
  button(ui.controls,"Racemate",()=>change(r,50));button(ui.controls,"95:5 R:S",()=>change(r,95));button(ui.controls,"Reset",()=>change(r,70));
  function update(){const m=enantiomerModel(Number(r.value));
    const w=600*m.rPercent/100;
    ui.visual.innerHTML=svg("R and S fractions and their excess",
      `<rect x="60" y="100" width="${w}" height="95" fill="#176b91"/><rect x="${60+w}" y="100" width="${600-w}" height="95" fill="#23765d"/>`+
      text(360,65,`${fmt(m.rPercent)}% R + ${fmt(m.sPercent)}% S`,31)+
      text(360,245,`${fmt(m.ee)}% enantiomeric excess`,31)+text(360,300,"Composition does not determine (+)/(−) rotation sign.",23));
    say(ui,`${fmt(m.ee)}% ee; major enantiomer: ${m.major}.`,"Only R and S are counted and their percentages sum to 100. R/S describes configuration, not optical-rotation sign or biological potency. Other chemical impurities are outside this composition model.");
    state(ui,{...m,rotationSign:"not inferred"});}
  update();
}

function mount11() {
  const ui=initDemo("day11");if(!ui)return;
  const rx=rangeControl(ui.controls,"Assumed [RX]",0,.2,.01,.05,"mol/L",update);
  const nu=rangeControl(ui.controls,"Assumed [Nu]",0,.2,.01,.1,"mol/L",update);
  button(ui.controls,"No nucleophile",()=>change(nu,0));button(ui.controls,"Double both from default",()=>{change(rx,.1);change(nu,.2);});
  button(ui.controls,"Reset",()=>{change(rx,.05);change(nu,.1);});
  function update(){const m=rateModel(Number(rx.value),Number(nu.value));
    const xs=Array.from({length:41},(_,i)=>i*.005),ymax=Math.max(.0002,.1*m.substrate*.2)*1.15;
    ui.visual.innerHTML=plot({title:"Two assumed rate laws at fixed substrate concentration",xDomain:[0,.2],yDomain:[0,ymax],xLabel:"[Nu] / mol L⁻¹",yLabel:"Initial rate / mol L⁻¹ s⁻¹",series:[
      {label:"Assumed SN1",color:"#176b91",points:xs.map(x=>[x,m.k1*m.substrate])},{label:"Assumed SN2",color:"#23765d",points:xs.map(x=>[x,m.k2*m.substrate*x])}],
      markers:[{x:m.nucleophile,y:m.sn1,color:"#176b91"},{x:m.nucleophile,y:m.sn2,color:"#23765d"}]});
    say(ui,`Blue assumed SN1: ${fmt(m.sn1)} mol L⁻¹ s⁻¹. Green assumed SN2: ${fmt(m.sn2)} mol L⁻¹ s⁻¹.`,
      "Illustrative supplied constants: k₁ = 0.010 s⁻¹; k₂ = 0.100 L mol⁻¹ s⁻¹. Models: v₁=k₁[RX] and v₂=k₂[RX][Nu]. They are separate assumed pathways, not a mechanism selector or a measured competition for a specified real substrate. No nucleophile gives zero SN2 model rate; zero substrate gives zero in both.");
    state(ui,m);}
  update();
}

function mount12() {
  const ui=initDemo("day12");if(!ui)return;
  const h=selectControl(ui.controls,"Reacting front-carbon hydrogen",[["0","H1 (12 o’clock)"],["1","H2 (4 o’clock)"],["2","H3 (8 o’clock)"]],"0",update);
  const rotation=rangeControl(ui.controls,"Rear Br rotation from 12 o’clock",0,360,10,60,"degrees",update);
  button(ui.controls,"Align selected H anti",()=>change(rotation,wrapAngle(Number(h.value)*120+180)));
  button(ui.controls,"Reset",()=>{h.value="0";change(rotation,60);});
  function update(){const m=e2Model(Number(rotation.value),Number(h.value));const cx=250,cy=185;
    let b=circle(cx,cy,44,"none","#153747")+circle(cx,cy,6,"#153747");
    const pos=(a,r)=>[cx+r*Math.sin(a*Math.PI/180),cy-r*Math.cos(a*Math.PI/180)];
    for(let j=0;j<3;j++){const a=m.rotation+j*120,p=pos(a,44),q=pos(a,126),t=pos(a,158);
      b+=line(...p,...q,j===0?"#176b91":"#153747")+text(t[0],t[1]+7,j===0?"Br":"H",26);}
    for(let j=0;j<3;j++){const q=pos(j*120,85),t=pos(j*120,110),col=j===m.hydrogen?"#af3b45":"#153747";
      b+=line(cx,cy,...q,col,4)+text(t[0],t[1]+7,`H${j+1}`,25,col);}
    b+=text(550,85,`${fmt(m.dihedral)}° dihedral`,28)+text(550,145,m.aligned?"Near anti":"Not near anti",29,m.aligned?"#23765d":"#af3b45")+
      text(550,215,"Front: β-CH₃",24)+text(550,255,"Rear: α-CH₂Br",24)+text(360,415,"CH₃CH₂Br → CH₂=CH₂",25);
    ui.visual.innerHTML=svg("Selectable bromoethane Newman geometry",b,720,450);
    say(ui,`${m.exactAnti?"Exactly anti":m.aligned?"Within the illustrative anti window":"Outside the illustrative anti window"}; H${m.hydrogen+1}–Cβ–Cα–Br = ${fmt(m.dihedral)}°.`,
      "Geometry-only model: the front carbon and its three labelled H atoms stay fixed; the rear CH₂Br rotates. A ±10° teaching window flags near 180°, not a reaction probability. Bromoethane can favour substitution with suitable nucleophiles; this explorer does not predict a dominant E2 pathway or rate.");
    state(ui,m);}
  update();
}

function mount13() {
  const ui=initDemo("day13");if(!ui)return;let checked=false;
  const compound=selectControl(ui.controls,"Structure to inspect",Object.keys(oxidationCases).map(k=>[k,oxidationCases[k].name]),"ethanol",()=>{checked=false;update();});
  const guess=selectControl(ui.controls,"My oxidation-state assignment",Array.from({length:9},(_,i)=>[String(i-4),String(i-4)]),"0",()=>{checked=false;update();});
  button(ui.controls,"Check my count",()=>{checked=true;update();});
  const advance=button(ui.controls,"Selected next oxidation",()=>{const m=oxidationModel(compound.value);if(m.next){compound.value=m.next;checked=false;}update();});
  button(ui.controls,"Reset",()=>{compound.value="ethanol";guess.value="0";checked=false;update();});
  function update(){const m=oxidationModel(compound.value);
    ui.visual.innerHTML=svg(`${m.name}: changing carbon bond inventory`,carbonSite({form:m.top==="O"?"carbonyl":"tetra",...m})+
      text(565,95,`${m.co} C–O bond orders`,25)+text(565,155,`${m.ch} C–H bonds`,25)+text(565,230,`${m.co} − ${m.ch} = ${m.value}`,32)+
      text(360,345,"C–O: +1 per bond order; C–H: −1; C–C: 0",23));
    let response=checked?(Number(guess.value)===m.value?"Your assignment agrees with the bond count.":`Recount: the assigned carbon value is ${m.value}.`):"Assign the highlighted carbon, then check your count.";
    say(ui,response,m.boundary+" These are oxidation states, not formal atomic charges. Next-oxidation maps summarise selected net transformations; actual reagents and conditions are not supplied as a protocol.");
    advance.disabled=!m.next;
    state(ui,{compound:compound.value,guess:Number(guess.value),checked,oxidationState:m.value,next:m.next,co:m.co,ch:m.ch});}
  update();
}

function mount14() {
  const ui=initDemo("day14");if(!ui)return;let stage=0,feedback="Choose the first electron-pair operation.";
  const action=selectControl(ui.controls,"Proposed next operation",[["ionise","C–Br pair → Br"],["waterAttack","Water lone pair → carbocation C"],["deprotonate","Water pair → oxonium H; O–H pair → O"],["neutralCapture","Water attacks and directly becomes neutral alcohol"],["eraseH","Delete the extra H without a recipient"],["methylCation","Replace the substrate with CH₃⁺"]],"ionise",update);
  const attempt=button(ui.controls,"Try this operation",()=>{const m=sn1Model(stage,action.value);feedback=(m.accepted?"Accepted. ":"Check again. ")+m.feedback;stage=m.nextStage;update();});
  button(ui.controls,"Reset mechanism",()=>{stage=0;action.value="ionise";feedback="Choose the first electron-pair operation.";update();});
  function update(){const m=sn1Model(stage);let b=text(250,185,m.centralCharge?"C⁺":"C",32);
    for(const [x,y] of [[100,175],[250,60],[250,290]])b+=line(250+(x-250)*.2,175+(y-175)*.2,250+(x-250)*.72,175+(y-175)*.72)+text(x,y+9,"CH₃",27);
    if(m.group)b+=line(275,175,370,175)+text(425,184,m.group,30,m.oxygenCharge?"#af3b45":"#153747");
    b+=text(570,70,m.br?"Br⁻ spectator":"Br still bonded",24)+text(570,120,`${m.waters} H₂O available`,24)+
      text(570,175,m.hydronium?"H₃O⁺ present":"No H₃O⁺ yet",24)+text(570,255,"Total charge = 0",27)+text(360,355,`State ${stage+1}/4: ${m.name}`,23);
    ui.visual.innerHTML=svg("SN1 species, charges and spectators",b,720,390);
    say(ui,feedback,"This is the selected tertiary aqueous substitution branch. Br⁻ remains after ionisation; the two explicit waters account for capture and deprotonation. E1 can compete in real conditions and is not calculated here. Only full-headed electron-pair operations are modelled.");
    attempt.disabled=stage===3;state(ui,{stage,action:action.value,feedback,totalCharge:m.totalCharge,brSpectator:m.br,waters:m.waters,hydronium:m.hydronium});}
  update();
}

function mount15() {
  const ui=initDemo("day15");if(!ui)return;
  const acid=numberControl(ui.controls,"Salicylic acid mass (g)",1.3812,0,20,.0001,update);
  const anh=numberControl(ui.controls,"Acetic anhydride mass (g)",1.53135,0,20,.00001,update);
  const isolate=numberControl(ui.controls,"Isolated material mass (g)",1.2,0,20,.01,update);
  button(ui.controls,"Anhydride-limited example",()=>{acid.value=1.3812;anh.value=.51045;isolate.value=.6;update();});
  button(ui.controls,"Zero acid",()=>{acid.value=0;update();});
  button(ui.controls,"Reset",()=>{acid.value=1.3812;anh.value=1.53135;isolate.value=1.2;update();});
  function update(){const inputs={acidMass:numerical(acid),anhydrideMass:numerical(anh),isolatedMass:numerical(isolate)},m=yieldModel(inputs.acidMass,inputs.anhydrideMass,inputs.isolatedMass);
    if(!m.valid){ui.visual.innerHTML=svg("Incomplete mass inputs",text(360,150,"Check the mass fields before calculating.",26));say(ui,m.error,"Inputs are masses in g. Empty and negative values are not silently converted into valid observations.");state(ui,{...inputs,valid:false,error:m.error});return;}
    const scale=Math.max(.001,m.acidMol,m.anhydrideMol),bar=(x,n,label)=>`<rect x="${x-45}" y="${280-n/scale*180}" width="90" height="${n/scale*180}" fill="#176b91"/>`+text(x,315,label,21)+text(x,75,`${fmt(n*1000,4)} mmol`,25);
    ui.visual.innerHTML=svg("Reactant amounts and maximum product amount",bar(145,m.acidMol,"acid")+bar(360,m.anhydrideMol,"anhydride")+bar(575,m.extent,"product max")+line(60,280,660,280));
    const warning=m.yieldPercent===null?"Yield is undefined because the theoretical product mass is zero.":`${fmt(m.yieldPercent,4)}% apparent isolated mass yield.`;
    say(ui,`Limiting reagent: ${m.limiting}. Theoretical aspirin mass: ${fmt(m.theoreticalMass,5)} g. ${warning}`,
      (m.aboveMaximum?"Above 100%: investigate residual solvent/water, impurities or input errors; this is not extra pure product beyond stoichiometry. ":"")+"Supplied balanced teaching equation: salicylic acid + acetic anhydride → aspirin + ethanoic acid (all coefficients 1). Molar masses: 138.12, 102.09 and 180.16 g/mol respectively. Yield is not conversion or purity. This is not a laboratory protocol.");
    state(ui,{...inputs,...m});}
  update();
}

function mount16() {
  const ui=initDemo("day16");if(!ui)return;let stage=0,checked=false;
  const route=selectControl(ui.controls,"Addition route",[["hydrate","Base-catalysed hydrate"],["hemiacetal","Acid-catalysed hemiacetal"]],"hydrate",()=>{stage=0;checked=false;update();});
  const target=selectControl(ui.controls,"Oxygen to inspect",[["a","A: original carbonyl oxygen"],["b","B: incoming nucleophile oxygen"],["c","C: proton-transfer partner oxygen"]],"a",()=>{checked=false;update();});
  const charge=selectControl(ui.controls,"My formal-charge assignment",[["-1","−1"],["0","0"],["1","+1"]],"0",()=>{checked=false;update();});
  const partner=selectControl(ui.controls,"Proton role in the operation entering this state",[["none","No proton transfer in this operation"],["waterDonates","Water C donates H⁺"],["hydroniumDonates","Hydronium C donates H⁺"],["waterAccepts","Water C accepts H⁺"]],"none",()=>{checked=false;update();});
  button(ui.controls,"Previous state",()=>{stage=Math.max(0,stage-1);checked=false;update();});
  button(ui.controls,"Next state",()=>{stage=Math.min(additionStates[route.value].length-1,stage+1);checked=false;update();});
  button(ui.controls,"Check labels and proton role",()=>{checked=true;update();});
  button(ui.controls,"Reset",()=>{route.value="hydrate";stage=0;target.value="a";charge.value="0";partner.value="none";checked=false;update();});
  const labelled=(letter,q,extra="")=>`O${letter}${q===-1?"⁻":q===1?"⁺":""}${extra}`;
  function update(){const m=additionModel(route.value,stage),isHemi=route.value==="hemiacetal";
    let top=labelled("ᴬ",m.a,m.form==="activated"||m.form==="tetra"&&m.a===0?"H":"");
    const bottom=m.form==="tetra"?(isHemi?labelled("ᴮ",m.b,m.b?"HCH₃":"CH₃"):"OᴮH"):null;
    let b=carbonSite({form:m.form,top,bottom})+text(560,80,`A charge: ${m.a}`,25)+text(560,130,`B charge: ${m.b}`,25)+text(560,180,`C charge: ${m.c}`,25)+text(560,245,m.side,24)+text(360,340,`State ${stage+1}: ${m.name}`,22);
    ui.visual.innerHTML=svg("Labelled carbonyl-addition species and proton partners",b,720,380);
    const ok=Number(charge.value)===m[target.value]&&partner.value===m.partner;
    const prefix=checked?(ok?"Both assignments agree. ":`Check: oxygen ${target.value.toUpperCase()} is ${m[target.value]}; the entering operation is ${m.partner}. `):"Inspect a labelled oxygen and choose its charge and the proton-transfer role. ";
    say(ui,prefix+m.explain,`Total charge of all explicitly tracked species is ${m.totalCharge}. Atom labels track actual oxygen identities; regenerating the same catalyst species does not require the incoming oxygen B to leave the product. Schemes are selected simplified pathways, not predictions of equilibrium proportions.`);
    state(ui,{route:route.value,stage,target:target.value,chargeGuess:Number(charge.value),partnerGuess:partner.value,checked,correct:checked&&ok,a:m.a,b:m.b,c:m.c,totalCharge:m.totalCharge});}
  update();
}

function mount17() {
  const ui=initDemo("day17");if(!ui)return;let stage=0,feedback="Choose the first operation in the supplied basic hydrolysis route.";
  const action=selectControl(ui.controls,"Proposed electron-pair operation",[["add","OH⁻ adds; C=O pair → O"],["collapse","O⁻ reforms C=O; OCH₃ departs"],["reverse","O⁻ reforms C=O; OH⁻ departs (reverse addition)"],["trap","Methoxide takes the acid OH proton"],["direct","Direct saturated-carbon-style replacement at acyl C"],["hydride","Collapse by expelling H⁻"],["methyl","Collapse by expelling CH₃⁻"],["stop","Stop with neutral acid + methoxide"]],"add",update);
  const attempt=button(ui.controls,"Try this operation",()=>{const m=acylModel(stage,action.value);feedback=(m.accepted?"Plausible step. ":"Reconsider. ")+m.feedback;stage=m.nextStage;update();});
  button(ui.controls,"Reset ester hydrolysis",()=>{stage=0;action.value="add";feedback="Choose the first operation in the supplied basic hydrolysis route.";update();});
  function update(){const m=acylModel(stage);let structure;
    if(stage===0)structure={form:"carbonyl",top:"O",right:"OCH₃"};
    else if(stage===1)structure={form:"tetra",top:"O⁻",right:"OCH₃",bottom:"OH"};
    else structure={form:"carbonyl",top:"O",right:stage===2?"OH":"O⁻"};
    ui.visual.innerHTML=svg("Acyl addition, collapse and proton-transfer trapping",carbonSite(structure)+text(570,105,m.side||"Still attached",25)+
      text(570,175,"Total charge −1",26)+text(570,245,stage===2?"Transient pair":stage===3?"Base-condition product":"Selected model state",23)+text(360,345,m.name,22));
    say(ui,feedback+" "+m.explain,"Supplied basic methyl-ethanoate hydrolysis model. Reverse addition from the tetrahedral intermediate is allowed. For the acid–base trap use compatible illustrative pKa values 4.8 for ethanoic acid and 15.5 for methanol: K≈10^(15.5−4.8)≈5×10¹⁰. Acid work-up, rates and real operating conditions are outside this explorer.");
    attempt.disabled=stage===3;state(ui,{stage,action:action.value,feedback,totalCharge:m.totalCharge,trapK:m.trapK});}
  update();
}

[mount10,mount11,mount12,mount13,mount14,mount15,mount16,mount17].forEach(mount=>mount());
})(typeof globalThis!=="undefined"?globalThis:window);
