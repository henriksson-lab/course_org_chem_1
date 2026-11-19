(function(global) {
"use strict";
// Shared browser-only interface helpers. No network, packages or build step.
let sequence = 0;
const esc = value => String(value).replace(/[&<>"']/g,
  c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
function initDemo(id) {
  if (typeof document === 'undefined') return null;
  const root = document.getElementById(id);
  if (!root) return null;
  root.innerHTML = '<div class="demo-controls"></div><div class="demo-visual"></div><div class="demo-result" role="status" aria-live="polite"></div><div class="demo-detail"></div>';
  root.dataset.ready = 'true';
  return {root, controls:root.querySelector('.demo-controls'), visual:root.querySelector('.demo-visual'),
    result:root.querySelector('.demo-result'), detail:root.querySelector('.demo-detail')};
}
function field(container, label) {
  const wrap=document.createElement('div'); wrap.className='demo-field';
  const lab=document.createElement('label'); lab.textContent=label;
  lab.htmlFor='demo-control-'+(++sequence); wrap.append(lab); container.append(wrap);
  return {wrap,lab,id:lab.htmlFor};
}
function selectControl(container, label, options, initial, update) {
  const {wrap,id}=field(container,label), input=document.createElement('select');
  input.id=id;
  for(const option of options) {
    const [value,text]=Array.isArray(option)?option:[option,option];
    const item=document.createElement('option'); item.value=value; item.textContent=text; input.append(item);
  }
  input.value=initial; input.addEventListener('change',()=>update(input.value)); wrap.append(input); return input;
}
function rangeControl(container, label, min, max, step, value, unit, update) {
  const {wrap,lab,id}=field(container,label), input=document.createElement('input'), output=document.createElement('output');
  input.type='range'; Object.assign(input,{id,min,max,step,value}); output.htmlFor=id;
  function refresh() { output.textContent=input.value+(unit?' '+unit:''); }
  input.addEventListener('input',()=>{refresh();update(Number(input.value));});
  refresh(); lab.append(' ',output); wrap.append(input); return input;
}
function numberControl(container, label, value, min, max, step, update) {
  const {wrap,id}=field(container,label), input=document.createElement('input');
  input.type='number'; Object.assign(input,{id,value,min,max,step});
  input.addEventListener('input',()=>update(input.value.trim()===''?NaN:Number(input.value))); wrap.append(input); return input;
}
function button(container, label, action) {
  const input=document.createElement('button'); input.type='button'; input.textContent=label;
  input.addEventListener('click',action); container.append(input); return input;
}
function state(ui, values) { ui.root.dataset.state=JSON.stringify(values); }
function svg(title, body, width=720, height=360) {
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" style="max-width:${width}px" role="img" aria-label="${esc(title)}"><title>${esc(title)}</title><defs><marker id="chem-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10Z" fill="#176b91"/></marker></defs>${body}</svg>`;
}
const text = (x,y,label,size=26,color='#153747',anchor='middle') => `<text x="${x}" y="${y}" font-family="DejaVu Sans,sans-serif" font-size="${size}" text-anchor="${anchor}" fill="${color}">${esc(label)}</text>`;
const line = (x1,y1,x2,y2,color='#153747',width=3,arrow=false,dash='') => `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${color}" stroke-width="${width}"${arrow?' marker-end="url(#chem-arrow)"':''}${dash?' stroke-dasharray="'+dash+'"':''}/>`;
const circle = (cx,cy,r=5,fill='#176b91',stroke='none') => `<circle cx="${cx}" cy="${cy}" r="${r}" fill="${fill}" stroke="${stroke}" stroke-width="2"/>`;
const fmt = (v,digits=3) => !Number.isFinite(v)?'undefined':v!==0&&(Math.abs(v)<.001||Math.abs(v)>=10000)?v.toExponential(digits-1).replace(/\.0+(?=e)/,''):Number(v.toPrecision(digits)).toString();
function plot({title='Calculated teaching model',series,xDomain,yDomain,xLabel,yLabel,markers=[]}) {
  const legend=series.length>1, extra=legend?Math.max(0,Math.ceil(series.length/2)-1)*30:0;
  const W=720,H=410+extra,L=85,R=680,T=45,B=310;
  const x=v=>L+(v-xDomain[0])/(xDomain[1]-xDomain[0])*(R-L);
  const y=v=>B-(v-yDomain[0])/(yDomain[1]-yDomain[0])*(B-T);
  let body=line(L,B,R,B)+line(L,B,L,T);
  for(let i=0;i<5;i++) {
    let xv=xDomain[0]+i*(xDomain[1]-xDomain[0])/4, yv=yDomain[0]+i*(yDomain[1]-yDomain[0])/4;
    body+=line(x(xv),B,x(xv),B+6)+text(x(xv),B+32,fmt(xv),23);
    body+=text(L-12,y(yv)+7,fmt(yv,2),21,'#153747','end');
  }
  body+=text((L+R)/2,390+extra,xLabel,25)+text(L,26,yLabel,24,'#153747','start');
  for(const [index,s] of series.entries()) {
    const pts=s.points.map(p=>`${x(p[0])},${y(p[1])}`).join(' ');
    body+=`<polyline points="${pts}" fill="none" stroke="${s.color||'#176b91'}" stroke-width="3"${s.dash?' stroke-dasharray="'+s.dash+'"':''}/>`;
    if(legend) {
      const lx=85+(index%2)*300, ly=355+Math.floor(index/2)*30;
      body+=line(lx,ly,lx+25,ly,s.color||'#176b91',3,false,s.dash||'')+text(lx+35,ly+7,s.label||'Series '+(index+1),20,s.color||'#176b91','start');
    }
  }
  for(const m of markers)body+=circle(x(m.x),y(m.y),7,m.color||'#af3b45');
  return svg(title,body,W,H);
}

const api={esc,initDemo,selectControl,rangeControl,numberControl,button,state,svg,text,line,circle,fmt,plot};
global.ChemDemo=api;
if(typeof module!=="undefined" && module.exports) module.exports=api;
})(typeof globalThis!=="undefined"?globalThis:window);
