"""Exercise all 25 demos offline, through real controls, on wide/narrow pages."""
import json
import time
from pathlib import Path
from browser_session import Browser

ROOT=Path(__file__).resolve().parents[1]
manifest=json.loads((ROOT/'plan/demos.json').read_text())
results=[];failures=[]

snapshot=r'''(()=>{
 const r=document.querySelector('.chem-demo');
 return {state:r.dataset.state||'',result:r.querySelector('.demo-result')?.textContent||'',
         visual:r.querySelector('.demo-visual')?.innerHTML||'',
         controls:Array.from(r.querySelectorAll('input,select')).map(e=>[e.tagName,e.type,e.value,e.checked])};
})()'''

layout=r'''(()=>{
 const root=document.querySelector('.chem-demo'), bad=[];
 for(const e of root.querySelectorAll('input,select,button')){
  const b=e.getBoundingClientRect();if(b.width>0&&(b.left< -2||b.right>innerWidth+2))bad.push(e.tagName+': control outside viewport');
  if(e.tagName!=='BUTTON'&&!e.labels?.length&&!e.getAttribute('aria-label')&&!e.getAttribute('aria-labelledby'))bad.push('Unlabelled form control');
 }
 for(const s of root.querySelectorAll('svg')){
  const box=s.viewBox.baseVal;if(!box.width)continue;
  for(const t of s.querySelectorAll('text')){
   const b=t.getBBox();if(b.x< -3||b.x+b.width>box.width+3||b.y< -3||b.y+b.height>box.height+3)
    bad.push('SVG label clipped: '+t.textContent);
  }
 }
 return {overflow:document.documentElement.scrollWidth>innerWidth+3,issues:bad};
})()'''

with Browser() as browser:
    for demo in manifest:
        browser.viewport(1440,1000)
        start_errors=len(browser.exceptions)
        url=(ROOT/'docs'/Path(demo['exercise_href']).with_suffix('.html')).as_uri()+'#explore'
        ready='document.getElementById('+json.dumps(demo['root_id'])+')?.dataset.ready === "true"'
        browser.navigate(url,ready)
        before=browser.evaluate(snapshot)
        if not before['visual'].strip() or not before['result'].strip():failures.append(f'Day {demo["day"]}: blank initial model')
        desktop=browser.evaluate(layout)
        changed=browser.evaluate(r'''(()=>{
          const r=document.querySelector('.chem-demo');
          for(const e of r.querySelectorAll('select,input')){
            if(e.disabled)continue;
            if(e.tagName==='SELECT' && e.options.length>1){e.selectedIndex=(e.selectedIndex+1)%e.options.length;e.dispatchEvent(new Event('change',{bubbles:true}));return true;}
            if(e.type==='range'||e.type==='number'){
              const v=Number(e.value),step=Number(e.step)||1,lo=e.min===''?0:Number(e.min),hi=e.max===''?v+step*3:Number(e.max);
              e.value=v+step<=hi?v+step:Math.max(lo,v-step);e.dispatchEvent(new Event('input',{bubbles:true}));return true;
            }
            if(e.type==='checkbox'){e.click();return true;}
          }
          const b=Array.from(r.querySelectorAll('button')).find(e=>!e.disabled&&!/reset/i.test(e.textContent));
          if(b){b.click();return true;}return false;
        })()''')
        after=browser.evaluate(snapshot)
        interactive=changed and any(before[k]!=after[k] for k in ['state','result','visual'])
        if changed and not interactive:
            # Some meaningful operations intentionally stage a choice before
            # applying it (for example, choosing which Fischer groups to swap).
            browser.evaluate('Array.from(document.querySelectorAll(".chem-demo button")).find(e=>!e.disabled&&!/reset/i.test(e.textContent))?.click();true')
            after=browser.evaluate(snapshot)
            interactive=any(before[k]!=after[k] for k in ['state','result','visual'])
        if not interactive:failures.append(f'Day {demo["day"]}: changing a control did not change the model/view')
        reset=browser.evaluate(r'''(()=>{
          const b=Array.from(document.querySelectorAll('.chem-demo button')).find(e=>/^(reset|defaults?)( |$)/i.test(e.textContent));
          if(b){b.click();return true;}return false;
        })()''')
        if not reset:failures.append(f'Day {demo["day"]}: no labelled reset/default control')
        restored=browser.evaluate(snapshot)
        if reset and restored['controls']!=before['controls']:failures.append(f'Day {demo["day"]}: reset did not restore control values')
        # Empty numeric inputs must not throw or expose a machine NaN/Infinity.
        numeric=browser.evaluate(r'''(()=>{
          const e=document.querySelector('.chem-demo input[type="number"]');
          if(!e)return false;e.value='';e.dispatchEvent(new Event('input',{bubbles:true}));return true;
        })()''')
        if numeric:
            text=browser.evaluate('document.querySelector(".chem-demo").innerText')
            if 'NaN' in text or 'Infinity' in text:failures.append(f'Day {demo["day"]}: invalid numeric input leaked a non-finite result')
            browser.evaluate('Array.from(document.querySelectorAll(".chem-demo button")).find(e=>/^(reset|defaults?)( |$)/i.test(e.textContent))?.click();true')
        # Inspect alternate cases and slider endpoints, where labels and geometry
        # often differ from the attractive default state.
        variants=[]
        select_count=browser.evaluate('document.querySelector(".chem-demo select")?.options.length||0')
        for index in sorted(set([0,max(0,select_count-1)])) if select_count else []:
            browser.evaluate(f'(()=>{{const e=document.querySelector(".chem-demo select");e.selectedIndex={index};e.dispatchEvent(new Event("change",{{bubbles:true}}));return true;}})()')
            variants.append(browser.evaluate(layout))
        for endpoint in ['min','max']:
            exists=browser.evaluate('Boolean(document.querySelector(".chem-demo input[type=range]"))')
            if exists:
                browser.evaluate('(()=>{const e=document.querySelector(".chem-demo input[type=range]");e.value=e.'+endpoint+';e.dispatchEvent(new Event("input",{bubbles:true}));return true;})()')
                variants.append(browser.evaluate(layout))
        browser.evaluate('Array.from(document.querySelectorAll(".chem-demo button")).find(e=>/^(reset|defaults?)( |$)/i.test(e.textContent))?.click();true')
        if demo['day']==3:
            key_before=browser.evaluate(snapshot)
            browser.evaluate('document.querySelector(".chem-demo input[type=range]").focus();true')
            browser.call('Input.dispatchKeyEvent',{'type':'keyDown','key':'ArrowRight','code':'ArrowRight','windowsVirtualKeyCode':39},session=True)
            browser.call('Input.dispatchKeyEvent',{'type':'keyUp','key':'ArrowRight','code':'ArrowRight','windowsVirtualKeyCode':39},session=True)
            if browser.evaluate(snapshot)['state']==key_before['state']:failures.append('Day 3: keyboard slider interaction failed')
        browser.viewport(390,900);time.sleep(.08)
        mobile=browser.evaluate(layout)
        if demo['day'] in [3,21]:
            browser.evaluate('document.querySelector(".chem-demo").scrollIntoView({block:"start"});true')
            browser.screenshot(f'/tmp/orgchem-demo-day{demo["day"]:02d}-mobile.png')
        for label,check in [('desktop',desktop),('mobile',mobile)]:
            if check['overflow'] or check['issues']:failures.append(f'Day {demo["day"]} {label}: '+json.dumps(check))
        for check in variants:
            if check['overflow'] or check['issues']:failures.append(f'Day {demo["day"]} alternate state: '+json.dumps(check))
        if len(browser.exceptions)>start_errors:failures.append(f'Day {demo["day"]}: browser script exception(s): '+str(browser.exceptions[start_errors:]))
        if demo['day'] in [3,6,12,18,21,22]:
            browser.viewport(1440,1000)
            browser.evaluate('document.querySelector(".chem-demo").scrollIntoView({block:"center"});true')
            browser.screenshot(f'/tmp/orgchem-demo-day{demo["day"]:02d}.png')
        results.append({'day':demo['day'],'initialised':True,'interaction_changed_view':interactive,
                        'reset_available':reset,'desktop':desktop,'mobile':mobile})
        print(f'Day {demo["day"]:02d}: loaded offline; interaction '+('passed' if interactive else 'FAILED'),flush=True)
    browser.viewport(1440,1000)
    browser.navigate((ROOT/'docs/index.html').as_uri())
    menu=browser.evaluate('Array.from(document.querySelectorAll("#quarto-sidebar a[href]"),e=>e.getAttribute("href"))')
    for demo in manifest:
        expected=demo['exercise_href'].replace('.qmd','.html')+'#explore'
        if not any(h.endswith(expected) for h in menu):failures.append(f'Day {demo["day"]}: missing rendered exercise/demo sidebar link')

    # Preserve already-published demo links while making Exercises canonical.
    for day in [1,15,25]:
        demo=manifest[day-1]
        alias=(ROOT/'docs'/Path(demo['href']).with_suffix('.html')).as_uri()
        expected=(ROOT/'docs'/Path(demo['exercise_href']).with_suffix('.html')).as_uri()
        browser.call('Page.navigate',{'url':alias},session=True)
        for _ in range(100):
            try:
                current=browser.evaluate('location.href') or ''
            except RuntimeError as error:
                if not any(s in str(error) for s in ['context','Context']):raise
                time.sleep(.05)
                continue
            if current==expected+'#explore':break
            time.sleep(.05)
        else:failures.append(f'Day {day}: old demo URL did not redirect to Exercises#explore')

(ROOT/'plan/demo-browser-check.json').write_text(json.dumps({'results':results,'failures':failures},indent=2)+'\n')
if failures:
    for failure in failures:print('ERROR:',failure)
    raise SystemExit(1)
print('All 25 demos respond, reset, fit desktop/mobile and have direct sidebar links with external requests blocked.')
