"""Optional browser check using installed Chrome and Python websocket-client.

Loads local course HTML, preloads every figure, checks slide bounds, and saves
representative screenshots in a chosen temporary output directory. No server
or external website is used. Run after make render.
"""
import base64
import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
import websocket

ROOT=Path(__file__).resolve().parents[1]
chrome=shutil.which('google-chrome') or shutil.which('chromium')
if not chrome: raise SystemExit('Install/use a supported local Chrome to run this optional check.')
units=json.loads((ROOT/'plan/units.json').read_text())
with tempfile.TemporaryDirectory(prefix='orgchem-browser-') as td:
    td=Path(td)
    with (td/'chrome.log').open('w') as log:
        proc=subprocess.Popen([chrome,'--headless','--no-sandbox','--disable-gpu',
            '--disable-dev-shm-usage','--disable-background-networking','--disable-extensions',
            '--disable-sync','--no-first-run','--no-default-browser-check',
            '--allow-file-access-from-files','--remote-debugging-port=0',
            '--window-size=1440,900','--user-data-dir='+str(td/'profile'),'about:blank'],
            stdout=log,stderr=log)
        try:
            address=None
            for _ in range(100):
                match=re.search(r'DevTools listening on (ws://\S+)',(td/'chrome.log').read_text())
                if match: address=match.group(1);break
                if proc.poll() is not None: break
                time.sleep(.1)
            if not address: raise RuntimeError('Chrome did not expose its local debugging connection')
            ws=websocket.create_connection(address,timeout=30,suppress_origin=True)
            sequence=0
            def call(method,params=None,session=None):
                global sequence
                sequence+=1
                req={'id':sequence,'method':method,'params':params or {}}
                if session: req['sessionId']=session
                ws.send(json.dumps(req))
                while True:
                    response=json.loads(ws.recv())
                    if response.get('id')==sequence:
                        if 'error' in response: raise RuntimeError(response['error'])
                        return response.get('result',{})
            target=call('Target.createTarget',{'url':'about:blank'})['targetId']
            session=call('Target.attachToTarget',{'targetId':target,'flatten':True})['sessionId']
            call('Page.enable',session=session)
            def evaluate(expression):
                result=call('Runtime.evaluate',{'expression':expression,'awaitPromise':True,
                      'returnByValue':True},session)
                if 'exceptionDetails' in result: raise RuntimeError(result['exceptionDetails'])
                return result.get('result',{}).get('value')
            results=[]
            for unit in units:
                call('Page.navigate',{'url':(ROOT/'docs'/unit['slug']/'lec.html').as_uri()},session)
                for _ in range(100):
                    if evaluate('Boolean(window.Reveal && Reveal.isReady())'): break
                    time.sleep(.1)
                else: raise RuntimeError('Reveal did not initialise: '+unit['slug'])
                result=evaluate('''(async()=>{
                  Reveal.configure({transition:'none',autoAnimate:false});
                  await Promise.all(Array.from(document.querySelectorAll('img')).map(img=>new Promise(resolve=>{
                    if(img.complete && img.naturalWidth>0){resolve();return;}
                    img.onload=resolve;img.onerror=resolve;
                    if(img.dataset.src) img.src=img.dataset.src;
                  })));
                  let overflows=[],missingImages=[];
                  for(const img of document.querySelectorAll('img')){
                    if(!img.naturalWidth)missingImages.push(img.src);
                  }
                  const slides=Reveal.getSlides();
                  for(let i=0;i<slides.length;i++){
                    Reveal.slide(i);Reveal.layout();
                    const slide=Reveal.getCurrentSlide();
                    const area=document.querySelector('.slides').getBoundingClientRect();
                    let bad=[];
                    for(const el of slide.querySelectorAll('h1,h2,h3,p,li,pre,table,img')){
                      if(el.closest('aside.notes') || !el.getClientRects().length)continue;
                      const b=el.getBoundingClientRect();
                      if(b.bottom>area.bottom+4 || b.right>area.right+4 || b.left<area.left-4){
                        bad.push({tag:el.tagName,text:(el.innerText||el.alt||'').slice(0,90)});
                      }
                    }
                    if(bad.length)overflows.push({slide:i+1,items:bad});
                  }
                  return {slides:slides.length,overflows,missingImages};
                })()''')
                result['unit']=unit['slug'];results.append(result)
                print(f'{unit["slug"]}: {result["slides"]} slides; '
                      f'{len(result["overflows"])} overflow(s); {len(result["missingImages"])} missing image(s)',flush=True)
                if unit['slug'] in ['acids-bases','stereochemistry','spectroscopy']:
                    evaluate('Reveal.slide(3);Reveal.layout();true')
                    time.sleep(.15)
                    shot=call('Page.captureScreenshot',{'format':'png'},session)
                    Path('/tmp/orgchem-'+unit['slug']+'-slide.png').write_bytes(base64.b64decode(shot['data']))
            ws.close()
            (ROOT/'plan/layout-check.json').write_text(json.dumps(results,indent=2)+'\n')
            if any(r['overflows'] or r['missingImages'] for r in results):
                raise SystemExit('Browser layout check failed; inspect plan/layout-check.json')
            print('All lecture slides fit their canvas with all figures loaded.')
        finally:
            proc.terminate()
            try: proc.wait(timeout=5)
            except subprocess.TimeoutExpired: proc.kill();proc.wait()
