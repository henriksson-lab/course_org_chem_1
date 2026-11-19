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
            history_results=[]
            for unit in units:
                expected_url=(ROOT/'docs'/unit['slug']/'lec.html').as_uri()
                call('Page.navigate',{'url':expected_url},session)
                for _ in range(100):
                    ready='location.href.split("#")[0] === '+json.dumps(expected_url)+' && Boolean(window.Reveal && Reveal.isReady())'
                    if evaluate(ready): break
                    time.sleep(.1)
                else: raise RuntimeError('Reveal did not initialise: '+unit['slug'])
                if evaluate('Reveal.getConfig().history') is not False:
                    raise RuntimeError('Slide changes must not add browser history entries: '+unit['slug'])
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
            # Reproduce the reported scaffold bug through the homepage's real
            # lecture link. Even after many slide changes, one browser Back must
            # return directly to that homepage, not to a previous slide.
            for slug in ['structures','substitution','spectroscopy']:
                home=(ROOT/'docs/index.html').as_uri()
                call('Page.navigate',{'url':home},session)
                for _ in range(100):
                    ready=evaluate('Boolean(document.querySelector(\'a[href$="'+slug+'/lec.html"]\'))')
                    if ready and evaluate('location.href')==home: break
                    time.sleep(.05)
                else: raise RuntimeError('Homepage lecture link unavailable: '+slug)
                before=call('Page.getNavigationHistory',session=session)['currentIndex']
                evaluate('document.querySelector(\'a[href$="'+slug+'/lec.html"]\').click(); true')
                for _ in range(100):
                    if slug+'/lec.html' in (evaluate('location.href') or '') and evaluate('Boolean(window.Reveal && Reveal.isReady())'): break
                    time.sleep(.05)
                else: raise RuntimeError('Could not enter lecture through homepage: '+slug)
                evaluate('Reveal.slide(2);Reveal.slide(7);Reveal.slide(12);true')
                time.sleep(.3)  # Allow Reveal's URL update to settle before Back.
                after=call('Page.getNavigationHistory',session=session)
                if after['currentIndex']!=before+1:
                    raise RuntimeError('Slide navigation added history entries: '+slug)
                call('Page.navigateToHistoryEntry',{'entryId':after['entries'][before]['id']},session)
                for _ in range(100):
                    if evaluate('location.href')==home: break
                    time.sleep(.05)
                else: raise RuntimeError('Browser Back did not return to course home: '+slug)
                history_results.append({'unit':slug,'browser_back_returns_home':True})
                print(f'{slug}: one browser Back returns to course home',flush=True)
            ws.close()
            (ROOT/'plan/layout-check.json').write_text(json.dumps(results,indent=2)+'\n')
            (ROOT/'plan/navigation-check.json').write_text(json.dumps(history_results,indent=2)+'\n')
            if any(r['overflows'] or r['missingImages'] for r in results):
                raise SystemExit('Browser layout check failed; inspect plan/layout-check.json')
            print('All lecture slides fit their canvas with all figures loaded.')
        finally:
            proc.terminate()
            try: proc.wait(timeout=5)
            except subprocess.TimeoutExpired: proc.kill();proc.wait()
