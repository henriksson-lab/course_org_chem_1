"""Screenshot every slide of a rendered Reveal deck and report overflowing elements.

Usage: python3 tools/slide_screenshots.py docs/<unit>/lec.html <output-dir>
Fragments are shown (each slide is captured in its final state).
"""
import base64,json,re,shutil,subprocess,sys,tempfile,time
from pathlib import Path
import websocket
html=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]); out.mkdir(exist_ok=True)
td=Path(tempfile.mkdtemp(dir=out))
log=(td/'chrome.log').open('w')
proc=subprocess.Popen(['google-chrome','--headless','--no-sandbox','--disable-gpu','--hide-scrollbars','--allow-file-access-from-files','--remote-debugging-port=0','--window-size=1280,720','--user-data-dir='+str(td/'p'),'about:blank'],stdout=log,stderr=log)
for _ in range(200):
    m=re.search(r'DevTools listening on (ws://\S+)',(td/'chrome.log').read_text())
    if m: break
    time.sleep(.1)
ws=websocket.create_connection(m.group(1),timeout=60,suppress_origin=True); seq=[0]
def call(method,params=None,session=None):
    seq[0]+=1; r={'id':seq[0],'method':method,'params':params or {}}
    if session: r['sessionId']=session
    ws.send(json.dumps(r))
    while True:
        x=json.loads(ws.recv())
        if x.get('id')==seq[0]: return x.get('result',{})
t=call('Target.createTarget',{'url':'about:blank'})['targetId']
s=call('Target.attachToTarget',{'targetId':t,'flatten':True})['sessionId']
call('Emulation.setDeviceMetricsOverride',{'width':1280,'height':720,'deviceScaleFactor':1,'mobile':False},s)
ev=lambda e: call('Runtime.evaluate',{'expression':e,'awaitPromise':True,'returnByValue':True},s).get('result',{}).get('value')
call('Page.navigate',{'url':html.as_uri()},s)
for _ in range(200):
    if ev('Boolean(window.Reveal && Reveal.isReady())'): break
    time.sleep(.1)
n=ev('Reveal.getSlides().length'); time.sleep(1)
for i in range(n):
    ev(f'Reveal.slide({i},0,99);Reveal.layout();true'); time.sleep(.4)
    bad=ev('''(()=>{const area=document.querySelector('.slides').getBoundingClientRect();let b=[];for(const el of Reveal.getCurrentSlide().querySelectorAll('h2,p,li,table,img')){if(el.closest('aside.notes')||!el.getClientRects().length)continue;const r=el.getBoundingClientRect();if(r.bottom>area.bottom+4||r.right>area.right+4)b.push((el.innerText||el.alt||el.tagName).slice(0,50));}return b;})()''')
    print(i,ev('Reveal.getCurrentSlide().querySelector("h1,h2")?.innerText'),'OVERFLOW' if bad else '',bad or '')
    Path(out/f'{i:02d}.png').write_bytes(base64.b64decode(call('Page.captureScreenshot',{'format':'png'},s)['data']))
proc.kill()
