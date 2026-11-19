"""Local Chrome/CDP helper for the interactive-demo checks."""
import base64
import json
import re
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
import websocket

class Browser:
    def __enter__(self):
        chrome=shutil.which('google-chrome') or shutil.which('chromium')
        if not chrome: raise RuntimeError('A local Chrome browser is needed for this optional check')
        self.temp=tempfile.TemporaryDirectory(prefix='orgchem-demo-browser-')
        self.folder=Path(self.temp.name);self.log=(self.folder/'chrome.log').open('w')
        self.proc=subprocess.Popen([chrome,'--headless','--no-sandbox','--disable-gpu',
            '--disable-dev-shm-usage','--disable-background-networking','--disable-extensions',
            '--disable-sync','--no-first-run','--no-default-browser-check',
            '--remote-debugging-port=0','--window-size=1440,1000',
            '--user-data-dir='+str(self.folder/'profile'),'about:blank'],stdout=self.log,stderr=self.log)
        address=None
        for _ in range(100):
            match=re.search(r'DevTools listening on (ws://\S+)',(self.folder/'chrome.log').read_text())
            if match:address=match.group(1);break
            if self.proc.poll() is not None:break
            time.sleep(.1)
        if not address:raise RuntimeError('Local Chrome did not start')
        self.ws=websocket.create_connection(address,timeout=30,suppress_origin=True)
        self.sequence=0;self.exceptions=[]
        target=self.call('Target.createTarget',{'url':'about:blank'})['targetId']
        self.session=self.call('Target.attachToTarget',{'targetId':target,'flatten':True})['sessionId']
        self.call('Page.enable',session=True);self.call('Runtime.enable',session=True)
        self.call('Network.enable',session=True)
        # Prove the daily workbenches do not need a CDN or external service.
        self.call('Network.setBlockedURLs',{'urls':['http://*','https://*']},session=True)
        return self
    def call(self,method,params=None,session=False):
        self.sequence+=1
        request={'id':self.sequence,'method':method,'params':params or {}}
        if session:request['sessionId']=self.session
        self.ws.send(json.dumps(request))
        while True:
            response=json.loads(self.ws.recv())
            if response.get('method')=='Runtime.exceptionThrown':
                self.exceptions.append(response['params']['exceptionDetails'])
            if response.get('id')==self.sequence:
                if 'error' in response:raise RuntimeError(response['error'])
                return response.get('result',{})
    def evaluate(self,expression):
        result=self.call('Runtime.evaluate',{'expression':expression,'awaitPromise':True,
                                            'returnByValue':True},session=True)
        if 'exceptionDetails' in result:raise RuntimeError(result['exceptionDetails'])
        return result.get('result',{}).get('value')
    def navigate(self,url,ready='document.readyState === "complete"'):
        self.call('Page.navigate',{'url':url},session=True)
        expression='location.href.split("#")[0] === '+json.dumps(url)+' && ('+ready+')'
        for _ in range(150):
            try:
                if self.evaluate(expression):return
            except RuntimeError as error:
                if not any(s in str(error) for s in ['context','Context']):raise
            time.sleep(.05)
        raise RuntimeError('Page did not become ready: '+url)
    def viewport(self,width,height):
        self.call('Emulation.setDeviceMetricsOverride',{'width':width,'height':height,
                  'deviceScaleFactor':1,'mobile':width<600},session=True)
    def screenshot(self,path):
        data=self.call('Page.captureScreenshot',{'format':'png'},session=True)
        Path(path).write_bytes(base64.b64decode(data['data']))
    def __exit__(self,*args):
        self.ws.close();self.proc.terminate()
        try:self.proc.wait(timeout=5)
        except subprocess.TimeoutExpired:self.proc.kill();self.proc.wait()
        self.log.close();self.temp.cleanup()
