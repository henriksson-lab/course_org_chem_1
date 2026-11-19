// Chemical/numerical anchors and boundary checks supplied with each demo group.
// No browser or external packages are needed for these pure-model checks.
const path = require('node:path');
global.ChemDemo = require('../interactives/core.js');
let count = 0;
for (const group of ['foundations','reactivity','applications']) {
  const modulePath = path.resolve(__dirname,'../interactives/'+group+'.js');
  const api = require(modulePath);
  if (typeof api.checks !== 'function') throw new Error(group+': missing model checks');
  const results = api.checks();
  if (!Array.isArray(results) || !results.length) throw new Error(group+': empty model checks');
  for (const check of results) {
    if (!check.name || check.pass !== true) throw new Error(group+': '+(check.name||'unnamed check')+' failed');
    count++;
  }
  console.log(group+': '+results.length+' model checks passed');
}
console.log(count+' chemical/numerical model checks passed.');
