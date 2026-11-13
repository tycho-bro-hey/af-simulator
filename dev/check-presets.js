const A = require('./sim-core.js');
const [freq,which]=[+process.argv[2],process.argv[3]];
const sim = A.createSim(freq);
const snap=()=>sim.detectPS().length+'/'+sim.countWaves()+'/'+sim.catheterCL().toFixed(0);
if (which==='cable') { for (const p of [{remodel:0},{remodel:1},{remodel:1,naBlock:1},{remodel:1,kBlock:1},{remodel:0.5},{naBlock:1},{kBlock:1}]) { Object.assign(sim.params,{remodel:0,naBlock:0,kBlock:0},p); const m=sim.measure(400); console.log(freq,'spm',sim.spm,JSON.stringify(p),'cv',m.cv.toFixed(2),'apd',m.apd.toFixed(0),'wl',m.wavelength.toFixed(0)); } }
else if (which==='drugK'||which==='drugNa') { sim.preset('trigger'); let out=[]; for(let s=0;s<8;s++){ sim.advance(500); out.push(snap()); } sim.params[which==='drugK'?'kBlock':'naBlock']=1; sim.applySubstrate(); out.push('|drug|'); for(let s=0;s<12;s++){ sim.advance(500); out.push(snap()); } console.log(freq,which,out.join(' ')); }
else { sim.preset(which); let out=[]; for(let s=0;s<12;s++){ sim.advance(500); out.push(snap()); } console.log(freq,which,out.join(' ')); }
