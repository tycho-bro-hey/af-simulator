# Assembles ../af-simulator.html from template.html and sim-core.js
import pathlib
here = pathlib.Path(__file__).parent
core = (here / 'sim-core.js').read_text()
assert '</script' not in core
out = (here / 'template.html').read_text().replace('/*__SIM_CORE__*/', core)
(here.parent / 'af-simulator.html').write_text(out)
print('wrote', here.parent / 'af-simulator.html', len(out), 'bytes')
