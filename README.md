# AF Simulator

An interactive teaching tool for electrophysiologists. Atrial fibrillation is not animated by hand: it emerges from a live excitable-tissue simulation running on the surface of a 3D left atrium, so triggers, substrate, signals and therapy can all be explored directly.

Everything is in one self-contained file, `af-simulator.html`. There is no build step for users, no dependencies and no network access.

## Running it

Open `af-simulator.html` in a recent Chrome, Edge or Firefox (WebGL2 is required).

Optional URL parameters:

| Parameter | Effect |
|---|---|
| `scenario=NAME` | Start in `sinus`, `trigger`, `rotor`, `wavelets`, `pvi` or `gap` |
| `warm=MS` | Advance the tissue by this many milliseconds before showing it |
| `worker=0` | Run the simulation on the main thread instead of a Web Worker |

## What you can do

**Scenarios** (top bar): sinus rhythm, PV trigger to AF, single rotor, multi-wavelet AF, after PVI, PVI with a gap. Each one sets the tissue and starts live; nothing is replayed.

**Tools** (under the atrium): rotate, pace anywhere, move the mapping catheter, paint ablation lesions, seed a rotor.

**Panels** (right):

- **Rhythm:** playback speed, sinus drive, premature beats and bursts from any pulmonary vein.
- **Substrate:** electrical remodelling, repolarisation heterogeneity, fibrosis, with live conduction velocity, APD and wavelength readouts.
- **Therapy:** one-click lesion sets (PVI, roof, inferior, mitral isthmus, posterior box) with an optional gap, class I and class III drug effects, cardioversion.
- **Display:** voltage, phase and activation-time maps, rotor core markers, mesh resolution.

**Signals** (bottom): a far-field pseudo-ECG lead with ventricular complexes from a simple AV node rule, plus unipolar and bipolar electrograms from the catheter, with cycle length and dominant frequency.

## The model

- Monodomain reaction-diffusion with the Fenton-Karma three-variable cell model, explicit Euler, on a triangulated surface with a cotangent Laplacian.
- Procedural left atrium: an ellipsoidal body with four pulmonary vein stubs, the appendage, and openings for the vein ostia and mitral annulus.
- Pulmonary vein sleeves have a shorter action potential and slower conduction than the body.
- Electrograms are computed from transmembrane current over distance to each electrode.

Indicative tissue values on the standard mesh:

| State | Conduction velocity | APD | Wavelength |
|---|---|---|---|
| Healthy | 80 cm/s | 200 ms | 16 cm |
| Fully remodelled | 38 cm/s | 110 ms | 4.2 cm |

## Limitations

This is a teaching model, not a patient-specific or validated one.

- The wall has no thickness or fibre direction, and there is no right atrium.
- The action potential upstroke is broader than in real tissue so that it resolves on the mesh.
- Fibrillation cycle lengths (about 105 to 140 ms) are shorter than typical clinical AF.
- Catheter traces are auto-gained and compensated for reduced coupling, so they do not show low voltage in remodelled tissue.
- Scenarios are tuned on the standard mesh; on the coarse mesh AF induction does not sustain.
- The Web Worker path and hands-on interaction (painting lesions, dragging, clicking to pace) have only been checked through headless runs and screenshots.

## Development

The HTML file is generated. Edit the sources and rebuild:

```
dev/sim-core.js      model: mesh, cell model, pacing, electrodes, lesions, presets (no DOM)
dev/template.html    interface: WebGL renderer, controls, signal strips
dev/build.py         writes af-simulator.html from the two files above
```

```
python3 dev/build.py
```

`dev/check-presets.js` runs the model under Node for quick checks:

```
node dev/check-presets.js 48 cable      # CV, APD and wavelength for several drug/remodelling states
node dev/check-presets.js 48 trigger    # rotor cores / wavefronts / catheter cycle length every 500 ms
node dev/check-presets.js 48 drugNa     # induce AF, then apply full class I block
```

The first argument is the mesh frequency (32, 48 or 64); the second is a scenario name, `cable`, `drugK` or `drugNa`.
