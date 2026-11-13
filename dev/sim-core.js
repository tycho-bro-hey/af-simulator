// AF simulator core: left atrial surface mesh + monodomain Fenton-Karma model.
// Pure JS, no DOM. Runs in a Web Worker, on the main thread (fallback) or under Node (tests).
const AFSim = (() => {
'use strict';

// ---------- small helpers ----------
function mulberry32(a) {
  return function () {
    a |= 0; a = a + 0x6D2B79F5 | 0;
    let t = Math.imul(a ^ a >>> 15, 1 | a);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
const norm3 = v => { const l = Math.hypot(v[0], v[1], v[2]); return [v[0] / l, v[1] / l, v[2] / l]; };
const cross3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot3 = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const clamp = (x, a, b) => x < a ? a : x > b ? b : x;
const smooth = x => { x = clamp(x, 0, 1); return x * x * x * (x * (6 * x - 15) + 10); };

// ---------- anatomy (mm; x = patient left, y = superior, z = anterior) ----------
const ANAT = {
  semi: [27, 25, 21],
  pv: { H: 11, tin: 0.18, tout: 0.40, hole: 0.155 },
  pvs: [
    { name: 'LSPV', axis: norm3([0.74, 0.42, -0.52]) },
    { name: 'LIPV', axis: norm3([0.74, -0.40, -0.54]) },
    { name: 'RSPV', axis: norm3([-0.72, 0.44, -0.54]) },
    { name: 'RIPV', axis: norm3([-0.72, -0.40, -0.56]) },
  ],
  laa: { name: 'LAA', axis: norm3([0.70, 0.30, 0.65]), H: 13, tin: 0.10, tout: 0.40 },
  mv: { name: 'Mitral annulus', axis: norm3([0.28, -0.62, 0.73]), hole: 0.52 },
  bachmann: norm3([-0.45, 0.55, 0.70]),
  avn: norm3([-0.55, -0.45, 0.50]),
  R: 24, // nominal radius used to convert angles to mm
};

function radiusAt(d) {
  const s = ANAT.semi;
  let r = 1 / Math.sqrt((d[0] / s[0]) ** 2 + (d[1] / s[1]) ** 2 + (d[2] / s[2]) ** 2);
  for (const p of ANAT.pvs) {
    const th = Math.acos(clamp(dot3(d, p.axis), -1, 1));
    r += ANAT.pv.H * (1 - smooth((th - ANAT.pv.tin) / (ANAT.pv.tout - ANAT.pv.tin)));
  }
  const th = Math.acos(clamp(dot3(d, ANAT.laa.axis), -1, 1));
  r += ANAT.laa.H * (1 - smooth((th - ANAT.laa.tin) / (ANAT.laa.tout - ANAT.laa.tin)));
  return r;
}

// ---------- mesh ----------
// Geodesic sphere; freq (edge subdivisions of the icosahedron) must be 2^k or 3 * 2^k.
function icosphere(freq) {
  let level = 0, three = freq % 3 === 0;
  for (let f = three ? freq / 3 : freq; f > 1; f >>= 1) level++;
  const t = (1 + Math.sqrt(5)) / 2;
  let verts = [[-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0], [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
    [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1]].map(norm3);
  let tris = [0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
    3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1];
  if (three) {
    const cache = new Map(), next = [];
    const third = (a, b) => { // point one third of the way from a to b
      const key = a * 1048576 + b;
      let m = cache.get(key);
      if (m === undefined) {
        const A = verts[a], B = verts[b];
        m = verts.length; verts.push(norm3([2 * A[0] + B[0], 2 * A[1] + B[1], 2 * A[2] + B[2]]));
        cache.set(key, m);
      }
      return m;
    };
    for (let i = 0; i < tris.length; i += 3) {
      const a = tris[i], b = tris[i + 1], c = tris[i + 2], A = verts[a], B = verts[b], C = verts[c];
      const ab1 = third(a, b), ab2 = third(b, a), bc1 = third(b, c), bc2 = third(c, b), ca1 = third(c, a), ca2 = third(a, c);
      const m = verts.length; verts.push(norm3([A[0] + B[0] + C[0], A[1] + B[1] + C[1], A[2] + B[2] + C[2]]));
      next.push(a, ab1, ca2, ab1, ab2, m, ab1, m, ca2, ab2, b, bc1, ab2, bc1, m, bc1, bc2, m, bc2, c, ca1, bc2, ca1, m, ca1, ca2, m);
    }
    tris = next;
  }
  for (let l = 0; l < level; l++) {
    const cache = new Map(), next = [];
    const mid = (a, b) => {
      const key = a < b ? a * 1048576 + b : b * 1048576 + a;
      let m = cache.get(key);
      if (m === undefined) {
        const A = verts[a], B = verts[b];
        m = verts.length; verts.push(norm3([A[0] + B[0], A[1] + B[1], A[2] + B[2]]));
        cache.set(key, m);
      }
      return m;
    };
    for (let i = 0; i < tris.length; i += 3) {
      const a = tris[i], b = tris[i + 1], c = tris[i + 2];
      const ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
      next.push(a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca);
    }
    tris = next;
  }
  return { verts, tris };
}

function buildMesh(freq) {
  const ico = icosphere(freq);
  const holes = ANAT.pvs.map(p => ({ axis: p.axis, cos: Math.cos(ANAT.pv.hole) }));
  holes.push({ axis: ANAT.mv.axis, cos: Math.cos(ANAT.mv.hole) });
  const inHole = d => holes.some(h => dot3(d, h.axis) > h.cos);

  // drop triangles inside the PV ostia and mitral annulus, then compact the vertex list
  const keepTri = [], remap = new Int32Array(ico.verts.length).fill(-1);
  for (let i = 0; i < ico.tris.length; i += 3) {
    const A = ico.verts[ico.tris[i]], B = ico.verts[ico.tris[i + 1]], C = ico.verts[ico.tris[i + 2]];
    if (inHole(norm3([A[0] + B[0] + C[0], A[1] + B[1] + C[1], A[2] + B[2] + C[2]]))) continue;
    keepTri.push(ico.tris[i], ico.tris[i + 1], ico.tris[i + 2]);
  }
  let n = 0;
  for (const v of keepTri) if (remap[v] < 0) remap[v] = n++;
  const dirs = new Float64Array(3 * n), pos = new Float32Array(3 * n);
  for (let o = 0; o < remap.length; o++) {
    const i = remap[o]; if (i < 0) continue;
    const d = ico.verts[o], r = radiusAt(d);
    dirs.set(d, 3 * i);
    pos[3 * i] = d[0] * r; pos[3 * i + 1] = d[1] * r; pos[3 * i + 2] = d[2] * r;
  }
  const idx = new Uint32Array(keepTri.length);
  for (let i = 0; i < idx.length; i++) idx[i] = remap[keepTri[i]];

  // cotangent Laplacian, lumped areas, vertex normals
  const area = new Float64Array(n), nrm = new Float32Array(3 * n), edge = new Map();
  let edgeSum = 0, edgeCnt = 0;
  const addEdge = (a, b, w) => {
    const key = a < b ? a * 1048576 + b : b * 1048576 + a;
    edge.set(key, (edge.get(key) || 0) + w);
  };
  for (let t = 0; t < idx.length; t += 3) {
    const ia = idx[t], ib = idx[t + 1], ic = idx[t + 2];
    const P = k => [pos[3 * k], pos[3 * k + 1], pos[3 * k + 2]];
    const A = P(ia), B = P(ib), C = P(ic);
    const ab = [B[0] - A[0], B[1] - A[1], B[2] - A[2]], ac = [C[0] - A[0], C[1] - A[1], C[2] - A[2]];
    const bc = [C[0] - B[0], C[1] - B[1], C[2] - B[2]];
    const nv = cross3(ab, ac), a2 = Math.hypot(nv[0], nv[1], nv[2]); // 2 * area
    for (const k of [ia, ib, ic]) { area[k] += a2 / 6; nrm[3 * k] += nv[0]; nrm[3 * k + 1] += nv[1]; nrm[3 * k + 2] += nv[2]; }
    // cot of the angle at each corner weights the opposite edge
    addEdge(ib, ic, 0.5 * dot3(ab, ac) / a2);
    addEdge(ia, ic, 0.5 * -dot3(ab, bc) / a2);
    addEdge(ia, ib, 0.5 * dot3(ac, bc) / a2);
    edgeSum += Math.hypot(...ab) + Math.hypot(...ac) + Math.hypot(...bc); edgeCnt += 3;
  }
  for (let i = 0; i < n; i++) {
    const l = Math.hypot(nrm[3 * i], nrm[3 * i + 1], nrm[3 * i + 2]) || 1;
    nrm[3 * i] /= l; nrm[3 * i + 1] /= l; nrm[3 * i + 2] /= l;
  }
  const deg = new Int32Array(n + 1);
  for (const key of edge.keys()) { deg[Math.floor(key / 1048576) + 1]++; deg[key % 1048576 + 1]++; }
  const rs = new Int32Array(n + 1);
  for (let i = 0; i < n; i++) rs[i + 1] = rs[i] + deg[i + 1];
  const col = new Int32Array(rs[n]), w0 = new Float64Array(rs[n]), fill = rs.slice(0, n);
  for (const [key, wRaw] of edge) {
    const a = Math.floor(key / 1048576), b = key % 1048576, w = Math.max(wRaw, 0.02);
    col[fill[a]] = b; w0[fill[a]++] = w; col[fill[b]] = a; w0[fill[b]++] = w;
  }

  // per-vertex anatomy: PV sleeve weight (0..1), which PV, LAA weight
  const sleeve = new Float32Array(n), pvId = new Int8Array(n).fill(-1), laa = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const d = [dirs[3 * i], dirs[3 * i + 1], dirs[3 * i + 2]];
    ANAT.pvs.forEach((p, k) => {
      const th = Math.acos(clamp(dot3(d, p.axis), -1, 1));
      const s = 1 - smooth((th - ANAT.pv.tin) / (ANAT.pv.tout - ANAT.pv.tin));
      if (s > sleeve[i]) { sleeve[i] = s; pvId[i] = k; }
    });
    const th = Math.acos(clamp(dot3(d, ANAT.laa.axis), -1, 1));
    laa[i] = 1 - smooth((th - ANAT.laa.tin) / (ANAT.laa.tout - ANAT.laa.tin));
  }

  const labelAt = (name, axis, lift) => {
    const r = radiusAt(axis) + lift;
    return { name, pos: [axis[0] * r, axis[1] * r, axis[2] * r], dir: axis };
  };
  const labels = ANAT.pvs.map(p => labelAt(p.name, p.axis, 3));
  labels.push(labelAt('LAA', ANAT.laa.axis, 3), labelAt('Mitral annulus', ANAT.mv.axis, -4),
    labelAt('Posterior wall', norm3([0, 0, -1]), 2), labelAt('Roof', norm3([0, 1, -0.15]), 2),
    labelAt('Septum', norm3([-1, -0.05, 0.25]), 2), labelAt("Bachmann's bundle", ANAT.bachmann, 2));

  return { n, pos, nrm, idx, dirs, area, rs, col, w0, sleeve, pvId, laa, labels, meanEdge: edgeSum / edgeCnt };
}

// ---------- tissue (generic graph: the atrial mesh or a 1D cable) ----------
// Fenton-Karma 3-variable model, atrial-like parameter set (times in ms).
const FK = {
  uc: 0.13, uv: 0.055, ucsi: 0.85, k: 10,
  tvp: 7, tvLow: 19.6, tvMid: 1000, twp: 360, twm: 70, // tvp is scaled with mesh spacing
  td: 0.27, t0: 8.3, tr: 50, tsi: 38,
  D: 0.33, // mm^2/ms at the reference spacing
};
const ACT_THR = 0.5;
// (1 + tanh(k (u - ucsi))) / 2 tabulated over u in [-0.25, 1.75)
let sgTab = null, sgKey = '';
function sigmoidTable() {
  const key = FK.k + ':' + FK.ucsi;
  if (key !== sgKey) {
    sgTab = new Float64Array(1024); sgKey = key;
    for (let i = 0; i < 1024; i++) sgTab[i] = 0.5 * (1 + Math.tanh(FK.k * ((i + 0.5) / 512 - 0.25 - FK.ucsi)));
  }
  return sgTab;
}

// Neighbours are stored padded to 6 per node (unused slots point at the node itself with weight 0)
// so the diffusion gather can be unrolled; slot of CSR entry k in row i is 6*i + (k - rs[i]).
function makeTissue(n, rs, col) {
  const f = () => new Float64Array(n), nb = new Int32Array(6 * n);
  for (let i = 0; i < n; i++) for (let s = 0; s < 6; s++) nb[6 * i + s] = rs[i] + s < rs[i + 1] ? col[rs[i] + s] : i;
  return {
    n, rs, col, nb, wgt: new Float64Array(6 * n), dead: new Uint8Array(n),
    gfi: f(), iwp: f(), u: f(), v: f().fill(1), w: f().fill(1), diff: f(),
    act: f().fill(-1e9), actPrev: f().fill(-1e9), t: 0, stims: [], itvp: 1 / FK.tvp,
  };
}

function restTissue(T) {
  T.u.fill(0); T.v.fill(1); T.w.fill(1); T.diff.fill(0);
  T.act.fill(-1e9); T.actPrev.fill(-1e9); T.t = 0; T.stims.length = 0;
}

function stepTissue(T, dt, nSteps) {
  const { n, nb, wgt, dead, gfi, iwp, u, v, w, diff, act, actPrev } = T;
  const uc = FK.uc, uv = FK.uv, sg = sigmoidTable();
  const itvp = T.itvp, itvLow = 1 / FK.tvLow, itvMid = 1 / FK.tvMid, itwm = 1 / FK.twm, it0 = 1 / FK.t0, itr = 1 / FK.tr, itsi = 1 / FK.tsi;
  for (let s = 0; s < nSteps; s++) {
    const t = T.t;
    if (T.stims.length) {
      for (let q = T.stims.length - 1; q >= 0; q--) {
        const st = T.stims[q];
        if (t >= st.until) { T.stims.splice(q, 1); continue; }
        const nodes = st.nodes, a = st.amp * dt;
        for (let j = 0; j < nodes.length; j++) { const i = nodes[j]; if (!dead[i] && u[i] < 1) u[i] += a; }
      }
    }
    for (let i = 0; i < n; i++) {
      if (dead[i]) continue;
      const ui = u[i], k = 6 * i;
      diff[i] = wgt[k] * (u[nb[k]] - ui) + wgt[k + 1] * (u[nb[k + 1]] - ui) + wgt[k + 2] * (u[nb[k + 2]] - ui)
        + wgt[k + 3] * (u[nb[k + 3]] - ui) + wgt[k + 4] * (u[nb[k + 4]] - ui) + wgt[k + 5] * (u[nb[k + 5]] - ui);
    }
    for (let i = 0; i < n; i++) {
      if (dead[i]) continue;
      const ui = u[i], vi = v[i], wi = w[i];
      let cur, dv, dw;
      if (ui >= uc) {
        cur = itr - vi * (1 - ui) * (ui - uc) * gfi[i];
        dv = -vi * itvp; dw = -wi * iwp[i];
      } else {
        cur = ui * it0;
        dv = (1 - vi) * (ui >= uv ? itvMid : itvLow); dw = (1 - wi) * itwm;
      }
      const si = (ui + 0.25) * 512;
      cur -= wi * sg[si < 0 ? 0 : si > 1023 ? 1023 : si | 0] * itsi;
      const un = ui + dt * (diff[i] - cur);
      if (ui < ACT_THR && un >= ACT_THR) { actPrev[i] = act[i]; act[i] = t; }
      u[i] = un; v[i] = vi + dt * dv; w[i] = wi + dt * dw;
    }
    T.t = t + dt;
  }
}

// Drug / remodelling multipliers shared by the mesh and the 1D reference cable.
// fs = mesh spacing / reference spacing: a coarser mesh gets a proportionally broader wavefront
// (td and D both scaled by fs), which keeps conduction velocity the same.
const H_REF = 0.535;
const KN = { na: 0.4, remodel: 0.6, k: 0.4, cx: 0.75, apdFreq: 3.5, fibDead: 0.4 };
function kinetics(p, fs) {
  return {
    gfi: (1 / (FK.td * fs)) * (1 - KN.na * p.naBlock),
    // plateau length follows the decay of the slow inward gate, so twp sets APD
    twpScale: (1 - KN.remodel * p.remodel) * (1 + KN.k * p.kBlock),
    // remodelling also uncouples cells (connexin loss), slowing conduction
    dScale: 1 - KN.cx * p.remodel,
  };
}

// 1D reference cable with the mesh's kinetics (fs) on a finer grid: reports CV, APD and wavelength
// for the current drug/remodelling state.
function measureCable(p, dx, dt, paceCL, fs) {
  const n = 220, rs = new Int32Array(n + 1), colA = [];
  for (let i = 0; i < n; i++) { if (i > 0) colA.push(i - 1); if (i < n - 1) colA.push(i + 1); rs[i + 1] = colA.length; }
  const T = makeTissue(n, rs, Int32Array.from(colA));
  const kin = kinetics(p, fs);
  for (let i = 0; i < n; i++) for (let k = rs[i]; k < rs[i + 1]; k++) T.wgt[6 * i + k - rs[i]] = FK.D * fs * kin.dScale / (dx * dx);
  T.itvp = 1 / (FK.tvp * fs); T.gfi.fill(kin.gfi); T.iwp.fill(1 / (FK.twp * kin.twpScale));
  const nodes = Int32Array.from([0, 1, 2, 3]), beats = 3, a = 70, b = 150, spm = Math.round(1 / dt);
  let cv = 0, apd = 0, up = -1;
  for (let beat = 0; beat < beats; beat++) {
    T.stims.push({ nodes, until: T.t + 2, amp: 0.5 });
    for (let ms = 0; ms < paceCL; ms++) {
      const before = T.u[110];
      stepTissue(T, dt, spm);
      if (beat === beats - 1) {
        if (before < ACT_THR && T.u[110] >= ACT_THR) up = T.t;
        if (up > 0 && !apd && before >= 0.1 && T.u[110] < 0.1) apd = T.t - up;
      }
    }
    if (beat === beats - 1) {
      const t0 = beat * paceCL;
      if (T.act[a] >= t0 && T.act[b] > T.act[a]) cv = (b - a) * dx / (T.act[b] - T.act[a]);
    }
  }
  return { cv, apd, wavelength: cv * apd }; // mm/ms, ms, mm
}

// ---------- the atrial simulation ----------
function createSim(freq) {
  const M = buildMesh(freq || 48);
  const n = M.n, T = makeTissue(n, M.rs, M.col);
  const dirOf = i => [M.dirs[3 * i], M.dirs[3 * i + 1], M.dirs[3 * i + 2]];
  const nearestToDir = d => { let best = 0, bd = -2; for (let i = 0; i < n; i++) { const c = M.dirs[3 * i] * d[0] + M.dirs[3 * i + 1] * d[1] + M.dirs[3 * i + 2] * d[2]; if (c > bd) { bd = c; best = i; } } return best; };
  const nodesNear = (p, radius) => {
    const out = [], r2 = radius * radius;
    for (let i = 0; i < n; i++) { const dx = M.pos[3 * i] - p[0], dy = M.pos[3 * i + 1] - p[1], dz = M.pos[3 * i + 2] - p[2]; if (dx * dx + dy * dy + dz * dz <= r2) out.push(i); }
    return Int32Array.from(out);
  };
  const posOf = i => [M.pos[3 * i], M.pos[3 * i + 1], M.pos[3 * i + 2]];

  const fs = M.meanEdge / H_REF;
  // explicit Euler: time step set by the stiffest node of the mesh
  let maxRate = 0;
  for (let i = 0; i < n; i++) { let r = 0; for (let k = M.rs[i]; k < M.rs[i + 1]; k++) r += M.w0[k]; r *= FK.D * fs / M.area[i]; if (r > maxRate) maxRate = r; }
  const spm = Math.max(4, Math.ceil(maxRate / 0.8)), dt = 1 / spm;
  const params = {
    naBlock: 0, kBlock: 0, remodel: 0, fibrosis: 0, fibSeed: 7, dispersion: 0.1,
    sinusOn: true, sinusCL: 800, qrs: true,
  };
  const ablated = new Uint8Array(n), fibrotic = new Uint8Array(n); // fibrotic: 1 = patch (slow), 2 = dead
  const undoStack = [];

  // smooth pseudo-random fields on the sphere (sum of plane waves)
  function noiseField(seed, freq) {
    const rnd = mulberry32(seed), K = 9, ws = [];
    for (let k = 0; k < K; k++) {
      const d = norm3([rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1]), f = freq * (0.6 + 0.8 * rnd());
      ws.push([d[0] * f, d[1] * f, d[2] * f, rnd() * 6.2832]);
    }
    const out = new Float32Array(n);
    let mx = 0;
    for (let i = 0; i < n; i++) {
      let s = 0;
      for (const q of ws) s += Math.cos(q[0] * M.dirs[3 * i] + q[1] * M.dirs[3 * i + 1] + q[2] * M.dirs[3 * i + 2] + q[3]);
      out[i] = s; if (Math.abs(s) > mx) mx = Math.abs(s);
    }
    for (let i = 0; i < n; i++) out[i] /= mx;
    return out;
  }
  const apdNoise = noiseField(1234, KN.apdFreq);
  let fibNoise = null, fibNoiseSeed = -1;

  function applySubstrate() {
    const kin = kinetics(params, fs);
    T.itvp = 1 / (FK.tvp * fs);
    // fibrosis: patches cover `fibrosis` of the surface; inside them ~40% of nodes are non-conducting
    fibrotic.fill(0);
    if (params.fibrosis > 0) {
      if (fibNoiseSeed !== params.fibSeed) { fibNoise = noiseField(params.fibSeed * 7919 + 13, 7); fibNoiseSeed = params.fibSeed; }
      const sorted = Float32Array.from(fibNoise).sort();
      const thr = sorted[Math.min(n - 1, Math.floor(n * (1 - params.fibrosis)))];
      const rnd = mulberry32(params.fibSeed * 31 + 5);
      for (let i = 0; i < n; i++) { const r = rnd(); if (fibNoise[i] >= thr) fibrotic[i] = r < KN.fibDead ? 2 : 1; }
    }
    for (let i = 0; i < n; i++) {
      T.dead[i] = (ablated[i] || fibrotic[i] === 2) ? 1 : 0;
      if (T.dead[i]) { T.u[i] = 0; T.v[i] = 1; T.w[i] = 1; T.diff[i] = 0; }
      const apdF = (1 - 0.30 * M.sleeve[i]) * (1 + params.dispersion * apdNoise[i]) * kin.twpScale;
      T.gfi[i] = kin.gfi * (1 - 0.15 * M.sleeve[i]);
      T.iwp[i] = 1 / (FK.twp * apdF);
    }
    for (let i = 0; i < n; i++) {
      const ci = fibrotic[i] === 1 ? 0.7 : 1, s = FK.D * fs * kin.dScale / M.area[i];
      for (let k = M.rs[i]; k < M.rs[i + 1]; k++) {
        const j = M.col[k];
        T.wgt[6 * i + k - M.rs[i]] = (T.dead[i] || T.dead[j]) ? 0 : s * M.w0[k] * ci * (fibrotic[j] === 1 ? 0.7 : 1);
      }
    }
    sigScale = 1 / kin.dScale;
    tissueDirty = true;
  }
  // electrode signals use the source strength without the coupling reduction, so remodelled tissue stays readable
  let tissueDirty = true, sigScale = 1;

  // ---------- pacing ----------
  const bachNode = nearestToDir(ANAT.bachmann), bachNodes = nodesNear(posOf(bachNode), 3);
  const avNode = nearestToDir(ANAT.avn);
  let lastSinus = -1e9;
  function stim(nodes, dur, amp) { T.stims.push({ nodes, until: T.t + (dur || 2), amp: amp || 0.5 }); }
  function paceAt(p, radius) { stim(nodesNear(p, radius || 2.5)); }
  // ectopic focus on the sleeve of each PV (posterior aspect of the tube)
  const pvFocus = ANAT.pvs.map(p => {
    const side = norm3(cross3(p.axis, [0, 1, 0])), a = 0.25;
    const d = norm3([p.axis[0] * Math.cos(a) + side[0] * Math.sin(a), p.axis[1] * Math.cos(a) + side[1] * Math.sin(a), p.axis[2] * Math.cos(a) + side[2] * Math.sin(a)]);
    return nodesNear(posOf(nearestToDir(d)), 2.5);
  });
  // one premature beat, timed from the last local activation of the sleeve (coupling interval, ms)
  function premature(pv, ci) {
    const last = T.act[pvFocus[pv][0]];
    schedule(Math.max(0, ci - (T.t - last)), () => stim(pvFocus[pv]));
  }
  const events = []; // {t, fn}
  function schedule(delay, fn) { events.push({ t: T.t + delay, fn }); events.sort((a, b) => a.t - b.t); }
  function ectopy(pv, count, cl, delay) {
    for (let b = 0; b < count; b++) schedule((delay || 0) + b * cl, () => stim(pvFocus[pv]));
  }

  // ---------- electrodes ----------
  // unipolar potential ~ sum over nodes of transmembrane current / distance
  function electrodeWeights(p) {
    const wts = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      const dx = M.pos[3 * i] - p[0], dy = M.pos[3 * i + 1] - p[1], dz = M.pos[3 * i + 2] - p[2];
      wts[i] = M.area[i] / Math.sqrt(dx * dx + dy * dy + dz * dz + 0.25);
    }
    return wts;
  }
  const leadAxis = norm3([0.5, -0.8, 0.2]);
  const ecgPos = electrodeWeights(leadAxis.map(c => c * 160)), ecgNeg = electrodeWeights(leadAxis.map(c => -c * 160));
  const ecgW = new Float64Array(n);
  for (let i = 0; i < n; i++) ecgW[i] = ecgPos[i] - ecgNeg[i];
  let cathW1 = null, cathW2 = null, cathNode = -1, cath = null;
  function setCatheter(p, nv) {
    const up = Math.abs(nv[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0], tg = norm3(cross3(nv, up));
    const e1 = [p[0] + nv[0] - 1.5 * tg[0], p[1] + nv[1] - 1.5 * tg[1], p[2] + nv[2] - 1.5 * tg[2]];
    const e2 = [p[0] + nv[0] + 1.5 * tg[0], p[1] + nv[1] + 1.5 * tg[1], p[2] + nv[2] + 1.5 * tg[2]];
    cathW1 = electrodeWeights(e1); cathW2 = electrodeWeights(e2);
    let best = 0, bd = 1e9;
    for (let i = 0; i < n; i++) { const d = (M.pos[3 * i] - p[0]) ** 2 + (M.pos[3 * i + 1] - p[1]) ** 2 + (M.pos[3 * i + 2] - p[2]) ** 2; if (d < bd && !T.dead[i]) { bd = d; best = i; } }
    cathNode = best; cath = { pos: p.slice(), nrm: nv.slice(), e1, e2 };
  }
  { const i = nearestToDir(norm3([0, 0.1, -1])); setCatheter(posOf(i), [M.nrm[3 * i], M.nrm[3 * i + 1], M.nrm[3 * i + 2]]); }

  // ---------- AV node + ventricular complexes for the pseudo-ECG ----------
  const avRnd = mulberry32(99);
  let avRefUntil = -1e9, avSeen = -1e9; const qrsTimes = [];
  function avInput(t) {
    if (t >= avRefUntil) {
      const early = Math.exp(-(t - avRefUntil) / 80);
      qrsTimes.push(t + 110 + 90 * early);
      avRefUntil = t + 330 + 120 * avRnd();
    } else avRefUntil = Math.min(avRefUntil + 28, t + 520); // concealed conduction
    while (qrsTimes.length && qrsTimes[0] < t - 700) qrsTimes.shift();
  }
  const g = (x, m, s) => Math.exp(-0.5 * ((x - m) / s) ** 2);
  function qrsAt(t) {
    let y = 0;
    for (const q of qrsTimes) {
      const x = t - q; if (x < -40 || x > 520) continue;
      y += -0.12 * g(x, -14, 6) + 1.0 * g(x, 0, 9) - 0.22 * g(x, 16, 8) + 0.26 * g(x, 250, 48);
    }
    return y;
  }

  // ---------- signals (1 kHz) ----------
  const ECG_GAIN = 60, EGM_GAIN = 3.5;
  const sig = { ecg: [], uni: [], bip: [] };
  function sample() {
    const diff = T.diff; let e = 0, a = 0, b = 0;
    for (let i = 0; i < n; i++) { const d = diff[i]; e += d * ecgW[i]; a += d * cathW1[i]; b += d * cathW2[i]; }
    sig.ecg.push(e * ECG_GAIN * sigScale + (params.qrs ? qrsAt(T.t) : 0)); sig.uni.push(a * EGM_GAIN * sigScale); sig.bip.push((a - b) * EGM_GAIN * sigScale);
  }

  function advance(ms) {
    for (let m = 0; m < ms; m++) {
      const t = T.t;
      if (params.sinusOn && t - Math.max(T.act[bachNode], lastSinus) >= params.sinusCL) { stim(bachNodes); lastSinus = t; }
      while (events.length && events[0].t <= t) events.shift().fn();
      stepTissue(T, dt, spm);
      if (T.act[avNode] !== avSeen) { avSeen = T.act[avNode]; avInput(avSeen); }
      sample();
    }
  }

  // ---------- state helpers ----------
  function reset() {
    restTissue(T); events.length = 0; lastSinus = -1e9; avRefUntil = -1e9; avSeen = -1e9; qrsTimes.length = 0;
    sig.ecg.length = sig.uni.length = sig.bip.length = 0;
    for (let i = 0; i < n; i++) T.dead[i] = (ablated[i] || fibrotic[i] === 2) ? 1 : 0;
  }
  function shock() { for (let i = 0; i < n; i++) if (!T.dead[i]) T.u[i] = 1; T.stims.length = 0; events.length = 0; }

  // great-arc geometry on the unit sphere of vertex directions
  function arcDist(d, A, B) { // angular distance from d to the arc A->B, plus position along it (0..1)
    const nrmv = norm3(cross3(A, B)), tot = Math.acos(clamp(dot3(A, B), -1, 1));
    const side = dot3(d, nrmv), inA = dot3(cross3(A, d), nrmv), inB = dot3(cross3(d, B), nrmv);
    if (inA >= 0 && inB >= 0) {
      const proj = norm3([d[0] - side * nrmv[0], d[1] - side * nrmv[1], d[2] - side * nrmv[2]]);
      return { dist: Math.abs(Math.asin(clamp(side, -1, 1))), side, along: Math.acos(clamp(dot3(A, proj), -1, 1)) / tot, inside: true };
    }
    const dA = Math.acos(clamp(dot3(d, A), -1, 1)), dB = Math.acos(clamp(dot3(d, B), -1, 1));
    return { dist: Math.min(dA, dB), side, along: dA < dB ? 0 : 1, inside: false };
  }

  // Seed a rotor: a broken wavefront along the arc tip->anchor, moving toward +side, refractory behind it.
  function seedRotor(tip, anchor) {
    for (let i = 0; i < n; i++) {
      if (T.dead[i]) continue;
      const a = arcDist(dirOf(i), tip, anchor); if (!a.inside) continue;
      const s = Math.asin(clamp(a.side, -1, 1)) * ANAT.R; // mm
      if (s >= 0 && s < 5) { T.u[i] = 1; T.v[i] = 0.9; T.w[i] = 0.9; }
      else if (s < 0 && s > -35) { const r = -s / 35; T.u[i] = 0; T.v[i] = r * r; T.w[i] = 0.3 + 0.7 * r; }
    }
  }

  // ---------- ablation ----------
  const HALF = 0.085; // lesion half-width (rad): lines about 4 mm wide
  const pairs = { left: [0, 1], right: [2, 3] };
  const ringOffset = ANAT.pv.tout + 0.07;
  function lesionNodes(name, gap) {
    const out = [], ax = k => ANAT.pvs[k].axis;
    const ringD = (d, side) => arcDist(d, ax(pairs[side][0]), ax(pairs[side][1])).dist;
    const insideRings = d => ringD(d, 'left') < ringOffset - HALF || ringD(d, 'right') < ringOffset - HALF;
    const line = (A, B, gapAt) => {
      for (let i = 0; i < n; i++) {
        const d = dirOf(i), a = arcDist(d, A, B);
        if (!a.inside || a.dist > HALF || insideRings(d)) continue;
        if (gap && gapAt !== undefined && Math.abs(a.along - gapAt) < 0.05) continue;
        out.push(i);
      }
    };
    const ring = (side, withGap) => {
      const A = ax(pairs[side][0]), B = ax(pairs[side][1]);
      // the gap sits on the posterior-wall side of the ring, level with the superior vein
      const toward = norm3([0, 0.25, -1]);
      for (let i = 0; i < n; i++) {
        const d = dirOf(i);
        if (Math.abs(arcDist(d, A, B).dist - ringOffset) > HALF) continue;
        if (withGap) { const gp = norm3([A[0] + 0.62 * toward[0], A[1] + 0.62 * toward[1], A[2] + 0.62 * toward[2]]); if (Math.acos(clamp(dot3(d, gp), -1, 1)) < 0.20) continue; }
        out.push(i);
      }
    };
    if (name === 'pviLeft') ring('left', gap);
    else if (name === 'pviRight') ring('right', gap);
    else if (name === 'roof') line(ax(0), ax(2), 0.5);
    else if (name === 'floor') line(ax(1), ax(3), 0.5);
    else if (name === 'mitral') line(ax(1), ANAT.mv.axis, 0.6);
    else if (name === 'box') { ring('left', false); ring('right', false); line(ax(0), ax(2), 0.5); line(ax(1), ax(3), 0.5); }
    return out;
  }
  function ablate(nodes) {
    const fresh = [];
    for (const i of nodes) if (!ablated[i]) { ablated[i] = 1; fresh.push(i); }
    if (fresh.length) { undoStack.push(fresh); applySubstrate(); }
    return fresh.length;
  }
  // extends the most recent undo entry so one brush stroke is one undo step
  function ablateStroke(p, radius, newStroke) {
    const fresh = [];
    for (const i of nodesNear(p, radius)) if (!ablated[i]) { ablated[i] = 1; fresh.push(i); }
    if (!fresh.length) return;
    if (newStroke || !undoStack.length) undoStack.push(fresh); else undoStack[undoStack.length - 1].push(...fresh);
    applySubstrate();
  }
  function undoAblation() { const last = undoStack.pop(); if (!last) return; for (const i of last) ablated[i] = 0; applySubstrate(); }
  function clearAblation() { ablated.fill(0); undoStack.length = 0; applySubstrate(); }

  // ---------- display ----------
  const phase = new Float32Array(n);
  function computePhase() { for (let i = 0; i < n; i++) phase[i] = Math.atan2(T.v[i] - 0.4, T.u[i] - 0.3); }
  function scalars(mode, out) {
    if (mode === 1) { computePhase(); for (let i = 0; i < n; i++) out[i] = phase[i] / 6.2832 + 0.5; }
    else if (mode === 2) { const t = T.t; for (let i = 0; i < n; i++) out[i] = t - T.act[i]; }
    else for (let i = 0; i < n; i++) out[i] = T.u[i];
    return out;
  }
  // phase singularities: triangles around which the phase winds by +-2pi; nearby opposite pairs cancel
  function detectPS() {
    computePhase();
    const idx = M.idx, raw = [], TWO = 6.2832, PI = Math.PI;
    const wrap = x => x > PI ? x - TWO : x < -PI ? x + TWO : x;
    for (let t = 0; t < idx.length; t += 3) {
      const a = idx[t], b = idx[t + 1], c = idx[t + 2];
      if (T.dead[a] || T.dead[b] || T.dead[c]) continue;
      const s = wrap(phase[b] - phase[a]) + wrap(phase[c] - phase[b]) + wrap(phase[a] - phase[c]);
      if (s > 3 || s < -3) raw.push([(M.pos[3 * a] + M.pos[3 * b] + M.pos[3 * c]) / 3, (M.pos[3 * a + 1] + M.pos[3 * b + 1] + M.pos[3 * c + 1]) / 3, (M.pos[3 * a + 2] + M.pos[3 * b + 2] + M.pos[3 * c + 2]) / 3, s > 0 ? 1 : -1]);
    }
    const kept = [];
    for (const p of raw) {
      let merged = false;
      for (let q = 0; q < kept.length; q++) {
        const k = kept[q];
        if ((k[0] - p[0]) ** 2 + (k[1] - p[1]) ** 2 + (k[2] - p[2]) ** 2 < 16) {
          k[3] += p[3]; merged = true; break;
        }
      }
      if (!merged) kept.push(p.slice());
    }
    return kept.filter(k => k[3] !== 0).map(k => [k[0], k[1], k[2], Math.sign(k[3])]);
  }
  // separate excited regions = wavefront count
  const seen = new Int32Array(n); let seenGen = 0; const stack = new Int32Array(n);
  function countWaves() {
    seenGen++; let count = 0;
    for (let s = 0; s < n; s++) {
      if (seen[s] === seenGen || T.dead[s] || T.u[s] < 0.6) continue;
      let top = 0, size = 0; stack[top++] = s; seen[s] = seenGen;
      while (top) {
        const i = stack[--top]; size++;
        for (let k = M.rs[i]; k < M.rs[i + 1]; k++) { const j = M.col[k]; if (seen[j] !== seenGen && !T.dead[j] && T.u[j] >= 0.6) { seen[j] = seenGen; stack[top++] = j; } }
      }
      if (size >= 8) count++;
    }
    return count;
  }
  function tissueFlags(out) { for (let i = 0; i < n; i++) out[i] = ablated[i] ? 3 : fibrotic[i]; return out; }
  function catheterCL() {
    if (cathNode < 0) return 0;
    const a = T.act[cathNode], b = T.actPrev[cathNode];
    return (b > -1e8 && T.t - a < 1500) ? a - b : 0;
  }

  // ---------- scenario presets ----------
  const PRESETS = {
    sinus: { p: { naBlock: 0, kBlock: 0, remodel: 0, fibrosis: 0, dispersion: 0.1, sinusOn: true, sinusCL: 800 }, lesions: [] },
    trigger: { p: { naBlock: 0, kBlock: 0, remodel: 1, fibrosis: 0.1, dispersion: 0.3, sinusOn: true, sinusCL: 700 }, lesions: [], run: () => ectopy(0, 12, 120, 1150) },
    rotor: { p: { naBlock: 0, kBlock: 0, remodel: 1, fibrosis: 0, dispersion: 0, sinusOn: false, sinusCL: 800 }, lesions: [], run: () => seedRotor(norm3([0, 0.1, -1]), ANAT.mv.axis) },
    wavelets: { p: { naBlock: 0, kBlock: 0, remodel: 1, fibrosis: 0, dispersion: 0.5, sinusOn: false, sinusCL: 800 }, lesions: [], run: () => { seedRotor(norm3([0, 0.1, -1]), ANAT.mv.axis); seedRotor(norm3([-0.3, 0.6, 0.7]), ANAT.pvs[2].axis); } },
    pvi: { p: { naBlock: 0, kBlock: 0, remodel: 1, fibrosis: 0.1, dispersion: 0.3, sinusOn: true, sinusCL: 700 }, lesions: [['pviLeft', false], ['pviRight', false]], run: () => ectopy(0, 12, 120, 1150) },
    gap: { p: { naBlock: 0, kBlock: 0, remodel: 1, fibrosis: 0.1, dispersion: 0.3, sinusOn: true, sinusCL: 700 }, lesions: [['pviLeft', true], ['pviRight', false]], run: () => ectopy(0, 12, 120, 1150) },
  };
  function preset(name) {
    const P = PRESETS[name]; if (!P) return;
    Object.assign(params, P.p);
    ablated.fill(0); undoStack.length = 0;
    for (const [l, gp] of P.lesions) for (const i of lesionNodes(l, gp)) ablated[i] = 1;
    applySubstrate(); reset();
    if (P.run) P.run();
  }

  applySubstrate();
  return {
    mesh: M, tissue: T, params, sig, dt,
    applySubstrate, advance, reset, shock, paceAt, ectopy, premature, stim, pvFocus, seedRotor, setCatheter,
    lesionNodes, ablate, ablateStroke, undoAblation, clearAblation,
    scalars, detectPS, countWaves, tissueFlags, catheterCL, preset, nodesNear, nearestToDir, posOf,
    get catheter() { return cath; },
    takeTissueDirty() { const d = tissueDirty; tissueDirty = false; return d; },
    measure(cl) { return measureCable(params, M.meanEdge / 2, dt, cl || 500, fs); },
    reseedFibrosis() { params.fibSeed = params.fibSeed % 97 + 1; applySubstrate(); },
    spm,
    lesionCount() { let c = 0; for (let i = 0; i < n; i++) c += ablated[i]; return c; },
  };
}

// ---------- message host: the same handler serves a Web Worker or the main thread ----------
function createHost(post, budgetMs) {
  let sim = null, carry = 0, frame = 0, measureAt = 0, lastPS = [], lastWaves = 0;
  const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());
  const SUBSTRATE = ['naBlock', 'kBlock', 'remodel', 'fibrosis', 'dispersion'];
  function init(freq) {
    sim = createSim(freq); carry = 0; measureAt = now();
    const M = sim.mesh;
    post({ type: 'mesh', n: M.n, pos: M.pos, nrm: M.nrm, idx: M.idx, labels: M.labels, meanEdge: M.meanEdge, freq, pvNames: ANAT.pvs.map(p => p.name) });
    post({ type: 'params', params: Object.assign({}, sim.params) });
  }
  const cmds = {
    preset(a) { sim.preset(a.name); post({ type: 'params', params: Object.assign({}, sim.params) }); measureAt = now(); },
    params(a) {
      let sub = false;
      for (const k in a.values) { sim.params[k] = a.values[k]; if (SUBSTRATE.includes(k)) sub = true; }
      if (sub) { sim.applySubstrate(); measureAt = now() + 250; }
    },
    pace(a) { sim.paceAt(a.p, 2.5); },
    catheter(a) { sim.setCatheter(a.p, a.n); },
    ectopy(a) { sim.ectopy(a.pv, a.count, a.cl, 0); },
    premature(a) { sim.premature(a.pv, a.ci); },
    ablateStroke(a) { sim.ablateStroke(a.p, a.r, a.fresh); },
    lesion(a) { sim.ablate(sim.lesionNodes(a.name, a.gap)); },
    undo() { sim.undoAblation(); },
    clear() { sim.clearAblation(); },
    shock() { sim.shock(); },
    rest() { sim.reset(); },
    reseed() { sim.reseedFibrosis(); },
    rotor(a) { sim.seedRotor(norm3(a.p), ANAT.mv.axis); },
    warm(a) { sim.advance(a.ms); },
  };
  return function handle(m) {
    if (m.type === 'init') return init(m.freq);
    if (!sim) return;
    if (m.type === 'cmd') return cmds[m.name](m.args || {});
    if (m.type !== 'tick') return;
    const t0 = now();
    carry = Math.min(carry + m.ms, 60);
    let done = 0;
    while (carry >= 1) { sim.advance(1); carry -= 1; done++; if (now() - t0 > budgetMs) { carry = 0; break; } }
    const out = {
      type: 'frame', t: sim.tissue.t, done, mode: m.mode, cost: now() - t0,
      val: sim.scalars(m.mode, new Float32Array(sim.mesh.n)),
      ecg: Float32Array.from(sim.sig.ecg), uni: Float32Array.from(sim.sig.uni), bip: Float32Array.from(sim.sig.bip),
      cl: sim.catheterCL(), catheter: sim.catheter, lesions: sim.lesionCount(),
    };
    sim.sig.ecg.length = sim.sig.uni.length = sim.sig.bip.length = 0;
    if (frame++ % 4 === 0) { lastPS = sim.detectPS(); lastWaves = sim.countWaves(); }
    out.ps = lastPS; out.waves = lastWaves;
    if (sim.takeTissueDirty()) out.flags = sim.tissueFlags(new Uint8Array(sim.mesh.n));
    if (measureAt && now() >= measureAt) { measureAt = 0; out.measure = sim.measure(400); }
    post(out);
  };
}

return { createSim, createHost, buildMesh, measureCable, ANAT, FK, KN, norm3 };
})();

if (typeof WorkerGlobalScope !== 'undefined' && self instanceof WorkerGlobalScope) {
  const handle = AFSim.createHost(m => self.postMessage(m), 45);
  self.onmessage = e => handle(e.data);
}

if (typeof module !== 'undefined' && module.exports) module.exports = AFSim;
