/* =======================================================================
   Phyto Poseidon — 3D hydrotropism simulator
   Two potted plants in transparent pots on either side of an opaque
   partition. Pot A has evenly moist soil (control); pot B has a water
   source on one side (treatment). Each root system bends toward the
   wetter side in the vertex shader.
   ======================================================================= */

const $ = (id) => document.getElementById(id);
const statusEl = $("stageStatus");
let loadDone = 0;
const LOAD_STEPS = 7;
function setProgress(pct, step) {
  const bar = $("loaderBar"), lab = $("loaderStep"), box = $("loader");
  if (bar) bar.style.width = pct + "%";
  if (lab && step) lab.textContent = step;
  if (box) box.setAttribute("aria-valuenow", String(Math.round(pct)));
}
function track(promise, name) {
  return promise.then((v) => { loadDone++; setProgress(16 + (loadDone / LOAD_STEPS) * 80, `${name} ready`); return v; });
}
const buzz = (ms = 10) => { try { navigator.vibrate && navigator.vibrate(ms); } catch (e) { /* not supported */ } };
function fail(err) {
  console.error("[PhytoPoseidon]", err);
  statusEl.textContent = location.protocol === "file:"
    ? "Open this page through a web server (for example GitHub Pages), not as a local file."
    : "The 3D scene needs WebGL and an internet connection.";
}

/* ---------------- Simulation model ----------------
   Hydro-gravitropic model (same form as the photo-gravitropic model of
   Bastien, Douady & Moulia 2015, with water in place of light).
   θ is the root tip angle from vertical (+ right), φ the direction of the
   water source (±90°, horizontal), G the moisture-gradient gain.
   u     = G·sin(φ − θ)              lateral moisture signal across the tip,
                                      from the gradient sensed one lag earlier
   ds/dt = (u − s) / τs              cortical growth asymmetry builds and decays
   dθ/dt = Kh·s − Kg(G)·sin θ        dry flank outgrows the wet flank; gravity
                                      pulls back toward vertical (sine law)
   A moisture gradient also breaks down the statolith amyloplasts of the root
   cap columella, so gravisensing weakens while the gradient lasts:
   dA/dt = (A* − A)/τ,  A* = 1 − 0.85·G,  Kg(t) = Kg·A
   (A falls to < 60 % in 30 min and ~15 % in 4 h under a strong gradient,
   Takahashi et al. 2003; τ = 0.8 h for loss, 6 h for regrowth, illustrative). */
const LAG       = 0.5;     // h; hydrotropic curvature starts within ~30 min
const MAX_BEND  = 90;      // degrees
const MAX_HOURS = 120;     // 5 days
const PLAY_DURATION_MS = 60000;  // full 120 h run takes 60 s (2 h per second)
const PSI_WET   = -0.05;   // MPa, water potential at the wet side
const DPSI_FULL = 1.0;     // MPa; Δψ at which the response saturates

const toRad = (d) => (d * Math.PI) / 180;
const smoothstep = (x) => { const t = Math.max(0, Math.min(1, x)); return t * t * (3 - 2 * t); };
const gradGain = (dpsi) => Math.max(0, Math.min(1, dpsi / DPSI_FULL));
const curLag = () => LAG;
const KH = 0.25, KG = 0.05, DT = 0.05, TAU_S = 1.0;   // rad h⁻¹, rad h⁻¹, h, h
const AMYLO_LOSS = 0.85, TAU_AMY_LOSS = 0.8, TAU_AMY_REGROW = 6;   // fraction, h, h
const amyTarget = (G) => 1 - AMYLO_LOSS * G;
// Flank growth from the turning rate (kinematics of a bending cylinder):
// the two flanks differ in elongation by ω·d, so each is ±ω·d/(2v) of the mean.
const ROOT_D = 0.12, ROOT_V = 0.2;        // mm between cortex flanks; mm h⁻¹ elongation (Arabidopsis primary root)
const flankPct = (rateDegH) => (100 * toRad(rateDegH) * ROOT_D) / (2 * ROOT_V);
const equilibriumDeg = (phiDeg, dpsi) => {
  const G = gradGain(dpsi);
  if (G <= 0) return 0;
  const p = KH * G, f = toRad(phiDeg), kg = KG * amyTarget(G);
  return (Math.atan2(p * Math.sin(f), kg + p * Math.cos(f)) * 180) / Math.PI;
};
// Settings are logged against simulated time, so the history is replayed exactly.
const events = [];
function settingsAt(h) {
  let lo = 0, hi = events.length - 1;
  while (lo < hi) { const mid = (lo + hi + 1) >> 1; if (events[mid].t <= h + 1e-9) lo = mid; else hi = mid - 1; }
  return events[lo];
}
// Five seedlings per room with natural variation (fixed, so runs are repeatable).
const SEEDLINGS = [
  { kp: 1.00, kg: 1.00, lag: 1.00 },
  { kp: 1.24, kg: 0.86, lag: 0.80 },
  { kp: 0.76, kg: 1.18, lag: 1.30 },
  { kp: 1.10, kg: 1.22, lag: 0.90 },
  { kp: 0.88, kg: 0.82, lag: 1.15 },
];
const traj = { A: null, B: null, Bset: [], Aset: [], Bmean: null, Bsd: null, Amean: null, Asd: null, rateB: null, amyB: null, dirty: true };
function runSeedling(stimOf, sp, amyOut = null) {
  const N = Math.round(MAX_HOURS / DT), out = new Float32Array(N + 1);
  let th = 0, s = 0, A = 1;
  if (amyOut) amyOut[0] = 1;
  for (let i = 1; i <= N; i++) {
    const h = i * DT;
    const lag = LAG * sp.lag;
    let u = 0, G = 0;
    if (h > lag) {
      const seen = stimOf(settingsAt(h - lag));
      G = gradGain(seen.dpsi);
      u = G * Math.sin(toRad(seen.phi) - th);
    }
    const At = amyTarget(G);
    A += ((At - A) / (At < A ? TAU_AMY_LOSS : TAU_AMY_REGROW)) * DT;
    s += ((u - s) / TAU_S) * DT;
    th += (KH * sp.kp * s - KG * A * sp.kg * Math.sin(th)) * DT;
    out[i] = (th * 180) / Math.PI;
    if (amyOut) amyOut[i] = A;
  }
  return out;
}
function stats(set) {
  const n = set.length, N = set[0].length, mean = new Float32Array(N), sd = new Float32Array(N);
  for (let i = 0; i < N; i++) {
    let m = 0; for (const a of set) m += a[i]; m /= n;
    let v = 0; for (const a of set) v += (a[i] - m) ** 2;
    mean[i] = m; sd[i] = Math.sqrt(v / (n - 1));
  }
  return { mean, sd };
}
function integrate() {
  traj.Aset = SEEDLINGS.map((sp) => runSeedling(() => ({ phi: 0, dpsi: 0 }), sp));   // even moisture: no gradient
  traj.amyB = new Float32Array(Math.round(MAX_HOURS / DT) + 1);
  traj.Bset = SEEDLINGS.map((sp, i) => runSeedling((ev) => ({ phi: 90 * ev.side, dpsi: ev.dpsi }), sp, i === 0 ? traj.amyB : null));
  const a = stats(traj.Aset), b = stats(traj.Bset);
  traj.Amean = a.mean; traj.Asd = a.sd; traj.Bmean = b.mean; traj.Bsd = b.sd;
  traj.A = traj.Aset[0]; traj.B = traj.Bset[0];
  // Turning rate of B (° h⁻¹), for the flank-growth view.
  const src = sim.reps === 5 ? traj.Bmean : traj.B;
  const r = new Float32Array(src.length);
  const w = Math.round(1 / DT);   // ±1 h window
  for (let i = 0; i < src.length; i++) { const a = Math.max(0, i - w), b = Math.min(src.length - 1, i + w); r[i] = b > a ? (src[b] - src[a]) / ((b - a) * DT) : 0; }
  traj.rateB = r;
  traj.dirty = false;
}
const angleAt = (arr, h) => {
  const x = Math.max(0, Math.min(MAX_HOURS, h)) / DT, i = Math.floor(x), f = x - i;
  return i >= arr.length - 1 ? arr[arr.length - 1] : arr[i] * (1 - f) + arr[i + 1] * f;
};
// The trajectory a room shows: seedling 1, or the mean of five.
const showA = () => (sim.reps === 5 ? traj.Amean : traj.A);
const showB = () => (sim.reps === 5 ? traj.Bmean : traj.B);
function logSettings() {
  const h = sim.hours;
  while (events.length > 1 && events[events.length - 1].t > h + 1e-9) events.pop(); // a change rewrites the future
  const last = events[events.length - 1];
  const rec = { side: sim.side, dpsi: sim.dpsi };
  if (events.length > 1 && h - last.t < 1.0) Object.assign(last, rec);   // one drag = one change
  else if (Math.abs(last.t - h) < 1e-6) Object.assign(last, rec);
  else events.push({ t: h, ...rec });
  traj.dirty = true;
}
// Roots grow down through the pot over the run; the shoot grows a little.
const rootGrowth  = (h) => 0.3 + 0.7 * Math.pow(Math.max(0, Math.min(1, h / MAX_HOURS)), 0.7);
const shootGrowth = (h) => 0.82 + 0.18 * Math.max(0, Math.min(1, h / MAX_HOURS));

const sim = {
  hours: 0,
  playing: false,
  playStart: 0,
  side: 1,                       // water source of pot B: −1 left, +1 right
  dpsi: 0.6,                     // MPa, effective water potential difference across pot B (0 when the water is off)
  dpsiSet: 0.6,                  // MPa, slider setting
  water: true,                   // water source on/off
  reps: 1,                       // 1 plant or 5 seedlings per room
};

events.push({ t: 0, side: sim.side, dpsi: sim.dpsi });

/* ---------------- Scene layout (world units) ---------------- */
const MODEL_SCALE  = 2.0;
const PLANT_ORIGIN = [0.125, 0.236, 0.065];  // stem base of the avocado, model units
const POT_X   = 1.16;      // pot centres at ±POT_X
const POT_H   = 1.5;      // round pot, scaled from the model (191 units tall)
const POT_K   = POT_H / 191.3;
const POT_RIM = 117.4 * POT_K;               // outer rim radius
const POT_FLOOR = 0.1 * POT_H;               // inside floor height
// Inside wall radius, linear in height (model: 79.6 at 10 %, 101.8 at 90 %).
const innerR = (y) => (79.6 + (101.8 - 79.6) * ((y / POT_H - 0.1) / 0.8)) * POT_K;
const SOIL_Y  = 1.36;      // soil surface (root crown) height
const ROOT_DEPTH = 1.12;    // full-grown root system depth
const BEND_ZONE_START = 0.2;   // fraction of the root depth that stays straight

/* ---------------- 3D ---------------- */
let THREE, GLTFLoader, OrbitControls, RoundedBoxGeometry, RoomEnvironment;

// Model files, relative to the site root. A host can override any of them
// by setting window.PHYTO_MODELS before this script loads.
const MODELS = Object.assign({
  glbPlant: "assets/models/avocado.glb",
  glbRoots: "assets/models/roots.glb",
  glbPot: "assets/models/pot.glb",
  glbWater: "assets/models/water.glb",
  glbMag: "assets/models/magnifier.glb",
  glbTip: "assets/models/root_tip.glb",
  glbMag: "assets/models/magnifier.glb",
}, window.PHYTO_MODELS || {});
const modelURL = (path) => (/^(\/|https?:)/.test(path) ? path : new URL(`../${path}`, import.meta.url).href);
const modelCache = {};
function b64ToBuffer(b64) {
  const bin = atob(b64), out = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) out[i] = bin.charCodeAt(i);
  return out.buffer;
}
function loadBuffer(key) {
  const inline = window.PHYTO_B64 && window.PHYTO_B64[key];
  if (!modelCache[key] && inline) modelCache[key] = Promise.resolve(b64ToBuffer(inline));
  if (!modelCache[key]) {
    modelCache[key] = fetch(modelURL(MODELS[key])).then((r) => {
      if (!r.ok) throw new Error(`Could not load ${MODELS[key]} (${r.status})`);
      return r.arrayBuffer();
    });
  }
  return modelCache[key];
}
// Each call gets its own copy, since a parsed buffer can't be shared by two scenes.
const parseGLB = (key) => loadBuffer(key).then((buf) => new Promise((res, rej) => new GLTFLoader().parse(buf.slice(0), "", res, rej)));

function toFloatAttr(attr) {
  const n = attr.count, s = attr.itemSize, arr = new Float32Array(n * s);
  for (let i = 0; i < n; i++) {
    arr[i * s] = attr.getX(i);
    if (s > 1) arr[i * s + 1] = attr.getY(i);
    if (s > 2) arr[i * s + 2] = attr.getZ(i);
    if (s > 3) arr[i * s + 3] = attr.getW(i);
  }
  return new THREE.BufferAttribute(arr, s);
}

// Shoot only (the model's own pot and soil are replaced by the glass pot):
// stem base at the origin, y up.
function bakeShoot(root) {
  root.updateMatrixWorld(true);
  const shoot = [];
  root.traverse((o) => {
    if (!o.isMesh) return;
    const name = (o.material && o.material.name) || "";
    if (/gorshok|material_0/i.test(name)) return;
    const g = new THREE.BufferGeometry();
    ["position", "normal", "uv", "uv1", "color"].forEach((k) => { if (o.geometry.attributes[k]) g.setAttribute(k, toFloatAttr(o.geometry.attributes[k])); });
    if (o.geometry.index) g.setIndex(o.geometry.index);
    g.applyMatrix4(o.matrixWorld);
    g.translate(-PLANT_ORIGIN[0], -PLANT_ORIGIN[1], -PLANT_ORIGIN[2]);
    g.scale(MODEL_SCALE, MODEL_SCALE, MODEL_SCALE);
    if (/stebel|vetvi/i.test(name)) g.computeVertexNormals(); // simplified bark: smooth shading
    shoot.push({ geometry: g, material: o.material });
  });
  return shoot;
}

function makeShoot(parts) {
  const g = new THREE.Group();
  parts.forEach((p) => {
    const mat = p.material.clone();
    mat.side = THREE.DoubleSide;
    const mesh = new THREE.Mesh(p.geometry, mat);
    mesh.castShadow = true; mesh.receiveShadow = true;
    g.add(mesh);
  });
  return g;
}

/* Root system: crown at the origin, depth 1 (scaled at runtime). Below the
   straight zone each horizontal slice of the root system is carried along a
   circular arc (the path of the primary root), whose tip angle is uBend
   (+ toward +x). Slices keep their shape, so the laterals sweep with the
   primary root toward the water instead of rotating around it. */
// Natural bend: the turning angle grows with depth, a(s) = θ·(s/L)^1.6, so
// the upper root stays nearly straight and the younger, deeper parts (the
// tips) curve toward the water. Each slice is carried along that path.
const BEND_POW = 1.6, BEND_STEPS = 12;
const BEND_GLSL_HEAD = "uniform float uBend;\nuniform float uY0;\nuniform float uTime;\nuniform float uPulse;\nvarying float vD;\n";
const BEND_NORMAL = `#include <beginnormal_vertex>
  {
    float nu = clamp((-position.y - uY0) / (1.0 - uY0), 0.0, 1.0);
    float na = 0.6 * uBend * pow(nu, ${BEND_POW.toFixed(2)});
    objectNormal = vec3(objectNormal.x * cos(na) - objectNormal.y * sin(na),
                        objectNormal.x * sin(na) + objectNormal.y * cos(na),
                        objectNormal.z);
  }`;
const BEND_POSITION = `#include <begin_vertex>
  {
    float bd = -transformed.y;
    vD = bd;
    if (abs(uBend) > 1e-5 && bd > uY0) {
      float L = 1.0 - uY0, u = bd - uY0, ds = u / ${BEND_STEPS}.0, px = 0.0, py = 0.0;
      for (int i = 0; i < ${BEND_STEPS}; i++) {
        float sm = (float(i) + 0.5) * ds;
        float a = uBend * pow(sm / L, ${BEND_POW.toFixed(2)});
        px += sin(a) * ds; py += cos(a) * ds;
      }
      transformed.x += px;
      transformed.y = -(uY0 + py);
    }
  }`;
// Same path in JS, for the root tip position.
function bendOffset(depth, theta) {
  const y0 = BEND_ZONE_START, L = 1 - y0;
  if (depth <= y0 || Math.abs(theta) < 1e-5) return [0, depth];
  const u = depth - y0, ds = u / BEND_STEPS;
  let px = 0, py = 0;
  for (let i = 0; i < BEND_STEPS; i++) { const a = theta * Math.pow(((i + 0.5) * ds) / L, BEND_POW); px += Math.sin(a) * ds; py += Math.cos(a) * ds; }
  return [px, y0 + py];
}
function injectBend(material, uniforms, withNormals) {
  material.onBeforeCompile = (sh) => {
    Object.assign(sh.uniforms, uniforms);
    let vs = BEND_GLSL_HEAD + sh.vertexShader.replace("#include <begin_vertex>", BEND_POSITION);
    if (withNormals) vs = vs.replace("#include <beginnormal_vertex>", BEND_NORMAL);
    sh.vertexShader = vs;
    if (withNormals) {
      // Homing pulse: waves of light run up from the tips while water flows.
      sh.fragmentShader = "uniform float uTime;\nuniform float uPulse;\nvarying float vD;\n" + sh.fragmentShader.replace("#include <emissivemap_fragment>",
        "#include <emissivemap_fragment>\n  float wave = pow(0.5 + 0.5 * sin(18.0 * vD + uTime * 5.0), 6.0);\n  totalEmissiveRadiance += vec3(0.25, 0.72, 1.0) * uPulse * (0.18 + 0.82 * wave) * smoothstep(0.15, 0.6, vD);");
    }
  };
  material.customProgramCacheKey = () => "phytoposeidon-bend-" + (withNormals ? "n" : "d");
}

function rootGeometry(gltf) {
  let geo = null;
  gltf.scene.updateMatrixWorld(true);
  gltf.scene.traverse((o) => {
    if (!o.isMesh || geo) return;
    geo = new THREE.BufferGeometry();
    geo.setAttribute("position", toFloatAttr(o.geometry.attributes.position));
    if (o.geometry.attributes.normal) geo.setAttribute("normal", toFloatAttr(o.geometry.attributes.normal));
    geo.setIndex(o.geometry.index);
    geo.applyMatrix4(o.matrixWorld);
  });
  // The file carries smoothed normals (baked offline); recompute only if missing.
  if (!geo.attributes.normal) geo.computeVertexNormals();
  else geo.attributes.normal.normalized = false;
  geo.computeBoundingSphere();
  return geo;
}

function makeRoots(geo, clipPlanes) {
  const uniforms = { uBend: { value: 0 }, uY0: { value: BEND_ZONE_START }, uTime: { value: 0 }, uPulse: { value: 0 } };
  const mat = new THREE.MeshStandardMaterial({ color: 0xeadcc0, roughness: 0.62, metalness: 0, clippingPlanes: clipPlanes, side: THREE.DoubleSide });
  injectBend(mat, uniforms, true);
  const mesh = new THREE.Mesh(geo, mat);
  // Roots cast shadows but don't receive them: thin geometry gives speckled shadow acne.
  mesh.castShadow = true; mesh.receiveShadow = false; mesh.frustumCulled = false;
  const depth = new THREE.MeshDepthMaterial({ depthPacking: THREE.RGBADepthPacking, clippingPlanes: clipPlanes });
  injectBend(depth, uniforms, false);
  mesh.customDepthMaterial = depth;
  const g = new THREE.Group();
  g.add(mesh);
  return { group: g, mesh, uniforms };
}

// Fine soil grain, multiplied over the moisture colours.
function soilTexture() {
  // Soft, low-contrast grain: fine speckle reads as noise once it is minified.
  const c = document.createElement("canvas");
  c.width = c.height = 512;
  const g = c.getContext("2d");
  g.fillStyle = "#f2f2f2"; g.fillRect(0, 0, 512, 512);
  let seed = 7;
  const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
  for (let i = 0; i < 2600; i++) {
    const v = 200 + Math.floor(rnd() * 55);
    g.fillStyle = `rgba(${v},${v - 6},${v - 14},${0.12 + rnd() * 0.22})`;
    const r = 2 + rnd() * 5;
    g.beginPath(); g.arc(rnd() * 512, rnd() * 512, r, 0, Math.PI * 2); g.fill();
  }
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(2, 2);
  t.anisotropy = 8;
  return t;
}

const SOIL_DRY = [0.86, 0.74, 0.55], SOIL_WET = [0.16, 0.11, 0.08], SOIL_EVEN = [0.40, 0.29, 0.20];
// Glass with a defined outline: nearly clear where the wall faces the
// viewer, bright where it turns away (Fresnel rim), so the silhouette reads.
function outlineGlass(color) {
  return new THREE.ShaderMaterial({
    uniforms: { uColor: { value: new THREE.Color(color) }, uBase: { value: 0.035 }, uRim: { value: 0.85 } },
    // Outline from the wall's silhouette only: the Fresnel term uses the
    // horizontal parts of the normal and view, so the floor and rim tops stay clear.
    vertexShader: `varying vec3 vN; varying vec3 vW;
      void main() {
        vec4 w = modelMatrix * vec4(position, 1.0);
        vW = w.xyz; vN = normalize(mat3(modelMatrix) * normal);
        gl_Position = projectionMatrix * viewMatrix * w;
      }`,
    fragmentShader: `uniform vec3 uColor; uniform float uBase; uniform float uRim; varying vec3 vN; varying vec3 vW;
      void main() {
        vec3 v = cameraPosition - vW;
        vec2 nh = vN.xz, vh = v.xz;
        float wall = clamp(length(nh), 0.0, 1.0);
        float f = 1.0 - abs(dot(normalize(nh + 1e-5), normalize(vh + 1e-5)));
        float a = uBase + uRim * pow(f, 4.0) * smoothstep(0.55, 0.95, wall);
        gl_FragColor = vec4(uColor, clamp(a, 0.0, 0.9));
        #include <colorspace_fragment>
      }`,
    transparent: true, depthWrite: false, side: THREE.DoubleSide,
  });
}

// Transparent round pot: see-through soil, roots inside.
function makePot(potGeo, soilTex, soilClip = null) {
  const g = new THREE.Group();
  const glassMat = outlineGlass(0xd8f1ff);
  const glass = new THREE.Mesh(potGeo, glassMat);
  glass.renderOrder = 3;
  g.add(glass);
  // Rim and base creases as crisp lines.
  // Keep the rim and outer base rings; drop the broken ring on the inside floor.
  const edgeGeo = new THREE.EdgesGeometry(potGeo, 35);
  {
    const p = edgeGeo.attributes.position, keep = [];
    for (let i = 0; i < p.count; i += 2) {
      const y0 = p.getY(i), y1 = p.getY(i + 1);
      if (Math.max(y0, y1) > 0.3 || Math.max(y0, y1) < 0.03) keep.push(p.getX(i), y0, p.getZ(i), p.getX(i + 1), y1, p.getZ(i + 1));
    }
    edgeGeo.setAttribute("position", new THREE.Float32BufferAttribute(keep, 3));
  }
  const edges = new THREE.LineSegments(edgeGeo,
    new THREE.LineBasicMaterial({ color: 0xe4f4ff, transparent: true, opacity: 0.7, depthWrite: false }));
  edges.renderOrder = 4;
  g.add(edges);
  // Soil: a cone inside the wall; translucent so the roots show.
  const SH = SOIL_Y - POT_FLOOR;
  const soilGeo = new THREE.CylinderGeometry(innerR(SOIL_Y) - 0.01, innerR(POT_FLOOR) - 0.01, SH, 96, 1);
  soilGeo.setAttribute("color", new THREE.BufferAttribute(new Float32Array(soilGeo.attributes.position.count * 3), 3));
  const side = new THREE.MeshStandardMaterial({ vertexColors: true, map: soilTex, color: 0xb8b8b8, transparent: true, opacity: 0.3, roughness: 1, depthWrite: false, side: THREE.DoubleSide, clippingPlanes: soilClip });
  const top = new THREE.MeshStandardMaterial({ vertexColors: true, map: soilTex, color: 0xb8b8b8, transparent: true, opacity: 0.93, roughness: 1, depthWrite: false, side: THREE.DoubleSide, clippingPlanes: soilClip });
  side.side = THREE.FrontSide; side.opacity = 0.22;
  // The far half of the soil body is denser, so the soil reads as a round
  // fill behind the roots; the near half stays light so the roots show.
  const back = new THREE.MeshStandardMaterial({ vertexColors: true, map: soilTex, color: 0xb8b8b8, transparent: true, opacity: 0.78, roughness: 1, depthWrite: false, side: THREE.BackSide, clippingPlanes: soilClip });
  const soilBack = new THREE.Mesh(soilGeo, back);
  soilBack.position.y = POT_FLOOR + SH / 2;
  soilBack.renderOrder = 1.5;
  soilBack.receiveShadow = true;
  g.add(soilBack);
  const soil = new THREE.Mesh(soilGeo, [side, top, side]);
  soil.position.y = POT_FLOOR + SH / 2;
  soil.renderOrder = 2;
  soil.receiveShadow = true;
  g.add(soil);
  const SW = innerR(SOIL_Y) * 2;
  function setTheme(light) {
    glassMat.uniforms.uColor.value.set(light ? 0x2c4a60 : 0xd8f1ff);
    edges.material.color.set(light ? 0x2c4a60 : 0xe4f4ff);
  }

  // Moisture colours: even (A) or a gradient toward the water side (B).
  const _c = new THREE.Color();
  function setMoisture(sideSign, dpsi) {
    const pos = soilGeo.attributes.position, col = soilGeo.attributes.color;
    const G = gradGain(dpsi);
    for (let i = 0; i < pos.count; i++) {
      let c;
      if (sideSign === 0) c = SOIL_EVEN;
      else {
        const f = (pos.getX(i) / SW + 0.5);                 // 0 left … 1 right
        const wet = sideSign > 0 ? f : 1 - f;               // 1 at the water
        const w = Math.pow(wet, 1.4) * G + (1 - G) * 0.55;  // flat when there is no gradient
        c = [0, 1, 2].map((k) => SOIL_DRY[k] + (SOIL_WET[k] - SOIL_DRY[k]) * w);
      }
      _c.setRGB(c[0], c[1], c[2], THREE.SRGBColorSpace);
      col.setXYZ(i, _c.r, _c.g, _c.b);
    }
    col.needsUpdate = true;
  }
  return { group: g, setMoisture, setTheme, soilTop: SOIL_Y };
}

/* Water reservoir for pot B: a crescent of water behind a porous membrane
   on the chosen side, with the animated water surface on top. Water moves
   from the reservoir into the soil, down the water-potential gradient
   (wet → dry); the drifting particles show that flow. */
const RES_CHORD = 0.5;          // membrane distance from the pot axis
const RES_TOP = SOIL_Y - 0.05;  // reservoir water level
function membraneTexture() {
  const c = document.createElement("canvas");
  c.width = c.height = 128;
  const g = c.getContext("2d");
  g.fillStyle = "rgba(255,255,255,0.8)"; g.fillRect(0, 0, 128, 128);
  g.globalCompositeOperation = "destination-out";
  for (let y = 4; y < 128; y += 8) for (let x = (y / 8) % 2 ? 8 : 4; x < 128; x += 8) { g.beginPath(); g.arc(x, y, 2.2, 0, Math.PI * 2); g.fill(); }
  const t = new THREE.CanvasTexture(c);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(5, 5);
  return t;
}
// A water droplet: round base, pointed top, with a highlight.
function dropletPath(g, x, y, r) {
  g.beginPath();
  g.moveTo(x, y - r * 1.55);
  g.bezierCurveTo(x + r * 0.35, y - r * 0.95, x + r, y - r * 0.35, x + r, y + r * 0.15);
  g.arc(x, y + r * 0.15, r, 0, Math.PI);
  g.bezierCurveTo(x - r, y - r * 0.35, x - r * 0.35, y - r * 0.95, x, y - r * 1.55);
  g.closePath();
}
function drawDroplet(g, x, y, r, fill = "rgba(110,195,255,0.92)") {
  dropletPath(g, x, y, r);
  g.fillStyle = fill; g.fill();
  g.strokeStyle = "rgba(225,245,255,0.85)"; g.lineWidth = Math.max(1, r * 0.16); g.stroke();
  g.fillStyle = "rgba(255,255,255,0.85)";
  g.beginPath(); g.ellipse(x - r * 0.38, y + r * 0.02, r * 0.18, r * 0.32, -0.35, 0, Math.PI * 2); g.fill();
}
function dotTexture() {
  const c = document.createElement("canvas");
  c.width = c.height = 64;
  const g = c.getContext("2d");
  drawDroplet(g, 32, 36, 17, "rgba(255,255,255,0.95)");
  return new THREE.CanvasTexture(c);
}
function ringTexture() {
  const c = document.createElement("canvas");
  c.width = c.height = 128;
  const g = c.getContext("2d");
  const grd = g.createRadialGradient(64, 64, 40, 64, 64, 62);
  grd.addColorStop(0, "rgba(120,210,255,0)"); grd.addColorStop(0.55, "rgba(150,225,255,0.95)"); grd.addColorStop(1, "rgba(120,210,255,0)");
  g.fillStyle = grd; g.beginPath(); g.arc(64, 64, 63, 0, Math.PI * 2); g.fill();
  return new THREE.CanvasTexture(c);
}
function makeReservoir(waterGltf, envTex, potX, coneClip) {
  const keepRes = new THREE.Plane(new THREE.Vector3(1, 0, 0), 0);    // reservoir side of the membrane
  const keepSoil = new THREE.Plane(new THREE.Vector3(-1, 0, 0), 1e3); // soil side (disabled when no gradient)
  const g = new THREE.Group();
  g.position.set(potX, 0, 0);
  const H = RES_TOP - POT_FLOOR;
  const waterMat = new THREE.MeshPhysicalMaterial({
    color: 0x2f9ee0, transparent: true, opacity: 0.5, roughness: 0.05, metalness: 0, envMap: envTex, envMapIntensity: 1.1,
    emissive: 0x0a3a60, emissiveIntensity: 0.45, depthWrite: false, side: THREE.DoubleSide, clippingPlanes: [keepRes],
  });
  const body = new THREE.Mesh(new THREE.CylinderGeometry(innerR(RES_TOP) - 0.012, innerR(POT_FLOOR) - 0.012, H, 96, 1), waterMat);
  body.position.y = POT_FLOOR + H / 2;
  body.renderOrder = 2.6;
  g.add(body);
  // Animated surface ("Water Animation", morph-target ripples), trimmed to the crescent.
  const surf = waterGltf.scene;
  const surfMat = new THREE.MeshPhysicalMaterial({
    color: 0x6cc6f5, transparent: true, opacity: 0.88, roughness: 0.03, metalness: 0, envMap: envTex, envMapIntensity: 1.6,
    clearcoat: 1, clearcoatRoughness: 0.03, side: THREE.DoubleSide, depthWrite: false, clippingPlanes: [keepRes, ...coneClip],
  });
  surf.traverse((o) => { if (o.isMesh) { o.material = surfMat; o.renderOrder = 2.7; o.frustumCulled = false; } });
  const surfHolder = new THREE.Group();
  surfHolder.add(surf);
  surfHolder.scale.setScalar(1.8 / 100);
  surfHolder.position.y = RES_TOP;
  g.add(surfHolder);
  const mixer = new THREE.AnimationMixer(surf);
  if (waterGltf.animations[0]) mixer.clipAction(waterGltf.animations[0]).play();
  // Porous membrane between the water and the soil (a trapezoid in the y–z plane).
  const zt = Math.sqrt(innerR(RES_TOP) ** 2 - RES_CHORD ** 2), zb = Math.sqrt(Math.max(0.01, innerR(POT_FLOOR) ** 2 - RES_CHORD ** 2));
  const shape = new THREE.Shape([new THREE.Vector2(-zb, POT_FLOOR), new THREE.Vector2(zb, POT_FLOOR), new THREE.Vector2(zt, SOIL_Y), new THREE.Vector2(-zt, SOIL_Y)]);
  const memGeo = new THREE.ShapeGeometry(shape);
  const uv = memGeo.attributes.uv, mp = memGeo.attributes.position;
  for (let i = 0; i < uv.count; i++) uv.setXY(i, mp.getX(i) * 0.8, mp.getY(i) * 0.8);
  const membrane = new THREE.Mesh(memGeo, new THREE.MeshStandardMaterial({ color: 0xd9e4ea, map: membraneTexture(), transparent: true, roughness: 0.6, side: THREE.DoubleSide, depthWrite: false }));
  membrane.renderOrder = 2.8;
  const memHolder = new THREE.Group();
  memHolder.add(membrane);
  g.add(memHolder);
  // Flow particles: water moving from the membrane into the drier soil.
  const N = 240;
  const pos = new Float32Array(N * 3), col = new Float32Array(N * 4);
  const pg = new THREE.BufferGeometry();
  pg.setAttribute("position", new THREE.BufferAttribute(pos, 3));
  pg.setAttribute("color", new THREE.BufferAttribute(col, 4));
  const pts = new THREE.Points(pg, new THREE.PointsMaterial({ size: 0.05, map: dotTexture(), vertexColors: true, transparent: true, depthWrite: false, clippingPlanes: coneClip }));
  pts.renderOrder = 2.9;
  pts.frustumCulled = false;
  g.add(pts);
  const P = Array.from({ length: N }, () => ({ d: Math.random() * 1.6, y: 0, z: 0, v: 0.7 + Math.random() * 0.6 }));
  const respawn = (p) => {
    p.d = 0;
    p.y = POT_FLOOR + 0.06 + Math.random() * (SOIL_Y - POT_FLOOR - 0.12);
    const zr = Math.sqrt(Math.max(0, innerR(p.y) ** 2 - RES_CHORD ** 2));
    p.z = (Math.random() * 2 - 1) * zr * 0.95;
    p.v = 0.7 + Math.random() * 0.6;
  };
  P.forEach((p) => { respawn(p); p.d = Math.random() * 1.6; });
  let side = 1, G = 0;
  function set(sideSign, gain) {
    side = sideSign; G = gain;
    const on = G > 0;
    g.visible = on;
    keepRes.normal.set(side, 0, 0); keepRes.constant = -RES_CHORD - side * potX;
    keepSoil.normal.set(-side, 0, 0); keepSoil.constant = on ? RES_CHORD + side * potX : 1e3;
    memHolder.position.x = side * RES_CHORD;
    memHolder.rotation.y = Math.PI / 2;
    waterMat.opacity = 0.35 + 0.25 * G;
  }
  function update(dt, light) {
    if (!g.visible) return;
    mixer.update(dt);
    // Flow speed and particle density scale with the gradient.
    const speed = 0.16 * (0.4 + G);
    const shown = Math.round(N * (0.25 + 0.75 * G));
    const c = light ? [0.05, 0.42, 0.75] : [0.55, 0.85, 1.0];
    for (let i = 0; i < N; i++) {
      const p = P[i];
      p.d += speed * p.v * dt;
      if (p.d > 1.55) respawn(p);
      const x = side * (RES_CHORD - p.d);
      const wob = Math.sin((p.d * 9) + i) * 0.012;
      pos[i * 3] = x; pos[i * 3 + 1] = p.y + wob; pos[i * 3 + 2] = p.z + wob;
      const a = i < shown ? Math.max(0, 1 - p.d / 1.55) ** 1.2 : 0;   // fewer water molecules toward the dry side
      col.set([c[0], c[1], c[2], a], i * 4);
    }
    pg.attributes.position.needsUpdate = true;
    pg.attributes.color.needsUpdate = true;
  }
  return { group: g, set, update, keepSoil };
}

function makeMagnifier(gltf) {
  const LENS_N = new THREE.Vector3(-0.46275, 0.37837, 0.80168).normalize();
  const HANDLE = new THREE.Vector3(0.47965, -0.65367, 0.58537).normalize();
  const LENS_C = new THREE.Vector3(-0.0062, 0.0085, -0.0076);
  const LENS_R = 1.0;
  const E = new THREE.Matrix4().makeBasis(HANDLE, new THREE.Vector3().crossVectors(LENS_N, HANDLE), LENS_N);
  const f1 = new THREE.Vector3(Math.SQRT1_2, -Math.SQRT1_2, 0), f3 = new THREE.Vector3(0, 0, 1);
  const F = new THREE.Matrix4().makeBasis(f1, new THREE.Vector3().crossVectors(f3, f1), f3);
  const R = F.multiply(E.transpose());
  const model = gltf.scene;
  const inner = new THREE.Group();
  inner.add(model);
  model.position.copy(LENS_C).multiplyScalar(-1);
  inner.applyMatrix4(R);
  inner.scale.multiplyScalar(1 / LENS_R);
  model.traverse((o) => { if (o.isMesh) { o.renderOrder = 10; o.material.depthTest = false; } });
  const holder = new THREE.Group();
  holder.add(inner);
  return holder;
}

/* ---------------- Root tip inside the lens ----------------
   The "Root Structure" model (a cut-away root tip), rendered into the lens.
   It bends with the model's tip angle: the curve sits in the elongation
   zone, as a true arc, so the outer (drier) flank is longer than the inner
   one. The meristem and cap stay straight and turn with the tip. Cortex on
   the drier flank of the elongation zone glows teal (MIZ1/SnRK2.2 activity).
   Behind it, soil runs from wet to dry and water molecules drift down the
   water-potential gradient; labels and amyloplasts are drawn on top. */
const LENS_PX = 768;
const RT_TOP = 0.040, RT_TIP = -0.148, RT_H = RT_TOP - RT_TIP;       // model units (tip down)
const RT_D0 = 0.12 * RT_H, RT_D1 = 0.78 * RT_H;                       // bending region, from the top
// Natural bend: the turning angle follows a smooth S-curve with depth,
// a(d) = θ·smoothstep(d0, d1, d), so curvature peaks mid-zone and fades at
// both ends. Each point sits at x along the local normal of the bent axis,
// so the outer flank is genuinely longer than the inner flank.
const RT_STEPS = 16;
const rtSmooth = (d) => { const t = Math.max(0, Math.min(1, (d - RT_D0) / (RT_D1 - RT_D0))); return t * t * (3 - 2 * t); };
function bendPoint(x, y, theta) {
  const d = Math.max(0, RT_TOP - y), ds = d / RT_STEPS;
  let cx = 0, cy = RT_TOP;
  for (let i = 0; i < RT_STEPS; i++) { const a = theta * rtSmooth((i + 0.5) * ds); cx += Math.sin(a) * ds; cy -= Math.cos(a) * ds; }
  const a = theta * rtSmooth(d);
  return [cx + x * Math.cos(a), cy + x * Math.sin(a)];
}
const TIP_GLSL_HEAD = "uniform float uTheta;\nuniform float uSide;\nuniform float uG;\nvarying float vEZ;\nvarying float vX;\n";
const RT_GLSL_FN = `float rtS(float d) { float t = clamp((d - ${RT_D0.toFixed(5)}) / ${(RT_D1 - RT_D0).toFixed(5)}, 0.0, 1.0); return t * t * (3.0 - 2.0 * t); }
`;
const TIP_BEND = `#include <begin_vertex>
  {
    float d = max(${RT_TOP.toFixed(4)} - transformed.y, 0.0);
    vX = transformed.x;
    float sd = rtS(d);
    vEZ = sd * (1.0 - sd) * 4.0;
    if (abs(uTheta) > 1e-4) {
      float ds = d / ${RT_STEPS}.0, cx = 0.0, cy = ${RT_TOP.toFixed(4)};
      for (int i = 0; i < ${RT_STEPS}; i++) { float a = uTheta * rtS((float(i) + 0.5) * ds); cx += sin(a) * ds; cy -= cos(a) * ds; }
      float a = uTheta * sd, x = transformed.x;
      transformed.x = cx + x * cos(a);
      transformed.y = cy + x * sin(a);
    }
  }`;
const TIP_NORMAL = `#include <beginnormal_vertex>
  {
    float a = uTheta * rtS(max(${RT_TOP.toFixed(4)} - position.y, 0.0));
    objectNormal = vec3(cos(a) * objectNormal.x - sin(a) * objectNormal.y, sin(a) * objectNormal.x + cos(a) * objectNormal.y, objectNormal.z);
  }`;
function makeLens(tipGltf) {
  const scene = new THREE.Scene();
  const VIEW = 0.145;
  const camera = new THREE.OrthographicCamera(-VIEW, VIEW, VIEW, -VIEW, 0.01, 10);
  camera.position.set(0, 0, 2);
  scene.add(new THREE.HemisphereLight(0xf2f7ff, 0x2a2016, 1.5));
  const key = new THREE.DirectionalLight(0xfff3dd, 2.2);
  key.position.set(1.2, 1.6, 2.5);
  scene.add(key);
  // Mipmapped, so the lens stays smooth when it is small on screen.
  const rt = new THREE.WebGLRenderTarget(LENS_PX, LENS_PX, { samples: 4, generateMipmaps: true, minFilter: THREE.LinearMipmapLinearFilter });
  rt.texture.anisotropy = 4;
  rt.texture.colorSpace = THREE.SRGBColorSpace;
  // Background (soil, water molecules) and overlay (labels, amyloplasts) canvases.
  const mkCanvas = () => { const c = document.createElement("canvas"); c.width = c.height = LENS_PX; const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace; return { c, g: c.getContext("2d"), t }; };
  const bg = mkCanvas(), ov = mkCanvas();
  const plane = (tex, z, order) => {
    // The background is opaque (drawn first); the overlay is transparent (drawn last).
    const m = new THREE.Mesh(new THREE.PlaneGeometry(VIEW * 2, VIEW * 2), new THREE.MeshBasicMaterial({ map: tex, transparent: order > 0, depthTest: false, depthWrite: false, toneMapped: false }));
    m.position.z = z; m.renderOrder = order; scene.add(m); return m;
  };
  plane(bg.t, -1, -1);
  plane(ov.t, 1, 10);
  // The model, without its baked text labels.
  const U = { uTheta: { value: 0 }, uSide: { value: 1 }, uG: { value: 0 } };
  const model = new THREE.Group();
  tipGltf.scene.updateMatrixWorld(true);
  tipGltf.scene.traverse((o) => {
    if (!o.isMesh) return;
    const g = new THREE.BufferGeometry();
    ["position", "normal", "uv"].forEach((k) => { if (o.geometry.attributes[k]) g.setAttribute(k, toFloatAttr(o.geometry.attributes[k])); });
    if (o.geometry.index) g.setIndex(o.geometry.index);
    g.applyMatrix4(o.matrixWorld);
    if (!g.attributes.normal) g.computeVertexNormals();
    const mat = o.material.clone();
    mat.side = THREE.DoubleSide;
    mat.onBeforeCompile = (sh) => {
      Object.assign(sh.uniforms, U);
      sh.vertexShader = TIP_GLSL_HEAD + RT_GLSL_FN + sh.vertexShader.replace("#include <begin_vertex>", TIP_BEND).replace("#include <beginnormal_vertex>", TIP_NORMAL);
      sh.fragmentShader = "uniform float uSide;\nuniform float uG;\nvarying float vEZ;\nvarying float vX;\n" + sh.fragmentShader.replace("#include <map_fragment>",
        // Recolour the model's green tissues to root cream, keeping their light/dark banding.
        "#include <map_fragment>\n  { float lum = dot(diffuseColor.rgb, vec3(0.3, 0.59, 0.11)); diffuseColor.rgb = mix(vec3(0.55, 0.42, 0.28), vec3(1.0, 0.95, 0.84), clamp(0.3 + lum * 1.5, 0.0, 1.0)); }").replace("#include <emissivemap_fragment>",
        "#include <emissivemap_fragment>\n  totalEmissiveRadiance += vec3(0.05, 0.55, 0.5) * uG * vEZ * (1.0 - smoothstep(-0.004, 0.004, uSide * vX)) * 0.9;");
    };
    mat.customProgramCacheKey = () => "phytoposeidon-tip";
    model.add(new THREE.Mesh(g, mat));
  });
  scene.add(model);
  const water = Array.from({ length: 150 }, (_, i) => ({ x: Math.random() * 2 - 1, y: Math.random() * 2 - 1, v: 0.5 + Math.random(), p: i }));
  let last = 0;
  // Model units → lens canvas pixels (after the model is centred).
  let off = [0, 0];
  const toPx = (p) => [((p[0] + off[0]) / VIEW * 0.5 + 0.5) * LENS_PX, (0.5 - (p[1] + off[1]) / VIEW * 0.5) * LENS_PX];
  function drawBackground(side, G, now, dt) {
    const g = bg.g;
    const grd = g.createLinearGradient(0, 0, LENS_PX, 0);
    const wetC = "#33291f", dryC = "#a88a63", evenC = "#5e4838";
    if (G > 0) { grd.addColorStop(0, side > 0 ? dryC : wetC); grd.addColorStop(1, side > 0 ? wetC : dryC); }
    else { grd.addColorStop(0, evenC); grd.addColorStop(1, evenC); }
    g.fillStyle = grd; g.fillRect(0, 0, LENS_PX, LENS_PX);
    g.fillStyle = "rgba(255,255,255,0.05)";
    for (let i = 0; i < 160; i++) { const r = ((i * 7919) % 97) / 97; g.beginPath(); g.arc((i * 1543) % LENS_PX, (i * 2711) % LENS_PX, 3 + r * 7, 0, Math.PI * 2); g.fill(); }
    for (const w of water) {
      if (G > 0) { w.x += -side * w.v * 0.12 * (0.4 + G) * dt; if (side * w.x < -1.05) { w.x = side * (0.7 + Math.random() * 0.35); w.y = Math.random() * 2 - 1; } }
      else { w.x += Math.sin(now / 700 + w.p) * 0.02 * dt; w.y += Math.cos(now / 900 + w.p) * 0.02 * dt; }
      const wet = G > 0 ? Math.max(0, Math.min(1, (side * w.x + 1) / 2)) : 0.55;
      if (G > 0 && ((w.p * 0.618) % 1) > 0.3 + 0.7 * wet) continue;   // fewer molecules on the dry side
      const x = LENS_PX / 2 + w.x * LENS_PX * 0.5, y = LENS_PX / 2 + w.y * LENS_PX * 0.5;
      drawDroplet(g, x, y, 8.5);
    }
    bg.t.needsUpdate = true;
  }
  function drawOverlay(theta, side, G, amy, light) {
    const g = ov.g;
    g.clearRect(0, 0, LENS_PX, LENS_PX);
    g.strokeStyle = light ? "rgba(20,40,60,0.5)" : "rgba(255,255,255,0.35)"; g.lineWidth = 6;
    g.beginPath(); g.arc(LENS_PX / 2, LENS_PX / 2, LENS_PX / 2 - 3, 0, Math.PI * 2); g.stroke();
    ov.t.needsUpdate = true;
  }
  function drawEmpty(light) {
    bg.g.fillStyle = light ? "#e9eef3" : "#11161d"; bg.g.fillRect(0, 0, LENS_PX, LENS_PX); bg.t.needsUpdate = true;
    const g = ov.g;
    g.clearRect(0, 0, LENS_PX, LENS_PX);
    g.fillStyle = light ? "#2b3a48" : "#dbe6ee"; g.font = "600 34px 'IBM Plex Sans', Roboto, sans-serif"; g.textAlign = "center"; g.textBaseline = "middle";
    g.fillText("Drag the lens", LENS_PX / 2, LENS_PX / 2 - 24); g.fillText("over a root tip", LENS_PX / 2, LENS_PX / 2 + 24);
    ov.t.needsUpdate = true;
  }
  function render(renderer, st, now) {
    const dt = last ? Math.min(0.1, (now - last) / 1000) : 0.016; last = now;
    if (!st) { drawEmpty(document.documentElement.dataset.theme === "light"); model.visible = false; }
    else {
      const th = toRad(st.theta);
      U.uTheta.value = th; U.uSide.value = st.side; U.uG.value = st.G;
      // Centre the bent model between its top and its tip.
      const top = [0, RT_TOP], tip = bendPoint(0, RT_TIP, th);
      // Centre on the middle of the bent root (fixed regardless of the water side).
      const mid = bendPoint(0, RT_TOP - RT_H * 0.5, th);
      off = [-(mid[0] * 0.5 + (top[0] + tip[0]) * 0.25), -(mid[1] * 0.5 + (top[1] + tip[1]) * 0.25)];
      model.position.set(off[0], off[1], 0);
      model.visible = true;
      drawBackground(st.side, st.G, now, dt);
      drawOverlay(th, st.side, st.G, st.amy, st.light);
    }
    renderer.setRenderTarget(rt);
    renderer.render(scene, camera);
    renderer.setRenderTarget(null);
  }
  return { tex: rt.texture, render };
}

async function init() {
  const [plantGltf, rootsA, rootsB, potGltf, waterGltf, magGltf, tipGltf] = await Promise.all([
    track(parseGLB("glbPlant"), "Plant"), track(parseGLB("glbRoots"), "Roots A"), track(parseGLB("glbRoots"), "Roots B"), track(parseGLB("glbPot"), "Pot"), track(parseGLB("glbWater"), "Water"), track(parseGLB("glbMag"), "Magnifier"), track(parseGLB("glbTip"), "Root tip"),
  ]);

  const stage = $("stage");
  const canvas = $("glCanvas");
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
  renderer.setClearColor(0x000000, 0);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.05;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.shadowMap.autoUpdate = false;   // refreshed only when the roots change
  renderer.localClippingEnabled = true;    // roots stay inside their pots

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x15171c);
  scene.fog = new THREE.Fog(0x15171c, 8, 16);
  scene.add(new THREE.HemisphereLight(0xd6deff, 0x2a1f14, 1.35));
  const fill = new THREE.DirectionalLight(0xbfcfff, 0.55);
  fill.position.set(0, 3, 8);
  scene.add(fill);
  // Overhead lab light (the only shadow caster).
  const key = new THREE.DirectionalLight(0xfff4e6, 1.9);
  key.position.set(-2.2, 5.5, 3.2);
  key.castShadow = true;
  key.shadow.mapSize.set(2048, 1024);
  key.shadow.bias = -0.0004;
  key.shadow.normalBias = 0.02;
  Object.assign(key.shadow.camera, { left: -3, right: 3, top: 3, bottom: -1.5, near: 1, far: 14 });
  scene.add(key);

  // Floor and partition.
  const floor = new THREE.Mesh(new THREE.CircleGeometry(30, 96), new THREE.MeshStandardMaterial({ color: 0x22252c, roughness: 0.95 }));
  const hemi = scene.children.find((o) => o.isHemisphereLight);
  floor.rotation.x = -Math.PI / 2;
  floor.receiveShadow = true;
  scene.add(floor);
  const pmrem = new THREE.PMREMGenerator(renderer);
  const envTex = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
  scene.add(makeDivider(envTex));

  // Pots, roots and shoots.
  const soilTex = soilTexture();
  const potGeo = (() => {
    let geo = null;
    potGltf.scene.updateMatrixWorld(true);
    potGltf.scene.traverse((o) => {
      if (!o.isMesh || geo) return;
      geo = new THREE.BufferGeometry();
      geo.setAttribute("position", toFloatAttr(o.geometry.attributes.position));
      geo.setIndex(o.geometry.index);
      geo.applyMatrix4(o.matrixWorld);
    });
    geo.scale(POT_K, POT_K, POT_K);
    // Keep only the outer skin (outer wall, rim face and base), so the glass
    // shows a single clean outline instead of inner and outer walls.
    const src = geo.toNonIndexed(), p = src.attributes.position, keep = [];
    const a = new THREE.Vector3(), b = new THREE.Vector3(), c = new THREE.Vector3(), n = new THREE.Vector3(), m = new THREE.Vector3();
    for (let i = 0; i < p.count; i += 3) {
      a.fromBufferAttribute(p, i); b.fromBufferAttribute(p, i + 1); c.fromBufferAttribute(p, i + 2);
      n.subVectors(c, b).cross(m.subVectors(a, b)).normalize();
      m.addVectors(a, b).add(c).divideScalar(3);
      const rad = Math.hypot(m.x, m.z) || 1;
      const outward = (n.x * m.x + n.z * m.z) / rad;
      const outerR = m.y > 0.85 * POT_H ? POT_RIM * 0.9 : innerR(Math.max(m.y, POT_FLOOR)) + 0.012;
      if ((outward > 0.25 && rad > outerR) || (n.y < -0.6 && m.y < 0.05)) keep.push(a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
    }
    const out = new THREE.BufferGeometry();
    out.setAttribute("position", new THREE.Float32BufferAttribute(keep, 3));
    out.computeVertexNormals();
    return out;
  })();
  const shootParts = bakeShoot(plantGltf.scene);
  const rootGeo = { left: rootGeometry(rootsA), right: rootGeometry(rootsB) };
  const plants = {};
  for (const side of ["left", "right"]) {
    const x = side === "left" ? -POT_X : POT_X;
    const chord = side === "right" ? new THREE.Plane(new THREE.Vector3(-1, 0, 0), 1e3) : null;  // membrane cut, pot B only
    const pot = makePot(potGeo, soilTex, chord ? [chord] : null);
    pot.group.position.set(x, 0, 0);
    // Roots stay inside the pot: a 16-sided cone just inside the wall, plus soil surface and floor.
    const yA = POT_FLOOR, yB = SOIL_Y, rA = innerR(yA) - 0.02, rB = innerR(yB) - 0.02;
    const slope = (rB - rA) / (yB - yA), r0 = rA - slope * yA;      // r(y) = r0 + slope·y
    const clip = [];
    for (let i = 0; i < 16; i++) {
      const f = (i / 16) * Math.PI * 2, c = Math.cos(f), sn = Math.sin(f);
      const n = new THREE.Vector3(-c, slope, -sn), len = n.length();
      clip.push(new THREE.Plane(n.divideScalar(len), (r0 + x * c) / len));
    }
    clip.push(new THREE.Plane(new THREE.Vector3(0, -1, 0), SOIL_Y + 0.005), new THREE.Plane(new THREE.Vector3(0, 1, 0), -yA));
    const roots = makeRoots(rootGeo[side], chord ? [...clip, chord] : clip);
    roots.group.position.set(x, SOIL_Y, 0);
    const shoot = makeShoot(shootParts);
    shoot.position.set(x, SOIL_Y - 0.01, 0);
    scene.add(pot.group, roots.group, shoot);
    plants[side] = { pot, roots, shoot, x, clip, chord };
  }
  const reservoir = makeReservoir(waterGltf, envTex, POT_X, plants.right.clip);
  scene.add(reservoir.group);
  // Homing beacon at root tip B: rings ping out while water is flowing.
  const ringTex = ringTexture();
  const pings = [0, 1].map(() => {
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: ringTex, color: 0x8fd8ff, transparent: true, depthTest: false, depthWrite: false, blending: THREE.AdditiveBlending }));
    sp.renderOrder = 6; sp.visible = false; scene.add(sp); return sp;
  });

  // Magnifier: a root-tip section drawn into the lens.
  const mag = makeMagnifier(magGltf);
  const lens = makeLens(tipGltf);
  const lensDisc = new THREE.Mesh(new THREE.CircleGeometry(0.94, 72), new THREE.MeshBasicMaterial({ map: lens.tex, depthTest: false, toneMapped: false, transparent: true }));
  lensDisc.renderOrder = 11;   // over the glass, so the view isn't tinted
  mag.add(lensDisc);
  mag.visible = false;
  scene.add(mag);
  const lensDefault = () => (stage.clientWidth < 560 ? 0.55 : 0.42);
  const MAG_HOME = new THREE.Vector3(POT_X - 0.62, 0.72, 0.95);
  const magState = { on: false, pos: MAG_HOME.clone(), radius: 0.42, dragging: false, locked: null };
  plants.left.pot.setMoisture(0, 0);
  plants.right.pot.setMoisture(sim.side, sim.dpsi);

  // Camera.
  const camera = new THREE.PerspectiveCamera(32, 1, 0.1, 50);
  const TARGET = new THREE.Vector3(0, 1.18, 0);
  const controls = new OrbitControls(camera, canvas);
  controls.target.copy(TARGET);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.enablePan = false;
  controls.minPolarAngle = toRad(52);
  controls.maxPolarAngle = toRad(92);
  controls.minAzimuthAngle = toRad(-55);
  controls.maxAzimuthAngle = toRad(55);

  // Views: whole bench, or one room. Distances fit the current aspect ratio.
  let view = null;           // null (whole bench) | "left" | "right"
  let camAnim = null;
  const VIEW_DIR = new THREE.Vector3(0, 0.16, 1).normalize();
  function viewPose(which) {
    const vHalf = toRad(camera.fov / 2);
    const hHalf = Math.atan(Math.tan(vHalf) * camera.aspect);
    if (!which) {
      const d = Math.max(2.2 / Math.tan(hHalf), 1.32 / Math.tan(vHalf)) + 0.6;
      return { target: TARGET.clone(), pos: TARGET.clone().addScaledVector(VIEW_DIR, d), dist: d };
    }
    const tgt = new THREE.Vector3(which === "left" ? -POT_X : POT_X, 1.18, 0);
    const d = Math.max(1.0 / Math.tan(hHalf), 1.3 / Math.tan(vHalf)) + 0.6;
    return { target: tgt, pos: tgt.clone().addScaledVector(VIEW_DIR, d), dist: d };
  }
  function goTo(which, animate = true) {
    const pose = viewPose(which);
    const overview = viewPose(null);
    controls.minDistance = Math.min(pose.dist, overview.dist) * 0.55;
    controls.maxDistance = overview.dist * 1.35;
    if (!animate) {
      camera.position.copy(pose.pos);
      controls.target.copy(pose.target);
      controls.update();
      camAnim = null;
      return;
    }
    camAnim = { t0: performance.now(), dur: 950, p0: camera.position.clone(), t0v: controls.target.clone(), p1: pose.pos, t1v: pose.target };
  }
  const easeOutBack = (x) => { const c1 = 1.3, c3 = c1 + 1; return 1 + c3 * Math.pow(x - 1, 3) + c1 * Math.pow(x - 1, 2); };
  function stepCamAnim(now) {
    if (!camAnim) return;
    const k = Math.min(1, (now - camAnim.t0) / camAnim.dur);
    const e = easeOutBack(k);
    camera.position.lerpVectors(camAnim.p0, camAnim.p1, e);
    controls.target.lerpVectors(camAnim.t0v, camAnim.t1v, e);
    if (k >= 1) camAnim = null;
  }

  function resize() {
    const w = stage.clientWidth, h = stage.clientHeight;
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    if (!userMoved) goTo(view, false);
  }
  let userMoved = false;
  controls.addEventListener("start", () => { userMoved = true; camAnim = null; });
  const ro = new ResizeObserver(resize);
  ro.observe(stage);
  resize();

  /* ---------- Room switch ---------- */
  const roomSwitch = $("roomSwitch"), roomThumb = $("roomThumb");
  function setRoom(which, animate = true) {
    const prev = view;
    view = which;
    roomSwitch.dataset.room = which || "";
    roomSwitch.querySelectorAll(".room-opt").forEach((b) => b.setAttribute("aria-checked", String(b.dataset.room === which)));
    if (which && prev !== which) buzz(12);
    if (which && prev !== which && roomThumb.animate) {
      roomThumb.animate(
        [{ scale: "1 1" }, { scale: "1.14 0.86", offset: 0.35 }, { scale: "0.95 1.05", offset: 0.7 }, { scale: "1 1" }],
        { duration: 620, easing: "ease-out" }
      );
    }
    userMoved = false;
    goTo(which, animate);
    // The water source belongs to pot B: in room A it dims and says so.
    const wp = document.querySelector(".water-panel");
    wp.classList.toggle("for-b", which === "left");
    $("potHint").hidden = which !== "left";
    if (which && typeof magState !== "undefined") {
      magState.pos.set((which === "left" ? -POT_X : POT_X) - 0.62, 0.72, 0.95);
      magState.locked = null;
    }
  }
  roomSwitch.querySelectorAll(".room-opt").forEach((b) => b.addEventListener("click", () => setRoom(view === b.dataset.room ? null : b.dataset.room)));

  /* ---------- Navigation pad + arrow keys ---------- */
  let navHold = null;
  const sph = new THREE.Spherical();
  function orbitStep(dir, amount = 1) {
    camAnim = null;
    userMoved = true;
    const off = camera.position.clone().sub(controls.target);
    sph.setFromVector3(off);
    if (dir === "left") sph.theta -= toRad(1.6) * amount;
    if (dir === "right") sph.theta += toRad(1.6) * amount;
    if (dir === "up") sph.phi -= toRad(1.1) * amount;
    if (dir === "down") sph.phi += toRad(1.1) * amount;
    sph.theta = Math.max(controls.minAzimuthAngle, Math.min(controls.maxAzimuthAngle, sph.theta));
    sph.phi = Math.max(controls.minPolarAngle, Math.min(controls.maxPolarAngle, sph.phi));
    camera.position.copy(controls.target).add(off.setFromSpherical(sph));
  }
  $("navPad").querySelectorAll(".nav").forEach((b) => {
    const stop = () => { navHold = null; b.classList.remove("held"); };
    b.addEventListener("pointerdown", (e) => { e.preventDefault(); navHold = b.dataset.dir; b.classList.add("held"); b.setPointerCapture(e.pointerId); });
    b.addEventListener("pointerup", stop);
    b.addEventListener("pointercancel", stop);
    b.addEventListener("lostpointercapture", stop);
    b.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); orbitStep(b.dataset.dir, 4); } });
  });
  stage.addEventListener("keydown", (e) => {
    const map = { ArrowLeft: "left", ArrowRight: "right", ArrowUp: "up", ArrowDown: "down" };
    if (map[e.key] && (e.target === stage || e.target === canvas)) { e.preventDefault(); orbitStep(map[e.key], 3); }
  });

  /* ---------- Reset layout ---------- */
  $("layoutBtn").addEventListener("click", () => {
    magState.pos.copy(MAG_HOME);
    magState.radius = lensDefault();
    magState.locked = null;
    setRoom(null);
    const lb = $("layoutBtn");
    lb.classList.add("pulse");
    setTimeout(() => lb.classList.remove("pulse"), 650);
  });

  /* ---------- Screen projection helpers ---------- */
  const tmp = new THREE.Vector3();
  function toScreen(v) {
    tmp.copy(v).project(camera);
    return { x: (tmp.x * 0.5 + 0.5) * stage.clientWidth, y: (-tmp.y * 0.5 + 0.5) * stage.clientHeight, behind: tmp.z > 1 };
  }
  // Screen-space tags stay inside the stage and step aside for the nav pad.
  const navPad = $("navPad");
  let layout = null;
  const tagSize = new Map();
  function measureLayout() {
    layout = { W: stage.clientWidth, H: stage.clientHeight, nx: navPad.offsetLeft, ny: navPad.offsetTop, nw: navPad.offsetWidth, nh: navPad.offsetHeight };
    tagSize.clear();
  }
  new ResizeObserver(() => { layout = null; }).observe(stage);
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => { layout = null; });
  function placeTag(el, v, dy = 0) {
    if (!layout) measureLayout();
    if (!tagSize.has(el)) tagSize.set(el, { w: el.offsetWidth, h: el.offsetHeight });
    const p = toScreen(v);
    const { w, h } = tagSize.get(el);
    const { W, H } = layout;
    let x = Math.max(w / 2 + 4, Math.min(W - w / 2 - 4, p.x));
    let y = Math.max(h / 2 + 4, Math.min(H - h / 2 - 4, p.y + dy));
    const { nx, ny, nw, nh } = layout;
    const overlaps = x + w / 2 > nx - 4 && x - w / 2 < nx + nw + 4 && y + h / 2 > ny - 4 && y - h / 2 < ny + nh + 4;
    if (overlaps) x = nx - w / 2 - 6;
    el.style.transform = `translate(${x.toFixed(1)}px, ${y.toFixed(1)}px) translate(-50%, -50%)`;
  }

  /* ---------- Root tip in world space (matches the shader) ---------- */
  function tipWorld(side, bendDeg) {
    const P = plants[side];
    const s = P.roots.group.scale.x;
    const [dx, dy] = bendOffset(1, toRad(bendDeg));
    return new THREE.Vector3(P.x + dx * s, SOIL_Y - dy * s, 0);
  }

  /* ---------- Magnifier dragging ---------- */
  const raycaster = new THREE.Raycaster();
  const dragPlane = new THREE.Plane();
  const hit = $("magHit");
  function pointerToWorld(e) {
    const r = stage.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    raycaster.setFromCamera(ndc, camera);
    const out = new THREE.Vector3();
    return raycaster.ray.intersectPlane(dragPlane, out) ? out : null;
  }
  hit.addEventListener("pointerdown", (e) => {
    e.preventDefault();
    magState.dragging = true;
    hit.classList.add("dragging");
    hit.setPointerCapture(e.pointerId);
    dragPlane.setFromNormalAndCoplanarPoint(camera.getWorldDirection(new THREE.Vector3()).negate(), magState.pos);
    magState.grab = magState.pos.clone().sub(pointerToWorld(e) || magState.pos);
  });
  hit.addEventListener("pointermove", (e) => {
    if (!magState.dragging) return;
    const p = pointerToWorld(e);
    if (!p) return;
    p.add(magState.grab);
    p.x = Math.max(-2.9, Math.min(2.9, p.x));
    p.y = Math.max(0.2, Math.min(2.6, p.y));
    p.z = Math.max(-1, Math.min(2, p.z));
    magState.pos.copy(p);
  });
  const endDrag = (e) => {
    magState.dragging = false;
    hit.classList.remove("dragging");
    if (e && hit.hasPointerCapture(e.pointerId)) hit.releasePointerCapture(e.pointerId);
  };
  hit.addEventListener("pointerup", endDrag);
  hit.addEventListener("pointercancel", endDrag);

  function setMagnifier(on) {
    if (on && !magState.on && !magState.sized) { magState.radius = lensDefault(); magState.sized = true; }
    magState.on = on;
    mag.visible = on;
    $("magLayer").style.display = on ? "block" : "none";
    $("callCard").classList.toggle("show", on);
    const btn = $("magBtn");
    btn.classList.toggle("active", on);
    btn.title = on ? "Hide magnifier" : "Show magnifier";
    btn.setAttribute("aria-label", btn.title);
  }
  $("magBtn").addEventListener("click", () => setMagnifier(!magState.on));
  $("magScaleDownBtn").addEventListener("click", () => { magState.radius = Math.max(0.24, magState.radius * 0.85); });
  $("magScaleUpBtn").addEventListener("click", () => { magState.radius = Math.min(0.85, magState.radius * 1.15); });

  /* ---------- Controls ---------- */
  function startPlay() {
    if (sim.hours >= MAX_HOURS) sim.hours = 0;
    sim.playing = true;
    sim.playStart = performance.now() - (sim.hours / MAX_HOURS) * PLAY_DURATION_MS;
    setPlayIcon(true);
    lastClockH = -1;
  }
  function setPlayIcon(on) {
    const b = $("playBtn");
    b.classList.toggle("active", on);
    b.title = on ? "Pause" : "Play";
    b.setAttribute("aria-label", b.title);
    b.innerHTML = on ? '<svg viewBox="0 0 24 24" fill="none"><rect x="6.5" y="5" width="4" height="14" rx="1.2" fill="currentColor"></rect><rect x="13.5" y="5" width="4" height="14" rx="1.2" fill="currentColor"></rect></svg>' : '<svg viewBox="0 0 24 24" fill="none"><path d="M8 5v14l11-7L8 5Z" fill="currentColor"></path></svg>';
  }
  function stopPlay() { sim.playing = false; setPlayIcon(false); lastClockH = -1; }
  $("playBtn").addEventListener("click", () => (sim.playing ? stopPlay() : startPlay()));
  /* ---------- Time clock: one full turn = 120 h, ticks every 12 h ---------- */
  {
    const ticks = [];
    for (let h = 0; h < MAX_HOURS; h += 12) {
      const a = (h / MAX_HOURS) * Math.PI * 2;
      const major = h % 24 === 0;
      const r1 = major ? 23 : 26, r2 = 31;
      ticks.push(`<line x1="${(50 + Math.sin(a) * r1).toFixed(1)}" y1="${(50 - Math.cos(a) * r1).toFixed(1)}" x2="${(50 + Math.sin(a) * r2).toFixed(1)}" y2="${(50 - Math.cos(a) * r2).toFixed(1)}" stroke="rgba(255,255,255,${major ? 0.55 : 0.28})" stroke-width="${major ? 1.8 : 1.2}" stroke-linecap="round"></line>`);
    }
    $("clockTicks").innerHTML = ticks.join("");
    const clock = $("timeClock");
    let dragging = false;
    const setFromPointer = (e) => {
      const r = clock.getBoundingClientRect();
      let f = Math.atan2(e.clientX - (r.left + r.width / 2), -(e.clientY - (r.top + r.height / 2))) / (Math.PI * 2);
      if (f < 0) f += 1;
      const prev = sim.hours / MAX_HOURS;
      if (prev > 0.75 && f < 0.25) f = 1;        // no wrap-around past the end
      else if (prev < 0.25 && f > 0.75) f = 0;   // or back past the start
      sim.hours = f * MAX_HOURS;
    };
    clock.addEventListener("pointerdown", (e) => { stopPlay(); dragging = true; clock.setPointerCapture(e.pointerId); setFromPointer(e); });
    clock.addEventListener("pointermove", (e) => { if (dragging) setFromPointer(e); });
    clock.addEventListener("pointerup", () => { dragging = false; });
    clock.addEventListener("pointercancel", () => { dragging = false; });
    clock.addEventListener("keydown", (e) => {
      const step = e.shiftKey ? 24 : 1;
      if (e.key === "ArrowRight" || e.key === "ArrowUp") { e.preventDefault(); stopPlay(); sim.hours = Math.min(MAX_HOURS, sim.hours + step); }
      if (e.key === "ArrowLeft" || e.key === "ArrowDown") { e.preventDefault(); stopPlay(); sim.hours = Math.max(0, sim.hours - step); }
      if (e.key === "Home") { e.preventDefault(); stopPlay(); sim.hours = 0; }
      if (e.key === "End") { e.preventDefault(); stopPlay(); sim.hours = MAX_HOURS; }
    });
  }
  let lastClockH = -1;
  function updateClock(h) {
    if (Math.abs(h - lastClockH) < 0.05) return;
    lastClockH = h;
    const f = h / MAX_HOURS;
    $("clockHand").setAttribute("transform", `rotate(${(f * 360).toFixed(2)} 50 50)`);
    $("clockArc").setAttribute("stroke-dashoffset", (238.76 * (1 - f)).toFixed(2));
    $("timeClock").setAttribute("aria-valuenow", h.toFixed(0));
    $("dayReadout").textContent = `Day ${Math.min(MAX_HOURS / 24, Math.floor(h / 24) + 1)} of ${MAX_HOURS / 24} · ${sim.playing ? "running" : "drag the hand"}`;
  }
  $("resetBtn").addEventListener("click", () => {
    stopPlay();
    sim.hours = 0;
    events.length = 0;
    events.push({ t: 0, side: sim.side, dpsi: sim.dpsi });
    traj.dirty = true;
    const rb = $("resetBtn");
    rb.classList.add("pulse");
    setTimeout(() => rb.classList.remove("pulse"), 650);
  });

  ["mechanism", "credits"].forEach((name) => {
    const card = $(`${name}Card`), toggle = $(`${name}CardToggle`);
    toggle.addEventListener("click", () => {
      const ex = card.getAttribute("data-expanded") === "true";
      card.setAttribute("data-expanded", String(!ex));
      toggle.setAttribute("aria-expanded", String(!ex));
    });
  });

  setProgress(100, "Ready");
  const loader = $("loader");
  loader.style.opacity = "0";
  setTimeout(() => loader.remove(), 480);

  /* ---------- Frame loop ---------- */
  const prevText = {};
  const setText = (id, t) => { if (prevText[id] !== t) { prevText[id] = t; $(id).textContent = t; } };
  let lastBend = { left: 0, right: 0 };
  let chartH = -1;

  /* ---------- Theme ---------- */
  function applyTheme(theme) {
    const light = theme === "light";
    document.documentElement.dataset.theme = light ? "light" : "dark";
    const btn = $("themeBtn");
    btn.title = light ? "Switch to dark theme" : "Switch to light theme";
    btn.setAttribute("aria-label", btn.title);
    btn.innerHTML = light
      ? '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.9" stroke-linecap="round"><circle cx="12" cy="12" r="4"></circle><path d="M12 2.5v2.2M12 19.3v2.2M2.5 12h2.2M19.3 12h2.2M5.3 5.3l1.6 1.6M17.1 17.1l1.6 1.6M18.7 5.3l-1.6 1.6M6.9 17.1l-1.6 1.6"></path></svg>'
      : '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.9" stroke-linecap="round"><path d="M20 14.5A8 8 0 0 1 9.5 4a8 8 0 1 0 10.5 10.5Z"></path></svg>';
    const bg = light ? 0xe4e7ec : 0x15171c;
    scene.background.set(bg);
    scene.fog.color.set(bg);
    floor.material.color.set(light ? 0xd9d8d3 : 0x22252c);
    plants.left.pot.setTheme(light); plants.right.pot.setTheme(light);
    if (hemi) hemi.intensity = light ? 1.7 : 1.35;
    chartDirty = true;
    layout = null;
    try { localStorage.setItem("phytoposeidon-theme", light ? "light" : "dark"); } catch (e) { /* storage unavailable */ }
  }
  $("themeBtn").addEventListener("click", () => applyTheme(document.documentElement.dataset.theme === "light" ? "dark" : "light"));

  /* ---------- Water source, seedlings, chart view ---------- */
  function setupSeg(id, onPick) {
    const seg = $(id);
    seg.querySelectorAll("button").forEach((b) => b.addEventListener("click", () => {
      seg.querySelectorAll("button").forEach((x) => x.setAttribute("aria-checked", String(x === b)));
      onPick(b.dataset.v);
    }));
  }
  function paintWater() {
    const G = gradGain(sim.dpsi);
    const where = sim.side > 0 ? "Right" : "Left";
    $("waterTitle").textContent = !sim.water ? "Water source · Off" : G > 0 ? `Water source · ${where} side of pot B` : "Water source · No gradient in pot B";
    document.querySelector(".water-panel").classList.toggle("off", !sim.water);
    const wb = $("waterBtn");
    wb.classList.toggle("active", sim.water);
    wb.setAttribute("aria-checked", String(sim.water));
    wb.title = sim.water ? "Turn water off" : "Turn water on";
    wb.setAttribute("aria-label", sim.water ? "Water source on" : "Water source off");
    $("tagExpSub").textContent = !sim.water ? "Water off" : G > 0 ? `Water on the ${where.toLowerCase()}` : "Even moisture";
    $("gradReadout").textContent = `Δψ ${sim.dpsiSet.toFixed(2)} MPa`;
    $("psiWet").textContent = `${fmtPsi(PSI_WET)} MPa`;
    $("psiDry").textContent = `${fmtPsi(PSI_WET - sim.dpsi)} MPa`;
    plants.right.pot.setMoisture(sim.side, sim.dpsi);
    reservoir.set(sim.side, G);
    const ch = plants.right.chord;
    ch.normal.set(-sim.side, 0, 0); ch.constant = G > 0 ? RES_CHORD + sim.side * POT_X : 1e3;
  }
  const fmtPsi = (v) => (v < 0 ? "−" : "") + Math.abs(v).toFixed(2);
  $("goBBtn").addEventListener("click", () => setRoom("right"));
  $("waterBtn").addEventListener("click", () => {
    sim.water = !sim.water;
    sim.dpsi = sim.water ? sim.dpsiSet : 0;
    logSettings(); paintWater(); chartDirty = true; buzz(12);
  });
  // Chart card collapses to its title.
  $("dataToggle").addEventListener("click", () => {
    const card = $("dataToggle").closest(".data-card");
    const collapsed = card.classList.toggle("collapsed");
    $("dataToggle").setAttribute("aria-expanded", String(!collapsed));
    $("dataToggle").title = collapsed ? "Expand chart" : "Collapse chart";
    if (!collapsed) chartDirty = true;
  });
  setupSeg("sideSeg", (v) => { sim.side = +v; logSettings(); paintWater(); chartDirty = true; buzz(10); });
  setupSeg("repSeg", (v) => { sim.reps = +v; traj.dirty = true; setLegend(); });
  setupSeg("viewSeg", (v) => { chartView = v; setLegend(); chartDirty = true; });
  $("gradSlider").addEventListener("input", (e) => {
    sim.dpsiSet = parseFloat(e.target.value);
    sim.dpsi = sim.water ? sim.dpsiSet : 0;
    logSettings();
    paintWater();
    chartDirty = true;
  });
  paintWater();

  /* ---------- Protractor ---------- */
  const _crown = new THREE.Vector3();
  let protOn = false, protKey = "";
  $("protBtn").addEventListener("click", () => {
    protOn = !protOn;
    const b = $("protBtn");
    b.classList.toggle("active", protOn);
    b.title = protOn ? "Hide protractor" : "Show protractor";
    b.setAttribute("aria-label", b.title);
    $("protractor").style.display = protOn ? "block" : "none";
    protKey = "";
  });
  function drawProtractor(side, bendDeg) {
    const base = toScreen(_crown.set(plants[side].x, SOIL_Y, 0));
    const tip = toScreen(tipWorld(side, bendDeg));
    const R = Math.max(46, Math.hypot(tip.x - base.x, tip.y - base.y) * 1.08);
    const P = (deg, r) => [base.x + Math.sin(toRad(deg)) * r, base.y + Math.cos(toRad(deg)) * r];   // 0° = straight down
    const gold = getComputedStyle(document.documentElement).getPropertyValue("--sun").trim() || "#ffd966";
    const ink = document.documentElement.dataset.theme === "light" ? "16,19,26" : "255,255,255";
    let s = "";
    const [ax0, ay0] = P(-90, R), [ax1, ay1] = P(90, R);
    s += `<path d="M ${ax0} ${ay0} A ${R} ${R} 0 0 0 ${ax1} ${ay1} Z" fill="rgba(${ink},0.05)" stroke="rgba(${ink},0.55)" stroke-width="1.2"></path>`;
    for (let d = -90; d <= 90; d += 5) {
      const len = d % 30 === 0 ? 10 : d % 10 === 0 ? 6 : 3;
      const [x1, y1] = P(d, R), [x2, y2] = P(d, R - len);
      s += `<line x1="${x1.toFixed(1)}" y1="${y1.toFixed(1)}" x2="${x2.toFixed(1)}" y2="${y2.toFixed(1)}" stroke="rgba(${ink},0.7)" stroke-width="1"></line>`;
      if (d % 30 === 0) {
        const [tx, ty] = P(d, R - 20);
        s += `<text x="${tx.toFixed(1)}" y="${ty.toFixed(1)}" font-size="9.5" fill="rgba(${ink},0.85)" text-anchor="middle" dominant-baseline="middle">${Math.abs(d)}°</text>`;
      }
    }
    const [vx, vy] = P(0, R);
    s += `<line x1="${base.x.toFixed(1)}" y1="${base.y.toFixed(1)}" x2="${vx.toFixed(1)}" y2="${vy.toFixed(1)}" stroke="rgba(${ink},0.45)" stroke-dasharray="4 4"></line>`;
    const a = bendDeg;
    if (Math.abs(a) > 0.5) {
      const r2 = R * 0.42;
      const [wx, wy] = P(a, r2), [w0x, w0y] = P(0, r2);
      s += `<path d="M ${base.x} ${base.y} L ${w0x} ${w0y} A ${r2} ${r2} 0 0 ${a > 0 ? 0 : 1} ${wx} ${wy} Z" fill="${gold}" fill-opacity="0.28"></path>`;
    }
    const [mx, my] = P(a, R);
    s += `<line x1="${base.x.toFixed(1)}" y1="${base.y.toFixed(1)}" x2="${mx.toFixed(1)}" y2="${my.toFixed(1)}" stroke="${gold}" stroke-width="2"></line>`;
    s += `<circle cx="${base.x.toFixed(1)}" cy="${base.y.toFixed(1)}" r="3" fill="${gold}"></circle>`;
    const [lx, ly] = P(a >= 0 ? a + 14 : a - 14, R + 16);
    s += `<g transform="translate(${lx.toFixed(1)} ${ly.toFixed(1)})"><rect x="-30" y="-11" width="60" height="22" rx="11" fill="rgba(12,14,18,0.78)" stroke="rgba(${ink},0.2)"></rect><text x="0" y="1" font-size="12" font-weight="600" fill="#fff" text-anchor="middle" dominant-baseline="middle">θ ${Math.abs(a).toFixed(1)}°</text></g>`;
    $("protractor").innerHTML = s;
  }

  /* ---------- Data panel ---------- */
  const readings = [];
  let chartDirty = true, chartLast = 0;
  const chart = $("chart");
  const cssVar = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  let chartView = "bend";
  const GOLD = "#4fc0ee", GREY = "#9aa0aa", ORANGE = "#e0a252";
  function setLegend() {
    const L = $("legend");
    L.innerHTML = chartView === "bend"
      ? `<span style="--c:${GREY}">A · Control</span><span style="--c:${GOLD}">B · Experimental</span>${sim.reps === 5 ? `<span class="band" style="--c:${GOLD}">± 1 SD (n = 5)</span>` : ""}<span class="tri">Setting changed</span>`
      : `<span style="--c:${GOLD}">B · left flank</span><span style="--c:${ORANGE}">B · right flank</span><span style="--c:${GREY}">A · both flanks</span>`;
    $("dataTitle").innerHTML = chartView === "bend" ? "Root tip angle over time <small>(+ right, − left)</small>" : "Cortex elongation rate <small>(% of control)</small>";
  }
  function drawChart(target = chart, scale = Math.min(window.devicePixelRatio || 1, 2), withBg = false) {
    const W = target === chart ? chart.clientWidth : target.width / scale;
    const H = target === chart ? chart.clientHeight : target.height / scale;
    if (target === chart) { chart.width = Math.round(W * scale); chart.height = Math.round(H * scale); }
    const g = target.getContext("2d");
    g.setTransform(scale, 0, 0, scale, 0, 0);
    g.clearRect(0, 0, W, H);
    if (withBg) { g.fillStyle = cssVar("--surface") || "#131519"; g.fillRect(0, 0, W, H); }
    const ink3 = cssVar("--ink-3"), line = cssVar("--line"), font = "11px 'IBM Plex Sans', Roboto, sans-serif";
    const m = { l: 42, r: 10, t: 10, b: 30 };
    const pw = W - m.l - m.r, ph = H - m.t - m.b;
    const growth = chartView === "growth";
    const yMin = growth ? 90 : -90, yMax = growth ? 110 : 90;
    const X = (h) => m.l + (h / MAX_HOURS) * pw;
    const Y = (v) => m.t + ph - ((Math.max(yMin, Math.min(yMax, v)) - yMin) / (yMax - yMin)) * ph;
    g.font = font; g.lineWidth = 1;
    const ticks = growth ? [90, 95, 100, 105, 110] : [-90, -45, 0, 45, 90];
    ticks.forEach((d) => {
      g.strokeStyle = (growth ? d === 100 : d === 0) ? ink3 : line;
      g.beginPath(); g.moveTo(m.l, Y(d) + 0.5); g.lineTo(W - m.r, Y(d) + 0.5); g.stroke();
      g.fillStyle = ink3; g.textAlign = "right"; g.textBaseline = "middle";
      g.fillText(growth ? d + "%" : (d > 0 ? "+" : "") + d + "°", m.l - 6, Y(d));
    });
    g.textBaseline = "top";
    for (let h = 0; h <= MAX_HOURS; h += (W < 520 ? 24 : 12)) {
      g.textAlign = h === 0 ? "left" : h >= MAX_HOURS ? "right" : "center";
      g.fillText(h + " h", X(h), H - m.b + 10);
    }
    // Lag after the latest change (shaded).
    const lastEv = events[events.length - 1];
    const lagEnd = Math.min(MAX_HOURS, lastEv.t + LAG);
    g.fillStyle = cssVar("--hover"); g.fillRect(X(lastEv.t), m.t, Math.max(1, X(lagEnd) - X(lastEv.t)), ph);
    g.fillStyle = ink3; g.textAlign = "left"; g.textBaseline = "top"; g.fillText("lag", X(lagEnd) + 3, m.t + 2);
    events.forEach((ev, i) => {
      if (i === 0) return;
      const x = X(ev.t);
      g.fillStyle = GOLD;
      g.beginPath(); g.moveTo(x, m.t + ph + 1); g.lineTo(x - 4, m.t + ph + 7); g.lineTo(x + 4, m.t + ph + 7); g.closePath(); g.fill();
    });
    const path = (fn, from, to, c, dashed, w = 2.4) => {
      g.strokeStyle = c; g.lineWidth = dashed ? 1.4 : w; g.globalAlpha = dashed ? 0.45 : 1; g.setLineDash(dashed ? [5, 5] : []);
      g.beginPath();
      for (let h = from, k = 0; h <= to + 1e-6; h += 0.5, k++) { const x = X(h), y = Y(fn(h)); k ? g.lineTo(x, y) : g.moveTo(x, y); }
      g.stroke(); g.globalAlpha = 1; g.setLineDash([]);
    };
    const dot = (v, c) => { g.fillStyle = c; g.beginPath(); g.arc(X(sim.hours), Y(v), 4, 0, Math.PI * 2); g.fill(); };
    if (!growth) {
      const A = showA(), B = showB();
      if (sim.reps === 5) {
        // ±1 SD band for B (A stays at 0 for every seedling).
        g.fillStyle = GOLD; g.globalAlpha = 0.16; g.beginPath();
        for (let h = 0; h <= MAX_HOURS; h += 1) g.lineTo(X(h), Y(angleAt(B, h) + angleAt(traj.Bsd, h)));
        for (let h = MAX_HOURS; h >= 0; h -= 1) g.lineTo(X(h), Y(angleAt(B, h) - angleAt(traj.Bsd, h)));
        g.closePath(); g.fill(); g.globalAlpha = 1;
        traj.Bset.forEach((arr) => path((h) => angleAt(arr, h), 0, sim.hours, GOLD, false, 0.8));
      }
      const eq = equilibriumDeg(90 * sim.side, sim.dpsi);
      {
        g.setLineDash([2, 4]); g.strokeStyle = GOLD; g.globalAlpha = 0.6; g.lineWidth = 1.2;
        g.beginPath(); g.moveTo(X(sim.hours), Y(eq)); g.lineTo(W - m.r, Y(eq)); g.stroke();
        g.globalAlpha = 1; g.setLineDash([]);
        g.fillStyle = GOLD; g.textAlign = "right"; g.textBaseline = eq >= 0 ? "bottom" : "top";
        g.fillText(`eq. ${eq >= 0 ? "+" : "−"}${Math.abs(eq).toFixed(1)}°`, W - m.r, Y(eq) + (eq >= 0 ? -3 : 3));
      }
      path((h) => angleAt(A, h), 0, sim.hours, GREY, false);
      path((h) => angleAt(B, h), sim.hours, MAX_HOURS, GOLD, true);
      path((h) => angleAt(B, h), 0, sim.hours, GOLD, false);
      dot(angleAt(A, sim.hours), GREY); dot(angleAt(B, sim.hours), GOLD);
      readings.forEach((r, i) => {
        g.strokeStyle = GOLD; g.lineWidth = 1.6; g.fillStyle = cssVar("--surface") || "#131519";
        if (r.bSd) { g.beginPath(); g.moveTo(X(r.t), Y(r.b - r.bSd)); g.lineTo(X(r.t), Y(r.b + r.bSd)); g.stroke(); }
        g.beginPath(); g.arc(X(r.t), Y(r.b), 5, 0, Math.PI * 2); g.fill(); g.stroke();
        g.fillStyle = ink3; g.font = "600 10px 'IBM Plex Sans', Roboto, sans-serif"; g.textAlign = "center"; g.textBaseline = r.b >= 0 ? "bottom" : "top";
        g.fillText(String(i + 1), X(r.t), Y(r.b) + (r.b >= 0 ? -8 - (r.bSd ? 6 : 0) : 8 + (r.bSd ? 6 : 0)));
        g.font = font;
      });
    } else {
      // Differential growth: turning right means the left flank is outgrowing the right.
      const R = traj.rateB;
      const left = (h) => 100 + flankPct(angleAt(R, h)), right = (h) => 100 - flankPct(angleAt(R, h));
      path(() => 100, 0, sim.hours, GREY, false);
      path(left, sim.hours, MAX_HOURS, GOLD, true); path(right, sim.hours, MAX_HOURS, ORANGE, true);
      path(left, 0, sim.hours, GOLD, false); path(right, 0, sim.hours, ORANGE, false);
      dot(100, GREY); dot(left(sim.hours), GOLD); dot(right(sim.hours), ORANGE);
      g.fillStyle = ink3; g.textAlign = "right"; g.textBaseline = "bottom";
      g.fillText("faster", W - m.r - 2, Y(yMax) + 12); g.textBaseline = "top"; g.fillText("slower", W - m.r - 2, Y(yMin) - 12);
    }
  }
  new ResizeObserver(() => { chartDirty = true; }).observe(chart);
  setLegend();

  function updateReasoning() {
    const h = sim.hours, th = lastBend.right;
    const seen = settingsAt(h - LAG);
    const G = gradGain(seen.dpsi), phi = 90 * seen.side;
    const eq = equilibriumDeg(90 * sim.side, sim.dpsi);
    const fmt = (x) => (Math.abs(x) < 0.05 ? "0°" : `${Math.abs(x).toFixed(1)}° ${x > 0 ? "right" : "left"}`);
    const sg = (x, d, paren = false) => { const s = (x < 0 ? "−" : "") + Math.abs(x).toFixed(d) + "°"; return paren && x < 0 ? `(${s})` : s; };
    const rate = angleAt(traj.rateB, h);
    const P = [];
    P.push(`<p><strong>Why root B is at ${fmt(th)}${sim.reps === 5 ? ` (mean of 5, SD ${sdNow("B").toFixed(1)}°)` : ""}</strong></p>`);
    if (G <= 0) {
      P.push(`<p>There is no moisture gradient across pot B, so both flanks of the root tip sit at the same water potential. No lateral signal forms and gravitropism keeps the root growing straight down. ${Math.abs(th) > 0.5 ? "The bend left from the earlier gradient is being straightened by gravitropism." : ""}</p>`);
    } else {
      P.push(`<p>The water is on the ${seen.side > 0 ? "right" : "left"}, 90° from the root's vertical axis, and the water potential falls by ${seen.dpsi.toFixed(2)} MPa across the pot. The hydrotropic signal is proportional to ${(G).toFixed(2)} × sin(${sg(phi, 0)} − ${sg(th, 1, true)}) = ${(G * Math.sin(toRad(phi - th))).toFixed(2)}, the gradient's sideways component across the tip. It shrinks as the root turns toward the water.</p>`);
      const amy = angleAt(traj.amyB, h);
      P.push(`<p>The gradient also breaks down the starch-filled amyloplasts that sense gravity in the root cap columella. They are now at ${(amy * 100).toFixed(0)}% of normal${amy > amyTarget(G) + 0.02 ? `, falling toward ${(amyTarget(G) * 100).toFixed(0)}%` : ""}, so the gravitropic pull is weaker by the same amount. Gravitropism pulls back in proportion to sin(${sg(th, 1)}) = ${Math.sin(toRad(th)).toFixed(2)}. The root settles where the two balance, here ${fmt(eq)}.</p>`);
    }
    let state;
    if (h <= LAG) state = `Lag phase (${(LAG * 60).toFixed(0)} min): the root tip has sensed the moisture gradient, but the cortex has not yet started to grow unevenly.`;
    else if (Math.abs(rate) < 0.15) state = "Not turning now: both flanks of the elongation zone are growing at the same rate.";
    else {
      const outer = rate > 0 ? "left" : "right", toWater = G > 0 && Math.sign(rate) === Math.sign(seen.side);
      const pct = flankPct(Math.abs(rate));
      state = `Turning ${Math.abs(rate).toFixed(1)}° per hour to the ${rate > 0 ? "right" : "left"}. To bend a root about ${(ROOT_D * 1000).toFixed(0)} µm thick that grows ${ROOT_V.toFixed(1)} mm per hour, the ${outer}${toWater ? " (drier)" : ""} flank must elongate ${pct.toFixed(1)}% faster than the mean and the other flank ${pct.toFixed(1)}% slower (see Cortex growth). In hydrotropism this differential growth happens in the cortex of the elongation zone.`;
    }
    P.push(`<p>${state}</p>`);
    if (sim.reps === 5) P.push(`<p>Five seedlings differ slightly in sensitivity and lag, as real seedlings do. The line is their mean and the band is ±1 standard deviation.</p>`);
    P.push(`<p>Pot A is watered evenly, so there is no gradient across its roots and they grow straight down by gravitropism alone.</p>`);
    $("reasoning").innerHTML = P.join("");
  }

  function renderReadings() {
    const wrap = $("resultsWrap");
    if (!readings.length) { wrap.innerHTML = '<p class="results-empty">No readings yet. Pause at any time and tap the record button.</p>'; return; }
    wrap.innerHTML = `<div class="scroll"><table class="results"><thead><tr><th>#</th><th>Time (h)</th><th>A (°)</th><th>B (°)</th><th>B eq. (°)</th><th>Water</th><th>Δψ (MPa)</th></tr></thead><tbody>${
      readings.map((r, i) => `<tr><td>${i + 1}</td><td>${r.t.toFixed(1)}</td><td>${r.a.toFixed(1)}</td><td>${r.b.toFixed(1)}${r.n === 5 ? ` ± ${r.bSd.toFixed(1)}` : ""}</td><td>${r.eq.toFixed(1)}</td><td>${r.dpsi > 0 ? (r.side > 0 ? "Right" : "Left") : "—"}</td><td>${r.dpsi.toFixed(2)}</td></tr>`).join("")
    }</tbody></table></div>`;
  }
  const csv = () => ["reading,time_h,n_seedlings,tip_angle_A_deg,tip_angle_B_mean_deg,tip_angle_B_sd_deg,equilibrium_B_deg,water_side,delta_psi_MPa"]
    .concat(readings.map((r, i) => [i + 1, r.t.toFixed(2), r.n, r.a.toFixed(2), r.b.toFixed(2), r.n === 5 ? r.bSd.toFixed(2) : "", r.eq.toFixed(2), r.side > 0 ? "right" : "left", r.dpsi.toFixed(2)].join(","))).join("\n");
  const note = (msg) => { $("dataNote").textContent = msg; };
  $("freezeBtn").addEventListener("click", () => {
    const prev = readings[readings.length - 1];
    if (prev && Math.abs(prev.t - sim.hours) < 0.05 && prev.side === sim.side && prev.dpsi === sim.dpsi && prev.n === sim.reps) {
      return note(`Already recorded at t = ${sim.hours.toFixed(1)} h. Move the clock or change a setting first.`);
    }
    if (traj.dirty) integrate();   // read from the model at this exact hour, not the last drawn frame
    readings.push({ t: sim.hours, a: angleAt(showA(), sim.hours), b: angleAt(showB(), sim.hours), bSd: sim.reps === 5 ? sdNow("B") : 0, n: sim.reps,
      side: sim.side, dpsi: sim.dpsi, eq: equilibriumDeg(90 * sim.side, sim.dpsi) });
    renderReadings(); chartDirty = true; buzz(15);
    note(`Reading ${readings.length} recorded at t = ${sim.hours.toFixed(1)} h.`);
  });
  $("clearBtn").addEventListener("click", () => { readings.length = 0; renderReadings(); chartDirty = true; note("Readings cleared."); });
  $("copyCsvBtn").addEventListener("click", async () => {
    if (!readings.length) return note("Record at least one reading first.");
    try { await navigator.clipboard.writeText(csv()); note("CSV copied. Paste it into a spreadsheet."); }
    catch (e) { note("Copy was blocked here. Use “Save CSV”, or select the table and copy it."); }
  });
  // Published page: the viewer's download capability. Standalone file: a plain link.
  const dlReady = (window.claude && window.claude.use) ? window.claude.use("downloads").catch(() => null) : Promise.resolve(null);
  async function download(name, blob) {
    const dl = await dlReady;
    if (dl) {
      try { await dl.save({ filename: name, data: blob }); note(`${name} saved.`); }
      catch (e) {
        const code = e && e.code;
        note(code === "declined" ? "Save cancelled." : code === "host_unavailable" ? "This app can't save files. Open the page in a browser to save, or use Copy CSV." : "Saving isn't available here. Use Copy CSV instead.");
      }
      return;
    }
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob); a.download = name; document.body.appendChild(a); a.click(); a.remove();
    note(`${name} saved.`);
  }
  $("saveCsvBtn").addEventListener("click", () => {
    if (!readings.length) return note("Record at least one reading first.");
    download("phyto-poseidon-readings.csv", new Blob([csv()], { type: "text/csv" }));
  });
  $("savePngBtn").addEventListener("click", () => {
    const sc = 2, W = Math.max(640, chart.clientWidth), H = 300;
    const c = document.createElement("canvas"); c.width = W * sc; c.height = H * sc;
    drawChart(c, sc, true);
    c.toBlob((blob) => { if (blob) download("phyto-poseidon-chart.png", blob); }, "image/png");
  });

  /* ---------- Performance: pause when hidden, adapt resolution ---------- */
  let onScreen = true;
  new IntersectionObserver((es) => { onScreen = es[0].isIntersecting; }).observe(stage);
  let perfT = 0, perfN = 0, perfLast = performance.now(), prCap = Math.min(window.devicePixelRatio || 1, 2);
  function adaptQuality(now) {
    perfT += now - perfLast; perfLast = now; perfN++;
    if (perfN < 90) return;
    const avg = perfT / perfN; perfT = 0; perfN = 0;
    const pr = renderer.getPixelRatio();
    const floorPR = Math.min(window.devicePixelRatio || 1, 1.5);   // below this, thin roots break up
    if (avg > 30 && pr > floorPR) { prCap = Math.max(floorPR, pr - 0.25); renderer.setPixelRatio(prCap); renderer.setSize(stage.clientWidth, stage.clientHeight, false); }
  }

  /* ---------- Text description for screen readers ---------- */
  let summaryLast = 0, summaryText = "";
  function updateSummary(now) {
    if (now - summaryLast < 2500) return;
    summaryLast = now;
    const side = (d) => (Math.abs(d) < 0.5 ? "growing straight down" : `turned ${Math.abs(d).toFixed(0)}° to the ${d > 0 ? "right" : "left"}`);
    const s = `At ${sim.hours.toFixed(0)} hours, the root of plant A, in evenly moist soil, is ${side(lastBend.left)}. The root of plant B is ${side(lastBend.right)}, with water on the ${sim.side > 0 ? "right" : "left"} and a gradient of ${sim.dpsi.toFixed(2)} megapascals.`;
    if (s !== summaryText) { summaryText = s; $("sceneSummary").textContent = s; }
  }

  let savedTheme = "dark";
  try { savedTheme = localStorage.getItem("phytoposeidon-theme") || "dark"; } catch (e) { /* storage unavailable */ }
  applyTheme(savedTheme);

  const _tagCtl = new THREE.Vector3(-POT_X, 0, POT_RIM * 0.75), _tagExp = new THREE.Vector3(POT_X, 0, POT_RIM * 0.75), _tagWater = new THREE.Vector3();
  let shadowKey = "", hiddenAt = 0, lastNow = 0, lensLast = 0;
  const _right = new THREE.Vector3(), _edge = new THREE.Vector3();
  function frame(now) {
    // A hidden tab freezes the clock. When the 3D view is only scrolled out of
    // sight, time, chart and readouts keep updating; just the scene isn't drawn.
    if (document.hidden) {
      if (!hiddenAt) hiddenAt = now;
      perfLast = now; requestAnimationFrame(frame); return;
    }
    if (hiddenAt) { if (sim.playing) sim.playStart += now - hiddenAt; hiddenAt = 0; }
    const drawScene = onScreen;
    if (drawScene) adaptQuality(now); else perfLast = now;
    // Time and the eased plant response.
    if (sim.playing) {
      sim.hours = Math.min(MAX_HOURS, ((now - sim.playStart) / PLAY_DURATION_MS) * MAX_HOURS);
      if (sim.hours >= MAX_HOURS) stopPlay();
    }
    if (traj.dirty) { integrate(); chartDirty = true; }

    const bend = { left: angleAt(showA(), sim.hours), right: angleAt(showB(), sim.hours) };
    lastBend = bend;
    const rg = rootGrowth(sim.hours) * ROOT_DEPTH, sg = shootGrowth(sim.hours);
    for (const side of ["left", "right"]) {
      plants[side].roots.group.scale.setScalar(rg);
      plants[side].roots.uniforms.uBend.value = toRad(bend[side]);
      plants[side].shoot.scale.setScalar(sg);
    }
    // Homing signal: while water flows, root B pulses and its tip pings.
    {
      const ev = settingsAt(sim.hours - LAG), G = sim.hours > LAG ? gradGain(ev.dpsi) : 0;
      const pulse = G > 0 ? 0.35 + 0.65 * G : 0;
      const tt = now / 1000;
      plants.right.roots.uniforms.uTime.value = tt;
      plants.right.roots.uniforms.uPulse.value = pulse * 0.95;
      const tip = tipWorld("right", bend.right);
      pings.forEach((sp, i) => {
        const ph = ((tt / 1.6) + i * 0.5) % 1;
        sp.visible = pulse > 0 && drawScene;
        sp.position.copy(tip);
        sp.scale.setScalar(0.06 + ph * 0.42 * (0.6 + 0.4 * G));
        sp.material.opacity = (1 - ph) * 0.9 * pulse;
      });
    }
    // Reservoir ripples and the water flowing into the soil.
    if (drawScene) reservoir.update(Math.min(0.05, (now - (lastNow || now)) / 1000), document.documentElement.dataset.theme === "light");
    lastNow = now;

    if (drawScene) {
      if (navHold) orbitStep(navHold);
      stepCamAnim(now);
      controls.update();
      // Magnifier faces the camera; the lens redraws at ~25 fps.
      if (magState.on) {
        mag.position.copy(magState.pos);
        mag.quaternion.copy(camera.quaternion);
        mag.scale.setScalar(magState.radius);
        const lensC = toScreen(magState.pos);
        _right.set(1, 0, 0).applyQuaternion(camera.quaternion).multiplyScalar(magState.radius);
        const edge = toScreen(_edge.copy(magState.pos).add(_right));
        const lensPx = Math.hypot(edge.x - lensC.x, edge.y - lensC.y);
        updateCellView(lensC, lensPx, bend);
        if (now - lensLast > 40) { lensLast = now; drawLens(bend, now); }
      }
      // Shadows: only when the roots or shoots change.
      const sk = `${rg.toFixed(3)}|${bend.left.toFixed(2)}|${bend.right.toFixed(2)}`;
      if (sk !== shadowKey) { shadowKey = sk; renderer.shadowMap.needsUpdate = true; }
      renderer.render(scene, camera);
    }

    // Screen-space labels.
    const wet = gradGain(sim.dpsi) > 0;
    $("tagWater").style.display = wet ? "" : "none";
    if (wet) placeTag($("tagWater"), _tagWater.set(POT_X + sim.side * (POT_RIM - 0.06), POT_H + 0.04, 0), -14);
    placeTag($("tagControl"), _tagCtl, 18);
    placeTag($("tagExperimental"), _tagExp, 18);

    // Readouts.
    setText("hoursReadout", `t = ${sim.hours.toFixed(1)} h`);
    updateClock(sim.hours);
    setText("controlStat", statText("A", bend.left));
    setText("experimentalStat", statText("B", bend.right));

    if (protOn) {
      const ps = view === "left" ? "left" : "right";
      const pk = `${ps}|${bend[ps].toFixed(2)}|${rg.toFixed(3)}|${camera.position.x.toFixed(3)}|${camera.position.y.toFixed(3)}|${camera.position.z.toFixed(3)}|${stage.clientWidth}|${document.documentElement.dataset.theme}`;
      if (pk !== protKey) { protKey = pk; drawProtractor(ps, bend[ps]); }
    }
    if (chartDirty || (now - chartLast > 150 && chartH !== sim.hours)) {
      drawChart(); chartDirty = false; chartLast = now; chartH = sim.hours;
      updateReasoning();
    }
    updateSummary(now);

    requestAnimationFrame(frame);
  }

  const sdNow = (room) => angleAt(room === "A" ? traj.Asd : traj.Bsd, sim.hours);
  /* ---------- Cellular view: lock onto the nearest root tip ---------- */
  function updateCellView(c, rPx, bend) {
    hit.setAttribute("cx", c.x.toFixed(1));
    hit.setAttribute("cy", c.y.toFixed(1));
    hit.setAttribute("r", Math.max(24, rPx).toFixed(1));
    const scr = { left: toScreen(tipWorld("left", bend.left)), right: toScreen(tipWorld("right", bend.right)) };
    const range = rPx * 2.4 + 60;
    const dist = (sd) => Math.hypot(scr[sd].x - c.x, scr[sd].y - c.y);
    if (!magState.locked || dist(magState.locked) > range) {
      magState.locked = null;
      const best = view || (dist("left") < dist("right") ? "left" : "right");
      if (dist(best) < range) magState.locked = best;
    }
    const lead = $("magLeaderLine"), dot = $("magApexDot"), sd = magState.locked;
    if (sd) {
      const a = scr[sd], d = Math.hypot(a.x - c.x, a.y - c.y) || 1;
      lead.setAttribute("x1", (c.x + ((a.x - c.x) / d) * rPx).toFixed(1));
      lead.setAttribute("y1", (c.y + ((a.y - c.y) / d) * rPx).toFixed(1));
      lead.setAttribute("x2", a.x.toFixed(1)); lead.setAttribute("y2", a.y.toFixed(1));
      dot.setAttribute("cx", a.x.toFixed(1)); dot.setAttribute("cy", a.y.toFixed(1));
      const show = d > rPx ? "1" : "0";
      lead.setAttribute("opacity", show); dot.setAttribute("opacity", show);
    } else { lead.setAttribute("opacity", "0"); dot.setAttribute("opacity", "0"); }
    // Call card.
    let t = "Magnifier", l1 = "Drag the lens over a root tip", l2 = "It locks onto the nearest one";
    if (sd === "left") { t = "Root tip · Pot A"; l1 = `${Math.abs(bend.left) < 0.05 ? "0" : Math.abs(bend.left).toFixed(1)}° · even moisture`; l2 = "Cortex grows evenly · amyloplasts 100%"; }
    if (sd === "right") {
      const ev = settingsAt(sim.hours - LAG), G = gradGain(ev.dpsi), amy = angleAt(traj.amyB, sim.hours), rate = angleAt(traj.rateB, sim.hours);
      t = "Root tip · Pot B";
      const ang = Math.abs(bend.right) < 0.05 ? "0°" : `${Math.abs(bend.right).toFixed(1)}° ${bend.right > 0 ? "right" : "left"}`;
      l1 = G > 0 ? `${ang} · water on the ${ev.side > 0 ? "right" : "left"}` : `${ang} · no gradient sensed`;
      l2 = sim.hours > LAG && Math.abs(rate) >= 0.15 ? `${rate > 0 ? "Left" : "Right"} flank +${flankPct(Math.abs(rate)).toFixed(1)}% · amyloplasts ${(amy * 100).toFixed(0)}%` : `Flanks even · amyloplasts ${(amy * 100).toFixed(0)}%`;
    }
    setText("callTitle", t); setText("callLine1", l1); setText("callLine2", l2);
  }
  function drawLens(bend, now) {
    const light = document.documentElement.dataset.theme === "light";
    const sd = magState.locked;
    if (!sd) return lens.render(renderer, null, now);
    if (sd === "left") return lens.render(renderer, { theta: bend.left, side: 1, G: 0, amy: 1, light }, now);
    const ev = settingsAt(sim.hours - LAG);
    lens.render(renderer, { theta: bend.right, side: ev.side, G: gradGain(ev.dpsi), amy: angleAt(traj.amyB, sim.hours), light }, now);
  }
  let lensEmpty = null;

  function statText(room, tipDeg) {
    if (room === "A") return Math.abs(tipDeg) < 0.5 ? "0° · even moisture, growing straight down" : `${Math.abs(tipDeg).toFixed(1)}° · straightening by gravitropism`;
    const ev = settingsAt(sim.hours), G = gradGain(ev.dpsi);
    const where = ev.side > 0 ? "right" : "left";
    if (G > 0 && sim.hours - ev.t <= LAG && Math.abs(tipDeg) < 0.5) return `0° · lag phase, water on the ${where}`;
    const eq = equilibriumDeg(90 * ev.side, ev.dpsi);
    const pm = sim.reps === 5 ? ` ± ${sdNow("B").toFixed(1)}° (n = 5)` : "";
    if (Math.abs(tipDeg) < 0.5 && Math.abs(eq) < 0.5) return "0° · no gradient, growing straight down";
    if (Math.abs(eq) < 0.5) return `${Math.abs(tipDeg).toFixed(1)}° ${tipDeg > 0 ? "right" : "left"}${sim.reps === 5 ? ` ± ${sdNow("B").toFixed(1)}° (n = 5)` : ""} · no gradient, straightening by gravitropism`;
    const lean = tipDeg > 0 ? "right" : "left";
    return `${Math.abs(tipDeg).toFixed(1)}° ${lean}${pm} · equilibrium ${Math.abs(eq).toFixed(1)}° ${eq >= 0 ? "right" : "left"}`;
  }

  requestAnimationFrame(frame);
}

/* ---------------- Light-tight partition ----------------
   Matte blackout panel in an anodised aluminium frame on two T-feet. */
function boardTexture() {
  const c = document.createElement("canvas");
  c.width = c.height = 256;
  const g = c.getContext("2d");
  const img = g.createImageData(256, 256);
  for (let y = 0; y < 256; y++) {
    for (let x = 0; x < 256; x++) {
      const grain = Math.sin(y * 0.9 + Math.sin(x * 0.05) * 3) * 6; // faint vertical fibre
      const v = 128 + (Math.random() - 0.5) * 46 + grain;
      const i = (y * 256 + x) * 4;
      img.data[i] = img.data[i + 1] = img.data[i + 2] = v;
      img.data[i + 3] = 255;
    }
  }
  g.putImageData(img, 0, 0);
  const t = new THREE.CanvasTexture(c);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.repeat.set(3, 3);
  return t;
}

function makeDivider(envTex) {
  const H = 2.4, D = 2.4, T = 0.045;
  const g = new THREE.Group();
  const add = (geo, mat, x, y, z) => {
    const m = new THREE.Mesh(geo, mat);
    m.position.set(x, y, z);
    m.castShadow = true; m.receiveShadow = true;
    g.add(m);
    return m;
  };
  const tex = boardTexture();
  const board = new THREE.MeshStandardMaterial({ color: 0x24272c, roughness: 0.92, roughnessMap: tex, bumpMap: tex, bumpScale: 0.35 });
  const alu = new THREE.MeshStandardMaterial({ color: 0xc3c7cd, metalness: 1, roughness: 0.34, envMap: envTex, envMapIntensity: 0.85 });
  const dark = new THREE.MeshStandardMaterial({ color: 0x15171a, roughness: 0.7, envMap: envTex, envMapIntensity: 0.3 });
  const steel = new THREE.MeshStandardMaterial({ color: 0x9aa0a8, metalness: 1, roughness: 0.25, envMap: envTex, envMapIntensity: 1 });

  const base = 0.06;
  add(new RoundedBoxGeometry(T, H - 0.1, D - 0.1, 2, 0.01), board, 0, base + H / 2, 0);            // panel
  add(new RoundedBoxGeometry(0.075, 0.06, D, 3, 0.012), alu, 0, base + H, 0);                       // top rail
  add(new RoundedBoxGeometry(0.075, 0.06, D, 3, 0.012), alu, 0, base + 0.03, 0);                    // bottom rail
  for (const z of [-D / 2, D / 2]) {
    add(new RoundedBoxGeometry(0.075, H + 0.06, 0.06, 3, 0.012), alu, 0, base + H / 2, z);          // posts
    add(new RoundedBoxGeometry(0.06, 0.006, 0.062, 2, 0.002), dark, 0, base + H * 0.5, z);          // shadow-gap groove
    add(new RoundedBoxGeometry(0.58, 0.035, 0.1, 3, 0.012), alu, 0, 0.0175, z);                     // T-foot
    for (const x of [-0.26, 0.26]) add(new THREE.CylinderGeometry(0.03, 0.034, 0.012, 20), dark, x, 0.006, z); // rubber pads
    for (const y of [0.35, H * 0.5 + 0.06, H - 0.25]) {                                              // fixing screws
      for (const sx of [-1, 1]) {
        const s = add(new THREE.CylinderGeometry(0.009, 0.009, 0.006, 16), steel, sx * 0.0405, base + y, z);
        s.rotation.z = Math.PI / 2;
      }
    }
  }
  // Brush seal where the panel meets the frame.
  add(new THREE.BoxGeometry(0.05, 0.012, D - 0.1), dark, 0, base + 0.065, 0);
  add(new THREE.BoxGeometry(0.05, 0.012, D - 0.1), dark, 0, base + H - 0.035, 0);
  return g;
}

/* ---------------- Start ---------------- */
try {
  setProgress(4, "Loading 3D engine");
  THREE = await import("three");
  setProgress(16, "Loading 3D engine");
  ({ GLTFLoader } = await import("three/addons/loaders/GLTFLoader.js"));
  ({ OrbitControls } = await import("three/addons/controls/OrbitControls.js"));
  ({ RoundedBoxGeometry } = await import("three/addons/geometries/RoundedBoxGeometry.js"));
  ({ RoomEnvironment } = await import("three/addons/environments/RoomEnvironment.js"));
  await init();
} catch (e) { fail(e); }
