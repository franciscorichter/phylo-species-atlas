"use strict";
// Species Phylogenies Catalogue: every described species as dots around the tree of life of
// the 62 groups. Reads data.json (the paper's numbers, written by build.py) and
// data/tree_of_groups.json (the order of the groups in the Open Tree of Life synthesis).
// ---------------------------------------------------------------- constants
const TREES_BASE = "../standardized/trees/"; // patched by CI for Pages deploy
const NAMES_BASE = "../standardized/names/"; // patched by CI for Pages deploy
const TREE_MAX = 6000;             // a group's tree shows at most this many tips
const LEVELS = ["named", "molecular", "dated", "uncertainty"];
const fmt = new Intl.NumberFormat("en-US");
const reduce = matchMedia("(prefers-reduced-motion: reduce)").matches;
const hash = (n) => { const x = Math.sin(n * 127.1 + 311.7) * 43758.5453; return x - Math.floor(x); };
const ease = (x) => (x < 0.5 ? 4 * x * x * x : 1 - (-2 * x + 2) ** 3 / 2);


// ---------------------------------------------------------------- load
// The dots of each level go to the groups by the largest-remainder method, so a level's
// dots sum to round(level total / SCALE); a group's dots at one level never exceed its dots
// at the level before. Bacteria and Archaea get their own dots, outside every eukaryotic total.
const SCALE = 500;
function largestRemainder(values, scale) {
  const target = Math.round(values.reduce((a, b) => a + b, 0) / scale);
  const raw = values.map(v => v / scale), base = raw.map(Math.floor);
  let rest = target - base.reduce((a, b) => a + b, 0);
  raw.map((x, i) => [x - base[i], i]).sort((a, b) => b[0] - a[0]).slice(0, rest).forEach(([, i]) => base[i]++);
  return base;
}
async function loadData() {
  const [d, t] = await Promise.all([
    fetch("data.json").then(r => { if (!r.ok) throw new Error("data.json " + r.status); return r.json(); }),
    fetch("data/tree_of_groups.json").then(r => { if (!r.ok) throw new Error("tree_of_groups.json " + r.status); return r.json(); }),
  ]);
  const trees = Object.fromEntries(d.trees.map(x => [x.filename, x]));
  const cov = d.coverage;
  const dots = Object.fromEntries(cov.map(r => [r.group, {}]));
  for (const part of [cov.filter(r => !r.is_prokaryotic), cov.filter(r => r.is_prokaryotic)]) {
    largestRemainder(part.map(r => r.described_species), SCALE).forEach((n, i) => { dots[part[i].group].all = n; });
    let prev = "all";
    for (const lv of LEVELS) {
      const vals = part.map(r => (r.levels && r.levels[lv] && r.levels[lv].n) || 0);
      largestRemainder(vals, SCALE).forEach((n, i) => { const g = dots[part[i].group]; g[lv] = Math.min(n, g[prev]); });
      prev = lv;
    }
  }
  const groups = cov.map(r => {
    const tr = trees[r.filename] || {};
    const L = r.levels || {};
    return {
      g: r.group, cat: r.category, prok: r.is_prokaryotic, tree: !r.is_unrepresented,
      described: r.described_species, described_text: r.described_text, described_source: r.described_source,
      levels: Object.fromEntries(LEVELS.filter(k => L[k]).map(k => [k, { n: L[k].n, n_text: L[k].n_text, pct: L[k].pct_text }])),
      dated: !!r.dated, est_low: r.estimated_low, est_high: r.estimated_high, est_total: r.estimated_total,
      est_kind: r.estimate_kind || null,
      file: r.filename, tips: r.tips, year: tr.year, study: tr.study, doi: tr.doi,
      in_package: !!tr.in_package, names_file: tr.names_file, dots: dots[r.group],
      withheld: tr.withheld || null,      // not redistributed: source link, never fetched
    };
  });
  const euk = cov.filter(r => !r.is_prokaryotic);
  return {
    scale: SCALE, paper: d.paper || {}, generated_from: d.generated_from || [],
    described_euk: euk.reduce((a, r) => a + r.described_species, 0),
    described_euk_text: d.summary.described_eukaryotic_text,
    described_prok_text: d.summary.described_prokaryotic_text,
    levels: d.levels, groups, tree: t.tree, tree_source: t.source, bins: t.bins,
  };
}

loadData().then(start).catch(e => {
  document.getElementById("titulo").textContent = "The data could not be loaded (" + e.message + ").";
});

function start(DATA) {
// ---------------------------------------------------------------- data
const G = Object.fromEntries(DATA.groups.map(g => [g.g, g]));
const order = [];                  // groups in tree order (leaves of the backbone)
(function walk(n) { if (n.g) order.push(n.g); (n.c || []).forEach(walk); })(DATA.tree);
const euk = DATA.groups.filter(g => !g.prok);

// Every dot of the ring: its group and its rank inside the group.
const ringDots = [];
for (const name of order) {
  const g = G[name];
  for (let k = 0; k < g.dots.all; k++) ringDots.push({ g: name, k });
}

// Estimated totals, summed over the eukaryotic groups (the release's values per group).
const EST = {
  low: euk.reduce((a, g) => a + (g.est_low ?? g.described), 0),
  total: euk.reduce((a, g) => a + (g.est_total ?? g.described), 0),
  high: euk.reduce((a, g) => a + (g.est_high ?? g.described), 0),
  published: euk.filter(g => g.est_kind && g.est_kind !== "none-published").length,
};
EST.lowDots = Math.round(EST.low / DATA.scale);
const describedDots = euk.reduce((a, g) => a + g.dots.all, 0);
const namedDots = euk.reduce((a, g) => a + g.dots.named, 0);
// all canvas dots: the ring, the faint dots up to the low estimate, and room for a group's tree
const POOL = Math.max(TREE_MAX, ringDots.length + (EST.lowDots - describedDots));
const mio = (n) => (n / 1e6).toFixed(n < 1e7 ? 2 : 1).replace(/\.?0+$/, "") + " million";

// ---------------------------------------------------------------- steps
const lv = Object.fromEntries(DATA.levels.map(l => [l.key, l]));
const pctOf = (n) => (100 * n / DATA.described_euk);
const STEPS = [
  { key: "campo", corto: "All",
    titulo: "Every described species of eukaryote",
    cifra: () => DATA.described_euk_text, sub: () => "described species",
    frase: () => `Each dot stands for ${fmt.format(DATA.scale)} described species: ${fmt.format(euk.reduce((a, g) => a + g.dots.all, 0))} dots. The open dots are Bacteria and Archaea, ${DATA.described_prok_text} GTDB species clusters, counted apart from every eukaryotic number.` },
  { key: "arbol", corto: "Tree",
    titulo: `${DATA.groups.length} groups, in the order of the tree of life`,
    cifra: () => String(DATA.groups.length), sub: () => "groups",
    frase: () => `Each arc of the ring holds the dots of one group, as wide as its described species. The tree inside is the branching order of the groups in the Open Tree of Life synthesis. Click or tap a group.`,
    nota: () => `${Object.keys(DATA.bins).length} catch-all groups (${Object.keys(DATA.bins).join(", ")}) are not clades; each sits at the common ancestor of its members. Groups smaller than one dot are drawn as a ring.` },
  ...LEVELS.map((k, i) => ({
    key: k, corto: String(i + 1), nivel: k,
    titulo: lv[k].label,
    cifra: () => `<span class="pct">${lv[k].pct_text}%</span>`,
    sub: () => `${lv[k].n_text} of ${DATA.described_euk_text} described species · ${lv[k].n_groups} groups`,
    frase: () => ({
      named: "Lit: the described species named in the reference tree of their group, one species-level phylogeny per group.",
      molecular: "Lit: the named species that the tree places with molecular data, not by taxonomy alone.",
      dated: "Lit: the species placed with molecular data in a chronogram, a tree whose branch lengths are absolute time.",
      uncertainty: "Lit: the dated species for which the source deposits the uncertainty of the dates (posterior trees or node intervals).",
    }[k]),
  })),
  { key: "final", corto: "Total",
    titulo: "Described species against the estimated total",
    cifra: () => `${(EST.low / 1e6).toFixed(1)}–${(EST.high / 1e6).toFixed(1)}<span style="font-size:.5em"> million</span>`,
    sub: () => `estimated species of eukaryotes; ${DATA.described_euk_text} described, ${(100 * DATA.described_euk / EST.high).toFixed(0)}–${(100 * DATA.described_euk / EST.low).toFixed(0)}% of the estimate`,
    frase: () => `At the same scale, one dot for ${fmt.format(DATA.scale)} species. The dark disc is the described species; its brass centre, the ${lv.named.n_text} named in a reference tree (${(100 * lv.named.n / EST.high).toFixed(1)}–${(100 * lv.named.n / EST.low).toFixed(1)}% of the estimate). The faint dots fill out to the low estimate, ${mio(EST.low)}; the thin circle encloses the central estimates, ${mio(EST.total)}; the dotted circle, the high estimates, ${mio(EST.high)}.`,
    nota: () => `Each total is the sum over the ${euk.length} eukaryotic groups of the group's low, central or high estimate in the release. For ${EST.published} of the ${euk.length} groups a published source states an estimate of the total; for the others the estimate is the described count, and a bound no source states is the described count too. Bacteria and Archaea are left out, as in every eukaryotic number.` },
];

const state = { step: 0, group: null, treeView: null, hover: null };

// ---------------------------------------------------------------- canvas
const caja = document.getElementById("lienzo");
const canvas = document.createElement("canvas");
caja.appendChild(canvas);
const ctx = canvas.getContext("2d");
let W = 0, H = 0, geo = null;

function resize() {
  const dpr = Math.min(devicePixelRatio || 1, 2);
  W = caja.clientWidth; H = caja.clientHeight;
  canvas.width = W * dpr; canvas.height = H * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  geo = layoutRing();
  retarget(true);
}

// The ring: one arc per group, columns of K dots, in tree order.
function layoutRing() {
  const cx = W / 2, cy = H / 2;
  const narrow = W < 520;
  const Ro = Math.min(W, H) / 2 - (narrow ? 34 : 70);
  const open = 0.05 * 2 * Math.PI;                     // a gap at the top, for the start of the ring
  const gapCols = 0.7;
  let best = null;
  for (let K = 4; K <= 12; K++) {
    const cols = order.reduce((a, n) => a + Math.max(1, Math.ceil(G[n].dots.all / K)), 0) + gapCols * order.length;
    // pitch p: ring of K dots thick must fit the arc length at mid radius
    let p = 0.3 * Ro / K;
    for (let it = 0; it < 30; it++) {
      const Rm = Ro - K * p / 2;
      const pMax = (2 * Math.PI - open) * Rm / cols;
      p = Math.min(0.3 * Ro / K, pMax);
    }
    if (!best || p > best.p) best = { K, p, cols };
  }
  const { K, p, cols } = best;
  const Ri = Ro - K * p;
  const dth = (2 * Math.PI - open) / cols;
  const start = -Math.PI / 2 + open / 2;
  const arcs = {};
  let col = 0;
  const pos = new Map();
  for (const n of order) {
    const g = G[n];
    const c = Math.max(1, Math.ceil(g.dots.all / K));
    const a0 = start + col * dth, a1 = a0 + c * dth;
    arcs[n] = { a0, a1, mid: (a0 + a1) / 2, cols: c };
    // fill row by row, from the tree outward: the lit share is a band on the inside
    for (let k = 0; k < g.dots.all; k++) {
      const cc = k % c, rr = Math.floor(k / c);
      const a = a0 + (cc + 0.5) * dth;
      const r = Ri + (rr + 0.5) * p;
      pos.set(n + "#" + k, { x: cx + Math.cos(a) * r, y: cy + Math.sin(a) * r });
    }
    col += c + gapCols;
  }
  // Backbone: a radial cladogram inside the ring, tips at the middle of each arc.
  const Rt = Ri - Math.max(6, p * 1.5);
  const root = d3.hierarchy(DATA.tree, d => d.c);
  const depth = root.height;
  root.eachAfter(d => {
    d.ang = d.children ? (d3.min(d.children, c => c.ang) + d3.max(d.children, c => c.ang)) / 2 : arcs[d.data.g].mid;
  });
  root.each(d => { d.rad = d.children ? Rt * (d.depth / depth) * 0.97 : Rt; });
  return { cx, cy, Ro, Ri, Rt, K, p, dth, start, arcs, pos, root, dotR: Math.max(0.9, p * 0.36), narrow };
}

// Field: a sunflower disc; dots keep their angular order so the opening is a radial move.
function fieldTargets() {
  const n = ringDots.length;
  const R = geo.Ro * 0.92, gold = Math.PI * (3 - Math.sqrt(5));
  const pts = [];
  for (let i = 0; i < n; i++) {
    const r = R * Math.sqrt((i + 0.5) / n), a = i * gold;
    pts.push({ x: geo.cx + Math.cos(a) * r, y: geo.cy + Math.sin(a) * r, ang: norm(a - geo.start) });
  }
  pts.sort((u, v) => u.ang - v.ang);
  return pts;
}
const norm = (a) => ((a % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
const prokCount = () => DATA.groups.filter(g => g.prok).reduce((a, g) => a + g.dots.all, 0);

// Finale: nested discs at one density. Sunflower order from the centre out: named (brass),
// described (ink), then faint dots up to the low estimate; circles for the central and high sums.
function finaleLayout() {
  const narrow = geo.narrow;
  const Rhigh = Math.min(W, H) / 2 - (narrow ? 26 : 56);
  const area = Math.PI * Rhigh * Rhigh / (EST.high / DATA.scale);   // px² per dot
  const rad = (dotsN) => Math.sqrt(dotsN * area / Math.PI);
  const pitch = Math.sqrt(area);
  const gold = Math.PI * (3 - Math.sqrt(5));
  const pt = (i) => { const r = rad(i + 0.5); const a = i * gold; return { x: geo.cx + Math.cos(a) * r, y: geo.cy + Math.sin(a) * r }; };
  return { pt, Rhigh, Rlow: rad(EST.lowDots), Rtot: rad(EST.total / DATA.scale), Rdesc: rad(describedDots), dotR: Math.max(0.5, pitch * 0.36) };
}
function finaleTargets() {
  const F = finaleLayout();
  geo.finale = F;
  const t = new Array(POOL);
  let lit = 0, ink = namedDots;
  ringDots.forEach((d, i) => {
    const g = G[d.g];
    if (g.prok) { t[i] = { ...geo.pos.get(d.g + "#" + d.k), r: geo.dotR, a: 0, rol: 2 }; return; }
    const isNamed = d.k < g.dots.named;
    const q = F.pt(isNamed ? lit++ : ink++);
    t[i] = { x: q.x, y: q.y, r: F.dotR, a: 0.9, rol: isNamed ? 1 : 0 };
  });
  for (let i = ringDots.length, j = describedDots; i < POOL && j < EST.lowDots; i++, j++) {
    const q = F.pt(j);
    t[i] = { x: q.x, y: q.y, r: F.dotR, a: 0.55, rol: 3 };
  }
  return t;
}

// ---------------------------------------------------------------- dots
// A dot: x, y, r, a (alpha), rol (0 ink, 1 brass, 2 open ink)
const dots = Array.from({ length: POOL }, () => ({ x: 0, y: 0, r: 0, a: 0, rol: 0 }));
let targets = dots.map(d => ({ ...d }));
let from = null, t0 = 0;
// ink → brass in 21 steps, for the colour fade
const MIX = Array.from({ length: 21 }, (_, k) => {
  const t = k / 20, a = [31, 35, 40], b = [176, 141, 79];
  return `rgb(${a.map((v, j) => Math.round(v + (b[j] - v) * t)).join(",")})`;
});

function litCount(g) {
  const s = STEPS[state.step];
  if (g.prok) return 0;
  if (s.nivel) return g.dots[s.nivel];
  return 0;
}

function ringTargets() {
  const s = STEPS[state.step];
  const out = [];
  const levelStep = !!s.nivel;
  for (const d of ringDots) {
    const g = G[d.g];
    const p = geo.pos.get(d.g + "#" + d.k);
    const lit = d.k < litCount(g);
    let rol = g.prok ? 2 : (lit ? 1 : 0);
    let a = g.prok ? 0.55 : (levelStep ? (lit ? 1 : 0.2) : 0.7);
    if (state.hover && state.hover !== d.g) a *= 0.45;
    out.push({ x: p.x, y: p.y, r: geo.dotR * (lit ? 1.15 : 1), a, rol });
  }
  return out;
}

function retarget(instant) {
  if (!geo) return;
  let t;
  if (state.treeView) {
    // the tips leave from the group's own arc: its ring dots go first
    const own = [], rest = [];
    ringDots.forEach((d, i) => (d.g === state.group ? own : rest).push(i));
    for (let i = ringDots.length; i < POOL; i++) rest.push(i);
    const idx = own.concat(rest);
    t = new Array(POOL);
    state.treeView.targets.forEach((m, j) => { t[idx[j]] = m; });
  }
  else if (STEPS[state.step].key === "final") {
    t = finaleTargets();
    // the faint dots appear where they stand: they only fade in
    for (let i = ringDots.length; i < POOL; i++) if (t[i] && dots[i].a < 0.02) { dots[i].x = t[i].x; dots[i].y = t[i].y; dots[i].r = t[i].r; dots[i].rol = 3; }
  }
  else if (STEPS[state.step].key === "campo") {
    const f = fieldTargets();
    t = ringTargets().map((d, i) => ({ ...d, x: f[i].x, y: f[i].y, a: d.rol === 2 ? 0.55 : 0.75, rol: d.rol === 2 ? 2 : 0 }));
  } else t = ringTargets();
  const full = [];
  for (let i = 0; i < POOL; i++) {
    full.push(t[i] || { x: dots[i].x || geo.cx, y: dots[i].y || geo.cy, r: dots[i].r || geo.dotR, a: 0, rol: dots[i].rol });
  }
  if (instant || reduce) { targets = full; for (let i = 0; i < POOL; i++) { Object.assign(dots[i], full[i]); dots[i].mix = full[i].rol === 1 ? 1 : 0; } from = null; return; }
  from = dots.map(d => ({ ...d }));
  targets = full;
  t0 = performance.now();
}

// ---------------------------------------------------------------- draw
let linesAlpha = 0, finaleAlpha = 0;
function drawFinale(alpha) {
  const F = geo.finale, { cx, cy } = geo;
  ctx.save();
  ctx.strokeStyle = "#1f2328"; ctx.fillStyle = "#1f2328";
  ctx.font = `${geo.narrow ? 11 : 13}px Libertinus, Georgia, serif`;
  ctx.textBaseline = "middle";
  const ring = (r, dash, w) => { ctx.setLineDash(dash); ctx.lineWidth = w; ctx.beginPath(); ctx.arc(cx, cy, r, 0, 2 * Math.PI); ctx.stroke(); };
  ctx.globalAlpha = alpha * 0.8;
  ring(F.Rtot, [], 0.7);
  ring(F.Rhigh, [1, 4], 1);
  ctx.setLineDash([]);
  // labels on a line to the upper right, each at its circle
  const ang = -Math.PI / 4.2;
  const lab = (r, text, dy = 0) => {
    const x = cx + Math.cos(ang) * r, y = cy + Math.sin(ang) * r;
    ctx.globalAlpha = alpha; ctx.beginPath(); ctx.arc(x, y, 1.8, 0, 2 * Math.PI); ctx.fill();
    ctx.textAlign = "left"; ctx.fillText(text, x + 6, y + dy);
  };
  const short = (n) => (n / 1e6).toFixed(1) + " M";
  lab(F.Rdesc, `${short(DATA.described_euk)} described`);
  lab(F.Rlow, `${short(EST.low)} low`);
  lab(F.Rtot, `${short(EST.total)} central`);
  lab(F.Rhigh, `${short(EST.high)} high`);
  ctx.restore();
}
function frame() {
  requestAnimationFrame(frame);
  if (!geo) return;
  const now = performance.now();
  const dt = (now - t0) / 1000;
  if (from) {
    for (let i = 0; i < POOL; i++) {
      const f = from[i], m = targets[i], d = dots[i];
      // a dot that keeps its place only changes colour: a short fade, no arc
      const moves = Math.abs(m.x - f.x) + Math.abs(m.y - f.y) > 1;
      const p = moves
        ? ease(Math.min(1, Math.max(0, (dt - hash(i + 99) * 0.45) / 0.9)))
        : ease(Math.min(1, Math.max(0, (dt - hash(i + 99) * 0.12) / 0.5)));
      const arco = moves ? Math.sin(p * Math.PI) * 26 * (hash(i + 3) - 0.5) : 0;
      d.x = f.x + (m.x - f.x) * p + arco; d.y = f.y + (m.y - f.y) * p - Math.abs(arco);
      d.r = f.r + (m.r - f.r) * p; d.a = f.a + (m.a - f.a) * p;
      d.mix = (f.mix ?? (f.rol === 1 ? 1 : 0)) * (1 - p) + (m.rol === 1 ? 1 : 0) * p;
      d.rol = m.rol === 2 || f.rol === 2 ? (p < 0.5 ? f.rol : m.rol) : m.rol;
    }
    if (dt > 1.5) { from = null; for (let i = 0; i < POOL; i++) dots[i].mix = targets[i].rol === 1 ? 1 : 0; }
  }
  const wantLines = state.treeView ? 0 : (["campo", "final"].includes(STEPS[state.step].key) ? 0 : 1);
  linesAlpha += (wantLines - linesAlpha) * (reduce ? 1 : 0.08);
  ctx.clearRect(0, 0, W, H);
  drawBackbone(linesAlpha);
  if (state.treeView) drawGroupTree(state.treeView, Math.min(1, dt / 1.2));
  const col = ["#1f2328", "#b08d4f", "#1f2328"];
  for (const d of dots) {
    if (d.a <= 0.01 || d.r <= 0) continue;
    ctx.globalAlpha = d.a;
    ctx.beginPath();
    ctx.arc(d.x, d.y, d.r, 0, Math.PI * 2);
    if (d.rol === 2) { ctx.strokeStyle = col[2]; ctx.lineWidth = 0.7; ctx.stroke(); }
    else { ctx.fillStyle = d.rol === 3 ? "#a3a8b0" : MIX[Math.round((d.mix ?? (d.rol === 1 ? 1 : 0)) * 20)]; ctx.fill(); }
  }
  ctx.globalAlpha = 1;
  if (linesAlpha > 0.02) drawLabels(linesAlpha);
  finaleAlpha += ((STEPS[state.step].key === "final" && !state.treeView ? 1 : 0) - finaleAlpha) * (reduce ? 1 : 0.05);
  if (finaleAlpha > 0.02 && geo.finale) drawFinale(finaleAlpha);
}

function drawBackbone(alpha) {
  if (alpha < 0.02) return;
  const { cx, cy, root } = geo;
  ctx.save();
  ctx.globalAlpha = alpha * 0.55;
  ctx.strokeStyle = "#1f2328";
  ctx.lineWidth = 0.8;
  ctx.beginPath();
  root.each(d => {
    if (d.children) {
      const a0 = d3.min(d.children, c => c.ang), a1 = d3.max(d.children, c => c.ang);
      ctx.moveTo(cx + Math.cos(a0) * d.rad, cy + Math.sin(a0) * d.rad);
      ctx.arc(cx, cy, d.rad, a0, a1);
      for (const c of d.children) {
        ctx.moveTo(cx + Math.cos(c.ang) * d.rad, cy + Math.sin(c.ang) * d.rad);
        ctx.lineTo(cx + Math.cos(c.ang) * c.rad, cy + Math.sin(c.ang) * c.rad);
      }
    }
  });
  ctx.stroke();
  // groups smaller than one dot: a ring at the tip
  for (const n of order) {
    if (G[n].dots.all > 0) continue;
    const a = geo.arcs[n].mid, r = (geo.Ri + geo.Ro) / 2;
    ctx.globalAlpha = alpha * (state.hover && state.hover !== n ? 0.35 : 0.9);
    ctx.beginPath(); ctx.arc(cx + Math.cos(a) * geo.Ri + Math.cos(a) * geo.p * 0.6, cy + Math.sin(a) * geo.Ri + Math.sin(a) * geo.p * 0.6, Math.max(1.6, geo.dotR * 1.2), 0, Math.PI * 2);
    ctx.stroke();
  }
  ctx.restore();
}

function drawLabels(alpha) {
  const { cx, cy, Ro, arcs, narrow } = geo;
  const minCols = narrow ? 14 : 6;
  ctx.save();
  ctx.font = `${narrow ? 11 : 13}px Libertinus, Georgia, serif`;
  ctx.fillStyle = "#1f2328";
  for (const n of order) {
    const a = arcs[n];
    const show = a.cols >= minCols || state.hover === n || state.group === n;
    if (!show) continue;
    ctx.globalAlpha = alpha * (state.hover && state.hover !== n ? 0.35 : 1);
    const r = Ro + 8;
    const x = cx + Math.cos(a.mid) * r, y = cy + Math.sin(a.mid) * r;
    const right = Math.cos(a.mid) >= 0;
    ctx.textAlign = right ? "left" : "right";
    ctx.textBaseline = "middle";
    ctx.fillText(n, x, y);
  }
  ctx.restore();
}

// ---------------------------------------------------------------- a group's own tree
function parseNewick(s) {
  let i = 0;
  s = s.replace(/\s/g, "").replace(/;$/, "");
  function node() {
    const n = { name: "", len: 0, c: null };
    if (s[i] === "(") {
      i++; n.c = [];
      for (;;) { n.c.push(node()); if (s[i] === ",") { i++; continue; } if (s[i] === ")") { i++; break; } throw new Error("Newick parse error at " + i); }
    }
    let j = i; while (j < s.length && !",():;".includes(s[j])) j++;
    n.name = s.slice(i, j); i = j;
    if (s[i] === ":") { i++; j = i; while (j < s.length && !",();".includes(s[j])) j++; n.len = parseFloat(s.slice(i, j)) || 0; i = j; }
    return n;
  }
  return node();
}
function leavesOf(n, out = []) { if (!n.c) out.push(n); else n.c.forEach(c => leavesOf(c, out)); return out; }
function prune(n, keep) {
  if (!n.c) return keep.has(n) ? n : null;
  const kids = n.c.map(c => prune(c, keep)).filter(Boolean);
  if (!kids.length) return null;
  if (kids.length === 1) { const k = kids[0]; return { ...k, len: (k.len || 0) + (n.len || 0) }; }
  return { ...n, c: kids };
}

async function openGroup(name) {
  const g = G[name];
  state.group = name;
  showCard(name);
  if (!g.file || g.withheld) { state.treeView = null; retarget(); return; }
  const status = document.getElementById("estado");
  if (status) status.textContent = `Loading the tree of ${name}…`;
  try {
    const [txt, names] = await Promise.all([
      fetch(TREES_BASE + g.file).then(r => { if (!r.ok) throw new Error("HTTP " + r.status); return r.text(); }),
      fetch(NAMES_BASE + (g.names_file || g.file) + ".json").then(r => r.ok ? r.json() : {}).catch(() => ({})),
    ]);
    if (state.group !== name) return;
    let tree = parseNewick(txt);
    const all = leavesOf(tree);
    let shown = all.length;
    if (all.length > TREE_MAX) {
      const keep = new Set(d3.shuffle(all.slice()).slice(0, TREE_MAX));
      tree = prune(tree, keep);
      shown = TREE_MAX;
    }
    state.treeView = buildTreeView(tree, g, names);
    if (status) status.textContent = shown < all.length
      ? `Here one dot is one species: ${fmt.format(shown)} of the ${fmt.format(all.length)} tips, drawn at random.`
      : `Here one dot is one species: all ${fmt.format(all.length)} tips of the tree.`;
    retarget(); render();
  } catch (e) {
    if (status) status.textContent = `The tree could not be loaded (${e.message}).`;
  }
}

function buildTreeView(tree, g, names) {
  const { cx, cy } = geo;
  const R = Math.min(W, H) / 2 - (geo.narrow ? 16 : 40);
  const root = d3.hierarchy(tree, d => d.c);
  d3.cluster().size([2 * Math.PI, R]).separation(() => 1)(root);
  let maxD = 0;
  root.each(d => { d.dist = (d.parent ? d.parent.dist : 0) + Math.max(0, d.data.len || 0); maxD = Math.max(maxD, d.dist); });
  const timed = g.dated && maxD > 0;
  if (timed) root.each(d => { d.y = R * d.dist / maxD; });
  const leaves = root.leaves();
  // ring position of the group's dots: tips leave from there
  const targets = leaves.map(l => ({
    x: cx + Math.cos(l.x - Math.PI / 2) * l.y, y: cy + Math.sin(l.x - Math.PI / 2) * l.y,
    r: Math.max(0.8, Math.min(2.6, 900 / leaves.length)), a: 0.95, rol: 1,
  }));
  return { root, leaves, targets, names, R, timed, maxD };
}

function drawGroupTree(v, p) {
  const { cx, cy } = geo;
  ctx.save();
  ctx.globalAlpha = 0.5 * p;
  ctx.strokeStyle = "#1f2328";
  ctx.lineWidth = v.leaves.length > 2000 ? 0.35 : 0.6;
  ctx.beginPath();
  v.root.each(d => {
    if (!d.children) return;
    const a0 = d3.min(d.children, c => c.x) - Math.PI / 2, a1 = d3.max(d.children, c => c.x) - Math.PI / 2;
    ctx.moveTo(cx + Math.cos(a0) * d.y, cy + Math.sin(a0) * d.y);
    ctx.arc(cx, cy, d.y, a0, a1);
    for (const c of d.children) {
      const a = c.x - Math.PI / 2;
      ctx.moveTo(cx + Math.cos(a) * d.y, cy + Math.sin(a) * d.y);
      ctx.lineTo(cx + Math.cos(a) * c.y, cy + Math.sin(a) * c.y);
    }
  });
  ctx.stroke();
  if (v.timed && v.maxD > 30) {                       // time rings, in Ma
    ctx.globalAlpha = 0.6 * p;
    ctx.strokeStyle = "#b08d4f"; ctx.fillStyle = "#8b9098";
    ctx.setLineDash([1, 3]);
    ctx.font = "11px Libertinus, Georgia, serif"; ctx.textAlign = "left";
    const step = [10, 20, 25, 50, 100, 200, 250, 500, 1000].find(s => v.maxD / s <= 5) || 1000;
    for (let age = step; age < v.maxD; age += step) {
      const r = v.R * (1 - age / v.maxD);
      ctx.beginPath(); ctx.arc(cx, cy, r, 0, 2 * Math.PI); ctx.stroke();
      ctx.fillText(`${age} Ma`, cx + 3, cy - r - 3);
    }
    ctx.setLineDash([]);
  }
  ctx.restore();
}

// ---------------------------------------------------------------- text
function render() {
  const s = STEPS[state.step];
  document.getElementById("escala").innerHTML = state.treeView
    ? "one dot = <b>one species</b>"
    : `one dot = <b>${fmt.format(DATA.scale)}</b> ${s.key === "final" ? "" : "described "}species`;
  document.getElementById("paso-n").textContent = `${state.step + 1} of ${STEPS.length}`;
  document.getElementById("titulo").textContent = s.titulo;
  document.getElementById("cifra").innerHTML = s.cifra();
  const sub = s.sub();
  document.getElementById("cifra-sub").textContent = sub;
  document.getElementById("frase").textContent = s.frase();
  document.getElementById("nota").textContent = s.nota ? s.nota() : "";
  const lit = !!s.nivel;
  document.getElementById("clave").innerHTML = s.key === "final"
    ? `<span><i style="background:#b08d4f"></i>named in a reference tree</span><span><i style="background:#1f2328"></i>described</span><span><i style="background:#a3a8b0"></i>estimated, up to the low sum</span>`
    : s.key === "campo" ? "" :
    (lit ? `<span><i style="background:#b08d4f"></i>${lv[s.nivel].label.toLowerCase()}</span><span><i style="background:#1f2328;opacity:.25"></i>not</span>` :
           `<span><i style="background:#1f2328"></i>described species</span>`) +
    `<span><i style="border:1px solid #1f2328"></i>Bacteria, Archaea (apart)</span>`;
  document.querySelectorAll("#nav .pt").forEach((b, i) => b.classList.toggle("on", i === state.step));
}

function showCard(name) {
  const c = document.getElementById("ficha");
  if (!name) { c.classList.remove("on"); c.innerHTML = ""; return; }
  const g = G[name];
  const L = g.levels || {};
  const rows = LEVELS.filter(k => L[k]).map(k =>
    `<tr><td>${lv[k].label}</td><td class="n">${L[k].n_text}</td><td class="p">${L[k].pct}%</td></tr>`).join("");
  const est = (g.est_low != null && g.est_high != null && g.est_high > g.est_low)
    ? `<tr><td>Estimated total</td><td class="n">${fmt.format(g.est_low)}–${fmt.format(g.est_high)}</td><td></td></tr>` : "";
  const ref = g.file
    ? `Reference tree: ${g.study || g.file}${g.year && !(g.study || "").includes(g.year) ? " " + g.year : ""}${g.doi ? ` · <a href="https://doi.org/${g.doi}" target="_blank" rel="noopener">doi</a>` : ""}${g.tips ? ` · ${fmt.format(g.tips)} tips` : ""}${g.dated ? " · chronogram" : ""}`
    : (g.prok ? "" : "No species-level reference tree for this group.");
  const dotsTxt = g.dots.all ? `${fmt.format(g.dots.all)} dot${g.dots.all > 1 ? "s" : ""}` : "smaller than one dot";
  c.innerHTML = `
    <div class="cat">${g.cat}${g.prok ? " · GTDB species clusters, counted apart" : ""}</div>
    <h3>${name}</h3>
    <table>
      <tr><td>Described species</td><td class="n">${g.described_text}</td><td class="p" style="color:#8b9098">${dotsTxt}</td></tr>
      ${rows}${est}
    </table>
    <div class="ref">${ref}</div>
    ${g.file && g.in_package ? `<div class="ref dl"><a href="${TREES_BASE + g.file}" download>Download the tree (.nwk)</a> · in R:</div><div class="r">tree &lt;- phyloatlas::load_atlas_tree("${g.file.replace(/\.nwk$/, "")}")</div>` : ""}
    ${g.file && g.withheld ? `<div class="ref dl">Tree file ${g.withheld.text}. Source: ${srcLink(g.withheld.source_url)}</div>` : ""}
    <div class="ref" id="estado"></div>
    ${state.treeView || state.group === name && g.file ? `<button type="button" id="cerrar">Back to the ring</button>` : ""}`;
  c.classList.add("on");
  const b = document.getElementById("cerrar");
  if (b) b.onclick = closeGroup;
}

// the source of a withheld tree: the first URL of the register's text is the link
function srcLink(s) {
  const esc = (x) => String(x).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/"/g, "&quot;");
  const m = /https?:\/\/[^\s;)]+/.exec(s || "");
  return m ? `<a href="${esc(m[0])}" target="_blank" rel="noopener">${esc(m[0])}</a>` : esc(s || "");
}

function closeGroup() {
  state.group = null; state.treeView = null; showCard(null); retarget(); render();
}

function go(i) {
  state.step = Math.max(0, Math.min(STEPS.length - 1, i));
  state.treeView = null; state.group = null; showCard(null);
  retarget(); render();
}

// ---------------------------------------------------------------- pointer
const tip = document.getElementById("tip");
function groupAt(px, py) {
  if (!geo || ["campo", "final"].includes(STEPS[state.step].key)) return null;
  const dx = px - geo.cx, dy = py - geo.cy, r = Math.hypot(dx, dy);
  if (r < geo.Ri - 2 * geo.p || r > geo.Ro + 40) return null;
  const a = Math.atan2(dy, dx);
  let best = null, bd = Infinity;
  for (const n of order) {
    const arc = geo.arcs[n];
    const d = Math.abs(norm(a - arc.mid + Math.PI) - Math.PI) - (arc.a1 - arc.a0) / 2;
    if (d < bd) { bd = d; best = n; }
  }
  return bd <= geo.dth * 1.2 ? best : null;
}
function tipAt(px, py) {
  const v = state.treeView; if (!v) return null;
  let best = null, bd = 64;
  v.targets.forEach((t, i) => { const d = (t.x - px) ** 2 + (t.y - py) ** 2; if (d < bd) { bd = d; best = i; } });
  return best == null ? null : v.leaves[best];
}
canvas.addEventListener("pointermove", (e) => {
  const b = canvas.getBoundingClientRect(), px = e.clientX - b.left, py = e.clientY - b.top;
  if (state.treeView) {
    const l = tipAt(px, py);
    if (l) { tip.style.display = "block"; tip.textContent = (state.treeView.names[l.data.name] || l.data.name).replace(/_/g, " "); tip.style.left = e.clientX + 12 + "px"; tip.style.top = e.clientY + 12 + "px"; }
    else tip.style.display = "none";
    return;
  }
  const n = groupAt(px, py);
  canvas.style.cursor = n ? "pointer" : "default";
  if (n !== state.hover) { state.hover = n; retarget(); }
  if (n) { tip.style.display = "block"; tip.textContent = n; tip.style.left = e.clientX + 12 + "px"; tip.style.top = e.clientY + 12 + "px"; }
  else tip.style.display = "none";
});
canvas.addEventListener("pointerleave", () => { tip.style.display = "none"; if (state.hover) { state.hover = null; retarget(); } });
canvas.addEventListener("click", (e) => {
  if (state.treeView) return;
  const b = canvas.getBoundingClientRect();
  const n = groupAt(e.clientX - b.left, e.clientY - b.top);
  if (!n) { if (state.group) closeGroup(); return; }
  state.hover = null;
  openGroup(n);
});
document.addEventListener("keydown", (e) => {
  if (e.key === "ArrowRight") go(state.step + 1);
  else if (e.key === "ArrowLeft") go(state.step - 1);
  else if (e.key === "Escape" && state.group) closeGroup();
});

// ---------------------------------------------------------------- start
const nav = document.getElementById("nav");
nav.innerHTML = `<button type="button" id="prev">←</button>` +
  STEPS.map((s, i) => `<button type="button" class="pt" data-i="${i}" title="${s.titulo}">${i + 1}</button>`).join("") +
  `<button type="button" id="next">→</button>`;
nav.querySelectorAll(".pt").forEach(b => b.onclick = () => go(+b.dataset.i));
document.getElementById("prev").onclick = () => go(state.step - 1);
document.getElementById("next").onclick = () => go(state.step + 1);
document.getElementById("pie").innerHTML =
  `Companion to Richter (2026), <a href="https://doi.org/${DATA.paper.doi}" target="_blank" rel="noopener"><em>${DATA.paper.title}</em></a>. ` +
  `Every number is read from the release's <code>data.json</code> (${DATA.generated_from.slice(0, 2).join(", ")}). ` +
  `Dots of each level are allocated to the groups by largest remainder, so every level's dots sum to its total ÷ ${DATA.scale}. ` +
  `Order of the groups: ${DATA.tree_source}. ` +
  `<a href="https://github.com/franciscorichter/phylo-species-atlas" target="_blank" rel="noopener">Data and code</a> · ` +
  `<a href="https://doi.org/10.5281/zenodo.20127157" target="_blank" rel="noopener">Zenodo</a> · ` +
  `CC-BY 4.0; cite the original studies for the trees.`;
new ResizeObserver(resize).observe(caja);
document.fonts.ready.then(() => { resize(); render(); });
resize(); render();
requestAnimationFrame(frame);
}
