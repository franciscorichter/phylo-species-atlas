// Phylo-Species Atlas browser — vanilla JS, no build step.
// Loads data.json, renders summary + plots, and lazy-loads tree files.

const STATE = {
  data: null,
  filtered: [],
  selected: null,
  layout: "linear", // or "radial"
  showLabels: false,
  currentTree: null,
};

const TREES_BASE = "../standardized/trees/"; // served from repo root on Pages
const TIP_RENDER_LIMIT = 5000; // above this, we sample for rendering

const fmt = new Intl.NumberFormat("en-US");

document.addEventListener("DOMContentLoaded", async () => {
  try {
    const res = await fetch("data.json");
    if (!res.ok) throw new Error(`data.json: ${res.status}`);
    STATE.data = await res.json();
  } catch (e) {
    document.body.innerHTML = `<p style="padding:40px;color:#c00">Failed to load data.json — ${e.message}</p>`;
    return;
  }

  renderSummary();
  renderCoveragePlot();
  renderSizePlot();
  wireFilters();
  applyFilters();
});

function renderSummary() {
  const s = STATE.data.summary;
  const items = [
    ["Trees", s.n_trees],
    ["Datasets", s.n_datasets],
    ["Dated", s.n_dated],
    ["Undated", s.n_undated],
    ["Total tips", s.total_tips],
  ];
  const root = document.getElementById("summary-stats");
  root.innerHTML = items.map(([label, value]) =>
    `<div class="stat"><span class="value">${fmt.format(value)}</span><span class="label">${label}</span></div>`
  ).join("");
}

function renderCoveragePlot() {
  const rows = STATE.data.coverage;
  if (!rows.length) return;

  const labels = rows.map(r => r.group);
  const values = rows.map(r => r.coverage_pct);
  const described = rows.map(r => r.described_species);
  const tips = rows.map(r => r.tips);
  const colors = rows.map(r => r.dated ? "#2a6fbf" : "#c97a2a");
  const customdata = rows.map(r => [r.study, r.year || "", r.tips, r.described_species]);

  const trace = {
    type: "bar",
    orientation: "h",
    x: values,
    y: labels,
    marker: { color: colors },
    customdata,
    hovertemplate:
      "<b>%{y}</b><br>" +
      "Coverage: %{x:.1f}%<br>" +
      "Tips: %{customdata[2]:,}<br>" +
      "Described: %{customdata[3]:,}<br>" +
      "%{customdata[0]} (%{customdata[1]})<extra></extra>",
  };

  const layout = {
    margin: { l: 130, r: 24, t: 12, b: 36 },
    xaxis: { title: "Coverage (%)", range: [0, 105], gridcolor: "#eee" },
    yaxis: { automargin: true, tickfont: { size: 11 } },
    height: Math.max(420, labels.length * 18),
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    showlegend: false,
  };

  Plotly.newPlot("coverage-plot", [trace], layout, { displayModeBar: false, responsive: true });

  document.getElementById("coverage-plot").on("plotly_click", (ev) => {
    const point = ev.points && ev.points[0];
    if (!point) return;
    const group = point.y;
    const target = STATE.data.trees.find(t => t.group === group);
    if (target) selectTree(target.filename, true);
  });
}

function renderSizePlot() {
  const trees = STATE.data.trees;
  const tips = trees.map(t => t.ntips || 0).filter(n => n > 0);

  const trace = {
    type: "histogram",
    x: tips,
    xbins: { start: 0, end: Math.log10(Math.max(...tips)) + 0.5, size: 0.25 },
    marker: { color: "#2a6fbf" },
    hovertemplate: "%{y} trees<br>%{x} tips<extra></extra>",
  };

  // Transform x to log scale
  trace.x = tips.map(n => Math.log10(n));
  trace.hovertemplate = "%{y} trees<extra></extra>";

  const layout = {
    margin: { l: 40, r: 16, t: 8, b: 36 },
    xaxis: {
      title: "Tips (log10)",
      tickvals: [1, 2, 3, 4, 5, 6],
      ticktext: ["10", "100", "1K", "10K", "100K", "1M"],
      gridcolor: "#eee",
    },
    yaxis: { title: "Trees", gridcolor: "#eee" },
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    bargap: 0.05,
  };

  Plotly.newPlot("size-plot", [trace], layout, { displayModeBar: false, responsive: true });
}

function wireFilters() {
  document.getElementById("tree-search").addEventListener("input", applyFilters);
  document.getElementById("filter-dated").addEventListener("change", applyFilters);
  document.getElementById("sort-by").addEventListener("change", applyFilters);
  document.getElementById("show-labels").addEventListener("change", (e) => {
    STATE.showLabels = e.target.checked;
    redrawTree();
  });
  document.getElementById("layout-toggle").addEventListener("click", () => {
    STATE.layout = STATE.layout === "linear" ? "radial" : "linear";
    document.getElementById("layout-toggle").textContent =
      STATE.layout === "linear" ? "Radial" : "Linear";
    redrawTree();
  });
  document.getElementById("reset-zoom").addEventListener("click", redrawTree);
}

function applyFilters() {
  const q = document.getElementById("tree-search").value.trim().toLowerCase();
  const datedOnly = document.getElementById("filter-dated").checked;
  const sort = document.getElementById("sort-by").value;

  let list = STATE.data.trees.slice();
  if (q) {
    list = list.filter(t =>
      (t.filename || "").toLowerCase().includes(q) ||
      (t.group || "").toLowerCase().includes(q) ||
      (t.study || "").toLowerCase().includes(q)
    );
  }
  if (datedOnly) list = list.filter(t => t.dated);

  list.sort((a, b) => {
    switch (sort) {
      case "ntips-asc": return (a.ntips || 0) - (b.ntips || 0);
      case "filename-asc": return a.filename.localeCompare(b.filename);
      case "year-desc": return (b.year || 0) - (a.year || 0);
      case "ntips-desc":
      default: return (b.ntips || 0) - (a.ntips || 0);
    }
  });

  STATE.filtered = list;
  renderTreeList();
}

function renderTreeList() {
  const ul = document.getElementById("tree-list");
  ul.innerHTML = STATE.filtered.map(t => {
    const cls = STATE.selected === t.filename ? "active" : "";
    const badge = t.dated ? "dated" : "undated";
    const name = t.filename.replace(/\.nwk$/, "");
    return `<li data-filename="${t.filename}" class="${cls}">
      <span class="tree-name"><span class="badge ${badge}"></span>${name}</span>
      <span class="tree-meta">${fmt.format(t.ntips || 0)}</span>
    </li>`;
  }).join("");
  ul.querySelectorAll("li").forEach(li => {
    li.addEventListener("click", () => selectTree(li.dataset.filename, false));
  });
}

async function selectTree(filename, scrollIntoView) {
  STATE.selected = filename;
  renderTreeList();
  if (scrollIntoView) {
    const li = document.querySelector(`#tree-list li[data-filename="${filename}"]`);
    if (li) li.scrollIntoView({ block: "nearest", behavior: "smooth" });
  }

  const t = STATE.data.trees.find(x => x.filename === filename);
  if (!t) return;

  document.getElementById("detail-empty").hidden = true;
  document.getElementById("detail").hidden = false;
  document.getElementById("detail-name").textContent = t.filename.replace(/\.nwk$/, "");

  const fields = [
    ["Group", t.group],
    ["Study", t.study],
    ["Year", t.year],
    ["Journal", t.journal],
    ["Tips", fmt.format(t.ntips || 0)],
    ["Dated", t.dated ? "yes" : "no"],
    ["Crown age", t.crown_ma ? `${t.crown_ma} Ma` : "—"],
    ["Described species", t.described_species ? fmt.format(t.described_species) : "—"],
    ["Coverage", t.coverage_pct != null ? `${t.coverage_pct.toFixed(1)}%` : "—"],
    ["DOI", t.doi ? `<a href="https://doi.org/${t.doi}" target="_blank" rel="noopener">${t.doi}</a>` : "—"],
    ["Source", t.data_source],
    ["File", `${(t.size_bytes / 1024).toFixed(1)} KB`],
  ];
  document.getElementById("detail-meta").innerHTML = fields
    .filter(([_, v]) => v !== null && v !== undefined && v !== "")
    .map(([k, v]) => `<div class="field"><span class="key">${k}</span><span class="val">${v}</span></div>`)
    .join("");

  const dl = document.getElementById("download-newick");
  dl.href = TREES_BASE + filename;
  dl.setAttribute("download", filename);

  await loadAndRenderTree(t);
}

async function loadAndRenderTree(t) {
  const status = document.getElementById("tree-status");
  const container = document.getElementById("tree-container");
  container.innerHTML = "";
  STATE.currentTree = null;

  status.classList.remove("warning");
  status.textContent = `Loading ${t.filename}…`;

  try {
    const res = await fetch(TREES_BASE + t.filename);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    let newick = await res.text();

    if ((t.ntips || 0) > TIP_RENDER_LIMIT) {
      const sampled = sampleNewick(newick, TIP_RENDER_LIMIT);
      newick = sampled.newick;
      status.classList.add("warning");
      status.textContent = `Tree has ${fmt.format(t.ntips)} tips — showing a random sample of ${fmt.format(sampled.kept)}. Download for the full tree.`;
    } else {
      status.textContent = `${fmt.format(t.ntips || 0)} tips — full tree shown.`;
    }

    STATE.currentTree = newick;
    drawTree(newick);
  } catch (e) {
    status.classList.add("warning");
    status.textContent = `Could not render: ${e.message}. Use the download link above.`;
  }
}

function redrawTree() {
  if (STATE.currentTree) drawTree(STATE.currentTree);
}

function drawTree(newick) {
  const container = document.getElementById("tree-container");
  container.innerHTML = "";

  const width = Math.max(640, container.clientWidth - 8);
  const tree = new phylotree.phylotree(newick);
  tree.render({
    container: "#tree-container",
    "show-scale": false,
    "align-tips": false,
    "show-labels": STATE.showLabels,
    "is-radial": STATE.layout === "radial",
    "left-right-spacing": "fit-to-size",
    "top-bottom-spacing": "fit-to-size",
    width,
    height: 520,
    "internal-names": false,
    selectable: false,
    "draw-size-bubbles": false,
  });
  tree.placenodes();
  tree.update();
}

// Cheap Newick subsampling: randomly drop tips until we hit `keep`.
// Re-parses the string, identifying leaf positions, and removes leaves
// (with surrounding comma/parens cleanup) until target count is reached.
function sampleNewick(newick, keep) {
  // Match leaf tokens: a label followed by optional :branchlength,
  // preceded by '(' or ',' — and NOT followed by '('.
  const tokens = [];
  const re = /([(,])\s*([A-Za-z0-9_.+\-]+)(:[-0-9.eE]+)?(?=[,)])/g;
  let m;
  while ((m = re.exec(newick)) !== null) {
    tokens.push({ start: m.index + m[1].length, end: m.index + m[0].length, raw: m[0].slice(1) });
  }
  if (tokens.length <= keep) return { newick, kept: tokens.length };

  // Randomly choose tips to drop
  const dropCount = tokens.length - keep;
  const indices = new Set();
  while (indices.size < dropCount) indices.add(Math.floor(Math.random() * tokens.length));

  // Build new string by removing chosen leaves (and their adjacent comma)
  // Sort indices descending to splice without shifting earlier offsets.
  const toRemove = [...indices].map(i => tokens[i]).sort((a, b) => b.start - a.start);
  let s = newick;
  for (const t of toRemove) {
    // Remove `<leaf>` plus a neighboring comma. If preceding char is ',', drop it.
    let start = t.start;
    let end = t.end;
    if (s[start - 1] === ",") start -= 1;
    else if (s[end] === ",") end += 1;
    s = s.slice(0, start) + s.slice(end);
  }
  // Collapse degenerate clades: "(X)" -> "X", "()" -> ""
  let prev;
  do {
    prev = s;
    s = s.replace(/\(([A-Za-z0-9_.+\-:]+(?::[-0-9.eE]+)?)\)/g, "$1");
    s = s.replace(/\(\)/g, "");
  } while (s !== prev);
  return { newick: s, kept: keep };
}
