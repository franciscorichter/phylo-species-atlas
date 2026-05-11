// Phylo-Species Atlas browser — vanilla JS, no build step.
// Loads data.json, renders summary + plots, and lazy-loads tree files.

const STATE = {
  data: null,
  filtered: [],
  selected: null,
  layout: "linear", // or "radial"
  showLabels: false,
  currentNewick: null,
};

const TREES_BASE = "../standardized/trees/"; // patched by CI for Pages deploy
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
  document.getElementById("summary-stats").innerHTML = items.map(([label, value]) =>
    `<div class="stat"><span class="value">${fmt.format(value)}</span><span class="label">${label}</span></div>`
  ).join("");
}

function renderCoveragePlot() {
  const rows = STATE.data.coverage;
  if (!rows.length) return;

  const labels = rows.map(r => r.group);
  const values = rows.map(r => r.coverage_pct);
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

  const chartHeight = Math.max(420, labels.length * 18 + 60);
  const container = document.getElementById("coverage-plot");
  container.style.height = chartHeight + "px";

  const layout = {
    margin: { l: 140, r: 24, t: 12, b: 40 },
    xaxis: { title: "Coverage (%)", range: [0, 105], gridcolor: "#eee" },
    yaxis: { automargin: false, tickfont: { size: 11 }, autorange: "reversed" },
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    showlegend: false,
  };

  Plotly.newPlot(container, [trace], layout, { displayModeBar: false, responsive: true });

  container.on("plotly_click", (ev) => {
    const point = ev.points && ev.points[0];
    if (!point) return;
    const target = STATE.data.trees.find(t => t.group === point.y);
    if (target) selectTree(target.filename, true);
  });
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
    ["File", t.size_bytes ? `${(t.size_bytes / 1024).toFixed(1)} KB` : "—"],
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
  STATE.currentNewick = null;

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

    STATE.currentNewick = newick;
    drawTree(newick);
  } catch (e) {
    status.classList.add("warning");
    status.textContent = `Could not render: ${e.message}. Use the download link above.`;
  }
}

function redrawTree() {
  if (STATE.currentNewick) drawTree(STATE.currentNewick);
}

// ---------- Newick parser ----------
// Produces {name, length, children}. Handles branch lengths, no support values.
function parseNewick(str) {
  let i = 0;
  function parseNode() {
    const node = { name: "", length: null, children: null };
    if (str[i] === "(") {
      i++;
      node.children = [];
      while (true) {
        node.children.push(parseNode());
        if (str[i] === ",") { i++; continue; }
        if (str[i] === ")") { i++; break; }
        throw new Error(`Newick parse error at ${i}: expected , or )`);
      }
    }
    let name = "";
    while (i < str.length && !",():;".includes(str[i])) {
      if (str[i] === ":") break;
      name += str[i++];
    }
    node.name = name.trim();
    if (str[i] === ":") {
      i++;
      let len = "";
      while (i < str.length && !",();".includes(str[i])) len += str[i++];
      node.length = parseFloat(len);
    }
    return node;
  }
  // Strip whitespace
  str = str.replace(/\s/g, "");
  const root = parseNode();
  return root;
}

// ---------- D3 tree renderer ----------
function drawTree(newick) {
  const container = document.getElementById("tree-container");
  container.innerHTML = "";

  let parsed;
  try { parsed = parseNewick(newick); }
  catch (e) {
    container.innerHTML = `<p style="padding:20px;color:#c00">Newick parse error: ${e.message}</p>`;
    return;
  }

  const root = d3.hierarchy(parsed, d => d.children);
  const nLeaves = root.leaves().length;

  const isRadial = STATE.layout === "radial";
  const w = Math.max(640, container.clientWidth - 4);

  if (isRadial) {
    drawRadial(root, w);
  } else {
    drawLinear(root, w, nLeaves);
  }
}

function drawLinear(root, width, nLeaves) {
  const labelWidth = STATE.showLabels ? 140 : 8;
  const margin = { top: 12, right: labelWidth, bottom: 12, left: 12 };
  const innerW = width - margin.left - margin.right;
  const tipSpacing = Math.max(1.2, Math.min(16, 600 / Math.max(40, nLeaves)));
  const innerH = Math.max(400, nLeaves * tipSpacing);

  const cluster = d3.cluster()
    .size([innerH, innerW])
    .separation(() => 1);
  cluster(root);

  const svg = d3.select("#tree-container").append("svg")
    .attr("viewBox", `0 0 ${width} ${innerH + margin.top + margin.bottom}`)
    .attr("preserveAspectRatio", "xMidYMid meet")
    .attr("width", "100%")
    .attr("height", innerH + margin.top + margin.bottom);

  const g = svg.append("g").attr("transform", `translate(${margin.left},${margin.top})`);

  // Step paths: horizontal then vertical (right-angle dendrogram)
  g.selectAll("path.branch")
    .data(root.descendants().slice(1))
    .join("path")
      .attr("class", "branch")
      .attr("d", d => `M${d.parent.y},${d.parent.x} V${d.x} H${d.y}`)
      .attr("fill", "none")
      .attr("stroke", "#1a1a1a")
      .attr("stroke-width", 0.8);

  if (STATE.showLabels) {
    g.selectAll("text.label")
      .data(root.leaves())
      .join("text")
        .attr("class", "label")
        .attr("x", d => d.y + 4)
        .attr("y", d => d.x)
        .attr("dy", "0.32em")
        .attr("font-size", Math.max(7, Math.min(11, tipSpacing * 0.9)))
        .attr("fill", "#1a1a1a")
        .text(d => d.data.name);
  }
}

function drawRadial(root, width) {
  const size = Math.min(width, 720);
  const radius = size / 2 - (STATE.showLabels ? 80 : 12);

  const cluster = d3.cluster()
    .size([2 * Math.PI, radius])
    .separation(() => 1);
  cluster(root);

  const svg = d3.select("#tree-container").append("svg")
    .attr("viewBox", `${-size / 2} ${-size / 2} ${size} ${size}`)
    .attr("preserveAspectRatio", "xMidYMid meet")
    .attr("width", "100%")
    .attr("height", size);

  const linkRadial = d3.linkRadial()
    .angle(d => d.x)
    .radius(d => d.y);

  svg.append("g")
    .selectAll("path")
    .data(root.links())
    .join("path")
      .attr("d", linkRadial)
      .attr("fill", "none")
      .attr("stroke", "#1a1a1a")
      .attr("stroke-width", 0.7);

  if (STATE.showLabels) {
    svg.append("g")
      .selectAll("text")
      .data(root.leaves())
      .join("text")
        .attr("transform", d => `rotate(${(d.x * 180 / Math.PI) - 90}) translate(${d.y + 4},0)${d.x >= Math.PI ? " rotate(180)" : ""}`)
        .attr("text-anchor", d => d.x < Math.PI ? "start" : "end")
        .attr("dy", "0.32em")
        .attr("font-size", 8)
        .attr("fill", "#1a1a1a")
        .text(d => d.data.name);
  }
}

// Cheap Newick subsampling: randomly drop tips until we hit `keep`.
function sampleNewick(newick, keep) {
  const tokens = [];
  const re = /([(,])\s*([A-Za-z0-9_.+\-]+)(:[-0-9.eE]+)?(?=[,)])/g;
  let m;
  while ((m = re.exec(newick)) !== null) {
    tokens.push({ start: m.index + m[1].length, end: m.index + m[0].length });
  }
  if (tokens.length <= keep) return { newick, kept: tokens.length };

  const dropCount = tokens.length - keep;
  const indices = new Set();
  while (indices.size < dropCount) indices.add(Math.floor(Math.random() * tokens.length));

  const toRemove = [...indices].map(i => tokens[i]).sort((a, b) => b.start - a.start);
  let s = newick;
  for (const t of toRemove) {
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
