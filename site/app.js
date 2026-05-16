// Phylo-Species Atlas browser — vanilla JS, no build step.
// Loads data.json, renders summary + plots, and lazy-loads tree files.

const STATE = {
  data: null,
  filtered: [],
  selected: null,
  layout: "linear", // or "radial"
  showLabels: false,
  currentNewick: null,
  coverageMode: "described", // or "estimated"
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
  wireCoverageToggle();
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

const FONT_STACK = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif";

function renderCoveragePlot() {
  const allRows = STATE.data.coverage;
  if (!allRows.length) return;
  const mode = STATE.coverageMode;
  const estimatedAvailable = (r) => r.coverage_pct_estimated != null;

  // In estimated mode, drop sub-clade datasets — their tips compared to the
  // parent partition's estimated total produce technically-correct but
  // misleading bars (e.g. "parrots: 3.7% of Birds estimate").
  const rows = mode === "estimated"
    ? allRows.filter(r => estimatedAvailable(r) && r.is_partition_anchor)
             .slice().sort((a, b) => b.coverage_pct_estimated - a.coverage_pct_estimated)
    : allRows;

  const labels = rows.map(r => r.group);
  const colors = rows.map(r => r.dated ? "#2a6fbf" : "#c97a2a");

  let values, valueText, customdata, hovertemplate, errorBars, xTitle;

  if (mode === "estimated") {
    values = rows.map(r => r.coverage_pct_estimated);
    valueText = rows.map(r => `${r.coverage_pct_estimated.toFixed(1)}%`);
    errorBars = {
      type: "data",
      symmetric: false,
      array: rows.map(r => Math.max(0, r.coverage_pct_estimated_high - r.coverage_pct_estimated)),
      arrayminus: rows.map(r => Math.max(0, r.coverage_pct_estimated - r.coverage_pct_estimated_low)),
      color: "#999",
      thickness: 1.2,
      width: 4,
    };
    customdata = rows.map(r => [
      r.study, r.year || "", r.tips, r.estimated_total, r.estimated_low, r.estimated_high,
      r.estimate_source || "—", r.estimate_confidence || "—",
      r.coverage_pct_estimated_low, r.coverage_pct_estimated_high,
    ]);
    hovertemplate =
      "<b>%{y}</b><br>" +
      "Coverage of estimate: %{x:.1f}% " +
      "<span style='color:#999'>(%{customdata[8]:.1f}–%{customdata[9]:.1f}%)</span><br>" +
      "Tips: %{customdata[2]:,}<br>" +
      "Estimated total: %{customdata[3]:,} (range %{customdata[4]:,}–%{customdata[5]:,})<br>" +
      "Estimate source: %{customdata[6]} · confidence %{customdata[7]}<br>" +
      "<i>%{customdata[0]}</i> (%{customdata[1]})<extra></extra>";
    xTitle = "Coverage of estimated true diversity (%) — error bars = low/high estimate range";
  } else {
    values = rows.map(r => r.coverage_pct);
    valueText = rows.map(r => `${r.coverage_pct.toFixed(1)}%`);
    errorBars = undefined;
    customdata = rows.map(r => [r.study, r.year || "", r.tips, r.described_species, r.described_source || "—"]);
    hovertemplate =
      "<b>%{y}</b><br>" +
      "Coverage: %{x:.1f}%<br>" +
      "Tips: %{customdata[2]:,}<br>" +
      "Described: %{customdata[3]:,} (source: %{customdata[4]})<br>" +
      "<i>%{customdata[0]}</i> (%{customdata[1]})<extra></extra>";
    xTitle = "Coverage of described species (%)";
  }

  const trace = {
    type: "bar",
    orientation: "h",
    x: values,
    y: labels,
    marker: { color: colors, line: { width: 0 } },
    text: valueText,
    textposition: "outside",
    textfont: { size: 11, color: "#1a1a1a", family: FONT_STACK },
    cliponaxis: false,
    customdata,
    hovertemplate,
    ...(errorBars ? { error_x: errorBars } : {}),
  };

  const chartHeight = Math.max(480, labels.length * 26 + 80);
  const container = document.getElementById("coverage-plot");
  container.style.height = chartHeight + "px";

  // x-axis range: 0–112 normally; expand for estimated mode if upper error bar pokes past it
  const maxX = mode === "estimated"
    ? Math.max(112, Math.ceil(Math.max(...rows.map(r => r.coverage_pct_estimated_high)) / 10) * 10 + 5)
    : 112;

  const layout = {
    margin: { l: 170, r: 70, t: 24, b: 52 },
    font: { family: FONT_STACK, size: 12, color: "#1a1a1a" },
    xaxis: {
      title: { text: xTitle, font: { size: 12, color: "#6b6b6b" }, standoff: 12 },
      range: [0, maxX],
      gridcolor: "#eee",
      zerolinecolor: "#ddd",
      ticksuffix: "%",
      tickfont: { size: 11, color: "#6b6b6b" },
      ticks: "outside",
      tickcolor: "#ddd",
    },
    yaxis: {
      automargin: false,
      tickfont: { size: 11.5, color: "#1a1a1a" },
      autorange: "reversed",
      ticks: "",
      ticksuffix: "   ",
    },
    bargap: 0.38,
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    showlegend: false,
    hoverlabel: {
      bgcolor: "#fff",
      bordercolor: "#e5e5e5",
      font: { size: 12, family: FONT_STACK, color: "#1a1a1a" },
      align: "left",
    },
  };

  Plotly.purge(container);
  Plotly.newPlot(container, [trace], layout, { displayModeBar: false, responsive: true });

  container.on("plotly_click", (ev) => {
    const point = ev.points && ev.points[0];
    if (!point) return;
    const row = rows[point.pointIndex];
    if (!row || !row.tree_name) return;
    const target = STATE.data.trees.find(t => t.filename === row.tree_name)
                || STATE.data.trees.find(t => t.group === row.group);
    if (target) selectTree(target.filename, true);
  });

  // Update the sub-line copy under the title to match the mode
  const sub = document.getElementById("coverage-sub");
  if (sub) {
    if (mode === "estimated") {
      sub.innerHTML = "Each bar is <strong>one published dataset</strong>. Coverage = tips in tree &divide; <strong>estimated true diversity</strong>. Error bars span the low–high estimate range from the cited source.";
    } else {
      sub.innerHTML = "Each bar is <strong>one published dataset</strong>. Coverage = tips in tree &divide; <strong>described</strong> species &times; 100. Click a bar to inspect that tree.";
    }
  }
}

function wireCoverageToggle() {
  document.querySelectorAll(".toggle-opt").forEach(btn => {
    btn.addEventListener("click", () => {
      const mode = btn.dataset.mode;
      if (mode === STATE.coverageMode) return;
      STATE.coverageMode = mode;
      document.querySelectorAll(".toggle-opt").forEach(b => {
        const on = b.dataset.mode === mode;
        b.classList.toggle("active", on);
        b.setAttribute("aria-pressed", on ? "true" : "false");
      });
      renderCoveragePlot();
    });
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

  const est = t.estimate;
  const describedVal = t.described_species
    ? `${fmt.format(t.described_species)}${t.described_source ? ` <span class="src">(${t.described_source})</span>` : ""}`
    : "—";
  const estimatedVal = est && est.estimated_total
    ? `${fmt.format(est.estimated_total)} <span class="range">(${fmt.format(est.estimated_low || est.estimated_total)}–${fmt.format(est.estimated_high || est.estimated_total)})</span>` +
      (est.estimate_source ? ` <span class="src">${est.estimate_source}${est.estimate_confidence ? `, ${est.estimate_confidence.toLowerCase()} conf.` : ""}</span>` : "")
    : null;
  const covDescribed = t.coverage_pct != null ? `${t.coverage_pct.toFixed(1)}% <span class="src">of described</span>` : "—";
  let covEstimated = null;
  if (est && est.estimated_total && t.ntips) {
    const main = (100 * t.ntips / est.estimated_total).toFixed(1);
    const lo = est.estimated_high ? (100 * t.ntips / est.estimated_high).toFixed(1) : null;
    const hi = est.estimated_low ? (100 * t.ntips / est.estimated_low).toFixed(1) : null;
    const subnote = t.is_partition_anchor
      ? `<span class="src">of ${t.partition_group || "partition"} estimate</span>`
      : `<span class="src">of ${t.partition_group} estimate &mdash; this tree is a sub-clade of ${t.partition_group}, so the % refers to the whole partition</span>`;
    covEstimated = `${main}% <span class="range">(${lo}&ndash;${hi}%)</span> ${subnote}`;
  }

  const fields = [
    ["Group", t.group + (t.partition_group && t.partition_group !== t.group ? ` <span class="src">→ ${t.partition_group}</span>` : "")],
    ["Study", t.study],
    ["Year", t.year],
    ["Journal", t.journal],
    ["Tips", fmt.format(t.ntips || 0)],
    ["Dated", t.dated ? "yes" : "no"],
    ["Crown age", t.crown_ma ? `${t.crown_ma} Ma` : "—"],
    ["Described species", describedVal],
    ["Estimated species", estimatedVal],
    ["Coverage (described)", covDescribed],
    ["Coverage (estimated)", covEstimated],
    ["DOI", t.doi ? `<a href="https://doi.org/${t.doi}" target="_blank" rel="noopener">${t.doi}</a>` : "—"],
    ["Source", t.data_source],
    ["File", t.size_bytes ? `${(t.size_bytes / 1024).toFixed(1)} KB` : "—"],
  ];
  document.getElementById("detail-meta").innerHTML = fields
    .filter(([_, v]) => v !== null && v !== undefined && v !== "")
    .map(([k, v]) => `<div class="field"><span class="key">${k}</span><span class="val">${v}</span></div>`)
    .join("");

  renderDatasetChooser(t);

  const dl = document.getElementById("download-newick");
  dl.href = TREES_BASE + filename;
  dl.setAttribute("download", filename);

  renderRSnippet(t);

  await loadAndRenderTree(t);
}

function renderDatasetChooser(t) {
  const host = document.getElementById("dataset-chooser");
  if (!host) return;
  const partition = t.partition_group;
  const peers = partition && STATE.data.datasets_by_partition
    ? (STATE.data.datasets_by_partition[partition] || [])
    : [];
  // For each peer dataset, the corresponding tree file is in STATE.data.trees.
  // A peer is interesting only if its tree exists and isn't the current tree.
  const items = peers
    .map(p => ({ peer: p, tree: STATE.data.trees.find(x => x.filename === p.tree_name) }))
    .filter(x => x.tree && x.tree.filename !== t.filename);

  if (!partition || items.length === 0) {
    host.hidden = true;
    host.innerHTML = "";
    return;
  }
  host.hidden = false;
  host.innerHTML = `
    <div class="chooser-label">Other datasets for <strong>${partition}</strong>:</div>
    <div class="chooser-buttons">
      ${items.map(({ peer, tree }) => {
        const label = peer.group;
        const meta = `${peer.tips ? fmt.format(peer.tips) + " tips" : ""}${peer.year ? " · " + peer.year : ""}`;
        return `<button type="button" class="chooser-btn" data-filename="${tree.filename}">
          <span class="chooser-name">${label}</span>
          <span class="chooser-meta">${meta}</span>
        </button>`;
      }).join("")}
    </div>`;
  host.querySelectorAll(".chooser-btn").forEach(btn => {
    btn.addEventListener("click", () => selectTree(btn.dataset.filename, false));
  });
}

function renderRSnippet(t) {
  const name = t.filename.replace(/\.nwk$/, "");
  const big = (t.ntips || 0) > 20000;
  const loadLine = big
    ? `tree <- load_atlas_tree("${name}", resolve_labels = FALSE)  # ${fmt.format(t.ntips)} tips`
    : `tree <- load_atlas_tree("${name}")`;
  const snippet =
    `# one-time install\n` +
    `# pak::pkg_install("franciscorichter/phylo-species-atlas/phyloatlas")\n\n` +
    `library(phyloatlas)\n` +
    `${loadLine}\n` +
    `plot(tree, show.tip.label = FALSE)`;
  const code = document.getElementById("r-snippet");
  code.textContent = snippet;

  const btn = document.getElementById("copy-r-snippet");
  btn.textContent = "Copy";
  btn.onclick = async () => {
    try {
      await navigator.clipboard.writeText(snippet);
      btn.textContent = "Copied";
      setTimeout(() => { btn.textContent = "Copy"; }, 1500);
    } catch {
      btn.textContent = "Copy failed";
    }
  };
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
