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
  query: "",
  datedOnly: false,
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

function sourceLinkHTML(info, fallbackLabel, confidence) {
  const label = (info && info.label) || fallbackLabel;
  if (!label) return "";
  const confSuffix = confidence ? `, ${confidence.toLowerCase()} conf.` : "";
  if (!info) return ` <span class="src">${label}${confSuffix}</span>`;
  const tipBits = [];
  if (info.citation) tipBits.push(info.citation);
  if (info.type === "paper" && info.doi) tipBits.push(`DOI: ${info.doi}`);
  const titleAttr = tipBits.length ? ` title="${tipBits.join(' — ').replace(/"/g, '&quot;')}"` : "";
  let href = null;
  if (info.doi) href = `https://doi.org/${info.doi}`;
  else if (info.url) href = info.url;
  if (!href) return ` <span class="src"${titleAttr}>${label}${confSuffix}</span>`;
  const typeTag = info.type === "paper" ? "paper" : "DB";
  return ` <a class="src-link" href="${href}" target="_blank" rel="noopener"${titleAttr}><span class="src-type ${info.type}">${typeTag}</span>${label}${confSuffix}</a>`;
}

function renderCoveragePlot() {
  const allRows = STATE.data.coverage;
  if (!allRows.length) return;
  const mode = STATE.coverageMode;
  const estimatedAvailable = (r) => r.coverage_pct_estimated != null;

  // Coverage chart shows one bar per partition — the canonical (largest)
  // tree. Sub-clade trees (parrots within Birds, primates within Mammals,
  // etc.) appear only when the user expands the partition in the tree
  // list on the left.
  const canonical = allRows.filter(r => r.is_partition_canonical);
  const rows = mode === "estimated"
    ? canonical.filter(estimatedAvailable)
               .slice().sort((a, b) => b.coverage_pct_estimated - a.coverage_pct_estimated)
    : canonical.slice().sort((a, b) => b.coverage_pct - a.coverage_pct);

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
    if (!row || !row.filename) return;
    selectTree(row.filename, true);
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
  STATE.query = q;
  STATE.datedOnly = datedOnly;

  let list = STATE.data.trees.slice();
  if (q) {
    list = list.filter(t =>
      (t.filename || "").toLowerCase().includes(q) ||
      (t.group || "").toLowerCase().includes(q) ||
      (t.partition_group || "").toLowerCase().includes(q) ||
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

const STATE_EXPANDED = new Set(); // partitions currently expanded in the sidebar

function treeRowHTML(t, { indent = false } = {}) {
  const cls = (STATE.selected === t.filename ? "active " : "") + (indent ? "child " : "");
  const badge = t.dated ? "dated" : "undated";
  const name = t.filename.replace(/\.nwk$/, "");
  return `<li data-filename="${t.filename}" class="tree-item ${cls.trim()}">
    <span class="tree-name"><span class="badge ${badge}"></span>${name}</span>
    <span class="tree-meta">${fmt.format(t.ntips || 0)}</span>
  </li>`;
}

function renderTreeList() {
  const ul = document.getElementById("tree-list");
  const list = STATE.filtered;
  const searching = STATE.query.length > 0 || STATE.datedOnly;

  let html;
  if (searching) {
    // Flat view when filtering — every match visible.
    html = list.map(t => treeRowHTML(t)).join("");
  } else {
    // Hierarchical view: partition headers → canonical → sub-clades (collapsed).
    // Bucket trees by partition; trees without a partition go to an "Other" bucket.
    const buckets = new Map();
    for (const t of list) {
      const p = t.partition_group || "Other";
      if (!buckets.has(p)) buckets.set(p, []);
      buckets.get(p).push(t);
    }

    // Order partitions: by canonical tips desc; then "Condamine 2019 families"
    // and "Other" pinned to the bottom.
    const partitionOrder = Array.from(buckets.keys())
      .map(p => {
        const canonical = buckets.get(p).find(t => t.is_partition_canonical)
                       || buckets.get(p)[0];
        return { partition: p, canonicalTips: canonical?.ntips || 0 };
      })
      .sort((a, b) => {
        const pinA = a.partition === "Condamine 2019 families" || a.partition === "Other" ? 1 : 0;
        const pinB = b.partition === "Condamine 2019 families" || b.partition === "Other" ? 1 : 0;
        if (pinA !== pinB) return pinA - pinB;
        return b.canonicalTips - a.canonicalTips;
      })
      .map(x => x.partition);

    html = partitionOrder.map(partition => {
      const items = buckets.get(partition);
      const canonical = items.find(t => t.is_partition_canonical);
      const subclades = items.filter(t => !t.is_partition_canonical);
      const hasSub = subclades.length > 0;
      const expanded = STATE_EXPANDED.has(partition);
      const chevron = hasSub ? `<span class="chev">${expanded ? "▾" : "▸"}</span>` : `<span class="chev empty"></span>`;
      const countBadge = hasSub ? `<span class="sub-count">+${subclades.length}</span>` : "";

      const headerRow = `<li class="partition-header${hasSub ? ' has-children' : ''}${expanded ? ' open' : ''}" data-partition="${partition}">
        ${chevron}<span class="partition-name">${partition}</span>${countBadge}
      </li>`;

      // For partitions without a canonical (e.g. Condamine 218 families), show
      // children as a generic list under the header. For others, show canonical
      // immediately; children appear when expanded.
      let canonicalRow = "";
      let childRows = "";
      if (canonical) {
        canonicalRow = treeRowHTML(canonical);
      }
      if (hasSub && expanded) {
        childRows = subclades.map(t => treeRowHTML(t, { indent: true })).join("");
      }
      return headerRow + canonicalRow + childRows;
    }).join("");

    // Append the dark-matter section: partitions with an estimate but no
    // shipped tree. Collapsible like the other partition headers.
    const unrepresented = STATE.data.unrepresented || [];
    if (unrepresented.length) {
      const SECTION_KEY = "__unrepresented__";
      const open = STATE_EXPANDED.has(SECTION_KEY);
      const chev = `<span class="chev">${open ? "▾" : "▸"}</span>`;
      let body = "";
      if (open) {
        // Sub-group by category.
        const byCat = new Map();
        for (const u of unrepresented) {
          const c = u.category || "Other";
          if (!byCat.has(c)) byCat.set(c, []);
          byCat.get(c).push(u);
        }
        for (const [cat, items] of byCat) {
          body += `<li class="dark-cat-header">${cat} <span class="sub-count">${items.length}</span></li>`;
          for (const u of items) {
            const desc = u.described ? fmt.format(u.described) + " described" : "—";
            const link = sourceLinkHTML(u.estimate_source_info, u.estimate_source);
            body += `<li class="tree-item dark child">
              <span class="tree-name"><span class="badge dark"></span>${u.group}${link}</span>
              <span class="tree-meta">${desc}</span>
            </li>`;
          }
        }
      }
      html += `<li class="partition-header dark-header has-children${open ? ' open' : ''}" data-partition="${SECTION_KEY}">
        ${chev}<span class="partition-name">Not yet represented</span><span class="sub-count">${unrepresented.length}</span>
      </li>${body}`;
    }
  }

  ul.innerHTML = html;

  ul.querySelectorAll(".partition-header.has-children").forEach(li => {
    li.addEventListener("click", () => {
      const p = li.dataset.partition;
      if (STATE_EXPANDED.has(p)) STATE_EXPANDED.delete(p);
      else STATE_EXPANDED.add(p);
      renderTreeList();
    });
  });
  ul.querySelectorAll(".tree-item").forEach(li => {
    if (!li.dataset.filename) return; // dark-matter rows have no tree to load
    li.addEventListener("click", () => selectTree(li.dataset.filename, false));
  });
}

async function selectTree(filename, scrollIntoView) {
  STATE.selected = filename;
  // If the selected tree is a sub-clade, make sure its partition is expanded
  // so the user can see it in the sidebar.
  const t0 = STATE.data.trees.find(x => x.filename === filename);
  if (t0 && t0.partition_group && !t0.is_partition_canonical) {
    STATE_EXPANDED.add(t0.partition_group);
  }
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
    ? `${fmt.format(t.described_species)}${sourceLinkHTML(t.described_source_info, t.described_source)}`
    : "—";
  const estimatedVal = est && est.estimated_total
    ? `${fmt.format(est.estimated_total)} <span class="range">(${fmt.format(est.estimated_low || est.estimated_total)}–${fmt.format(est.estimated_high || est.estimated_total)})</span>` +
      sourceLinkHTML(t.estimate_source_info, est.estimate_source, est.estimate_confidence)
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

  const treePaperHTML = t.doi
    ? `<a class="src-link" href="https://doi.org/${t.doi}" target="_blank" rel="noopener" title="${(t.study || "") + (t.year ? " (" + t.year + ")" : "") + (t.journal ? " — " + t.journal : "")}"><span class="src-type paper">paper</span>${(t.study || "").split(",")[0].split(" et")[0] || t.doi}${t.year ? " " + t.year : ""}</a>`
    : (t.study ? `<span class="src">${t.study}${t.year ? " (" + t.year + ")" : ""}</span>` : "—");

  const fields = [
    ["Group", t.group + (t.partition_group && t.partition_group !== t.group ? ` <span class="src">→ ${t.partition_group}</span>` : "")],
    ["Tree paper", treePaperHTML],
    ["Journal", t.journal],
    ["Tips", fmt.format(t.ntips || 0)],
    ["Dated", t.dated ? "yes" : "no"],
    ["Crown age", t.crown_ma ? `${t.crown_ma} Ma` : "—"],
    ["Described species", describedVal],
    ["Estimated species", estimatedVal],
    ["Coverage (described)", covDescribed],
    ["Coverage (estimated)", covEstimated],
    ["Tree file source", t.data_source],
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

  const items = peers
    .map(p => ({ peer: p, tree: STATE.data.trees.find(x => x.filename === p.filename) }))
    .filter(x => x.tree);

  // Hide when there's no partition or it has only one tree (nothing to switch).
  if (!partition || items.length <= 1) {
    host.hidden = true;
    host.innerHTML = "";
    return;
  }
  host.hidden = false;

  const canonical = items.find(x => x.peer.is_partition_canonical) || items[0];
  const subclades = items.filter(x => !x.peer.is_partition_canonical);

  const renderBtn = ({ peer, tree }, role) => {
    const isActive = tree.filename === t.filename;
    const tips = peer.tips ? `${fmt.format(peer.tips)} tips` : "";
    const year = peer.year ? `${peer.year}` : "";
    const datedBadge = peer.dated
      ? `<span class="chooser-badge dated" title="dated"></span>`
      : `<span class="chooser-badge undated" title="undated"></span>`;
    const studyShort = (peer.study || "").split(",")[0].split(" et")[0];
    const meta = [tips, year, studyShort].filter(Boolean).join(" · ");
    return `<button type="button" class="chooser-btn ${role} ${isActive ? "active" : ""}" data-filename="${tree.filename}" aria-pressed="${isActive}">
        <span class="chooser-name">${datedBadge}${peer.group}</span>
        <span class="chooser-meta">${meta}</span>
      </button>`;
  };

  host.innerHTML = `
    <div class="chooser-label">Partition <strong>${partition}</strong></div>
    <div class="chooser-tier">
      <div class="tier-label">Whole-partition tree</div>
      <div class="chooser-buttons">${renderBtn(canonical, "canonical")}</div>
    </div>
    ${subclades.length ? `
    <div class="chooser-tier">
      <div class="tier-label">Sub-clade trees <span class="tier-count">${subclades.length}</span></div>
      <div class="chooser-buttons">${subclades.map(x => renderBtn(x, "subclade")).join("")}</div>
    </div>` : ""}
  `;

  host.querySelectorAll(".chooser-btn").forEach(btn => {
    btn.addEventListener("click", () => {
      if (btn.dataset.filename === STATE.selected) return;
      selectTree(btn.dataset.filename, false);
    });
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
