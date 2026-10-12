// Phylo-Species Atlas browser — vanilla JS, no build step.
// Loads data.json, renders summary + plots, and lazy-loads tree files.

const STATE = {
  data: null,
  selected: null,
  layout: "linear", // or "radial"
  showLabels: false,
  currentNewick: null,
  currentNames: {}, // {id: standardized_name} for the active tree
};

const TREES_BASE = "../standardized/trees/"; // patched by CI for Pages deploy
const NAMES_BASE = "../standardized/names/"; // patched by CI for Pages deploy
const TIP_RENDER_LIMIT = 12000; // above this, we sample for rendering. Covers fish/birds/squamates/amphibians/ferns/mammals/frogs at full resolution; only seed_plants (~342k) still gets sampled.

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
  wireTreeViewerControls();
});

function renderSummary() {
  const s = STATE.data.summary;
  const levels = STATE.data.levels || [];

  // Per-category group counts (with a reference tree / all groups) for the
  // mini-bar legend: the same rows as the chart.
  const catStats = new Map();
  for (const r of (STATE.data.coverage || [])) {
    const c = r.category || "Other";
    if (!catStats.has(c)) catStats.set(c, { total: 0, withTree: 0 });
    const entry = catStats.get(c);
    entry.total += 1;
    if (!r.is_unrepresented) entry.withTree += 1;
  }

  // The four coverage levels of the paper, each over the described
  // eukaryotic species; values are copied from data.json as printed.
  const levelStats = levels.map(l => `
        <div class="stat">
          <span class="value">${l.pct_text}%</span>
          <span class="label">${escapeHTML(l.label)}</span>
          <span class="sub">${l.n_text} species · ${l.n_groups} groups</span>
        </div>`).join("");
  const coverageRow = `
    <div class="stats-group">
      <div class="stats-group-label">Coverage of the ${s.described_eukaryotic_text} described eukaryotic species (${s.n_partitions_eukaryotic} groups; ${s.n_partitions_with_tree_eukaryotic} with a reference tree)</div>
      <div class="stats-row">${levelStats}
      </div>
      <div class="stats-row">
        <div class="stat">
          <span class="value">${s.n_partitions_with_tree} / ${s.n_partitions}</span>
          <span class="label">Groups with a reference tree</span>
        </div>
        <div class="stat">
          <span class="value">${fmt.format(s.n_trees)}</span>
          <span class="label">Trees redistributed (reference + other trees)</span>
        </div>
        <div class="stat stat-stack">
          <span class="value">${s.described_prokaryotic_text}</span>
          <span class="label">Bacteria and Archaea: GTDB species clusters</span>
          <span class="sub">counted on their own, never in the eukaryotic levels</span>
        </div>
      </div>
    </div>`;

  const catOrder = ["Vertebrates", "Plants", "Arthropods", "Other animals", "Microbes & protists"];
  const maxTotal = Math.max(...catOrder.map(c => catStats.get(c)?.total || 0));
  const catRow = `
    <div class="stats-group stats-group-cat">
      <div class="stats-group-label">Groups with a reference tree, by category</div>
      <div class="cat-bars">
        ${catOrder.map(c => {
          const e = catStats.get(c) || { total: 0, withTree: 0 };
          const widthPct = maxTotal ? Math.round(100 * e.total / maxTotal) : 0;
          const fillPct = e.total ? Math.round(100 * e.withTree / e.total) : 0;
          return `
            <div class="cat-bar-row">
              <span class="cat-bar-label">${c}</span>
              <div class="cat-bar" style="width:${widthPct}%">
                <div class="cat-bar-fill" style="width:${fillPct}%"></div>
              </div>
              <span class="cat-bar-count">${e.withTree}/${e.total}</span>
            </div>`;
        }).join("")}
      </div>
    </div>`;

  document.getElementById("summary-stats").innerHTML = coverageRow + catRow;

  const p = STATE.data.paper || {};
  const foot = document.getElementById("paper-ref");
  if (foot && p.title) {
    foot.innerHTML = `Companion to Richter (2026), <a href="https://doi.org/${escapeHTML(p.doi || "")}" target="_blank" rel="noopener"><em>${escapeHTML(p.title)}</em></a>. Every number on this page is read from the paper's files in this release (<code>registers/coverage_basis.csv</code>, <code>paper/</code>).`;
  }
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

// Macro-category ordering matches the paper's Figure 1 narrative.
const CATEGORY_ORDER = ["Vertebrates", "Plants", "Arthropods", "Other animals", "Microbes & protists"];  // 5 high-level groups
const CATEGORY_RANK = Object.fromEntries(CATEGORY_ORDER.map((c, i) => [c, i]));

function levelsText(r) {
  // Four levels of one group, as the paper's Table 1 prints them.
  const L = r.levels || {};
  const line = (k, label) => L[k] ? `${label}: ${L[k].pct_text}% (${L[k].n_text})` : "";
  return [line("named", "Named"), line("molecular", "Molecular"),
          line("dated", "Dated"), line("uncertainty", "Date uncertainty")]
    .filter(Boolean).join("<br>");
}

function renderCoveragePlot() {
  const allRows = STATE.data.coverage;
  if (!allRows.length) return;

  // One row per group (62). Markers, on a log scale of species:
  //   ◯ described species (the paper's denominator)
  //   ◇ named in the reference tree (accepted species; blue = chronogram in
  //     absolute time, orange = not)
  //   ● placed with molecular data
  //   ━ estimated total of the group (low–high), from the paper's Figure 1 data
  const rows = allRows.slice().sort((a, b) => {
    const ra = CATEGORY_RANK[a.category] ?? 99;
    const rb = CATEGORY_RANK[b.category] ?? 99;
    if (ra !== rb) return ra - rb;
    return (b.estimated_total ?? b.described_species ?? 0) - (a.estimated_total ?? a.described_species ?? 0);
  });

  const labels = rows.map(r => r.group);
  const tickText = rows.map(r => r.group);

  const intervalClass = (r) => {
    const ip = r.interval_provenance;
    return (ip && ip.overall_classification) || "unaudited";
  };
  const classLabel = (cls) => ({
    "paper-published-range": "paper-published range",
    "partly-heuristic": "partly heuristic (one bound is the described count or atlas-derived)",
    "fully-heuristic": "heuristic — no true-diversity estimate published",
    "derived": "derived from described count",
    "unaudited": "not yet audited for interval provenance",
  }[cls] || cls);
  const intervalColor = (r) => {
    const cls = intervalClass(r);
    return (cls === "fully-heuristic" || cls === "partly-heuristic")
      ? "rgba(180,140,90,0.55)"
      : "#5a5a5a";
  };

  const intervalShapes = rows.map(r => {
    if (r.estimated_low == null || r.estimated_high == null) return null;
    if (r.estimated_low >= r.estimated_high) return null;
    return {
      type: "line",
      xref: "x", yref: "y",
      x0: r.estimated_low, x1: r.estimated_high,
      y0: r.group, y1: r.group,
      line: { color: intervalColor(r), width: 3 },
      layer: "above",
    };
  }).filter(Boolean);
  const intervalAnchor = (r) => {
    if (r.estimated_low && r.estimated_high) {
      return Math.sqrt(r.estimated_low * r.estimated_high);
    }
    return r.estimated_total ?? null;
  };
  const estimatedHoverTrace = {
    type: "scatter",
    mode: "markers",
    name: "Estimated total (low–high)",
    x: rows.map(intervalAnchor),
    y: labels,
    marker: {
      symbol: "line-ns-open",
      size: 12,
      color: "rgba(0,0,0,0)",
      line: { color: "rgba(0,0,0,0)", width: 0 },
    },
    hoverinfo: "text",
    customdata: rows.map(r => [
      r.estimated_low, r.estimated_high, r.estimate_source || "—",
      r.estimate_confidence || "—", classLabel(intervalClass(r)),
      r.estimated_total,
    ]),
    hovertemplate:
      "<b>%{y}</b><br>" +
      "Estimated total: %{customdata[5]:,} (range %{customdata[0]:,}–%{customdata[1]:,})<br>" +
      "<i style='color:#7a7a7a'>%{customdata[4]}</i><br>" +
      "Source: %{customdata[2]} · confidence %{customdata[3]}<extra></extra>",
  };

  const describedMarkerColor = rows.map(r => r.is_unrepresented ? "#9a9a9a" : "#1a1a1a");
  const describedTrace = {
    type: "scatter",
    mode: "markers",
    name: "Described species",
    x: rows.map(r => r.described_species ?? null),
    y: labels,
    marker: {
      symbol: "circle-open",
      size: 11,
      color: describedMarkerColor,
      line: { color: describedMarkerColor, width: 1.5 },
    },
    customdata: rows.map(r => [
      r.described_text || "—", r.described_source || "—",
      r.is_unrepresented ? " · no reference tree" : "",
      r.is_prokaryotic ? " · GTDB species clusters, not in the eukaryotic levels" : "",
    ]),
    hovertemplate:
      "<b>%{y}</b>%{customdata[2]}%{customdata[3]}<br>" +
      "Described species: %{customdata[0]}<br>" +
      "Source: %{customdata[1]}<extra></extra>",
  };

  const withTree = rows.filter(r => !r.is_unrepresented);
  const namedTrace = {
    type: "scatter",
    mode: "markers",
    name: "Named in the reference tree (blue: chronogram in absolute time)",
    x: rows.map(r => r.levels ? r.levels.named.n : null),
    y: labels,
    marker: {
      symbol: "diamond",
      size: 10,
      color: rows.map(r => r.dated ? "#2a6fbf" : "#c97a2a"),
      line: { color: "#1a1a1a", width: 0.5 },
    },
    customdata: rows.map(r => [
      r.levels ? levelsText(r) : "",
      r.filename || "", r.time_scale_text || "",
      r.interval_kind_text ? ` · ${r.interval_kind_text}` : "",
    ]),
    hovertemplate:
      "<b>%{y}</b> · %{customdata[1]}<br>" +
      "%{customdata[0]}<br>" +
      "<i style='color:#7a7a7a'>%{customdata[2]}%{customdata[3]}</i><extra></extra>",
  };
  const molecularTrace = {
    type: "scatter",
    mode: "markers",
    name: "Placed with molecular data",
    x: rows.map(r => r.levels ? r.levels.molecular.n : null),
    y: labels,
    marker: { symbol: "circle", size: 6, color: "#1a1a1a" },
    customdata: rows.map(r => [r.levels ? r.levels.molecular.n_text : "", r.levels ? r.levels.molecular.pct_text : ""]),
    hovertemplate: "<b>%{y}</b><br>Placed with molecular data: %{customdata[0]} (%{customdata[1]}%)<extra></extra>",
  };

  const bands = [];
  for (let i = 0; i < rows.length; i++) {
    const cat = rows[i].category || "Other";
    const last = bands[bands.length - 1];
    if (!last || last.category !== cat) bands.push({ category: cat, start: i, end: i });
    else last.end = i;
  }
  const LEFT_MARGIN = 230;
  const annotations = bands.map(b => ({
    xref: "paper", yref: "y",
    x: 0, xanchor: "left", xshift: -LEFT_MARGIN + 6,
    y: rows[b.start].group, yanchor: "middle",
    text: `<b>${b.category.toUpperCase()}</b>`,
    showarrow: false, align: "left",
    font: { family: FONT_STACK, size: 10.5, color: "#7a7a7a" },
  }));
  const shapes = [
    ...bands.slice(1).map(b => ({
      type: "line", xref: "paper", yref: "y",
      x0: 0, x1: 1, y0: rows[b.start].group, y1: rows[b.start].group,
      line: { color: "#e4e4e4", width: 1, dash: "dot" }, layer: "below",
    })),
    ...intervalShapes,
  ];

  const container = document.getElementById("coverage-plot");
  container.style.height = Math.max(560, labels.length * 28 + bands.length * 8 + 110) + "px";

  const allXVals = rows.flatMap(r => [r.levels && r.levels.molecular.n, r.described_species, r.estimated_high])
    .filter(v => v != null && v > 0);
  const maxX = Math.pow(10, Math.ceil(Math.log10(Math.max(...allXVals))));
  const minX = 10;

  const layout = {
    margin: { l: LEFT_MARGIN, r: 30, t: 18, b: 64 },
    annotations,
    shapes,
    font: { family: FONT_STACK, size: 12, color: "#1a1a1a" },
    xaxis: {
      title: { text: "Species (log scale) — ◇ named in the reference tree · ● placed with molecular data · ◯ described · ━ estimated total (low–high)", font: { size: 12, color: "#6b6b6b" }, standoff: 14 },
      type: "log",
      range: [Math.log10(minX), Math.log10(maxX)],
      gridcolor: "#eee",
      zerolinecolor: "#ddd",
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
      tickmode: "array",
      tickvals: labels,
      ticktext: tickText,
    },
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    showlegend: true,
    legend: {
      orientation: "h", y: 1.04, x: 0, xanchor: "left",
      font: { size: 11, family: FONT_STACK },
      bgcolor: "rgba(255,255,255,0)",
    },
    hoverlabel: {
      bgcolor: "#fff", bordercolor: "#e5e5e5",
      font: { size: 12, family: FONT_STACK, color: "#1a1a1a" },
      align: "left",
    },
  };

  Plotly.purge(container);
  Plotly.newPlot(container, [estimatedHoverTrace, describedTrace, namedTrace, molecularTrace], layout,
                 { displayModeBar: false, responsive: true });

  container.on("plotly_click", (ev) => {
    const point = ev.points && ev.points[0];
    if (!point) return;
    const row = rows[point.pointIndex];
    if (!row) return;
    if (row.is_unrepresented) {
      selectDarkMatter(row.partition_group, true);
      return;
    }
    selectTree(row.filename, true);
  });

  const sub = document.getElementById("coverage-sub");
  if (sub) {
    sub.innerHTML = `Each row = one group (${withTree.length} with a reference tree, ${rows.length - withTree.length} without). <strong>◇</strong> = species named in the reference tree (accepted names; chronogram in absolute time <span style='color:#2a6fbf'>■</span> / not <span style='color:#c97a2a'>■</span>), <strong>●</strong> = placed with molecular data, <strong>◯</strong> = described species. Hover a row for its four coverage levels (Table 1 of the paper). <strong>━━━━</strong> = estimated total of the group (low–high; <em>solid grey</em> = range published by a study, <em>faded amber</em> = heuristic). Click a row with a tree to inspect it; click a row without one to see its notes.`;
  }
}

function wireTreeViewerControls() {
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

// Select a dark-matter (no-tree) partition: show the partition's audit
// info — interval, evidence-paper hyperlinks, constituents — instead of
// trying to load a tree that doesn't exist.
function selectDarkMatter(partitionName, scrollIntoView) {
  STATE.selected = null;
  document.getElementById("detail-empty").hidden = true;
  document.getElementById("detail").hidden = false;
  document.getElementById("detail-name").textContent = partitionName;

  // Hide the tree-related cards.
  document.getElementById("dataset-chooser").hidden = true;
  document.getElementById("detail-meta").innerHTML = "";
  document.getElementById("uncertainty-panel").hidden = true;
  document.getElementById("r-snippet-card").hidden = true;
  document.getElementById("tree-viewer-card").hidden = true;
  document.getElementById("dark-matter-panel").hidden = false;

  const audit = (STATE.data.audits || {})[partitionName] || {};
  const est = audit.estimate || {};
  const src = est.source || {};
  const ip = est.interval_provenance || {};
  const canTree = audit.canonical_tree || {};
  const constituents = audit.constituents || [];

  const classLabel = (cls) => ({
    "paper-published-range": "paper-published range",
    "partly-heuristic": "partly heuristic",
    "fully-heuristic": "heuristic — no true-diversity estimate published",
    "derived": "derived from described count",
    "unaudited": "not yet audited for interval provenance",
  }[cls] || cls);
  const clsClass = ip.overall_classification || "unaudited";

  // Paper hyperlinks for the interval evidence.
  let srcLinks = "";
  if (src.doi) {
    srcLinks += `<a href="https://doi.org/${escapeHTML(src.doi)}" target="_blank" rel="noopener" class="ev-link"><span class="src-type paper">paper</span>${escapeHTML(src.key || "source")} · DOI ${escapeHTML(src.doi)}</a>`;
  } else if (src.url) {
    srcLinks += `<a href="${escapeHTML(src.url)}" target="_blank" rel="noopener" class="ev-link"><span class="src-type database">database</span>${escapeHTML(src.key || "source")}</a>`;
  } else if (src.key) {
    srcLinks += `<span class="src">${escapeHTML(src.key)}</span>`;
  }

  // Range row.
  // Numbers come from the group row (the paper's files); the note supplies text.
  const row = (STATE.data.coverage || []).find(r => r.partition_group === partitionName) || {};
  const rangeText = (row.estimated_low != null && row.estimated_high != null)
    ? `${fmt.format(row.estimated_low)}–${fmt.format(row.estimated_high)} <span class="src">estimated total ${fmt.format(row.estimated_total)}${row.estimate_source ? ", " + escapeHTML(row.estimate_source) : ""}</span>`
    : (row.estimated_total != null ? `${fmt.format(row.estimated_total)} (single point)` : "—");
  const siteTrees = (STATE.data.trees || []).filter(t => t.partition_group === partitionName);
  const siteTreesHTML = siteTrees.length ? `
      <section class="dm-section">
        <h3>Trees of a part of this group (${siteTrees.length})</h3>
        <p class="dm-bnote">Not a reference tree: not counted in the group's coverage.</p>
        <div class="chooser-buttons">${siteTrees.map(t => `<button type="button" class="chooser-btn subclade" data-filename="${escapeHTML(t.filename)}">
          <span class="chooser-name">${escapeHTML(t.group)}</span><span class="chooser-meta">${t.ntips_text} tips · ${escapeHTML((t.study || "").split(",")[0].split(" -")[0].split(" —")[0])}</span></button>`).join("")}</div>
      </section>` : "";

  // Constituents table (only when present, e.g. for "Others").
  let constituentsHTML = "";
  if (constituents.length) {
    constituentsHTML = `
      <section class="dm-section">
        <h3>Constituent clades (${constituents.length})</h3>
        <table class="dm-table">
          <thead><tr>
            <th>Clade</th><th>Described</th><th>Estimated</th><th>Notes</th>
          </tr></thead>
          <tbody>
            ${constituents.map(c => `<tr>
              <td><strong>${escapeHTML(c.name)}</strong></td>
              <td>${fmt.format(c.described || 0)}</td>
              <td>${fmt.format(c.estimated || 0)}</td>
              <td class="dm-note">${c.note ? escapeHTML(c.note) : ""}</td>
            </tr>`).join("")}
          </tbody>
        </table>
      </section>`;
  }

  // Interval provenance per-bound, with verbatim quotes when present.
  const renderBound = (label, bound) => {
    if (!bound) return "";
    const note = bound.note ? `<div class="dm-bnote">${escapeHTML(bound.note)}</div>` : "";
    const quote = bound.verbatim_quote ? `<div class="dm-quote">“${escapeHTML(bound.verbatim_quote)}”</div>` : "";
    return `<tr>
      <th>${escapeHTML(label)}</th>
      <td><span class="ip-type ip-type-${(bound.type || "").replace(/-/g, "_")}">${escapeHTML(bound.type || "—")}</span></td>
      <td>${escapeHTML(bound.source || "—")}${quote}${note}</td>
    </tr>`;
  };

  const intervalProvHTML = (ip && Object.keys(ip).length > 0) ? `
    <section class="dm-section">
      <h3>Interval provenance <span class="ip-cls ip-cls-${clsClass}">${escapeHTML(classLabel(clsClass))}</span></h3>
      <table class="ip-table">
        ${renderBound("Described", ip.described)}
        ${renderBound("Estimated total", ip.estimated_total)}
        ${renderBound("Low bound", ip.estimated_low)}
        ${renderBound("High bound", ip.estimated_high)}
      </table>
      ${ip.auditor_note ? `<p class="dm-bnote">${escapeHTML(ip.auditor_note)}</p>` : ""}
    </section>` : "";

  const host = document.getElementById("dark-matter-panel");
  host.innerHTML = `
    <div class="card dm-card">
      <div class="card-header">
        <h3>No reference tree for this group</h3>
        <span class="audit-badge audit-deferred">${escapeHTML((audit.audit || {}).status || "deferred")}</span>
      </div>
      <p class="card-sub">${escapeHTML(audit.category || "")} · The search found no species-level tree of this group; its described species count in the denominator of every coverage level and in no numerator.</p>

      ${canTree.reason ? `<section class="dm-section">
        <h3>Why no tree</h3>
        <p class="dm-reason">${escapeHTML(canTree.reason)}</p>
      </section>` : ""}

      <section class="dm-section">
        <h3>Diversity estimate</h3>
        <dl class="dm-dl">
          <dt>Described species</dt>
          <dd>${row.described_text || "—"} <span class="src">${escapeHTML(row.described_source || "")}</span></dd>
          <dt>Estimated range</dt>
          <dd>${rangeText}</dd>
          ${est.confidence ? `<dt>Confidence</dt><dd>${escapeHTML(est.confidence)}</dd>` : ""}
          <dt>Evidence</dt>
          <dd>${srcLinks || '<span class="src">—</span>'}</dd>
        </dl>
      </section>

      ${siteTreesHTML}
      ${constituentsHTML}
      ${intervalProvHTML}

      <p class="dm-footer"><a href="audit.html#${encodeURIComponent(partitionName)}" target="_blank" rel="noopener">Full audit entry →</a></p>
    </div>
  `;

  host.querySelectorAll(".chooser-btn").forEach(btn => {
    btn.addEventListener("click", () => selectTree(btn.dataset.filename, false));
  });

  if (scrollIntoView) {
    document.getElementById("detail").scrollIntoView({ block: "start", behavior: "smooth" });
  }
}

async function selectTree(filename, scrollIntoView) {
  STATE.selected = filename;
  const t = STATE.data.trees.find(x => x.filename === filename);
  if (!t) return;
  if (scrollIntoView) {
    document.getElementById("detail").scrollIntoView({ block: "start", behavior: "smooth" });
  }

  document.getElementById("detail-empty").hidden = true;
  document.getElementById("detail").hidden = false;
  document.getElementById("detail-name").textContent = t.filename.replace(/\.nwk$/, "");
  document.getElementById("dark-matter-panel").hidden = true;
  document.getElementById("r-snippet-card").hidden = !t.in_package;
  document.getElementById("tree-viewer-card").hidden = false;

  // The group row: the paper's values for the group whose reference tree this is.
  const row = t.is_partition_canonical
    ? (STATE.data.coverage || []).find(r => r.filename === t.filename)
    : null;

  const treePaperHTML = t.doi
    ? `<a class="src-link" href="https://doi.org/${t.doi}" target="_blank" rel="noopener" title="${escapeHTML((t.study || "") + (t.year ? " (" + t.year + ")" : "") + (t.journal ? " — " + t.journal : ""))}"><span class="src-type paper">paper</span>${escapeHTML((t.study || "").split(",")[0].split(" et")[0] || t.doi)}${t.year ? " " + t.year : ""}</a>`
    : (t.study ? `<span class="src">${escapeHTML(t.study)}${t.year ? " (" + t.year + ")" : ""}</span>` : "—");

  const audit = STATE.data.audits?.[t.partition_group];
  const auditMethods = t.is_partition_canonical ? audit?.tree?.methods : t.methods;

  let roleHTML;
  if (t.is_partition_canonical) roleHTML = `reference tree of ${escapeHTML(t.partition_group)}`;
  else if (t.site_tree === "derived") roleHTML = `derived in this site from <code>${escapeHTML(t.derived_from)}</code>, one tip per species; not redistributed in <code>standardized/trees/</code>`;
  else if (t.site_tree === "site_only") roleHTML = `held by this site only; not redistributed in <code>standardized/trees/</code> and not used by the paper`;
  else roleHTML = t.partition_group ? `other tree of ${escapeHTML(t.partition_group)}; not counted in its coverage` : "other tree; not counted in any group's coverage";

  const fields = [
    ["Group", escapeHTML(t.partition_group || "—")],
    ["Role", roleHTML],
    ["Tree paper", treePaperHTML],
    ["Tips", `${t.ntips_text} tips`],
  ];

  if (row) {
    const L = row.levels;
    const lv = (k) => `${L[k].pct_text}% <span class="src">${L[k].n_text} of ${row.described_text} described</span>`;
    fields.push(
      ["Described species", `${row.described_text} <span class="src">${escapeHTML(row.described_source || "")}</span>`],
      ["Named in the tree", lv("named") + (row.is_prokaryotic ? ` <span class="src">GTDB species clusters</span>` : ` <span class="src">accepted names (Catalogue of Life COL26.9)</span>`)],
      ["Placed with molecular data", lv("molecular")],
      ["In a chronogram in absolute time", row.dated ? lv("dated") : `no <span class="src">${escapeHTML(row.time_scale_text || "")}</span>`],
      ["Date uncertainty available", row.uncertainty_available
        ? `${lv("uncertainty")} <span class="src">${escapeHTML(row.interval_kind_text || "")}</span>`
        : `no${row.interval_kind_text ? ` <span class="src">${escapeHTML(row.interval_kind_text)}</span>` : ""}`],
    );
    if (row.crown) {
      const c = row.crown;
      const iv = c.interval_lo_ma_text
        ? ` <span class="range">(${c.interval_lo_ma_text}–${c.interval_hi_ma_text} Ma${c.n_trees > 1 ? ", over " + fmt.format(c.n_trees) + " trees" : ""})</span>`
        : "";
      fields.push(["Crown age", `${c.age_ma_text} Ma${iv} <span class="src">${escapeHTML(c.clade)}; ${escapeHTML(c.age_statistic)}</span>`]);
    } else {
      fields.push(["Crown age", "—"]);
    }
    if (row.estimated_total) {
      fields.push(["Estimated total", `${fmt.format(row.estimated_total)} <span class="range">(${fmt.format(row.estimated_low ?? row.estimated_total)}–${fmt.format(row.estimated_high ?? row.estimated_total)})</span> <span class="src">${escapeHTML(row.estimate_source || "")}${row.estimate_confidence ? ", " + escapeHTML(row.estimate_confidence.toLowerCase()) + " conf." : ""}</span>`]);
    }
  } else {
    fields.push(["Dated", t.dated_text ? escapeHTML(t.dated_text) : (t.dated ? "yes" : "no")]);
    if (t.paper) {
      const p = t.paper;
      fields.push(
        ["Accepted names", p.accepted_names_text || "—"],
        ["Described species", p.described_text ? `${p.described_text} <span class="src">${escapeHTML(p.described_source || "")}</span>` : "—"],
        ["Named", p.named_pct_text ? `${p.named_pct_text}% <span class="src">of the described species of this clade (Supplement table of the trees)</span>` : "—"],
      );
    }
    if (t.derivation) fields.push(["Derivation", escapeHTML(t.derivation)]);
  }

  if (auditMethods) {
    const m = auditMethods;
    const methodFields = [
      ["Inference", m.inference],
      ["Software", Array.isArray(m.software) ? m.software.join(", ") : m.software],
      ["Data type", m.data_type],
      ["Loci / genes", m.n_loci || m.n_genes],
      ["Total bp", m.total_bp ? fmt.format(m.total_bp) : null],
      ["Substitution model", m.substitution_model],
      ["Dating method", m.dating_method],
      ["Fossil calibrations", m.fossil_calibrations],
      ["Support metric", m.support_metric],
    ];
    for (const [k, v] of methodFields) {
      if (v != null && v !== "") fields.push([k, escapeHTML(String(v))]);
    }
  }

  fields.push(
    ["Tree file source", t.data_source ? escapeHTML(t.data_source) : null],
    ["File", `<code>${escapeHTML(t.tree_url)}</code> · ${(t.size_bytes / 1024).toFixed(1)} KB`],
  );

  document.getElementById("detail-meta").innerHTML = fields
    .filter(([_, v]) => v !== null && v !== undefined && v !== "")
    .map(([k, v]) => `<div class="field"><span class="key">${k}</span><span class="val">${v}</span></div>`)
    .join("");

  renderDatasetChooser(t);
  renderAuditPanel(t);
  renderUncertaintyPanel(t);

  const dl = document.getElementById("download-newick");
  dl.href = treeURL(t);
  dl.setAttribute("download", t.filename);

  if (t.in_package) renderRSnippet(t);

  await loadAndRenderTree(t);
}

// The file the viewer shows and the download link offers: the redistributed
// file in standardized/trees/, or, for a tree only this site holds, its file
// under data/partitions/.
function treeURL(t) {
  return t.in_package ? TREES_BASE + t.filename : t.tree_url_override;
}

// Render the audit findings panel for the partition the selected tree belongs
// to. Pulls from STATE.data.audits[partition] which is populated from
// site/data/partitions/<slug>/info.yaml at build time. Methods + uncertainty
// from the audit are folded inline into detail-meta (above), so this panel
// only carries the status badge and the DOI-correction block.
function renderAuditPanel(t) {
  const host = document.getElementById("audit-panel");
  if (!host) return;
  const audits = STATE.data.audits || {};
  const partition = t.partition_group;
  const audit = partition && audits[partition];
  if (!audit) {
    host.hidden = true;
    host.innerHTML = "";
    return;
  }
  const status = (audit.audit && audit.audit.status) || "stub";
  const lastAudited = audit.audit && audit.audit.last_audited;
  const estSrc = (audit.estimate && audit.estimate.source) || {};

  let estSrcHTML = "";
  if (estSrc.paper_doi_status === "broken" && estSrc.live_doi) {
    estSrcHTML = `
      <div class="audit-correction">
        <strong>Estimate source DOI was broken</strong> and has been replaced for the web display.
        <div class="audit-doi-pair">
          <span class="doi-broken"><span class="src">paper cites:</span>
            <a href="https://doi.org/${estSrc.paper_doi}" target="_blank" rel="noopener">${estSrc.paper_key} · ${estSrc.paper_doi}</a>
            <span class="badge-broken">404</span>
          </span>
          <span class="doi-live"><span class="src">live:</span>
            <a href="https://doi.org/${estSrc.live_doi}" target="_blank" rel="noopener">${estSrc.live_key || "live"} · ${estSrc.live_doi}</a>
            ${estSrc.live_cited_value ? `<span class="src"> · ${fmt.format(estSrc.live_cited_value)} species</span>` : ""}
          </span>
        </div>
        ${estSrc.note ? `<p class="audit-note">${escapeHTML(estSrc.note).replace(/\n/g, "<br>")}</p>` : ""}
      </div>`;
  }

  // Caveats — known data-integrity issues with the shipped tree
  // (e.g., Crustaceans: Decapoda-only; Bryozoa: chronogram not deposited).
  const caveats = audit.caveats || [];
  const caveatsHTML = caveats.length ? `
    <div class="audit-block audit-caveats-block">
      <h4>Caveats</h4>
      <ul class="audit-caveat-list">
        ${caveats.map(c => `
          <li class="audit-caveat sev-${c.severity || 'medium'}">
            <span class="caveat-sev sev-${c.severity || 'medium'}">${escapeHTML(c.severity || 'medium')}</span>
            <span class="caveat-summary">${escapeHTML(c.summary || '')}</span>
            ${c.detail ? `<div class="caveat-detail">${escapeHTML(c.detail).replace(/\n/g, "<br>")}</div>` : ""}
          </li>`).join("")}
      </ul>
    </div>` : "";

  // Dated alternatives — for any partition whose canonical is undated, or
  // where dated alternatives exist as sub-clades / candidates. Includes the
  // 'reason_canonical_undated' field that explains why dating wasn't done.
  const da = audit.dated_alternatives || {};
  const daCands = da.candidates || [];
  const daHTML = (da.reason_canonical_undated || daCands.length) ? `
    <div class="audit-block audit-dated-alt">
      <h4>Dating context${daCands.length ? ` · ${daCands.length} alternative${daCands.length===1?'':'s'}` : ''}</h4>
      ${da.reason_canonical_undated ? `<p class="da-reason">${escapeHTML(da.reason_canonical_undated)}</p>` : ""}
      ${daCands.length ? `
        <ul class="da-cands">
          ${daCands.map(c => `
            <li>
              <div class="da-cand-head">
                <strong>${escapeHTML(c.key || '')}</strong>
                ${c.scope ? `<span class="src">${escapeHTML(c.scope)}</span>` : ""}
                ${c.tips ? `<span class="src">${escapeHTML(String(c.tips))} tips</span>` : ""}
              </div>
              ${c.study ? `<div class="da-study">${escapeHTML(c.study)}</div>` : ""}
              ${c.note ? `<div class="da-note">${escapeHTML(c.note)}</div>` : ""}
              ${c.doi ? `<div class="da-doi"><a href="https://doi.org/${escapeHTML(c.doi)}" target="_blank" rel="noopener">doi.org/${escapeHTML(c.doi)}</a></div>` : ""}
            </li>`).join("")}
        </ul>` : ""}
    </div>` : "";

  // Superseded canonical — when an override replaced an earlier shipped tree
  // (e.g., diatoms: Nakov 2018 replaced Parks 2018). Show the demoted study
  // so users see the lineage decision.
  const sup = (audit.tree && audit.tree.superseded_canonical) || null;
  const supHTML = sup ? `
    <div class="audit-block audit-superseded">
      <h4>Superseded by current canonical</h4>
      <div class="da-cand-head">
        <strong>${escapeHTML(sup.key || '')}</strong>
        ${sup.doi ? `<a class="da-doi-inline" href="https://doi.org/${escapeHTML(sup.doi)}" target="_blank" rel="noopener">doi.org/${escapeHTML(sup.doi)}</a>` : ""}
      </div>
      ${sup.study ? `<div class="da-study">${escapeHTML(sup.study)}</div>` : ""}
      ${sup.note ? `<div class="da-note">${escapeHTML(sup.note)}</div>` : ""}
    </div>` : "";

  // Resolutions — issues that were closed during auditing.
  const resolutions = audit.resolutions || [];
  const resHTML = resolutions.length ? `
    <div class="audit-block audit-resolutions">
      <h4>Resolved issues (${resolutions.length})</h4>
      <ul class="audit-res-list">
        ${resolutions.map(r => `
          <li>
            <div><strong>Issue:</strong> ${escapeHTML(r.issue || '')}</div>
            <div class="res-fix"><strong>Resolved:</strong> ${escapeHTML(r.resolved_by || '')}</div>
          </li>`).join("")}
      </ul>
    </div>` : "";

  // Show the panel when there's ANY notable audit content to surface — not
  // just DOI corrections. Verified status alone is no reason to hide caveats
  // or dating context that the curator went to the trouble of documenting.
  const hasContent = !!(estSrcHTML || caveatsHTML || daHTML || supHTML || resHTML);
  if (!hasContent && status === "verified") {
    host.hidden = true;
    host.innerHTML = "";
    return;
  }

  host.hidden = false;
  host.innerHTML = `
    <div class="audit-header">
      <span class="audit-badge audit-${status}">${status.replace(/_/g, " ")}</span>
      <span class="src">audit · last reviewed ${lastAudited || "—"}</span>
    </div>
    ${estSrcHTML}
    ${caveatsHTML}
    ${supHTML}
    ${daHTML}
    ${resHTML}
  `;
}

// Posterior samples / HPD trees / BEAST XML — surfaced as a card with
// click-through links to the original data host (Dryad, GitHub, VertLife).
// The atlas does not mirror these files; the host is canonical.
function renderUncertaintyPanel(t) {
  const host = document.getElementById("uncertainty-panel");
  if (!host) return;
  const audit = STATE.data.audits?.[t.partition_group];
  const unc = audit?.uncertainty;
  if (!unc || !t.is_partition_canonical) {
    host.hidden = true;
    host.innerHTML = "";
    return;
  }
  const files = unc.files || [];
  const fileRowsHTML = files.map(f => {
    const sizeText = f.size_bytes
      ? `${f.size_bytes >= 1048576 ? (f.size_bytes / 1048576).toFixed(1) + " MB" : (f.size_bytes / 1024).toFixed(0) + " KB"}`
      : f.format && /large/i.test(f.format) ? "large" : "—";
    const fmtBadge = f.format ? `<span class="unc-fmt">${escapeHTML(f.format)}</span>` : "";
    return `
      <li class="unc-file">
        <div class="unc-file-head">
          <a class="unc-link" href="${f.url}" target="_blank" rel="noopener">${escapeHTML(f.label || f.name)}</a>
          ${fmtBadge}
          <span class="unc-size">${sizeText}</span>
        </div>
        ${f.description ? `<div class="unc-desc">${escapeHTML(f.description).replace(/\n/g, "<br>")}</div>` : ""}
        ${f.name && f.name !== f.label ? `<div class="unc-fname"><code>${escapeHTML(f.name)}</code></div>` : ""}
      </li>`;
  }).join("");
  host.hidden = false;
  host.innerHTML = `
    <div class="card">
      <div class="card-header">
        <h3>Posterior &amp; uncertainty data</h3>
        ${unc.host_url ? `<a class="btn-secondary" href="${unc.host_url}" target="_blank" rel="noopener">Browse at ${escapeHTML(unc.host || "host")}</a>` : ""}
      </div>
      <p class="card-sub">
        The release redistributes the reference tree. The source's tree sets, node intervals and dating input stay at the data host below.
      </p>
      <ul class="unc-list">${fileRowsHTML}</ul>
      ${unc.note ? `<p class="unc-note">${escapeHTML(unc.note).replace(/\n/g, "<br>")}</p>` : ""}
    </div>`;
}

function escapeHTML(s) {
  return String(s).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[c]));
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

  const canonical = items.find(x => x.peer.is_partition_canonical) || null;
  const subclades = items.filter(x => !x.peer.is_partition_canonical);

  const renderBtn = ({ peer, tree }, role) => {
    const isActive = tree.filename === t.filename;
    const tips = peer.ntips ? `${fmt.format(peer.ntips)} tips` : "";
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
    <div class="chooser-label">Group <strong>${partition}</strong></div>
    ${canonical ? `
    <div class="chooser-tier">
      <div class="tier-label">Reference tree</div>
      <div class="chooser-buttons">${renderBtn(canonical, "canonical")}</div>
    </div>` : ""}
    ${subclades.length ? `
    <div class="chooser-tier">
      <div class="tier-label">Other trees (not counted in the coverage) <span class="tier-count">${subclades.length}</span></div>
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
  STATE.currentNames = {};

  status.classList.remove("warning");
  status.textContent = `Loading ${t.filename}…`;

  try {
    const treeUrl = treeURL(t);
    const nTips = t.ntips || 0;
    const [treeRes, namesRes] = await Promise.all([
      fetch(treeUrl),
      t.names_file ? fetch(NAMES_BASE + t.names_file + ".json").catch(() => null) : Promise.resolve(null),
    ]);
    if (!treeRes.ok) throw new Error(`HTTP ${treeRes.status}`);
    let newick = await treeRes.text();
    if (namesRes && namesRes.ok) {
      try { STATE.currentNames = await namesRes.json(); } catch { STATE.currentNames = {}; }
    }

    if (nTips > TIP_RENDER_LIMIT) {
      const sampled = sampleNewick(newick, TIP_RENDER_LIMIT);
      newick = sampled.newick;
      status.classList.add("warning");
      status.textContent = `Tree has ${t.ntips_text} tips — showing a random sample of ${fmt.format(sampled.kept)}. Download for the full tree.`;
    } else {
      status.textContent = `${t.ntips_text} tips — full tree shown.`;
    }

    STATE.currentNewick = newick;
    STATE.currentDated = !!t.dated;
    drawTree(newick);
  } catch (e) {
    status.classList.add("warning");
    status.textContent = `Could not render: ${e.message}. Use the download link above.`;
  }
}

// Resolve a tip's numeric ID to its standardized name. Returns the raw label
// when no shard entry exists.
function tipDisplayName(rawName) {
  const n = STATE.currentNames[rawName];
  return n || rawName;
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
  const axisHeight = STATE.currentDated ? 30 : 0;
  const margin = { top: 12, right: labelWidth, bottom: 12 + axisHeight, left: 12 };
  const innerW = width - margin.left - margin.right;
  const tipSpacing = Math.max(1.2, Math.min(16, 600 / Math.max(40, nLeaves)));
  const innerH = Math.max(400, nLeaves * tipSpacing);

  // Use d3.cluster() to assign x (vertical position, evenly-spaced leaves).
  // Then OVERRIDE the y coordinate (horizontal, depth) using cumulative
  // branch lengths when the tree is dated — that's what turns a cladogram
  // into a time-tree.
  const cluster = d3.cluster()
    .size([innerH, innerW])
    .separation(() => 1);
  cluster(root);

  // Compute root-to-node cumulative branch length. Honor only when the
  // tree is flagged dated AND has non-trivial branch lengths.
  let maxRootDist = 0;
  root.each(d => {
    const blen = (d.data && d.data.length != null && !Number.isNaN(d.data.length))
      ? Math.max(0, d.data.length) : 0;
    d.rootDist = (d.parent ? d.parent.rootDist : 0) + blen;
    if (d.rootDist > maxRootDist) maxRootDist = d.rootDist;
  });
  const honorBranchLengths = STATE.currentDated && maxRootDist > 0;
  if (honorBranchLengths) {
    const yScale = innerW / maxRootDist;
    root.each(d => { d.y = d.rootDist * yScale; });
  }

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

  // Invisible hover targets at each tip so the cursor can pick up species
  // names even when labels are off. Native SVG <title> drives the tooltip.
  const tipHit = g.append("g").attr("class", "tip-hits");
  tipHit.selectAll("circle")
    .data(root.leaves())
    .join("circle")
      .attr("cx", d => d.y)
      .attr("cy", d => d.x)
      .attr("r", Math.max(2.5, tipSpacing * 0.6))
      .attr("fill", "transparent")
      .style("pointer-events", "all")
    .append("title")
      .text(d => tipDisplayName(d.data.name));

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
        .text(d => tipDisplayName(d.data.name))
      .append("title")
        .text(d => tipDisplayName(d.data.name));
  }

  // Time axis — only for dated trees with real branch lengths.
  if (honorBranchLengths) {
    const axisY = innerH + 4;
    const axisG = g.append("g").attr("class", "tree-time-axis")
                                .attr("transform", `translate(0, ${axisY})`);
    // x maps node.y back to root distance. Root is at y=0 (left, oldest),
    // tips at y=innerW (right, present). Age-from-tip axis: tip = 0,
    // root = maxRootDist (decreasing rightward).
    const xToAge = (y) => maxRootDist - (y / innerW) * maxRootDist;
    const niceTicks = (max) => {
      // Pick ~6 ticks at human-friendly steps.
      const steps = [1, 2, 5, 10, 20, 25, 50, 100, 200, 250, 500, 1000];
      const target = max / 6;
      const step = steps.find(s => s >= target) || Math.ceil(target / 100) * 100;
      const out = [];
      for (let v = 0; v <= max; v += step) out.push(v);
      return out;
    };
    const ticks = niceTicks(maxRootDist);
    // Axis line
    axisG.append("line")
      .attr("x1", 0).attr("x2", innerW)
      .attr("y1", 0).attr("y2", 0)
      .attr("stroke", "#888").attr("stroke-width", 0.8);
    // Tick marks + labels
    const tickG = axisG.selectAll("g.tick")
      .data(ticks)
      .join("g")
        .attr("class", "tick")
        .attr("transform", age => `translate(${innerW - (age / maxRootDist) * innerW}, 0)`);
    tickG.append("line")
      .attr("y1", 0).attr("y2", 5)
      .attr("stroke", "#888").attr("stroke-width", 0.8);
    tickG.append("text")
      .attr("y", 16)
      .attr("text-anchor", "middle")
      .attr("font-size", 10)
      .attr("fill", "#555")
      .text(age => age === 0 ? "0" : age);
    // Axis caption. Branch-length units aren't recorded per tree in the
    // atlas pipeline — most chronograms are in Ma but some standardized
    // files (e.g. turtles, max ~10) appear to carry rate-normalized
    // distances. Heuristic: max > 30 → likely Ma; otherwise generic.
    const captionText = maxRootDist > 30
      ? "Time before present (Ma)"
      : "Root-to-tip distance (branch-length units; not all source trees publish in Ma)";
    axisG.append("text")
      .attr("class", "tree-time-axis-caption")
      .attr("x", innerW / 2)
      .attr("y", 28)
      .attr("text-anchor", "middle")
      .attr("font-size", 10)
      .attr("fill", "#777")
      .text(captionText);
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

  // Hover targets — circles at each leaf so users can scrub for names.
  const hits = svg.append("g").attr("class", "tip-hits");
  hits.selectAll("circle")
    .data(root.leaves())
    .join("circle")
      .attr("cx", d => Math.cos(d.x - Math.PI / 2) * d.y)
      .attr("cy", d => Math.sin(d.x - Math.PI / 2) * d.y)
      .attr("r", 3)
      .attr("fill", "transparent")
      .style("pointer-events", "all")
    .append("title")
      .text(d => tipDisplayName(d.data.name));

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
        .text(d => tipDisplayName(d.data.name))
      .append("title")
        .text(d => tipDisplayName(d.data.name));
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
