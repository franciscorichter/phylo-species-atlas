// Phylo-Species Atlas browser — vanilla JS, no build step.
// Loads data.json, renders summary + plots, and lazy-loads tree files.

const STATE = {
  data: null,
  filtered: [],
  selected: null,
  layout: "linear", // or "radial"
  showLabels: false,
  currentNewick: null,
  currentNames: {}, // {id: standardized_name} for the active tree
  query: "",
  datedOnly: false,
};

const TREES_BASE = "../standardized/trees/"; // patched by CI for Pages deploy
const NAMES_BASE = "../standardized/names/"; // patched by CI for Pages deploy
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
  const auditedFrac = s.n_audited && s.n_groups
    ? `${s.n_audits_verified || 0}/${s.n_audited}`
    : null;
  const items = [
    ["Trees", fmt.format(s.n_trees)],
    ["Datasets", fmt.format(s.n_datasets)],
    ["Dated", fmt.format(s.n_dated)],
    ["Undated", fmt.format(s.n_undated)],
    ["Total tips", fmt.format(s.total_tips)],
  ];
  if (auditedFrac) {
    items.push([`<a href="audit.html" class="audit-link">Audits verified</a>`, auditedFrac]);
  }
  document.getElementById("summary-stats").innerHTML = items.map(([label, value]) =>
    `<div class="stat"><span class="value">${value}</span><span class="label">${label}</span></div>`
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

// Macro-category ordering matches the paper's Figure 1 narrative.
const CATEGORY_ORDER = ["Vertebrates", "Plants", "Arthropods", "Other animals", "Microbes & protists"];
const CATEGORY_RANK = Object.fromEntries(CATEGORY_ORDER.map((c, i) => [c, i]));

function renderCoveragePlot() {
  const allRows = STATE.data.coverage;
  if (!allRows.length) return;

  // Single chart per partition (canonical tree row). Three visual elements:
  //   ▸ filled marker at tree-tip count (what the atlas actually ships)
  //   ▸ open marker at described-species count (what taxonomists have catalogued)
  //   ▸ a line interval from estimated_low to estimated_high with a marker at
  //     estimated_total (what diversity studies estimate is really out there)
  // X-axis is log-scale species counts because the range across partitions
  // spans 6+ orders of magnitude (60 echinoderm tips to 4 million bacteria).
  const canonical = allRows.filter(r => r.is_partition_canonical);
  const rows = canonical.slice().sort((a, b) => {
    const ra = CATEGORY_RANK[a.category] ?? 99;
    const rb = CATEGORY_RANK[b.category] ?? 99;
    if (ra !== rb) return ra - rb;
    return (b.estimated_total ?? b.described_species ?? 0) - (a.estimated_total ?? a.described_species ?? 0);
  });

  const labels = rows.map(r => r.group);

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

  // Interval line + estimated-total marker. Heuristic ranges get a muted
  // amber colour so the eye reads them as less authoritative than the
  // paper-published ones (Stork 2018 insects, Mora 2011 amphibians, etc).
  const intervalColor = (r) => {
    const cls = intervalClass(r);
    return (cls === "fully-heuristic" || cls === "partly-heuristic")
      ? "rgba(180,140,90,0.55)"
      : "#5a5a5a";
  };

  // Estimate is a RANGE, not a point. Render as a line segment from low to
  // high for each row — implemented via Plotly shapes (no markers at all).
  // A separate invisible scatter trace at the midpoint provides a hover
  // target so users can still inspect interval details.
  const intervalShapes = rows.map(r => {
    if (r.estimated_low == null || r.estimated_high == null) return null;
    if (r.estimated_low >= r.estimated_high) return null;  // degenerate
    return {
      type: "line",
      xref: "x", yref: "y",
      x0: r.estimated_low, x1: r.estimated_high,
      y0: r.group, y1: r.group,
      line: { color: intervalColor(r), width: 3 },
      layer: "above",
    };
  }).filter(Boolean);
  // Hover-only trace: invisible scatter at the geometric midpoint so users
  // can still hover the interval to see source/range info.
  const intervalAnchor = (r) => {
    if (r.estimated_low && r.estimated_high) {
      return Math.sqrt(r.estimated_low * r.estimated_high);
    }
    return r.estimated_total ?? null;
  };
  const estimatedHoverTrace = {
    type: "scatter",
    mode: "markers",
    name: "Estimated true diversity (interval)",
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
      "Estimated range: %{customdata[0]:,}–%{customdata[1]:,}<br>" +
      (`<i style='color:#7a7a7a'>cited central: %{customdata[5]:,} · %{customdata[4]}</i><br>` +
       "Source: %{customdata[2]} · confidence %{customdata[3]}<extra></extra>"),
  };

  // Mute described/estimated markers for unrepresented partitions so the
  // missing-tree story reads at a glance (no ◇ tip + slightly faded ◯ ●).
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
      r.described_species ?? 0, r.described_source || "—",
      r.is_unrepresented ? " · no shipped tree" : "",
    ]),
    hovertemplate:
      "<b>%{y}</b>%{customdata[2]}<br>" +
      "Described species: %{x:,}<br>" +
      "Source: %{customdata[1]}<extra></extra>",
  };

  // Tip marker: ◇ for paper-shipped trees, ★ for atlas-derived species-level
  // pruning (Turtles, Bryozoa, Seed plants, Fungi, ...). The atlas-derived
  // ones are NOT what the source paper published verbatim — they're one-rep-
  // per-species reductions of multi-individual MCCs, transparent in info.yaml.
  // Build.py exposes this via the `tree_url_override` field per partition;
  // the chart picks up the override status by joining via filename.
  const overrideByFilename = {};
  for (const t of STATE.data.trees) {
    if (t.tree_url_override) overrideByFilename[t.filename] = true;
  }
  const isDerived = (r) => !!overrideByFilename[r.filename];

  const tipsTrace = {
    type: "scatter",
    mode: "markers",
    name: "Tips in shipped tree (◇ paper-shipped / ★ atlas-derived)",
    x: rows.map(r => r.tips ?? null),
    y: labels,
    marker: {
      symbol: rows.map(r => isDerived(r) ? "star" : "diamond"),
      size: rows.map(r => isDerived(r) ? 12 : 10),
      color: rows.map(r => r.dated ? "#2a6fbf" : "#c97a2a"),
      line: { color: "#1a1a1a", width: 0.5 },
    },
    customdata: rows.map(r => [
      r.tips ?? 0, r.described_species ?? 0,
      r.described_species ? (100 * r.tips / r.described_species).toFixed(1) : "—",
      r.study || "", r.year || "",
      isDerived(r) ? `atlas-derived species-level pruning (raw paper tree had ${fmt.format(r.raw_tips || r.tips)} tips)` : "paper-shipped MCC",
    ]),
    hovertemplate:
      "<b>%{y}</b><br>" +
      "Tips in tree: %{x:,}<br>" +
      "Coverage: %{customdata[2]}% of described<br>" +
      "<i style='color:#7a7a7a'>%{customdata[5]}</i><br>" +
      "<i>%{customdata[3]}</i> (%{customdata[4]})<extra></extra>",
  };

  // Category bands → labels in the left gutter + dotted dividers between groups.
  const bands = [];
  for (let i = 0; i < rows.length; i++) {
    const cat = rows[i].category || "Other";
    const last = bands[bands.length - 1];
    if (!last || last.category !== cat) bands.push({ category: cat, start: i, end: i });
    else last.end = i;
  }
  const LEFT_MARGIN = 230;
  const annotations = bands.map(b => {
    const midRow = Math.floor((b.start + b.end) / 2);
    return {
      xref: "paper", yref: "y",
      x: 0, xanchor: "left", xshift: -LEFT_MARGIN + 6,
      y: rows[midRow].group, yanchor: "middle",
      text: `<b>${b.category.toUpperCase()}</b>`,
      showarrow: false, align: "left",
      font: { family: FONT_STACK, size: 10.5, color: "#7a7a7a" },
    };
  });
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

  // Log-scale x-axis: range from a sensible floor to just above the highest
  // upper bound across all partitions.
  const allXVals = rows.flatMap(r => [r.tips, r.described_species, r.estimated_high])
    .filter(v => v != null && v > 0);
  const maxX = Math.pow(10, Math.ceil(Math.log10(Math.max(...allXVals))));
  const minX = 10;

  const layout = {
    margin: { l: LEFT_MARGIN, r: 30, t: 18, b: 64 },
    annotations,
    shapes,
    font: { family: FONT_STACK, size: 12, color: "#1a1a1a" },
    xaxis: {
      title: { text: "Species count (log scale) — ◇ tips · ◯ described · ● estimated total (with low–high whiskers)", font: { size: 12, color: "#6b6b6b" }, standoff: 14 },
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
  Plotly.newPlot(container, [estimatedHoverTrace, describedTrace, tipsTrace], layout,
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
    if (!row.filename) return;
    selectTree(row.filename, true);
  });

  const sub = document.getElementById("coverage-sub");
  if (sub) {
    const nShipped = rows.filter(r => !r.is_unrepresented).length;
    const nDark = rows.filter(r => r.is_unrepresented).length;
    sub.innerHTML = `Each row = one partition (${nShipped} with a shipped tree, ${nDark} 'dark matter' — no tree yet). <strong>◇</strong> = paper-shipped tree tips, <strong>★</strong> = atlas-derived species-level pruning (dated <span style='color:#2a6fbf'>■</span> / undated <span style='color:#c97a2a'>■</span>). <strong>◯</strong> = described species (taxonomic catalogue). <strong>━━━━</strong> = estimated true diversity interval (low–high); <em>solid grey</em> = paper-published range, <em>faded amber</em> = atlas-derived heuristic. Click a tree row to inspect it; click a dark-matter row to see the evidence behind its interval.`;
  }
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

    // Order partitions: by macro category (CATEGORY_ORDER), then by canonical
    // tips desc within each category. "Condamine 2019 families" and "Other"
    // are pinned to the bottom of their category.
    const partitionInfo = Array.from(buckets.keys()).map(p => {
      const items = buckets.get(p);
      const canonical = items.find(t => t.is_partition_canonical) || items[0];
      const category = items.find(t => t.category)?.category || "Other";
      return { partition: p, canonicalTips: canonical?.ntips || 0, category };
    });
    partitionInfo.sort((a, b) => {
      const ra = CATEGORY_RANK[a.category] ?? 99;
      const rb = CATEGORY_RANK[b.category] ?? 99;
      if (ra !== rb) return ra - rb;
      const pinA = a.partition === "Condamine 2019 families" || a.partition === "Other" ? 1 : 0;
      const pinB = b.partition === "Condamine 2019 families" || b.partition === "Other" ? 1 : 0;
      if (pinA !== pinB) return pinA - pinB;
      return b.canonicalTips - a.canonicalTips;
    });

    let lastCategory = null;
    html = partitionInfo.map(({ partition, category }) => {
      const items = buckets.get(partition);
      const canonical = items.find(t => t.is_partition_canonical);
      const subclades = items.filter(t => !t.is_partition_canonical);
      const hasSub = subclades.length > 0;
      const expanded = STATE_EXPANDED.has(partition);
      const chevron = hasSub ? `<span class="chev">${expanded ? "▾" : "▸"}</span>` : `<span class="chev empty"></span>`;
      const countBadge = hasSub ? `<span class="sub-count">+${subclades.length}</span>` : "";

      let categoryRow = "";
      if (category !== lastCategory) {
        categoryRow = `<li class="category-header">${category}</li>`;
        lastCategory = category;
      }

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
      return categoryRow + headerRow + canonicalRow + childRows;
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
  const rangeText = (est.low != null && est.high != null)
    ? `${fmt.format(est.low)}–${fmt.format(est.high)}`
    : (est.total != null ? `${fmt.format(est.total)} (single point)` : "—");

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
        <h3>No shipped tree for this partition</h3>
        <span class="audit-badge audit-deferred">${escapeHTML((audit.audit || {}).status || "deferred")}</span>
      </div>
      <p class="card-sub">${escapeHTML(audit.category || "")} · This partition lacks a canonical phylogenetic tree in the atlas. The information below describes the diversity estimate and its evidence.</p>

      ${canTree.reason ? `<section class="dm-section">
        <h3>Why no tree</h3>
        <p class="dm-reason">${escapeHTML(canTree.reason)}</p>
      </section>` : ""}

      <section class="dm-section">
        <h3>Diversity estimate</h3>
        <dl class="dm-dl">
          <dt>Described species</dt>
          <dd>${est.described != null ? fmt.format(est.described) : "—"}</dd>
          <dt>Estimated range</dt>
          <dd>${rangeText}</dd>
          ${est.confidence ? `<dt>Confidence</dt><dd>${escapeHTML(est.confidence)}</dd>` : ""}
          <dt>Evidence</dt>
          <dd>${srcLinks || '<span class="src">—</span>'}</dd>
        </dl>
      </section>

      ${constituentsHTML}
      ${intervalProvHTML}

      <p class="dm-footer"><a href="audit.html#${encodeURIComponent(partitionName)}" target="_blank" rel="noopener">Full audit entry →</a></p>
    </div>
  `;

  if (scrollIntoView) {
    document.getElementById("detail").scrollIntoView({ block: "start", behavior: "smooth" });
  }
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
  // Hide dark-matter panel; show tree-related cards.
  document.getElementById("dark-matter-panel").hidden = true;
  document.getElementById("r-snippet-card").hidden = false;
  document.getElementById("tree-viewer-card").hidden = false;

  const est = t.estimate;
  const describedVal = t.described_species
    ? `${fmt.format(t.described_species)}${sourceLinkHTML(t.described_source_info, t.described_source)}`
    : "—";
  const estimatedVal = est && est.estimated_total
    ? `${fmt.format(est.estimated_total)} <span class="range">(${fmt.format(est.estimated_low || est.estimated_total)}–${fmt.format(est.estimated_high || est.estimated_total)})</span>` +
      sourceLinkHTML(t.estimate_source_info, est.estimate_source, est.estimate_confidence)
    : null;
  // Coverage (described). Use the override's species count when an atlas-derived
  // canonical tree exists; otherwise check tip_taxonomy from the original shard.
  // Coverage uses the species-level effective count when available
  // (tip_taxonomy_override.unique_species for atlas-derived canonicals like
  // turtles; tip_taxonomy.unique_species when the tree has redundancy). Falls
  // back to raw ntips when neither applies. This matches build.py's coverage_pct.
  const txForCoverage = t.tip_taxonomy_override || t.tip_taxonomy;
  const effectiveCount = (txForCoverage?.unique_species) ?? t.ntips;

  let covDescribed = t.coverage_pct != null
    ? `${t.coverage_pct.toFixed(1)}% <span class="src">of described</span>`
    : "—";
  if (effectiveCount !== t.ntips && t.described_species) {
    covDescribed = `${t.coverage_pct.toFixed(1)}% <span class="src">at species level (${fmt.format(effectiveCount)} / ${fmt.format(t.described_species)} described)</span>`;
  }

  let covEstimated = null;
  if (est && est.estimated_total && effectiveCount) {
    const main = (100 * effectiveCount / est.estimated_total).toFixed(1);
    const lo = est.estimated_high ? (100 * effectiveCount / est.estimated_high).toFixed(1) : null;
    const hi = est.estimated_low ? (100 * effectiveCount / est.estimated_low).toFixed(1) : null;
    const subnote = t.is_partition_anchor
      ? `<span class="src">of ${t.partition_group || "partition"} estimate</span>`
      : `<span class="src">of ${t.partition_group} estimate &mdash; this tree is a sub-clade of ${t.partition_group}, so the % refers to the whole partition</span>`;
    covEstimated = `${main}% <span class="range">(${lo}&ndash;${hi}%)</span> ${subnote}`;
  }

  const treePaperHTML = t.doi
    ? `<a class="src-link" href="https://doi.org/${t.doi}" target="_blank" rel="noopener" title="${(t.study || "") + (t.year ? " (" + t.year + ")" : "") + (t.journal ? " — " + t.journal : "")}"><span class="src-type paper">paper</span>${(t.study || "").split(",")[0].split(" et")[0] || t.doi}${t.year ? " " + t.year : ""}</a>`
    : (t.study ? `<span class="src">${t.study}${t.year ? " (" + t.year + ")" : ""}</span>` : "—");

  // When an audit exists for the partition, prefer its authoritative figures
  // (n_loci, species_represented, multiple_individuals_per_species) over the
  // CSV-derived ones, and fold the methods block inline.
  const audit = STATE.data.audits?.[t.partition_group];
  const auditTree = audit?.tree;
  const auditMethods = auditTree?.methods;
  const useAuditMethods = !!auditMethods && t.is_partition_canonical;

  // Tips field. Layered fallback:
  //   1) tip_taxonomy_override (atlas-derived species-level tree present)
  //   2) tip_taxonomy (cross-partition, from the full name shard)
  //   3) audit-confirmed species_represented (from info.yaml)
  //   4) plain ntips
  let tipsVal = fmt.format(t.ntips || 0) + " tips";
  const txOv = t.tip_taxonomy_override;
  const tx = t.tip_taxonomy;
  if (txOv) {
    // Atlas-derived: report the override tree's actual content + flag the
    // derivation, with the original multi-individual count in parens.
    const parts = [`${fmt.format(txOv.total_tips)} tips`];
    if (txOv.unique_species != null) parts.push(`${fmt.format(txOv.unique_species)} species`);
    parts.push(`<span class="src">atlas-derived from ${fmt.format(t.ntips || 0)}-tip multi-individual MCC</span>`);
    tipsVal = parts.join(" · ");
  } else if (tx) {
    const parts = [`${fmt.format(t.ntips || 0)} tips`];
    if (tx.unique_ids != null && tx.unique_ids !== t.ntips) {
      parts.push(`${fmt.format(tx.unique_ids)} unique IDs`);
    }
    if (tx.unique_species != null && tx.unique_species !== (tx.unique_ids ?? t.ntips)) {
      parts.push(`${fmt.format(tx.unique_species)} species`);
    }
    if (tx.subspecies_tips > 0) {
      parts.push(`<span class="src">${fmt.format(tx.subspecies_tips)} subspecies-level</span>`);
    }
    tipsVal = parts.join(" · ");
  }
  if (useAuditMethods && auditTree.species_represented && !txOv && (!tx || tx.unique_species == null)) {
    tipsVal = `${fmt.format(t.ntips || 0)} tips · ${fmt.format(auditTree.species_represented)} species`;
  }

  const fields = [
    ["Group", t.group + (t.partition_group && t.partition_group !== t.group ? ` <span class="src">→ ${t.partition_group}</span>` : "")],
    ["Tree paper", treePaperHTML],
    ["Tips", tipsVal],
    ["Dated", t.dated ? "yes" : "no"],
    ["Crown age", t.crown_ma ? `${t.crown_ma} Ma` : "—"],
    ["Described species", describedVal],
    ["Estimated species", estimatedVal],
    ["Coverage (described)", covDescribed],
    ["Coverage (estimated)", covEstimated],
  ];

  // Inline methods + uncertainty block when audit available.
  if (useAuditMethods) {
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
      ["Support summary", m.support_typical],
      ["Divergence CIs", m.divergence_ci_available ? (m.divergence_ci_note || "yes") : "no"],
      ["Posterior trees", m.posterior_n ? fmt.format(m.posterior_n) : null],
    ];
    for (const [k, v] of methodFields) {
      if (v != null && v !== "") fields.push([k, escapeHTML(String(v))]);
    }
  }

  fields.push(
    ["Tree file source", t.data_source],
    ["File", t.size_bytes ? `${(t.size_bytes / 1024).toFixed(1)} KB` : "—"],
  );

  document.getElementById("detail-meta").innerHTML = fields
    .filter(([_, v]) => v !== null && v !== undefined && v !== "")
    .map(([k, v]) => `<div class="field"><span class="key">${k}</span><span class="val">${v}</span></div>`)
    .join("");

  renderDatasetChooser(t);
  renderUncertaintyPanel(t);

  const dl = document.getElementById("download-newick");
  dl.href = t.tree_url_override || (TREES_BASE + filename);
  dl.setAttribute("download", filename);

  renderRSnippet(t);

  await loadAndRenderTree(t);
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

  // If there's nothing notable beyond status (no DOI correction), keep the
  // panel hidden — the methods data is already in detail-meta.
  if (!estSrcHTML && status === "verified") {
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
        The atlas ships only the MCC consensus tree. The original posterior trees, HPD intervals, and BEAST input live at the data host below.
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
  STATE.currentNames = {};

  status.classList.remove("warning");
  status.textContent = `Loading ${t.filename}…`;

  try {
    // Prefer the atlas-derived override tree when build.py flags one (e.g.
    // turtles' species-level pruning). The name shard still resolves IDs
    // against the full dictionary, so it's the same shard.
    const treeUrl = t.tree_url_override
      ? t.tree_url_override
      : TREES_BASE + t.filename;
    // The effective tip count: if there's an override, use its taxonomy's
    // total_tips; otherwise the CSV-recorded ntips.
    const txOverride = t.tip_taxonomy_override;
    const effectiveTips = txOverride && txOverride.total_tips != null
      ? txOverride.total_tips
      : (t.ntips || 0);
    const [treeRes, namesRes] = await Promise.all([
      fetch(treeUrl),
      fetch(NAMES_BASE + t.filename + ".json").catch(() => null),
    ]);
    if (!treeRes.ok) throw new Error(`HTTP ${treeRes.status}`);
    let newick = await treeRes.text();
    if (namesRes && namesRes.ok) {
      try { STATE.currentNames = await namesRes.json(); } catch { STATE.currentNames = {}; }
    }

    if (effectiveTips > TIP_RENDER_LIMIT) {
      const sampled = sampleNewick(newick, TIP_RENDER_LIMIT);
      newick = sampled.newick;
      status.classList.add("warning");
      status.textContent = `Tree has ${fmt.format(effectiveTips)} tips — showing a random sample of ${fmt.format(sampled.kept)}. Download for the full tree.`;
    } else {
      const suffix = t.tree_url_override
        ? " — atlas-derived species-level tree (one tip per species)"
        : " — full tree shown.";
      status.textContent = `${fmt.format(effectiveTips)} tips${suffix}`;
    }

    STATE.currentNewick = newick;
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
