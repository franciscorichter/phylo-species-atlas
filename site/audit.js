// Audit log page — renders site/data.json's `audits` block as a per-partition
// verification trail. Reads the same data.json the main page uses.

const fmt = new Intl.NumberFormat("en-US");

function escapeHTML(s) {
  return String(s ?? "").replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;","'":"&#39;"}[c]));
}

const FILTERS = { status: "all", category: "all" };
let LOADED_DATA = null;

async function init() {
  const res = await fetch("data.json");
  if (!res.ok) {
    document.getElementById("audit-list").innerHTML = `<p class="empty">Could not load data.json (HTTP ${res.status}).</p>`;
    return;
  }
  LOADED_DATA = await res.json();
  renderSummary(LOADED_DATA);
  wireFilters();
  renderList(LOADED_DATA);
}

function wireFilters() {
  document.querySelectorAll("[data-filter-status]").forEach(btn => {
    btn.addEventListener("click", () => {
      FILTERS.status = btn.dataset.filterStatus;
      document.querySelectorAll("[data-filter-status]").forEach(b =>
        b.classList.toggle("active", b === btn));
      renderList(LOADED_DATA);
    });
  });
  document.querySelectorAll("[data-filter-cat]").forEach(btn => {
    btn.addEventListener("click", () => {
      FILTERS.category = btn.dataset.filterCat;
      document.querySelectorAll("[data-filter-cat]").forEach(b =>
        b.classList.toggle("active", b === btn));
      renderList(LOADED_DATA);
    });
  });
}

function renderSummary(data) {
  const audits = data.audits || {};
  const list = Object.values(audits);
  const stats = data.summary || {};
  const verified = list.filter(a => (a.audit || {}).status === "verified").length;
  const inProgress = list.filter(a => (a.audit || {}).status === "in_progress").length;
  const totalPartitions = (stats.n_groups || 0) + (data.unrepresented || []).length;
  const items = [
    ["Audited partitions", `${list.length} / ${totalPartitions || "—"}`],
    ["Verified", verified],
    ["In progress", inProgress],
    ["Total resolutions", list.reduce((acc, a) => acc + ((a.resolutions || []).length), 0)],
  ];
  document.getElementById("audit-summary").innerHTML = items
    .map(([k, v]) => `<div class="stat"><span class="stat-num">${v}</span><span class="stat-key">${k}</span></div>`)
    .join("");
}

function renderList(data) {
  const audits = data.audits || {};
  const allEntries = Object.entries(audits);
  if (!allEntries.length) {
    document.getElementById("audit-empty").hidden = false;
    return;
  }
  // Apply filters.
  let entries = allEntries.filter(([_, a]) => {
    if (FILTERS.status !== "all" && (a.audit || {}).status !== FILTERS.status) return false;
    if (FILTERS.category !== "all" && a.category !== FILTERS.category) return false;
    return true;
  });
  // Order: verified first, then in_progress, then alphabetic.
  const order = { verified: 0, in_progress: 1, shipped: 2, stub: 3, deferred: 4 };
  entries.sort(([na, a], [nb, b]) => {
    const sa = order[(a.audit || {}).status] ?? 9;
    const sb = order[(b.audit || {}).status] ?? 9;
    if (sa !== sb) return sa - sb;
    return na.localeCompare(nb);
  });

  const countEl = document.getElementById("audit-count");
  if (countEl) {
    countEl.textContent = entries.length === allEntries.length
      ? `Showing all ${allEntries.length} partitions`
      : `Showing ${entries.length} of ${allEntries.length} partitions`;
  }

  document.getElementById("audit-list").innerHTML = entries
    .map(([partition, a]) => renderPartition(partition, a))
    .join("");
}

function renderPartition(partition, audit) {
  const auditMeta = audit.audit || {};
  const status = auditMeta.status || "stub";
  const lastAudited = auditMeta.last_audited || "—";
  const category = audit.category || "";

  const estObj = audit.estimate || {};
  const estSrc = estObj.source || {};
  const ip = estObj.interval_provenance;

  let intervalHTML = "";
  if (ip) {
    const clsLabel = {
      "paper-published-range": "paper-published range",
      "partly-heuristic": "partly heuristic",
      "fully-heuristic": "fully heuristic",
      "derived": "derived from described count",
    }[ip.overall_classification] || ip.overall_classification;
    const clsClass = ip.overall_classification || "unaudited";
    const fieldRow = (label, sub) => {
      if (!sub) return "";
      const note = sub.note ? `<div class="ip-note">${escapeHTML(sub.note)}</div>` : "";
      const quote = sub.verbatim_quote ? `<div class="ip-quote">“${escapeHTML(sub.verbatim_quote)}”</div>` : "";
      return `<tr>
        <th>${escapeHTML(label)}</th>
        <td><span class="ip-type ip-type-${(sub.type || "").replace(/-/g, "_")}">${escapeHTML(sub.type || "—")}</span></td>
        <td>${escapeHTML(sub.source || "—")}${quote}${note}</td>
      </tr>`;
    };
    intervalHTML = `
      <section class="aud-block">
        <h3>Interval provenance <span class="ip-cls ip-cls-${clsClass}">${escapeHTML(clsLabel)}</span></h3>
        <table class="ip-table">
          ${fieldRow("Described", ip.described)}
          ${fieldRow("Estimated total", ip.estimated_total)}
          ${fieldRow("Low bound", ip.estimated_low)}
          ${fieldRow("High bound", ip.estimated_high)}
        </table>
        ${ip.auditor_note ? `<p class="aud-note">${escapeHTML(ip.auditor_note)}</p>` : ""}
      </section>`;
  }

  // Deferred reason: why this partition has no shipped tree (from canonical_tree.reason).
  let deferredHTML = "";
  const ct = audit.canonical_tree || {};
  if (auditMeta.status === "deferred" && ct.reason) {
    deferredHTML = `
      <section class="aud-block aud-deferred">
        <h3>No shipped canonical tree</h3>
        <p class="aud-deferred-reason">${escapeHTML(ct.reason)}</p>
      </section>`;
  }

  let doiFixHTML = "";
  if (estSrc.paper_doi_status === "broken" && estSrc.live_doi) {
    doiFixHTML = `
      <section class="aud-block">
        <h3>Estimate-source DOI correction</h3>
        <div class="aud-doi">
          <div class="doi-row doi-broken">
            <span class="doi-label">Paper cites</span>
            <a href="https://doi.org/${escapeHTML(estSrc.paper_doi)}" target="_blank" rel="noopener">${escapeHTML(estSrc.paper_key || "")} · ${escapeHTML(estSrc.paper_doi)}</a>
            <span class="badge-broken">404</span>
          </div>
          <div class="doi-row doi-live">
            <span class="doi-label">Live (replacement)</span>
            <a href="https://doi.org/${escapeHTML(estSrc.live_doi)}" target="_blank" rel="noopener">${escapeHTML(estSrc.live_key || "live")} · ${escapeHTML(estSrc.live_doi)}</a>
            ${estSrc.live_year ? `<span class="muted">${estSrc.live_year}</span>` : ""}
            ${estSrc.live_cited_value ? `<span class="muted">${fmt.format(estSrc.live_cited_value)} species</span>` : ""}
          </div>
        </div>
        ${estSrc.note ? `<p class="aud-note">${escapeHTML(estSrc.note).replace(/\n/g, "<br>")}</p>` : ""}
      </section>`;
  }

  const resolutions = audit.resolutions || [];
  const resHTML = resolutions.length
    ? `
      <section class="aud-block">
        <h3>Resolutions (${resolutions.length})</h3>
        <ol class="aud-resolutions">
          ${resolutions.map(r => `
            <li>
              <div class="res-issue"><strong>Issue:</strong> ${escapeHTML(r.issue)}</div>
              <div class="res-fix"><strong>Resolved by:</strong> ${escapeHTML(r.resolved_by)}</div>
            </li>`).join("")}
        </ol>
      </section>`
    : "";

  const cands = audit.candidate_subclades || [];
  const candHTML = cands.length
    ? `
      <section class="aud-block">
        <h3>Candidate sub-clades (${cands.length})</h3>
        <ul class="aud-cands">
          ${cands.map(c => `
            <li>
              <div class="cand-head">
                <strong>${escapeHTML(c.taxon || "")}</strong>
                <span class="cand-priority cand-${c.priority || "low"}">${escapeHTML(c.priority || "low")}</span>
                <span class="cand-status">${escapeHTML(c.status || "proposed")}</span>
              </div>
              ${c.rationale ? `<p class="cand-rationale">${escapeHTML(c.rationale)}</p>` : ""}
              ${(c.candidate_papers || []).map(p => `
                <p class="cand-paper">
                  ${p.doi ? `<a href="https://doi.org/${escapeHTML(p.doi)}" target="_blank" rel="noopener">${escapeHTML(p.doi)}</a>` : ""}
                  ${p.year ? `<span class="muted">${p.year}</span>` : ""}
                  · ${escapeHTML(p.study || "")}
                  ${p.ntips ? `<span class="muted">${fmt.format(p.ntips)} tips</span>` : ""}
                </p>`).join("")}
            </li>`).join("")}
        </ul>
      </section>`
    : "";

  const shipped = audit.shipped_subclades || [];
  const shipHTML = shipped.length
    ? `
      <section class="aud-block">
        <h3>Shipped sub-clades (${shipped.length})</h3>
        <ul class="aud-cands">
          ${shipped.map(s => `
            <li>
              <strong>${escapeHTML(s.taxon || "")}</strong>
              <span class="muted">${escapeHTML(s.study || "")}</span>
              ${s.ntips ? `<span class="muted">${fmt.format(s.ntips)} tips</span>` : ""}
            </li>`).join("")}
        </ul>
      </section>`
    : "";

  return `
    <article class="aud-card">
      <header class="aud-card-head">
        <h2>${escapeHTML(partition)}</h2>
        <div class="aud-meta">
          <span class="audit-badge audit-${status}">${status.replace(/_/g, " ")}</span>
          ${category ? `<span class="muted">${escapeHTML(category)}</span>` : ""}
          <span class="muted">last reviewed ${escapeHTML(lastAudited)}</span>
        </div>
      </header>
      ${deferredHTML}
      ${intervalHTML}
      ${doiFixHTML}
      ${resHTML}
      ${shipHTML}
      ${candHTML}
    </article>`;
}

init();
