"""Regenerate the GitHub wiki pages from the canonical data sources.

One markdown page per atlas partition (62 total: 26 represented, 36
unrepresented), plus a Home.md index. Each page synthesises:

  - site/data/partitions/<slug>/info.yaml       (partition metadata)
  - audits/partitions/<slug>.md                 (audit verification frontmatter)
  - audits/table_s5_canonical_tree_provenance.csv (current atlas canonical)
  - audits/table_s7_canonical_succession.csv     (historical succession)
  - data_provenance.csv                         (download URLs, sub-clades)

The "Notes & open questions" section in each page is enclosed by HTML sentinel
comments and is preserved across regenerations: existing free-form human notes
stay in place, while the structural blocks above are overwritten with the latest
data.

Usage:
  python3 scripts/regen_wiki.py             # writes pages to wiki/
  python3 scripts/regen_wiki.py --target /tmp/wiki  # writes elsewhere

The default target is a `wiki/` sibling directory at the repo root. To deploy
to the live GitHub wiki, clone the wiki repo, point this script at it, and
push:

  git clone https://github.com/franciscorichter/phylo-species-atlas.wiki.git wiki
  python3 scripts/regen_wiki.py --target wiki
  cd wiki && git add . && git commit -m 'wiki sync from atlas vX.Y.Z' && git push
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is required (pip install pyyaml)", file=sys.stderr)
    sys.exit(1)

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent
PARTITIONS_DIR = REPO_ROOT / "site" / "data" / "partitions"
AUDITS_DIR = REPO_ROOT / "audits" / "partitions"
S5_PATH = REPO_ROOT / "audits" / "table_s5_canonical_tree_provenance.csv"
S7_PATH = REPO_ROOT / "audits" / "table_s7_canonical_succession.csv"
DATA_PROV_PATH = REPO_ROOT / "data_provenance.csv"

RELEASE = "v1.0.4"
RELEASE_DATE = "2026-05-31"

SENTINEL_START = "<!-- HUMAN-NOTES-START -->"
SENTINEL_END = "<!-- HUMAN-NOTES-END -->"
DEFAULT_NOTES_BLOCK = (
    "## Notes & open questions\n\n"
    "_This section is human-editable and preserved across page regenerations._\n\n"
    "_(no notes yet — add observations, flag stale references, or propose canonical updates here)_\n"
)


# ─────────────────────────── helpers ───────────────────────────

def load_audit_frontmatter(slug: str) -> dict:
    path = AUDITS_DIR / f"{slug}.md"
    if not path.exists():
        return {}
    text = path.read_text()
    m = re.match(r"^---\n(.*?)\n---\n", text, re.DOTALL)
    if not m:
        return {}
    try:
        return yaml.safe_load(m.group(1)) or {}
    except yaml.YAMLError:
        return {}


def load_info_yaml(slug: str) -> dict:
    path = PARTITIONS_DIR / slug / "info.yaml"
    if not path.exists():
        return {}
    try:
        return yaml.safe_load(path.read_text()) or {}
    except yaml.YAMLError:
        return {}


def load_csv_rows(path: pathlib.Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def s5_for_slug(s5_rows: list[dict], slug: str) -> dict | None:
    """The S5 row whose partition_slug matches; returns None if unrepresented."""
    for r in s5_rows:
        if r.get("partition_slug") == slug and slug != "(condamine_bulk)":
            return r
    return None


def s7_for_slug(s7_rows: list[dict], slug: str) -> list[dict]:
    return [r for r in s7_rows if r.get("partition_slug") == slug]


def page_filename(name: str) -> str:
    """GitHub wiki: spaces become dashes in URL; filenames use dashes."""
    safe = re.sub(r"[^A-Za-z0-9 \-]", "", name).strip()
    return safe.replace(" ", "-") + ".md"


def preserve_notes(existing_path: pathlib.Path) -> str:
    """Extract the human-editable Notes block from an existing wiki page."""
    if not existing_path.exists():
        return DEFAULT_NOTES_BLOCK
    text = existing_path.read_text()
    m = re.search(
        rf"{re.escape(SENTINEL_START)}(.*?){re.escape(SENTINEL_END)}",
        text, re.DOTALL,
    )
    if m:
        return m.group(1).strip("\n") + "\n"
    return DEFAULT_NOTES_BLOCK


def fmt_float(x, ndigits=2):
    try:
        v = float(x)
        return f"{v:.{ndigits}f}"
    except (TypeError, ValueError):
        return "—"


def fmt_int(x):
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return "—"


def doi_link(doi: str) -> str:
    if not doi:
        return "—"
    doi_clean = doi.strip()
    return f"[{doi_clean}](https://doi.org/{doi_clean})"


def safe(x, default="—"):
    if x is None or x == "":
        return default
    return str(x)


def inline(x, default="—"):
    """Collapse multi-line strings into a single line for markdown table cells."""
    if x is None or x == "":
        return default
    return re.sub(r"\s+", " ", str(x)).strip()


# ─────────────────────────── page templates ───────────────────────────

def render_represented(slug: str, info: dict, audit: dict,
                       s5_row: dict, s7_rows: list[dict],
                       data_prov: list[dict], notes_block: str) -> str:
    name = info.get("partition", slug.replace("_", " ").title())
    category = info.get("category", "—")
    estimate = info.get("estimate", {}) or {}
    interval = estimate.get("interval_provenance", {}) or {}
    audit_status = audit.get("status", "—")
    audit_date = audit.get("last_audited", "—")
    audit_who = audit.get("auditor", "—")

    # Coverage from S5
    s5_ntips = fmt_int(s5_row.get("ntips"))
    s5_dated = s5_row.get("dated", "no")
    s5_mol = s5_row.get("molecular_tip_fraction", "—")
    s5_doi = (s5_row.get("doi") or "").strip()
    s5_study = (s5_row.get("tree_label") or s5_row.get("citation_key") or "—").strip()
    s5_dating_code = s5_row.get("dating_code", "—")
    s5_otl = s5_row.get("otl_pipeline_involvement", "—")

    described = estimate.get("described")
    coverage = "—"
    try:
        if described and s5_row.get("ntips"):
            pct = 100.0 * int(s5_row["ntips"]) / int(described)
            pct = min(pct, 100.0)
            coverage = f"{pct:.1f}%"
    except (TypeError, ValueError, ZeroDivisionError):
        pass

    # ─── header ───
    out = [f"# {name}", ""]
    out.append(f"> **Atlas canonical** · {s5_study} · {s5_ntips} tips · "
               f"{'dated' if s5_dated.lower() == 'yes' else 'undated'} ({s5_dating_code}) · {coverage} coverage")
    out.append(f"> **Category** · {category}")
    out.append(f"> **Audit** · {audit_status} · last reviewed {audit_date} by {audit_who}")
    out.append(f"> **Release** · {RELEASE} ([{RELEASE_DATE}]("
               f"https://github.com/franciscorichter/phylo-species-atlas/releases/tag/{RELEASE}))")
    out.append("")
    out.append("This page is the living review surface for the *" + name + "* "
               "partition of the Phylo-Species Atlas. The structured blocks below "
               "are regenerated from `info.yaml`, the per-partition audit, and "
               "Supplementary Tables S5/S7 at every release; the **Notes & open "
               "questions** section at the bottom is human-editable and preserved "
               "across regenerations.")
    out.append("")
    out.append("---")
    out.append("")

    # ─── diversity estimate ───
    out.append("## Diversity estimate")
    out.append("")
    src = estimate.get("source", {}) or {}
    src_key = src.get("key", "—")
    src_doi = src.get("doi", "")
    out.append("| Field | Value | Provenance |")
    out.append("|---|---|---|")
    out.append(f"| Described species | **{fmt_int(estimate.get('described'))}** | "
               f"{inline((interval.get('described', {}) or {}).get('source'))} |")
    out.append(f"| Range (low–high) | {fmt_int(estimate.get('low'))} – "
               f"{fmt_int(estimate.get('high'))} | "
               f"{inline((interval.get('estimated_total', {}) or {}).get('source'))} |")
    out.append(f"| Total (central) | {fmt_int(estimate.get('total'))} | — |")
    out.append(f"| Confidence | {safe(estimate.get('confidence'))} | "
               f"{inline((interval.get('overall_classification') or '—').replace('_',' '))} |")
    src_note = src.get("note") or ("DOI verified per audit"
        if (audit.get("estimate_source") or {}).get("doi_verified") else "—")
    out.append(f"| Source | {inline(src_key)} {(' · ' + doi_link(src_doi)) if src_doi else ''} | "
               f"{inline(src_note)} |")
    auditor_note = inline(interval.get("auditor_note"), default="")
    if auditor_note:
        out.append("")
        out.append(f"> **Audit caveat.** {auditor_note}")
    out.append("")
    out.append("---")
    out.append("")

    # ─── atlas-canonical tree (S5) ───
    out.append("## Atlas-canonical tree (Supplementary Table S5)")
    out.append("")
    out.append("| Field | Value |")
    out.append("|---|---|")
    out.append(f"| Study | **{s5_study}** |")
    out.append(f"| DOI | {doi_link(s5_doi)} |")
    out.append(f"| Tips | {s5_ntips} |")
    out.append(f"| Dated | {s5_dated} |")
    out.append(f"| `molecular_tip_fraction` | **{s5_mol}** |")
    out.append(f"| `dating_code` | `{s5_dating_code}` |")
    out.append(f"| `otl_pipeline_involvement` | `{s5_otl}` |")
    out.append(f"| Atlas coverage (C0) | {coverage} |")

    # data_provenance hints (download URL)
    dp_rows = [r for r in data_prov if r.get("group") == slug or r.get("group") == slug + "_jetz"]
    if dp_rows:
        for dp in dp_rows[:2]:
            durl = dp.get("download_url", "")
            if durl:
                out.append(f"| Repository | [{dp.get('data_source', 'source')}]({durl}) |")
                break
    out.append("")
    out.append("---")
    out.append("")

    # ─── canonical succession (S7) ───
    if s7_rows:
        out.append("## Canonical succession (Supplementary Table S7)")
        out.append("")
        out.append("The historical sequence of trees that have served as the "
                   "canonical for this partition between 2006 and the present. "
                   "The **atlas canonical** is the row flagged "
                   "`canonical_in_atlas=yes`; the others are predecessors that "
                   "were superseded, or recent literature that the atlas has "
                   "not (yet) adopted.")
        out.append("")
        out.append("| Year | Study | DOI | Tips | Dated | Mol. frac | Dating | Conf. |")
        out.append("|---|---|---|---:|:-:|---:|:-:|:-:|")
        for r in sorted(s7_rows, key=lambda x: int(x.get("year") or 0)):
            year = r.get("year", "?")
            label = r.get("study_label", "?")
            doi = r.get("doi", "")
            tips = fmt_int(r.get("ntips"))
            dated = "✓" if (r.get("dated") or "").lower() == "yes" else "—"
            mol = fmt_float(r.get("molecular_tip_fraction"))
            dc = r.get("dating_code", "—")
            conf = r.get("confidence", "—")
            is_atlas = (r.get("canonical_in_atlas") or "").lower() == "yes"
            if is_atlas:
                label = f"**{label} — atlas canonical**"
            year_str = f"**{year}**" if is_atlas else str(year)
            out.append(f"| {year_str} | {label} | {doi_link(doi)} | {tips} | "
                       f"{dated} | {mol} | {dc} | {conf} |")
        out.append("")
        out.append("---")
        out.append("")

    # ─── audit status ───
    out.append("## Audit status")
    out.append("")
    out.append(f"- **Audit version**: {safe(audit.get('audit_version'))}")
    out.append(f"- **Status**: {audit_status}")
    out.append(f"- **Last audited**: {audit_date} by {audit_who}")
    es = audit.get("estimate_source", {}) or {}
    if es.get("doi_verified"):
        out.append("- **Estimate source DOI verified**: ✓")
    if es.get("newer_authoritative_known"):
        nc = es.get("newer_candidate", {}) or {}
        rec = nc.get("recommendation", "watch")
        out.append(f"- **Newer estimate source**: {nc.get('study', '—')} "
                   f"(action: {rec})")
    ct = audit.get("canonical_tree", {}) or {}
    if ct.get("doi_verified"):
        out.append("- **Canonical tree DOI verified**: ✓")
    out.append(f"- **Full audit notes**: "
               f"[`audits/partitions/{slug}.md`]("
               f"https://github.com/franciscorichter/phylo-species-atlas/blob/"
               f"{RELEASE}/audits/partitions/{slug}.md)")
    out.append("")

    # ─── caveats from info.yaml ───
    caveats = info.get("caveats", []) or []
    if caveats:
        out.append("---")
        out.append("")
        out.append("## Open data caveats")
        out.append("")
        for c in caveats:
            sev = c.get("severity", "info").upper()
            summary = c.get("summary", "—")
            detail = (c.get("detail") or "").strip()
            out.append(f"- **[{sev}]** {summary}")
            if detail:
                for line in detail.split("\n"):
                    out.append(f"  > {line.strip()}")
        out.append("")

    # ─── dated alternatives (if canonical is undated) ───
    dated_alts = info.get("dated_alternatives") or {}
    if dated_alts:
        out.append("---")
        out.append("")
        out.append("## Dated-tree alternatives")
        out.append("")
        reason = (dated_alts.get("reason_canonical_undated") or "").strip()
        if reason:
            out.append(f"_{reason}_")
            out.append("")
        for cand in dated_alts.get("candidates", []) or []:
            study = cand.get("study", cand.get("key", "—"))
            cdoi = cand.get("doi", "")
            out.append(f"- **{cand.get('key', '—')}** · "
                       f"{doi_link(cdoi)} · scope: {cand.get('scope', '—')}, "
                       f"~{cand.get('tips', '?')} tips")
            note = (cand.get("note") or "").strip()
            if note:
                for line in note.split("\n"):
                    out.append(f"  > {line.strip()}")
        out.append("")

    # ─── notes (preserved) ───
    out.append("---")
    out.append("")
    out.append(SENTINEL_START)
    out.append("")
    out.append(notes_block.strip())
    out.append("")
    out.append(SENTINEL_END)
    out.append("")
    out.append("---")
    out.append("")
    out.append(
        f"*Page auto-generated from `site/data/partitions/{slug}/info.yaml`, "
        f"`audits/partitions/{slug}.md`, S5, and S7. Last regenerated: "
        f"{RELEASE_DATE} from release {RELEASE}.*"
    )
    out.append("")
    return "\n".join(out)


def render_unrepresented(slug: str, info: dict, audit: dict,
                         data_prov: list[dict], notes_block: str) -> str:
    name = info.get("partition", slug.replace("_", " ").title())
    category = info.get("category", "—")
    estimate = info.get("estimate", {}) or {}
    audit_status = audit.get("status", "—")
    audit_date = audit.get("last_audited", "—")

    out = [f"# {name} (unrepresented)", ""]
    out.append(f"> **Status** · No broadly representative species-level "
               f"phylogeny identified as of {RELEASE_DATE}")
    out.append(f"> **Category** · {category}")
    out.append(f"> **Audit** · {audit_status} · last reviewed {audit_date}")
    out.append(f"> **Release** · {RELEASE}")
    out.append("")
    out.append("This partition is part of the atlas's partition framework but "
               "currently has no canonical species-level phylogenetic tree "
               "ingested. The page documents diversity estimates, any partial "
               "sub-clade trees that exist in the literature, and what would "
               "be required for an atlas-canonical to land here.")
    out.append("")
    out.append("---")
    out.append("")

    out.append("## Diversity estimate")
    out.append("")
    out.append(f"- **Described species**: {fmt_int(estimate.get('described'))}")
    out.append(f"- **Range**: {fmt_int(estimate.get('low'))} – "
               f"{fmt_int(estimate.get('high'))}")
    src = estimate.get("source", {}) or {}
    if src:
        out.append(f"- **Source**: {src.get('key', '—')} {(' · ' + doi_link(src.get('doi', ''))) if src.get('doi') else ''}")
    out.append("")

    caveats = info.get("caveats", []) or []
    if caveats:
        out.append("## Why unrepresented")
        out.append("")
        for c in caveats:
            summary = c.get("summary", "—")
            detail = (c.get("detail") or "").strip()
            out.append(f"- {summary}")
            if detail:
                for line in detail.split("\n"):
                    out.append(f"  > {line.strip()}")
        out.append("")

    out.append("## Sub-clade or partial trees in the literature")
    out.append("")
    out.append("_Documented partial coverage. To propose elevating any of these "
               "to atlas-canonical, open an issue or PR against "
               "`site/data/partitions/" + slug + "/info.yaml`._")
    out.append("")

    out.append("---")
    out.append("")
    out.append(SENTINEL_START)
    out.append("")
    out.append(notes_block.strip())
    out.append("")
    out.append(SENTINEL_END)
    out.append("")
    out.append("---")
    out.append("")
    out.append(f"*Page auto-generated from `site/data/partitions/{slug}/info.yaml` "
               f"and `audits/partitions/{slug}.md`. Last regenerated: "
               f"{RELEASE_DATE} from release {RELEASE}.*")
    out.append("")
    return "\n".join(out)


def render_home(rep: list[dict], unrep: list[dict]) -> str:
    out = ["# Phylo-Species Atlas — Wiki", ""]
    out.append(
        "Living per-partition review surface for the "
        "[Phylo-Species Atlas](https://github.com/franciscorichter/phylo-species-atlas) "
        "(release " + RELEASE + ", " + RELEASE_DATE + ")."
    )
    out.append("")
    out.append("Each page below documents one *partition* — a taxonomic unit "
               "for which the atlas tracks the species-level phylogenetic "
               "tree (or its absence). Structured blocks on each page are "
               "regenerated from the repo's canonical data files at every "
               "release; the **Notes & open questions** section is human-"
               "editable and preserved across regenerations.")
    out.append("")
    out.append("## Represented partitions")
    out.append("")
    out.append("These partitions ship a canonical species-level tree in the atlas.")
    out.append("")
    out.append("| Partition | Atlas canonical | Tips | Dated | Coverage | Audit |")
    out.append("|---|---|---:|:-:|---:|:-:|")
    for r in sorted(rep, key=lambda x: x["name"]):
        out.append(
            f"| [{r['name']}]({r['link']}) | {r['study']} | {r['tips']} | "
            f"{r['dated']} | {r['coverage']} | {r['audit']} |"
        )
    out.append("")
    out.append("## Unrepresented partitions")
    out.append("")
    out.append("Partitions in the atlas framework that currently lack a "
               "canonical species-level tree. Their pages document diversity "
               "estimates, sub-clade or partial trees in the literature, and "
               "what's needed for a canonical to land.")
    out.append("")
    out.append("| Partition | Category | Described | Range |")
    out.append("|---|---|---:|---|")
    for r in sorted(unrep, key=lambda x: x["name"]):
        out.append(f"| [{r['name']}]({r['link']}) | {r['category']} | "
                   f"{r['described']} | {r['range']} |")
    out.append("")
    out.append("---")
    out.append("")
    out.append(f"*Wiki regenerated from release {RELEASE} on {RELEASE_DATE}.* "
               "See [scripts/regen_wiki.py](https://github.com/franciscorichter/"
               "phylo-species-atlas/blob/" + RELEASE + "/scripts/regen_wiki.py) "
               "for the generator.")
    return "\n".join(out)


# ─────────────────────────── driver ───────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default=str(REPO_ROOT.parent / "phylo-atlas-wiki-out"),
                    help="Directory to write pages into (defaults to a sibling dir)")
    args = ap.parse_args()

    target = pathlib.Path(args.target)
    target.mkdir(parents=True, exist_ok=True)

    s5_rows = load_csv_rows(S5_PATH)
    s7_rows = load_csv_rows(S7_PATH)
    data_prov = load_csv_rows(DATA_PROV_PATH)

    slugs = sorted(p.name for p in PARTITIONS_DIR.iterdir() if (p / "info.yaml").exists())

    rep_index, unrep_index = [], []
    n_rep = n_unrep = 0

    for slug in slugs:
        info = load_info_yaml(slug)
        audit = load_audit_frontmatter(slug)
        s5_row = s5_for_slug(s5_rows, slug)
        s7_for = s7_for_slug(s7_rows, slug)

        name = info.get("partition", slug.replace("_", " ").title())
        fname = page_filename(name)
        existing = target / fname
        notes = preserve_notes(existing)

        if s5_row:
            md = render_represented(slug, info, audit, s5_row, s7_for, data_prov, notes)
            n_rep += 1
            tips = fmt_int(s5_row.get("ntips"))
            dated = "✓" if (s5_row.get("dated") or "").lower() == "yes" else "—"
            described = (info.get("estimate") or {}).get("described")
            cov = "—"
            try:
                if described and s5_row.get("ntips"):
                    pct = 100.0 * int(s5_row["ntips"]) / int(described)
                    cov = f"{min(pct, 100.0):.1f}%"
            except (TypeError, ValueError, ZeroDivisionError):
                pass
            rep_index.append({
                "name": name,
                "link": fname.replace(".md", ""),  # GitHub wiki links don't use .md
                "study": (s5_row.get("tree_label") or "—").split("(")[0].strip().rstrip(","),
                "tips": tips,
                "dated": dated,
                "coverage": cov,
                "audit": "✓" if (audit.get("status") == "verified") else audit.get("status", "?"),
            })
        else:
            md = render_unrepresented(slug, info, audit, data_prov, notes)
            n_unrep += 1
            estimate = info.get("estimate", {}) or {}
            unrep_index.append({
                "name": name,
                "link": fname.replace(".md", ""),
                "category": info.get("category", "—"),
                "described": fmt_int(estimate.get("described")),
                "range": f"{fmt_int(estimate.get('low'))} – {fmt_int(estimate.get('high'))}",
            })

        (target / fname).write_text(md)

    # Home
    (target / "Home.md").write_text(render_home(rep_index, unrep_index))

    print(f"wrote {n_rep} represented + {n_unrep} unrepresented + 1 Home page to {target}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
