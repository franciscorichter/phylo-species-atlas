#!/usr/bin/env python3
"""Scaffold the per-partition and per-category expert agent folders.

Internal documentation system. Not coupled to the parent paper or the live
website. Creates one folder per partition / category under
``Data/experts/`` with seven artifacts:

  AGENT.md         research-protocol instructions (from _shared template)
  report.html      scholarly wiki (skeleton from _shared template; agent fills)
  timeline.csv     year-by-year list of every notable phylogeny (schema in _shared)
  subclades.yaml   internal taxonomic structure of the partition
  discoveries.yaml external resources / preprints / in-prep not yet in atlas
  gaps.md          narrative on what's still dark at sub-partition resolution
  references.bib   BibTeX bibliography accumulated by the agent
  state.yaml       maintenance log

The HTML report embeds the timeline plot via Vega-Lite (loaded from CDN at
view time) and the reference list rendered from BibTeX at scaffold time.

Idempotency: report.html is overwritten only when it still carries the
template's stub placeholders (``<em>Empty section`` text); once the agent
has filled in narrative the file is preserved (unless ``--force``). Refresh
flow: agent edits report.html, timeline.csv, references.bib in place; re-run
this script to refresh the embedded timeline plot data and the rendered
reference list while preserving the agent's prose.

Usage:
    cd Data
    python3 scripts/build_experts.py [--force]
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import html
import json
import pathlib
import re
import sys
from typing import Any

try:
    import yaml
except ImportError:
    sys.stderr.write("pyyaml is required. Install with: pip install pyyaml\n")
    sys.exit(1)

# ─── Paths ─────────────────────────────────────────────────────────────────
DATA = pathlib.Path(__file__).resolve().parent.parent  # Data/
PARTITION_INFOS = DATA / "site" / "data" / "partitions"

EXPERTS = DATA / "experts"
SHARED = EXPERTS / "_shared"
CATEGORY_DIR = EXPERTS / "_category"
PARTITION_DIR = EXPERTS / "partitions"

AGENT_TEMPLATE = SHARED / "AGENT_TEMPLATE.md"
REPORT_TEMPLATE = SHARED / "report_template.html"

CATEGORY_SLUG_MAP = {
    "Vertebrates": "vertebrates",
    "Plants": "plants",
    "Arthropods": "arthropods",
    "Other animals": "other_animals",
    "Microbes & protists": "microbes_protists",
    "Microbes and protists": "microbes_protists",
    "Not yet represented": "not_yet_represented",
}


# ─── Helpers ───────────────────────────────────────────────────────────────

def load_partition_info(slug: str) -> dict[str, Any]:
    info_path = PARTITION_INFOS / slug / "info.yaml"
    if not info_path.exists():
        return {}
    with info_path.open() as f:
        return yaml.safe_load(f) or {}


def discover_partition_slugs() -> list[str]:
    if not PARTITION_INFOS.exists():
        return []
    return sorted(
        d.name for d in PARTITION_INFOS.iterdir()
        if d.is_dir() and (d / "info.yaml").exists()
    )


def render_template(template_text: str, fields: dict[str, str]) -> str:
    def sub(m: re.Match[str]) -> str:
        return fields.get(m.group(1).strip(), m.group(0))
    return re.sub(r"\{\{([A-Z_]+)\}\}", sub, template_text)


def write_if_safe(path: pathlib.Path, content: str, force: bool, *,
                  hand_edit_markers: list[str] | None = None,
                  rendered_fields: dict[str, str] | None = None) -> str:
    """Write content, preserving hand-edited files.

    When ``hand_edit_markers`` is provided and NONE of the markers are present
    in the existing file (i.e. the agent has filled in the narrative), the
    file is preserved EXCEPT for ``{{...}}`` placeholders, which are refreshed
    in place from ``rendered_fields``. This lets the scaffold refresh embedded
    reference lists and timeline JSON without disturbing hand-written prose.
    """
    if not path.exists():
        path.write_text(content)
        return "created"
    if force:
        path.write_text(content)
        return "overwritten"
    existing = path.read_text()
    if not existing.strip():
        path.write_text(content)
        return "overwritten"

    # Has hand-edit markers semantics? (narrative file like report.html)
    if hand_edit_markers is not None:
        if any(m in existing for m in hand_edit_markers):
            # Still skeleton-shaped → overwrite with fresh template
            path.write_text(content)
            return "overwritten"
        # Narrative has been written. Refresh only the {{...}} placeholders.
        if "{{" in existing and rendered_fields:
            refreshed = render_template(existing, rendered_fields)
            if refreshed != existing:
                path.write_text(refreshed)
                return "refreshed"
        return "preserved"

    # No hand-edit semantics (source files): preserve once written, unless
    # the file is still in template state (contains {{}} or matches initial
    # comment-only stub). We use {{}} as the indicator of template state.
    if "{{" in existing:
        path.write_text(content)
        return "overwritten"
    return "preserved"


# ─── BibTeX → numbered reference list ──────────────────────────────────────

_BIB_ENTRY_RE = re.compile(r"@(\w+)\s*\{\s*([^,\s]+)\s*,(.*?)\n\}", re.DOTALL)
_BIB_FIELD_RE = re.compile(r"(\w+)\s*=\s*[{\"](.*?)[}\"]\s*,?\s*$", re.MULTILINE)


def parse_bib(bib_text: str) -> list[dict[str, str]]:
    entries: list[dict[str, str]] = []
    for m in _BIB_ENTRY_RE.finditer(bib_text):
        etype, ekey, ebody = m.group(1), m.group(2), m.group(3)
        fields = {"_type": etype, "_key": ekey}
        for fm in _BIB_FIELD_RE.finditer(ebody):
            fields[fm.group(1).lower()] = fm.group(2).strip()
        entries.append(fields)
    return entries


def format_reference_html(e: dict[str, str]) -> str:
    parts = []
    if e.get("author"):
        parts.append(html.escape(e["author"]))
    if e.get("year"):
        parts.append(f"({html.escape(e['year'])})")
    if e.get("title"):
        parts.append(html.escape(e["title"]).rstrip(".") + ".")
    if e.get("journal"):
        parts.append(f"<em>{html.escape(e['journal'])}</em>")
    bits = []
    if e.get("volume"):
        bits.append(html.escape(e["volume"]))
    if e.get("pages"):
        bits.append(html.escape(e["pages"]))
    if bits:
        parts.append(", ".join(bits) + ".")
    text = " ".join(parts)
    if e.get("doi"):
        text += f' <a href="https://doi.org/{html.escape(e["doi"])}">doi:{html.escape(e["doi"])}</a>'
    return f'<li id="ref-{html.escape(e.get("_key",""))}">{text}</li>'


def render_references(bib_path: pathlib.Path) -> str:
    if not bib_path.exists():
        return '    <li><em>No references yet.</em></li>'
    entries = parse_bib(bib_path.read_text())
    if not entries:
        return '    <li><em>No references yet — add BibTeX entries to <code>references.bib</code>.</em></li>'
    # Sort by year for the rendered list (ascending); agents typically cite chronologically.
    entries.sort(key=lambda e: (e.get("year", "9999"), e.get("_key", "")))
    return "\n".join("    " + format_reference_html(e) for e in entries)


# ─── Timeline.csv → embedded JSON ──────────────────────────────────────────

def _split_timeline_rows(timeline_path: pathlib.Path) -> tuple[list[dict], list[dict]]:
    """Split timeline rows into tree-bearing (n_tips > 0) and non-tree
    (reviews, dissertations, books, methods papers, databases). Non-tree
    rows should not appear in the scatter plot — they would be misleadingly
    plotted at y=1 — and instead are listed separately below the figure."""
    if not timeline_path.exists():
        return [], []
    with timeline_path.open() as f:
        rows = list(csv.DictReader(f))
    with_tree: list[dict] = []
    no_tree: list[dict] = []
    for r in rows:
        try:
            r["year"] = int(r["year"])
        except (KeyError, ValueError):
            continue
        try:
            n = int(r["n_tips"] or 0)
        except (KeyError, ValueError):
            n = 0
        r["n_tips"] = n
        if n > 0:
            with_tree.append(r)
        else:
            no_tree.append(r)
    return with_tree, no_tree


def render_timeline_json(timeline_path: pathlib.Path) -> str:
    with_tree, _ = _split_timeline_rows(timeline_path)
    return json.dumps(with_tree, ensure_ascii=False)


def render_timeline_non_tree_html(timeline_path: pathlib.Path) -> str:
    """Render the non-tree studies (reviews, books, dissertations, methods,
    catalogues, databases) as a compact dated list below the scatter plot."""
    _, no_tree = _split_timeline_rows(timeline_path)
    if not no_tree:
        return '    <li class="muted"><em>No non-tree studies recorded.</em></li>'
    no_tree.sort(key=lambda r: (r.get("year", 0), r.get("citekey", "")))
    items = []
    for r in no_tree:
        role = r.get("atlas_role", "").strip()
        study = html.escape(r.get("study", ""))
        journal = html.escape(r.get("journal", ""))
        year = r.get("year", "")
        role_suffix = f" — {html.escape(role)}" if role else ""
        items.append(
            f'    <li><strong>{year}</strong> · {study} '
            f'<span class="muted">({journal})</span>{role_suffix}</li>'
        )
    return "\n".join(items)


# ─── Per-partition scaffold ────────────────────────────────────────────────

DEFAULT_SUBCLADES_YAML = """# Internal taxonomic structure of the partition.
# Schema (per node):
#   - name: <clade name>
#     rank: <order | family | etc.>
#     described: <int or null>
#     phylo_status: <comprehensive | partial | sub-clade-only | dark>
#     atlas_resource: <slug-of-tree | null>
#     children:
#       - ... (recursive)
#
# Populated by the expert agent during research.
clade: {NAME}
rank: partition
children: []
"""

DEFAULT_DISCOVERIES_YAML = """# External phylogenetic resources and frontier work the expert has surfaced
# but the atlas does not yet incorporate. Internal-only; surfacing only.
#
# Per entry:
#   - kind: <preprint | deposited_tree | sequencing_initiative | review | other>
#     citekey: <bib key, optional>
#     title: <string>
#     year: <int>
#     url: <doi or URL>
#     why_relevant: <one-liner>
#     status: <noted | candidate_for_atlas | rejected>
#     notes: <free text>
items: []
"""

DEFAULT_GAPS_MD = """# Sub-partition gaps — {NAME}

Internal narrative on what remains phylogenetically dark **inside** the
{NAME} partition, even where the atlas already ships a partition-canonical
tree. Populated by the expert agent during research.

The atlas-level audit reports this partition's overall status; this file
records the finer-grained picture (which orders/families/lineages within the
partition are well-resolved, partial, or unsampled).
"""


def write_partition_scaffold(slug: str, info: dict[str, Any], force: bool) -> dict[str, str]:
    out_dir = PARTITION_DIR / slug
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(exist_ok=True)

    name = info.get("partition", slug.replace("_", " ").title())
    category = info.get("category", "Unknown")

    fields = {
        "NAME": name,
        "SLUG": slug,
        "TYPE": "partition",
        "CATEGORY": category,
        "SNAPSHOT_DATE": dt.date.today().isoformat(),
    }

    timeline_csv = out_dir / "timeline.csv"
    bib_path = out_dir / "references.bib"
    fields["TIMELINE_JSON"] = render_timeline_json(timeline_csv)
    fields["TIMELINE_NON_TREE_HTML"] = render_timeline_non_tree_html(timeline_csv)
    fields["REFERENCES_HTML"] = render_references(bib_path)

    actions: dict[str, str] = {}

    actions["AGENT.md"] = write_if_safe(
        out_dir / "AGENT.md", render_template(AGENT_TEMPLATE.read_text(), fields), force,
    )
    actions["report.html"] = write_if_safe(
        out_dir / "report.html", render_template(REPORT_TEMPLATE.read_text(), fields), force,
        hand_edit_markers=["<em>Empty section"],
        rendered_fields=fields,
    )
    actions["references.bib"] = write_if_safe(
        bib_path,
        f"% Bibliography for the {name} expert wiki.\n"
        f"% Internal documentation under the Phylogenetic Atlas project.\n"
        f"% Add BibTeX entries with verified DOIs as research proceeds.\n",
        force,
    )
    actions["timeline.csv"] = write_if_safe(
        timeline_csv,
        "year,citekey,study,journal,doi,n_tips,species_in_clade_at_time,dated,"
        "crown_ma,methods,calibration,imputation,in_atlas,atlas_role,notes\n",
        force,
    )
    actions["subclades.yaml"] = write_if_safe(
        out_dir / "subclades.yaml",
        DEFAULT_SUBCLADES_YAML.format(NAME=name),
        force,
    )
    actions["discoveries.yaml"] = write_if_safe(
        out_dir / "discoveries.yaml",
        DEFAULT_DISCOVERIES_YAML,
        force,
    )
    actions["gaps.md"] = write_if_safe(
        out_dir / "gaps.md",
        DEFAULT_GAPS_MD.format(NAME=name),
        force,
    )
    state_yaml = (
        f"slug: {slug}\n"
        f"name: {name}\n"
        f"type: partition\n"
        f"category: {category}\n"
        f"last_refresh: null\n"
        f"refresh_cadence_days: 90\n"
        f"references_count: 0\n"
        f"timeline_rows: 0\n"
        f"report_status: skeleton\n"
        f"sections_completed: []\n"
        f"pending_atlas_updates: []\n"
        f"open_questions: []\n"
    )
    actions["state.yaml"] = write_if_safe(out_dir / "state.yaml", state_yaml, force)
    return actions


def write_category_scaffold(category: str, member_slugs: list[str], force: bool) -> dict[str, str]:
    cat_slug = CATEGORY_SLUG_MAP.get(category, category.lower().replace(" ", "_").replace("&", "and"))
    out_dir = CATEGORY_DIR / cat_slug
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "figures").mkdir(exist_ok=True)

    bib_path = out_dir / "references.bib"

    fields = {
        "NAME": category,
        "SLUG": cat_slug,
        "TYPE": "category",
        "CATEGORY": category,
        "SNAPSHOT_DATE": dt.date.today().isoformat(),
        "TIMELINE_JSON": "[]",
        "TIMELINE_NON_TREE_HTML": '    <li class="muted"><em>No non-tree studies recorded.</em></li>',
        "REFERENCES_HTML": render_references(bib_path),
    }

    actions: dict[str, str] = {}
    actions["AGENT.md"] = write_if_safe(
        out_dir / "AGENT.md", render_template(AGENT_TEMPLATE.read_text(), fields), force,
    )
    actions["report.html"] = write_if_safe(
        out_dir / "report.html", render_template(REPORT_TEMPLATE.read_text(), fields), force,
        hand_edit_markers=["<em>Empty section"],
    )
    actions["references.bib"] = write_if_safe(
        bib_path,
        f"% Bibliography for the {category} category-coordinator wiki.\n"
        f"% Internal documentation.\n",
        force,
    )
    state_yaml = (
        f"slug: {cat_slug}\n"
        f"name: {category}\n"
        f"type: category\n"
        f"member_partitions: {member_slugs}\n"
        f"last_refresh: null\n"
        f"refresh_cadence_days: 180\n"
        f"report_status: skeleton\n"
    )
    actions["state.yaml"] = write_if_safe(out_dir / "state.yaml", state_yaml, force)
    return actions


# ─── Entry point ───────────────────────────────────────────────────────────

def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="Overwrite hand-edited files. Use with care.")
    args = ap.parse_args()

    slugs = discover_partition_slugs()
    if not slugs:
        sys.stderr.write(f"No partitions found under {PARTITION_INFOS}\n")
        return 1

    summary: dict[str, dict[str, str]] = {}
    cat_members: dict[str, list[str]] = {}

    for slug in slugs:
        info = load_partition_info(slug)
        summary[slug] = write_partition_scaffold(slug, info, args.force)
        cat_members.setdefault(info.get("category", "Unknown"), []).append(slug)

    for cat, members in cat_members.items():
        if cat == "Not yet represented":
            continue
        cat_slug_path = f"_category/{CATEGORY_SLUG_MAP.get(cat, cat)}"
        summary[cat_slug_path] = write_category_scaffold(cat, sorted(members), args.force)

    n_partitions = len(slugs)
    n_categories = len([c for c in cat_members if c != "Not yet represented"])
    print(f"Scaffolded {n_partitions} partition experts + {n_categories} category coordinators.")
    by_action: dict[str, int] = {}
    for files in summary.values():
        for a in files.values():
            by_action[a] = by_action.get(a, 0) + 1
    print(f"Files: {by_action}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
