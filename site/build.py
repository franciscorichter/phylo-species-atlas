"""Build site/data.json from standardized/ CSVs and data_provenance.csv.

Reads from the repo root (one level up) and writes a single bundle the static
page can fetch on load. Tree files themselves are served directly from
standardized/trees/ — we don't duplicate them here.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STD = ROOT / "standardized"
PROVENANCE = ROOT / "data_provenance.csv"
OUT = Path(__file__).resolve().parent / "data.json"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def parse_bool(v: str) -> bool:
    return str(v).strip().upper() == "TRUE"


def parse_int(v) -> int | None:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def parse_float(v) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


TIP_RE = re.compile(r"[(,]\s*(\d+|[A-Za-z_][\w.]*)")


def count_tips(newick_path: Path) -> int:
    """Cheap tip count: number of leaf labels in the Newick string.

    A leaf is a token that appears immediately after '(' or ',' and is not
    itself followed by '('. We approximate by counting all `[(,]label` matches
    where label is a numeric ID or alphanumeric token — this is exact for the
    standardized trees (numeric tip labels)."""
    text = newick_path.read_text(encoding="utf-8", errors="ignore")
    return len(TIP_RE.findall(text))


def main() -> None:
    metadata = read_csv(STD / "metadata.csv")
    provenance = read_csv(PROVENANCE)
    prov_by_group = {row["group"]: row for row in provenance}

    trees = []
    by_group_aggregate = defaultdict(lambda: {"trees": 0, "tips": 0})
    dated_count = 0

    for row in metadata:
        filename = row["filename"]
        group = row["group"]
        study = row["study"]
        ntips_meta = parse_int(row["ntips"])
        dated = parse_bool(row["dated"])
        dated_count += int(dated)

        tree_path = STD / "trees" / filename
        size_bytes = tree_path.stat().st_size if tree_path.exists() else None

        prov = prov_by_group.get(group, {})
        described = parse_int(prov.get("described_species"))
        coverage_pct = parse_float(prov.get("coverage_pct"))

        trees.append({
            "filename": filename,
            "group": group,
            "study": study,
            "ntips": ntips_meta,
            "dated": dated,
            "size_bytes": size_bytes,
            "year": parse_int(prov.get("year")),
            "journal": prov.get("journal") or None,
            "doi": prov.get("doi") or None,
            "crown_ma": parse_float(prov.get("crown_ma")),
            "described_species": described,
            "coverage_pct": coverage_pct,
            "data_source": prov.get("data_source") or None,
            "download_url": prov.get("download_url") or None,
            "methods_brief": prov.get("methods_brief") or None,
            "notes": prov.get("notes") or None,
        })

        agg = by_group_aggregate[group]
        agg["trees"] += 1
        agg["tips"] += ntips_meta or 0

    # Sort by tips desc — better default for tree picker
    trees.sort(key=lambda r: (-(r["ntips"] or 0), r["filename"]))

    # Coverage by group, for the bar chart. Use the provenance row when present.
    coverage_rows = []
    for prov_row in provenance:
        described = parse_int(prov_row.get("described_species"))
        cov = parse_float(prov_row.get("coverage_pct"))
        tips = parse_int(prov_row.get("tips"))
        if described and cov is not None:
            coverage_rows.append({
                "group": prov_row["group"],
                "study": prov_row.get("study") or "",
                "year": parse_int(prov_row.get("year")),
                "tips": tips,
                "described_species": described,
                "coverage_pct": cov,
                "dated": parse_bool(prov_row.get("dated", "FALSE")),
            })
    coverage_rows.sort(key=lambda r: -r["coverage_pct"])

    summary = {
        "n_trees": len(trees),
        "n_dated": dated_count,
        "n_undated": len(trees) - dated_count,
        "n_groups": len({r["group"] for r in trees}),
        "n_datasets": len(provenance),
        "total_tips": sum((r["ntips"] or 0) for r in trees),
    }

    bundle = {
        "schema_version": 1,
        "generated_from": "standardized/metadata.csv + data_provenance.csv",
        "summary": summary,
        "trees": trees,
        "coverage": coverage_rows,
    }

    OUT.write_text(json.dumps(bundle, indent=2, ensure_ascii=False))
    size_kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT.name}: {len(trees)} trees, {len(coverage_rows)} coverage rows, {size_kb:.1f} KB", file=sys.stderr)


if __name__ == "__main__":
    main()
