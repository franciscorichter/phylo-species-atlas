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
ESTIMATES = ROOT / "data_estimates.csv"
OUT = Path(__file__).resolve().parent / "data.json"

# Map provenance group (dataset-level slug) → partition group in data_estimates.csv
# (the user-facing Figure 1 category). Sub-clade datasets inherit their parent's
# estimate. Cross-clade synthesis trees (eukaryotes, timetree, condamine) have no
# single partition home and are left unmapped.
PROVENANCE_TO_PARTITION = {
    "mammals": "Mammals", "birds_jetz": "Birds", "birds_mctavish": "Birds",
    "amphibians": "Amphibians", "frogs": "Amphibians", "salamanders": "Amphibians",
    "fish": "Fish", "neotropical_fish": "Fish",
    "squamates": "Squamates", "sharks": "Sharks", "turtles": "Turtles",
    "primates": "Mammals", "parrots": "Birds", "carnivora": "Mammals", "cetaceans": "Mammals",
    "seed_plants": "Seed plants", "ferns": "Ferns", "bryophytes": "Bryophytes",
    "grasses": "Seed plants", "orchids": "Seed plants", "asteraceae": "Seed plants",
    "solanaceae": "Seed plants", "cacti": "Seed plants", "eucalyptus": "Seed plants",
    "conifers": "Conifers",
    "insects": "Insects", "insects_chesters": "Insects", "butterflies": "Insects",
    "bees": "Insects", "hemiptera": "Insects", "termites": "Insects", "ants": "Insects",
    "spiders": "Spiders", "crustaceans": "Crustaceans",
    "corals": "Cnidaria", "cephalopods": "Cephalopods", "gastropods": "Gastropods",
    "sponges": "Sponges", "echinoderms": "Echinoderms", "nematodes": "Nematodes",
    "diatoms": "Diatoms", "bacteria": "Bacteria", "archaea": "Archaea", "fungi": "Fungi",
    "bivalves": "Bivalves", "bryozoa": "Bryozoa",
}


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


def _est_pct(tips, denominator):
    if tips is None or denominator is None or denominator <= 0:
        return None
    return min(100.0, round(100.0 * tips / denominator, 1))


def main() -> None:
    metadata = read_csv(STD / "metadata.csv")
    provenance = read_csv(PROVENANCE)
    estimates_rows = read_csv(ESTIMATES) if ESTIMATES.exists() else []
    prov_by_group = {row["group"]: row for row in provenance}
    estimates_by_partition = {row["group"]: row for row in estimates_rows}

    def estimate_for(prov_group: str | None):
        partition = PROVENANCE_TO_PARTITION.get(prov_group or "")
        if not partition:
            return None, None
        est = estimates_by_partition.get(partition)
        if not est:
            return partition, None
        return partition, {
            "estimated_total": parse_int(est.get("estimated_total")),
            "estimated_low": parse_int(est.get("estimated_low")),
            "estimated_high": parse_int(est.get("estimated_high")),
            "estimate_source": est.get("estimate_source") or None,
            "estimate_confidence": est.get("confidence") or None,
            "category": est.get("category") or None,
        }

    # Group provenance rows by their partition so the UI can offer dataset switches.
    datasets_by_partition: dict[str, list[dict]] = defaultdict(list)
    for prov_row in provenance:
        partition = PROVENANCE_TO_PARTITION.get(prov_row.get("group") or "")
        if not partition:
            continue
        datasets_by_partition[partition].append({
            "group": prov_row["group"],
            "tree_name": prov_row.get("tree_name") or None,
            "study": prov_row.get("study") or "",
            "year": parse_int(prov_row.get("year")),
            "tips": parse_int(prov_row.get("tips")),
            "described_species": parse_int(prov_row.get("described_species")),
            "dated": parse_bool(prov_row.get("dated", "FALSE")),
        })

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
        partition, est = estimate_for(group)

        is_anchor = bool(est and est.get("estimated_total") and described
                         and described >= 0.5 * est["estimated_total"])
        trees.append({
            "filename": filename,
            "group": group,
            "partition_group": partition,
            "is_partition_anchor": is_anchor,
            "study": study,
            "ntips": ntips_meta,
            "dated": dated,
            "size_bytes": size_bytes,
            "year": parse_int(prov.get("year")),
            "journal": prov.get("journal") or None,
            "doi": prov.get("doi") or None,
            "crown_ma": parse_float(prov.get("crown_ma")),
            "described_species": described,
            "described_source": prov.get("species_count_source") or None,
            "coverage_pct": coverage_pct,
            "data_source": prov.get("data_source") or None,
            "download_url": prov.get("download_url") or None,
            "methods_brief": prov.get("methods_brief") or None,
            "notes": prov.get("notes") or None,
            "estimate": est,
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
        partition, est = estimate_for(prov_row.get("group"))
        if described and cov is not None:
            row_out = {
                "group": prov_row["group"],
                "partition_group": partition,
                "tree_name": prov_row.get("tree_name") or None,
                "study": prov_row.get("study") or "",
                "year": parse_int(prov_row.get("year")),
                "tips": tips,
                "described_species": described,
                "described_source": prov_row.get("species_count_source") or None,
                "coverage_pct": cov,
                "dated": parse_bool(prov_row.get("dated", "FALSE")),
            }
            if est and est.get("estimated_total"):
                row_out["estimated_total"] = est["estimated_total"]
                row_out["estimated_low"] = est["estimated_low"]
                row_out["estimated_high"] = est["estimated_high"]
                row_out["estimate_source"] = est["estimate_source"]
                row_out["estimate_confidence"] = est["estimate_confidence"]
                row_out["coverage_pct_estimated"] = _est_pct(tips, est["estimated_total"])
                row_out["coverage_pct_estimated_low"] = _est_pct(tips, est["estimated_high"])  # high denom = low cov
                row_out["coverage_pct_estimated_high"] = _est_pct(tips, est["estimated_low"])  # low denom = high cov
                # Anchor = dataset's described count covers most of the partition's described.
                # Sub-clade trees (parrots within Birds, primates within Mammals) are not anchors.
                row_out["is_partition_anchor"] = bool(described and described >= 0.5 * est["estimated_total"])
            coverage_rows.append(row_out)
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
        "schema_version": 2,
        "generated_from": "standardized/metadata.csv + data_provenance.csv + data_estimates.csv",
        "summary": summary,
        "trees": trees,
        "coverage": coverage_rows,
        "datasets_by_partition": dict(datasets_by_partition),
    }

    OUT.write_text(json.dumps(bundle, indent=2, ensure_ascii=False))
    size_kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT.name}: {len(trees)} trees, {len(coverage_rows)} coverage rows, {size_kb:.1f} KB", file=sys.stderr)


if __name__ == "__main__":
    main()
