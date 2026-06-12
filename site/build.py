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
SOURCES = ROOT / "data_sources.csv"
DICTIONARY = STD / "dictionary.csv"
NAMES_DIR = STD / "names"
SITE_DATA = Path(__file__).resolve().parent / "data" / "partitions"  # web-data overrides
OUT = Path(__file__).resolve().parent / "data.json"

try:
    import yaml  # type: ignore
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

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

# Captures numeric tip IDs from a Newick string. Tips appear after '(' or ',',
# end at ':' or ',' or ')'. This is exact for the standardized trees (all tips
# are integer IDs from standardized/dictionary.csv).
TIP_ID_RE = re.compile(r"([(,])(-?\d+)(?=[:,)])")


def _json_default(o):
    """JSON encoder hook for types PyYAML produces (date, datetime)."""
    import datetime
    if isinstance(o, (datetime.date, datetime.datetime)):
        return o.isoformat()
    raise TypeError(f"Object of type {o.__class__.__name__} is not JSON serializable")


_AUDIT_FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?\n)---\s*$", re.DOTALL | re.MULTILINE)


def _load_audit_frontmatter(audits_dir: Path) -> dict[str, dict]:
    """Parse YAML frontmatter from every audits/partitions/<slug>.md file.

    Keyed by the `partition:` field. These records carry the audit's working
    metadata (shipped/candidate sub-clades, full estimate-source verification
    block, etc.) that the audit page renders. info.yaml carries only the
    verified subset the data view needs.
    """
    if not HAS_YAML or not audits_dir.exists():
        return {}
    out: dict[str, dict] = {}
    for md_path in sorted(audits_dir.glob("*.md")):
        text = md_path.read_text(encoding="utf-8")
        m = _AUDIT_FRONTMATTER_RE.match(text)
        if not m:
            continue
        try:
            data = yaml.safe_load(m.group(1))
        except yaml.YAMLError as e:
            print(f"WARN: bad frontmatter in {md_path}: {e}", file=sys.stderr)
            continue
        if not isinstance(data, dict):
            continue
        partition = data.get("partition")
        if not partition:
            continue
        out[partition] = data
    return out


def load_audit_overrides() -> dict[str, dict]:
    """Read site/data/partitions/<slug>/info.yaml AND audits/partitions/<slug>.md
    frontmatter, merging them into a single dict keyed by partition name. The
    info.yaml provides web-facing facts; the .md frontmatter provides the
    audit's process metadata (candidate/shipped sub-clades, audit trail).
    """
    if not HAS_YAML:
        print("PyYAML not installed — skipping audit overrides", file=sys.stderr)
        return {}
    out: dict[str, dict] = {}

    audits_md = _load_audit_frontmatter(ROOT / "audits" / "partitions")

    if SITE_DATA.exists():
        for info_path in sorted(SITE_DATA.glob("*/info.yaml")):
            try:
                data = yaml.safe_load(info_path.read_text(encoding="utf-8"))
            except yaml.YAMLError as e:
                print(f"WARN: skipping {info_path}: {e}", file=sys.stderr)
                continue
            if not isinstance(data, dict):
                continue
            partition = data.get("partition")
            if not partition:
                continue
            sub_dir = info_path.parent / "sub-phylos"
            subs = []
            if sub_dir.exists():
                for sub_info in sorted(sub_dir.glob("*/info.yaml")):
                    try:
                        sub_data = yaml.safe_load(sub_info.read_text(encoding="utf-8"))
                        if isinstance(sub_data, dict):
                            subs.append(sub_data)
                    except yaml.YAMLError:
                        continue
            data["_sub_phylos_loaded"] = subs
            # Merge audit .md frontmatter — adds shipped_subclades,
            # candidate_subclades, full estimate_source/canonical_tree blocks,
            # api_calls trail. info.yaml fields take precedence on conflicts.
            md = audits_md.pop(partition, None)
            if md:
                for k, v in md.items():
                    data.setdefault(k, v)
            out[partition] = data

    # Audits that haven't been migrated to site/data/ yet — still surface.
    for partition, md in audits_md.items():
        out[partition] = md
    return out


def _find_override_path(partition: str | None, group: str, filename: str) -> Path | None:
    """Locate an atlas-derived tree override under site/data/.

    Selection rule: each `group` (from standardized/metadata.csv) maps to
    exactly one possible override location.
      - Partition canonical (group == partition_slug):
        site/data/partitions/<partition_slug>/tree.nwk
      - Sub-clade (group != partition_slug):
        site/data/partitions/<partition_slug>/sub-phylos/<group>/tree.nwk

    The candidate is considered an override only when its file size differs
    from standardized/trees/<filename>.nwk (same size ⇒ just a copy, no
    derivation intent).
    """
    if not partition:
        return None
    pslug = partition.lower().replace(" ", "_")
    standard = STD / "trees" / filename
    if not standard.exists():
        return None
    if group == pslug:
        candidate = SITE_DATA / pslug / "tree.nwk"
    else:
        candidate = SITE_DATA / pslug / "sub-phylos" / group / "tree.nwk"
    if candidate.exists() and candidate.stat().st_size != standard.stat().st_size:
        return candidate
    return None


def _override_url(partition: str | None, group: str, filename: str) -> str | None:
    """Same selection as _find_override_path, but returns the URL path the
    viewer fetches (relative to the site root). The viewer combines this
    with its base URL (../  in dev, ./ on Pages) — both resolve correctly
    so long as we emit `data/partitions/...`."""
    path = _find_override_path(partition, group, filename)
    if not path:
        return None
    # SITE_DATA.parent.parent == site/, so this yields "data/partitions/..."
    rel = path.relative_to(SITE_DATA.parent.parent)
    return str(rel).replace("\\", "/")


def _override_taxonomy(partition: str | None, group: str, filename: str,
                       full_names: dict[str, str]) -> dict | None:
    if not partition or not full_names:
        return None
    candidate = _find_override_path(partition, group, filename)
    if not candidate:
        return None
    return tip_taxonomy_for_override(candidate, full_names)


def _override_tree_meta(partition: str | None, group: str, filename: str) -> dict:
    """When an override is active, read info.yaml's tree block to capture
    the curator-asserted ntips/dated/source — those reflect the override
    file's reality, not the upstream CSV row (which still describes the
    superseded shipped tree)."""
    if not HAS_YAML or not partition:
        return {}
    path = _find_override_path(partition, group, filename)
    if not path:
        return {}
    info_yaml = path.parent / "info.yaml"
    if not info_yaml.exists():
        return {}
    try:
        data = yaml.safe_load(info_yaml.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError:
        return {}
    tree = data.get("tree") or {}
    if not isinstance(tree, dict):
        return {}
    out = {}
    for k in ("ntips", "dated"):
        if k in tree and tree[k] is not None:
            out[k] = tree[k]
    src = tree.get("source") or {}
    if src.get("key"): out["override_source_key"] = src["key"]
    if src.get("doi"): out["override_source_doi"] = src["doi"]
    if src.get("study"): out["override_source_study"] = src["study"]
    return out


def tip_taxonomy_for_override(override_path: Path, names: dict[str, str]) -> dict | None:
    """Like tip_taxonomy_from_names but for an atlas-derived override tree:
    parse the override Newick directly, extract its tip IDs, and build the
    taxonomy stats from THAT subset (not from the full shard, which still
    reflects the original multi-individual tree)."""
    if not override_path.exists():
        return None
    text = override_path.read_text(encoding="utf-8", errors="ignore")
    ids = [m.group(2) for m in TIP_ID_RE.finditer(text)]
    if not ids:
        return None
    subset = {i: names[i] for i in set(ids) if i in names}
    tx = tip_taxonomy_from_names(subset)
    tx["total_tips"] = len(ids)  # actual tip positions in the override tree
    return tx


def tip_taxonomy_from_names(names: dict[str, str]) -> dict:
    """Break a tree's tip names down by taxonomic level.

    Tries to parse `Genus_species[_subspecies_or_code]` (Latin binomial)
    naming, where Genus is capitalized + lowercase-tail and species is
    all-lowercase. This applies cleanly to vertebrates, most plants, fungi,
    and many invertebrate trees.

    BUT many partitions use non-binomial labels:
      - bacteria/archaea use GTDB assembly IDs (`GB_GCA_024640485.1`)
      - bryophytes use voucher-prefixed names (`UFG_393202_P031.Buxbaumia_aphylla`)
      - crustaceans use single-word names (`Sympagurus`)
      - diatoms use concatenated genus+species (`Pseudo-nitzschiaheimii`)

    For non-binomial trees, treating each unique ID as a distinct species
    is the safe interpretation (each tip is a distinct sequence/specimen).
    `binomial_format` records which case we hit.
    """
    species: set[str] = set()
    subspecies_tips = 0
    binomial_count = 0
    n = 0
    for name in names.values():
        if not name:
            continue
        n += 1
        parts = name.split("_")
        if len(parts) >= 2:
            first, second = parts[0], parts[1]
            is_binomial = (
                first[:1].isupper() and first[1:].islower() and len(first) >= 2
                and second[:1].islower()
                and second.replace(".", "").replace("-", "").isalpha()
            )
            if is_binomial:
                binomial_count += 1
                species.add(f"{first}_{second}")
                if len(parts) > 2:
                    third = parts[2]
                    if third and third[0].islower() and third.isalpha():
                        subspecies_tips += 1
                continue
        # Not binomial — treat the full name as a unique species label.
        species.add(name)
    binomial_format = n > 0 and binomial_count > n * 0.5
    if binomial_format:
        unique_species = len(species)
    else:
        # Don't trust the parsing — fall back to unique-IDs.
        unique_species = len(names)
        subspecies_tips = 0
    return {
        "unique_ids": len(names),
        "unique_species": unique_species,
        "subspecies_tips": subspecies_tips,
        "redundancy": len(names) - unique_species,
        "binomial_format": binomial_format,
    }


def build_name_shards(metadata: list[dict]) -> tuple[int, dict[str, dict[str, str]]]:
    """Generate one `<filename>.json` per tree with {id: name} for every tip
    that tree uses. Stored in standardized/names/, loaded on demand by the
    web viewer to populate hover tooltips with species/specimen names.

    Also returns the in-memory {filename: {id: name}} mapping so callers
    (compute_tip_taxonomy) can derive per-tree taxonomy stats without
    re-reading from disk — this matters on CI runs where the shards don't
    exist yet when the main build loop runs.

    Returns (total_bytes_written, names_by_filename).
    """
    names_by_filename: dict[str, dict[str, str]] = {}
    if not DICTIONARY.exists():
        print("dictionary.csv missing — skipping name shards", file=sys.stderr)
        return 0, names_by_filename
    NAMES_DIR.mkdir(exist_ok=True)
    print("loading dictionary…", file=sys.stderr)
    name_by_id: dict[str, str] = {}
    with DICTIONARY.open(newline="", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            name_by_id[row["id"]] = row["standardized_name"]
    total = 0
    for row in metadata:
        filename = row["filename"]
        tree_path = STD / "trees" / filename
        if not tree_path.exists():
            continue
        text = tree_path.read_text(encoding="utf-8", errors="ignore")
        ids = {m.group(2) for m in TIP_ID_RE.finditer(text)}
        shard = {i: name_by_id[i] for i in ids if i in name_by_id}
        names_by_filename[filename] = shard
        out_path = NAMES_DIR / f"{filename}.json"
        out_path.write_text(json.dumps(shard, separators=(",", ":"), ensure_ascii=False))
        total += out_path.stat().st_size
    return total, names_by_filename


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


def find_provenance_match(tree_group: str, tree_ntips, provenance: list[dict]):
    """Resolve a standardized tree to its data_provenance.csv row.

    Standardized groups are simplified (`birds`) while provenance can have
    dataset-level rows (`birds_jetz`, `birds_mctavish`). We try exact match
    first, then prefix-with-underscore variants, preferring the row whose
    `tips` matches the standardized tree's tip count.
    """
    for p in provenance:
        if p.get("group") == tree_group:
            return p
    candidates = [p for p in provenance if (p.get("group") or "").startswith(tree_group + "_")]
    if not candidates:
        return None
    if tree_ntips is not None:
        for p in candidates:
            try:
                if int(p.get("tips") or "") == tree_ntips:
                    return p
            except (TypeError, ValueError):
                pass
    return candidates[0]


def main() -> None:
    metadata = read_csv(STD / "metadata.csv")
    provenance = read_csv(PROVENANCE)
    estimates_rows = read_csv(ESTIMATES) if ESTIMATES.exists() else []
    estimates_by_partition = {row["group"]: row for row in estimates_rows}
    sources_rows = read_csv(SOURCES) if SOURCES.exists() else []
    sources_by_key = {row["source_key"]: row for row in sources_rows}

    def source_info(key: str | None):
        if not key:
            return None
        row = sources_by_key.get(key)
        if not row:
            return {"label": key, "type": None, "doi": None, "url": None, "citation": None}
        return {
            "label": key,
            "type": row.get("type") or None,
            "doi": row.get("doi") or None,
            "url": row.get("url") or None,
            "citation": row.get("citation") or None,
        }

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

    # Build name shards FIRST so the in-memory names map is available
    # when the main loop attaches tip_taxonomy. Otherwise on a fresh CI
    # checkout (no standardized/names/ on disk yet) every tree gets
    # tip_taxonomy = None and the web app misses species-level coverage.
    shard_bytes, names_by_filename = build_name_shards(metadata)

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

        prov = find_provenance_match(group, ntips_meta, provenance) or {}
        prov_group = prov.get("group") or group
        described = parse_int(prov.get("described_species"))
        coverage_pct = parse_float(prov.get("coverage_pct"))
        partition, est = estimate_for(prov_group)

        # Special case: Condamine 2019 Crocodylia is the family-set member
        # that maps to its own partition (Crocodilians) — promote it BEFORE
        # the described-species fallback so it inherits the partition estimate.
        if group == "condamine_Crocodylia":
            partition = "Crocodilians"
            est_row = estimates_by_partition.get("Crocodilians") or {}
            est = {
                "category": "Vertebrates",
                "estimated_total": parse_int(est_row.get("total")),
                "estimated_low": parse_int(est_row.get("low")),
                "estimated_high": parse_int(est_row.get("high")),
                "estimate_source": est_row.get("estimate_source"),
                "estimate_confidence": est_row.get("confidence"),
            }

        # Fallback: when data_provenance.csv has no described_species for
        # this tree (e.g. Bryozoa, where the original CSV row was empty),
        # use the partition-level described count from data_estimates.csv
        # so the tree still appears in the coverage chart.
        if not described and partition:
            est_row = estimates_by_partition.get(partition) or {}
            described = parse_int(est_row.get("described"))
        if described and coverage_pct is None and ntips_meta:
            coverage_pct = round(100 * ntips_meta / described, 1)

        # Condamine 2019 family trees are 218 small family-level trees — group
        # them under one synthetic partition so they collapse cleanly in the UI.
        if not partition and group.startswith("condamine_"):
            partition = "Condamine 2019 families"
            # Family-level insect trees from Condamine 2019 → Arthropods.
            if not est:
                est = {"category": "Arthropods", "estimated_total": None,
                       "estimated_low": None, "estimated_high": None,
                       "estimate_source": None, "estimate_confidence": None}

        category = est.get("category") if est else None

        is_anchor = bool(est and est.get("estimated_total") and described
                         and described >= 0.5 * est["estimated_total"])
        trees.append({
            "filename": filename,
            "group": group,
            "provenance_group": prov_group,
            "partition_group": partition,
            "category": category,
            "is_partition_anchor": is_anchor,
            "study": study or prov.get("study") or "",
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
            "described_source_info": source_info(prov.get("species_count_source")),
            "estimate_source_info": source_info(est.get("estimate_source")) if est else None,
            "tip_taxonomy": (tip_taxonomy_from_names(names_by_filename[filename])
                             if filename in names_by_filename else None),
            # When site/data/partitions/<slug>/tree.nwk exists and differs from
            # standardized/trees/<filename>.nwk, it's an atlas-derived override
            # (e.g., turtles species-level pruning). The viewer prefers it; the
            # original stays under uncertainty/ as the multi-individual version.
            "tree_url_override": _override_url(partition, group, filename),
            "tip_taxonomy_override": _override_taxonomy(partition, group, filename,
                                                       names_by_filename.get(filename, {})),
        })
        # Apply override-side metadata when the curator has installed a
        # replacement tree (e.g., Nakov 2018 swap-in for diatoms). The
        # info.yaml tree block then carries the truth about ntips/dated/source
        # for the override file, since the upstream CSV row still describes
        # the superseded shipped tree.
        if trees[-1].get("tree_url_override"):
            meta = _override_tree_meta(partition, group, filename)
            if "ntips" in meta:
                trees[-1]["ntips"] = meta["ntips"]
            if "dated" in meta:
                trees[-1]["dated"] = bool(meta["dated"])
            if meta.get("override_source_key"):
                trees[-1]["study"] = meta.get("override_source_study") or meta["override_source_key"]
            if meta.get("override_source_doi"):
                trees[-1]["doi"] = meta["override_source_doi"]

        agg = by_group_aggregate[group]
        agg["trees"] += 1
        agg["tips"] += ntips_meta or 0

    # Reclassify the Conifers tree (Ran 2018 Pinaceae) as a sub-clade of
    # Seed plants. Conifers were in data_estimates.csv as a top-level
    # partition but Conifer species (~615) double-count within Seed plants
    # (~383k, since Smith Brown 2018 ALLMB includes gymnosperms). Curated
    # the tree's partition_group manually here so the chart's partition
    # filter sees Conifers under Seed plants → sub-phylos → conifers/.
    for t in trees:
        if t.get("group") == "conifers":
            t["partition_group"] = "Seed plants"

    # Identify the canonical tree per partition: the one with the most tips.
    # All other trees in that partition are sub-clades (Layer 2).
    partition_max_tips: dict[str, int] = {}
    for t in trees:
        p = t.get("partition_group")
        if not p:
            continue
        n = t.get("ntips") or 0
        if n > partition_max_tips.get(p, -1):
            partition_max_tips[p] = n
    for t in trees:
        p = t.get("partition_group")
        n = t.get("ntips") or 0
        # Condamine 218 family trees are all peer family-level — no one of them
        # is "the canonical Condamine tree", so leave them all as sub-clades.
        if p == "Condamine 2019 families":
            t["is_partition_canonical"] = False
        elif p == "Crocodilians" and t.get("group") == "condamine_Crocodylia":
            t["is_partition_canonical"] = True
        else:
            t["is_partition_canonical"] = bool(p and partition_max_tips.get(p) == n)

    # Sort by tips desc — better default for tree picker
    trees.sort(key=lambda r: (-(r["ntips"] or 0), r["filename"]))

    # Coverage by tree — one bar per shipped standardized tree that has
    # described_species and coverage_pct attached from its provenance row.
    coverage_rows = []
    for t in trees:
        # When the tree has an atlas-derived species-level override (e.g.,
        # turtles' 287-tip pruning), use the override's species count as the
        # numerator for both 'described' and 'estimated' coverage. The CSV
        # coverage_pct came from the original multi-individual tip count and
        # is misleading (turtles shows 100% tip-level but only 79.7% species-level).
        described = t["described_species"]
        if not described:
            continue
        tx_ov = t.get("tip_taxonomy_override") or {}
        tx = t.get("tip_taxonomy") or {}
        # When an override file is installed, trust the curator-asserted
        # ntips (which we synced from info.yaml above) — overrides may use
        # non-atlas tip labels that the dictionary-based taxonomy lookup
        # can't parse (e.g., Nakov 2018 diatom binomials not in atlas
        # dictionary). Otherwise fall back to per-species taxonomy stats.
        if t.get("tree_url_override"):
            effective_count = t["ntips"]
        else:
            effective_count = (tx_ov.get("unique_species")
                               or tx.get("unique_species")
                               or t["ntips"])
        if effective_count == t["ntips"] and t["coverage_pct"] is not None:
            # No taxonomic correction needed — keep CSV's pre-computed coverage.
            cov = t["coverage_pct"]
        else:
            # Recompute from species-level effective count.
            cov = round(100 * effective_count / described, 1)
        row_out = {
            "filename": t["filename"],
            "group": t["group"],
            "provenance_group": t["provenance_group"],
            "partition_group": t["partition_group"],
            "category": t.get("category"),
            "study": t["study"],
            "year": t["year"],
            "tips": effective_count,            # what users see as the species/tip count
            "raw_tips": t["ntips"],             # CSV-recorded raw tip positions
            "described_species": described,
            "described_source": t["described_source"],
            "coverage_pct": cov,
            "dated": t["dated"],
            "is_partition_anchor": t["is_partition_anchor"],
            "is_partition_canonical": t["is_partition_canonical"],
        }
        est = t["estimate"]
        if est and est.get("estimated_total"):
            row_out["estimated_total"] = est["estimated_total"]
            row_out["estimated_low"] = est["estimated_low"]
            row_out["estimated_high"] = est["estimated_high"]
            row_out["estimate_source"] = est["estimate_source"]
            row_out["estimate_confidence"] = est["estimate_confidence"]
            row_out["coverage_pct_estimated"] = _est_pct(effective_count, est["estimated_total"])
            row_out["coverage_pct_estimated_low"] = _est_pct(effective_count, est["estimated_high"])  # high denom = low cov
            row_out["coverage_pct_estimated_high"] = _est_pct(effective_count, est["estimated_low"])  # low denom = high cov
        coverage_rows.append(row_out)
    coverage_rows.sort(key=lambda r: -r["coverage_pct"])

    # Datasets per partition — sourced from shipped trees, so every entry maps
    # to a real .nwk the chooser can switch to. (For partitions where the
    # release ships a single canonical tree but provenance documents multiple
    # alternative source studies, only the shipped tree appears here.)
    datasets_by_partition: dict[str, list[dict]] = defaultdict(list)
    for t in trees:
        partition = t.get("partition_group")
        if not partition:
            continue
        datasets_by_partition[partition].append({
            "filename": t["filename"],
            "group": t["group"],
            "provenance_group": t["provenance_group"],
            "category": t.get("category"),
            "study": t["study"],
            "year": t["year"],
            "tips": t["ntips"],
            "dated": t["dated"],
            "is_partition_anchor": t["is_partition_anchor"],
            "is_partition_canonical": t["is_partition_canonical"],
        })
    # Sort each partition's datasets: canonical first, then sub-clades by tips desc.
    for partition in datasets_by_partition:
        datasets_by_partition[partition].sort(
            key=lambda x: (not x["is_partition_canonical"], -(x["tips"] or 0))
        )

    # "Not yet represented" lineages — partitions present in data_estimates.csv
    # but with no shipped tree. These are the paper's dark-matter clades.
    # Conifers is excluded: it's a sub-clade of Seed plants (the 615 Conifer
    # species double-count within Seed plants' 383k), so the audit reclassifies
    # Conifers as a sub-phylo rather than a top-level partition.
    mapped_partitions = {t["partition_group"] for t in trees if t.get("partition_group")}
    PARTITION_EXCLUDED_FROM_CHART = {"Conifers"}
    unrepresented = []
    for r in estimates_rows:
        g = r.get("group")
        if not g or g in mapped_partitions or g in PARTITION_EXCLUDED_FROM_CHART:
            continue
        unrepresented.append({
            "group": g,
            "category": r.get("category") or None,
            "described": parse_int(r.get("described")),
            "estimated_total": parse_int(r.get("estimated_total")),
            "estimated_low": parse_int(r.get("estimated_low")),
            "estimated_high": parse_int(r.get("estimated_high")),
            "estimate_source": r.get("estimate_source") or None,
            "estimate_source_info": source_info(r.get("estimate_source")),
            "estimate_confidence": r.get("confidence") or None,
        })
    # Stable order: by category then by described count desc.
    CAT_ORDER = {"Vertebrates": 0, "Plants": 1, "Arthropods": 2, "Other animals": 3,
                 "Microbes & protists": 4}
    unrepresented.sort(key=lambda x: (CAT_ORDER.get(x["category"], 9),
                                      -(x["described"] or 0)))

    # Synthetic atlas-only partitions: appear in audits/ but NOT in
    # data_estimates.csv. Currently 'Others' — aggregates Tuatara, Coelacanths,
    # Lungfishes, Onychophora, Placozoa, Loricifera, Cycliophora, Chaetognatha
    # so the atlas's partition is collectively exhaustive at the species level.
    synthetic_partitions = []
    existing_in_csv = {r["group"] for r in estimates_rows}
    existing_unrep = {u["group"] for u in unrepresented}
    # We need the audits dict; load it now (was below, but coverage_rows
    # need it here). Idempotent if already loaded.
    _audits_for_synth = load_audit_overrides()
    for partition_name, audit in _audits_for_synth.items():
        if partition_name in existing_in_csv or partition_name in existing_unrep:
            continue
        est = audit.get("estimate") or {}
        if not est.get("described") and not est.get("total"):
            continue  # no numeric data
        synthetic_partitions.append({
            "group": partition_name,
            "category": audit.get("category") or "Not yet represented",
            "described": est.get("described"),
            "estimated_total": est.get("total"),
            "estimated_low": est.get("low"),
            "estimated_high": est.get("high"),
            "estimate_source": (est.get("source") or {}).get("key"),
            "estimate_source_info": None,
            "estimate_confidence": est.get("confidence"),
            "synthetic": True,
        })

    # Surface the unrepresented + synthetic partitions on the main coverage
    # chart — gives users the full "dark matter" story alongside the shipped
    # trees. Each gets tips=null (no shipped tree → no ◇/★ marker), but the
    # ◯ described and ● estimated whiskers all render. Category is taken
    # from the partition's audit info.yaml (if present), allowing the
    # 'Not yet represented' CSV category to be overridden into the row's
    # natural taxonomic category (Acari → Arthropods, etc.).
    audit_categories = {p: a.get("category") for p, a in _audits_for_synth.items()}
    for u in list(unrepresented) + synthetic_partitions:
        u_cat = audit_categories.get(u["group"]) or u["category"]
        # If the synthetic partition has its own info.yaml tree (e.g.,
        # Xiphosura — hand-built 4-species dated chronogram with no upstream
        # standardized counterpart), promote it from "unrepresented" to a
        # full canonical row with tip count + dated flag from info.yaml.
        partition_audit = _audits_for_synth.get(u["group"]) or {}
        tree_meta = partition_audit.get("tree") or {}
        has_synth_tree = isinstance(tree_meta, dict) and tree_meta.get("file")
        if has_synth_tree:
            partition_slug = u["group"].lower().replace(" ", "_")
            tree_path = SITE_DATA / partition_slug / tree_meta["file"]
            if tree_path.exists():
                row_out = {
                    "filename": f"{partition_slug}.nwk",
                    "group": partition_slug,
                    "provenance_group": partition_slug,
                    "partition_group": u["group"],
                    "category": u_cat,
                    "study": (tree_meta.get("source") or {}).get("study") or (tree_meta.get("source") or {}).get("key"),
                    "year": None,
                    "tips": tree_meta.get("ntips"),
                    "raw_tips": tree_meta.get("ntips"),
                    "described_species": u["described"],
                    "described_source": None,
                    "coverage_pct": tree_meta.get("coverage_pct_described"),
                    "dated": bool(tree_meta.get("dated")),
                    "is_partition_anchor": False,
                    "is_partition_canonical": True,
                    "is_unrepresented": False,
                    "tree_url_override": str(tree_path.relative_to(SITE_DATA.parent.parent)).replace("\\", "/"),
                    "estimated_total": u["estimated_total"],
                    "estimated_low": u["estimated_low"],
                    "estimated_high": u["estimated_high"],
                    "estimate_source": u["estimate_source"],
                    "estimate_confidence": u["estimate_confidence"],
                }
                if u["estimated_total"]:
                    row_out["coverage_pct_estimated"] = _est_pct(tree_meta["ntips"], u["estimated_total"])
                coverage_rows.append(row_out)
                # Also push into the trees[] list so selectTree() in the viewer
                # can locate the tree by filename (the viewer reads from trees,
                # not coverage). Keys mirror what build_trees creates above.
                source = (tree_meta.get("source") or {})
                trees.append({
                    "filename": f"{partition_slug}.nwk",
                    "group": partition_slug,
                    "provenance_group": partition_slug,
                    "partition_group": u["group"],
                    "category": u_cat,
                    "is_partition_anchor": False,
                    "study": source.get("study") or source.get("key") or "",
                    "ntips": tree_meta.get("ntips"),
                    "dated": bool(tree_meta.get("dated")),
                    "size_bytes": tree_path.stat().st_size,
                    "year": None,
                    "journal": None,
                    "doi": source.get("doi"),
                    "crown_ma": None,
                    "described_species": u["described"],
                    "described_source": None,
                    "coverage_pct": tree_meta.get("coverage_pct_described"),
                    "data_source": None,
                    "download_url": source.get("url"),
                    "methods_brief": (tree_meta.get("methods") or {}).get("note"),
                    "notes": None,
                    "estimate": {
                        "category": u_cat,
                        "estimated_total": u["estimated_total"],
                        "estimated_low": u["estimated_low"],
                        "estimated_high": u["estimated_high"],
                        "estimate_source": u["estimate_source"],
                        "estimate_confidence": u["estimate_confidence"],
                    },
                    "described_source_info": None,
                    "estimate_source_info": None,
                    "tip_taxonomy": None,
                    "tree_url_override": row_out["tree_url_override"],
                    "tip_taxonomy_override": None,
                    "is_partition_canonical": True,
                })
                continue
        # Default: no tree shipped — render as dark-matter bar.
        row_out = {
            "filename": None,
            "group": u["group"].lower().replace(" ", "_"),
            "provenance_group": None,
            "partition_group": u["group"],
            "category": u_cat,
            "study": None,
            "year": None,
            "tips": None,
            "raw_tips": None,
            "described_species": u["described"],
            "described_source": None,
            "coverage_pct": None,
            "dated": False,
            "is_partition_anchor": False,
            "is_partition_canonical": True,
            "is_unrepresented": True,
            "estimated_total": u["estimated_total"],
            "estimated_low": u["estimated_low"],
            "estimated_high": u["estimated_high"],
            "estimate_source": u["estimate_source"],
            "estimate_confidence": u["estimate_confidence"],
        }
        coverage_rows.append(row_out)

    summary = {
        "n_trees": len(trees),
        "n_dated": dated_count,
        "n_undated": len(trees) - dated_count,
        "n_groups": len({r["group"] for r in trees}),
        "n_datasets": len(provenance),
        "total_tips": sum((r["ntips"] or 0) for r in trees),
    }

    # Load per-partition audit overrides from site/data/partitions/<slug>/info.yaml.
    # Each loaded record contributes:
    #   - methods + caveats + audit linkage → audits[partition_name]
    #   - corrected estimate source when the paper-cited DOI is broken
    audits = load_audit_overrides()
    summary["n_audited"] = len(audits)
    summary["n_audits_in_progress"] = sum(
        1 for a in audits.values()
        if (a.get("audit") or {}).get("status") == "in_progress"
    )
    summary["n_audits_verified"] = sum(
        1 for a in audits.values()
        if (a.get("audit") or {}).get("status") == "verified"
    )
    summary["n_audits_deferred"] = sum(
        1 for a in audits.values()
        if (a.get("audit") or {}).get("status") == "deferred"
    )
    # Partitions (= chart rows where is_partition_canonical=true)
    summary["n_partitions"] = sum(1 for r in coverage_rows if r.get("is_partition_canonical"))
    # Atlas-derived canonical trees (Turtles, Bryozoa, Primates, Seed plants, Fungi)
    summary["n_atlas_derived"] = sum(
        1 for t in trees if t.get("tree_url_override")
    )
    # Interval-class distribution
    from collections import Counter
    cls_counts = Counter()
    for a in audits.values():
        ip = (a.get("estimate") or {}).get("interval_provenance") or {}
        cls_counts[ip.get("overall_classification") or "unaudited"] += 1
    summary["n_paper_published_intervals"] = cls_counts.get("paper-published-range", 0)
    summary["n_partly_heuristic_intervals"] = cls_counts.get("partly-heuristic", 0)
    summary["n_fully_heuristic_intervals"] = cls_counts.get("fully-heuristic", 0)
    # Total described species (summed across partitions, excluding overlap)
    summary["total_described_species"] = sum(
        (r.get("described_species") or 0)
        for r in coverage_rows if r.get("is_partition_canonical")
    )
    # Coverage: partitions that ship a tree (chart rows that are NOT
    # is_unrepresented). The complementary set is the "dark matter" partitions
    # where only an estimate is shown.
    summary["n_partitions_with_tree"] = sum(
        1 for r in coverage_rows
        if r.get("is_partition_canonical") and not r.get("is_unrepresented")
    )
    # Unique species represented across all shipped partition canonicals.
    # Sum of tips on canonical (with-tree) rows — partitions are mutually
    # exclusive, so this is an upper bound on actual unique species; close
    # enough for the dashboard headline.
    summary["unique_species_in_trees"] = sum(
        (r.get("tips") or 0)
        for r in coverage_rows
        if r.get("is_partition_canonical") and not r.get("is_unrepresented")
    )
    # Estimated true-diversity range across all partitions
    summary["estimated_total_low"] = sum(
        (r.get("estimated_low") or 0)
        for r in coverage_rows if r.get("is_partition_canonical")
    )
    summary["estimated_total_high"] = sum(
        (r.get("estimated_high") or 0)
        for r in coverage_rows if r.get("is_partition_canonical")
    )
    summary["estimated_total_central"] = sum(
        (r.get("estimated_total") or 0)
        for r in coverage_rows if r.get("is_partition_canonical")
    )

    # Apply broken-DOI substitution AND attach interval_provenance on coverage
    # rows when the audit has them.
    for row in coverage_rows:
        partition = row.get("partition_group")
        if not partition or partition not in audits:
            continue
        est_obj = audits[partition].get("estimate") or {}
        src = est_obj.get("source") or {}
        if src.get("paper_doi_status") == "broken" and src.get("live_doi"):
            row["estimate_source_corrected"] = {
                "key": src.get("live_key"),
                "doi": src.get("live_doi"),
                "year": src.get("live_year"),
                "url": src.get("url"),
                "note": src.get("note"),
            }
            row["estimate_source_paper_doi_status"] = "broken"
        if est_obj.get("interval_provenance"):
            row["interval_provenance"] = est_obj["interval_provenance"]

    bundle = {
        "schema_version": 4,
        "generated_from": "standardized/metadata.csv + data_provenance.csv + data_estimates.csv + data_sources.csv + site/data/partitions/*/info.yaml",
        "summary": summary,
        "trees": trees,
        "coverage": coverage_rows,
        "datasets_by_partition": dict(datasets_by_partition),
        "unrepresented": unrepresented,
        "audits": audits,
    }

    OUT.write_text(json.dumps(bundle, indent=2, ensure_ascii=False, default=_json_default))
    size_kb = OUT.stat().st_size / 1024
    print(f"wrote {OUT.name}: {len(trees)} trees, {len(coverage_rows)} coverage rows, "
          f"{len(unrepresented)} unrepresented partitions, {size_kb:.1f} KB", file=sys.stderr)

    # Name shards were already built (and counted) at the top of main().
    if shard_bytes:
        print(f"wrote name shards: {len(names_by_filename)} files, {shard_bytes/1024/1024:.1f} MB total",
              file=sys.stderr)


if __name__ == "__main__":
    main()
