"""Build site/data.json, the data of the browser site, from this release.

    python3 site/build.py            # write site/data.json and the name shards
    python3 site/build.py --check    # compare with the committed data.json

Every number the site shows is read from the files of the paper in this
release and copied, never typed:

    registers/coverage_basis.csv     groups, described species (the paper's
                                     denominators), the reference tree of each
                                     group, named and molecular species, time
                                     scale and date uncertainty of each tree
    paper/headline_numbers.json      the totals of the four coverage levels
    paper/crown_ages.csv             crown ages and intervals of the reference
                                     trees (the crown register)
    paper/table_datasets.csv         tips, accepted names, described species
                                     and named percentage of the other trees
    paper/fig1_coverage_bars.csv     estimated totals of each group
    registers/dataset_table_registry.csv, registers/tree_provenance.csv
                                     year, journal and DOI of each tree
    standardized/metadata.csv        the redistributed trees; each tree file is
                                     read from standardized/trees/

Percentages are formatted as the paper prints them: one decimal, two below 1.
The per-group notes in site/data/partitions/*/info.yaml supply text only
(methods, caveats, the provenance of each estimate) and the trees that exist
only in the site (sub-clade trees that are not redistributed, and a tree
derived from a redistributed one).  data_provenance.csv supplies only the
download link and DOI of the trees that are not reference trees.

The page reads every reference tree and every other redistributed tree from
standardized/trees/; site/data/partitions holds no copy of a reference tree.
Requires PyYAML.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
SITE = Path(__file__).resolve().parent
STD = ROOT / "standardized"
TREES = STD / "trees"
NAMES_DIR = STD / "names"
PARTITIONS = SITE / "data" / "partitions"
OUT = SITE / "data.json"

INPUTS = [
    "registers/coverage_basis.csv", "paper/headline_numbers.json",
    "paper/crown_ages.csv", "paper/table_datasets.csv",
    "paper/fig1_coverage_bars.csv", "registers/dataset_table_registry.csv",
    "registers/tree_provenance.csv", "standardized/metadata.csv",
    "standardized/withdrawn.csv", "data_provenance.csv", "CITATION.cff",
    "site/data/partitions/*/info.yaml",
]

LEVELS = [
    # key in headline_numbers.json rungs, key in data.json, label on the page
    ("C0", "named", "Named in a reference tree"),
    ("C1", "molecular", "Placed with molecular data"),
    ("C1d", "dated", "In a chronogram in absolute time"),
    ("C1u", "uncertainty", "Date uncertainty available"),
]

TIME_SCALE_TEXT = {
    "absolute": "chronogram in absolute time",
    "relative": "ultrametric, relative time only",
    "not_ultrametric": "not a chronogram",
}
INTERVAL_KIND_TEXT = {
    "posterior_sample": "posterior tree set",
    "bootstrap_dating_sample": "bootstrap dating tree set",
    "date_assignment_sample": "tree set with dates drawn at random from published estimates",
    "calibration_prior_sample": "tree set dated under root ages drawn from the calibration prior",
    "node_hpd_annotation": "node age intervals on the summary tree",
    "node_label_interval": "node age intervals on the summary tree",
    "taxon_addition_sample": "taxon-addition tree set: placement of added tips only, no date uncertainty",
}

TIP_ID_RE = re.compile(r"([(,])(-?\d+)(?=[:,)])")
TIP_RE = re.compile(r"[(,]\s*([^(),:;\s]+)")


def read_csv(rel: str) -> list[dict]:
    with (ROOT / rel).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def as_int(v):
    try:
        return int(round(float(v)))
    except (TypeError, ValueError):
        return None


def thousands(n) -> str | None:
    return None if n is None else "{:,}".format(int(round(float(n))))


def pct_text(x) -> str | None:
    """The paper's rule: one decimal; values below 1 get two decimals."""
    if x is None:
        return None
    return "%.1f" % x if abs(x) >= 1.0 else "%.2f" % x


def pct(num, den):
    if num is None or not den:
        return None
    return 100.0 * float(num) / float(den)


def ma_text(v) -> str | None:
    f = None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return "%.1f" % f


def tree_tip_labels(path: Path) -> list[str]:
    return TIP_RE.findall(path.read_text(encoding="utf-8", errors="ignore"))


def load_yaml(path: Path) -> dict:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


_FRONTMATTER = re.compile(r"^---\s*\n(.*?\n)---\s*$", re.DOTALL | re.MULTILINE)


def load_notes() -> tuple[dict, dict, list]:
    """Per-group notes: {partition: info.yaml (+ audit frontmatter)},
    {slug: partition}, and the sub-clade folders [(partition, slug, dir, info)]."""
    audits_md = {}
    adir = ROOT / "audits" / "partitions"
    if adir.exists():
        for md in sorted(adir.glob("*.md")):
            m = _FRONTMATTER.match(md.read_text(encoding="utf-8"))
            if not m:
                continue
            d = yaml.safe_load(m.group(1))
            if isinstance(d, dict) and d.get("partition"):
                audits_md[d["partition"]] = d
    notes, slugs, subs = {}, {}, []
    for info in sorted(PARTITIONS.glob("*/info.yaml")):
        d = load_yaml(info)
        part = d.get("partition")
        if not part:
            continue
        slugs[info.parent.name] = part
        for k, v in (audits_md.pop(part, None) or {}).items():
            d.setdefault(k, v)
        # A note may point at a file the site no longer holds (a former copy
        # of a reference tree); such entries are not shown.
        unc = d.get("uncertainty")
        if isinstance(unc, dict) and unc.get("files"):
            unc["files"] = [f for f in unc["files"]
                            if not (f.get("local") and not (SITE / f.get("url", "")).exists())]
        # The tree block of a group note described the site's former copy of
        # the reference tree; the reference tree is described by the paper.
        tree = d.pop("tree", None)
        if isinstance(tree, dict):
            d["tree"] = {k: tree[k] for k in ("methods", "source", "superseded_canonical")
                         if k in tree}
        d.pop("species_count", None)
        notes[part] = d
        for sub in sorted((info.parent / "sub-phylos").glob("*/info.yaml")):
            subs.append((part, sub.parent.name, sub.parent, load_yaml(sub)))
    return notes, slugs, subs


def build() -> dict:
    basis = read_csv("registers/coverage_basis.csv")
    head = json.loads((ROOT / "paper" / "headline_numbers.json").read_text(encoding="utf-8"))
    crowns = {r["partition"]: r for r in read_csv("paper/crown_ages.csv")}
    bars = {r["group"]: r for r in read_csv("paper/fig1_coverage_bars.csv")}
    dsets = {r["filename"]: r for r in read_csv("paper/table_datasets.csv") if r["filename"]}
    registry = {r["filename"]: r for r in read_csv("registers/dataset_table_registry.csv")}
    tprov = {r["filename"]: r for r in read_csv("registers/tree_provenance.csv")}
    metadata = read_csv("standardized/metadata.csv")
    withdrawn = {r["filename"] for r in read_csv("standardized/withdrawn.csv")}
    prov = defaultdict(list)
    for r in read_csv("data_provenance.csv"):
        prov[r.get("group") or ""].append(r)
    notes, slugs, subs = load_notes()
    slug_of = {v: k for k, v in slugs.items()}

    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    m = re.search(r"preferred-citation:.*?\n\s+title:\s*\"([^\"]*)\"", cff, re.S)
    paper_title = m.group(1) if m else None
    m = re.search(r"preferred-citation:.*?\n\s+doi:\s*\"?([^\"\s]+)", cff, re.S)
    paper_doi = m.group(1) if m else None

    # ---------------------------------------------------------- groups ---
    rows = []
    canon_of = {}
    for b in basis:
        part = b["partition"]
        described = as_int(b["described"])
        ref = b["canonical_tree"] or None
        prok = b["is_prokaryotic"] == "yes"
        row = {
            "partition_group": part,
            "group": part,
            "slug": slug_of.get(part),
            "category": b["category"],
            "is_prokaryotic": prok,
            "described_species": described,
            "described_text": thousands(described),
            "described_source": b["described_source"] or None,
            "filename": ref,
            "is_partition_canonical": True,
            "is_unrepresented": ref is None,
        }
        if ref:
            canon_of[ref] = part
            c0, c1 = as_int(b["c0_units"]), as_int(b["c1_units"])
            dated = b["is_chronogram"] == "yes"
            unc = dated and b["uncertainty_available"] == "yes"
            levels = {
                "named": c0, "molecular": c1,
                "dated": c1 if dated else 0, "uncertainty": c1 if unc else 0,
            }
            row.update({
                "tips": as_int(b["labels_unique"]),
                "tips_text": thousands(b["labels_unique"]),
                "time_scale": b["time_scale"] or None,
                "time_scale_text": TIME_SCALE_TEXT.get(b["time_scale"]),
                "dated": dated,
                "uncertainty_available": unc,
                "interval_kind": b["interval_kind"] or None,
                "interval_kind_text": INTERVAL_KIND_TEXT.get(b["interval_kind"]),
                "molecular_basis": b["c1_basis"] or None,
                "levels": {k: {"n": v, "n_text": thousands(v),
                               "pct_text": pct_text(pct(v, described))}
                           for k, v in levels.items()},
            })
            cr = crowns.get(part)
            if cr:
                row["crown"] = {
                    "clade": cr["crown_clade"],
                    "age_ma_text": ma_text(cr["crown_age_ma"]),
                    "age_statistic": cr["age_statistic"],
                    "interval_lo_ma_text": ma_text(cr["interval_lo_ma"]) if unc else None,
                    "interval_hi_ma_text": ma_text(cr["interval_hi_ma"]) if unc else None,
                    "interval_statistic": cr["center_statistic"] if unc else None,
                    "n_trees": as_int(cr["n_trees"]) if unc else None,
                }
        est = bars.get(part) or {}
        if as_int(est.get("estimated_total")):
            row.update({
                "estimated_total": as_int(est["estimated_total"]),
                "estimated_low": as_int(est.get("estimated_low")),
                "estimated_high": as_int(est.get("estimated_high")),
                "estimate_source": est.get("estimate_source") or None,
                "estimate_confidence": est.get("confidence") or None,
            })
        ip = ((notes.get(part) or {}).get("estimate") or {}).get("interval_provenance")
        if ip:
            row["interval_provenance"] = ip
        rows.append(row)
    by_part = {r["partition_group"]: r for r in rows}

    # ----------------------------------------------------------- trees ---
    sub_of = {}            # metadata group -> (partition, folder, info)
    for part, slug, folder, info in subs:
        sub_of[slug] = (part, folder, info)

    def study_fields(fn, group):
        reg = registry.get(fn) or {}
        tp = tprov.get(fn) or {}
        pv = (prov.get(group) or [{}])[0]
        doi = tp.get("doi") or pv.get("doi") or None
        return {
            "year": as_int(reg.get("year")) or as_int(pv.get("year")),
            "journal": reg.get("journal") or None,
            "doi": doi,
            "download_url": pv.get("download_url") or None,
            "data_source": pv.get("data_source") or None,
        }

    trees = []
    for md in metadata:
        fn, group = md["filename"], md["group"]
        if fn in withdrawn:
            raise SystemExit("standardized/metadata.csv lists withdrawn %s" % fn)
        path = TREES / fn
        labels = tree_tip_labels(path)
        t = {
            "filename": fn, "group": group, "study": md["study"],
            "ntips": len(labels), "ntips_text": thousands(len(labels)),
            "size_bytes": path.stat().st_size,
            "tree_url": "standardized/trees/" + fn,
            "names_file": fn,
            "in_package": True,
        }
        t.update(study_fields(fn, group))
        if fn in canon_of:
            r = by_part[canon_of[fn]]
            t.update({"partition_group": r["partition_group"], "category": r["category"],
                      "is_partition_canonical": True, "dated": r["dated"],
                      "dated_text": r["time_scale_text"]})
        else:
            part = None
            if group in sub_of:
                part = sub_of[group][0]
            elif group.startswith("condamine_"):
                part = "Tetrapod families (Condamine et al. 2019)"
            ds = dsets.get(fn) or {}
            t.update({
                "partition_group": part,
                "category": (by_part.get(part) or {}).get("category"),
                "is_partition_canonical": False,
                "dated": md["dated"].strip().upper() == "TRUE",
                "dated_text": ds.get("dated") or None,
            })
            if ds:
                t["paper"] = {
                    "accepted_names_text": thousands(ds["accepted_names"]) if ds["accepted_names"] else None,
                    "described_text": thousands(ds["described"]) if ds["described"] else None,
                    "named_pct_text": ds["coverage_pct"] or None,
                    "described_source": ds["described_source"] or None,
                }
        trees.append(t)

    # Trees that only the site holds: sub-clade trees that are not
    # redistributed, and a tree derived from a redistributed one.
    std_names = {md["filename"] for md in metadata}
    for part, slug, folder, info in subs:
        f = folder / "tree.nwk"
        if not f.exists():
            continue
        std = TREES / (slug + ".nwk")
        rel = f.relative_to(SITE).as_posix()
        tb = info.get("tree") or {}
        src = tb.get("source") or {}
        if (slug + ".nwk") in std_names:
            if f.read_bytes() == std.read_bytes():
                continue          # the same file as the redistributed tree
            if not tb.get("derivation"):
                raise SystemExit("%s differs from %s and states no derivation"
                                 % (rel, std.relative_to(ROOT)))
            fn, names_file = slug + "_one_tip_per_species.nwk", slug + ".nwk"
            kind, label = "derived", slug + " (one tip per species)"
        else:
            fn, names_file, kind, label = slug + ".nwk", None, "site_only", slug
        labels = tree_tip_labels(f)
        trees.append({
            "filename": fn, "group": label, "study": src.get("key") or src.get("study") or "",
            "ntips": len(labels), "ntips_text": thousands(len(labels)),
            "size_bytes": f.stat().st_size,
            "tree_url": "site/" + rel, "tree_url_override": rel,
            "names_file": names_file, "in_package": False, "site_tree": kind,
            "derivation": tb.get("derivation"), "derived_from": names_file,
            "year": None, "journal": None, "doi": src.get("doi"),
            "download_url": src.get("url"), "data_source": None,
            "partition_group": part, "category": (by_part.get(part) or {}).get("category"),
            "is_partition_canonical": False, "dated": bool(tb.get("dated")),
            "dated_text": None, "methods": tb.get("methods"),
        })

    trees.sort(key=lambda t: (-(t["ntips"] or 0), t["filename"]))
    datasets = defaultdict(list)
    for t in trees:
        if t.get("partition_group"):
            datasets[t["partition_group"]].append(
                {k: t.get(k) for k in ("filename", "group", "study", "year", "ntips",
                                       "dated", "is_partition_canonical")})
    for v in datasets.values():
        v.sort(key=lambda x: (not x["is_partition_canonical"], -(x["ntips"] or 0)))

    # --------------------------------------------------------- summary ---
    den = head["denominators"]
    euk_den = den["eukaryotic"]["described"]
    levels = []
    for rung, key, label in LEVELS:
        r = head["rungs"][rung]
        levels.append({
            "key": key, "label": label, "n": r["numerator_union"],
            "n_text": thousands(r["numerator_union"]),
            "pct_text": pct_text(r["percent_union"]),
            "n_groups": r["n_partitions_contributing"],
        })
    summary = {
        "n_trees": len(metadata),
        "n_site_trees": sum(1 for t in trees if not t["in_package"]),
        "n_partitions": len(rows),
        "n_partitions_eukaryotic": den["eukaryotic"]["n_partitions"],
        "n_partitions_prokaryotic": den["prokaryotic"]["n_partitions"],
        "n_partitions_with_tree": sum(1 for r in rows if not r["is_unrepresented"]),
        "n_partitions_with_tree_eukaryotic": sum(
            1 for r in rows if not r["is_unrepresented"] and not r["is_prokaryotic"]),
        "described_eukaryotic": euk_den,
        "described_eukaryotic_text": thousands(euk_den),
        "described_prokaryotic_text": thousands(den["prokaryotic"]["described"]),
        "estimated_total_low_text": thousands(den["all_life"]["estimated_low"]),
        "estimated_total_high_text": thousands(den["all_life"]["estimated_high"]),
        "n_audited": len(notes),
    }
    return {
        "schema_version": 5,
        "generated_from": INPUTS,
        "paper": {"title": paper_title, "doi": paper_doi},
        "summary": summary,
        "levels": levels,
        "coverage": rows,
        "trees": trees,
        "datasets_by_partition": dict(datasets),
        "audits": notes,
    }


def _json_default(o):
    import datetime
    if isinstance(o, (datetime.date, datetime.datetime)):
        return o.isoformat()
    raise TypeError("not JSON serializable: %r" % type(o))


def dump(bundle) -> str:
    return json.dumps(bundle, indent=1, ensure_ascii=False, default=_json_default) + "\n"


def write_name_shards(trees) -> int:
    """standardized/names/<file>.json: {tip id: standardized name} per tree,
    read by the tree viewer (not committed; rebuilt here)."""
    dictionary = STD / "dictionary.csv"
    names = {}
    with dictionary.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            names[r["id"]] = r["standardized_name"]
    NAMES_DIR.mkdir(exist_ok=True)
    n = 0
    for t in trees:
        if not t["in_package"]:
            continue
        text = (TREES / t["filename"]).read_text(encoding="utf-8", errors="ignore")
        ids = {m.group(2) for m in TIP_ID_RE.finditer(text)}
        shard = {i: names[i] for i in sorted(ids) if i in names}
        (NAMES_DIR / (t["filename"] + ".json")).write_text(
            json.dumps(shard, separators=(",", ":"), ensure_ascii=False), encoding="utf-8")
        n += 1
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="exit 1 if site/data.json differs from a fresh build")
    ap.add_argument("--no-names", action="store_true",
                    help="do not rebuild standardized/names/")
    a = ap.parse_args(argv)
    bundle = build()
    text = dump(bundle)
    if a.check:
        old = OUT.read_text(encoding="utf-8") if OUT.exists() else None
        if old != text:
            print("site/data.json differs from what site/build.py builds", file=sys.stderr)
            return 1
        print("site/data.json is current", file=sys.stderr)
        return 0
    OUT.write_text(text, encoding="utf-8")
    print("wrote site/data.json: %d groups, %d trees (%d only in the site)"
          % (len(bundle["coverage"]), len(bundle["trees"]),
             bundle["summary"]["n_site_trees"]), file=sys.stderr)
    if not a.no_names:
        print("wrote %d name shards in standardized/names/" % write_name_shards(bundle["trees"]),
              file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
