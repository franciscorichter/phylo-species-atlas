#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""match_col.py -- resolve every eukaryotic canonical tip to a Catalogue of
Life accepted name (finding S05; decisions D4, D9, D13, D18).

Authority: the Catalogue of Life BASE release COL26.9 (ChecklistBank dataset
316321, doi:10.48580/dgz5n, released 2026-09-11), read from its DwC-A export
https://download.checklistbank.org/col/monthly/2026-09-11_dwca.zip
(sha256 d50f75ae601e6657f08f45144226f3187d893f42aa74d7a99836e345294787f5).
The release is not distributed with the atlas; this script reads a compact
index of its Taxon.tsv (every row with a specific epithet), built by --index.

RULE, applied to every tip of every eukaryotic canonical tree:
  1. Query.  The species key of the standardized name (D4: the first
     Genus_epithet pair; scripts/canonical_trees.py species_key), or the
     query in Data/standardized/name_overrides.csv for a label whose format
     carries no binomial (each row cites its evidence).  A label that is not a
     species label (D13: sp., cf., aff., hybrid) is 'unresolved_label'.
  2. Accepted.  Species-rank rows of the release with that genus and epithet
     and status accepted or provisionally accepted.  If several (homonyms in
     different kingdoms), the one in the partition's kingdom is taken; the
     partition's kingdom is the majority kingdom of its unambiguous accepted
     matches (measured, not declared).
  3. Synonym.  Otherwise, synonym and ambiguous-synonym rows with that genus
     and epithet; each points to an accepted name (an infraspecific target is
     taken at its species).  Targets outside the partition's kingdom are
     dropped.  One target: resolved to it.  More than one: 'ambiguous'.
  4. A comparison that ignores hyphens and case is tried when 2-3 find
     nothing ('Pseudonitzschia' / 'Pseudo-nitzschia').
  5. Otherwise 'unmatched'.
A resolved tip counts under its accepted name; tips resolving to one accepted
name count once (C0, decision on S05).

OUTPUT: Data/standardized/col_resolution.csv, one row per (partition, tip);
Data/standardized/col_resolution_noncanonical.csv, the same rule applied to
every non-canonical eukaryotic tree, one row per (file, tip) (finding S09).

USAGE
    python scripts/match_col.py --index DWCA_ZIP     # build the index once
    python scripts/match_col.py                       # resolve
    python scripts/match_col.py --family-counts       # species per family
    python scripts/match_col.py --opentree            # P01: Open Tree names

--opentree applies the same rule to the Open Tree of Life species placed by
phylogeny (Data/external/opentree16.1_phylo_species.tsv.gz, written by
opentree_extract.py), one call per atlas clade so that the homonym kingdom is
measured per clade as for a tree, and writes
Data/external/opentree16.1_phylo_species_col.csv (ott_id, clade, query,
match_type, col_id).  rebuild_atlas.py then compares the two sources on
accepted names as well as on exact spellings (finding R1-14).
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import sys
import zipfile
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, HERE)
from canonical_trees import species_key  # noqa: E402

DATA = os.path.join(ROOT, "Data")
STD = os.path.join(DATA, "standardized")
INDEX = os.path.join(DATA, "raw", "col", "col269_species_index.tsv")
OUT = os.path.join(STD, "col_resolution.csv")
OUT_NONCANON = os.path.join(STD, "col_resolution_noncanonical.csv")
OVERRIDES = os.path.join(STD, "name_overrides.csv")
OT_SPECIES = os.path.join(DATA, "external", "opentree16.1_phylo_species.tsv.gz")
OT_OUT = os.path.join(DATA, "external", "opentree16.1_phylo_species_col.csv")
RELEASE = "COL26.9 base release (ChecklistBank 316321, doi:10.48580/dgz5n)"

ACCEPTED = ("accepted", "provisionally accepted")
SYNONYM = ("synonym", "ambiguous synonym")

csv.field_size_limit(1 << 30)


def build_index(zip_path):
    os.makedirs(os.path.dirname(INDEX), exist_ok=True)
    z = zipfile.ZipFile(zip_path)
    n = 0
    with z.open("Taxon.tsv") as fh, open(INDEX, "w", encoding="utf-8") as out:
        t = io.TextIOWrapper(fh, encoding="utf-8", newline="")
        t.readline()
        for line in t:
            f = line.rstrip("\n").split("\t")
            if len(f) < 26 or not f[13]:
                continue
            # taxonID acceptedID status rank genus epithet infra kingdom
            # family scientificName
            out.write("\t".join([f[0], f[2], f[6], f[7], f[11], f[13], f[14],
                                 f[20], f[25], f[8]]) + "\n")
            n += 1
    print("indexed %d rows into %s" % (n, os.path.relpath(INDEX, ROOT)))


def suprageneric_genera():
    """Genus names of the release's species-rank rows (any status) that end in
    a suprageneric suffix of canonical_trees.SUPRAGENERIC."""
    from canonical_trees import SUPRAGENERIC
    out = set()
    with open(INDEX, encoding="utf-8") as fh:
        for line in fh:
            g = line.split("\t", 5)[4]
            if g.lower().endswith(SUPRAGENERIC):
                out.add(g)
    return sorted(out)


def write_suprageneric_genera():
    from canonical_trees import SUPRAGENERIC_GENERA_FILE
    gen = suprageneric_genera()
    with open(SUPRAGENERIC_GENERA_FILE, "w", encoding="utf-8") as fh:
        fh.write("# Genera of %s whose names end in a suprageneric suffix;\n"
                 "# written by scripts/match_col.py --genera from %s\n"
                 % (RELEASE, os.path.relpath(INDEX, ROOT)))
        fh.write("\n".join(gen) + "\n")
    print("wrote %d genera to %s"
          % (len(gen), os.path.relpath(SUPRAGENERIC_GENERA_FILE, ROOT)))


FAMILY_COUNTS = os.path.join(STD, "col269_family_species_counts.csv")


def write_family_counts():
    """Accepted species per (kingdom, family) in the release: species-rank
    rows with status accepted or provisionally accepted, no infraspecific
    epithet, counted once per (genus, epithet).  The denominators of the
    family-specific sampling fractions of the worked example (audit A15)."""
    keys = defaultdict(set)
    with open(INDEX, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[2] in ACCEPTED and f[3] == "species" and not f[6] and f[8]:
                keys[(f[7], f[8])].add((f[4], f[5]))
    with open(FAMILY_COUNTS, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["kingdom", "family", "n_accepted_species"])
        for (k, fam), v in sorted(keys.items()):
            w.writerow([k, fam, len(v)])
    print("wrote %d families to %s"
          % (len(keys), os.path.relpath(FAMILY_COUNTS, ROOT)))


def load_index():
    by_key = defaultdict(list)       # (genus, epithet) -> rows
    by_loose = defaultdict(list)     # hyphen/case-insensitive key -> rows
    by_id = {}
    with open(INDEX, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            row = dict(id=f[0], acc=f[1], status=f[2], rank=f[3], genus=f[4],
                       epithet=f[5], infra=f[6], kingdom=f[7], family=f[8],
                       name=f[9])
            by_id[row["id"]] = row
            if row["infra"]:
                continue
            by_key[(row["genus"], row["epithet"])].append(row)
            by_loose[loose(row["genus"], row["epithet"])].append(row)
    return by_key, by_loose, by_id


def loose(g, e):
    return (g.replace("-", "").lower(), e.replace("-", "").lower())


def species_of(row, by_key, by_id):
    """The accepted species-rank row a (possibly infraspecific) accepted row
    belongs to."""
    if row["rank"] == "species":
        return row
    for r in by_key.get((row["genus"], row["epithet"]), []):
        if r["status"] in ACCEPTED and r["rank"] == "species":
            return r
    return None


def candidates(rows, by_key, by_id):
    acc = [r for r in rows if r["status"] in ACCEPTED and r["rank"] == "species"]
    syn = []
    for r in rows:
        if r["status"] in SYNONYM and r["acc"] in by_id:
            t = species_of(by_id[r["acc"]], by_key, by_id)
            if t is not None:
                syn.append(t)
    return acc, syn


def resolve_ids(ids, names, override, by_key, by_loose, by_id):
    """Resolve one tree's ids by the rule in the module docstring.  The
    kingdom used to separate homonyms is the majority kingdom of the tree's
    unambiguous accepted matches.  Returns (rows, kingdom)."""
    queries = {}
    for sid in ids:
        if sid in override:
            queries[sid] = (override[sid], "override")
            continue
        k, ok = species_key(names[sid], sid)
        queries[sid] = (k, "species_key") if ok else (None, "")
    kcount = Counter()
    cand = {}
    for sid, (q, _how) in queries.items():
        if q is None:
            continue
        g, e = q.split("_", 1)
        rows = by_key.get((g, e)) or []
        acc, syn = candidates(rows, by_key, by_id)
        loosed = False
        if not acc and not syn:
            acc, syn = candidates(by_loose.get(loose(g, e)) or [],
                                  by_key, by_id)
            loosed = bool(acc or syn)
        cand[sid] = (acc, syn, loosed)
        if len({r["id"] for r in acc}) == 1:
            kcount[acc[0]["kingdom"]] += 1
    kingdom = kcount.most_common(1)[0][0] if kcount else ""
    out = []
    for sid in ids:
        q, how = queries[sid]
        row = dict(standardized_id=sid, query=q or "", query_source=how,
                   match_type="", col_id="", col_name="", col_kingdom="",
                   col_family="")
        if q is None:
            row["match_type"] = "unresolved_label"
            out.append(row)
            continue
        acc, syn, loosed = cand[sid]
        pick, mt = None, None
        accs = {r["id"]: r for r in acc}
        if len(accs) > 1:
            accs = {i: r for i, r in accs.items() if r["kingdom"] == kingdom}
        if len(accs) == 1:
            pick = list(accs.values())[0]
            mt = ("provisional" if pick["status"] != "accepted"
                  else "accepted")
        elif not accs:
            tg = {r["id"]: r for r in syn if r["kingdom"] == kingdom
                  or not kingdom}
            if len(tg) == 1:
                pick, mt = list(tg.values())[0], "synonym"
            elif len(tg) > 1:
                mt = "ambiguous"
            else:
                mt = "unmatched"
        else:
            mt = "ambiguous"
        if pick is not None and pick["kingdom"] != kingdom and kingdom:
            pick, mt = None, "other_kingdom"
        if pick is not None:
            if loosed:
                mt += "_orthographic"
            row.update(col_id=pick["id"], col_name=pick["name"],
                       col_kingdom=pick["kingdom"],
                       col_family=pick["family"])
        row["match_type"] = mt
        out.append(row)
    return out, kingdom


def write_opentree():
    """R1-14: resolve the Open Tree phylogeny-placed species to COL26.9
    accepted names by the rule of this module (no overrides: OTT names are
    binomials)."""
    import gzip
    by_clade = defaultdict(list)
    with gzip.open(OT_SPECIES, "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            ott, nm, clade = line.rstrip("\n").split("\t")
            by_clade[clade].append((ott, nm.replace(" ", "_")))
    by_key, by_loose, by_id = load_index()
    out = []
    for clade in sorted(by_clade):
        ids = [o for o, _n in by_clade[clade]]
        names = dict(by_clade[clade])
        rows, kingdom = resolve_ids(ids, names, {}, by_key, by_loose, by_id)
        for row in rows:
            out.append(dict(ott_id=row["standardized_id"], clade=clade,
                            query=row["query"], match_type=row["match_type"],
                            col_id=row["col_id"]))
        c = Counter(r["match_type"] for r in rows)
        print("%-14s kingdom %-9s %s" % (clade, kingdom, dict(sorted(c.items()))))
    with open(OT_OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, ["ott_id", "clade", "query", "match_type",
                                "col_id"], lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print("wrote %s (%d rows); authority %s"
          % (os.path.relpath(OT_OUT, ROOT), len(out), RELEASE))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--opentree", action="store_true",
                    help="resolve the Open Tree phylogeny-placed species "
                         "(Data/external/opentree16.1_phylo_species_col.csv)")
    ap.add_argument("--index", metavar="DWCA_ZIP")
    ap.add_argument("--genera", action="store_true",
                    help="write the CoL genera whose names end in a "
                         "suprageneric suffix (canonical_trees.py)")
    ap.add_argument("--family-counts", action="store_true",
                    help="write accepted species per family "
                         "(Data/standardized/col269_family_species_counts.csv)")
    a = ap.parse_args(argv)
    if a.opentree:
        write_opentree()
        return 0
    if a.family_counts:
        write_family_counts()
        return 0
    if a.index:
        build_index(a.index)
    if a.index or a.genera:
        write_suprageneric_genera()
        return 0
    if not os.path.exists(INDEX):
        sys.exit("no index at %s; run with --index <COL DwC-A zip>"
                 % os.path.relpath(INDEX, ROOT))
    names = {}
    with open(os.path.join(STD, "dictionary.csv"), newline="",
              encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            names[r["id"]] = r["standardized_name"]
    override = {}
    if os.path.exists(OVERRIDES):
        with open(OVERRIDES, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                override[r["standardized_id"].strip()] = r["query"].strip()
    prok = set()
    with open(os.path.join(DATA, "_rebuild", "truth_species.csv"), newline="",
              encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["is_prokaryotic"].strip().lower() in ("true", "yes", "1"):
                prok.add(r["partition"])
    with open(os.path.join(DATA, "canonical_sources.csv"), newline="",
              encoding="utf-8") as fh:
        canon = [(r["partition"], r["filename"]) for r in csv.DictReader(fh)
                 if r["partition"] not in prok]

    by_key, by_loose, by_id = load_index()
    out = []
    for part, fn in sorted(canon):
        with open(os.path.join(STD, "names", fn + ".json"),
                  encoding="utf-8") as fh:
            ids = sorted(json.load(fh), key=int)
        rows, kingdom = resolve_ids(ids, names, override, by_key, by_loose,
                                    by_id)
        for row in rows:
            row["partition"] = part
        out += rows
        c = Counter(r["match_type"] for r in rows)
        print("%-14s kingdom %-9s %s" % (part, kingdom, dict(sorted(c.items()))))
    hdr = ["partition", "standardized_id", "query", "query_source",
           "match_type", "col_id", "col_name", "col_kingdom", "col_family"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, hdr, lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print("wrote %s (%d rows); authority %s"
          % (os.path.relpath(OUT, ROOT), len(out), RELEASE))

    # S09: the same rule over every NON-canonical eukaryotic tree, one row per
    # (file, tip).  rebuild_atlas.py forms the union bound on C0 from it.
    noncanon = []
    with open(os.path.join(DATA, "_rebuild", "truth_dated.csv"), newline="",
              encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["is_canonical"].strip().lower() != "yes" \
                    and r["partition_name"] not in prok:
                noncanon.append(r["filename"])
    out2 = []
    for fn in sorted(noncanon):
        with open(os.path.join(STD, "names", fn + ".json"),
                  encoding="utf-8") as fh:
            ids = sorted(json.load(fh), key=int)
        rows, _k = resolve_ids(ids, names, override, by_key, by_loose, by_id)
        for row in rows:
            row["filename"] = fn
        out2 += rows
    hdr2 = ["filename"] + hdr[1:]
    with open(OUT_NONCANON, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, hdr2, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(out2)
    print("wrote %s (%d rows from %d non-canonical trees)"
          % (os.path.relpath(OUT_NONCANON, ROOT), len(out2), len(noncanon)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
