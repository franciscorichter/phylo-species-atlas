#!/usr/bin/env python3
"""curate.py -- the measured half of curating one candidate tree.

For a tree file a watch run flagged, it writes living/curation/<group>/<name>.json
and .md with what can be measured without reading the source study:

  tips, tip labels with a species key (canonical_trees.species_key, D4/D13);
  each species key matched in the Catalogue of Life COL26.9 base release
  (ChecklistBank 316321): accepted species -> itself, synonym -> its accepted
  species, anything else (genus-level match, no match) -> unresolved; the
  accepted species counted once (the C0 rule of the paper);
  the group of each accepted species (Data/col_denominators.csv), so the
  tips outside the group are listed by name (D7 removes them);
  the root-to-tip coefficient of variation and root height (CV <= 0.05 is the
  paper's ultrametric rule; whether the time axis is absolute is read from
  the source, not measured);
  the number of trees in the file;
  and the deposit's own metadata (title, creators, description, related
  identifiers), which is where the source study is named.

What it cannot measure, and the dossier leaves to the curator with the
evidence: which tips were placed with molecular data, whether the tree is
dated in absolute time and how, what date uncertainty the archive carries,
and whether another file of the deposit is the one to use.

USAGE
  python3 living/curate.py --group Scorpiones --source zenodo --record 20973872 \\
      --file living/cache/downloads/<file> [--member path/inside/archive]
"""
from __future__ import annotations

import argparse
import io
import json
import os
import re
import sys
import tarfile
import urllib.parse
import zipfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import watch as W  # noqa: E402

OUT = os.path.join(HERE, "curation")


def read_member(path, member):
    blob = open(path, "rb").read()
    if not member:
        return blob
    if path.endswith(".zip"):
        return zipfile.ZipFile(io.BytesIO(blob)).read(member)
    t = tarfile.open(fileobj=io.BytesIO(blob))
    return t.extractfile(member).read()


def match_species(key):
    """(accepted species name or '', status, family, classification names)."""
    q = key.replace("_", " ")
    url = ("https://api.checklistbank.org/dataset/%d/match/nameusage?" %
           W.COL_DATASET + urllib.parse.urlencode({"q": q, "rank": "species"}))
    try:
        d = W.http(url)
    except Exception as e:
        return "", "error: %s" % e, "", []
    u = d.get("usage") or {}
    cls = u.get("classification", []) or []
    names = [u.get("name", "")] + [c.get("name", "") for c in cls]
    fam = next((c["name"] for c in cls if c.get("rank") == "family"), "")
    if u.get("rank") == "species" and u.get("status") in (
            "accepted", "provisionally accepted"):
        return u["name"], "accepted", fam, names
    if u.get("status") in ("synonym", "ambiguous synonym"):
        sp = next((c for c in cls if c.get("rank") == "species"), None)
        if sp:
            return sp["name"], "synonym", fam, names
    return "", d.get("type", "none"), fam, names


def record_meta(source, record):
    if source == "zenodo":
        d = W.http("https://zenodo.org/api/records/%s" % record)
        md = d.get("metadata", {})
        return {"title": md.get("title", ""),
                "creators": [c.get("name", "") for c in md.get("creators", [])],
                "doi": d.get("doi", ""),
                "publication_date": md.get("publication_date", ""),
                "description": re.sub(r"<[^>]+>", " ", md.get("description", ""))[:3000],
                "related": md.get("related_identifiers", []),
                "files": [(f.get("key"), f.get("size")) for f in d.get("files", [])]}
    if source == "figshare":
        d = W.http("https://api.figshare.com/v2/articles/%s" % record)
        return {"title": d.get("title", ""),
                "creators": [a.get("full_name", "") for a in d.get("authors", [])],
                "doi": d.get("doi", ""), "publication_date": d.get("published_date", ""),
                "description": re.sub(r"<[^>]+>", " ", d.get("description", ""))[:3000],
                "related": [d.get("resource_doi", "")],
                "files": [(f.get("name"), f.get("size")) for f in d.get("files", [])]}
    if source == "dryad":
        ident = urllib.parse.quote(record, safe="")
        d = W.http("https://datadryad.org/api/v2/datasets/%s" % ident)
        return {"title": d.get("title", ""),
                "creators": ["%s %s" % (a.get("firstName", ""), a.get("lastName", ""))
                             for a in d.get("authors", [])],
                "doi": d.get("identifier", ""),
                "publication_date": d.get("publicationDate", ""),
                "description": re.sub(r"<[^>]+>", " ", (d.get("abstract") or "") + " " +
                                      (d.get("methods") or ""))[:3000],
                "related": d.get("relatedWorks", []), "files": []}
    return {}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--record", required=True)
    ap.add_argument("--file", required=True)
    ap.add_argument("--member", default="")
    a = ap.parse_args(argv)

    raw = read_member(a.file, a.member)
    got, n_trees = W.tree_texts(raw)
    nwk, trans = got
    root = W.parse(nwk)
    labels = []
    for leaf in W.leaves(root):
        lab = leaf.label
        labels.append(trans.get(lab, lab).strip("'\"") if trans else lab)
    keyed = {}
    for i, lab in enumerate(labels):
        k, ok = W.species_key(lab.replace("'", ""), str(i))
        if ok:
            keyed.setdefault(k, []).append(lab)
    m = W.measure_tree(raw)

    keys = sorted(keyed)
    with ThreadPoolExecutor(8) as pool:
        res = dict(zip(keys, pool.map(match_species, keys)))
    groups = W.load_groups()
    accepted, unresolved, off = {}, [], Counter()
    off_names = {}
    for k in keys:
        name, status, fam, names = res[k]
        if not name:
            unresolved.append((k, status))
            continue
        g = W.group_of(names, groups) or "(no group)"
        accepted.setdefault(name, (g, status, fam, k))
    in_group = {n: v for n, v in accepted.items() if v[0] == a.group}
    for n, v in accepted.items():
        if v[0] != a.group:
            off[v[0]] += 1
            off_names.setdefault(v[0], []).append(n)
    basis = W.load_basis().get(a.group, {})
    meta = record_meta(a.source, a.record)
    out = {
        "group": a.group, "source": a.source, "record": a.record,
        "file": os.path.basename(a.file), "member": a.member,
        "record_meta": meta,
        "tips": len(labels), "tips_with_species_key": sum(len(v) for v in keyed.values()),
        "species_keys": len(keys), "n_trees_in_file": n_trees,
        "root_to_tip_cv": m["cv"] if m else "", "root_height": m["root_height"] if m else "",
        "ultrametric": (m["cv"] != "" and m["cv"] <= W.CHRONOGRAM_TOLERANCE) if m else None,
        "accepted_species": len(accepted),
        "accepted_in_group": len(in_group),
        "via_synonym_in_group": sum(1 for v in in_group.values() if v[1] == "synonym"),
        "unresolved_keys": unresolved,
        "off_group": {g: sorted(v) for g, v in off_names.items()},
        "reference_tree": basis.get("canonical_tree", ""),
        "reference_accepted": int(float(basis.get("c0_units") or 0)),
        "reference_is_chronogram": basis.get("is_chronogram", ""),
        "reference_uncertainty": basis.get("uncertainty_available", ""),
        "families_in_group": Counter(v[2] for v in in_group.values()).most_common(),
    }
    d = os.path.join(OUT, a.group.replace(" ", "_"))
    os.makedirs(d, exist_ok=True)
    stem = "%s_%s_%s" % (a.source, re.sub(r"[^A-Za-z0-9.]", "_", a.record),
                         re.sub(r"[^A-Za-z0-9.]", "_",
                                os.path.basename(a.member or a.file))[-60:])
    json.dump(out, open(os.path.join(d, stem + ".json"), "w"), indent=1)
    L = ["# %s: %s" % (a.group, meta.get("title", "")[:150]), "",
         "Deposit: %s %s (doi %s, %s). Creators: %s." % (
             a.source, a.record, meta.get("doi", ""), meta.get("publication_date", ""),
             ", ".join(meta.get("creators", [])[:8])),
         "File: `%s`%s." % (os.path.basename(a.file), " member `%s`" % a.member if a.member else ""),
         "", "## Measured", "",
         "| Quantity | Value |", "|---|---|",
         "| Tips | %d |" % out["tips"],
         "| Tips with a species key | %d |" % out["tips_with_species_key"],
         "| Species keys | %d |" % out["species_keys"],
         "| Accepted species (COL26.9) | %d |" % out["accepted_species"],
         "| Accepted species in %s | %d (%d through a synonym) |" % (
             a.group, out["accepted_in_group"], out["via_synonym_in_group"]),
         "| Unresolved species keys | %d |" % len(unresolved),
         "| Trees in the file | %d |" % n_trees,
         "| Root-to-tip CV | %s |" % (round(out["root_to_tip_cv"], 5) if out["root_to_tip_cv"] != "" else "no branch lengths"),
         "| Root height | %s |" % (round(out["root_height"], 3) if out["root_height"] != "" else ""),
         "| Reference tree | %s: %d accepted species, chronogram %s, date uncertainty %s |" % (
             out["reference_tree"] or "none", out["reference_accepted"],
             out["reference_is_chronogram"] or "-", out["reference_uncertainty"] or "-"),
         ""]
    if off_names:
        L += ["Tips outside %s (removed by D7):" % a.group, ""]
        for g, v in sorted(off_names.items()):
            L.append("- %s (%d): %s" % (g, len(v), ", ".join(sorted(v)[:12]) +
                                        (" ..." if len(v) > 12 else "")))
        L.append("")
    if unresolved:
        L += ["Unresolved keys (%d): %s" % (len(unresolved), ", ".join(
            k for k, _ in unresolved[:30]) + (" ..." if len(unresolved) > 30 else "")), ""]
    L += ["Families: " + ", ".join("%s %d" % fc for fc in out["families_in_group"][:20]), "",
          "## From the deposit", "", meta.get("description", "")[:1500], "",
          "Related: " + "; ".join(str(r.get("identifier", r) if isinstance(r, dict) else r)
                                  for r in (meta.get("related") or [])[:6]), "",
          "## For the curator", "",
          "- Molecular tips: which tips the source placed with sequence data.",
          "- Time scale: absolute (calibrated, in Ma) or relative; how it was dated.",
          "- Date uncertainty archived: tree set, intervals, or none.",
          "- The file to use: is another file of the deposit the reference candidate?",
          "- Decision, with its rule and evidence."]
    open(os.path.join(d, stem + ".md"), "w").write("\n".join(L) + "\n")
    print(os.path.join(d, stem + ".md"))
    print("tips %d, accepted %d, in %s %d (ref %d), unresolved %d, off-group %s, cv %s" % (
        out["tips"], out["accepted_species"], a.group, out["accepted_in_group"],
        out["reference_accepted"], len(unresolved),
        {g: len(v) for g, v in off_names.items()}, out["root_to_tip_cv"]))


if __name__ == "__main__":
    main()
