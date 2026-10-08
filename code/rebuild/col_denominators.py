#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""col_denominators.py -- single-authority denominators from the Catalogue of
Life, for the denominator sensitivity of finding S08 (decisions D9, D18).

For every partition: the number of accepted, extant species-rank names under
the taxa that define it (DEF below, from the partition scopes in
Data/site/data/partitions/*/info.yaml), in the Catalogue of Life base release
(ChecklistBank alias 3LR; COL26.9 when written), read through the
ChecklistBank API.

    count(taxon) = total(status=accepted, rank=species, TAXON_ID=taxon)
                 - total(same, extinct=true)

A partition's count is the sum over its plus-terms minus the sum over its
minus-terms.  A plus-term whose classification already contains another
plus-term of the same partition is skipped (the releases place Chilopoda
differently: under Arthropoda in the base release, under Myriapoda in XR).

COVERAGE HOLES.  The base release covers some classes thinly or not at all
(algae and protists; Thecostraca, Branchiura, Pentastomida, Diplura).  The
rule, applied to every partition: within each plus-term, every class is
compared between the base release and the extended release XR (alias 3LXR) of
the same Catalogue of Life version; where base < HOLE_FRACTION x XR, the class
is a hole and its XR count replaces its base count.  A plus-term that is
itself a class or lower rank is one unit; a plus-term absent from the base
release is a hole in full.  Classes that contain a minus-term are never
filled.  The filled count is the sensitivity denominator; the base count and
the filled classes are reported beside it.  A partition whose filled count is
0 has no Catalogue of Life denominator; the sensitivity keeps its primary one
and says so.

COMPLETENESS.  The animal partitions are summed and compared with the Animalia
total of the same release; the residual is written to the output (a partition
scope that leaves a phylum out shows up there).

OUTPUT: Data/col_denominators.csv, one row per partition.

USAGE
    python scripts/col_denominators.py            # needs network access
"""

from __future__ import annotations

import csv
import datetime
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
BASIS = os.path.join(ROOT, "Data", "coverage_basis.csv")
OUT = os.path.join(ROOT, "Data", "col_denominators.csv")
API = "https://api.checklistbank.org"
BASE, XR = "3LR", "3LXR"
HOLE_FRACTION = 0.5          # rule parameter: base < 0.5 x XR is a hole
S = requests.Session()


def get(path, params=None):
    err = None
    for i in range(6):
        try:
            r = S.get(API + path, params=params, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as e:                              # noqa: BLE001
            err = e
            time.sleep(3 * (i + 1))
    raise err


def T(n, r, a=None, opt=False):
    return dict(name=n, rank=r, anc=a, opt=opt)


A, C, P, PR = "Animalia", "Chromista", "Plantae", "Protozoa"
DEF = {
 "Acanthocephala": ([T("Acanthocephala", "phylum", A)], []),
 "Acari": ([T("Parasitiformes", "superorder", "Arachnida"),
            T("Acariformes", "superorder", "Arachnida")], []),
 "Amoebozoa": ([T("Amoebozoa", "phylum", PR)], []),
 "Amphibians": ([T("Amphibia", "class", "Chordata")], []),
 "Annelids": ([T("Annelida", "phylum", A)],
              [T("Echiura", "subclass", "Annelida", True)]),
 "Apicomplexa": ([T("Miozoa", "phylum", C), T("Apicomplexa", "subphylum", C, True)],
                 [T("Dinophyceae", "class", "Miozoa")]),
 "Archaea": ([T("Archaea", "domain")], []),
 "Bacteria": ([T("Bacteria", "domain")], []),
 "Birds": ([T("Aves", "class", "Chordata")], []),
 "Bivalves": ([T("Bivalvia", "class", "Mollusca")], []),
 "Brown algae": ([T("Phaeophyceae", "class", C, True)], []),
 "Bryophytes": ([T("Bryophyta", "phylum", P), T("Marchantiophyta", "phylum", P),
                 T("Anthocerotophyta", "phylum", P)], []),
 "Bryozoa": ([T("Bryozoa", "phylum", A)], []),
 "Cephalopods": ([T("Cephalopoda", "class", "Mollusca")], []),
 "Cercozoa": ([T("Cercozoa", "phylum", C)], []),
 "Choanoflagellates": ([T("Choanoflagellatea", "class", PR)], []),
 "Chrysophyta": ([T("Chrysophyceae", "class", C, True),
                  T("Synurophyceae", "class", C, True)], []),
 "Ciliates": ([T("Ciliophora", "phylum", C)], []),
 "Cnidaria": ([T("Cnidaria", "phylum", A)],
              [T("Myxozoa", "class", "Cnidaria", True)]),
 "Crocodilians": ([T("Crocodylia", "order", "Reptilia")], []),
 "Crustaceans": ([T("Crustacea", "subphylum", "Arthropoda")], []),
 "Ctenophores": ([T("Ctenophora", "phylum", A)], []),
 "Diatoms": ([T("Bacillariophyta", "phylum", C, True),
              T("Bacillariophyceae", "class", C, True),
              T("Coscinodiscophyceae", "class", C, True),
              T("Mediophyceae", "class", C, True)], []),
 "Dinoflagellates": ([T("Dinophyceae", "class", C),
                      T("Syndiniophyceae", "class", C, True),
                      T("Noctilucophyceae", "class", C, True)], []),
 "Echinoderms": ([T("Echinodermata", "phylum", A)], []),
 "Euglenozoa": ([T("Euglenozoa", "phylum", PR)], []),
 "Ferns": ([T("Polypodiopsida", "class", "Tracheophyta")], []),
 "Fish": ([T("Actinopterygii", "gigaclass", "Chordata")], []),
 "Flatworms": ([T("Platyhelminthes", "phylum", A)], []),
 "Foraminifera": ([T("Foraminifera", "phylum", C)], []),
 "Fungi": ([T("Fungi", "kingdom", "Eukaryota")], []),
 "Gastropods": ([T("Gastropoda", "class", "Mollusca")], []),
 "Green algae": ([T("Chlorophyta", "phylum", P), T("Charophyta", "phylum", P)], []),
 "Hagfish": ([T("Myxini", "class", "Chordata")], []),
 "Haptophytes": ([T("Haptophyta", "phylum", C, True)], []),
 "Insects": ([T("Insecta", "class", "Arthropoda")], []),
 "Lampreys": ([T("Petromyzonti", "class", "Chordata")], []),
 "Lycophytes": ([T("Lycopodiopsida", "class", "Tracheophyta")], []),
 "Mammals": ([T("Mammalia", "class", "Chordata")], []),
 "Minor animal phyla": ([T(n, "phylum", A) for n in
                         ("Onychophora", "Placozoa", "Loricifera",
                          "Cycliophora", "Chaetognatha")], []),
 # Chilopoda: a class under Arthropoda in the base release, inside Myriapoda
 # in XR (skipped there as nested).
 "Myriapods": ([T("Myriapoda", "subphylum", "Arthropoda"),
                T("Chilopoda", "class", "Arthropoda")], []),
 "Myxozoa": ([T("Myxozoa", "class", A)], []),
 "Nematodes": ([T("Nematoda", "phylum", A)], []),
 "Nemerteans": ([T("Nemertea", "phylum", A)], []),
 "Oomycetes": ([T("Oomycota", "phylum", C)], []),
 "Other Mollusca": ([T("Mollusca", "phylum", A)],
                    [T("Cephalopoda", "class", "Mollusca"),
                     T("Bivalvia", "class", "Mollusca"),
                     T("Gastropoda", "class", "Mollusca")]),
 "Other arachnids": ([T("Arachnida", "class", "Chelicerata"),
                      T("Xiphosura", "order", "Chelicerata")],
                     [T("Araneae", "order", "Arachnida"),
                      T("Scorpiones", "order", "Arachnida"),
                      T("Parasitiformes", "superorder", "Arachnida"),
                      T("Acariformes", "superorder", "Arachnida")]),
 "Other hexapods": ([T("Collembola", "class", "Arthropoda"),
                     T("Protura", "class", "Arthropoda"),
                     T("Diplura", "class", "Arthropoda")], []),
 # The info.yaml scope lists ten phyla 'primarily'; the animal phyla that no
 # other partition takes (measured by the completeness residual) join here.
 "Other phyla": ([T(n, "phylum", A) for n in
                  ("Gastrotricha", "Kinorhyncha", "Priapulida", "Sipuncula",
                   "Phoronida", "Hemichordata", "Brachiopoda", "Entoprocta",
                   "Gnathostomulida", "Nematomorpha", "Xenacoelomorpha",
                   "Dicyemida", "Orthonectida")]
                 + [T("Cephalochordata", "subphylum", "Chordata"),
                    T("Echiura", "subclass", "Annelida", True)], []),
 "Other vertebrates": ([T("Actinistia", "superclass", "Chordata"),
                        T("Dipnoi", "superclass", "Chordata"),
                        T("Sphenodontida", "suborder", "Reptilia")], []),
 "Pycnogonida": ([T("Pycnogonida", "class", "Chelicerata")], []),
 "Rhodophyta": ([T("Rhodophyta", "phylum", P)], []),
 "Rotifera": ([T("Rotifera", "phylum", A)], []),
 "Scorpiones": ([T("Scorpiones", "order", "Arachnida")], []),
 "Seed plants": ([T("Tracheophyta", "phylum", P)],
                 [T("Polypodiopsida", "class", "Tracheophyta"),
                  T("Lycopodiopsida", "class", "Tracheophyta")]),
 "Sharks": ([T("Chondrichthyes", "parvphylum", "Chordata")], []),
 "Spiders": ([T("Araneae", "order", "Arachnida")], []),
 "Sponges": ([T("Porifera", "phylum", A)], []),
 "Squamates": ([T("Squamata", "order", "Reptilia")], []),
 "Tardigrada": ([T("Tardigrada", "phylum", A)], []),
 "Tunicates": ([T("Tunicata", "subphylum", "Chordata")], []),
 "Turtles": ([T("Testudines", "order", "Reptilia")], []),
}
RANK_ORDER = ["domain", "kingdom", "phylum", "subphylum", "parvphylum",
              "gigaclass", "superclass", "class", "subclass", "superorder",
              "order", "suborder"]


def resolve(ds, t):
    j = get("/dataset/%s/nameusage/search" % ds,
            dict(q=t["name"], content="SCIENTIFIC_NAME", type="EXACT",
                 limit=50))
    hits = []
    for r in j.get("result", []):
        u = r["usage"]
        if u["name"]["scientificName"] != t["name"]:
            continue
        path = [c["name"] for c in r["classification"][:-1]]
        if u["status"] == "accepted" and u["name"].get("rank") == t["rank"] \
                and (t["anc"] is None or t["anc"] in path):
            hits.append(dict(id=r["id"], path=path))
    if len(hits) > 1:
        raise RuntimeError("ambiguous %s in %s: %s" % (t["name"], ds, hits))
    return hits[0] if hits else None


def n(ds, tid, **kw):
    return get("/dataset/%s/nameusage/search" % ds,
               dict(TAXON_ID=tid, rank="species", limit=0, **kw))["total"]


def extant(ds, tid):
    return n(ds, tid, status="accepted") - n(ds, tid, status="accepted",
                                             extinct="true")


def classes(ds, tid):
    """{class name: id} of accepted classes under a taxon."""
    out, off = {}, 0
    while True:
        j = get("/dataset/%s/nameusage/search" % ds,
                dict(TAXON_ID=tid, rank="class", status="accepted",
                     limit=500, offset=off))
        for r in j.get("result", []):
            out[r["usage"]["name"]["scientificName"]] = r["id"]
        off += 500
        if off >= j.get("total", 0):
            return out


def is_class_or_lower(rank):
    return rank in RANK_ORDER and RANK_ORDER.index(rank) >= \
        RANK_ORDER.index("class")


def units_of(ds, hp, hm):
    """{unit: extant count} over a release's plus-terms.  A plus-term of class
    rank or lower is one unit; a higher one contributes its classes plus a
    '<term> (no class)' remainder.  Classes that are or contain a minus-term
    are left out of the comparison (they are never filled)."""
    minus_names = {t["name"] for t, h in hm if h}
    minus_paths = [set(h["path"]) for _t, h in hm if h]
    units = {}
    for t, h in hp:
        if h is None:
            continue
        total = extant(ds, h["id"])
        if is_class_or_lower(t["rank"]):
            units[t["name"]] = units.get(t["name"], 0) + total
            continue
        rest = total
        for cname, cid in sorted(classes(ds, h["id"]).items()):
            c = extant(ds, cid)
            rest -= c
            if cname in minus_names or any(cname in mp for mp in minus_paths):
                continue
            units[cname] = units.get(cname, 0) + c
        units["%s (no class)" % t["name"]] = rest
    return units


def evaluate(part, keys):
    plus, minus = DEF[part]
    res, counts = {}, {}
    for ds in (BASE, XR):
        ks = keys[ds]
        hp = [(t, resolve(ks, t)) for t in plus]
        names = {t["name"] for t, h in hp if h}
        hp = [(t, h) for t, h in hp
              if h is None or not (names & set(h["path"]))]
        hm = [(t, resolve(ks, t)) for t in minus]
        res[ds] = (hp, hm)
        counts[ds] = (sum(extant(ks, h["id"]) for _t, h in hp if h)
                      - sum(extant(ks, h["id"]) for _t, h in hm if h))
    ub = units_of(keys[BASE], *res[BASE])
    ux = units_of(keys[XR], *res[XR])
    fills = [(u, ub.get(u, 0), x) for u, x in sorted(ux.items())
             if x > 0 and ub.get(u, 0) < HOLE_FRACTION * x]
    filled = counts[BASE] + sum(x - b for _u, b, x in fills)
    return dict(
        partition=part,
        col_taxa=" + ".join(t["name"] for t in plus)
        + "".join(" - " + t["name"] for t in minus),
        col_base=counts[BASE], col_xr=counts[XR], col_filled=filled,
        filled_classes="; ".join("%s base %d < %.2f x XR %d"
                                 % (u, b, HOLE_FRACTION, x)
                                 for u, b, x in fills),
        absent_in_base="; ".join(t["name"] for t, h in res[BASE][0]
                                 if h is None))


def main():
    with open(BASIS, newline="", encoding="utf-8") as fh:
        parts = [r["partition"] for r in csv.DictReader(fh)]
    if set(parts) != set(DEF):
        sys.exit("partitions differ from DEF: %s" % (set(parts) ^ set(DEF)))
    meta = {ds: get("/dataset/%s" % ds) for ds in (BASE, XR)}
    keys = {ds: meta[ds]["key"] for ds in (BASE, XR)}
    with ThreadPoolExecutor(6) as ex:
        rows = list(ex.map(lambda p: evaluate(p, keys), sorted(parts)))
    # completeness: animal partitions vs Animalia
    anim = [p for p in parts if any(t["anc"] == A or t["anc"] in
            ("Arachnida", "Chordata", "Mollusca", "Arthropoda", "Annelida",
             "Cnidaria", "Chelicerata", "Reptilia")
            for t in DEF[p][0])]
    tot = {}
    for ds in (BASE, XR):
        h = resolve(keys[ds], T("Animalia", "kingdom"))
        tot[ds] = extant(keys[ds], h["id"])
    s_base = sum(r["col_base"] for r in rows if r["partition"] in anim)
    s_xr = sum(r["col_xr"] for r in rows if r["partition"] in anim)
    stamp = datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="seconds")
    rel = lambda ds: "%s (ChecklistBank %s, doi:%s, issued %s)" % (
        meta[ds]["alias"], meta[ds]["key"], meta[ds]["doi"], meta[ds]["issued"])
    for r in rows:
        r.update(release_base=rel(BASE), release_xr=rel(XR), retrieved=stamp,
                 animal_partition="yes" if r["partition"] in anim else "no")
    hdr = ["partition", "col_taxa", "col_base", "col_filled", "col_xr",
           "filled_classes", "absent_in_base", "animal_partition",
           "release_base", "release_xr", "retrieved"]
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, hdr, lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
        fh.write("# completeness: Animalia extant accepted species %s %d, sum "
                 "of animal partitions %d, residual %d; %s %d, %d, %d\n"
                 % (BASE, tot[BASE], s_base, tot[BASE] - s_base, XR, tot[XR],
                    s_xr, tot[XR] - s_xr))
    print("wrote %s; Animalia residual base %d, XR %d"
          % (os.path.relpath(OUT, ROOT), tot[BASE] - s_base, tot[XR] - s_xr))


if __name__ == "__main__":
    main()
