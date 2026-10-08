#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""gbif_unmatched.py -- the tips that the Catalogue of Life base release does
not match, queried against the GBIF backbone (audit A10, a C0 sensitivity).

C0 counts tips that resolve to a COL26.9 base-release accepted name (D18).
Valid species absent from that release fall into 'unmatched' and do not count,
while each partition's denominator comes from its own authority.  This script
asks the GBIF backbone taxonomy about every 'unmatched' query of
Data/standardized/col_resolution.csv, so that the rebuild can report C0 with
the GBIF-resolved tips added, as a sensitivity.  It never changes C0.

RULE.  One request per distinct (query, kingdom):
    GET https://api.gbif.org/v1/species/match?name=<Genus epithet>
        &kingdom=<partition kingdom>&strict=true
The partition kingdom is the majority col_kingdom of its resolved rows.  A
query counts as resolved when matchType is EXACT or FUZZY, rank is SPECIES and
status is ACCEPTED, SYNONYM, HETEROTYPIC_SYNONYM, HOMOTYPIC_SYNONYM,
PROPARTE_SYNONYM or DOUBTFUL; the counting unit is the accepted species
(speciesKey and its canonical name `species`).

OUTPUT  Data/external/gbif_unmatched.csv, one row per (query, kingdom),
        with the date of the request.  Rows already present are reused, so the
        file is the cache; --refresh re-queries everything.

USAGE   python scripts/gbif_unmatched.py [--jobs N] [--refresh]
Needs the network and `requests`.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

import requests

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
COL_RESOLUTION = os.path.join(ROOT, "Data", "standardized", "col_resolution.csv")
OUT = os.path.join(ROOT, "Data", "external", "gbif_unmatched.csv")
API = "https://api.gbif.org/v1/species/match"
HEADER = ["query", "kingdom", "matchType", "rank", "status", "confidence",
          "usageKey", "speciesKey", "species", "gbif_kingdom", "family",
          "queried_on"]
RESOLVED_MATCH = ("EXACT", "FUZZY")
RESOLVED_STATUS = ("ACCEPTED", "SYNONYM", "HETEROTYPIC_SYNONYM",
                   "HOMOTYPIC_SYNONYM", "PROPARTE_SYNONYM", "DOUBTFUL")

csv.field_size_limit(1 << 30)


def kingdoms(rows):
    by = defaultdict(Counter)
    for r in rows:
        if r["col_kingdom"]:
            by[r["partition"]][r["col_kingdom"]] += 1
    return {p: c.most_common(1)[0][0] for p, c in by.items()}


def ask(args):
    query, kingdom = args
    params = {"name": query.replace("_", " "), "kingdom": kingdom,
              "strict": "true"}
    for attempt in range(6):
        try:
            r = requests.get(API, params=params, timeout=60)
            if r.status_code == 200:
                j = r.json()
                return {
                    "query": query, "kingdom": kingdom,
                    "matchType": j.get("matchType", ""),
                    "rank": j.get("rank", ""), "status": j.get("status", ""),
                    "confidence": j.get("confidence", ""),
                    "usageKey": j.get("usageKey", ""),
                    "speciesKey": j.get("speciesKey", ""),
                    "species": j.get("species", ""),
                    "gbif_kingdom": j.get("kingdom", ""),
                    "family": j.get("family", ""),
                    "queried_on": datetime.date.today().isoformat()}
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    return None


def resolved(row):
    return (row["matchType"] in RESOLVED_MATCH and row["rank"] == "SPECIES"
            and row["status"] in RESOLVED_STATUS and row["speciesKey"] != "")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--refresh", action="store_true")
    a = ap.parse_args(argv)
    with open(COL_RESOLUTION, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    king = kingdoms(rows)
    want = sorted({(r["query"], king.get(r["partition"], ""))
                   for r in rows if r["match_type"] == "unmatched"
                   and r["query"]})
    have = {}
    if os.path.exists(OUT) and not a.refresh:
        with open(OUT, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                have[(r["query"], r["kingdom"])] = r
    todo = [k for k in want if k not in have]
    print("%d distinct unmatched queries; %d cached; %d to ask"
          % (len(want), len(want) - len(todo), len(todo)))
    failed = 0
    with ThreadPoolExecutor(a.jobs) as ex:
        for i, got in enumerate(ex.map(ask, todo), 1):
            if got is None:
                failed += 1
                continue
            have[(got["query"], got["kingdom"])] = got
            if i % 2000 == 0:
                print("  %d / %d" % (i, len(todo)))
    out = [have[k] for k in want if k in have]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, HEADER, lineterminator="\n",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    n_res = sum(1 for r in out if resolved(r))
    print("wrote %s (%d rows, %d resolved to a GBIF species; %d requests "
          "failed)" % (os.path.relpath(OUT, ROOT), len(out), n_res, failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
