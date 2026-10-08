#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""prepare_canonical_sources.py -- relabel a deposited source tree to atlas ids.

Data/canonical_source_prep.csv names, per canonical file whose source is a
deposited tree other than the current shipped file (decisions D8, D15, D16,
D17 and the birds provenance correction), the raw deposit and how its labels
map to the name layer.  This script reads that deposit (Data/raw, which is not
distributed), and writes the id-labelled, UNPRUNED source to
Data/_archive/2026-10-revision/canonical_sources/<filename>, which is
distributed and from which scripts/canonical_trees.py builds the shipped file.

Columns of the register:
    filename      the canonical file it becomes
    raw_source    path under Data/ of the deposited tree file
    tree_index    1-based index of the tree in that file
    branch_scale  multiplier applied to every branch length (decimal
                  arithmetic, so a length x scale is exact); 1 for none
    crosswalk     optional TSV under Data/ (raw_label <TAB> standardized_name)
                  for labels that differ from the name layer
    group         full_mapping.csv group for labels new to the dictionary
    capitalize    'yes' to capitalise the first letter of a new label's
                  standardized name (deposits written in lower case)
    evidence      provenance, checksums, the decision that applies

A raw label resolves, in order: crosswalk -> full_mapping.csv original_label of
the group -> dictionary.csv standardized_name (spaces and underscores
equivalent).  A label that resolves nowhere is NEW: it gets the next free id
(new names sorted, ids assigned in that order) and rows in dictionary.csv,
full_mapping.csv; scripts/match_col.py then resolves it against the Catalogue
of Life.  Comments ([&...]) are dropped; branch-length
text is kept unless scaled.

USAGE
    python scripts/prepare_canonical_sources.py [--only FILENAME ...] [--check]
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from decimal import Decimal

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
sys.path.insert(0, HERE)
import canonical_trees as CT  # noqa: E402

DATA = os.path.join(ROOT, "Data")
STD = os.path.join(DATA, "standardized")
REGISTER = os.path.join(DATA, "canonical_source_prep.csv")
OUTDIR = os.path.join(DATA, "_archive", "2026-10-revision", "canonical_sources")
DICT = os.path.join(STD, "dictionary.csv")
FULL = os.path.join(STD, "full_mapping.csv")

csv.field_size_limit(1 << 30)


def read_nth_tree(path, n):
    """The n-th tree (1-based) of a Newick or NEXUS file, with a NEXUS
    translate table applied.  Returns the Newick string."""
    translate = {}
    in_tr = False
    count = 0
    buf = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            low = s.lower()
            if low.startswith("translate"):
                in_tr = True
                continue
            if in_tr:
                for item in s.rstrip(";").split(","):
                    item = item.strip()
                    if item:
                        kv = item.split(None, 1)
                        if len(kv) == 2:
                            translate[kv[0]] = kv[1].strip().strip("'")
                if s.endswith(";"):
                    in_tr = False
                continue
            if buf is not None:
                buf += s
                if ";" in s:
                    break
                continue
            if "(" in s and (low.startswith("tree ") or low.startswith("utree ")
                             or s.startswith("(")):
                count += 1
                if count == n:
                    buf = s[s.index("("):] if "=" not in s.split("(")[0] \
                        else s[s.index("(", s.index("=")):]
                    if ";" in buf:
                        break
    if buf is None:
        raise SystemExit("%s has fewer than %d trees" % (path, n))
    buf = buf[:buf.index(";") + 1]
    if translate:
        def sub(m):
            return m.group(1) + translate.get(m.group(2), m.group(2))
        buf = re.sub(r"([(,])\s*([^(),:;\[\s]+)", sub, buf)
    return buf


def read_csv_rows(path):
    with open(path, newline="", encoding="utf-8") as fh:
        r = csv.reader(fh)
        hdr = next(r)
        return hdr, list(r)


def terminator(path):
    with open(path, "rb") as fh:
        return "\r\n" if fh.readline().endswith(b"\r\n") else "\n"


def append_rows(path, rows):
    term = terminator(path)
    with open(path, "a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator=term)
        w.writerows(rows)


def norm(s):
    return s.replace(" ", "_")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    with open(REGISTER, newline="", encoding="utf-8") as fh:
        reg = list(csv.DictReader(fh))
    _h, drows = read_csv_rows(DICT)
    by_name = {norm(r[1]): r[0] for r in drows}
    by_id = {r[0]: r[1] for r in drows}
    max_id = max(int(r[0]) for r in drows)
    fh_, frows = read_csv_rows(FULL)
    os.makedirs(OUTDIR, exist_ok=True)
    for row in reg:
        fn = row["filename"].strip()
        if a.only and fn not in a.only:
            continue
        group = row["group"].strip()
        orig = {}
        for r in frows:
            if r[3] == group:
                orig.setdefault(r[1], r[0])
        cw = {}
        if row.get("crosswalk", "").strip():
            with open(os.path.join(DATA, row["crosswalk"].strip()),
                      encoding="utf-8") as fh:
                for line in fh:
                    if line.strip() and not line.startswith("#"):
                        k, _t, v = line.rstrip("\n").partition("\t")
                        cw[k] = v
        text = read_nth_tree(os.path.join(DATA, row["raw_source"].strip()),
                             int(row["tree_index"] or 1))
        root = CT.parse(text)
        scale = Decimal(row.get("branch_scale") or "1")
        new_names = {}
        add_f = []
        stack = [root]
        while stack:
            nd = stack.pop()
            stack.extend(nd.kids)
            if nd.bl is not None and scale != 1:
                nd.bl = str(Decimal(str(nd.bl)) * scale)
            if nd.kids:
                continue
            lab = nd.label
            std = cw.get(lab)
            sid = None
            if std is None and lab in orig:
                sid = orig[lab]
            if sid is None:
                std = norm(std if std is not None else lab)
                sid = by_name.get(std)
            if sid is None:
                if row.get("capitalize", "").strip() == "yes":
                    std = std[:1].upper() + std[1:]
                sid = by_name.get(std)
            if sid is None:
                new_names.setdefault(std, []).append(nd)
                nd.raw = "\x00" + std
                continue
            if lab not in orig:
                add_f.append([sid, lab, by_id[sid], group])
                orig[lab] = sid
            nd.label = nd.raw = sid
        # internal labels: keep as written (support values, node names)
        add_d = []
        for std in sorted(new_names):
            max_id += 1
            sid = str(max_id)
            by_name[std] = sid
            by_id[sid] = std
            for nd in new_names[std]:
                add_f.append([sid, nd.label, std, group])
                nd.label = nd.raw = sid
            add_d.append([sid, std])
        leaves = CT.leaves(root)
        ids = [l.label for l in leaves]
        print("%-24s %d leaves, %d distinct ids, %d new names"
              % (fn, len(ids), len(set(ids)), len(add_d)))
        if a.check:
            continue
        out = CT.serialise(root)
        with open(os.path.join(OUTDIR, fn), "w", encoding="utf-8") as fh:
            fh.write(out + "\n")
        if add_d:
            append_rows(DICT, add_d)
        if add_f:
            append_rows(FULL, add_f)
            frows.extend(add_f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
