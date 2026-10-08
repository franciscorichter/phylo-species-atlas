#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""canonical_trees.py -- build the shipped canonical trees from their sources.

Decisions D4, D7 and D8 (Paper/revision-2026-10/DECISIONS.md) change what a
canonical tree file contains.  This script applies them by rule, so that no
shipped canonical file is edited by hand:

    Data/canonical_sources.csv   one row per canonical tree: the archived,
                                 id-labelled, UNPRUNED source file it is built
                                 from (Data/_archive/2026-10-revision/
                                 canonical_sources/<filename>) and where that
                                 source came from.
    Data/pruned_tips.csv         D7: tips outside the partition, by
                                 standardized id, per canonical file.

For each row the shipped file Data/standardized/trees/<filename> is

    1. the source tree, with every leaf listed in pruned_tips.csv removed (D7);
    2. then, walking the leaves left to right, every leaf whose SPECIES KEY was
       already seen removed, so one tip per distinct species label remains (D4);
    3. unary nodes left by the removals collapsed, their branch lengths summed;
       a root left with one child is replaced by that child (its stem branch is
       dropped), so the root of the shipped file is the MRCA of the retained
       tips.

Every branch length that is not a sum is copied as its original text.

SPECIES KEY (D4, D13).  From the standardized name, the first pair of adjacent
tokens (split on '_', '.' and space) of the form Genus epithet, Genus matching
[A-Z][a-z]+ and not ending in a suprageneric suffix, epithet matching
[a-z][a-z-]+.  If the epithet is an open-nomenclature or hybrid marker (sp,
spp, cf, aff, af, nr, indet, x, hybrid, ...), or such a marker occurs anywhere
in the name, the label is UNRESOLVED and its key is the id itself: it is never
merged with another tip.  Two tips merge only when their keys are equal.

The names map Data/standardized/names/<filename>.json is rewritten to the
retained ids.  Run with --check to report what would change without writing.

USAGE
    python scripts/canonical_trees.py [--check] [--only FILENAME ...]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
DATA = os.path.join(ROOT, "Data")
STD = os.path.join(DATA, "standardized")
SOURCES = os.path.join(DATA, "canonical_sources.csv")
PRUNED = os.path.join(DATA, "pruned_tips.csv")
DICTIONARY = os.path.join(STD, "dictionary.csv")

csv.field_size_limit(1 << 30)
sys.setrecursionlimit(10000)

OPEN_MARKERS = {"sp", "spp", "sp.", "spp.", "cf", "cf.", "aff", "aff.", "af",
                "nr", "nr.", "indet", "indet.", "x", "hybrid", "nov", "n",
                "undescribed", "unidentified", "unknown", "gen", "trans",
                "genome", "sp1", "sp2", "sp3", "sp4", "sp5"}
SUPRAGENERIC = ("idae", "aceae", "inae", "oidea", "ales", "phyta", "mycota",
                "mycetes", "opsida", "ineae")
# Genera of the Catalogue of Life (COL26.9 base release, D18) whose names end
# in one of the SUPRAGENERIC suffixes (Xerophyta, Elapsoidea, Pimephales, ...).
# The suffix test marks a token as a family/order name only when the authority
# does not list it as a genus.  Written by scripts/match_col.py --genera from
# the release index; verify_atlas.py re-derives it when the index is present.
SUPRAGENERIC_GENERA_FILE = os.path.join(ROOT, "Data", "standardized",
                                        "col269_suprageneric_genera.txt")


def _load_suprageneric_genera():
    out = set()
    with open(SUPRAGENERIC_GENERA_FILE, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                out.add(line)
    return frozenset(out)


SUPRAGENERIC_GENERA = _load_suprageneric_genera()
GENUS_RE = re.compile(r"^[A-Z][a-z]+$")


def is_genus_token(g):
    """A capitalised token is a genus unless it ends in a suprageneric suffix
    AND the authority lists no genus of that name."""
    return bool(GENUS_RE.match(g)) and (
        g in SUPRAGENERIC_GENERA or not g.lower().endswith(SUPRAGENERIC))
EPITHET_RE = re.compile(r"^[a-z][a-z-]+$")
SPLIT_RE = re.compile(r"[_.\s]+")


def species_key(name, sid):
    """(key, resolved).  resolved=False means the label is not a species
    label (D13) and the key is the id, so it merges with nothing."""
    toks = [t for t in SPLIT_RE.split(name or "") if t]
    low = {t.lower() for t in toks}
    if low & {"x", "hybrid", "cf", "cf.", "aff", "aff.", "af"}:
        return "id:" + sid, False
    for i in range(len(toks) - 1):
        g, e = toks[i], toks[i + 1]
        if is_genus_token(g) and EPITHET_RE.match(e):
            if e.lower() in OPEN_MARKERS or re.match(r"^sp\d*$", e):
                return "id:" + sid, False
            return g + "_" + e, True
        if GENUS_RE.match(g) and e.lower() in OPEN_MARKERS:
            return "id:" + sid, False
    return "id:" + sid, False


# --------------------------------------------------------------------------
# Newick: parse keeping branch-length text, prune, serialise
# --------------------------------------------------------------------------

TOKEN_RE = re.compile(r"\[[^\]]*\]|'(?:[^']|'')*'|[(),;:]|[^(),;:\[']+")


class Node(object):
    __slots__ = ("kids", "label", "raw", "bl", "parent")

    def __init__(self):
        self.kids = []
        self.label = ""
        self.raw = ""       # the label token as written (quotes kept)
        self.bl = None      # original text of the branch length, or a float
        self.parent = None


def parse(text):
    root = Node()
    cur = root
    expect_bl = False
    for m in TOKEN_RE.finditer(text):
        t = m.group(0)
        if t.startswith("["):
            continue
        if t == "(":
            n = Node()
            n.parent = cur
            cur.kids.append(n)
            cur = n
            expect_bl = False
        elif t == ",":
            n = Node()
            n.parent = cur.parent
            cur.parent.kids.append(n)
            cur = n
            expect_bl = False
        elif t == ")":
            cur = cur.parent
            expect_bl = False
        elif t == ":":
            expect_bl = True
        elif t == ";":
            break
        else:
            s = t.strip()
            if not s:
                continue
            if expect_bl:
                cur.bl = s
                expect_bl = False
            else:
                cur.raw = s
                if s.startswith("'") and s.endswith("'"):
                    s = s[1:-1].replace("''", "'")
                cur.label = s
    # A leading '(' makes `root` an empty wrapper around the real root.
    while not root.raw and root.bl is None and len(root.kids) == 1 \
            and root.kids[0].parent is root:
        root = root.kids[0]
        root.parent = None
    return root


def leaves(root):
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        if n.kids:
            stack.extend(reversed(n.kids))
        else:
            out.append(n)
    return out


def _blf(x):
    if x is None:
        return 0.0
    if isinstance(x, float):
        return x
    try:
        return float(x)
    except ValueError:
        return float("nan")


def remove_leaves(root, drop):
    """Remove the leaf nodes in `drop` (a set of Node objects); collapse unary
    nodes; return the new root."""
    for leaf in drop:
        p = leaf.parent
        p.kids.remove(leaf)
        leaf.parent = None
        while p is not None and not p.kids:
            q = p.parent
            if q is None:
                return None
            q.kids.remove(p)
            p = q
    # collapse unary nodes, bottom-up, iteratively
    order, stack = [], [root]
    while stack:
        n = stack.pop()
        order.append(n)
        stack.extend(n.kids)
    for n in reversed(order):
        if len(n.kids) == 1:
            c = n.kids[0]
            if n.parent is None:
                c.parent = None
                c.bl = None
                root = c
                continue
            if n.bl is None and c.bl is None:
                c.bl = None
            else:
                c.bl = _blf(n.bl) + _blf(c.bl)
            p = n.parent
            i = p.kids.index(n)
            p.kids[i] = c
            c.parent = p
    while len(root.kids) == 1:
        root = root.kids[0]
        root.parent = None
        root.bl = None
    return root


def _fmt_bl(bl):
    if bl is None:
        return ""
    if isinstance(bl, float):
        return ":" + repr(bl)
    return ":" + bl


def serialise(root):
    out = []
    stack = [(root, 0)]
    while stack:
        n, state = stack.pop()
        if n.kids and state < len(n.kids):
            out.append("(" if state == 0 else ",")
            stack.append((n, state + 1))
            stack.append((n.kids[state], 0))
        else:
            if n.kids:
                out.append(")")
            out.append(n.raw)
            out.append(_fmt_bl(n.bl))
    out.append(";")
    return "".join(out)


# --------------------------------------------------------------------------

def read_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_names():
    out = {}
    with open(DICTIONARY, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            out[r["id"]] = r["standardized_name"]
    return out


def build_one(src_text, pruned_ids, names):
    """Return (newick, retained ids in order, stats)."""
    root = parse(src_text)
    lv = leaves(root)
    drop, seen, n_d4, n_d7 = set(), set(), 0, 0
    for leaf in lv:
        sid = leaf.label
        if sid in pruned_ids:
            drop.add(leaf)
            n_d7 += 1
            continue
        key, _resolved = species_key(names.get(sid, ""), sid)
        if key in seen:
            drop.add(leaf)
            n_d4 += 1
            continue
        seen.add(key)
    root = remove_leaves(root, drop)
    kept = [l.label for l in leaves(root)]
    return serialise(root), kept, {"source_leaves": len(lv), "d7_pruned": n_d7,
                                   "d4_collapsed": n_d4, "shipped": len(kept)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args(argv)
    names = load_names()
    pruned = {}
    for r in read_csv(PRUNED):
        pruned.setdefault(r["filename"].strip(), set()).add(
            r["standardized_id"].strip())
    rc = 0
    for r in read_csv(SOURCES):
        fn = r["filename"].strip()
        if a.only and fn not in a.only:
            continue
        src = os.path.join(ROOT, r["source"].strip())
        with open(src, encoding="utf-8") as fh:
            text = fh.read()
        missing = [i for i in pruned.get(fn, set())
                   if not re.search(r"[(,]%s[:,);]" % re.escape(i), text)]
        if missing:
            print("REFUSED %s: pruned ids not in the source: %s"
                  % (fn, missing[:5]))
            rc = 1
            continue
        nwk, kept, st = build_one(text, pruned.get(fn, set()), names)
        unknown = [k for k in kept if k not in names]
        if unknown:
            print("REFUSED %s: %d tip id(s) not in dictionary.csv, e.g. %s"
                  % (fn, len(unknown), unknown[:5]))
            rc = 1
            continue
        dest = os.path.join(STD, "trees", fn)
        old = open(dest, encoding="utf-8").read() if os.path.exists(dest) else ""
        status = "unchanged" if old == nwk + "\n" or old == nwk else "changed"
        print("%-28s %s  source %d  D7 -%d  D4 -%d  shipped %d"
              % (fn, status, st["source_leaves"], st["d7_pruned"],
                 st["d4_collapsed"], st["shipped"]))
        if a.check or status == "unchanged":
            continue
        with open(dest, "w", encoding="utf-8") as fh:
            fh.write(nwk + "\n")
        nm = {k: names[k] for k in kept}
        with open(os.path.join(STD, "names", fn + ".json"), "w",
                  encoding="utf-8") as fh:
            json.dump(nm, fh, ensure_ascii=False, separators=(",", ":"))
    return rc


if __name__ == "__main__":
    sys.exit(main())
