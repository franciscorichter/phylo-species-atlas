#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""summary_tree_node_ages.py -- where the node ages of a date-assignment
summary tree come from (F4-27; decision D60).

The worked example tests the restriction only on a tree whose model
assumptions hold for it (D59, D60): a binary tree whose node ages are
estimated by a dating analysis of the data.  A synthesis that dates its nodes
by assigning published ages to the nodes those studies share and
interpolating the others does not estimate its node ages.  For each row of

    Data/summary_tree_node_ages.csv   one row per eligible group of the
        worked example.  Curated columns:
        partition    the group
        node_ages    dating_analysis (the summary tree's node ages are
                     estimated by the source's dating analysis) or
                     date_assignment (published ages assigned to nodes,
                     the other nodes interpolated or dated at taxon addition)
        summary      for date_assignment: the summary tree (NEXUS, under Data/)
        ages         for date_assignment: the source's JSON of published ages
                     per node label ('node_ages' key)
        placed       for date_assignment: the source's tree of the tips an
                     input phylogeny places, before taxon addition (Newick)
        evidence     why node_ages is what it is

this script measures, for the date_assignment rows, and writes back:

        n_tips, n_internal               the summary tree
        n_internal_published             its internal nodes whose label has
                                         at least one published age
        n_tips_pruned, n_internal_pruned the summary tree pruned to the tips
                                         of `placed` (the sequenced tips)
        n_internal_pruned_published      nodes of the pruned tree whose
                                         corresponding summary-tree node has
                                         a published age
        n_tips_placed, n_internal_placed the tree of placed tips
        n_polytomy_events_placed         branching events at its nodes with
                                         more than two children (a node with
                                         k children is k - 1 events)

The dating_analysis rows are left blank in these columns: their node ages are
the output of the analysis, which the worked example fits.
"""
import csv
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REGISTER = os.path.join(ROOT, "Data", "summary_tree_node_ages.csv")
CURATED = ["partition", "node_ages", "summary", "ages", "placed", "evidence"]
MEASURED = ["n_tips", "n_internal", "n_internal_published", "n_tips_pruned",
            "n_internal_pruned", "n_internal_pruned_published",
            "n_tips_placed", "n_internal_placed", "n_polytomy_events_placed"]
KINDS = ("dating_analysis", "date_assignment")


def parse_newick(s):
    """Nodes as [label, children]; comments in [] removed, branch lengths
    dropped.  Iterative, so a 10,000-tip caterpillar does not recurse."""
    s = re.sub(r"\[[^\]]*\]", "", s).strip().rstrip(";")
    stack, root, i, n = [], None, 0, len(s)
    cur = None
    while i < n:
        c = s[i]
        if c == "(":
            node = ["", []]
            if stack:
                stack[-1][1].append(node)
            stack.append(node)
            i += 1
        elif c == ",":
            i += 1
        elif c == ")":
            cur = stack.pop()
            i += 1
            j = i
            while j < n and s[j] not in ",();":
                j += 1
            cur[0] = s[i:j].split(":")[0].strip().strip("'")
            i = j
            if not stack:
                root = cur
        else:
            j = i
            while j < n and s[j] not in ",();":
                j += 1
            leaf = [s[i:j].split(":")[0].strip().strip("'"), None]
            stack[-1][1].append(leaf)
            i = j
    return root


def nexus_tree(path):
    t = open(path, encoding="utf-8").read()
    k = t.upper().index("TREE ", t.upper().index("BEGIN TREES"))
    t = t[k:]
    return t[t.index("("):t.index(";", t.index("(")) + 1]


def walk(root):
    out, st = [], [root]
    while st:
        x = st.pop()
        out.append(x)
        if x[1]:
            st.extend(x[1])
    return out


def tips(root):
    return [x[0] for x in walk(root) if x[1] is None]


def internal(root):
    return [x for x in walk(root) if x[1] is not None]


def prune(root, keep):
    """The tree on the tips in `keep`, unary nodes suppressed; each kept
    internal node carries the label of the summary-tree node it is."""
    post, st = [], [(root, False)]
    while st:
        x, done = st.pop()
        if x[1] is None or done:
            post.append(x)
        else:
            st.append((x, True))
            st.extend((c, False) for c in x[1])
    new = {}
    for x in post:
        if x[1] is None:
            new[id(x)] = [x[0], None] if x[0] in keep else None
            continue
        kids = [new[id(c)] for c in x[1] if new[id(c)] is not None]
        new[id(x)] = (None if not kids else kids[0] if len(kids) == 1
                      else [x[0], kids])
    return new[id(root)]


def measure(row):
    summ = parse_newick(nexus_tree(os.path.join(ROOT, "Data", row["summary"])))
    ages = json.load(open(os.path.join(ROOT, "Data", row["ages"]),
                          encoding="utf-8"))["node_ages"]
    placed = parse_newick(open(os.path.join(ROOT, "Data", row["placed"]),
                               encoding="utf-8").read())
    keep = set(tips(placed))
    if len(keep) != len(tips(placed)):
        sys.exit("%s: duplicate tip labels" % row["placed"])
    missing = keep - set(tips(summ))
    if missing:
        sys.exit("%s: %d placed tips are not in the summary tree"
                 % (row["partition"], len(missing)))
    pr = prune(summ, keep)
    si, pi, qi = internal(summ), internal(pr), internal(placed)
    return {
        "n_tips": len(tips(summ)), "n_internal": len(si),
        "n_internal_published": sum(1 for x in si if x[0] in ages),
        "n_tips_pruned": len(tips(pr)), "n_internal_pruned": len(pi),
        "n_internal_pruned_published": sum(1 for x in pi if x[0] in ages),
        "n_tips_placed": len(keep), "n_internal_placed": len(qi),
        "n_polytomy_events_placed": sum(len(x[1]) - 1 for x in qi
                                        if len(x[1]) > 2),
    }


def main():
    with open(REGISTER, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        if r["node_ages"] not in KINDS:
            sys.exit("%s: unknown node_ages %r" % (r["partition"], r["node_ages"]))
        m = measure(r) if r["node_ages"] == "date_assignment" else {}
        for k in MEASURED:
            r[k] = m.get(k, "")
        if m:
            print("%s: %s" % (r["partition"], m))
    with open(REGISTER, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, CURATED + MEASURED, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
