#!/usr/bin/env python3
"""Derive a species-level tree from a standardized multi-individual tree.

Reads a Newick tree where tips are numeric IDs (resolved against
standardized/dictionary.csv to specimen-level names like
`Genus_species_specimen_code`) and prunes it to one representative tip
per species (defined as the first two underscore-separated tokens of
the standardized name).

Branch-length policy: when collapsing siblings of the same species, the
representative inherits its own original branch length. Internal unary
nodes are collapsed by summing branch lengths into the surviving child.

Usage:
    python scripts/derive_species_tree.py \\
        --in standardized/trees/turtles.nwk \\
        --dictionary standardized/dictionary.csv \\
        --out site/data/partitions/turtles/tree.nwk
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path


def parse_newick(s: str) -> dict:
    """Parse Newick into nested dict: {name, length, children?}"""
    s = s.strip().rstrip(";").strip()
    i = [0]
    def node():
        n = {"name": "", "length": None, "children": None}
        if i[0] < len(s) and s[i[0]] == "(":
            i[0] += 1
            n["children"] = []
            while True:
                n["children"].append(node())
                if i[0] < len(s) and s[i[0]] == ",":
                    i[0] += 1
                    continue
                if i[0] < len(s) and s[i[0]] == ")":
                    i[0] += 1
                    break
                raise ValueError(f"newick parse error at {i[0]}: expected , or )")
        name = ""
        while i[0] < len(s) and s[i[0]] not in ",():":
            name += s[i[0]]
            i[0] += 1
        n["name"] = name.strip()
        if i[0] < len(s) and s[i[0]] == ":":
            i[0] += 1
            ln = ""
            while i[0] < len(s) and s[i[0]] not in ",()":
                ln += s[i[0]]
                i[0] += 1
            try:
                n["length"] = float(ln)
            except ValueError:
                n["length"] = None
        return n
    return node()


def write_newick(node: dict) -> str:
    if node.get("children"):
        inside = ",".join(write_newick(c) for c in node["children"])
        body = f"({inside})"
    else:
        body = node["name"]
    if node["length"] is not None:
        body += f":{node['length']:g}"
    elif node.get("children"):
        body += node["name"]
    return body


def prune_keep_refs(node: dict, keep: set[int]) -> dict | None:
    """Recursively prune the tree to retain only leaf NODES whose Python id() is
    in `keep`. Keying by node identity (not name) is required because the
    standardized turtle tree contains duplicate numeric IDs at multiple tip
    positions — the redundancy that prompted this whole derivation.

    Collapses internal nodes with a single remaining child by summing branch
    lengths into the surviving child.
    """
    if not node.get("children"):
        return node if id(node) in keep else None

    kept_children = []
    for c in node["children"]:
        pruned = prune_keep_refs(c, keep)
        if pruned is not None:
            kept_children.append(pruned)

    if not kept_children:
        return None
    if len(kept_children) == 1:
        only = kept_children[0]
        only_len = only["length"] or 0
        node_len = node["length"] or 0
        if node["length"] is not None or only["length"] is not None:
            only["length"] = node_len + only_len
        return only
    node["children"] = kept_children
    return node


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_path", required=True)
    ap.add_argument("--dictionary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", action="store_true", help="Print summary of the derivation")
    args = ap.parse_args()

    # Load dictionary id → name
    name_by_id: dict[str, str] = {}
    with open(args.dictionary, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            name_by_id[row["id"]] = row["standardized_name"]

    # Read tree, parse, identify tip IDs.
    text = Path(args.in_path).read_text(encoding="utf-8")
    tree = parse_newick(text)

    # Walk to collect tips and group by species.
    tips: list[dict] = []
    def collect(n):
        if not n.get("children"):
            tips.append(n)
        else:
            for c in n["children"]:
                collect(c)
    collect(tree)

    # Group tip NODES (by Python identity) per species. Multiple tip nodes can
    # share the same numeric ID — that's the redundancy we're collapsing.
    species_to_tips: dict[str, list[dict]] = {}
    skipped = []
    for tip in tips:
        tid = tip["name"]
        sp_name = name_by_id.get(tid)
        if not sp_name:
            skipped.append(tid); continue
        parts = sp_name.split("_")
        if len(parts) < 2:
            skipped.append(tid); continue
        sp = "_".join(parts[:2])
        species_to_tips.setdefault(sp, []).append(tip)

    # Pick representative tip NODE per species — deterministic: smallest id,
    # ties broken by tree-traversal order. We keep node identities so the
    # same numeric ID at a different tree position is not accidentally kept.
    keep_node_ids: set[int] = set()
    for sp, sp_tips in species_to_tips.items():
        # Sort by (int(name), index_in_tips) for stable picking.
        ordered = sorted(sp_tips, key=lambda n: (int(n["name"]) if n["name"].isdigit() else n["name"], tips.index(n)))
        keep_node_ids.add(id(ordered[0]))

    # Prune.
    pruned = prune_keep_refs(tree, keep_node_ids)
    if pruned is None:
        sys.exit("ERROR: pruning removed every tip — check dictionary or input file")

    out_text = write_newick(pruned) + ";"
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(out_text + "\n", encoding="utf-8")

    if args.report:
        # Properly count tips in the output.
        out_tree = parse_newick(out_text)
        out_tips: list[dict] = []
        def collect_out(n):
            if not n.get("children"):
                out_tips.append(n)
            else:
                for c in n["children"]: collect_out(c)
        collect_out(out_tree)
        print(f"input:  {len(tips)} tip positions ({len(species_to_tips)} unique species)", file=sys.stderr)
        print(f"output: {len(out_tips)} tip positions ({len({t['name'] for t in out_tips})} unique IDs)", file=sys.stderr)
        if skipped:
            print(f"skipped {len(skipped)} tips (no dictionary match or no Genus_species)", file=sys.stderr)


if __name__ == "__main__":
    main()
