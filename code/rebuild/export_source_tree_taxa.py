#!/usr/bin/env python3
"""export_source_tree_taxa.py -- list the tip labels of the molecular source trees.

For seven groups the molecular tips are the tips of a tree of sequenced species
that the source study distributes (registers/molecular_only_trees.csv) but this
repository does not redistribute.  This script writes, for each, the labels of
that tree's first tree, as the rebuild reads them (scripts/rebuild_atlas.py
read_first_tree and tree_facts: NEXUS translate tables applied, comments
removed), one label per line, with the source file, its SHA-256 and the study
in the header.  scripts/reproduce_paper.py re-derives the molecular tips of
these groups from the lists.

It runs in the layout of the paper repository (forge:pancho/phylo-atlas),
where the source trees are under Data/raw/:

    python3 export_source_tree_taxa.py --root PATH/TO/phylo-atlas \
        --out registers/molecular_sources/source_trees
"""

import argparse
import csv
import hashlib
import importlib.util
import os
import sys


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="paper repository root")
    ap.add_argument("--out", required=True, help="output directory")
    a = ap.parse_args(argv)
    root = os.path.abspath(a.root)
    spec = importlib.util.spec_from_file_location(
        "rebuild_atlas", os.path.join(root, "scripts", "rebuild_atlas.py"))
    ra = importlib.util.module_from_spec(spec)
    sys.path.insert(0, os.path.join(root, "scripts"))
    spec.loader.exec_module(ra)
    ra.PARSE_NEWICK = ra.load_audit_parser()
    sys.setrecursionlimit(100000)

    with open(os.path.join(root, "Data", "tree_provenance.csv"), newline="",
              encoding="utf-8") as fh:
        prov = {r["partition_name"]: r for r in csv.DictReader(fh)
                if r["is_partition_canonical"] in ("yes", "TRUE", "true")}
    with open(os.path.join(root, "Data", "coverage_basis.csv"), newline="",
              encoding="utf-8") as fh:
        basis = list(csv.DictReader(fh))
    os.makedirs(a.out, exist_ok=True)
    for r in basis:
        src = r["molecular_source_file"]
        if r["f_basis"] != "source_tree_intersection" or not src.startswith("Data/raw/"):
            continue
        p = os.path.join(root, src)
        text = ra.read_first_tree(p)
        labels, _ = ra.tree_facts(text)
        with open(p, "rb") as fh:
            sha = hashlib.sha256(fh.read()).hexdigest()
        pv = prov.get(r["partition"], {})
        out = os.path.join(a.out, "%s.taxa.txt" % r["group"])
        with open(out, "w", encoding="utf-8") as fh:
            fh.write("# tip labels of the first tree of %s, the tree of sequenced "
                     "species behind the %s reference tree (%s, doi:%s)\n"
                     % (os.path.basename(src), r["partition"],
                        pv.get("citation_key", ""), pv.get("doi", "")))
            fh.write("# sha256 %s; %d labels; the tree itself is not "
                     "redistributed; its source is cited in "
                     "registers/tree_provenance.csv\n" % (sha, len(labels)))
            for lab in labels:
                fh.write(lab + "\n")
        print("%-12s %6d labels  %s" % (r["partition"], len(labels), out))


if __name__ == "__main__":
    main()
