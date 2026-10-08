#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""opentree_extract.py -- the Open Tree of Life side of finding P01.

Reads the synthesis release and its taxonomy (not distributed with the atlas):

    opentree16.1_tree.tgz  https://files.opentreeoflife.org/synthesis/
                           opentree16.1/opentree16.1_tree.tgz
                           (labelled_supertree/labelled_supertree.tre,
                            grafted_solution/grafted_solution.tre)
    ott3.7.3.tgz           https://files.opentreeoflife.org/ott/ott3.7.3/
                           ott3.7.3.tgz (taxonomy.tsv)

'Tips from phylogeny' (159,925 in the v16.1 release notes) are the leaves of
grafted_solution.tre.  Each is assigned to the clades of Data/opentree_clades.csv
by OTT ancestry (never by synthesis subtree: Crustacea and Crocodylia are
broken taxa in the synthesis).  Writes

    Data/external/opentree16.1_counts.csv
        per clade: synthesis tips, tips that are extant species, phylogeny-
        supported leaves, and those that are extant species-rank taxa;
        with the sha256 of both inputs.
    Data/external/opentree16.1_phylo_species.tsv.gz
        every phylogeny-supported extant species in an atlas clade:
        ott_id, OTT name, clade.

rebuild_atlas.py compares these with the atlas (Paper/inputs/
opentree_comparison.csv); no atlas number is computed here.

USAGE
    python3 scripts/opentree_extract.py TREE_DIR TAXONOMY_TSV TREE_TGZ OTT_TGZ
(TREE_DIR is the extracted opentree16.1_tree directory.)
"""

import csv
import gzip
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
CLADES = os.path.join(ROOT, "Data", "opentree_clades.csv")
OUT_COUNTS = os.path.join(ROOT, "Data", "external", "opentree16.1_counts.csv")
OUT_SPECIES = os.path.join(ROOT, "Data", "external",
                           "opentree16.1_phylo_species.tsv.gz")
LEAF_RE = re.compile(r"[(,](ott\d+)(?=[,)])")
RELEASE = ("opentree16.1 (synthesis 2025-12-20, OTT 3.7.3); API "
           "tree_of_life/about synth_id opentree16.1")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv):
    tree_dir, tax, tree_tgz, ott_tgz = argv[1:5]
    with open(os.path.join(tree_dir, "labelled_supertree",
                           "labelled_supertree.tre")) as fh:
        tips = [int(x[3:]) for x in LEAF_RE.findall(fh.read())]
    with open(os.path.join(tree_dir, "grafted_solution",
                           "grafted_solution.tre")) as fh:
        graft = [int(x[3:]) for x in LEAF_RE.findall(fh.read())]
    par, rank, flags, name = {}, {}, {}, {}
    with open(tax, encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            p = line.split("\t|\t")
            u = int(p[0])
            par[u] = int(p[1]) if p[1] else None
            name[u] = p[2]
            rank[u] = p[3]
            flags[u] = p[6]
    clades = []
    with open(CLADES, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            clades.append((r["partition"],
                           [int(x) for x in r["ott_ids"].split(";")]))
    target = {}
    for part, ids in clades:
        for i in ids:
            target.setdefault(i, []).append(part)

    def parts_of(u):
        out, v = [], u
        while v is not None:
            out += target.get(v, [])
            v = par.get(v)
        return out

    def extant_species(u):
        return rank.get(u) == "species" and "extinct" not in flags.get(u, "")

    keys = ("ot_tips", "ot_tips_extant_species", "ot_phylo",
            "ot_phylo_extant_species")
    cnt = {p: dict.fromkeys(keys, 0) for p, _ in clades}
    for u in tips:
        es = extant_species(u)
        for p in set(parts_of(u)):
            cnt[p]["ot_tips"] += 1
            cnt[p]["ot_tips_extant_species"] += es
    species = []
    for u in sorted(set(graft)):
        es = extant_species(u)
        ps = set(parts_of(u))
        for p in ps:
            cnt[p]["ot_phylo"] += 1
            cnt[p]["ot_phylo_extant_species"] += es
        if es:
            for p in sorted(ps):
                if not p.startswith("__"):
                    species.append((u, name[u], p))
    os.makedirs(os.path.dirname(OUT_COUNTS), exist_ok=True)
    sh_t, sh_o = sha256(tree_tgz), sha256(ott_tgz)
    with open(OUT_COUNTS, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["clade"] + list(keys) + ["release", "tree_tgz_sha256",
                                             "ott_tgz_sha256"])
        for p, _ in clades:
            w.writerow([p] + [cnt[p][k] for k in keys] + [RELEASE, sh_t, sh_o])
    with gzip.open(OUT_SPECIES, "wt", encoding="utf-8", newline="") as fh:
        fh.write("ott_id\tname\tclade\n")
        for u, nm, p in species:
            fh.write("%d\t%s\t%s\n" % (u, nm, p))
    print("grafted leaves %d, labelled tips %d; %d phylogeny-supported "
          "extant species rows" % (len(set(graft)), len(tips), len(species)))


if __name__ == "__main__":
    main(sys.argv)
