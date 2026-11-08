#!/usr/bin/env python3
"""partition_sensitivity.py -- how much the coverage numbers depend on the partition of described
species (decision D66; Supplement S7, Table S-partition).

    python3 scripts/partition_sensitivity.py        # writes Data/partition_sensitivity.csv
    import partition_sensitivity; partition_sensitivity.compute()   # the rows, for the rebuild

Two estimands, at three levels, under five partitions.

Estimands
    A  any tree:  the described eukaryotic species named in at least one tree of the catalogue
                  that is eligible at the level. No partition.
    B  best tree: each group of the partition takes the eligible tree that names the most of
                  its species, and the group counts only that tree's species of the group
                  (the paper's one-reference-tree-per-group rule, applied at the partition's
                  scale). B <= A for every partition; B -> A as the groups become finer.

Levels (read the same way for every tree, so A and B are comparable)
    named        a tip resolves to an accepted Catalogue of Life (COL26.9) species, directly or
                 through a synonym (Data/standardized/col_resolution*.csv), the paper's C0 match.
    dated        ... in a tree that is a chronogram in absolute time: Data/time_scale_checks.csv
                 where it has the tree, otherwise truth_dated.csv is_chronogram.
    uncertainty  ... in a reference tree whose archive carries date uncertainty
                 (Data/coverage_basis.csv uncertainty_available). Other trees are not assessed.
    The paper's molecular level is enumerated only for the reference trees and is not compared.
    These levels count named species, not the paper's nested levels: the dated level here is not
    restricted to species placed with molecular data.

Partitions
    groups       the paper's groups: a species in a reference tree takes that tree's group; any
                 other species takes the group whose reference trees carry its family most often,
                 else its tree's group (the rule of the union bound, rebuild_atlas.union_bound).
    phylum, class, order, family
                 the COL26.9 classification of the species' family
                 (Data/external/col269_family_classification.csv); a species without a family
                 or a family without the rank forms a group "(none) <kingdom>".
"""

import collections
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OK = {"accepted", "synonym", "provisional", "accepted_orthographic", "synonym_orthographic"}
PROK = {"Bacteria", "Archaea"}
PARTITIONS = ("groups", "phylum", "class", "order", "family")
LEVELS = ("named", "dated", "uncertainty")
HEADER = ["partition", "groups_with_species", "level", "A_any_tree", "B_best_tree",
          "A_minus_B", "B_over_A"]


def _csv(root, *p):
    with open(os.path.join(root, *p), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def compute(root=ROOT):
    basis = _csv(root, "Data", "coverage_basis.csv")
    euk = [r for r in basis if r["is_prokaryotic"] == "no"]
    ref_tree = {r["partition"]: r["canonical_tree"] for r in euk if r["canonical_tree"]}
    unc_trees = {r["canonical_tree"] for r in euk
                 if r["canonical_tree"] and r["uncertainty_available"] == "yes"}
    parts = {r["partition"] for r in euk}

    truth_dated = {r["filename"]: r for r in _csv(root, "Data", "_rebuild", "truth_dated.csv")}
    tree_dated = {fn: (r["is_chronogram"] == "yes") for fn, r in truth_dated.items()}
    for r in _csv(root, "Data", "time_scale_checks.csv"):
        tree_dated[r["filename"]] = (r["time_scale"] == "absolute")
    fam_cls = {r["family"]: r for r in _csv(root, "Data", "external",
                                              "col269_family_classification.csv")}

    names = collections.defaultdict(set)      # tree -> accepted COL ids (eukaryotes)
    sp_family, sp_kingdom, sp_ref = {}, {}, {}
    fam_votes = collections.defaultdict(collections.Counter)
    noncanon = collections.defaultdict(list)

    for r in _csv(root, "Data", "standardized", "col_resolution.csv"):
        if r["match_type"] not in OK or not r["col_id"] or r["col_kingdom"] in PROK:
            continue
        p = r["partition"]
        t = ref_tree.get(p)
        if t is None:
            continue
        i = r["col_id"]
        names[t].add(i)
        sp_family.setdefault(i, r["col_family"])
        sp_kingdom.setdefault(i, r["col_kingdom"])
        sp_ref.setdefault(i, p)
        if r["col_family"]:
            fam_votes[r["col_family"]][p] += 1
    fam_part = {f: sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                for f, c in fam_votes.items()}
    for r in _csv(root, "Data", "standardized", "col_resolution_noncanonical.csv"):
        if r["match_type"] not in OK or not r["col_id"] or r["col_kingdom"] in PROK:
            continue
        noncanon[r["filename"]].append(r)
    tree_part = {}
    for fn, rows in noncanon.items():
        tag = (truth_dated.get(fn, {}).get("partition_name") or "").strip()
        if tag in parts:
            tree_part[fn] = tag
        else:
            c = collections.Counter(fam_part[r["col_family"]] for r in rows
                                    if r["col_family"] in fam_part)
            tree_part[fn] = (sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                             if c else "(unfiled)")
        for r in rows:
            i = r["col_id"]
            names[fn].add(i)
            sp_family.setdefault(i, r["col_family"])
            sp_kingdom.setdefault(i, r["col_kingdom"])

    def group_of(i):
        if i in sp_ref:
            return sp_ref[i]
        return fam_part.get(sp_family.get(i) or "")

    first_tree = {}
    for fn in sorted(noncanon):
        for i in names[fn]:
            first_tree.setdefault(i, fn)

    def groups_part(i):
        return group_of(i) or tree_part.get(first_tree.get(i), "(unfiled)")

    def rank_group(i, rank):
        fam = sp_family.get(i) or ""
        if rank == "family":
            return fam or "(none) " + sp_kingdom.get(i, "")
        v = (fam_cls.get(fam) or {}).get(rank, "")
        return v or "(none) " + sp_kingdom.get(i, "")

    pf = {"groups": groups_part}
    for rank in PARTITIONS[1:]:
        pf[rank] = (lambda rk: (lambda i: rank_group(i, rk)))(rank)
    eligible = {
        "named": list(names),
        "dated": [t for t in names if tree_dated.get(t)],
        "uncertainty": [t for t in names if t in unc_trees],
    }
    out = []
    for pname in PARTITIONS:
        grp = {i: pf[pname](i) for i in sp_kingdom}
        per_tree = {t: collections.Counter(grp[i] for i in ids) for t, ids in names.items()}
        for lv in LEVELS:
            trees = eligible[lv]
            union = set().union(*(names[t] for t in trees)) if trees else set()
            best = collections.Counter()
            for t in trees:
                for g, n in per_tree[t].items():
                    if n > best[g]:
                        best[g] = n
            A, B = len(union), sum(best.values())
            out.append(collections.OrderedDict([
                ("partition", pname), ("groups_with_species", len(set(grp.values()))),
                ("level", lv), ("A_any_tree", A), ("B_best_tree", B), ("A_minus_B", A - B),
                ("B_over_A", "%.4f" % (B / A) if A else "")]))
    return out


def main():
    rows = compute()
    p = os.path.join(ROOT, "Data", "partition_sensitivity.csv")
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=HEADER, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("%-8s %7s %-12s A %9s  B %9s  B/A %s" % (
            r["partition"], r["groups_with_species"], r["level"], r["A_any_tree"],
            r["B_best_tree"], r["B_over_A"]))
    print("->", os.path.relpath(p, ROOT))


if __name__ == "__main__":
    main()
