#!/usr/bin/env python3
"""run_worked_example.py -- refit the worked example from this repository.

Builds, in a temporary directory, the layout that worked_example_bd.R reads
(Data/..., Paper/inputs/crown_ages.csv) from the files of this repository,
runs the script with Rscript (needs R with ape; phytools for its cross-check),
then (1) recomputes every printed worked-example number from the fresh fits
with scripts/reproduce_paper.py --fits, which must pass, and (2) reports, for
each output table, the largest relative difference of any number from the
committed copy in paper/ (the optimiser's last printed digit can differ
between machines) and any text column that differs.

The fits over the archived tree sets (worked_example_sets.R, written to
paper/worked_example_bd_sets.csv) are not rerun here: they read the sources'
posterior and date-assignment sets, which are not redistributed
(registers/worked_example_tree_sets.csv names each set and where it is
archived; README.md, "Rebuilding the fits").  reproduce_paper.py --fits
checks that block against the committed copy.

    python3 code/rebuild/run_worked_example.py [--keep DIR]

The fits take minutes per large tree; run it on a machine with 16 GB or more.
"""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir))
OUTPUTS = ("worked_example_bd.csv", "worked_example_bd_selection.csv",
           "worked_example_bd_family_fractions.csv",
           "worked_example_bd_crosscheck.csv", "worked_example_bd_polytomy.csv",
           "worked_example_bd_family_order.csv")


def link(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    os.symlink(os.path.abspath(src), dst)


NEAR_ZERO = 1e-6     # estimates on the mu = 0 boundary print as ~1e-12


def compare(a, b):
    """(largest relative difference of numbers away from zero, list of text
    differences).  Pairs where both values are below NEAR_ZERO in magnitude
    are boundary estimates and are not ratioed."""
    with open(a, newline="", encoding="utf-8") as fh:
        ra = list(csv.DictReader(fh))
    with open(b, newline="", encoding="utf-8") as fh:
        rb = list(csv.DictReader(fh))
    if len(ra) != len(rb):
        return None, ["%d rows, committed %d" % (len(ra), len(rb))]
    worst, text = 0.0, []
    for i, (x, y) in enumerate(zip(ra, rb)):
        if set(x) != set(y):
            return None, ["columns differ: %s" % sorted(set(x) ^ set(y))]
        for k in x:
            u, v = x[k], y[k]
            try:
                fu, fv = float(u), float(v)
                if fu != fv and max(abs(fu), abs(fv)) >= NEAR_ZERO:
                    worst = max(worst, abs(fu - fv) / max(abs(fu), abs(fv)))
            except ValueError:
                if u != v:
                    text.append("row %d %s: %r, committed %r" % (i + 1, k, u, v))
    return worst, text


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", help="build the layout here and keep it")
    a = ap.parse_args(argv)
    work = a.keep or tempfile.mkdtemp(prefix="worked_example_")
    D = os.path.join(work, "Data")
    link(os.path.join(ROOT, "registers", "coverage_basis.csv"),
         os.path.join(D, "coverage_basis.csv"))
    link(os.path.join(ROOT, "registers", "molecular_tips.csv"),
         os.path.join(D, "molecular_tips.csv"))
    link(os.path.join(ROOT, "standardized", "trees"),
         os.path.join(D, "standardized", "trees"))
    link(os.path.join(ROOT, "standardized", "col269_family_species_counts.csv"),
         os.path.join(D, "standardized", "col269_family_species_counts.csv"))
    with gzip.open(os.path.join(ROOT, "standardized", "col_resolution.csv.gz"),
                   "rb") as src, open(os.path.join(
                       D, "standardized", "col_resolution.csv"), "wb") as dst:
        shutil.copyfileobj(src, dst)
    link(os.path.join(ROOT, "paper", "crown_ages.csv"),
         os.path.join(work, "Paper", "inputs", "crown_ages.csv"))
    script = os.path.join(work, "Paper", "code", "analysis", "worked_example_bd.R")
    link(os.path.join(HERE, "worked_example_bd.R"), script)
    subprocess.check_call(["Rscript", script, work])
    failed = 0
    for f in OUTPUTS:
        new = os.path.join(work, "Paper", "inputs", f)
        if not os.path.exists(new):
            print("missing output %s" % f)
            failed += 1
            continue
        worst, text = compare(new, os.path.join(ROOT, "paper", f))
        print("%-42s largest relative difference %s%s" % (
            f, "--" if worst is None else "%.1e" % worst,
            "" if not text else ", %d text difference(s)" % len(text)))
        for t in text[:10]:
            print("    " + t)
        failed += worst is None or bool(text)
    print("worked_example_bd_sets.csv: not refitted (needs the archived tree "
          "sets of registers/worked_example_tree_sets.csv); checked from paper/")
    rc = subprocess.call([sys.executable, "-I",
                          os.path.join(ROOT, "scripts", "reproduce_paper.py"),
                          "--quiet", "--fits", os.path.join(work, "Paper", "inputs")])
    failed += rc != 0
    if not a.keep:
        shutil.rmtree(work)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
