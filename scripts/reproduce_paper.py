#!/usr/bin/env python3
"""reproduce_paper.py -- recompute the paper's reported numbers from this repository.

Python 3, standard library only.  Run from anywhere:

    python3 scripts/reproduce_paper.py            # full run, under a minute
    python3 scripts/reproduce_paper.py --quiet    # failures and summary only

The script reads the redistributed trees and registers of this repository and
nothing else.  It re-implements every counting rule the paper states, from the
rules, without importing the code that produced the paper (code/rebuild/), and
compares each recomputed quantity with the committed copy of the paper's
numbers:

    paper/headline_numbers.json     every reported quantity, as computed
    paper/numbers_macros.tex        the numbers as printed in the text
    paper/crown_ages_macros.tex     crown ages and intervals as printed
    paper/table_datasets_macros.tex counts printed with Table 1
    paper/table_*.csv               Tables 1 to 3, cell by cell
    paper/fig*_*.csv, crown_ages.csv   the data of Figures 1 to 4
    paper/supp_*.csv, supp_macros.tex  the tables of Supplement S1 to S10, cell
                                    by cell, and every number its text prints
                                    (each table's .tex twin is checked to
                                    print the counts of its .csv)
    paper/worked_example_bd*.csv    the birth-death fits (Figure 5, Table 3);
                                    the fits themselves are rerun by
                                    code/rebuild/run_worked_example.py (R)
    registers/coverage_basis.csv    the per-group table
    registers/molecular_tips.csv    the molecular tips, re-derived from the
                                    distributed lists

It exits 0 when every comparison agrees and every printed number has been
recomputed, and 1 otherwise.  A few rule parameters (the interval quantiles,
the Catalogue of Life hole fraction, the cross-check size limit and the list of
named checks) are read from the pipeline's source in code/rebuild/, and the
counts of restored branch lengths from registers/restored_branch_lengths.csv.

Inputs (all in this repository):

    standardized/trees/*.nwk                 redistributed trees, numeric tips
    standardized/dictionary.csv              tip id -> standardized name
    standardized/full_mapping.csv            source label -> standardized name
    standardized/col_resolution.csv.gz       Catalogue of Life (COL26.9) match
                                             of every reference-tree tip
    standardized/col_resolution_noncanonical.csv   the same for the other trees
    standardized/col269_suprageneric_genera.txt    genera with family-like names
    standardized/canonical_sources/*.nwk     the unpruned source of each
                                             reference tree
    standardized/superseded_canonicals/*.nwk the reference trees replaced in
                                             this release (dictionary check)
    standardized/metadata.csv, withdrawn.csv per-tree metadata; withheld files
    registers/canonical_sources.csv          which source each is built from
    registers/pruned_tips.csv                tips outside the group, removed
    registers/coverage_basis.csv             groups, described species, the
                                             reference tree of each group
    registers/molecular_tips.csv             molecular tips, by tip id
    registers/molecular_sources/*.taxa.txt   taxa with sequence data in the
                                             sources' alignments
    registers/molecular_sources/source_trees/*.taxa.txt
                                             tip labels of the sources' trees
                                             of sequenced species
    registers/molecular_tip_evidence.csv     per-tip molecular evidence
    registers/time_scale_checks.csv          absolute or relative time
    registers/calibration_register.csv       the calibration on each focal
                                             node and its prior range
    registers/crown_intervals.csv            crown-node intervals of the sources
    registers/col_denominators.csv           Catalogue of Life denominators
    registers/tree_provenance.csv            assembly fields of the sources
    registers/dating_deposit_checks.csv      what each source study dated
    registers/molecular_only_trees.csv       trees of molecular species
    registers/dataset_table_registry.csv     the rows of Table 1
    registers/external/gbif_unmatched.csv    GBIF matches of unmatched tips
    registers/external/opentree16.1_*        Open Tree synthesis 16.1 counts
                                             and phylogeny-placed species
    registers/_rebuild/truth_*.csv           per-file group and imputation
                                             fields read by the rebuild
    registers/reference_tree_sources.csv     the source file of each
                                             reference tree and the rule that
                                             chose it
    registers/restored_branch_lengths.csv    the trees whose branch lengths
                                             release v1.0.9 restored
    standardized/name_overrides.csv, name_rejected_targets.csv
                                             name-resolution registers
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import re
import sys
from statistics import NormalDist
from collections import Counter, OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))
csv.field_size_limit(1 << 30)


def path(*parts):
    return os.path.join(ROOT, *parts)


# ---------------------------------------------------------------------------
# Rules stated in the paper's Methods and Supplement (parameters, not results)
# ---------------------------------------------------------------------------

CHRONOGRAM_CV_MAX = 0.05            # root-to-tip coefficient of variation
COL_RESOLVED = {"accepted", "provisional", "synonym", "accepted_orthographic",
                "provisional_orthographic", "synonym_orthographic"}
INTERVAL_KIND_COUNTS = {            # does the spread estimate a divergence time
    "node_hpd_annotation": True, "node_label_interval": True,
    "posterior_sample": True, "bootstrap_dating_sample": True,
    "calibration_prior_sample": False, "taxon_addition_sample": False,
    "date_assignment_sample": False}
PRIOR_SAMPLE_KINDS = {"calibration_prior_sample"}
DEGENERATE_WIDTH_MYR = 1e-3
BOUND_TOL_MYR = 1e-2
PRIOR_CONTAINMENT_MIN = 0.9
PRIOR_WIDTH_RATIO_MIN = 0.5
INTERVAL_PCT = 95
CROWN_POINT_TOL_MYR = 1e-3
CROWN_CENTER_KINDS = {"node_hpd_annotation", "posterior_sample",
                      "bootstrap_dating_sample", "calibration_prior_sample"}
CROWN_SUMMARY_KINDS = {"node_hpd_annotation", "node_label_interval"}
GBIF_MATCH = {"EXACT", "FUZZY"}
GBIF_STATUS = {"ACCEPTED", "SYNONYM", "HETEROTYPIC_SYNONYM",
               "HOMOTYPIC_SYNONYM", "PROPARTE_SYNONYM", "DOUBTFUL"}
CATEGORY_ORDER = ["Vertebrates", "Plants", "Arthropods", "Other animals",
                  "Microbes & protists"]
LEVELS = ("C0", "C1", "C1d", "C1b")
# Table 1 marks a group whose tips exceed its named species by this fraction
TIPS_OVER_NAMED = 0.10
# interval kinds that come from a set of trees
SET_KINDS = ("posterior_sample", "bootstrap_dating_sample", "calibration_prior_sample",
             "taxon_addition_sample", "date_assignment_sample")
LOSS_REASONS = ("NotChron", "RelTime", "NoInterval", "Kind", "Calib")
UNREP_NAMED = 4                     # largest groups without a tree, named
# Agreement the Supplement states as tolerances (the differences themselves are
# rounding noise that changes between machines)
UNIFORM_TOL = 1e-6                  # family vs uniform log-likelihood at equal f
CROSSCHECK_RTOL = 1e-4              # rates against phytools::fit.bd, relative
TEST_LEVEL = 0.05                   # level of the test for mu = 0
# Under mu = 0 the statistic is an equal mixture of a point mass at zero and a
# chi-squared with one degree of freedom, so its critical value at level a is
# qchisq(1 - 2a, 1) = z_{1-a}^2 (2.71 at a = 0.05).
LR_CRITICAL = NormalDist().inv_cdf(1.0 - TEST_LEVEL) ** 2
# Species key (decisions D4, D13): see species_key() below.
OPEN_MARKERS = {"sp", "spp", "sp.", "spp.", "cf", "cf.", "aff", "aff.", "af",
                "nr", "nr.", "indet", "indet.", "x", "hybrid", "nov", "n",
                "undescribed", "unidentified", "unknown", "gen", "trans",
                "genome", "sp1", "sp2", "sp3", "sp4", "sp5"}
HYBRID_OR_OPEN = {"x", "hybrid", "cf", "cf.", "aff", "aff.", "af"}
SUPRAGENERIC_SUFFIXES = ("idae", "aceae", "inae", "oidea", "ales", "phyta",
                         "mycota", "mycetes", "opsida", "ineae")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_csv(p, comment=False):
    opener = gzip.open if p.endswith(".gz") else open
    with opener(p, "rt", newline="", encoding="utf-8") as fh:
        lines = (ln for ln in fh if not (comment and ln.startswith("#")))
        return list(csv.DictReader(lines))


def as_int(v, default=0):
    try:
        return int(float(str(v).replace(",", "")))
    except (TypeError, ValueError):
        return default


def as_float(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def tree_arrays(text):
    """parent index, branch length, label and child count of every node of a
    Newick string (nodes in the order they open, parents before children)."""
    s, quote = [], None
    for ch in text.strip():
        if quote:
            s.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            s.append(ch)
        elif not ch.isspace():
            s.append(ch)
    s = "".join(s)
    if s.endswith(";"):
        s = s[:-1]
    parent, bl, name, nkids = [-1], [0.0], [""], [0]
    cur, i, n = 0, 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            parent.append(cur); bl.append(0.0); name.append(""); nkids.append(0)
            nkids[cur] += 1
            cur = len(parent) - 1
            i += 1
        elif c == ",":
            p = parent[cur]
            parent.append(p); bl.append(0.0); name.append(""); nkids.append(0)
            nkids[p] += 1
            cur = len(parent) - 1
            i += 1
        elif c == ")":
            cur = parent[cur]
            i += 1
        else:
            j = i
            while j < n and s[j] not in "(),":
                j += 1
            lab, _, b = s[i:j].partition(":")
            name[cur] = lab.strip()
            v = as_float(b.strip(), 0.0) if b.strip() else 0.0
            bl[cur] = v if v == v else 0.0
            i = j
    return parent, bl, name, nkids


def polytomy_counts(text, keep, scale, tol=1e-6):
    """(branching events, events at nodes with more than two children) of a
    tree restricted to the tips in `keep` (None: all), after branches shorter
    than `tol` (after scaling) between two branching nodes are collapsed."""
    parent, bl, name, nkids = tree_arrays(text)
    kids = [[] for _ in parent]
    for k in range(1, len(parent)):
        kids[parent[k]].append(k)
    res = [None] * len(parent)
    tot = [0, 0]

    def fin(k_):
        tot[0] += k_ - 1
        if k_ > 2:
            tot[1] += k_ - 1
    for k in range(len(parent) - 1, -1, -1):
        b = bl[k] * scale
        if not kids[k]:
            res[k] = ("tip", 0, b) if (keep is None or name[k] in keep) else None
            continue
        sub = [res[c] for c in kids[k] if res[c] is not None]
        if not sub:
            continue
        if len(sub) == 1:
            res[k] = (sub[0][0], sub[0][1], sub[0][2] + b)
            continue
        k_eff = 0
        for kind, kk, ln in sub:
            if kind == "node" and ln < tol:
                k_eff += kk
            else:
                k_eff += 1
                if kind == "node":
                    fin(kk)
        res[k] = ("node", k_eff, b)
    if res[0] and res[0][0] == "node":
        fin(res[0][1])
    return tot[0], tot[1]


def leaf_depths(text):
    """Leaves of a Newick string as (label, root-to-tip depth) pairs, and the
    number of NaN branch lengths.  Whitespace outside quoted labels is
    ignored; the root's own branch length, if written, is included; a NaN
    branch length contributes zero depth."""
    s, quote = [], None
    for ch in text.strip():
        if quote:
            s.append(ch)
            if ch == quote:
                quote = None
        elif ch in "'\"":
            quote = ch
            s.append(ch)
        elif not ch.isspace():
            s.append(ch)
    s = "".join(s)
    if s.endswith(";"):
        s = s[:-1]
    parent, bl, name, nkids = [-1], [0.0], [""], [0]
    nan = 0
    cur, i, n = 0, 0, len(s)
    while i < n:
        c = s[i]
        if c == "(":
            parent.append(cur); bl.append(0.0); name.append(""); nkids.append(0)
            nkids[cur] += 1
            cur = len(parent) - 1
            i += 1
        elif c == ",":
            p = parent[cur]
            parent.append(p); bl.append(0.0); name.append(""); nkids.append(0)
            nkids[p] += 1
            cur = len(parent) - 1
            i += 1
        elif c == ")":
            cur = parent[cur]
            i += 1
        else:
            j = i
            while j < n and s[j] not in "(),":
                j += 1
            lab, _, b = s[i:j].partition(":")
            name[cur] = lab.strip()
            b = b.strip()
            if b.lower() in ("nan", "-nan"):
                nan += 1
                v = 0.0
            else:
                v = as_float(b, 0.0) if b else 0.0
                if v != v:
                    nan += 1
                    v = 0.0
            bl[cur] = v
            i = j
    depth = [0.0] * len(parent)
    out = []
    for k in range(len(parent)):        # parents are created before children
        depth[k] = (depth[parent[k]] if parent[k] >= 0 else 0.0) + bl[k]
    for k in range(len(parent)):
        if nkids[k] == 0:
            out.append((name[k], depth[k]))
    return out, nan


def cv(xs):
    if not xs:
        return None
    m = sum(xs) / len(xs)
    if m == 0:
        return None
    return math.sqrt(sum((x - m) ** 2 for x in xs) / len(xs)) / abs(m)


def load_suprageneric_genera():
    out = set()
    with open(path("standardized", "col269_suprageneric_genera.txt"),
              encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if ln and not ln.startswith("#"):
                out.add(ln)
    return out


GENUS_RE = re.compile(r"^[A-Z][a-z]+$")
EPITHET_RE = re.compile(r"^[a-z][a-z-]+$")
SPLIT_RE = re.compile(r"[_.\s]+")
SUPRA_GENERA = None


def species_key(name):
    """(key, resolved) for a standardized name (decisions D4, D13).  The key is
    the first adjacent 'Genus epithet' pair, the genus not ending in a
    family-level suffix unless the Catalogue of Life lists it as a genus.  A
    name with an open-nomenclature or hybrid marker is not a species label."""
    toks = [t for t in SPLIT_RE.split(name or "") if t]
    if {t.lower() for t in toks} & HYBRID_OR_OPEN:
        return None, False
    for a, b in zip(toks, toks[1:]):
        genus = bool(GENUS_RE.match(a)) and (
            a in SUPRA_GENERA or not a.lower().endswith(SUPRAGENERIC_SUFFIXES))
        if genus and EPITHET_RE.match(b):
            if b.lower() in OPEN_MARKERS or re.match(r"^sp\d*$", b):
                return None, False
            return a + "_" + b, True
        if GENUS_RE.match(a) and b.lower() in OPEN_MARKERS:
            return None, False
    return None, False


def key_of(name):
    k, ok = species_key(name)
    return k if ok else "label:" + name


def macro_name(s):
    return "".join(w[:1].upper() + w[1:] for w in re.split(r"[^A-Za-z]+", s) if w)


def fmt_int(n):
    return "{:,}".format(int(round(float(n)))).replace(",", "{,}")


def fmt_pct(x):
    return "%.1f" % x if abs(x) >= 1.0 else "%.2f" % x


def name_list(names):
    """'A, B and C', as group names are printed in running text."""
    names = list(names)
    return names[0] if len(names) == 1 else \
        ", ".join(names[:-1]) + " and " + names[-1]


def pct(a, b):
    return 100.0 * float(a) / float(b)


def rule_from_code(script, name, block=False, r_assign=False):
    """A rule parameter of the pipeline, read from its source in code/rebuild/
    (the scripts are not imported: several read the source deposits)."""
    import ast
    with open(path("code", "rebuild", script), encoding="utf-8") as fh:
        src = fh.read()
    if block:
        m = re.search(r"^%s = (\[.*?^\])" % name, src, re.M | re.S)
    elif r_assign:
        m = re.search(r"^%s <- ([^#\n]+)" % name, src, re.M)
    else:
        m = re.search(r"^%s = ([^#\n]+)" % name, src, re.M)
    if not m:
        sys.exit("cannot read %s from code/rebuild/%s" % (name, script))
    return ast.literal_eval(m.group(1).strip())


def read_macros(*files):
    out = {}
    rx = re.compile(r"\\newcommand\{\\([A-Za-z]+)\}\{(.*)\}\s*$")
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for ln in fh:
                m = rx.match(ln.strip())
                if m:
                    out[m.group(1)] = m.group(2)
    return out


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

class Report(object):
    def __init__(self, quiet):
        self.quiet = quiet
        self.n = 0
        self.fail = []
        self.macros_checked = set()
        self.supp_macros = {}

    def section(self, title):
        if not self.quiet:
            print("\n== %s" % title)

    def eq(self, label, got, want, tol=None):
        self.n += 1
        if tol is not None and isinstance(got, (int, float)) \
                and isinstance(want, (int, float)):
            ok = abs(got - want) <= tol
        else:
            ok = got == want
        if not ok:
            self.fail.append("%s: recomputed %r, paper %r" % (label, got, want))
            print("  FAIL %s: recomputed %r, paper %r" % (label, got, want))
        elif not self.quiet:
            print("  ok   %s = %s" % (label, got))
        return ok

    def macro(self, macros, name, got):
        self.macros_checked.add(name)
        if name not in macros:
            self.n += 1
            self.fail.append("macro \\%s is not in the paper's macro files" % name)
            print("  FAIL macro \\%s is not in the paper's macro files" % name)
            return
        self.eq("\\" + name, got, macros[name])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--quiet", action="store_true",
                    help="print failures and the summary only")
    ap.add_argument("--fits", default=path("paper"),
                    help="directory holding worked_example_bd*.csv (default "
                         "paper/; code/rebuild/run_worked_example.py passes "
                         "a fresh refit)")
    a = ap.parse_args(argv)
    R = Report(a.quiet)
    global SUPRA_GENERA
    SUPRA_GENERA = load_suprageneric_genera()

    head = json.load(open(path("paper", "headline_numbers.json"), encoding="utf-8"))
    macros = read_macros(path("paper", "numbers_macros.tex"),
                         path("paper", "crown_ages_macros.tex"),
                         path("paper", "table_datasets_macros.tex"))

    # ---------------------------------------------------------- inputs ---
    basis = read_csv(path("registers", "coverage_basis.csv"))
    bmap = {r["partition"]: r for r in basis}
    euk = [r for r in basis if r["is_prokaryotic"] == "no"]
    prok = [r for r in basis if r["is_prokaryotic"] == "yes"]
    canon = [r for r in basis if r["canonical_tree"]]
    canon_euk = [r for r in canon if r["is_prokaryotic"] == "no"]
    den = sum(as_int(r["described"]) for r in euk)

    dictionary = {}
    with open(path("standardized", "dictionary.csv"), newline="",
              encoding="utf-8") as fh:
        for row in csv.reader(fh):
            if len(row) >= 2 and row[0] != "id":
                dictionary[row[0]] = row[1]
    by_group = defaultdict(dict)
    for r in read_csv(path("standardized", "full_mapping.csv")):
        by_group[r["group"]][r["original_label"]] = r["standardized_name"]
    col_rows = defaultdict(dict)
    for r in read_csv(path("standardized", "col_resolution.csv.gz")):
        col_rows[r["partition"]][r["standardized_id"]] = r
    time_scale = {r["filename"]: r["time_scale"]
                  for r in read_csv(path("registers", "time_scale_checks.csv"))}
    intervals = {r["partition"]: r
                 for r in read_csv(path("registers", "crown_intervals.csv"))}
    truth_dated = {r["filename"]: r for r in read_csv(
        path("registers", "_rebuild", "truth_dated.csv"))}
    truth_prov = {r["partition_name"]: r for r in read_csv(
        path("registers", "_rebuild", "truth_provenance.csv"))
        if r["is_canonical"] == "yes"}
    mol_rows = defaultdict(set)
    for r in read_csv(path("registers", "molecular_tips.csv")):
        mol_rows[r["partition"]].add(r["standardized_id"])

    # ------------------------------------------------ redistributed trees ---
    R.section("Redistributed trees: tip counts and root-to-tip variation")
    tree_dir = path("standardized", "trees")
    files = sorted(f for f in os.listdir(tree_dir) if f.endswith(".nwk"))
    facts = {}
    for f in files:
        with open(os.path.join(tree_dir, f), encoding="utf-8",
                  errors="replace") as fh:
            leaves, nan = leaf_depths(fh.read())
        depths = [d for _, d in leaves]
        facts[f] = dict(labels=[l for l, _ in leaves], depths=depths, nan=nan,
                        cv=cv(depths), root=max(depths) if depths else None)
    meta = read_csv(path("standardized", "metadata.csv"))
    R.eq("trees in standardized/trees", len(files), head["trees"]["n_shipped"])
    R.eq("trees listed in standardized/metadata.csv", sorted(m["filename"] for m in meta), files)
    bad_ntips = [m["filename"] for m in meta
                 if as_int(m["ntips"]) != len(set(facts[m["filename"]]["labels"]))
                 and as_int(m["ntips"]) != len(facts[m["filename"]]["labels"])]
    R.eq("metadata.csv ntips that match no tip count of the file", bad_ntips, [])
    n_chron_all = sum(1 for f in files if facts[f]["cv"] is not None
                      and facts[f]["cv"] <= CHRONOGRAM_CV_MAX)
    R.eq("trees with root-to-tip CV <= %g" % CHRONOGRAM_CV_MAX, n_chron_all,
         head["trees"]["n_chronograms_all_shipped"])
    R.macro(macros, "NumTreesShipped", fmt_int(len(files)))
    R.macro(macros, "NumTreesChron", fmt_int(n_chron_all))
    cv_bad = [f for f in files if f in truth_dated and (
        ("%.6f" % facts[f]["cv"] if facts[f]["cv"] is not None else "")
        != truth_dated[f]["roottotip_cv"]
        or str(len(facts[f]["labels"])) != truth_dated[f]["leaves"])]
    R.eq("files whose CV or leaf count differs from registers/_rebuild/truth_dated.csv",
         cv_bad, [])
    withdrawn = [r["filename"] for r in read_csv(path("standardized", "withdrawn.csv"))]
    R.eq("withdrawn files present in standardized/trees",
         [f for f in withdrawn if f in files], [])

    # ------------------------------------- reference trees, built by rule ---
    R.section("Reference trees: sources minus off-group tips, one tip per species")
    pruned = defaultdict(set)
    src_tips = {}
    for r in read_csv(path("registers", "pruned_tips.csv")):
        pruned[r["filename"]].add(r["standardized_id"])
    for r in read_csv(path("registers", "canonical_sources.csv")):
        f = r["filename"]
        src = path("standardized", "canonical_sources", os.path.basename(r["source"]))
        with open(src, encoding="utf-8", errors="replace") as fh:
            sleaves, _ = leaf_depths(fh.read())
        src_tips[f] = len(sleaves)
        seen, kept = set(), []
        missing = pruned[f] - {l for l, _ in sleaves}
        for lab, d in sleaves:
            if lab in pruned[f]:
                continue
            k, ok = species_key(dictionary.get(lab, ""))
            k = k if ok else "id:" + lab
            if k in seen:
                continue
            seen.add(k)
            kept.append((lab, d))
        shipped = list(zip(facts[f]["labels"], facts[f]["depths"]))
        R.eq("%s: tips = source tips - pruned - repeated species, in order" % f,
             ([l for l, _ in kept] == [l for l, _ in shipped], sorted(missing)),
             (True, []))
        if [l for l, _ in kept] == [l for l, _ in shipped] and kept:
            offs = [ds - dk for (_, ds), (_, dk) in zip(kept, shipped)]
            scale = max(abs(d) for _, d in kept) or 1.0
            R.eq("%s: root-to-tip depths equal the source's less one offset" % f,
                 (max(offs) - min(offs)) <= 1e-6 * scale, True)

    srcdir = path("standardized", "canonical_sources")
    used = set()
    for f in files:
        used.update(facts[f]["labels"])
    for d_ in (srcdir, path("standardized", "superseded_canonicals")):
        for f in os.listdir(d_):
            with open(os.path.join(d_, f), encoding="utf-8", errors="replace") as fh:
                used.update(l for l, _ in leaf_depths(fh.read())[0])
    R.eq("tip ids missing from standardized/dictionary.csv",
         len(used - set(dictionary)), 0)
    R.eq("dictionary ids that occur in no distributed tree",
         len(set(dictionary) - used), 0)

    # ------------------------------------------------------ named (C0) ---
    R.section("Named coverage: tips resolving to Catalogue of Life accepted names")
    tips = {}           # partition -> {name: first id}
    unit_of = {}        # partition -> {id: counting unit}
    c0 = {}             # partition -> set of units
    unresolved = {}
    for r in canon:
        p, f = r["partition"], r["canonical_tree"]
        first = {}
        for lab in facts[f]["labels"]:
            nm = dictionary.get(lab, lab)
            if nm not in first or lab < first[nm]:
                first[nm] = lab
        tips[p] = first
        if r["is_prokaryotic"] == "yes":
            unit_of[p] = {sid: nm for nm, sid in first.items()}
            c0[p] = set(first)
            unresolved[p] = 0
            continue
        rows = col_rows.get(p, {})
        ids = set(first.values())
        R.eq("%s: CoL rows cover exactly the tree's tips" % p,
             sorted(ids ^ set(rows))[:5], [])
        u = {sid: rows[sid]["col_id"] for sid in ids if sid in rows
             and rows[sid]["match_type"] in COL_RESOLVED and rows[sid]["col_id"]}
        unit_of[p] = u
        c0[p] = set(u.values())
        unresolved[p] = len(ids) - len(u)

    # --------------------------------------------------- molecular (C1) ---
    R.section("Molecular coverage: tips in the molecular analysis behind each tree")
    c1, c1_basis, fmol = {}, {}, {}
    for r in canon:
        p = r["partition"]
        names_in_tree = tips[p]
        if mol_rows.get(p):
            ids = mol_rows[p]
            c1[p] = {unit_of[p][sid] for sid in ids if sid in unit_of[p]}
            c1_basis[p] = "enumerated"
            mol_names = {dictionary.get(sid, sid) for sid in ids} & set(names_in_tree)
            fmol[p] = len(mol_names) / float(len(names_in_tree))
        else:
            tp = truth_prov.get(p, {})
            asserted = (tp.get("imputation_status") == "none"
                        and (as_float(tp.get("molecular_tip_fraction")) or 0) >= 1.0)
            R.eq("%s: no molecular tip list, every tip asserted molecular "
                 "(imputation_status none, fraction 1)" % p, asserted, True)
            c1[p] = set(c0[p])
            c1_basis[p] = "all_tips"
            fmol[p] = 1.0

    # Independent re-derivation of the molecular tip lists from the
    # redistributed sources: alignment taxa lists and per-tip evidence.
    evidence = defaultdict(list)
    for e in read_csv(path("registers", "molecular_tip_evidence.csv")):
        evidence[e["partition"]].append(e)
    not_rederived = []
    exact_tips = {}     # verbatim source labels equal to a tip's name
    for r in canon:
        p = r["partition"]
        if c1_basis[p] != "enumerated":
            continue
        src = r["molecular_source_file"]
        names_in_tree = set(tips[p])
        listed = {dictionary.get(s, s) for s in mol_rows[p]} & names_in_tree
        if src.startswith("Data/molecular_sources/") or src.startswith("Data/raw/"):
            fn = (path("registers", "molecular_sources", os.path.basename(src))
                  if src.startswith("Data/molecular_sources/") else
                  path("registers", "molecular_sources", "source_trees",
                       r["group"] + ".taxa.txt"))
            with open(fn, encoding="utf-8") as fh:
                labels = [ln.strip() for ln in fh if ln.strip()
                          and not ln.startswith("#")]
            gm = by_group.get(r["group"], {})
            keys = {key_of(gm.get(l, l)) for l in labels}
            derived = {nm for nm in names_in_tree if key_of(nm) in keys}
            exact_tips[p] = len(set(labels) & names_in_tree)
            R.eq("%s: molecular tips re-derived from %s, differences from "
                 "registers/molecular_tips.csv" % (p, os.path.relpath(fn, ROOT)),
                 len(derived ^ listed), 0)
        elif src == "Data/molecular_tip_evidence.csv":
            derived = {dictionary.get(e["standardized_id"]) for e in evidence[p]}
            derived = {d for d in derived if d in names_in_tree}
            R.eq("%s: molecular tips re-derived from per-tip evidence, "
                 "differences from registers/molecular_tips.csv" % p,
                 len(derived ^ listed), 0)
        else:
            not_rederived.append((p, src))
            R.eq("%s: molecular tip source is redistributed" % p, src, "")

    # ---------------------------------------------------- dated (C1d) ---
    R.section("Dated coverage: chronograms in absolute time")
    chron, reltime = {}, set()
    for r in canon:
        p, f = r["partition"], r["canonical_tree"]
        ultra = facts[f]["cv"] is not None and facts[f]["cv"] <= CHRONOGRAM_CV_MAX
        ts = time_scale.get(f, "") if ultra else "not_ultrametric"
        if ultra and ts not in ("absolute", "relative"):
            R.eq("%s: ultrametric file has a time-scale row" % p, ts, "absolute|relative")
        chron[p] = ultra and ts == "absolute"
        if ultra and ts == "relative" and r["is_prokaryotic"] == "no":
            reltime.add(p)
        R.eq("%s: chronogram, time scale, CV" % p,
             ("yes" if chron[p] else "no", ts, "%.6f" % facts[f]["cv"]),
             (r["is_chronogram"], r["time_scale"], r["roottotip_cv"]))

    # ------------------------------------------------- interval (C1b) ---
    R.section("Interval coverage: crown intervals not set by a calibration")
    counts, why_of = {}, {}
    for p, iv in sorted(intervals.items()):
        lo, hi = float(iv["interval_lo_ma"]), float(iv["interval_hi_ma"])
        center = as_float(iv["interval_center_ma"])
        point = as_float(iv["source_point_height_ma"])
        c = center if center is not None else point
        w = hi - lo
        degenerate = w < DEGENERATE_WIDTH_MYR
        at_bound = (not degenerate and c is not None
                    and (c - lo < BOUND_TOL_MYR or hi - c < BOUND_TOL_MYR))
        repro = False
        if iv["calibration"] == "two_sided":
            plo, phi = float(iv["prior_lo_ma"]), float(iv["prior_hi_ma"])
            cont = max(0.0, min(hi, phi) - max(lo, plo)) / w if w > 0 else 1.0
            ratio = w / (phi - plo)
            repro = cont >= PRIOR_CONTAINMENT_MIN and ratio >= PRIOR_WIDTH_RATIO_MIN
            R.eq("%s: prior containment and width ratio" % p,
                 ("%.6f" % cont, "%.6f" % ratio),
                 (iv["prior_containment"], iv["prior_width_ratio"]))
        why = ("prior_reproduced" if iv["kind"] in PRIOR_SAMPLE_KINDS
               else "kind" if not INTERVAL_KIND_COUNTS[iv["kind"]]
               else "degenerate" if degenerate else "at_bound" if at_bound
               else "prior_reproduced" if repro else "")
        counts[p] = not why
        why_of[p] = why
        R.eq("%s: interval classification" % p,
             (why, "yes" if not why else "no"),
             (iv["c1b_exclusion"], iv["counts_for_c1b"]))
    calreg = {r["partition"]: r for r in read_csv(
        path("registers", "calibration_register.csv"))}
    for p, iv in sorted(intervals.items()):
        if INTERVAL_KIND_COUNTS[iv["kind"]]:
            cr = calreg.get(p, {})
            # the prior range enters the test only for a two-sided calibration
            two = cr.get("calibration") == "two_sided"
            f6 = lambda x: "%.6f" % float(x) if x and two else ""
            R.eq("%s: calibration and prior range as in the calibration register" % p,
                 (iv["calibration"], f6(iv["prior_lo_ma"]), f6(iv["prior_hi_ma"])),
                 (cr.get("calibration"), f6(cr.get("prior_lo_ma")), f6(cr.get("prior_hi_ma"))))
    bounded = {p: chron.get(p, False) and counts.get(p, False) for p in c0}

    # Per-group table columns
    R.section("Per-group table (registers/coverage_basis.csv)")
    for r in canon:
        p = r["partition"]
        R.eq("%s: named, molecular, basis, unresolved, tips, molecular fraction" % p,
             (len(c0[p]), len(c1[p]), c1_basis[p], unresolved[p], len(tips[p]),
              "%.6f" % fmol[p], "yes" if bounded[p] else "no"),
             (as_int(r["c0_units"]), as_int(r["c1_units"]), r["c1_basis"],
              as_int(r["c0_unresolved_labels"]), as_int(r["labels_unique"]),
              r["f_mol"], r["interval_recoverable"]))

    # ------------------------------------------------------- the levels ---
    R.section("The four levels over the %d eukaryotic groups" % len(euk))
    E = [r["partition"] for r in canon_euk]

    def level_sets(p, bset=None):
        b = bounded if bset is None else bset
        return (c0[p], c1[p], c1[p] if chron[p] else set(),
                c1[p] if chron[p] and b[p] else set())

    def ladder(parts, denominator, bset=None):
        out = OrderedDict()
        for i, lv in enumerate(LEVELS):
            sets = [level_sets(p, bset)[i] for p in parts]
            s = sum(len(x) for x in sets)
            u = len(set().union(*sets)) if sets else 0
            out[lv] = dict(sum=s, union=u, pct_sum=pct(s, denominator),
                           pct_union=pct(u, denominator),
                           n_part=sum(1 for x in sets if x))
        return out

    lad = ladder(E, den)
    R.eq("eukaryotic denominator", den, head["denominators"]["eukaryotic"]["described"])
    for lv, nm in zip(LEVELS, ("CZero", "COne", "COneD", "COneB")):
        h = head["rungs"][lv]
        R.eq("%s numerator, sum over groups" % lv, lad[lv]["sum"],
             int(round(h["numerator_sum"])))
        R.eq("%s numerator, union" % lv, lad[lv]["union"], h["numerator_union"])
        R.eq("%s percent (union)" % lv, lad[lv]["pct_union"], h["percent_union"], 1e-9)
        R.eq("%s contributing groups" % lv, lad[lv]["n_part"],
             h["n_partitions_contributing"])
        R.macro(macros, nm + "N", fmt_int(lad[lv]["union"]))
        R.macro(macros, nm + "Pct", fmt_pct(lad[lv]["pct_union"]))
        R.macro(macros, nm + "NPart", fmt_int(lad[lv]["n_part"]))
    R.macro(macros, "COneNEnumerated", fmt_int(sum(c1_basis[p] == "enumerated" for p in E)))
    R.macro(macros, "COneNAllTips", fmt_int(sum(c1_basis[p] == "all_tips" for p in E)))
    asserted = sum(len(c1[p]) for p in E if c1_basis[p] == "all_tips")
    R.macro(macros, "COneAssertedN", fmt_int(asserted))
    R.macro(macros, "COneAssertedPct", fmt_pct(pct(asserted, lad["C1"]["sum"])))
    for (a_, an), (b_, bn) in zip((("C1", "COne"), ("C1d", "COneD"), ("C1b", "COneB")),
                                  (("C0", "CZero"), ("C1", "COne"), ("C1d", "COneD"))):
        R.macro(macros, an + "Of" + bn, fmt_pct(pct(lad[a_]["union"], lad[b_]["union"])))

    # Counts of groups and trees
    R.section("Groups, reference trees and denominators")
    R.macro(macros, "NumPartitionsAll", fmt_int(len(basis)))
    R.macro(macros, "NumPartitionsEuk", fmt_int(len(euk)))
    R.macro(macros, "NumPartitionsProk", fmt_int(len(prok)))
    R.macro(macros, "NumCanonical", fmt_int(len(canon)))
    R.macro(macros, "NumCanonicalEuk", fmt_int(len(canon_euk)))
    R.macro(macros, "NumCanonicalChron", fmt_int(sum(chron[p] for p in E)))
    R.macro(macros, "NumCanonicalRelTime", fmt_int(len(reltime)))
    R.macro(macros, "DescribedEuk", fmt_int(den))
    den_prok = sum(as_int(r["described"]) for r in prok)
    R.macro(macros, "DescribedProk", fmt_int(den_prok))
    represented = [p for p in E if c0[p]]
    den_rep = sum(as_int(bmap[p]["described"]) for p in represented)
    R.macro(macros, "DescribedRepresented", fmt_int(den_rep))
    unrep = sorted((r for r in euk if not r["canonical_tree"]),
                   key=lambda r: (-as_int(r["described"]), r["partition"]))
    den_unrep = sum(as_int(r["described"]) for r in unrep)
    R.macro(macros, "UnrepresentedPct", fmt_pct(pct(den_unrep, den)))
    big = [r["partition"] for r in unrep[:UNREP_NAMED]]
    R.macro(macros, "UnrepLargestNames", name_list(big))
    reg = read_csv(path("registers", "dataset_table_registry.csv"))
    listed, n_rows, n_agg = set(), 0, 0
    for g in reg:
        flags = set((g["flags"] or "").split(";"))
        if "aggregate" in flags:
            rx = re.compile("^" + re.escape(g["filename"]).replace(r"\*", ".*") + "$")
            hit = [f for f in files if rx.match(f)]
            listed.update(hit)
            n_agg += len(hit)
            n_rows += 1
        elif g["filename"] in files:
            listed.add(g["filename"])
            n_rows += 1
    R.eq("redistributed trees not listed in Table 1's registry",
         sorted(set(files) - listed), [])
    R.macro(macros, "DTnRows", str(n_rows))
    R.macro(macros, "DTnAggregate", str(n_agg))
    R.macro(macros, "DTnCanonical", str(len(canon)))
    R.macro(macros, "DTnCanonicalEuk", str(len(canon_euk)))
    R.macro(macros, "DTnUnrepresented", str(len(unrep)))
    R.macro(macros, "DTdescribedUnrepresented", fmt_int(den_unrep))
    R.macro(macros, "DTdescribedEuk", fmt_int(den))
    pl = head["prokaryotic_line"]
    R.eq("prokaryotic line: GTDB clusters named / described",
         (sum(len(c0[r["partition"]]) for r in prok), den_prok),
         (pl["C0_sum"], pl["described"]))

    # Per group, as printed
    R.section("Per group, as printed")
    for p in E:
        n = macro_name(p)
        d = as_int(bmap[p]["described"])
        R.macro(macros, "Described" + n, fmt_int(d))
        R.macro(macros, "CZeroN" + n, fmt_int(len(c0[p])))
        R.macro(macros, "CZeroPct" + n, fmt_pct(pct(len(c0[p]), d)))
        R.macro(macros, "COneN" + n, fmt_int(len(c1[p])))
        R.macro(macros, "COnePct" + n, fmt_pct(pct(len(c1[p]), d)))
        R.macro(macros, "FMol" + n, "%.2f" % float("%.6f" % fmol[p]))
        R.macro(macros, "TipsNotNamed" + n, fmt_int(len(tips[p]) - len(c0[p])))
    R.macro(macros, "TipsNotNamedN", fmt_int(sum(len(tips[p]) - len(c0[p]) for p in E)))
    R.macro(macros, "SensTipsCZeroPct", fmt_pct(pct(sum(
        min(len(tips[p]), as_int(bmap[p]["described"])) for p in E), den)))
    labels_euk = sum(len(tips[p]) for p in E)
    unres = sum(unresolved[p] for p in E)
    R.macro(macros, "LabelsEukN", fmt_int(labels_euk))
    R.macro(macros, "UnresolvedLabelsN", fmt_int(unres))
    R.macro(macros, "UnresolvedLabelsPct", fmt_pct(pct(unres, labels_euk)))

    # Per category (sum over the groups of a category)
    R.section("Per category")
    for cat in CATEGORY_ORDER:
        rows = [r for r in euk if r["category"] == cat]
        d = sum(as_int(r["described"]) for r in rows)
        parts = [r["partition"] for r in rows if r["canonical_tree"]]
        v = [sum(len(level_sets(p)[i]) for p in parts) for i in range(3)]
        h = head["per_category_eukaryotic"][cat]
        R.eq("%s: described, named, molecular, dated" % cat, (d, v[0], v[1], v[2]),
             (h["described"], h["C0_sum"], h["C1_sum"], h["C1d_sum"]))
        for i, key in enumerate(("CZero", "COne", "COneD")):
            R.macro(macros, "Cat%s%sPct" % (macro_name(cat), key), fmt_pct(pct(v[i], d)))

    # Largest contributor
    top = max(E, key=lambda p: (len(c0[p]), p))
    for i, (lv, nm) in enumerate(zip(LEVELS, ("CZero", "COne", "COneD", "COneB"))):
        R.macro(macros, macro_name(top) + "Share" + nm,
                fmt_pct(pct(len(level_sets(top)[i]), lad[lv]["sum"])))

    # The levels summed over groups with the largest group's tree and its
    # described species both removed (Limitations)
    den_top = as_int(bmap[top]["described"])
    for i, (lv, nm) in enumerate(zip(LEVELS[:3], ("CZero", "COne", "COneD"))):
        R.macro(macros, "Excl" + macro_name(top) + nm + "Pct",
                fmt_pct(pct(lad[lv]["sum"] - len(level_sets(top)[i]), den - den_top)))

    # Union convention and the losses between levels
    R.section("Losses between levels (each species once, in its furthest group)")
    reason = {}
    for p in E:
        if not chron[p]:
            reason[p] = "RelTime" if p in reltime else "NotChron"
        elif bounded[p]:
            reason[p] = ""
        elif p not in intervals:
            reason[p] = "NoInterval"
        elif why_of[p] == "kind":
            reason[p] = "Kind"
        else:
            reason[p] = "Calib"
    best = {}
    for p in sorted(E):
        for sp in c0[p]:
            lvl = (3 if sp in c1[p] and chron[p] and bounded[p]
                   else 2 if sp in c1[p] and chron[p] else 1 if sp in c1[p] else 0)
            if sp not in best or lvl > best[sp][0]:
                best[sp] = (lvl, p)
    by_cat = {lv: Counter() for lv in LEVELS}
    lost = Counter()
    for sp, (lvl, p) in best.items():
        for i in range(lvl + 1):
            by_cat[LEVELS[i]][bmap[p]["category"]] += 1
        if lvl == 0:
            lost["TaxPlaced"] += 1
        elif lvl < 3:
            lost[reason[p]] += 1
    cu = head["coverage_union"]
    for lv in LEVELS:
        R.eq("%s by category (union)" % lv,
             {c: by_cat[lv][c] for c in CATEGORY_ORDER},
             {c: cu["levels_by_category"][lv][c] for c in CATEGORY_ORDER})
    R.eq("losses by reason", {k: lost[k] for k in cu["lost"]}, dict(cu["lost"]))
    for k in LOSS_REASONS:
        R.macro(macros, "Lost%sN" % k, fmt_int(lost[k]))
    R.macro(macros, "LostUndatedOfCOnePct",
            fmt_pct(pct(lost["NotChron"] + lost["RelTime"], lad["C1"]["union"])))
    R.macro(macros, "LostRelTimeOfUndatedPct",
            fmt_pct(pct(lost["RelTime"], lost["NotChron"] + lost["RelTime"])))
    for k in ("NoInterval", "Kind", "Calib"):
        R.macro(macros, "Lost%sOfCOneDPct" % k, fmt_pct(pct(lost[k], lad["C1d"]["union"])))
    R.eq("species named in two or more groups",
         sum(1 for v in Counter(sp for p in E for sp in c0[p]).values() if v > 1),
         cu["species_in_several_groups"])

    # ------------------------------------------------------ sensitivities ---
    R.section("Sensitivity analyses")
    s = head["sensitivities"]
    # Denominator restricted to groups with a reference tree
    rep = ladder(E, den_rep)
    for lv, nm in zip(LEVELS, ("CZero", "COne", "COneD", "COneB")):
        R.macro(macros, "SensRep%sPct" % nm, fmt_pct(rep[lv]["pct_union"]))
    # Catalogue of Life single-authority denominators
    coldn = {r["partition"]: r for r in read_csv(
        path("registers", "col_denominators.csv"), comment=True)}
    den_col = 0
    for r in euk:
        c = as_int(coldn[r["partition"]]["col_filled"])
        den_col += c if c > 0 else as_int(r["described"])
    R.macro(macros, "SensColDenominator", fmt_int(den_col))
    col = ladder(E, den_col)
    for lv, nm in zip(LEVELS, ("CZero", "COne", "COneD", "COneB")):
        R.macro(macros, "SensCol%sPct" % nm, fmt_pct(col[lv]["pct_union"]))
    # Molecular coverage as named coverage times the molecular fraction
    prod = {lv: 0.0 for lv in ("C1", "C1d")}
    for p in E:
        x = len(c0[p]) * float("%.6f" % fmol[p])
        prod["C1"] += x
        prod["C1d"] += x if chron[p] else 0.0
    R.macro(macros, "SensProductCOnePct", fmt_pct(pct(round(prod["C1"], 4), den)))
    R.macro(macros, "SensProductCOneDPct", fmt_pct(pct(round(prod["C1d"], 4), den)))
    # Every archived interval counted; calibration-set intervals counted
    alli = ladder(E, den, {p: p in intervals for p in c0})
    R.macro(macros, "SensAllIntCOneBN", fmt_int(alli["C1b"]["union"]))
    R.macro(macros, "SensAllIntCOneBPct", fmt_pct(alli["C1b"]["pct_union"]))
    calset = {p: counts.get(p, False) or why_of.get(p) in
              ("degenerate", "at_bound", "prior_reproduced") for p in c0}
    cal = ladder(E, den, calset)
    R.macro(macros, "SensCalibCOneBN", fmt_int(cal["C1b"]["union"]))
    R.macro(macros, "SensCalibCOneBPct", fmt_pct(cal["C1b"]["pct_union"]))
    # Stricter classifications: an interval at least as wide as the range of a
    # two-sided prior on its node is set by the calibration; intervals whose
    # calibration was not checked are left out as well
    WIDTH_RULE_RATIO = 1.0
    wide = sorted(p for p in c0 if bounded[p] and intervals[p]["calibration"] == "two_sided"
                  and float(intervals[p]["prior_width_ratio"]) >= WIDTH_RULE_RATIO)
    unch = sorted(p for p in c0 if bounded[p]
                  and intervals[p]["calibration"] == "not_checked")
    strict_sets = {}
    for nm, drop in (("Width", set(wide)), ("Strict", set(wide) | set(unch))):
        lad_ = ladder(E, den, {p: bounded[p] and p not in drop for p in c0})
        strict_sets[nm] = lad_
        R.macro(macros, "Sens%sCOneBN" % nm, fmt_int(lad_["C1b"]["union"]))
        R.macro(macros, "Sens%sCOneBPct" % nm, fmt_pct(lad_["C1b"]["pct_union"]))
    R.macro(macros, "SensWidthNames", name_list(wide))
    R.macro(macros, "SensUncheckedNames", name_list([p for p in unch if p not in wide]))
    R.macro(macros, "RuleWidthRatio", "%g" % WIDTH_RULE_RATIO)
    R.macro(macros, "RuleTipsOverNamedPct", "%.0f" % (100 * TIPS_OVER_NAMED))
    R.macro(macros, "RuleElasticityPct", "%.0f" % (100 * 0.2))
    for lv, nm in (("C0", "CZero"), ("C1", "COne")):
        R.macro(macros, "SensAllMinus%sPct" % nm,
                fmt_pct(pct(lad[lv]["union"], den * (1 - 0.2))))
        R.macro(macros, "SensAllPlus%sPct" % nm,
                fmt_pct(pct(lad[lv]["union"], den * (1 + 0.2))))
    # Without the automated-synthesis tree
    prov = {r["partition_name"]: r for r in read_csv(path("registers", "tree_provenance.csv"))
            if r["is_partition_canonical"] in ("yes", "TRUE", "True", "true")}

    def assembly(p):
        r = prov.get(p, {})
        if (r.get("otl_pipeline_involvement") or "").strip() == "full_pipeline":
            return "automated_synthesis"
        if (r.get("inference") or "").strip() in ("synthesis", "literature_compilation",
                                                 "supertree_re_dated"):
            return "compiled"
        return "primary_analysis"
    R.eq("topology assembly classes", {p: assembly(p) for p in E},
         {p: bmap[p]["topology_assembly"] for p in E})
    nob = ladder([p for p in E if assembly(p) != "automated_synthesis"], den)
    R.macro(macros, "SensNoBirdsCZeroPct", fmt_pct(nob["C0"]["pct_union"]))
    R.macro(macros, "SensNoBirdsCOnePct", fmt_pct(nob["C1"]["pct_union"]))
    # GBIF backbone for the tips the Catalogue of Life leaves unmatched
    gb = {(r["query"], r["kingdom"]): r for r in read_csv(
        path("registers", "external", "gbif_unmatched.csv"))}
    added = set()
    gbif_queried = gbif_resolved = 0
    for p in E:
        rows = col_rows[p]
        kc = Counter(x["col_kingdom"] for x in rows.values() if x["col_kingdom"])
        king = kc.most_common(1)[0][0] if kc else ""
        have = {" ".join(x["col_name"].split()[:2]) for x in rows.values()
                if x["match_type"] in COL_RESOLVED and x["col_id"]}
        for x in rows.values():
            if x["match_type"] != "unmatched" or not x["query"]:
                continue
            g = gb.get((x["query"], king))
            gbif_queried += 1
            if g is None:
                R.eq("%s: GBIF row for unmatched query %s" % (p, x["query"]), None, "row")
                continue
            if (g["matchType"] in GBIF_MATCH and g["rank"] == "SPECIES"
                    and g["status"] in GBIF_STATUS and g["species"]):
                gbif_resolved += 1
                if g["species"] not in have:
                    added.add((p, g["species"]))
    R.macro(macros, "SensGbifSpeciesAdded", fmt_int(len(added)))
    R.macro(macros, "SensGbifCZeroPct", fmt_pct(pct(lad["C0"]["sum"] + len(added), den)))
    # Union of all redistributed eukaryotic trees
    fam = defaultdict(Counter)
    for p in E:
        for x in col_rows[p].values():
            if x["match_type"] in COL_RESOLVED and x["col_family"]:
                fam[x["col_family"]][p] += 1
    fam_part = {f: sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))[0][0]
                for f, c in fam.items()}
    counted = set().union(*(c0[p] for p in E))
    nc_rows = defaultdict(list)
    for x in read_csv(path("standardized", "col_resolution_noncanonical.csv")):
        nc_rows[x["filename"]].append(x)
    prok_names = {r["partition"] for r in prok}
    noncanon = sorted(f for f, d in truth_dated.items() if d["is_canonical"] != "yes"
                      and d["partition_name"].strip() not in prok_names)
    union_added = defaultdict(set)
    for f in noncanon:
        res = [x for x in nc_rows[f] if x["match_type"] in COL_RESOLVED and x["col_id"]]
        tag = truth_dated[f]["partition_name"].strip()
        if tag in E:
            tree_part = tag
        else:
            cc = Counter(fam_part[x["col_family"]] for x in res if x["col_family"] in fam_part)
            tree_part = sorted(cc.items(), key=lambda kv: (-kv[1], kv[0]))[0][0] if cc else ""
        for x in res:
            part = fam_part.get(x["col_family"]) or tree_part
            if part and x["col_id"] not in counted:
                union_added[part].add(x["col_id"])
    n_add = len(set().union(*union_added.values())) if union_added else 0
    R.macro(macros, "SensUnionNTrees", fmt_int(len(noncanon)))
    R.macro(macros, "SensUnionNamesAdded", fmt_int(n_add))
    R.macro(macros, "SensUnionCZeroPct", fmt_pct(pct(lad["C0"]["sum"] + n_add, den)))
    R.macro(macros, "SensUnionCOneHighPct", fmt_pct(pct(lad["C1"]["sum"] + n_add, den)))

    # ---------------------------------------------- source studies, dating ---
    R.section("Source studies: molecular-species trees and what each study dated")
    taxo = sorted(p for p in E if c1_basis[p] == "enumerated"
                  and len(mol_rows[p]) < len(tips[p]))
    molonly = sorted(r["partition"] for r in read_csv(
        path("registers", "molecular_only_trees.csv")))
    R.eq("groups with taxonomically placed tips = groups listed in "
         "registers/molecular_only_trees.csv", taxo, molonly)
    R.macro(macros, "NGroupsTaxonomicTips", fmt_int(len(taxo)))
    R.macro(macros, "NMolOnlyTrees", fmt_int(len(molonly)))
    sd = {r["partition"]: r["study_dating"] for r in read_csv(
        path("registers", "dating_deposit_checks.csv"))}
    for key, nm in (("absolute_deposited_fewer_tips", "FewerTips"),
                    ("absolute_not_deposited", "NotDeposited"),
                    ("absolute_not_read", "NotRead"),
                    ("relative_only", "Relative"), ("not_dated", "NotDated")):
        R.macro(macros, "StudyDating%sN" % nm,
                fmt_int(sum(1 for p in E if sd.get(p) == key and not chron[p])))

    # ------------------------------------------------------- crown ages ---
    R.section("Crown ages and intervals (registers/crown_intervals.csv)")
    crown_csv = {r["partition"]: r for r in read_csv(path("paper", "crown_ages.csv"))}
    rows = [p for p in sorted(E) if chron[p]]
    R.eq("absolute chronograms with a crown-age row", rows, sorted(crown_csv))
    for p in rows:
        f = bmap[p]["canonical_tree"]
        h = facts[f]["root"]
        iv = intervals.get(p)
        if iv and iv["kind"] in CROWN_CENTER_KINDS and iv["interval_center_ma"]:
            age = float(iv["interval_center_ma"])
        elif iv and iv["kind"] in CROWN_SUMMARY_KINDS:
            age = float(iv["source_point_height_ma"])
        else:
            age = h
        if iv and iv["kind"] in CROWN_SUMMARY_KINDS:
            sp = float(iv["source_point_height_ma"])
            if as_int(iv["n_focal_outgroups"]) == 0:
                R.eq("%s: focal node of the source summary tree is the reference "
                     "tree's root (Myr)" % p, abs(sp - h) <= CROWN_POINT_TOL_MYR, True)
            else:
                sr = float(iv["source_root_height_ma"])
                R.eq("%s: source root = reference root, focal node below it" % p,
                     abs(sr - h) <= CROWN_POINT_TOL_MYR and sp < h, True)
        n = macro_name(p)
        R.eq("%s: crown age and root height (Ma)" % p, ("%.4f" % age, "%.4f" % h),
             (crown_csv[p]["crown_age_ma"], crown_csv[p]["geometry_root_ma"]))
        R.macro(macros, "CrownAge" + n, "%.1f" % age)
        if iv:
            lo, hi = float(iv["interval_lo_ma"]), float(iv["interval_hi_ma"])
            R.macro(macros, "CrownLo" + n, "%.1f" % lo)
            R.macro(macros, "CrownHi" + n, "%.1f" % hi)
            R.macro(macros, "CrownWidth" + n, "%.1f" % (hi - lo))
            R.macro(macros, "CrownWidthPct" + n, "%.1f" % (100.0 * (hi - lo) / age))
    R.macro(macros, "CrownNChron", fmt_int(len(rows)))
    R.macro(macros, "CrownNCounted", fmt_int(sum(bounded[p] for p in rows)))
    R.macro(macros, "CrownNamesCounted", name_list(sorted(p for p in rows if bounded[p])))
    excl = Counter(why_of[p] for p in rows if p in intervals and why_of[p])
    for k, v in excl.items():
        R.macro(macros, "CrownExcl" + macro_name(k), fmt_int(v))
    R.macro(macros, "CrownNWithInterval", fmt_int(sum(p in intervals for p in rows)))
    R.macro(macros, "CrownNNoInterval", fmt_int(sum(p not in intervals for p in rows)))
    cset = sorted(p for p in rows if bounded[p] and intervals[p]["kind"] in SET_KINDS)
    csum = sorted(p for p in rows if bounded[p] and intervals[p]["kind"] not in SET_KINDS)
    oset = sorted(p for p in rows if p in intervals and not bounded[p]
                  and intervals[p]["kind"] in SET_KINDS)
    for nm, ps in (("CountedSet", cset), ("CountedSummary", csum),
                   ("SetNotCounted", oset)):
        R.macro(macros, "CrownN" + nm, fmt_int(len(ps)))
        R.macro(macros, "CrownNames" + nm, name_list(ps))

    # ------------------------------------------------ Open Tree of Life ---
    R.section("Open Tree of Life synthesis 16.1: species placed by a phylogeny")
    ext = path("registers", "external")
    cnt = {r["clade"]: r for r in read_csv(os.path.join(ext, "opentree16.1_counts.csv"))}
    ot_keys = defaultdict(set)
    with gzip.open(os.path.join(ext, "opentree16.1_phylo_species.tsv.gz"), "rt",
                   encoding="utf-8") as fh:
        next(fh)
        for ln in fh:
            _u, nm, clade = ln.rstrip("\n").split("\t")
            k, ok = species_key(nm.replace(" ", "_"))
            if ok:
                ot_keys[clade].add(k)
    ot_acc = defaultdict(set)
    for r in read_csv(os.path.join(ext, "opentree16.1_phylo_species_col.csv")):
        if r["match_type"] in COL_RESOLVED and r["col_id"]:
            ot_acc[r["clade"]].add(r["col_id"])
    ot_rows = []
    for p in E:
        if p not in cnt or not c0[p]:
            continue
        rows = col_rows[p]
        q = lambda nm: rows[tips[p][nm]]["query"] or None
        akeys = {k for k in map(q, tips[p]) if k}
        o, oa = ot_keys.get(p, set()), ot_acc.get(p, set())
        row = dict(partition=p, c1=len(c1[p]),
                   otp=as_int(cnt[p]["ot_phylo_extant_species"]),
                   both=len(akeys & o), enum=c1_basis[p] == "enumerated")
        if row["enum"]:
            mnames = {dictionary.get(sid) for sid in mol_rows[p]} & set(tips[p])
            mkeys = {k for k in map(q, mnames) if k}
            macc = {unit_of[p][sid] for sid in mol_rows[p] if sid in unit_of[p]}
            row.update(mkeys=len(mkeys), mand=len(mkeys & o), acc=len(macc),
                       accand=len(macc & oa), accnot=len(macc - oa))
        ot_rows.append(row)
    euk_ot = as_int(cnt["__Eukaryota"]["ot_phylo_extant_species"])
    tot_c1 = sum(r["c1"] for r in ot_rows)
    tot_ot = sum(r["otp"] for r in ot_rows)
    enum = [r for r in ot_rows if r["enum"]]
    accnot = sum(r["accnot"] for r in enum)
    R.macro(macros, "OTPhyloEukN", fmt_int(euk_ot))
    R.macro(macros, "OTPhyloEukPct", fmt_pct(pct(euk_ot, den)))
    R.macro(macros, "OTCOneRatio", "%.2f" % round(tot_c1 / float(tot_ot), 4))
    R.macro(macros, "OTKeysBoth", fmt_int(sum(r["both"] for r in ot_rows)))
    R.macro(macros, "OTKeysMolecular", fmt_int(sum(r["mkeys"] for r in enum)))
    R.macro(macros, "OTKeysMolecularAndOT", fmt_int(sum(r["mand"] for r in enum)))
    R.macro(macros, "OTAccMolecular", fmt_int(sum(r["acc"] for r in enum)))
    R.macro(macros, "OTAccMolecularAndOT", fmt_int(sum(r["accand"] for r in enum)))
    R.macro(macros, "OTAccMolecularNotOT", fmt_int(accnot))
    R.macro(macros, "OTAccMolecularAndOTPct", fmt_pct(pct(
        sum(r["accand"] for r in enum), sum(r["acc"] for r in enum))))
    R.macro(macros, "OTNPartitionsCompared", fmt_int(len(enum)))
    for r in enum:
        R.macro(macros, "OTNotOTSharePct" + macro_name(r["partition"]),
                fmt_pct(pct(r["accnot"], accnot)))

    # --------------------------------------------------- worked example ---
    R.section("Worked example: eligible groups and the numbers printed from the fits")
    # The fits are made by code/rebuild/worked_example_bd.R (R, ape), which
    # code/rebuild/run_worked_example.py reruns from this repository and
    # compares with paper/worked_example_bd*.csv.  Here: the eligible set and
    # every input of the fits are recomputed, and every printed number is
    # formed from the committed fits by the rules the text states.
    WE_LL_TOL, WE_RATIO_RTOL = 1e-3, 1e-4      # same log-likelihood; ratio
    WE_LR_CRIT = LR_CRITICAL
    elig = [p for p in sorted(E) if chron[p] and c1_basis[p] == "enumerated"
            and len(mol_rows[p]) < len(tips[p])]
    sel = read_csv(os.path.join(a.fits, "worked_example_bd_selection.csv"))
    R.eq("eligible groups (absolute chronogram, molecular tips listed, fewer "
         "than all tips)", elig, [x["partition"] for x in sel])
    for x in sel:
        p = x["partition"]
        R.eq("%s: tips, molecular tips, named, described, crown age, root height" % p,
             (len(tips[p]), len(mol_rows[p]), len(c0[p]), as_int(bmap[p]["described"]),
              "%.4f" % float(crown_csv[p]["crown_age_ma"]), "%.4f" % facts[bmap[p]["canonical_tree"]]["root"]),
             (as_int(x["labels_unique"]), as_int(x["n_molecular"]), as_int(x["c0_units"]),
              as_int(x["described"]), "%.4f" % float(x["crown_age_ma"]),
              "%.4f" % float(x["geometry_root_ma"])))
    fits = {(f["partition"], f["tree"], f["rho_convention"]): f
            for f in read_csv(os.path.join(a.fits, "worked_example_bd.csv"))}

    def g(p, tree, conv, k):
        return float(fits[(p, tree, conv)][k])

    def overlap(p, a_, b_, v):
        return (g(p, a_[0], a_[1], v + "_lo") <= g(p, b_[0], b_[1], v + "_hi")
                and g(p, b_[0], b_[1], v + "_lo") <= g(p, a_[0], a_[1], v + "_hi"))
    M, MC, F, FM = (("molecular", "tips"), ("molecular", "full_misapplied"),
                    ("full", "tips"), ("molecular", "family"))
    n_ident = n_bound = n_rej = n_fg = n_fs = n_fo = n_famo = n_el = n_eh = 0
    for p in elig:
        n = macro_name(p)
        rho_t, rho_c = g(p, *M, k="rho"), g(p, *MC, k="rho")
        R.eq("%s: rho of the two same-tree fits = molecular tips, full-tree tips, "
             "over described" % p, ("%.6f" % rho_t, "%.6f" % rho_c),
             ("%.6f" % (len(mol_rows[p]) / float(bmap[p]["described"])),
              "%.6f" % (len(tips[p]) / float(bmap[p]["described"]))))
        lam_t, lam_c = g(p, *M, k="lambda"), g(p, *MC, k="lambda")
        ll_t, ll_c = g(p, *M, k="loglik"), g(p, *MC, k="loglik")
        pred, ratio, dll = rho_t / rho_c, lam_c / lam_t, ll_t - ll_c
        ident = abs(dll) <= WE_LL_TOL and abs(ratio - pred) <= WE_RATIO_RTOL * pred
        lr = round(2.0 * dll, 4)
        n_ident += ident
        n_bound += not ident
        n_rej += (not ident) and lr > WE_LR_CRIT
        rf, rm = g(p, *F, k="r"), g(p, *M, k="r")
        rr = round(rf / rm, 4)
        n_fg += rr > 1
        n_fs += rr < 1
        n_fo += overlap(p, F, M, "r")
        n_famo += overlap(p, M, FM, "epsilon")
        # family and uniform sampling give the same tips different sampling
        # probabilities, so their likelihoods are not compared: the text
        # compares the extinction fractions only
        n_el += g(p, *FM, k="epsilon") < g(p, *M, k="epsilon")
        n_eh += g(p, *FM, k="epsilon") > g(p, *M, k="epsilon")
        R.macro(macros, "WERhoTips" + n, "%.2f" % rho_t)
        R.macro(macros, "WERhoFull" + n, "%.2f" % rho_c)
        R.macro(macros, "WEEpsTips" + n, "%.2f" % g(p, *M, k="epsilon"))
        R.macro(macros, "WEEpsFull" + n, "%.2f" % g(p, *MC, k="epsilon"))
        R.macro(macros, "WELambdaRatio" + n, "%.2f" % round(ratio, 4))
        R.macro(macros, "WELR" + n, "%.2f" % lr)
        R.macro(macros, "WERFull" + n, "%.3f" % rf)
        R.macro(macros, "WERMol" + n, "%.3f" % rm)
        R.macro(macros, "WERRatio" + n, "%.2f" % rr)
        R.macro(macros, "WEEpsFamily" + n, "%.2f" % g(p, *FM, k="epsilon"))
        R.macro(macros, "WERFamily" + n, "%.3f" % g(p, *FM, k="r"))
    R.macro(macros, "WENEligible", fmt_int(len(elig)))
    R.macro(macros, "WENIdentity", fmt_int(n_ident))
    R.macro(macros, "WENBoundary", fmt_int(n_bound))
    R.macro(macros, "WENRejected", fmt_int(n_rej))
    R.macro(macros, "WENFullGreater", fmt_int(n_fg))
    R.macro(macros, "WENFullSmaller", fmt_int(n_fs))
    R.macro(macros, "WENFullOverlap", fmt_int(n_fo))
    R.macro(macros, "WENFamilyOverlap", fmt_int(n_famo))
    R.macro(macros, "WENFamilyEpsLower", fmt_int(n_el))
    R.macro(macros, "WENFamilyEpsHigher", fmt_int(n_eh))
    R.macro(macros, "WENamesEligible", name_list(elig))
    # Polytomies of the fitted trees and the fits with them resolved at random
    poly = defaultdict(lambda: defaultdict(dict))
    for r in read_csv(os.path.join(a.fits, "worked_example_bd_polytomy.csv")):
        poly[r["partition"]][r["seed"]][(r["tree"], r["rho_convention"])] = r
    sfac = {x["partition"]: float(x["time_scale_factor"]) for x in sel}
    with_poly, allr, nb, nrj = [], [], 0, 0
    poly_rows = OrderedDict()
    for p in elig:
        with open(os.path.join(tree_dir, bmap[p]["canonical_tree"]), encoding="utf-8",
                  errors="replace") as fh:
            text_ = fh.read()
        for tree, keep in (("full", None), ("molecular", mol_rows[p])):
            ne, npe = polytomy_counts(text_, keep, sfac[p])
            for conv in (("tips", "full_misapplied", "family") if tree == "molecular"
                         else ("tips",)):
                R.eq("%s %s/%s: branching events and those at polytomies, counted "
                     "on the tree" % (p, tree, conv), (ne, npe),
                     (as_int(fits[(p, tree, conv)]["n_events"]),
                      as_int(fits[(p, tree, conv)]["n_polytomy_events"])))
        ff_, fm_ = fits[(p,) + F], fits[(p,) + M]
        if as_int(ff_["n_polytomy_events"]) or as_int(fm_["n_polytomy_events"]):
            with_poly.append(p)
        rr, lrs, em, bnd, rej = [], [], [], [], []
        for sd_, f_ in sorted(poly[p].items()):
            R.eq("%s seed %s: the three uniform fits resolved" % (p, sd_), len(f_), 3)
            x_ = round(2.0 * (float(f_[M]["loglik"]) - float(f_[MC]["loglik"])), 4)
            e_ = float(f_[MC]["epsilon"])
            rr.append(round(float(f_[F]["r"]) / float(f_[M]["r"]), 4))
            lrs.append(x_)
            em.append(float(f_[M]["epsilon"]))
            bnd.append(e_ < 1e-6 and x_ > 2e-3)
            rej.append(e_ < 1e-6 and x_ > WE_LR_CRIT)
        allr += rr
        nb += all(bnd)
        nrj += all(rej)
        n = macro_name(p)
        if p in with_poly:
            R.macro(macros, "WEPolyPctFull" + n, fmt_pct(pct(
                as_int(ff_["n_polytomy_events"]), as_int(ff_["n_events"]))))
            R.macro(macros, "WEPolyPctPruned" + n, fmt_pct(pct(
                as_int(fm_["n_polytomy_events"]), as_int(fm_["n_events"]))))
            ll_t, ll_c = g(p, *M, k="loglik"), g(p, *MC, k="loglik")
            if g(p, *MC, k="epsilon") < 1e-6 and abs(ll_t - ll_c) > WE_LL_TOL:
                R.macro(macros, "WEPolyLRMin" + n, "%.1f" % min(lrs))
                R.macro(macros, "WEPolyLRMax" + n, "%.1f" % max(lrs))
            R.macro(macros, "WEPolyRRatioMin" + n, "%.2f" % min(rr))
            R.macro(macros, "WEPolyRRatioMax" + n, "%.2f" % max(rr))

        def rng(v, f="%.2f"):
            lo_, hi_ = f % min(v), f % max(v)
            return lo_ if lo_ == hi_ else "%s--%s" % (lo_, hi_)
        fc_ = fits[(p,) + MC]
        poly_rows[p] = OrderedDict([
            ("partition", p), ("events_full", as_int(ff_["n_events"])),
            ("polytomy_events_full", as_int(ff_["n_polytomy_events"])),
            ("events_pruned", as_int(fm_["n_events"])),
            ("polytomy_events_pruned", as_int(fm_["n_polytomy_events"])),
            ("r_ratio", "%.2f" % round(float(ff_["r"]) / float(fm_["r"]), 4)),
            ("r_ratio_resolved", rng(rr)),
            ("epsilon_matching", "%.2f" % float(fm_["epsilon"])),
            ("epsilon_matching_resolved", rng(em)),
            ("lr", "%.2f" % round(2.0 * (float(fm_["loglik"]) - float(fc_["loglik"])), 4)),
            ("lr_resolved", rng(lrs)),
            ("outcome_resolved", "mu0_rejected" if all(rej) else "mu0" if all(bnd)
             else "rescaled" if not any(bnd) else "mixed")])
    R.macro(macros, "WEPolyNGroups", fmt_int(len(with_poly)))
    R.macro(macros, "WEPolyNames", name_list(with_poly))
    R.macro(macros, "WEPolyRRatioMin", "%.2f" % min(allr))
    R.macro(macros, "WEPolyRRatioMax", "%.2f" % max(allr))
    R.macro(macros, "WEPolyNBoundary", fmt_int(nb))
    R.macro(macros, "WEPolyNRejected", fmt_int(nrj))
    sf = [float(x["time_scale_factor"]) for x in sel]
    R.macro(macros, "WEScaleMin", "%.2f" % min(sf))
    R.macro(macros, "WEScaleMax", "%.2f" % max(sf))
    for x in sel:
        p = x["partition"]
        R.eq("%s: time-scale factor = crown age / root height" % p,
             "%.6f" % (float(crown_csv[p]["crown_age_ma"]) / facts[bmap[p]["canonical_tree"]]["root"]),
             "%.6f" % float(x["time_scale_factor"]))

    # ------------------------------------------------- main-text tables ---
    R.section("Main-text tables (paper/table_*.csv), cell by cell")

    def rho_fmt(x):
        return "%.2f" % x if x >= 0.1 else "%#.2g" % x

    def table_eq(name, got, key):
        with open(path("paper", name), newline="", encoding="utf-8") as fh:
            want = list(csv.DictReader(fh))
        want = {w[key]: w for w in want}
        R.eq("%s: rows" % name, sorted(got), sorted(want))
        for k, row in got.items():
            w = want.get(k, {})
            bad = {c: (v, w.get(c)) for c, v in row.items() if str(v) != w.get(c)}
            R.eq("%s: %s" % (name, k), bad, {})

    reg_by_file = {g["filename"]: g for g in reg}
    t1 = OrderedDict()
    for p in E:
        b = bmap[p]
        d = as_int(b["described"])
        mt = len(mol_rows[p]) if c1_basis[p] == "enumerated" else len(tips[p])
        iv = intervals.get(p) if chron[p] else None
        kind = iv["kind"] if iv else ""
        ts = ("absolute" if chron[p] else "relative" if p in reltime
              else "not_ultrametric")
        sdp = sd.get(p) or ("absolute_deposited" if chron[p] else
                            "relative_only" if p in reltime else "")
        t1[p] = OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("bibkey", reg_by_file[b["canonical_tree"]]["bibkey"]),
            ("flags", reg_by_file[b["canonical_tree"]]["flags"]),
            ("tips", len(tips[p])), ("named_pct", fmt_pct(pct(len(c0[p]), d))),
            ("molecular_tips", mt), ("molecular_basis", c1_basis[p]),
            ("molecular_pct", fmt_pct(pct(len(c1[p]), d))), ("time_scale", ts),
            ("crown_age_ma", "%.1f" % float(crown_csv[p]["crown_age_ma"]) if chron[p] else ""),
            ("interval_lo_ma", "%.1f" % float(iv["interval_lo_ma"]) if kind else ""),
            ("interval_hi_ma", "%.1f" % float(iv["interval_hi_ma"]) if kind else ""),
            ("interval_kind", kind),
            ("interval_status", "" if not chron[p] else "none" if not kind
             else "counted" if counts[p] else why_of[p]),
            ("rho_full", rho_fmt(len(tips[p]) / float(d))),
            ("rho_molecular", rho_fmt(mt / float(d))), ("described", d),
            ("study_dating", sdp),
            ("molecular_only_tree", "yes" if p in molonly else "no"),
            ("molecular_share", "%.2f" % (mt / float(len(tips[p])))),
            ("tips_exceed_named", "yes" if len(tips[p]) > (1 + TIPS_OVER_NAMED)
             * len(c0[p]) else "no"),
        ])
    table_eq("table_which_tree.csv", t1, "partition")

    t2 = OrderedDict()
    for lv, label in zip(LEVELS, ("named", "molecular", "dated", "interval")):
        t2[label] = OrderedDict([
            ("level", label), ("species", fmt_int(lad[lv]["union"])),
            ("pct", fmt_pct(lad[lv]["pct_union"])),
            ("rep", fmt_pct(rep[lv]["pct_union"])), ("col", fmt_pct(col[lv]["pct_union"])),
            ("union", fmt_pct(pct(lad["C0"]["sum"] + n_add, den)) if lv == "C0"
             else "<= " + fmt_pct(pct(lad["C1"]["sum"] + n_add, den)) if lv == "C1" else ""),
            ("gbif", fmt_pct(pct(lad["C0"]["sum"] + len(added), den)) if lv == "C0" else ""),
            ("product", fmt_pct(pct(round(prod[lv], 4), den)) if lv in prod else ""),
            ("calib", fmt_pct(cal["C1b"]["pct_union"]) if lv == "C1b" else ""),
            ("allint", fmt_pct(alli["C1b"]["pct_union"]) if lv == "C1b" else ""),
            ("nobirds", fmt_pct(nob[lv]["pct_union"])),
        ])
    table_eq("table_rungs_sensitivity.csv", t2, "level")

    t3 = OrderedDict()
    for p in elig:
        ft, fc, ff, fm = (fits[(p,) + M], fits[(p,) + MC], fits[(p,) + F],
                          fits[(p,) + FM])
        f2 = lambda x: "%.2f" % float(x)
        f3 = lambda x: "%.3f" % float(x)
        lam = round(float(fc["lambda"]) / float(ft["lambda"]), 4)
        dll = float(ft["loglik"]) - float(fc["loglik"])
        ident = (abs(dll) <= WE_LL_TOL and abs(float(fc["lambda"]) / float(ft["lambda"])
                 - float(ft["rho"]) / float(fc["rho"]))
                 <= WE_RATIO_RTOL * float(ft["rho"]) / float(fc["rho"]))
        lr = round(2.0 * dll, 4)
        t3[p] = OrderedDict([
            ("partition", p), ("tips_full", as_int(ff["n_tips"])),
            ("tips_molecular", as_int(ft["n_tips"])),
            ("rho_molecular", f2(ft["rho"])), ("rho_fulltree", f2(fc["rho"])),
            ("lambda_ratio", "%.2f" % lam),
            ("eps_molecular", f2(ft["epsilon"])), ("eps_molecular_lo", f2(ft["epsilon_lo"])),
            ("eps_molecular_hi", f2(ft["epsilon_hi"])),
            ("eps_fullrho", f2(fc["epsilon"])), ("eps_fullrho_lo", f2(fc["epsilon_lo"])),
            ("eps_fullrho_hi", f2(fc["epsilon_hi"])), ("lr", "%.2f" % lr),
            ("classification", ("identity" if ident else "boundary")
             + (" (rejected)" if (not ident) and lr > WE_LR_CRIT else "")),
            ("r_full", f3(ff["r"])), ("r_full_lo", f3(ff["r_lo"])), ("r_full_hi", f3(ff["r_hi"])),
            ("r_molecular", f3(ft["r"])), ("r_molecular_lo", f3(ft["r_lo"])),
            ("r_molecular_hi", f3(ft["r_hi"])),
            ("eps_family", f2(fm["epsilon"])), ("eps_family_lo", f2(fm["epsilon_lo"])),
            ("eps_family_hi", f2(fm["epsilon_hi"])),
            ("r_family", f3(fm["r"])), ("r_family_lo", f3(fm["r_lo"])),
            ("r_family_hi", f3(fm["r_hi"])),
        ])
    table_eq("table_worked_example.csv", t3, "partition")

    R.section("Figure source data (paper/fig*_*.csv)")
    casc = read_csv(path("paper", "fig1_coverage_cascade.csv"))
    got = {("level", r["level"], r["group"]): by_cat[r["level"]][r["group"]]
           for r in casc if r["kind"] == "level"}
    got.update({("loss", r["level"], r["group"]): lost[r["group"]]
                for r in casc if r["kind"] == "loss"})
    R.eq("Figure 1: levels by category and losses by reason",
         {k: str(v) for k, v in got.items()},
         {(r["kind"], r["level"], r["group"]): r["species"] for r in casc})
    for fig in ("fig3_c0c1_bars.csv", "fig2_coverage_scatter.csv"):
        rows = read_csv(path("paper", fig))
        mine = {p: (str(as_int(bmap[p]["described"])), str(len(c0[p])),
                    "%.4f" % pct(len(c0[p]), as_int(bmap[p]["described"])),
                    "%.4f" % pct(len(c1[p]), as_int(bmap[p]["described"])),
                    "yes" if chron[p] else "no") for p in E}
        theirs = {r["partition"]: (r["described"], r["c0_units"], r["coverage_C0_pct"],
                                   r["coverage_C1_pct"], r["is_chronogram"]) for r in rows}
        R.eq("%s: described, named, named %%, molecular %%, chronogram" % fig,
             {k: v for k, v in mine.items() if theirs.get(k) != v}, {})
        if "f_mol" in rows[0]:
            R.eq("%s: molecular fraction" % fig,
                 {r["partition"]: r["f_mol"] for r in rows
                  if r["f_mol"] != "%.6f" % fmol[r["partition"]]}, {})

    # ------------------------------------------------------- Supplement ---
    # Supplement S1-S10: every table (paper/supp_*.csv) cell by cell and
    # every number its text prints (paper/supp_macros.tex).
    smac = read_macros(path("paper", "supp_macros.tex"))
    R.supp_macros = smac

    def supp_rows(name):
        return read_csv(path("paper", name + ".csv"))

    def supp_eq(name, got, key):
        want = OrderedDict((w[key], w) for w in supp_rows(name))
        R.eq("%s: rows, in order" % name, list(got), list(want))
        for k, row in got.items():
            w = want.get(k, {})
            bad = {c: (str(v), w.get(c)) for c, v in row.items()
                   if str(v) != w.get(c)}
            R.eq("%s: %s" % (name, k), bad, {})
        # the .tex table the Supplement prints carries every count of the csv
        tex = open(path("paper", name + ".tex"), encoding="utf-8").read()
        miss = []
        for k, row in want.items():
            for c, v in row.items():
                if re.fullmatch(r"\d{4,}", v or "") and fmt_int(int(v)) not in tex:
                    miss.append((k, c, v))
        R.eq("%s.tex prints every count of %s.csv" % (name, name), miss, [])

    def by_cat(parts):
        return sorted(parts, key=lambda p: (CATEGORY_ORDER.index(bmap[p]["category"]), p))

    def sm(name, got):
        R.macro(smac, name, got)

    def fix(x, dp):
        return "{0:.{1}f}".format(float(x), dp)

    # S1: reference trees and the rule that chose each
    R.section("Supplement S1: reference trees and their sources")
    rts = {r["partition"]: r for r in read_csv(path("registers", "reference_tree_sources.csv"))}
    prep = {r["filename"]: r for r in read_csv(path("registers", "canonical_source_prep.csv"))}
    prep_sel = {"D8": "study_chronogram", "D15": "more_branch_lengths",
                "D16": "first_of_set", "D17": "mcc_of_set"}
    got = OrderedDict()
    for p in by_cat(E):
        s_, fn = rts[p], bmap[p]["canonical_tree"]
        R.eq("%s: register names the reference tree file" % p, s_["filename"], fn)
        if fn in prep:
            # the last decision of a leading 'D6/D16' is the one that chose
            d_ = re.match(r"\s*(D\d+(?:/D\d+)*)", prep[fn]["evidence"])
            dec = d_.group(1).split("/")[-1] if d_ else ""
            R.eq("%s: selection = the decision the prepared source cites" % p,
                 s_["selection"], prep_sel.get(dec, "release"))
            R.eq("%s: source file = the prepared source" % p,
                 s_["source_file"], prep[fn]["raw_source"])
        got[p] = OrderedDict([
            ("partition", p), ("category", bmap[p]["category"]), ("filename", fn),
            ("bibkey", reg_by_file[fn]["bibkey"]),
            ("source_file", os.path.basename(s_["source_file"])),
            ("tree_index", s_["tree_index"]), ("selection", s_["selection"]),
            ("tree_used", s_["tree_used"]), ("archive", s_["archive"])])
    supp_eq("supp_reference_trees", got, "partition")
    for code in ("as_archived", "release", "study_chronogram",
                 "more_branch_lengths", "first_of_set", "mcc_of_set"):
        nm = macro_name(code.replace("_", " "))
        ps = sorted(p for p in E if rts[p]["selection"] == code)
        sm("SuppSel" + nm, fmt_int(len(ps)))
        sm("SuppSelNames" + nm, name_list(ps) if ps else "none")

    # S2: tips removed and merged
    R.section("Supplement S2: tips removed outside the group and merged per species")
    got = OrderedDict()
    for p in by_cat(E):
        fn = bmap[p]["canonical_tree"]
        n_src, n_rm = src_tips[fn], len(pruned[fn])
        n_open = sum(1 for nm in tips[p] if not species_key(nm)[1])
        got[p] = OrderedDict([
            ("partition", p), ("category", bmap[p]["category"]),
            ("source_tips", n_src), ("removed_outside", n_rm),
            ("merged_same_species", n_src - n_rm - len(tips[p])),
            ("tips", len(tips[p])), ("unresolved_labels_kept", n_open)])
    supp_eq("supp_pruned", got, "partition")
    for k, col_ in (("SuppSourceTipsN", "source_tips"), ("SuppRemovedN", "removed_outside"),
                    ("SuppMergedN", "merged_same_species"), ("SuppTipsN", "tips"),
                    ("SuppOpenLabelsN", "unresolved_labels_kept")):
        sm(k, fmt_int(sum(r[col_] for r in got.values())))
    sm("SuppRemovedTrees", fmt_int(sum(1 for r in got.values() if r["removed_outside"])))
    sm("SuppMergedTrees", fmt_int(sum(1 for r in got.values() if r["merged_same_species"])))

    # S3: name resolution
    R.section("Supplement S3: name resolution per group")
    gb_add = Counter(p_ for p_, _ in added)
    got = OrderedDict()
    for p in by_cat(E):
        mts = [col_rows[p][sid]["match_type"] for sid in tips[p].values()
               if sid in col_rows[p] and col_rows[p][sid]["match_type"] in COL_RESOLVED
               and col_rows[p][sid]["col_id"]]
        syn = sum(1 for m in mts if m.startswith("synonym"))
        acc = len(mts) - syn
        n_lab = len(tips[p])
        got[p] = OrderedDict([
            ("partition", p), ("category", bmap[p]["category"]), ("labels", n_lab),
            ("accepted", acc), ("synonym", syn), ("unresolved", n_lab - len(mts)),
            ("unresolved_pct", fix(pct(n_lab - len(mts), n_lab), 1)),
            ("named_species", len(c0[p])), ("merged", len(mts) - len(c0[p])),
            ("gbif_added", gb_add[p]),
            ("override", sum(1 for sid in tips[p].values() if sid in col_rows[p]
                             and col_rows[p][sid]["match_type"] in COL_RESOLVED
                             and col_rows[p][sid]["col_id"]
                             and col_rows[p][sid].get("query_source") == "override"))])
    supp_eq("supp_names", got, "partition")
    sm("SuppSynonymLabelsN", fmt_int(sum(r["synonym"] for r in got.values())))
    sm("SuppMergedNamesN", fmt_int(sum(r["merged"] for r in got.values())))
    sm("SuppOverridesN", fmt_int(len(read_csv(path("standardized", "name_overrides.csv")))))
    sm("SuppRejectedTargetsN", fmt_int(len(read_csv(
        path("standardized", "name_rejected_targets.csv")))))
    sm("SuppSupraGeneraN", fmt_int(len(SUPRA_GENERA)))
    sm("SuppGbifQueriedN", fmt_int(gbif_queried))
    sm("SuppGbifResolvedN", fmt_int(gbif_resolved))

    # S4: molecular tips
    R.section("Supplement S4: how the molecular tips of each tree are found")
    got = OrderedDict()
    for p in by_cat(E):
        b = bmap[p]
        mt = len(mol_rows[p]) if c1_basis[p] == "enumerated" else len(tips[p])
        ex = exact_tips.get(p, "")
        got[p] = OrderedDict([
            ("partition", p), ("category", b["category"]), ("f_basis", b["f_basis"]),
            ("source_file", os.path.basename(b["molecular_source_file"] or "")),
            ("tips", len(tips[p])), ("molecular_tips", mt),
            ("molecular_fraction", "%.3f" % (mt / float(len(tips[p])))),
            ("exact_label_tips", ex), ("molecular_species", len(c1[p])),
            ("product", fix(len(c0[p]) * float("%.6f" % fmol[p]), 1)),
            ("f_low", "%.3f" % float(b["f_low"])), ("f_high", "%.3f" % float(b["f_high"])),
            ("molecular_only_tree", "yes" if p in molonly else "no")])
    supp_eq("supp_molecular", got, "partition")
    for code in ("source_tree_intersection", "source_alignment_intersection",
                 "per_tip_evidence", "asserted_no_imputation"):
        sm("SuppMolBasis" + macro_name(code.replace("_", " ")),
           fmt_int(sum(1 for r in got.values() if r["f_basis"] == code)))
    sm("SuppExactZeroNames", name_list(sorted(p for p, r in got.items()
                                              if r["exact_label_tips"] == 0)))
    banded = [p for p in E if float(bmap[p]["f_high"]) > float(bmap[p]["f_low"])]
    R.eq("groups with a molecular-fraction band", banded, ["Birds"])
    sm("SuppBandLow", "%.3f" % float(bmap[banded[0]]["f_low"]))
    sm("SuppBandHigh", "%.3f" % float(bmap[banded[0]]["f_high"]))

    # S5: time scale
    R.section("Supplement S5: time scale of each reference tree")
    got = OrderedDict()
    for p in by_cat(E):
        f = bmap[p]["canonical_tree"]
        ultra = facts[f]["cv"] <= CHRONOGRAM_CV_MAX
        ts = time_scale.get(f, "") if ultra else "not_ultrametric"
        got[p] = OrderedDict([
            ("partition", p), ("category", bmap[p]["category"]),
            ("roottotip_cv", "%.4f" % facts[f]["cv"]),
            ("ultrametric", "yes" if ultra else "no"),
            ("root_height", "%.4g" % facts[f]["root"]), ("time_scale", ts),
            ("study_dating", sd.get(p) or ("absolute_deposited" if chron[p] else
                                           "relative_only" if p in reltime else ""))])
    supp_eq("supp_time_scale", got, "partition")
    ucv = [facts[bmap[p]["canonical_tree"]]["cv"] for p in E
           if got[p]["ultrametric"] == "yes"]
    ocv = [facts[bmap[p]["canonical_tree"]]["cv"] for p in E
           if got[p]["ultrametric"] == "no"]
    sm("SuppMaxCVUltrametric", "%.4f" % max(ucv))
    sm("SuppMinCVOther", "%.4f" % min(ocv))
    sm("SuppNUltrametric", fmt_int(len(ucv)))

    # S6: crown ages and intervals
    R.section("Supplement S6: crown ages, intervals and the rules that classify them")
    clades = {r["partition"]: r["crown_clade"] for r in read_csv(
        path("registers", "crown_clades.csv"))}
    got = OrderedDict()
    for p in by_cat([p for p in E if chron[p]]):
        iv = intervals.get(p) or {}
        kind = iv.get("kind", "")
        c_ = crown_csv[p]
        status = ("none" if not kind else "counted" if counts[p] else why_of[p])
        pc, pw = iv.get("prior_containment") or "", iv.get("prior_width_ratio") or ""
        got[p] = OrderedDict([
            ("partition", p), ("category", bmap[p]["category"]),
            ("crown_clade", clades[p]),
            ("crown_age_ma", "%.1f" % float(c_["crown_age_ma"])),
            ("age_statistic", c_["age_statistic"]),
            ("interval_lo_ma", "%.1f" % float(iv["interval_lo_ma"]) if kind else ""),
            ("interval_hi_ma", "%.1f" % float(iv["interval_hi_ma"]) if kind else ""),
            ("interval_width_myr", "%.3g" % float("%.4f" % float(iv["interval_width_myr"]))
             if kind else ""),                                   # Myr, to 4 decimals
            ("interval_kind", kind), ("n_trees", iv.get("n_trees", "") if kind else ""),
            ("calibration", iv.get("calibration", "") if kind else ""),
            ("prior_lo_ma", "%.1f" % float(iv["prior_lo_ma"]) if iv.get("prior_lo_ma") else ""),
            ("prior_hi_ma", "%.1f" % float(iv["prior_hi_ma"]) if iv.get("prior_hi_ma") else ""),
            ("prior_containment_pct", "%.0f" % (100 * float(pc)) if pc else ""),
            ("prior_width_ratio", "%.2f" % float(pw) if pw else ""),
            ("status", status)])
    supp_eq("supp_crown", got, "partition")
    R.eq("crown interval widths = hi - lo of registers/crown_intervals.csv",
         [p for p in got if intervals.get(p) and abs(
             float(intervals[p]["interval_width_myr"])
             - (float(intervals[p]["interval_hi_ma"]) - float(intervals[p]["interval_lo_ma"])))
          > 1e-4], [])
    for key, nm in (("degenerate", "Degenerate"), ("at_bound", "AtBound"),
                    ("prior_reproduced", "Prior"), ("kind", "Kind")):
        ps = sorted(p for p, r in got.items() if r["status"] == key)
        sm("SuppCrown" + nm, name_list(ps) if ps else "none")
    sm("SuppCrownNone", name_list(sorted(p for p, r in got.items() if r["status"] == "none")))
    qs = rule_from_code("crown_intervals.py", "QUANTILES")
    sm("SuppCrownIntervalPct", "%.0f" % (100 * (qs[-1] - qs[0])))

    # S7: sensitivity analyses and per-group denominators
    R.section("Supplement S7: sensitivity analyses")
    ELASTIC = 0.2
    sm("SuppElasticityPct", "%.0f" % (100 * ELASTIC))
    d_ins = as_int(bmap["Insects"]["described"])
    sm("SuppInsectSharePct", fmt_pct(pct(d_ins, den)))
    R.eq("the insect count is the largest single count of the denominator",
         max(euk, key=lambda r: as_int(r["described"]))["partition"], "Insects")
    den_base, keep_base = 0, []
    for r in euk:
        cb = as_int(coldn[r["partition"]]["col_base"])
        den_base += cb if cb > 0 else as_int(r["described"])
        if cb <= 0:
            keep_base.append(r["partition"])
    sm("SuppColBaseKeepN", fmt_int(len(keep_base)))
    cls = Counter(bmap[p]["topology_assembly"] for p in E)
    sm("SuppAssemblyPrimaryN", fmt_int(cls["primary_analysis"]))
    sm("SuppAssemblyCompiledN", fmt_int(cls["compiled"]))
    sm("SuppAssemblyAutomatedN", fmt_int(cls["automated_synthesis"]))
    sm("SuppAssemblyCompiledNames", name_list(sorted(
        p for p in E if bmap[p]["topology_assembly"] == "compiled")))
    sm("SuppAssemblyAutomatedNames", name_list(sorted(
        p for p in E if bmap[p]["topology_assembly"] == "automated_synthesis")))
    sm("SuppUnionTreesN", fmt_int(len(noncanon)))

    def ladder_pcts(parts, denominator):
        l_ = ladder(parts, denominator)
        return OrderedDict((lv, pct(l_[lv]["union"], denominator)) for lv in LEVELS)

    def srow(key, d_, pcts):
        return OrderedDict([("alternative", key), ("denominator", int(round(d_)))]
                           + [(lv, "" if pcts.get(lv) is None else fmt_pct(pcts[lv]))
                              for lv in LEVELS])
    got = OrderedDict()
    got["reported"] = srow("reported", den, ladder_pcts(E, den))
    for key, sgn in (("insects_minus", -1), ("insects_plus", 1)):
        d_ = den + sgn * int(round(ELASTIC * d_ins))
        got[key] = srow(key, d_, ladder_pcts(E, d_))
    for key, sgn in (("all_minus", -1), ("all_plus", 1)):
        d_ = den * (1 + sgn * ELASTIC)
        got[key] = srow(key, d_, ladder_pcts(E, d_))
    got["col_filled"] = srow("col_filled", den_col, ladder_pcts(E, den_col))
    got["col_base"] = srow("col_base", den_base, ladder_pcts(E, den_base))
    got["rep"] = srow("rep", den_rep, ladder_pcts(E, den_rep))
    got["no_seed"] = srow("no_seed", den - den_top, OrderedDict(
        (lv, pct(lad[lv]["sum"] - len(level_sets(top)[i]), den - den_top))
        for i, lv in enumerate(LEVELS)))
    for key in ("drop_automated_synthesis", "drop_compiled", "drop_primary_analysis",
                "retain_only_automated_synthesis", "retain_only_compiled",
                "retain_only_primary_analysis"):
        c_ = key[len("drop_"):] if key.startswith("drop_") else key[len("retain_only_"):]
        parts = [p for p in E if (bmap[p]["topology_assembly"] == c_)
                 == key.startswith("retain")]
        got[key] = srow(key, den, ladder_pcts(parts, den))
    bp = banded[0]
    scale = float(bmap[bp]["f_low"]) / float(bmap[bp]["f_mol"])
    low = OrderedDict()
    for i, lv in enumerate(LEVELS):
        n_ = lad[lv]["union"]
        if lv in ("C1", "C1d") and level_sets(bp)[i]:
            n_ = lad[lv]["sum"] - len(level_sets(bp)[i]) * (1 - scale)
        low[lv] = pct(n_, den)
    got["bird_band_low"] = srow("bird_band_low", den, low)
    got["union"] = srow("union", den, OrderedDict(
        [("C0", pct(lad["C0"]["sum"] + n_add, den)),
         ("C1", pct(lad["C1"]["sum"] + n_add, den)), ("C1d", None), ("C1b", None)]))
    for key, nm in (("width_rule", "Width"), ("strict", "Strict")):
        got[key] = srow(key, den, OrderedDict(
            [("C0", None), ("C1", None), ("C1d", None),
             ("C1b", strict_sets[nm]["C1b"]["pct_union"])]))
    supp_eq("supp_sensitivity", got, "alternative")

    got = OrderedDict()
    for r in sorted(euk, key=lambda r: (CATEGORY_ORDER.index(r["category"]),
                                        -as_int(r["described"]), r["partition"])):
        p, d_ = r["partition"], as_int(r["described"])
        n0 = len(c0.get(p, ())) if r["canonical_tree"] else 0
        cf = as_int(coldn[p]["col_filled"])
        nu = n0 + len(union_added.get(p, ()))
        got[p] = OrderedDict([
            ("partition", p), ("category", r["category"]), ("described", d_),
            ("described_source", r["described_source"]), ("col_count", cf),
            ("col_filled_from", coldn[p]["filled_classes"] or ""),
            ("named_pct", fmt_pct(pct(n0, d_)) if n0 else ""),
            ("named_col_pct", fmt_pct(pct(n0, cf)) if n0 and cf else ""),
            ("named_all_trees_pct", fmt_pct(pct(nu, d_)) if nu else "")])
    supp_eq("supp_denominators", got, "partition")
    sm("SuppColHolePct", "%.0f" % (100 * rule_from_code("col_denominators.py",
                                                        "HOLE_FRACTION")))

    # S8: Open Tree of Life
    R.section("Supplement S8: Open Tree of Life, per group")
    got = OrderedDict()
    for r in sorted(enum, key=lambda r: (CATEGORY_ORDER.index(bmap[r["partition"]]["category"]),
                                         r["partition"])):
        p = r["partition"]
        got[p] = OrderedDict([
            ("partition", p), ("category", bmap[p]["category"]),
            ("ot_phylo_species", r["otp"]), ("molecular_species", r["acc"]),
            ("both", r["accand"]), ("molecular_not_ot", r["accnot"]),
            ("both_pct", fix(pct(r["accand"], r["acc"]), 1) if r["acc"] else ""),
            ("molecular_keys", r["mkeys"]), ("keys_both", r["mand"]),
            ("keys_both_pct", fix(pct(r["mand"], r["mkeys"]), 1) if r["mkeys"] else "")])
    supp_eq("supp_opentree", got, "partition")

    # S9: every fit and the family sampling fractions
    R.section("Supplement S9: every fit of the worked example")

    def sig3(x):
        x = float(x)
        return "0" if x == 0 or abs(x) < 1e-9 else "%.3g" % x
    got = OrderedDict()
    for p in elig:
        for tree, conv in (F, M, MC, FM):
            f_ = fits[(p, tree, conv)]
            got["%s/%s/%s" % (p, tree, conv)] = OrderedDict([
                ("partition", p), ("tree", tree), ("rho_convention", conv),
                ("n_tips", as_int(f_["n_tips"])),
                ("rho", "" if conv == "family" else rho_fmt(float(f_["rho"]))),
                ("lambda", sig3(f_["lambda"])), ("mu", sig3(f_["mu"])),
                ("r", "%.4f" % float(f_["r"])), ("r_lo", "%.4f" % float(f_["r_lo"])),
                ("r_hi", "%.4f" % float(f_["r_hi"])),
                ("epsilon", "%.3f" % float(f_["epsilon"])),
                ("epsilon_lo", "%.3f" % float(f_["epsilon_lo"])),
                ("epsilon_hi", "%.3f" % float(f_["epsilon_hi"])),
                ("loglik", "%.2f" % float(f_["loglik"])),
                ("n_events", as_int(f_["n_events"])),
                ("n_polytomy_events", as_int(f_["n_polytomy_events"]))])
            R.eq("%s %s/%s: n_tips = the tree fitted" % (p, tree, conv),
                 as_int(f_["n_tips"]), len(tips[p]) if tree == "full" else len(mol_rows[p]))
            lam_, mu_ = float(f_["lambda"]), float(f_["mu"])
            R.eq("%s %s/%s: r = lambda - mu and mu = epsilon lambda, to the "
                 "precision of the file" % (p, tree, conv),
                 (abs(float(f_["r"]) - (lam_ - mu_)) <= 1e-5 * lam_,
                  abs(float(f_["epsilon"]) * lam_ - mu_) <= 1e-5 * lam_), (True, True))
    want = OrderedDict(("%s/%s/%s" % (w["partition"], w["tree"], w["rho_convention"]), w)
                       for w in supp_rows("supp_worked_fits"))
    R.eq("supp_worked_fits: rows, in order", list(got), list(want))
    for k, row in got.items():
        R.eq("supp_worked_fits: %s" % k,
             {c: (str(v), want[k].get(c)) for c, v in row.items()
              if str(v) != want[k].get(c)}, {})
    famf = read_csv(os.path.join(a.fits, "worked_example_bd_family_fractions.csv"))
    ford = read_csv(os.path.join(a.fits, "worked_example_bd_family_order.csv"))
    got = OrderedDict()
    for p in elig:
        fs = sorted(float(r["f"]) for r in famf if r["partition"] == p)
        n_ = len(fs)
        med = fs[n_ // 2] if n_ % 2 else 0.5 * (fs[n_ // 2 - 1] + fs[n_ // 2])
        fm, ft = fits[(p,) + FM], fits[(p,) + M]
        fr_ = [r for r in famf if r["partition"] == p and r["n_col_species"]]
        pooled = (sum(as_int(r["n_accepted_names"]) for r in fr_)
                  / float(sum(as_int(r["n_col_species"]) for r in fr_)))
        fo = {r["child"]: r for r in ford if r["partition"] == p}
        R.eq("%s: family fit with the value from the first child = the reported "
             "family fit" % p, (fo["first"]["r"], fo["first"]["loglik"]),
             (fm["r"], fm["loglik"]))
        got[p] = OrderedDict([
            ("partition", p), ("families", n_),
            ("tips_with_family", as_int(fm["n_tips_family_fraction"])),
            ("tips", as_int(ft["n_tips"])), ("f_min", "%.2f" % fs[0]),
            ("f_median", "%.2f" % med), ("f_max", "%.2f" % fs[-1]),
            ("f_pooled", rho_fmt(pooled)),
            ("rho_uniform", rho_fmt(float(ft["rho"]))),
            ("r_first", "%.4f" % float(fo["first"]["r"])),
            ("r_last", "%.4f" % float(fo["last"]["r"])),
            ("epsilon_first", "%.3f" % float(fo["first"]["epsilon"])),
            ("epsilon_last", "%.3f" % float(fo["last"]["epsilon"]))])
        R.eq("%s: family likelihood at equal fractions = uniform likelihood, "
             "within %g" % (p, UNIFORM_TOL),
             abs(float(fm["family_uniform_check"])) <= UNIFORM_TOL, True)
    supp_eq("supp_worked_families", got, "partition")
    sm("SuppFamilyMinTips", fmt_int(min(r["tips_with_family"] for r in got.values())))
    fo = defaultdict(dict)
    for r in ford:
        fo[r["partition"]][r["child"]] = r
    sm("SuppFamOrderMaxRDiffPct", "%.1f" % (100 * max(
        abs(float(v["last"]["r"]) - float(v["first"]["r"])) / float(v["first"]["r"])
        for v in fo.values())))
    sm("SuppFamOrderMaxEpsDiff", "%.3f" % max(
        abs(float(v["last"]["epsilon"]) - float(v["first"]["epsilon"]))
        for v in fo.values()))
    sm("SuppFamOrderMaxDLL", "%.1f" % max(
        abs(float(v["first"]["loglik_last_at_first_mle"]) - float(v["first"]["loglik"]))
        for v in fo.values()))
    supp_eq("supp_worked_polytomy", poly_rows, "partition")
    sm("SuppPolySeedsN", fmt_int(len({sd_ for p in poly for sd_ in poly[p]})))
    chk = read_csv(os.path.join(a.fits, "worked_example_bd_crosscheck.csv"))
    rd = 0.0
    for r in chk:
        lam = max(float(r["this_lambda"]), float(r["phytools_lambda"]))
        for x_, y_ in (("this_lambda", "phytools_lambda"), ("this_mu", "phytools_mu")):
            x, y = float(r[x_]), float(r[y_])
            if lam > 1e-6:
                rd = max(rd, abs(x - y) / lam)
    sm("SuppCrossCheckN", fmt_int(len(chk)))
    R.eq("rates agree with phytools within %g (relative)" % CROSSCHECK_RTOL,
         rd <= CROSSCHECK_RTOL, True)
    sm("SuppCrossCheckTolPct", "%g" % (100 * CROSSCHECK_RTOL))
    sm("SuppUniformCheckTol", r"$10^{%d}$" % round(math.log10(UNIFORM_TOL)))
    sm("SuppCrossCheckMaxTips", fmt_int(int(rule_from_code(
        "worked_example_bd.R", "CROSSCHECK_MAX_TIPS", r_assign=True))))
    sm("SuppProfileCutoff", "%.2f" % (NormalDist().inv_cdf(0.5 + INTERVAL_PCT / 200.0) ** 2 / 2))

    # S10: data and code
    R.section("Supplement S10: data and code")
    restored = [r for r in read_csv(path("registers", "restored_branch_lengths.csv"))
                if r["file"].startswith("Data/standardized/trees/") and r["changed"] == "yes"
                and os.path.basename(r["file"]) in files]
    sm("SuppRestoredTreesN", fmt_int(len(restored)))
    sm("SuppRestoredBranchLengthsN", fmt_int(sum(as_int(r["branch_lengths_restored"])
                                                 for r in restored)))
    sm("SuppShippedTreesN", fmt_int(len(files)))
    vr = Counter((r.get("rule") or "").strip()
                 for r in read_csv(path("standardized", "tree_references.csv")))
    sm("SuppRefBuiltN", fmt_int(vr["built_from_source"]))
    sm("SuppRefRelabelN", fmt_int(vr["tokens_equal_after_tip_relabel"]))
    sm("SuppRefEdgeN", fmt_int(vr["edge_length_multiset_equal"]))
    checks = rule_from_code("verify_atlas.py", "CHECKS", block=True)
    sm("SuppVerifyChecksN", fmt_int(len(checks)))
    R.eq("supp_checks: the named checks of code/rebuild/verify_atlas.py",
         [(r["check"], r["title"]) for r in supp_rows("supp_checks")],
         [tuple(c) for c in checks])
    bcols = set(basis[0])
    for r in supp_rows("supp_terms"):
        cols = [c.strip() for c in r["columns"].split(";") if c.strip()]
        R.eq("supp_terms: %s: columns in registers/coverage_basis.csv" % r["term"],
             [c for c in cols if c not in bcols], [])
        if r["file"]:
            R.eq("supp_terms: %s: %s is in this repository" % (r["term"], r["file"]),
                 os.path.exists(path(*r["file"].rstrip("/").split("/"))), True)

    # ------------------------------------------------------------- rules ---
    R.section("Rule parameters printed in the text")
    R.macro(macros, "RuleChronCV", "%g" % CHRONOGRAM_CV_MAX)
    R.macro(macros, "RuleDegenerateMyr", "%g" % DEGENERATE_WIDTH_MYR)
    R.macro(macros, "RuleBoundMyr", "%g" % BOUND_TOL_MYR)
    R.macro(macros, "RulePriorContainPct", "%g" % (100 * PRIOR_CONTAINMENT_MIN))
    R.macro(macros, "RulePriorWidthRatio", "%g" % PRIOR_WIDTH_RATIO_MIN)
    R.macro(macros, "RuleIntervalPct", "%g" % INTERVAL_PCT)
    R.macro(macros, "RuleLRCritical", "%.2f" % LR_CRITICAL)
    R.macro(macros, "RuleTestLevelPct", "%g" % (100 * TEST_LEVEL))

    # ----------------------------------------------------------- summary ---
    print("\n== Summary")
    print("  comparisons: %d, failures: %d" % (R.n, len(R.fail)))
    print("  printed macros recomputed: %d of %d"
          % (len(R.macros_checked & set(macros)), len(macros)))
    sm_all = set(R.supp_macros)
    print("  Supplement macros recomputed: %d of %d"
          % (len(R.macros_checked & sm_all), len(sm_all)))
    srest = sorted(sm_all - R.macros_checked)
    if srest:
        print("  Supplement macros not recomputed: %s" % ", ".join(srest))
        R.fail.append("Supplement macros not recomputed")
    rest = sorted(set(macros) - R.macros_checked)
    where = OrderedDict([
])
    if rest:
        R.fail.append("printed macros not recomputed")
        print("  not recomputed here (%d):" % len(rest))
        for pre, txt in where.items():
            k = [m for m in rest if m.startswith(pre)]
            if k:
                print("    %-4s %3d  %s" % (pre + "*", len(k), txt))
        other = [m for m in rest if not any(m.startswith(x) for x in where)]
        if other:
            print("    other: %s" % ", ".join(other))
    if not_rederived:
        print("  molecular tip lists taken as distributed (the source tree is not")
        print("  redistributed; registers/molecular_only_trees.csv cites it):")
        for p, src in not_rederived:
            print("    %-12s %s" % (p, src))
    if R.fail:
        print("\nFAILED: %d comparison(s) disagree." % len(R.fail))
        return 1
    print("\nPASSED: every recomputed number agrees with the paper.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
