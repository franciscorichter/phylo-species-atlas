#!/usr/bin/env python3
"""reproduce_paper.py -- recompute the paper's reported numbers from this repository.

Python 3, standard library only.  Run from anywhere:

    python3 scripts/reproduce_paper.py            # full run, about 15 seconds
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
    paper/worked_example_bd*.csv    the birth-death fits (Figure 5, Table 3);
                                    the fits themselves are rerun by
                                    code/rebuild/run_worked_example.py (R)
    registers/coverage_basis.csv    the per-group table
    registers/molecular_tips.csv    the molecular tips, re-derived from the
                                    distributed lists

It exits 0 when every comparison agrees and 1 otherwise.  The summary counts
the printed numbers recomputed and lists any that are not.

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
LOSS_REASONS = ("NotChron", "RelTime", "NoInterval", "Kind", "Calib")
UNREP_NAMED = 4                     # largest groups without a tree, named
LR_CRITICAL = 2.71                  # qchisq(0.90, 1), printed rule
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


def pct(a, b):
    return 100.0 * float(a) / float(b)


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
    for r in read_csv(path("registers", "pruned_tips.csv")):
        pruned[r["filename"]].add(r["standardized_id"])
    for r in read_csv(path("registers", "canonical_sources.csv")):
        f = r["filename"]
        src = path("standardized", "canonical_sources", os.path.basename(r["source"]))
        with open(src, encoding="utf-8", errors="replace") as fh:
            sleaves, _ = leaf_depths(fh.read())
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
    R.macro(macros, "UnrepLargestNames", ", ".join(big[:-1]) + " and " + big[-1])
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
            if g is None:
                R.eq("%s: GBIF row for unmatched query %s" % (p, x["query"]), None, "row")
                continue
            if (g["matchType"] in GBIF_MATCH and g["rank"] == "SPECIES"
                    and g["status"] in GBIF_STATUS and g["species"]
                    and g["species"] not in have):
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
    excl = Counter(why_of[p] for p in rows if p in intervals and why_of[p])
    for k, v in excl.items():
        R.macro(macros, "CrownExcl" + macro_name(k), fmt_int(v))
    R.macro(macros, "CrownNWithInterval", fmt_int(sum(p in intervals for p in rows)))
    R.macro(macros, "CrownNNoInterval", fmt_int(sum(p not in intervals for p in rows)))

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
    WE_LR_CRIT = 2.705543                      # qchisq(0.90, 1)
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
    n_ident = n_bound = n_rej = n_fg = n_fs = n_fo = n_famo = n_fb = n_fw = 0
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
        dfam = round(g(p, *FM, k="loglik") - ll_t, 4)
        n_fb += dfam > 0
        n_fw += dfam < 0
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
        R.macro(macros, "WEFamilyDLL" + n, "%.1f" % dfam)
    R.macro(macros, "WENEligible", fmt_int(len(elig)))
    R.macro(macros, "WENIdentity", fmt_int(n_ident))
    R.macro(macros, "WENBoundary", fmt_int(n_bound))
    R.macro(macros, "WENRejected", fmt_int(n_rej))
    R.macro(macros, "WENFullGreater", fmt_int(n_fg))
    R.macro(macros, "WENFullSmaller", fmt_int(n_fs))
    R.macro(macros, "WENFullOverlap", fmt_int(n_fo))
    R.macro(macros, "WENFamilyOverlap", fmt_int(n_famo))
    R.macro(macros, "WENFamilyBetter", fmt_int(n_fb))
    R.macro(macros, "WENFamilyWorse", fmt_int(n_fw))
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
        return "%.2f" % x if x >= 0.1 else "%.2g" % x

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
            ("dll_family", "%.1f" % round(float(fm["loglik"]) - float(ft["loglik"]), 4)),
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

    # ------------------------------------------------------------- rules ---
    R.section("Rule parameters printed in the text")
    R.macro(macros, "RuleChronCV", "%g" % CHRONOGRAM_CV_MAX)
    R.macro(macros, "RuleDegenerateMyr", "%g" % DEGENERATE_WIDTH_MYR)
    R.macro(macros, "RuleBoundMyr", "%g" % BOUND_TOL_MYR)
    R.macro(macros, "RulePriorContainPct", "%g" % (100 * PRIOR_CONTAINMENT_MIN))
    R.macro(macros, "RulePriorWidthRatio", "%g" % PRIOR_WIDTH_RATIO_MIN)
    R.macro(macros, "RuleIntervalPct", "%g" % INTERVAL_PCT)
    R.macro(macros, "RuleLRCritical", "%.2f" % LR_CRITICAL)

    # ----------------------------------------------------------- summary ---
    print("\n== Summary")
    print("  comparisons: %d, failures: %d" % (R.n, len(R.fail)))
    print("  printed macros recomputed: %d of %d"
          % (len(R.macros_checked & set(macros)), len(macros)))
    rest = sorted(set(macros) - R.macros_checked)
    where = OrderedDict([
])
    if rest:
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
