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

and checks the browser site: every number of site/data.json against the
macros, Table 1 (paper/table_which_tree.csv), paper/crown_ages.csv and
paper/fig1_coverage_bars.csv, and every tree file the site serves against its
file in standardized/trees/.
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
# D55: does the archive carry date uncertainty a user can propagate (level
# C1u), whatever its relation to the calibration prior.  A taxon-addition set
# varies only the placement of added tips (registers/taxon_addition_set_check.csv).
INTERVAL_KIND_AVAILABLE = {
    "node_hpd_annotation": True, "node_label_interval": True,
    "posterior_sample": True, "bootstrap_dating_sample": True,
    "calibration_prior_sample": True, "taxon_addition_sample": False,
    "date_assignment_sample": True}


def mixture_p(lr):
    """F3-15: P(T >= lr) for T the 50:50 mixture of a point mass at zero and
    chi2(1), i.e. P(Z >= sqrt(lr)) for a standard normal Z."""
    if lr <= 0:
        return 1.0
    return 0.5 * math.erfc(math.sqrt(lr) / math.sqrt(2.0))


def p_text(q):
    """A p-value as the paper prints it (math mode)."""
    if q >= 0.001:
        return "%.3f" % q
    e = int(math.floor(math.log10(q)))
    return "%.1f \\times 10^{%d}" % (q / 10 ** e, e)


def prior_quantile(spec, x):
    """F3-14: the quantile of x under the calibration density of
    registers/calibration_register.csv (column prior_density)."""
    form, _, rest = spec.partition("(")
    a = [float(v) for v in rest.strip().rstrip(")").split(",")]
    form = form.strip()
    if form == "uniform":
        q = min(1.0, max(0.0, (x - a[0]) / (a[1] - a[0])))
    elif form == "soft_uniform":         # 95% inside soft bounds
        if x < a[0]:
            return "<0.025"
        if x > a[1]:
            return ">0.975"
        q = 0.025 + 0.95 * (x - a[0]) / (a[1] - a[0])
    elif form == "offset_exp_q95":       # 95% quantile at the maximum
        q = 0.0 if x <= a[0] else 1.0 - math.exp(
            math.log(0.05) * (x - a[0]) / (a[1] - a[0]))
    elif form == "offset_gamma":         # integer shape: Erlang
        z = max(0.0, (x - a[0]) / a[2])
        q = 1.0 - sum(math.exp(-z) * z ** k / math.factorial(k)
                      for k in range(int(a[1])))
    else:
        raise SystemExit("unknown prior_density %r" % spec)
    return "%.3f" % q
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
LEVELS = ("C0", "C1", "C1d", "C1u")
LEVEL_MACRO = ("CZero", "COne", "COneD", "COneU")
# Table 1 marks a group whose tips exceed its named species by this fraction
TIPS_OVER_NAMED = 0.10
# interval kinds that come from a set of trees
SET_KINDS = ("posterior_sample", "bootstrap_dating_sample", "calibration_prior_sample",
             "taxon_addition_sample", "date_assignment_sample")
LOSS_REASONS = ("NotChron", "RelTime", "NoInterval", "Placement")
LOSS_STEP = {"TaxPlaced": "C0>C1", "NotChron": "C1>C1d", "RelTime": "C1>C1d",
             "NoInterval": "C1d>C1u", "Placement": "C1d>C1u"}
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
    if not names:
        return ""
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


# ---------------------------------------------------------------------------
# The browser site (site/): its numbers are the paper's, its trees the release's
# ---------------------------------------------------------------------------

SITE_RETIRED_FIELDS = ("coverage_pct", "crown_ma", "tip_taxonomy", "redistributed",
                       "unique_species_in_trees", "total_described_species")


def _nwk_parse(text):
    """Newick -> [label, length, children]; lengths kept as strings."""
    s = text.strip().rstrip(";")
    pos = [0]

    def node():
        kids = None
        if s[pos[0]] == "(":
            pos[0] += 1
            kids = []
            while True:
                kids.append(node())
                c = s[pos[0]]
                pos[0] += 1
                if c == ")":
                    break
        j = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ",():":
            pos[0] += 1
        label, length = s[j:pos[0]], None
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            j = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ",()":
                pos[0] += 1
            length = s[j:pos[0]]
        return [label, length, kids]
    return node()


def _nwk_tips(n):
    return [n[0]] if n[2] is None else [t for k in n[2] for t in _nwk_tips(k)]


def _nwk_write(n):
    out = "(" + ",".join(_nwk_write(k) for k in n[2]) + ")" if n[2] is not None else ""
    return out + n[0] + (":" + n[1] if n[1] is not None else "")


def one_tip_per_species(text, names):
    """The site's derived tree: per Genus_species the tip with the smallest
    id; other tips removed; a node left with one child joined to it, the
    branch lengths added as decimals."""
    from decimal import Decimal
    tree = _nwk_parse(text)
    best = {}
    for t in _nwk_tips(tree):
        nm = names.get(t, t)
        p = nm.split("_")
        key = nm
        if len(p) >= 2 and len(p[0]) >= 2 and p[0][:1].isupper() and p[0][1:].islower() \
                and p[1][:1].islower() and p[1].replace(".", "").replace("-", "").isalpha():
            key = p[0] + "_" + p[1]
        if key not in best or int(t) < int(best[key]):
            best[key] = t
    keep = set(best.values())

    def prune(n):
        if n[2] is None:
            return n if n[0] in keep else None
        kids = [c for c in (prune(k) for k in n[2]) if c is not None]
        if not kids:
            return None
        if len(kids) == 1:
            c = kids[0]
            ln = c[1] if n[1] is None else (n[1] if c[1] is None
                                             else str(Decimal(n[1]) + Decimal(c[1])))
            return [c[0], ln, c[2]]
        return [n[0], n[1], kids]
    return _nwk_write(prune(tree)) + ";\n"


def check_site(R, macros, basis, head):
    """site/data.json against the paper; the site's tree files against the
    release's files."""
    R.section("The browser site: site/data.json and site/data/partitions")
    fp = path("site", "data.json")
    R.eq("site/data.json exists", os.path.exists(fp), True)
    if not os.path.exists(fp):
        return
    raw = open(fp, encoding="utf-8").read()
    site = json.loads(raw)
    unmacro = lambda v: None if v is None else v.replace("{,}", ",")
    mac = lambda k: unmacro(macros.get(k))
    t1 = {r["partition"]: r for r in read_csv(path("paper", "table_which_tree.csv"))}
    crowns = {r["partition"]: r for r in read_csv(path("paper", "crown_ages.csv"))}
    bars = {r["group"]: r for r in read_csv(path("paper", "fig1_coverage_bars.csv"))}
    rows = {r["partition_group"]: r for r in site.get("coverage", [])}
    R.eq("site groups = registers/coverage_basis.csv groups",
         sorted(rows), sorted(r["partition"] for r in basis))

    for b in basis:
        p = b["partition"]
        r = rows.get(p)
        if r is None:
            continue
        mn = macro_name(p)
        d = int(b["described"])
        R.eq("site %s: described species = coverage_basis" % p, r["described_species"], d)
        if "Described" + mn in macros:
            R.eq("site %s: described = \\Described%s" % (p, mn), r["described_text"],
                 mac("Described" + mn))
        R.eq("site %s: reference tree = coverage_basis" % p, r["filename"],
             b["canonical_tree"] or None)
        R.eq("site %s: estimated total = paper/fig1_coverage_bars.csv" % p,
             (r.get("estimated_total"), r.get("estimated_low"), r.get("estimated_high")),
             tuple(as_int(bars.get(p, {}).get(k), None) for k in
                   ("estimated_total", "estimated_low", "estimated_high")))
        if not b["canonical_tree"]:
            R.eq("site %s: no levels without a reference tree" % p, r.get("levels"), None)
            continue
        L = r["levels"]
        c0, c1 = as_int(b["c0_units"]), as_int(b["c1_units"])
        dated = b["is_chronogram"] == "yes"
        unc = dated and b["uncertainty_available"] == "yes"
        R.eq("site %s: level counts = coverage_basis" % p,
             tuple(L[k]["n"] for k in ("named", "molecular", "dated", "uncertainty")),
             (c0, c1, c1 if dated else 0, c1 if unc else 0))
        R.eq("site %s: level percentages = counts over described" % p,
             tuple(L[k]["pct_text"] for k in ("named", "molecular", "dated", "uncertainty")),
             tuple(fmt_pct(pct(n, d)) for n in (c0, c1, c1 if dated else 0, c1 if unc else 0)))
        if b["is_prokaryotic"] == "yes":
            pl = head["prokaryotic_line"]["per_partition"][p]
            R.eq("site %s: named = prokaryotic line" % p, L["named"]["pct_text"],
                 fmt_pct(pl["C0_pct"]))
            continue
        for lev, key in (("named", "CZero"), ("molecular", "COne")):
            R.eq("site %s: %s %% = \\%sPct%s" % (p, lev, key, mn), L[lev]["pct_text"],
                 mac(key + "Pct" + mn))
            R.eq("site %s: %s = \\%sN%s" % (p, lev, key, mn), L[lev]["n_text"],
                 mac(key + "N" + mn))
        t = t1.get(p)
        R.eq("site %s: in Table 1" % p, t is not None, True)
        if t:
            R.eq("site %s: Table 1 tips, named %%, molecular %%" % p,
                 (r["tips_text"], L["named"]["pct_text"], L["molecular"]["pct_text"]),
                 ("{:,}".format(int(t["tips"])), t["named_pct"], t["molecular_pct"]))
            R.eq("site %s: Table 1 dates (A) and date uncertainty (s, i)" % p,
                 (r["dated"], r["uncertainty_available"]),
                 (t["time_scale"] == "absolute", t["uncertainty_status"] in ("set", "summary")))
        c = crowns.get(p)
        R.eq("site %s: crown age only from the crown register" % p,
             r.get("crown") is not None, c is not None)
        if c:
            R.eq("site %s: crown age = \\CrownAge%s" % (p, mn), r["crown"]["age_ma_text"],
                 mac("CrownAge" + mn))
            if unc:
                R.eq("site %s: crown interval = \\CrownLo%s, \\CrownHi%s" % (p, mn, mn),
                     (r["crown"]["interval_lo_ma_text"], r["crown"]["interval_hi_ma_text"]),
                     (mac("CrownLo" + mn), mac("CrownHi" + mn)))

    s = site["summary"]
    lv = {x["key"]: x for x in site["levels"]}
    for key, m in (("named", "CZero"), ("molecular", "COne"), ("dated", "COneD"),
                   ("uncertainty", "COneU")):
        R.eq("site total %s = \\%sN, \\%sPct, \\%sNPart" % (key, m, m, m),
             (lv[key]["n_text"], lv[key]["pct_text"], str(lv[key]["n_groups"])),
             (mac(m + "N"), mac(m + "Pct"), mac(m + "NPart")))
    R.eq("site totals: described, groups, reference trees, trees",
         (s["described_eukaryotic_text"], s["described_prokaryotic_text"],
          str(s["n_partitions"]), str(s["n_partitions_eukaryotic"]),
          str(s["n_partitions_with_tree"]), str(s["n_partitions_with_tree_eukaryotic"]),
          str(s["n_trees"])),
         (mac("DescribedEuk"), mac("DescribedProk"), mac("NumPartitionsAll"),
          mac("NumPartitionsEuk"), mac("NumCanonical"), mac("NumCanonicalEuk"),
          mac("NumTreesShipped")))
    al = head["denominators"]["all_life"]
    R.eq("site estimated totals = headline all-life estimates",
         (s["estimated_total_low_text"], s["estimated_total_high_text"]),
         ("{:,}".format(al["estimated_low"]), "{:,}".format(al["estimated_high"])))

    # Trees: every redistributed tree is served from standardized/trees/
    meta = [r["filename"] for r in read_csv(path("standardized", "metadata.csv"))]
    withdrawn = {r["filename"] for r in read_csv(path("standardized", "withdrawn.csv"))}
    trees = site["trees"]
    pkg = [t for t in trees if t["in_package"]]
    R.eq("site trees in the package = standardized/metadata.csv",
         sorted(t["filename"] for t in pkg), sorted(meta))
    R.eq("site lists no withdrawn tree (timetree.nwk)",
         sorted(w for w in withdrawn if w in raw or any(
             t["filename"] == w or t["group"] == w[:-4] for t in trees)), [])
    bad = [t["filename"] for t in pkg
           if t["tree_url"] != "standardized/trees/" + t["filename"]
           or not os.path.exists(path("standardized", "trees", t["filename"]))]
    R.eq("site serves every redistributed tree from standardized/trees/", bad, [])
    canon = {b["canonical_tree"]: b["partition"] for b in basis if b["canonical_tree"]}
    R.eq("site reference trees = coverage_basis reference trees, served from "
         "standardized/trees/",
         sorted((t["filename"], t["partition_group"]) for t in pkg if t["is_partition_canonical"]),
         sorted(canon.items()))
    tipre = re.compile(r"[(,]\s*([^(),:;\s]+)")
    bad = [t["filename"] for t in pkg if t["ntips"] != len(tipre.findall(
        open(path("standardized", "trees", t["filename"]), encoding="utf-8").read()))]
    R.eq("site tip counts = tips of the files", bad, [])

    # site/data/partitions: no copy of a reference tree; every other tree file
    # is identical to its release file, or derived from it by the stated rule,
    # or a tree the release does not hold (labels are names, never ids)
    names = None
    files, problems = [], []
    for d, _, fs in os.walk(path("site")):
        for f in fs:
            if f.endswith(".nwk"):
                files.append(os.path.relpath(os.path.join(d, f), ROOT).replace(os.sep, "/"))
    served = {t.get("tree_url"): t for t in trees}
    for rel in sorted(files):
        parts = rel.split("/")
        if len(parts) == 5 and parts[4] == "tree.nwk" or "uncertainty" in parts[:5] \
                and parts.index("uncertainty") == 4:
            problems.append("%s: a copy of a reference tree" % rel)
            continue
        if "sub-phylos" not in parts:
            problems.append("%s: not a sub-clade tree" % rel)
            continue
        slug = parts[parts.index("sub-phylos") + 1]
        std = path("standardized", "trees", slug + ".nwk")
        text = open(path(*parts), encoding="utf-8").read()
        if os.path.exists(std):
            stext = open(std, encoding="utf-8").read()
            if text == stext:
                continue
            t = served.get(rel)
            if t is None or t.get("site_tree") != "derived":
                problems.append("%s differs from standardized/trees/%s.nwk" % (rel, slug))
                continue
            if names is None:
                names = {}
                with open(path("standardized", "dictionary.csv"), newline="",
                          encoding="utf-8") as fh:
                    for r in csv.DictReader(fh):
                        names[r["id"]] = r["standardized_name"]
            if one_tip_per_species(stext, names) != text:
                problems.append("%s is not one tip per species of standardized/trees/%s.nwk"
                                % (rel, slug))
        else:
            labels = tipre.findall(text)
            if any(x.lstrip("-").isdigit() for x in labels):
                problems.append("%s: numeric tip ids in a tree the release does not hold" % rel)
            if served.get(rel, {}).get("site_tree") != "site_only":
                problems.append("%s: not listed as a site-only tree" % rel)
    R.eq("site tree files: no reference-tree copy; each other file = its release "
         "file, its stated derivation, or a site-only tree", problems, [])
    js = open(path("site", "app.js"), encoding="utf-8").read()
    R.eq("site/app.js reads no retired field",
         [f for f in SITE_RETIRED_FIELDS if re.search(r"\.%s\b" % f, js)], [])


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
    src_counts = {}     # T route: source labels, merged, unmatched (F4-24)
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
            if src.startswith("Data/raw/"):
                tk = {key_of(nm) for nm in names_in_tree}
                nl = len(set(labels))
                src_counts[p] = (nl, nl - len(keys), len(keys - tk))
                R.eq("%s: source labels - merged - unmatched = molecular tips"
                     % p, nl - src_counts[p][1] - src_counts[p][2], len(listed))
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

    # ------------------------------------------------- interval (C1u) ---
    R.section("Date uncertainty available (C1u) and the comparison of each crown "
              "interval with its calibrations (a description, D58)")
    counts, why_of, avail_of = {}, {}, {}
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
        # rule (p) compares the interval with the crown's own two-sided
        # calibration, or, for a crown without one, with the band that the
        # calibrations on the nodes around it allow (D52, prior_range band)
        if iv.get("prior_range") in ("calibration", "band"):
            plo, phi = float(iv["prior_lo_ma"]), float(iv["prior_hi_ma"])
            cont = max(0.0, min(hi, phi) - max(lo, plo)) / w if w > 0 else 1.0
            ratio = w / (phi - plo)
            repro = cont >= PRIOR_CONTAINMENT_MIN and ratio >= PRIOR_WIDTH_RATIO_MIN
            R.eq("%s: prior containment and width ratio" % p,
                 ("%.6f" % cont, "%.6f" % ratio),
                 (iv["prior_containment"], iv["prior_width_ratio"]))
        why = ("prior_sampled" if iv["kind"] in PRIOR_SAMPLE_KINDS
               else "kind" if not INTERVAL_KIND_COUNTS[iv["kind"]]
               else "degenerate" if degenerate else "at_bound" if at_bound
               else "prior_reproduced" if repro else "")
        counts[p] = not why
        why_of[p] = why
        avail_of[p] = INTERVAL_KIND_AVAILABLE[iv["kind"]]
        R.eq("%s: interval classification" % p,
             (why, "yes" if not why else "no",
              "yes" if avail_of[p] else "no"),
             (iv["rule_exclusion"], iv["not_excluded_by_rule"],
              iv["uncertainty_available"]))
        R.eq("%s: a crown interval no rule excludes has its date "
             "uncertainty available (D55)" % p,
             bool(avail_of[p] or not counts[p]), True)
    calreg = {r["partition"]: r for r in read_csv(
        path("registers", "calibration_register.csv"))}
    for p, iv in sorted(intervals.items()):
        if INTERVAL_KIND_COUNTS[iv["kind"]]:
            cr = calreg.get(p, {})
            # the range enters the test for a two-sided calibration on the
            # crown, otherwise for the band the register records (D52)
            two = cr.get("calibration") == "two_sided"
            band = (not two and bool(cr.get("band_lo_ma"))
                    and bool(cr.get("band_hi_ma")))
            src = "calibration" if two else "band" if band else ""
            f6 = lambda x: "%.6f" % float(x) if x else ""
            want = ((cr.get("prior_lo_ma"), cr.get("prior_hi_ma")) if two
                    else (cr.get("band_lo_ma"), cr.get("band_hi_ma")) if band
                    else ("", ""))
            got = ((iv["prior_lo_ma"], iv["prior_hi_ma"]) if src else ("", ""))
            R.eq("%s: calibration, range and its source as in the calibration register" % p,
                 (iv["calibration"], iv.get("prior_range", ""), f6(got[0]), f6(got[1])),
                 (cr.get("calibration"), src, f6(want[0]), f6(want[1])))
    available = {p: chron.get(p, False) and avail_of.get(p, False) for p in c0}
    # D55: a taxon-addition set is not date uncertainty: pruned to the
    # sequenced tips, every tree of it has the dates of the source's dated tree
    tac = {r["partition"]: r for r in read_csv(
        path("registers", "taxon_addition_set_check.csv"))}
    for p, iv in sorted(intervals.items()):
        if iv["kind"] == "taxon_addition_sample":
            t_ = tac.get(p, {})
            R.eq("%s: taxon-addition set leaves the dates of the sequenced "
                 "species unchanged (largest difference below %g Myr)"
                 % (p, DEGENERATE_WIDTH_MYR),
                 (t_.get("set"), float(t_.get("max_abs_branching_time_diff_myr")
                                       or "inf") < DEGENERATE_WIDTH_MYR),
                 (iv["source"], True))
    # Per-group table columns
    R.section("Per-group table (registers/coverage_basis.csv)")
    for r in canon:
        p = r["partition"]
        R.eq("%s: named, molecular, basis, unresolved, tips, molecular fraction" % p,
             (len(c0[p]), len(c1[p]), c1_basis[p], unresolved[p], len(tips[p]),
              "%.6f" % fmol[p], "yes" if available[p] else "no"),
             (as_int(r["c0_units"]), as_int(r["c1_units"]), r["c1_basis"],
              as_int(r["c0_unresolved_labels"]), as_int(r["labels_unique"]),
              r["f_mol"], r["uncertainty_available"]))

    # ------------------------------------------------------- the levels ---
    R.section("The four levels over the %d eukaryotic groups" % len(euk))
    E = [r["partition"] for r in canon_euk]

    def level_sets(p):
        return (c0[p], c1[p], c1[p] if chron[p] else set(),
                c1[p] if available[p] else set())

    def ladder(parts, denominator):
        out = OrderedDict()
        for i, lv in enumerate(LEVELS):
            sets = [level_sets(p)[i] for p in parts]
            s = sum(len(x) for x in sets)
            u = len(set().union(*sets)) if sets else 0
            out[lv] = dict(sum=s, union=u, pct_sum=pct(s, denominator),
                           pct_union=pct(u, denominator),
                           n_part=sum(1 for x in sets if x))
        return out

    lad = ladder(E, den)
    R.eq("eukaryotic denominator", den, head["denominators"]["eukaryotic"]["described"])
    for lv, nm in zip(LEVELS, LEVEL_MACRO):
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
    for (a_, an), (b_, bn) in zip(list(zip(LEVELS, LEVEL_MACRO))[1:],
                                  list(zip(LEVELS, LEVEL_MACRO))[:-1]):
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
    for i, (lv, nm) in enumerate(zip(LEVELS, LEVEL_MACRO)):
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
        elif available[p]:
            reason[p] = ""
        elif p not in intervals:
            reason[p] = "NoInterval"
        else:
            reason[p] = "Placement"
    best = {}
    for p in sorted(E):
        for sp in c0[p]:
            lvl = (3 if sp in c1[p] and available[p]
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
        elif lvl < len(LEVELS) - 1:
            R.eq("%s: loss reason %s is a loss from %s" % (p, reason[p], LEVELS[lvl]),
                 LOSS_STEP[reason[p]].split(">")[0], LEVELS[lvl])
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
    for k in ("NoInterval", "Placement"):
        R.macro(macros, "Lost%sOfCOneDPct" % k, fmt_pct(pct(lost[k], lad["C1d"]["union"])))
    R.eq("species named in two or more groups",
         sum(1 for v in Counter(sp for p in E for sp in c0[p]).values() if v > 1),
         cu["species_in_several_groups"])

    # ------------------------------------------------------ sensitivities ---
    R.section("Sensitivity analyses")
    s = head["sensitivities"]
    # Denominator restricted to groups with a reference tree
    rep = ladder(E, den_rep)
    for lv, nm in zip(LEVELS, LEVEL_MACRO):
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
    for lv, nm in zip(LEVELS, LEVEL_MACRO):
        R.macro(macros, "SensCol%sPct" % nm, fmt_pct(col[lv]["pct_union"]))
    # Molecular coverage as named coverage times the molecular fraction
    prod = {lv: 0.0 for lv in ("C1", "C1d")}
    for p in E:
        x = len(c0[p]) * float("%.6f" % fmol[p])
        prod["C1"] += x
        prod["C1d"] += x if chron[p] else 0.0
    R.macro(macros, "SensProductCOnePct", fmt_pct(pct(round(prod["C1"], 4), den)))
    R.macro(macros, "SensProductCOneDPct", fmt_pct(pct(round(prod["C1d"], 4), den)))
    # D58: the crowns whose interval no rule excludes, by the calibration on
    # the crown (a description in the Supplement; no level rests on it)
    bounded = {p: chron.get(p, False) and counts.get(p, False) for p in c0}
    prun = sorted(r["partition"] for r in read_csv(
        path("registers", "calibration_register.csv"))
        if r.get("prior_basis") == "effective_prior_run")
    R.macro(macros, "CrownNPriorRun", fmt_int(len(prun)))
    R.macro(macros, "CrownNamesPriorRun", name_list(prun))
    R.eq("no interval that no rule excludes sits on an unchecked crown calibration",
         sorted(p for p in c0 if bounded[p]
                and intervals[p]["calibration"] == "not_checked"), [])
    for k_, nm_ in (("none", "NoCal"), ("one_sided", "OneSided"),
                    ("two_sided", "TwoSided")):
        R.macro(macros, "CrownCounted%sNames" % nm_, name_list(sorted(
            p for p in c0 if bounded[p] and intervals[p]["calibration"] == k_)))
    byu = sorted(((len(c1[p]), p) for p in E if available[p]), reverse=True)
    R.macro(macros, "COneULargestName", byu[0][1])
    R.macro(macros, "COneULargestN", fmt_int(byu[0][0]))
    R.macro(macros, "COneULargestSharePct",
            "%.0f" % (100.0 * byu[0][0] / sum(x for x, _ in byu)))
    # D61: date uncertainty only where the archive holds dated trees of a
    # posterior or bootstrap dating of sequence data with a free crown
    seqd = sorted(p for p in E if available[p]
                  and intervals[p]["kind"] in ("posterior_sample",
                                               "bootstrap_dating_sample")
                  and intervals[p]["degenerate"] == "no"
                  and intervals[p]["at_bound"] == "no")
    seq_lad = ladder(seqd, den)
    R.macro(macros, "SensSeqDatingCOneUN", fmt_int(seq_lad["C1u"]["union"]))
    R.macro(macros, "SensSeqDatingCOneUPct", fmt_pct(seq_lad["C1u"]["pct_union"]))
    R.macro(macros, "SensSeqDatingOfCOneD", fmt_pct(pct(
        seq_lad["C1u"]["union"], lad["C1d"]["union"])))
    R.macro(macros, "SensSeqDatingNPart", fmt_int(len(seqd)))
    R.macro(macros, "SensSeqDatingNames", name_list(seqd))
    R.macro(macros, "DescribedUntracedNames", name_list(sorted(
        r["partition"] for r in euk
        if "source not traced" in (r.get("described_source") or ""))))
    R.macro(macros, "RuleTipsOverNamedPct", "%.0f" % (100 * TIPS_OVER_NAMED))
    R.macro(macros, "RuleElasticityPct", "%.0f" % (100 * 0.2))
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
    # F3-16: what the branch lengths of each tree of sequenced species are, as
    # the register records them (the trees are not redistributed)
    modat = {r["partition"]: r["branch_lengths"] for r in read_csv(
        path("registers", "molecular_only_trees.csv"))}
    R.eq("branch-length classes of the trees of sequenced species are known",
         sorted(set(modat.values()) - {"time", "substitutions", "none"}), [])
    R.macro(macros, "NMolOnlyDated", fmt_int(sum(v == "time" for v in modat.values())))
    R.macro(macros, "MolOnlyDatedNames", name_list(sorted(
        p for p, v in modat.items() if v == "time")))
    R.macro(macros, "MolOnlyPhylogramNames", name_list(sorted(
        p for p, v in modat.items() if v == "substitutions")))
    R.macro(macros, "MolOnlyNoLengthNames", name_list(sorted(
        p for p, v in modat.items() if v == "none")))
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
    width_pct = {}
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
            width_pct[p] = 100.0 * (hi - lo) / age
    # F3-14, F3-11: where each interval lies in its stated prior density, and
    # the containment and width ratio of every two-sided prior
    for p in rows:
        iv = intervals.get(p)
        if not iv:
            continue
        n = macro_name(p)
        spec = calreg.get(p, {}).get("prior_density") or ""
        if spec:
            cx = (float(iv["interval_center_ma"]) if iv["interval_center_ma"]
                  else float(iv["source_point_height_ma"]))
            for nm, x in (("Lo", float(iv["interval_lo_ma"])), ("Center", cx),
                          ("Hi", float(iv["interval_hi_ma"]))):
                q = prior_quantile(spec, x)
                R.macro(macros, "CrownPriorQ%s%s" % (nm, n),
                        "$%s$%s" % (q[0], q[1:]) if q[:1] in "<>" else "%.2f" % float(q))
        if iv.get("prior_containment"):
            R.macro(macros, "CrownPriorIn" + n,
                    "%.0f" % (100 * float(iv["prior_containment"])))
            R.macro(macros, "CrownPriorRatio" + n,
                    "%.2f" % float(iv["prior_width_ratio"]))
            R.macro(macros, "CrownPriorLo" + n, "%.1f" % float(iv["prior_lo_ma"]))
            R.macro(macros, "CrownPriorHi" + n, "%.1f" % float(iv["prior_hi_ma"]))
    wc_ = [width_pct[p] for p in rows if bounded[p]]
    R.macro(macros, "CrownWidthPctMinCounted", "%.1f" % min(wc_))
    R.macro(macros, "CrownWidthPctMaxCounted", "%.1f" % max(wc_))
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
    R.macro(macros, "CrownNamesNoInterval", name_list(sorted(p for p in rows
                                                             if p not in intervals)))
    R.macro(macros, "CrownNamesKind", name_list(sorted(
        p for p in rows if p in intervals and why_of[p] == "kind")))
    R.macro(macros, "CrownNamesCalib", name_list(sorted(
        p for p in rows if p in intervals and why_of[p] and why_of[p] != "kind")))
    # D55: date uncertainty available (C1u)
    av_ = sorted(p for p in rows if p in intervals and avail_of[p])
    R.macro(macros, "CrownNAvailable", fmt_int(len(av_)))
    R.macro(macros, "CrownNamesAvailable", name_list(av_))
    avs_ = [p for p in av_ if intervals[p]["kind"] in SET_KINDS]
    avm_ = [p for p in av_ if intervals[p]["kind"] not in SET_KINDS]
    R.macro(macros, "CrownNAvailableSet", fmt_int(len(avs_)))
    R.macro(macros, "CrownNamesAvailableSet", name_list(avs_))
    R.macro(macros, "CrownNAvailableSummary", fmt_int(len(avm_)))
    R.macro(macros, "CrownNamesAvailableSummary", name_list(avm_))
    R.macro(macros, "CrownNamesAvailableNotBeyond",
            name_list([p for p in av_ if not counts[p]]))
    R.macro(macros, "CrownNamesPlacement", name_list(sorted(
        p for p in rows if p in intervals and not avail_of[p])))
    R.macro(macros, "CrownNamesDateDraw", name_list(
        [p for p in av_ if why_of[p] == "kind"]))
    # how a source built its tree set, as the source states it, against the
    # measured set (registers/tree_set_designs.csv)
    for d_ in read_csv(path("registers", "tree_set_designs.csv")):
        p = d_["partition"]
        nd, dr, nt = (as_int(d_[k]) for k in ("n_dated_trees",
                                              "draws_per_dated_tree", "n_trees"))
        ns, na = as_int(d_["n_sequenced"]), as_int(d_["n_added"])
        R.eq("%s: set design, dated trees x draws = trees in the measured set" % p,
             (nd * dr, d_["set_file"]),
             (as_int(intervals[p]["n_trees"]), intervals[p]["source"]))
        R.eq("%s: set design, sequenced + added = tips; sequenced = molecular tips" % p,
             (ns + na, ns), (len(tips[p]), len(mol_rows[p])))
        n = macro_name(p)
        for k_, v_ in (("NDated", nd), ("NDraws", dr), ("NTrees", nt),
                       ("NSeq", ns), ("NAdded", na)):
            R.macro(macros, "Design%s%s" % (n, k_), fmt_int(v_))
    # F4-23: the shark set comparison (registers/shark_set_comparison.csv)
    for r_ in read_csv(path("registers", "shark_set_comparison.csv")):
        R.eq("shark set comparison: internal nodes = sequenced tips - 1",
             int(r_["n_internal_nodes"]), len(mol_rows[r_["partition"]]) - 1)
        for k_, c_ in (("NCompared", "n_compared"), ("NNodes", "n_internal_nodes"),
                       ("TopoEqual", "pruned_topology_equal"),
                       ("AbsentMedian", "absent_clades_median"),
                       ("AbsentMin", "absent_clades_min"),
                       ("AbsentMax", "absent_clades_max"),
                       ("AbsentGenusMedian", "absent_within_genus_median"),
                       ("EqualAgeMedian", "equal_age_median"),
                       ("EqualAgeMin", "equal_age_min"),
                       ("EqualAgeMax", "equal_age_max"),
                       ("YoungMyr", "young_myr"),
                       ("YoungComplete", "young_complete_median"),
                       ("YoungSeq", "young_sequenced_median")):
            R.macro(macros, "SharkCmp" + k_, fmt_int(r_[c_]))
        for k_, c_ in (("GammaComplete", "gamma_complete_median"),
                       ("GammaSeq", "gamma_sequenced_median"),
                       ("RootDiff", "root_diff_median_myr")):
            R.macro(macros, "SharkCmp" + k_, r_[c_].replace("-", "$-$"))
        R.macro(macros, "SharkCmpEqualAgePct", "%g" % (100 * float(r_["equal_age_tol"])))
    # tree sets excluded at the crown by a calibration rule (the classification
    # is of the crown only)
    R.macro(macros, "CrownNamesSetCalib", name_list(sorted(
        p for p in oset if why_of[p] and why_of[p] != "kind")))
    R.macro(macros, "CrownNamesPostCalib", name_list(sorted(
        p for p in oset if why_of[p] and why_of[p] != "kind"
        and intervals[p]["kind"] == "posterior_sample")))

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
            row.update(mkeys=len(mkeys), mand=len(mkeys & o), okeys=len(o),
                       acc=len(macc),
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
    R.macro(macros, "OTKeysMolecularAndOTPct", fmt_pct(pct(
        sum(r["mand"] for r in enum), sum(r["mkeys"] for r in enum))))
    R.macro(macros, "OTKeysOTPhyloCompared", fmt_int(sum(r["okeys"] for r in enum)))
    R.macro(macros, "OTKeysMolecularAndOTOfOTPct", fmt_pct(pct(
        sum(r["mand"] for r in enum), sum(r["okeys"] for r in enum))))
    notc = [r for r in ot_rows if not r["enum"]]
    R.macro(macros, "OTPhyloCompared", fmt_int(sum(r["otp"] for r in enum)))
    R.macro(macros, "OTPhyloNotComparedN", fmt_int(sum(r["otp"] for r in notc)))
    R.macro(macros, "OTPhyloNotComparedNames",
            name_list(sorted(r["partition"] for r in notc)))
    R.macro(macros, "OTPhyloOutsideGroups", fmt_int(euk_ot - tot_ot))
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
    # A test only on a binary summary tree (D59): the likelihood reads a node
    # with k children as k - 1 simultaneous speciation events, of probability
    # zero under the model.  The family is the eligible groups whose pruned
    # tree has no polytomy (counted on the tree below), each test at the
    # Bonferroni level over that family.
    # A test is reported only where its model assumptions hold for the tree
    # (D60): the summary tree is also dated by a dating analysis, not by
    # published ages assigned to some nodes and interpolated at the others
    # (registers/summary_tree_node_ages.csv, measured by
    # code/rebuild/summary_tree_node_ages.py from the source's files).
    nage = {r["partition"]: r for r in read_csv(
        path("registers", "summary_tree_node_ages.csv"))}
    R.eq("summary_tree_node_ages.csv: one row per eligible group",
         sorted(nage), sorted(elig))
    poly_ = [p for p in elig
             if as_int(fits[(p, "molecular", "tips")]["n_polytomy_events"]) != 0]
    tested = [p for p in elig if p not in poly_
              and nage[p]["node_ages"] == "dating_analysis"]
    R.macro(macros, "WENamesNotTestedPoly", name_list(poly_) or "none")
    R.macro(macros, "WENamesNotTestedAssigned", name_list(
        [p for p in elig if p not in tested and p not in poly_]) or "none")
    for p in sorted(nage):
        r_ = nage[p]
        if r_["node_ages"] != "date_assignment":
            continue
        v_ = {k: as_int(r_[k]) for k in (
            "n_tips", "n_internal", "n_internal_published", "n_tips_pruned",
            "n_internal_pruned", "n_internal_pruned_published", "n_tips_placed",
            "n_internal_placed", "n_polytomy_events_placed")}
        R.eq("%s: node-age register pruned and placed tips = the fit's "
             "sequenced tips" % p, (v_["n_tips_pruned"], v_["n_tips_placed"]),
             (as_int(fits[(p, "molecular", "tips")]["n_tips"]),) * 2)
        R.eq("%s: summary and pruned trees binary, placed-tip tree not" % p,
             (v_["n_internal"] == v_["n_tips"] - 1,
              v_["n_internal_pruned"] == v_["n_tips_pruned"] - 1,
              v_["n_internal_placed"] == v_["n_tips_placed"] - 1),
             (True, True, False))
        n_ = macro_name(p)
        R.macro(macros, "NodeAgesInternal" + n_, fmt_int(v_["n_internal"]))
        R.macro(macros, "NodeAgesUnpublished" + n_,
                fmt_int(v_["n_internal"] - v_["n_internal_published"]))
        R.macro(macros, "NodeAgesUnpublishedPct" + n_, fmt_pct(pct(
            v_["n_internal"] - v_["n_internal_published"], v_["n_internal"])))
        R.macro(macros, "NodeAgesPrunedInternal" + n_,
                fmt_int(v_["n_internal_pruned"]))
        R.macro(macros, "NodeAgesPrunedUnpublished" + n_, fmt_int(
            v_["n_internal_pruned"] - v_["n_internal_pruned_published"]))
        R.macro(macros, "NodeAgesPrunedUnpublishedPct" + n_, fmt_pct(pct(
            v_["n_internal_pruned"] - v_["n_internal_pruned_published"],
            v_["n_internal_pruned"])))
        R.macro(macros, "NodeAgesPlacedTips" + n_, fmt_int(v_["n_tips_placed"]))
        R.macro(macros, "NodeAgesPlacedInternal" + n_,
                fmt_int(v_["n_internal_placed"]))
        R.macro(macros, "NodeAgesPlacedPolyEvents" + n_,
                fmt_int(v_["n_polytomy_events_placed"]))
        R.macro(macros, "NodeAgesResolved" + n_,
                fmt_int(v_["n_internal_pruned"] - v_["n_internal_placed"]))
    WE_LEVEL_B = TEST_LEVEL / len(tested)
    WE_LR_CRIT = NormalDist().inv_cdf(1.0 - WE_LEVEL_B) ** 2
    R.macro(macros, "WENTests", fmt_int(len(tested)))
    R.macro(macros, "WENamesTested", name_list(tested))
    R.macro(macros, "WENamesNotTested", name_list([p for p in elig
                                                   if p not in tested]))

    def g(p, tree, conv, k):
        return float(fits[(p, tree, conv)][k])

    def overlap(p, a_, b_, v):
        return (g(p, a_[0], a_[1], v + "_lo") <= g(p, b_[0], b_[1], v + "_hi")
                and g(p, b_[0], b_[1], v + "_lo") <= g(p, a_[0], a_[1], v + "_hi"))
    M, MC, F, FM = (("molecular", "tips"), ("molecular", "full_misapplied"),
                    ("full", "tips"), ("molecular", "family"))
    n_ident = n_bound = n_rej = n_fg = n_fs = n_fo = n_famo = n_el = n_eh = 0
    wep = {}
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
        n_rej += p in tested and (not ident) and lr > WE_LR_CRIT
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
        if not ident and p in tested:
            # F3-15: p-value under the 50:50 mixture of zero and chi2(1)
            pv = mixture_p(lr)
            R.macro(macros, "WEP" + n, p_text(pv))
            wep[p] = pv
        R.macro(macros, "WERFull" + n, "%.3f" % rf)
        R.macro(macros, "WERMol" + n, "%.3f" % rm)
        R.macro(macros, "WERRatio" + n, "%.2f" % rr)
        R.macro(macros, "WEEpsFamily" + n, "%.2f" % g(p, *FM, k="epsilon"))
        R.macro(macros, "WERFamily" + n, "%.3f" % g(p, *FM, k="r"))
        # the edge of the restriction, 1 - rho_m/rho_f, and the profile
        # interval of the matching extinction fraction on the summary tree
        R.macro(macros, "WEEdge" + n, "%.2f" % (1.0 - rho_t / rho_c))
        R.macro(macros, "WEEpsTipsLo" + n, "%.2f" % g(p, *M, k="epsilon_lo"))
        R.macro(macros, "WEEpsTipsHi" + n, "%.2f" % g(p, *M, k="epsilon_hi"))
    R.macro(macros, "WENEligible", fmt_int(len(elig)))
    R.macro(macros, "WENIdentity", fmt_int(n_ident))
    R.macro(macros, "WENBoundary", fmt_int(n_bound))
    R.macro(macros, "WENRejected", fmt_int(n_rej))
    R.macro(macros, "RuleBonferroniLevel", "%.4f" % WE_LEVEL_B)
    R.macro(macros, "RuleBonferroniCritical", "%.2f" % WE_LR_CRIT)
    rej_names = sorted(p for p, v in wep.items() if v < WE_LEVEL_B)
    R.macro(macros, "WENamesRejected", name_list(rej_names) or "none")
    bind_ = sorted(wep)
    R.macro(macros, "WENamesTestedBinding", name_list(bind_) or "none")
    # F4-30 (D63): the family with the date-assignment summary trees added,
    # stated beside the family; the text says their statistic does not reach
    # its critical value
    asg_ = [p for p in elig if p not in tested and p not in poly_]
    n_wa = len(tested) + len(asg_)
    crit_wa = NormalDist().inv_cdf(1.0 - TEST_LEVEL / n_wa) ** 2
    rej_wa = sorted(p for p in asg_ if round(2.0 * (
        g(p, "molecular", "tips", "loglik")
        - g(p, "molecular", "full_misapplied", "loglik")), 4) > crit_wa)
    R.eq("no date-assignment summary tree rejects in the family with those "
         "trees added (as the text states)", rej_wa, [])
    R.macro(macros, "WENTestsWithAssigned", fmt_int(n_wa))
    R.macro(macros, "RuleBonferroniLevelWithAssigned", "%.4f" % (TEST_LEVEL / n_wa))
    R.macro(macros, "RuleBonferroniCriticalWithAssigned", "%.2f" % crit_wa)
    R.macro(macros, "WENamesRejectedWithAssigned", name_list(rej_wa) or "none")
    R.eq("no tested summary tree on which the restriction binds (the text "
         "says no test of the family can reject)", (bind_, rej_names), ([], []))
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
    with_poly, allr, nb = [], [], 0
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
        rr, lrs, em, bnd = [], [], [], []
        for sd_, f_ in sorted(poly[p].items()):
            R.eq("%s seed %s: the three uniform fits resolved" % (p, sd_), len(f_), 3)
            x_ = round(2.0 * (float(f_[M]["loglik"]) - float(f_[MC]["loglik"])), 4)
            e_ = float(f_[MC]["epsilon"])
            rr.append(round(float(f_[F]["r"]) / float(f_[M]["r"]), 4))
            lrs.append(x_)
            em.append(float(f_[M]["epsilon"]))
            bnd.append(e_ < 1e-6 and x_ > 0.0)
        allr += rr
        nb += all(bnd)
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
            ("outcome_resolved", "mu0" if all(bnd)
             else "rescaled" if not any(bnd) else "mixed")])
    R.macro(macros, "WEPolyNGroups", fmt_int(len(with_poly)))
    R.macro(macros, "WEPolyNames", name_list(with_poly))
    R.macro(macros, "WEPolyRRatioMin", "%.2f" % min(allr))
    R.macro(macros, "WEPolyRRatioMax", "%.2f" % max(allr))
    R.macro(macros, "WEPolyNBoundary", fmt_int(nb))
    # The fits repeated over the tree sets the sources archive
    # (registers/worked_example_tree_sets.csv; code/rebuild/worked_example_sets.R).
    # Each row is checked against the rules of the paper: the fractions are
    # tips over described species, a pruned tree keeps exactly the molecular
    # tips, the statistic is twice the log-likelihood difference, and the two
    # fits obey the likelihood's identity (same model inside the restriction,
    # mu = 0 only outside it).
    SET_TAG = {"posterior": ("Posterior", "posterior set"),
               "date_random": ("DateRandom", "date assignment, random dates"),
               "date_mean": ("DateMean", "date assignment, mean dates"),
               "prior": ("Prior", "root ages drawn from the prior"),
               "molecular_only": ("MolOnly", "source's tree of sequenced species")}
    set_reg = read_csv(path("registers", "worked_example_tree_sets.csv"))
    # The fits over archived tree sets (code/rebuild/worked_example_sets.R)
    # need the sources' posterior sets, which this repository does not
    # redistribute; a --fits directory without them is checked against the
    # committed copy in paper/.
    sets_csv = os.path.join(a.fits, "worked_example_bd_sets.csv")
    if not os.path.exists(sets_csv):
        print("  note: %s has no worked_example_bd_sets.csv; the tree-set fits "
              "are checked from paper/ (they need the archived tree sets of "
              "registers/worked_example_tree_sets.csv)" % a.fits)
        sets_csv = path("paper", "worked_example_bd_sets.csv")
    set_fit = read_csv(sets_csv)
    R.eq("tree-set fits: every row belongs to a register set",
         sorted({(r["partition"], r["set"]) for r in set_fit}),
         sorted({(x["partition"], x["set"]) for x in set_reg}))

    def med(v):
        v = sorted(v)
        k = len(v) // 2
        return v[k] if len(v) % 2 else (v[k - 1] + v[k]) / 2.0

    def rng2(lo_v, hi_v, f="%.2f"):
        lo_, hi_ = f % lo_v, f % hi_v
        return lo_ if lo_ == hi_ else "%s--%s" % (lo_, hi_)
    set_rows, set_sum = [], OrderedDict()
    for x in set_reg:
        p, st = x["partition"], x["set"]
        R.eq("tree set %s %s: eligible group, known set" % (p, st),
             p in elig and st in SET_TAG, True)
        d = float(bmap[p]["described"])
        rf_want = float(bmap[p]["labels_unique"]) / d
        rs = [r for r in set_fit if r["partition"] == p and r["set"] == st]
        R.eq("tree set %s %s: trees fitted" % (p, st), len(rs), int(x["n_select"]))
        bad = []
        for r in rs:
            n_ = int(r["n_tips"])
            rm_, rf_ = float(r["rho_matching"]), float(r["rho_full"])
            em_, ef_ = float(r["epsilon_matching"]), float(r["epsilon_full_rho"])
            lr_ = float(r["lr_statistic"])
            ok = (abs(rm_ - n_ / d) <= 1e-7 * rm_ and abs(rf_ - rf_want) <= 1e-7 * rf_
                  and abs(lr_ - max(0.0, 2.0 * (float(r["loglik_matching"])
                                                - float(r["loglik_full_rho"])))) <= 2e-3
                  and (x["tips"] != "molecular" or n_ == len(mol_rows[p])))
            edge = 1.0 - rm_ / rf_
            if em_ > edge + 1e-3:
                ok = ok and lr_ <= 2e-3 and abs((1 - ef_) - (1 - em_) * rf_ / rm_) <= 2e-3
            if ef_ < 1e-6 and lr_ > 2e-3:
                ok = ok and em_ < edge + 1e-3
            if not ok:
                bad.append(r["tree_index"])
        R.eq("tree set %s %s: fractions, tips, statistic and identity" % (p, st),
             bad, [])
        eps_ = [float(r["epsilon_matching"]) for r in rs]
        lrs_ = [float(r["lr_statistic"]) for r in rs]
        rr_ = [float(r["r_full_tree"]) / float(r["r_matching"])
               for r in rs if r["r_full_tree"]]
        # where the restriction binds is a property of each tree; the trees
        # are draws from one set, so no count above a critical value is formed
        # F4-26: binds where the fit with the full tree's fraction is on the
        # boundary and the statistic as written is positive, which is
        # exactly where the matching estimate lies below the edge
        nbnd_ = sum(1 for r in rs if float(r["epsilon_full_rho"]) < 1e-6
                    and float(r["lr_statistic"]) > 0.0)
        R.eq("tree set %s %s: binds exactly below the edge" % (p, st),
             sorted(r["tree_index"] for r in rs
                    if float(r["epsilon_full_rho"]) < 1e-6
                    and float(r["lr_statistic"]) > 0.0),
             sorted(r["tree_index"] for r in rs
                    if float(r["epsilon_matching"]) < 1.0
                    - float(r["rho_matching"]) / float(r["rho_full"])))
        set_sum[(p, st)] = (x["tips"], len(rs), lrs_, rr_)
        n = macro_name(p) + SET_TAG[st][0]
        R.macro(macros, "WESetN" + n, fmt_int(len(rs)))
        R.macro(macros, "WESetEpsMin" + n, "%.2f" % min(eps_))
        R.macro(macros, "WESetEpsMax" + n, "%.2f" % max(eps_))
        R.macro(macros, "WESetLRMin" + n, "%.1f" % min(lrs_))
        R.macro(macros, "WESetLRMax" + n, "%.1f" % max(lrs_))
        R.macro(macros, "WESetNBoundary" + n, fmt_int(nbnd_))
        R.macro(macros, "WESetLRMedian" + n, "%.2f" % med(lrs_))
        # the statistic over the set against each critical value: a
        # description of the set, not a test
        R.macro(macros, "WESetNAboveCrit" + n,
                fmt_int(sum(1 for x in lrs_ if x > LR_CRITICAL)))
        R.macro(macros, "WESetNAboveCorr" + n,
                fmt_int(sum(1 for x in lrs_ if x > WE_LR_CRIT)))
        if rr_:
            R.macro(macros, "WESetRRatioMin" + n, "%.2f" % min(rr_))
            R.macro(macros, "WESetRRatioMax" + n, "%.2f" % max(rr_))
        ft_, fc_ = fits[(p,) + M], fits[(p,) + MC]
        set_rows.append(OrderedDict([
            ("partition", p), ("set", SET_TAG[st][1]), ("n_trees", str(len(rs))),
            ("tips_pruned", "/".join(str(v) for v in sorted({int(r["n_tips"]) for r in rs}))),
            ("tips_full", "/".join(str(v) for v in sorted({int(r["n_tips_full_tree"])
                                                            for r in rs if r["n_tips_full_tree"]}))),
            ("epsilon_matching", "%.2f (%.2f--%.2f)" % (med(eps_), min(eps_), max(eps_))),
            ("lr", "%.2f (%s)" % (med(lrs_), rng2(min(lrs_), max(lrs_)))),
            ("n_boundary", str(nbnd_)),
            ("r_ratio", rng2(min(rr_), max(rr_)) if rr_ else ""),
            ("summary_epsilon_matching", "%.2f" % float(ft_["epsilon"])),
            ("summary_lr", "%.2f" % round(max(0.0, 2.0 * (float(ft_["loglik"])
                                                          - float(fc_["loglik"]))), 4))]))
    # posterior sets: the primary result for their groups, described by where
    # the statistic lies against the corrected critical value
    post_ = {p: v[2] for (p, st), v in set_sum.items() if st == "posterior"}
    R.macro(macros, "WEPostNames", name_list(sorted(post_)))
    R.macro(macros, "WESetNamesAboveAll", name_list(sorted(
        p for p, v in post_.items() if min(v) > WE_LR_CRIT)) or "none")
    R.macro(macros, "WESetNamesStraddle", name_list(sorted(
        p for p, v in post_.items() if min(v) <= WE_LR_CRIT < max(v))) or "none")
    R.macro(macros, "WESetNamesBelowAll", name_list(sorted(
        p for p, v in post_.items() if max(v) <= WE_LR_CRIT)) or "none")
    nopost_ = [p for p in elig if p not in post_]
    R.macro(macros, "WENoPostN", fmt_int(len(nopost_)))
    R.macro(macros, "WENoPostNames", name_list(sorted(nopost_)))
    R.macro(macros, "WENoPostNRejected", fmt_int(sum(p in rej_names for p in nopost_)))
    grp_ = sorted({p for p, _ in set_sum})
    R.macro(macros, "WESetNGroups", fmt_int(len(grp_)))
    R.macro(macros, "WESetNamesGroups", name_list(grp_))
    allrr_ = [v for t in set_sum.values() for v in t[3]]
    R.macro(macros, "WESetRRatioMin", "%.2f" % min(allrr_))
    R.macro(macros, "WESetRRatioMax", "%.2f" % max(allrr_))
    # crown age over the root height of the fitted tree: reported, not
    # applied; every tree is fitted with its archived branch lengths
    sf = [float(crown_csv[x["partition"]]["crown_age_ma"])
          / facts[bmap[x["partition"]]["canonical_tree"]]["root"] for x in sel]
    R.macro(macros, "WEScaleMin", "%.2f" % min(sf))
    R.macro(macros, "WEScaleMax", "%.2f" % max(sf))
    for x in sel:
        R.eq("%s: fitted with the archived branch lengths (factor 1)" % x["partition"],
             float(x["time_scale_factor"]), 1.0)

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
    # larger trees of a group or a subclade the search missed, not counted
    larger = defaultdict(list)
    for x in read_csv(path("registers", "larger_trees_not_counted.csv")):
        larger[x["partition"]].append(x)
    R.macro(macros, "TabWhichTreeLargerNote", "; ".join(
        "%s: %s" % (p, " and ".join(
            "\\citet{%s}, %s%s %s of %s%s" % (
                x["bibkey"], x["n_qualifier"] + " " if x["n_qualifier"] else "",
                fmt_int(x["n"]), x["unit"], x["clade"],
                ", dated" if x["dated"] == "yes" else "") for x in xs))
        for p, xs in sorted(larger.items())))
    # F4-07: sizes of a larger tree and of its molecular source tree
    for xs in larger.values():
        for x in xs:
            if x.get("n_molecular_source"):
                k_ = "".join(t[:1].upper() + t[1:] for t in
                             re.split(r"[^A-Za-z]+", x["bibkey"]) if t)
                R.macro(macros, "LargerN" + k_, fmt_int(x["n"]))
                R.macro(macros, "LargerNMol" + k_, fmt_int(x["n_molecular_source"]))
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
            ("uncertainty_status", "" if not chron[p] else "none" if not kind
             else "placement" if not avail_of[p] else "set"
             if kind in SET_KINDS else "summary"),
            ("rho_full", rho_fmt(len(tips[p]) / float(d))),
            ("rho_molecular", rho_fmt(mt / float(d))), ("described", d),
            ("study_dating", sdp),
            ("molecular_only_tree", "yes" if p in molonly else "no"),
            ("molecular_only_dated", ("yes" if modat[p] == "time" else "no")
             if p in molonly else ""),
            ("molecular_share", "%.2f" % (mt / float(len(tips[p])))),
            ("tips_exceed_named", "yes" if len(tips[p]) > (1 + TIPS_OVER_NAMED)
             * len(c0[p]) else "no"),
            ("larger_tree_not_counted", ";".join(x["bibkey"] for x in larger.get(p, []))),
        ])
    table_eq("table_which_tree.csv", t1, "partition")

    t2 = OrderedDict()
    for lv, label in zip(LEVELS, ("named", "molecular", "dated",
                                  "uncertainty available")):
        t2[label] = OrderedDict([
            ("level", label), ("species", fmt_int(lad[lv]["union"])),
            ("pct", fmt_pct(lad[lv]["pct_union"])),
            ("rep", fmt_pct(rep[lv]["pct_union"])), ("col", fmt_pct(col[lv]["pct_union"])),
            ("union", fmt_pct(pct(lad["C0"]["sum"] + n_add, den)) if lv == "C0"
             else "<= " + fmt_pct(pct(lad["C1"]["sum"] + n_add, den)) if lv == "C1" else ""),
            ("gbif", fmt_pct(pct(lad["C0"]["sum"] + len(added), den)) if lv == "C0" else ""),
            ("product", fmt_pct(pct(round(prod[lv], 4), den)) if lv in prod else ""),
            ("nobirds", fmt_pct(nob[lv]["pct_union"])),
            ("seqdating", fmt_pct(seq_lad[lv]["pct_union"]) if lv == "C1u" else ""),
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
             + (" (rejected)" if p in tested and (not ident) and lr > WE_LR_CRIT
                else "")),
            ("r_full", f3(ff["r"])), ("r_full_lo", f3(ff["r_lo"])), ("r_full_hi", f3(ff["r_hi"])),
            ("r_molecular", f3(ft["r"])), ("r_molecular_lo", f3(ft["r_lo"])),
            ("r_molecular_hi", f3(ft["r_hi"])),
            ("eps_family", f2(fm["epsilon"])), ("eps_family_lo", f2(fm["epsilon_lo"])),
            ("eps_family_hi", f2(fm["epsilon_hi"])),
            ("r_family", f3(fm["r"])), ("r_family_lo", f3(fm["r_lo"])),
            ("r_family_hi", f3(fm["r_hi"])),
            ("basis", "summary"), ("n_trees", ""), ("n_boundary", ""),
            ("lr_lo", ""), ("lr_hi", ""),
        ])
        # A group whose source archives a posterior set, else a set with
        # dates drawn at random, is reported in Table 3 as the median and
        # range over that set; its summary-tree fits are in the Supplement.
        st_ = next((x_ for x_ in ("posterior", "date_random")
                    if any(r_["partition"] == p and r_["set"] == x_
                           for r_ in set_fit)), None)
        if st_:
            rs_ = [r_ for r_ in set_fit if r_["partition"] == p and r_["set"] == st_]
            row_ = t3[p]
            for key_, col_, ff_ in (("eps_molecular", "epsilon_matching", f2),
                                    ("eps_fullrho", "epsilon_full_rho", f2),
                                    ("r_full", "r_full_tree", f3),
                                    ("r_molecular", "r_matching", f3),
                                    ("lr", "lr_statistic", f2)):
                xs_ = [float(r_[col_]) for r_ in rs_]
                row_[key_] = ff_(med(xs_))
                row_[key_ + "_lo"], row_[key_ + "_hi"] = ff_(min(xs_)), ff_(max(xs_))
            lam_ = [(float(r_["r_full_rho"]) / (1 - float(r_["epsilon_full_rho"])))
                    / (float(r_["r_matching"]) / (1 - float(r_["epsilon_matching"])))
                    for r_ in rs_]
            row_["lambda_ratio"] = "%.2f" % med(lam_)
            nb_ = sum(1 for r_ in rs_ if float(r_["epsilon_full_rho"]) < 1e-6
                      and float(r_["lr_statistic"]) > 0.0)
            row_["classification"] = ("boundary" if nb_ == len(rs_) else
                                      "identity" if nb_ == 0 else "mixed")
            row_["basis"], row_["n_trees"], row_["n_boundary"] = (
                st_, str(len(rs_)), str(nb_))
    table_eq("table_worked_example.csv", t3, "partition")

    R.section("Figure source data (paper/fig*_*.csv)")
    casc = read_csv(path("paper", "fig1_coverage_cascade.csv"))
    got = {("level", r["level"], r["group"]): by_cat[r["level"]][r["group"]]
           for r in casc if r["kind"] == "level"}
    got.update({("loss", r["level"], r["group"]): lost[r["group"]]
                for r in casc if r["kind"] == "loss"})
    R.eq("Figure 1: each loss sits between the levels its reason names",
         sorted((r["group"], r["level"]) for r in casc if r["kind"] == "loss"),
         sorted((r["group"], LOSS_STEP[r["group"]]) for r in casc
                if r["kind"] == "loss"))
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
    n_ovr = sum(r["override"] for r in got.values())
    sm("SuppOverridesResolvedN", fmt_int(n_ovr))
    sm("SuppOverridesUnresolvedN", fmt_int(
        len(read_csv(path("standardized", "name_overrides.csv"))) - n_ovr))
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
            ("exact_label_tips", ex),
            ("source_labels", src_counts[p][0] if p in src_counts else ""),
            ("source_merged", src_counts[p][1] if p in src_counts else ""),
            ("source_unmatched", src_counts[p][2] if p in src_counts else ""),
            ("molecular_species", len(c1[p])),
            ("product", fix(len(c0[p]) * float("%.6f" % fmol[p]), 1)),
            ("f_low", "%.3f" % float(b["f_low"])), ("f_high", "%.3f" % float(b["f_high"])),
            ("molecular_only_tree", "yes" if p in molonly else "no"),
            ("molecular_only_dated", ("yes" if modat[p] == "time" else "no")
             if p in molonly else "")])
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
        status = ("none" if not kind else "counted" if counts[p] else
                  "placement" if not avail_of[p] else why_of[p])
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
            ("prior_q_lo", prior_quantile(calreg[p]["prior_density"],
                                          float(iv["interval_lo_ma"]))
             if kind and calreg.get(p, {}).get("prior_density") else ""),
            ("prior_q_center", prior_quantile(
                calreg[p]["prior_density"], float(iv["interval_center_ma"])
                if iv["interval_center_ma"] else float(iv["source_point_height_ma"]))
             if kind and calreg.get(p, {}).get("prior_density") else ""),
            ("prior_q_hi", prior_quantile(calreg[p]["prior_density"],
                                          float(iv["interval_hi_ma"]))
             if kind and calreg.get(p, {}).get("prior_density") else ""),
            ("parent_calibration", calreg.get(p, {}).get("parent_calibration", "")
             if kind else ""),
            ("child_calibration", calreg.get(p, {}).get("child_calibration", "")
             if kind else ""),
            ("prior_range", ("calibration" if calreg.get(p, {}).get("calibration")
                             == "two_sided" else "band"
                             if calreg.get(p, {}).get("band_lo_ma") else "")
             if kind else ""),
            ("status", status)])
        for k_ in ("prior_q_lo", "prior_q_center", "prior_q_hi"):
            R.eq("%s: registers/crown_intervals.csv %s" % (p, k_),
                 iv.get(k_, "") if kind else "", got[p][k_])
    supp_eq("supp_crown", got, "partition")
    R.eq("crown interval widths = hi - lo of registers/crown_intervals.csv",
         [p for p in got if intervals.get(p) and abs(
             float(intervals[p]["interval_width_myr"])
             - (float(intervals[p]["interval_hi_ma"]) - float(intervals[p]["interval_lo_ma"])))
          > 1e-4], [])
    for key, nm in (("degenerate", "Degenerate"), ("at_bound", "AtBound"),
                    ("prior_reproduced", "Prior"), ("prior_sampled", "PriorSampled"),
                    ("kind", "Kind"), ("placement", "Placement")):
        ps = sorted(p for p, r in got.items() if r["status"] == key)
        sm("SuppCrown" + nm, name_list(ps) if ps else "none")
    sm("SuppCrownNone", name_list(sorted(p for p, r in got.items() if r["status"] == "none")))
    tac_ = read_csv(path("registers", "taxon_addition_set_check.csv"))
    R.eq("one taxon-addition set is checked", len(tac_), 1)
    sm("SuppTaxonAdditionNTrees", fmt_int(as_int(tac_[0]["n_trees"])))
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
        if lv != "C0" and level_sets(bp)[i]:
            n_ = lad[lv]["sum"] - len(level_sets(bp)[i]) * (1 - scale)
        low[lv] = pct(n_, den)
    got["bird_band_low"] = srow("bird_band_low", den, low)
    got["union"] = srow("union", den, OrderedDict(
        [("C0", pct(lad["C0"]["sum"] + n_add, den)),
         ("C1", pct(lad["C1"]["sum"] + n_add, den)), ("C1d", None), ("C1u", None)]))
    got["c1u_sequence_dating"] = srow("c1u_sequence_dating", den, OrderedDict(
        [("C0", None), ("C1", None), ("C1d", None),
         ("C1u", seq_lad["C1u"]["pct_union"])]))
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
            ("molecular_keys", r["mkeys"]), ("ot_keys", r["okeys"]),
            ("keys_both", r["mand"]),
            ("keys_both_pct", fix(pct(r["mand"], r["mkeys"]), 1) if r["mkeys"] else ""),
            ("keys_both_ot_pct", fix(pct(r["mand"], r["okeys"]), 1) if r["okeys"] else "")])
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
            ("r_first", "%.5f" % float(fo["first"]["r"])),
            ("r_last", "%.5f" % float(fo["last"]["r"])),
            ("epsilon_first", "%.4f" % float(fo["first"]["epsilon"])),
            ("epsilon_last", "%.4f" % float(fo["last"]["epsilon"]))])
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
    sm("SuppFamOrderMaxEpsDiff", "%.4f" % max(
        abs(float(v["last"]["epsilon"]) - float(v["first"]["epsilon"]))
        for v in fo.values()))
    sm("SuppFamOrderMaxDLL", "%.1f" % max(
        abs(float(v["first"]["loglik_last_at_first_mle"]) - float(v["first"]["loglik"]))
        for v in fo.values()))
    supp_eq("supp_worked_polytomy", poly_rows, "partition")
    want_sets = supp_rows("supp_worked_sets")
    R.eq("supp_worked_sets: rows, in order",
         [(r["partition"], r["set"]) for r in set_rows],
         [(r["partition"], r["set"]) for r in want_sets])
    for r, w_ in zip(set_rows, want_sets):
        R.eq("supp_worked_sets: %s %s" % (r["partition"], r["set"]),
             {c: (v, w_.get(c)) for c, v in r.items() if v != w_.get(c)}, {})
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
    # ERRATUM.md states the counts of the same register
    reg_ = read_csv(path("registers", "restored_branch_lengths.csv"))
    wd_ = {os.path.basename(r.get("filename") or r.get("file") or "")
           for r in read_csv(path("standardized", "withdrawn.csv"))}
    tr_ = [r for r in reg_ if r["file"].startswith("Data/standardized/trees/")]
    st_ = [r for r in reg_ if r["file"].startswith("Data/site/")]
    with open(path("ERRATUM.md"), encoding="utf-8") as fh:
        et_ = re.sub(r"\s+", " ", fh.read())
    m_ = re.search(r"Of its (\d+) rows, (\d+) have `changed` = `yes`: (\d+) are "
                   r"trees redistributed in `standardized/trees/`, (\d+) is "
                   r"`timetree.nwk`.*? and (\d+) are the copies of trees under "
                   r"`site/data/partitions/`, of which the browser site holds (\d+)\..*? "
                   r"The other (\d+) rows are files of "
                   r"`standardized/trees/` that were checked and found unchanged", et_)
    R.eq("ERRATUM.md: rows, changed, redistributed trees, withdrawn, site copies, "
         "site copies held, unchanged = registers/restored_branch_lengths.csv",
         tuple(int(x) for x in m_.groups()) if m_ else None,
         (len(reg_), sum(r["changed"] == "yes" for r in reg_),
          sum(r["changed"] == "yes" and os.path.basename(r["file"]) not in wd_
              and os.path.basename(r["file"]) in files for r in tr_),
          sum(r["changed"] == "yes" and os.path.basename(r["file"]) in wd_ for r in tr_),
          sum(r["changed"] == "yes" for r in st_),
          sum(r["changed"] == "yes" and os.path.exists(path(*r["file"].split("/")[1:]))
              for r in st_),
          sum(r["changed"] != "yes" and os.path.basename(r["file"]) in files for r in tr_)))
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

    check_site(R, macros, basis, head)

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
