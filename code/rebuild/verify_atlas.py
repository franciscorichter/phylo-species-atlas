#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verify_atlas.py -- THE ACCEPTANCE TEST FOR THE PHYLO-ATLAS DEPOSIT.

Exit status 0 iff every check passes.  Any failure exits non-zero.

    python3 scripts/verify_atlas.py               # full run, every shipped tree
    python3 scripts/verify_atlas.py --quick       # canonical trees only
    python3 scripts/verify_atlas.py --json r.json # machine-readable report
    python3 scripts/verify_atlas.py --list        # list the checks and exit


WHY THIS FILE EXISTS
====================
The published 23.7 % coverage headline was a hard-coded literal
(Paper/code/figures/fig2_sensitivity_bounds.py:61) reachable from no
numerator/denominator pair anywhere in the deposit; its implied denominator,
2,084,409, exists in no file.  It survived three submissions because nothing
ever tried to recompute it.  This script is the thing that tries.

THE INDEPENDENCE RULE.  This script NEVER imports scripts/rebuild_atlas.py and
never trusts a number it produced.  It re-derives every headline from

    Data/_rebuild/truth_{dated,tips,species,provenance}.csv   measured ground truth
    Data/standardized/trees/*.nwk                             the shipped files
    Data/standardized/{dictionary,full_mapping}.csv           the name layer
    Data/raw/**                                               the source analyses
    Data/tree_provenance.csv                                  pre-audit ledger,
                                                              annotation fields only

and then asserts that Paper/inputs/headline_numbers.json, Data/coverage_basis.csv,
Data/molecular_tips.csv and the figure-input CSVs agree with what it measured.
coverage_basis.csv is read as a CLAIM UNDER TEST, not as an input.

The ONE thing it does import is the audit's recursive Newick parser at
Paper/audit-2026-08-20/code/common.py (`parse_newick` / `leaf_depths`).  A naive
one-pass parser gets ultrametricity wrong, so the reference implementation is
reused rather than reinvented -- reimplementing it here would be the defect this
script exists to catch.

WHAT IT CHECKS  (each is a named check; --list prints them)
    INPUTS        every required file exists
    TREES         re-measure every shipped .nwk; truth_dated/truth_tips must match
    NAMES         S05: per partition, labels resolving to a Catalogue of Life
                  accepted name (col_resolution.csv, COL26.9 base), re-derived
    DECISIONS     D4 one tip per distinct species label; D7 no pruned tip
                  survives and the checked/unchecked record is re-derived; D8
                  each deposited chronogram's coverage re-measured
    BRANCHLEN     every shipped tree reproduces the reference declared for it in
                  Data/standardized/tree_references.csv (catches rewritten
                  branch lengths, which the CV filter cannot)
    FILE-TRUTH    dated <=> ultrametric in every ledger of record, both directions;
                  no NaN branch length in a canonical tree carried as dated
    PREAUDIT      the pre-audit ledgers disagree with geometry -- that is the
                  paper's FINDING, so it is reported, not failed; what IS failed
                  is truth_dated.csv mis-describing that disagreement
    INVARIANTS    described <= estimated_high; estimated_low <= estimated_high;
                  tips <= described; f_low <= f_mol <= f_high; f in [0,1]
    BASIS         coverage_basis.csv's 62-row deposit contract, column by column
    FMOL          every f_mol re-measured from the source tree; every f_basis
                  re-derived from the deposit's own decision rule
    HEADLINES     every headline recomputed independently and asserted equal
    PERCENTS      every percentage in headline_numbers.json is reproducible from
                  a numerator and a denominator printed beside it
    MONOTONICITY  a stricter setting never increases coverage, per partition,
                  per rung, and across every sensitivity variant
    MOLTIPS       molecular_tips.csv enumerates exactly the measured tip sets
    UNITS         "species" at C0 only; "tip-equivalents" at C1 and below
    CROSSLEDGER   no field holds two different values anywhere in the deposit,
                  unless the ground truth documents the disagreement
    NOLITERALS    no coverage number is hard-coded in a figure script
    NOABSPATHS    no absolute home-directory path under Paper/code/ or
                  scripts/ -- the deposit must run from a clone anywhere
    RETIRED       no retired figure (23.7 %, 8.65 %, 2,084,409, 494,005 as a
                  Table-1 numerator, "17/16/18 dated trees") survives in live
                  figure code or live manuscript sources, except inside an
                  explicit retirement disclosure
    SUPPLEMENT    the Supplement has one section per Supplement number the
                  article cites, in order of first citation; every supp_*.csv
                  table and supp_macros.tex value is re-derived; each reference
                  tree's archived source reproduces the branch lengths of the
                  source file Data/reference_tree_sources.csv names
    RELEASE       release/path_map.csv names an existing file for every file of
                  the public data release; every path the paper prints, and
                  every generated file it inputs, is in the release; the
                  release's citation metadata carries the paper's title and
                  the release named in the text
    BIORXIV       the bioRxiv wrappers read the revised body, metadata and
                  Supplement, are named, carry no review-only element, and
                  print the data and code statements

CONVENTIONS THIS SCRIPT ENFORCES, stated so a failure can be read
    * a chronogram is a tree whose root-to-tip CV is <= CHRONOGRAM_TOLERANCE;
      the ledger's `dated` column is never a source of truth
    * f_mol is published to FRAC_DP decimal places and the PUBLISHED value is
      the operative one, so re-measurements are rounded to FRAC_DP before
      comparison -- the deposited file must reproduce the published numbers
      exactly, not to within a rounding error
    * C0 is quoted on the UNION convention; the per-partition SUM is checked
      alongside so the two can never drift
    * prokaryotes never enter a eukaryotic numerator or denominator

NOT CHECKED HERE, ON PURPOSE
    * licensing.  timetree.nwk is not redistributed (decision D2); it is not
      among the shipped trees, and TREES fails if it reappears without a row
      in truth_dated.csv.
    * the MEE article type.
"""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.util
import json
import math
import os
import ast
import re
import shutil
import sys
import threading
import traceback
from collections import Counter, OrderedDict, defaultdict

# The reference parser is recursive and the largest shipped tree has ~356k
# leaves.  Same limit the generator uses; the work is additionally run on a
# thread with a large stack where the platform allows it.
sys.setrecursionlimit(300000)

# ==========================================================================
# PATHS -- all derived from this file's own location.  NO ABSOLUTE PATH may
# appear in this file: it is itself inside the scope of the NOABSPATHS check.
# ==========================================================================

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir))


def P(*parts):
    return os.path.join(ROOT, *parts)


def rel(path):
    try:
        return os.path.relpath(path, ROOT)
    except ValueError:
        return path


TRUTH = OrderedDict([
    ("dated", P("Data", "_rebuild", "truth_dated.csv")),
    ("tips", P("Data", "_rebuild", "truth_tips.csv")),
    ("species", P("Data", "_rebuild", "truth_species.csv")),
    ("provenance", P("Data", "_rebuild", "truth_provenance.csv")),
])

AUDIT_COMMON = P("Paper", "audit-2026-08-20", "code", "common.py")
STD_TREES = P("Data", "standardized", "trees")
DICTIONARY = P("Data", "standardized", "dictionary.csv")
FULL_MAPPING = P("Data", "standardized", "full_mapping.csv")
STD_METADATA = P("Data", "standardized", "metadata.csv")
TREE_PROVENANCE = P("Data", "tree_provenance.csv")
DATA_ESTIMATES = P("Data", "data_estimates.csv")
DATA_PROVENANCE = P("Data", "data_provenance.csv")
PARTITION_SUMMARY = P("Data", "site", "partition_summary.csv")

BASIS = P("Data", "coverage_basis.csv")
MOLTIPS = P("Data", "molecular_tips.csv")
HEADLINE = P("Paper", "inputs", "headline_numbers.json")
COVERAGE_PP = P("Paper", "inputs", "coverage_per_partition.csv")
FIG1_BARS = P("Paper", "inputs", "fig1_coverage_bars.csv")
FIG2_SCATTER = P("Paper", "inputs", "fig2_coverage_scatter.csv")
FIG3_BARS = P("Paper", "inputs", "fig3_c0c1_bars.csv")
FIG3_COUNTS = P("Paper", "inputs", "fig3_dated_undated_counts.csv")

PRE_AUDIT = P("Data", "_archive", "2026-08-20-pre-audit")

# Per-tip molecular evidence (f_basis per_tip_evidence).
MOL_TIP_EVIDENCE = P("Data", "molecular_tip_evidence.csv")
MOL_EVIDENCE_KINDS = ("source_tree_label", "reference_genome",
                      "analysis_terminal")

# Revision 2026-10 registers: off-target checks and pruned tips (D7), deposited
# chronograms (D8).
OFFTARGET_CHECKS = P("Data", "offtarget_checks.csv")
PRUNED_TIPS = P("Data", "pruned_tips.csv")
COL_RESOLUTION_NONCANON = P("Data", "standardized", "col_resolution_noncanonical.csv")
COL_DENOMINATORS = P("Data", "col_denominators.csv")
OT_COUNTS = P("Data", "external", "opentree16.1_counts.csv")
OT_SPECIES = P("Data", "external", "opentree16.1_phylo_species.tsv.gz")
OT_CSV = P("Paper", "inputs", "opentree_comparison.csv")
OT_SPECIES_COL = P("Data", "external", "opentree16.1_phylo_species_col.csv")
# R1-06 and R1-09 registers
MOLECULAR_ONLY_TREES = P("Data", "molecular_only_trees.csv")
STUDY_DATING_VALUES = {
    "absolute_deposited": "absolute chronogram deposited",
    "absolute_deposited_fewer_tips": "dated; chronogram lacks tips",
    "absolute_not_deposited": "dated; chronogram not archived",
    "absolute_not_read": "dated; archive not readable",
    "relative_only": "relative time only",
    "not_dated": "not dated"}
STUDY_DATING_MACRO = {"absolute_deposited_fewer_tips": "FewerTips",
                      "absolute_not_deposited": "NotDeposited",
                      "absolute_not_read": "NotRead",
                      "relative_only": "Relative", "not_dated": "NotDated"}
WORKED_FITS = P("Paper", "inputs", "worked_example_bd.csv")
WORKED_SEL = P("Paper", "inputs", "worked_example_bd_selection.csv")
WORKED_CHECK = P("Paper", "inputs", "worked_example_bd_crosscheck.csv")
WORKED_TOL = 1e-6      # relative, on sampling fractions written with 8 s.f.
WORKED_FAMILY_FRACTIONS = P("Paper", "inputs",
                            "worked_example_bd_family_fractions.csv")
COL_FAMILY_COUNTS = P("Data", "standardized", "col269_family_species_counts.csv")
WORKED_UNIFORM_TOL = 1e-6      # clade likelihood at equal f == uniform one
WORKED_CROSSCHECK_RTOL = 1e-4  # rates against phytools, relative (Supplement S9)
WORKED_EPS_BOUNDARY = 1e-6     # = EPS_BOUNDARY in worked_example_bd.R
WORKED_IDENTITY_LL_TOL = 1e-3  # log-likelihoods are written at 4 dp
WORKED_IDENTITY_RATIO_RTOL = 1e-4
WORKED_LR_CUTOFF = 2.705543    # qchisq(0.90, 1): 5% point of the 50:50 chi2(0):chi2(1) mixture (D27)
WORKED_TEST_LEVEL = 0.05
WORKED_RESOLVED = ("accepted", "provisional", "synonym", "accepted_orthographic",
                   "provisional_orthographic", "synonym_orthographic")
CROWN_SOURCES = P("Data", "crown_interval_sources.csv")
CROWN_INTERVALS = P("Data", "crown_intervals.csv")
CROWN_AGES = P("Paper", "inputs", "crown_ages.csv")
# S02.  Which interval kinds make an interval recoverable (rung C1b); the
# table is restated here, not imported, so the two scripts check each other.
CROWN_KIND_COUNTS = {
    "node_hpd_annotation": "yes", "node_label_interval": "yes",
    "posterior_sample": "yes", "bootstrap_dating_sample": "yes",
    "calibration_prior_sample": "no", "taxon_addition_sample": "no",
    "date_assignment_sample": "no"}
CROWN_SUMMARY_KINDS = ("node_hpd_annotation", "node_label_interval")
CROWN_POINT_TOL = 1e-3       # Myr
CROWN_DEGENERATE_MYR = 1e-3  # D5: an interval narrower than this is fixed
CROWN_BOUND_TOL_MYR = 1e-2   # D22: a centre this close to an end sits on a bound
CROWN_CENTER_KINDS = ("node_hpd_annotation", "posterior_sample",
                      "bootstrap_dating_sample", "calibration_prior_sample")
CROWN_CLADES = P("Data", "crown_clades.csv")
CROWN_VALUE_TOL = 1e-4       # Ma; crown_ages.csv is written at 4 dp
CALIBRATION_REGISTER = P("Data", "calibration_register.csv")
PRIOR_CONTAINMENT_MIN = 0.9  # D25: interval inside the prior's 95% range ...
PRIOR_WIDTH_RATIO_MIN = 0.5  # ... and not even halved by the data
CAL_CLASSES = ("two_sided", "one_sided", "none", "not_checked")


def d25_prior_test(lo, hi, cal):
    """D25, implemented independently of crown_intervals.py: (class, p_lo,
    p_hi, containment, width_ratio, reproduced)."""
    if not cal:
        return "", None, None, None, None, False
    k = cal.get("calibration")
    if k != "two_sided":
        return k, None, None, None, None, False
    plo, phi = float(cal["prior_lo_ma"]), float(cal["prior_hi_ma"])
    ov = min(hi, phi) - max(lo, plo)
    ov = ov if ov > 0 else 0.0
    cont = ov / (hi - lo) if hi > lo else 1.0
    ratio = (hi - lo) / (phi - plo)
    return k, plo, phi, cont, ratio, (cont >= PRIOR_CONTAINMENT_MIN
                                      and ratio >= PRIOR_WIDTH_RATIO_MIN)
DATING_DEPOSIT_CHECKS = P("Data", "dating_deposit_checks.csv")

# S05 single-authority name resolution (Catalogue of Life).
COL_RESOLVED_MATCH_TYPES = ("exact", "synonym")
C0_NAME_RESOLUTION = P("Paper", "inputs", "c0_name_resolution.csv")

# Files withdrawn from distribution (decision D2 and any later one).
WITHDRAWN = P("Data", "standardized", "withdrawn.csv")
# S05 / D18 and D14 registers.
COL_RESOLUTION = P("Data", "standardized", "col_resolution.csv")
COL_RESOLUTION_RESOLVED = ("accepted", "provisional", "synonym",
                           "accepted_orthographic",
                           "provisional_orthographic",
                           "synonym_orthographic")
NAME_OVERRIDES = P("Data", "standardized", "name_overrides.csv")
TIME_SCALE_CHECKS = P("Data", "time_scale_checks.csv")
# D14 plausibility: a tree declared in absolute time (Ma) must be older than
# this at its root; relative-time trees are conventionally scaled to 1.
MIN_ABSOLUTE_ROOT_MA = 1.0

# Every shipped tree declares the file its branch lengths must reproduce.
TREE_REFERENCES = P("Data", "standardized", "tree_references.csv")
PRE_REMAP_DICT = P("Paper", "archived", "wrong_root_standardized_20260512_012852",
                   "dictionary.csv")

# Scopes for the source-hygiene checks.
CODE_SCOPES = [P("Paper", "code"), P("scripts")]
FIGURE_SCOPES = [P("Paper", "code", "figures"), P("Paper", "scripts")]
LIVE_TEX_SCOPES = [P("Paper", "shared"), P("Paper", "MEE"), P("Paper", "bioRxiv")]

# Anything under a directory whose name matches this is history, not the
# deposit.  Never scanned, never modified, never "fixed".
EXCLUDED_DIR_RE = re.compile(
    r"^(_.*|\..*|__pycache__|archived?|backups?|old|deprecated|superseded"
    r"|audit-\d{4}-\d{2}-\d{2}|forks|papers|reports"
    r"|.*-\d{4}-\d{2}-\d{2}|data-audit.*|expert-panel.*|pnas-review.*"
    r"|methods-rewrite.*|system\d-.*|figure.*rework.*|figure\d.*-\d{4}.*)$",
    re.I)

# The generator is exempt from the TEXT scans and only from those: it is
# REQUIRED to name what it retired, and every number it emits is independently
# recomputed by the HEADLINES check, which is a far stronger guarantee than a
# grep could ever be.
TEXT_SCAN_EXEMPT = {os.path.abspath(P("scripts", "rebuild_atlas.py")),
                    os.path.abspath(__file__)}

# ==========================================================================
# RULE PARAMETERS.  Rules, not reported quantities.  No coverage figure
# appears anywhere in this file.
# ==========================================================================

YES = {"yes", "true", "1"}
NO = {"no", "false", "0"}

# An absolute home-directory prefix, assembled rather than written out, so
# that THIS FILE is inside the scope of its own NOABSPATHS check and passes
# it. A verifier exempt from its own rule is a verifier nobody trusts.
HOME_PATH_MARKERS = (os.sep + "Users" + os.sep, os.sep + "home" + os.sep)

CHRONOGRAM_TOLERANCE = 0.05      # root-to-tip CV at or below this => chronogram
FRAC_DP = 6                      # published precision of f_mol / roottotip_cv
PCT_DP = 4                       # published precision of percentages in CSVs
ELASTICITY_FRACTION = 0.20       # denominator elasticity demonstration
# Round F1 rule parameters (declared here as in rebuild_atlas.py)
C1B_WIDTH_RULE_RATIO = 1.0       # interval/prior width at which the prior sets it
TABLE1_TIPS_OVER_NAMED = 0.10    # Table 1 flag: tips exceed named species by this
ZERO_BRANCH_TOL = 1e-6           # = ZERO_BRANCH_TOL in worked_example_bd.R
WORKED_POLY = P("Paper", "inputs", "worked_example_bd_polytomy.csv")
WORKED_FORDER = P("Paper", "inputs", "worked_example_bd_family_order.csv")
CROWN_SET_KINDS = ("posterior_sample", "bootstrap_dating_sample",
                   "calibration_prior_sample", "taxon_addition_sample",
                   "date_assignment_sample")

SOURCE_NUMERATOR_TOL = 0.02      # |measured - recorded| / recorded
SOURCE_NONSUBSET_TOL = 0.05      # |source \ shipped| / |source|
SOURCE_SCAN_MAX_BYTES = 64 * 1024 * 1024
SOURCE_TREE_EXT = (".tre", ".nwk", ".newick", ".tree", ".taxa.txt")
TAXA_LIST_EXT = ".taxa.txt"

CV_RECHECK_TOL = 1e-3            # re-measured CV vs truth_dated.csv
ATOL_PCT = 1e-6                  # percentages
ATOL_SUM = 5e-4                  # weighted numerator sums
ATOL_FRAC = 5e-7                 # fractions published at FRAC_DP
BL_REL_TOL = 1e-9                # branch length vs its declared reference

CATEGORY_ORDER = ["Vertebrates", "Plants", "Arthropods",
                  "Other animals", "Microbes & protists"]

RUNGS = ("C0", "C1", "C1d", "C1b")

ASSEMBLY_AUTOMATED = "automated_synthesis"
ASSEMBLY_COMPILED = "compiled"
ASSEMBLY_PRIMARY = "primary_analysis"
ASSEMBLY_OTL_TRIGGER = "full_pipeline"
ASSEMBLY_INFERENCE_TRIGGERS = ("synthesis", "literature_compilation",
                               "supertree_re_dated")

RAW_PATH_RE = re.compile(r"(?:raw|molecular_sources)/[A-Za-z0-9_./()+\-]+")
RATIO_RE = re.compile(r"(?<![\d.])(\d[\d,]{1,9})\s*/\s*(\d[\d,]{1,9})(?![\d.])")
FRACTION_RE = re.compile(r"(?<![\d.])(0?\.\d+|1\.0+)(?![\d])")
SENTENCE_SPLIT_RE = re.compile(r"(?<=[.;])\s+")


# ==========================================================================
# small utilities
# ==========================================================================

def truthy(v):
    return str(v).strip().lower() in YES


def tri(v):
    """'yes' / 'no' / '' -- the deposit's ledgers spell booleans four ways."""
    s = str(v).strip().lower()
    if s in YES:
        return "yes"
    if s in NO:
        return "no"
    return ""


def as_int(v, default=None):
    try:
        return int(round(float(str(v).replace(",", "").strip())))
    except (TypeError, ValueError):
        return default


def as_float(v, default=None):
    try:
        return float(str(v).replace(",", "").strip())
    except (TypeError, ValueError):
        return default


def fnum(x, dp):
    """The generator's single formatter.  Reproduced because the PUBLISHED
    precision is the operative value: a re-measurement must round to the same
    string, not merely come close."""
    if x is None:
        return ""
    return "{0:.{1}f}".format(float(x), dp)


def pub(x, dp=FRAC_DP):
    """A measured fraction, taken to the precision at which it is published."""
    if x is None:
        return None
    return float(fnum(x, dp))


def pct(num, den):
    if den in (None, 0):
        return None
    return 100.0 * float(num) / float(den)


def close(a, b, atol):
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    try:
        return abs(float(a) - float(b)) <= atol
    except (TypeError, ValueError):
        return False


def commas(n):
    if n is None or n == "":
        return "--"
    try:
        return "{:,}".format(int(round(float(n))))
    except (TypeError, ValueError):
        return str(n)


def read_csv(path, required=True):
    if not os.path.exists(path):
        if required:
            raise IOError("missing required input: %s" % rel(path))
        return []
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh))


def csv_header(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh):
            return row
    return []


def walk_files(root, exts):
    """Every file under `root` with one of `exts`, skipping history."""
    out = []
    if not os.path.isdir(root):
        return out
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if not EXCLUDED_DIR_RE.match(d)]
        for fn in sorted(filenames):
            if fn.lower().endswith(exts):
                out.append(os.path.join(dirpath, fn))
    return sorted(out)


# ==========================================================================
# THE REPORT
# ==========================================================================

class Report(object):
    """Collects failures and notes.  Failures are fatal; notes are the things
    a reader must SEE but that are not defects -- above all the pre-audit
    ledger disagreements, which are the paper's finding, not its bug."""

    def __init__(self):
        self.checks = OrderedDict()
        self.current = None
        self.aborted = False

    def begin(self, cid, title):
        self.current = cid
        self.checks[cid] = OrderedDict([("id", cid), ("title", title),
                                        ("failures", []), ("notes", []),
                                        ("assertions", 0), ("skipped", None)])

    def bump(self, n=1):
        if self.current:
            self.checks[self.current]["assertions"] += n

    def fail(self, msg, detail=None):
        entry = OrderedDict([("message", msg)])
        if detail:
            entry["detail"] = detail if isinstance(detail, list) else [detail]
        self.checks[self.current]["failures"].append(entry)

    def note(self, msg, detail=None):
        entry = OrderedDict([("message", msg)])
        if detail:
            entry["detail"] = detail if isinstance(detail, list) else [detail]
        self.checks[self.current]["notes"].append(entry)

    def skip(self, reason):
        self.checks[self.current]["skipped"] = reason

    def eq(self, label, got, want, atol=0.0):
        """One assertion.  Numeric when BOTH sides parse as numbers, exact
        otherwise.  Every failure prints BOTH values -- a verifier that says
        only 'mismatch' is a verifier nobody uses."""
        self.bump()
        if got is None and want is None:
            return True
        gf, wf = as_float(got), as_float(want)
        if gf is not None and wf is not None:
            if abs(gf - wf) <= atol:
                return True
        elif got == want or str(got) == str(want):
            return True
        self.fail("%s: deposit says %r, this script measures %r"
                  % (label, got, want))
        return False

    @property
    def n_failures(self):
        return sum(len(c["failures"]) for c in self.checks.values())

    @property
    def n_notes(self):
        return sum(len(c["notes"]) for c in self.checks.values())


R = Report()


def run_check(cid, title, fn, *args, **kw):
    """Run one check.  An exception inside a check is itself a failure -- it is
    never allowed to abort the run silently or to be mistaken for a pass."""
    R.begin(cid, title)
    try:
        fn(*args, **kw)
    except Exception:                                      # noqa: BLE001
        R.fail("the check raised an exception (this is a FAILURE, not a skip)",
               traceback.format_exc().strip().splitlines())
    c = R.checks[cid]
    # FAIL outranks SKIP: a partially-scoped run that still found a defect
    # reports the defect, never the reduced scope.
    if c["failures"]:
        state = "FAIL"
    elif c["skipped"]:
        state = "SKIP"
    else:
        state = "ok"
    sys.stdout.write("  %-4s %-13s %-58s %d assertions\n"
                     % (state, cid, title[:58], c["assertions"]))
    sys.stdout.flush()
    for f in c["failures"]:
        sys.stdout.write("        ! %s\n" % f["message"])
        for d in f.get("detail", [])[:12]:
            sys.stdout.write("            %s\n" % d)
    return not c["failures"]


# ==========================================================================
# THE REFERENCE PARSER -- imported from the audit, never reimplemented
# ==========================================================================

def load_audit_common():
    if not os.path.exists(AUDIT_COMMON):
        raise IOError(
            "the audit's reference Newick parser is missing at %s. It is not "
            "optional: a naive one-pass parser gets ultrametricity wrong, and "
            "reimplementing it here would be the very defect this script "
            "exists to catch." % rel(AUDIT_COMMON))
    spec = importlib.util.spec_from_file_location("_atlas_audit_common",
                                                  AUDIT_COMMON)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    for fn in ("parse_newick", "leaf_depths"):
        if not hasattr(mod, fn):
            raise IOError("%s does not define %s()" % (rel(AUDIT_COMMON), fn))
    return mod


COMMON = None       # set in main()


def coefficient_of_variation(xs):
    if not xs:
        return None
    mean = sum(xs) / float(len(xs))
    if mean == 0:
        return None
    var = sum((x - mean) ** 2 for x in xs) / float(len(xs))
    return math.sqrt(var) / abs(mean)


# -- WHITESPACE NORMALISATION -------------------------------------------------
# The audit's reference parser does not skip whitespace at structural
# positions: after a comma it tests s[pos] == "(" directly, so a tree written
# as ", (" is read as an empty-named leaf and the rest of the file is
# mis-nested.  cetaceans.nwk is written that way and the reference parser reads
# its 85 leaves as 4; truth_tips.csv's 85 is correct and the PARSER is the
# defect.  The parser is still reused unmodified -- it is handed normalised
# text.  Normalisation is quote-aware and touches only whitespace adjacent to a
# structural character, so quoted labels such as 'Acanthella acuta' survive.
_NEWICK_STRUCTURAL = "(),;:"


def normalize_newick(text):
    """Drop whitespace at structural positions only; never inside a label."""
    out = []
    i, n, quote = 0, len(text), None
    while i < n:
        ch = text[i]
        if quote is not None:
            out.append(ch)
            if ch == quote:
                quote = None
            i += 1
            continue
        if ch in "'\"":
            quote = ch
            out.append(ch)
            i += 1
            continue
        if ch.isspace():
            j = i
            while j < n and text[j].isspace():
                j += 1
            prev = out[-1] if out else None
            nxt = text[j] if j < n else None
            if prev is None or nxt is None \
                    or prev in _NEWICK_STRUCTURAL or nxt in _NEWICK_STRUCTURAL:
                i = j
                continue
            out.append(" ")
            i = j
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def _has_structural_whitespace(text):
    """True iff whitespace sits INSIDE the tree between structural characters.

    A trailing newline after the final ';' is cosmetic and does not count --
    the reference parser strips it.  What counts is whitespace the reference
    parser will mistake for a label, e.g. the ', (' and ' : ' of cetaceans.nwk.
    """
    body = text.strip()
    quote = None
    for i, ch in enumerate(body):
        if quote is not None:
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            continue
        if ch.isspace():
            j = i
            while j < len(body) and body[j].isspace():
                j += 1
            prev = body[i - 1] if i else None
            nxt = body[j] if j < len(body) else None
            if (prev and prev in _NEWICK_STRUCTURAL) \
                    or (nxt and nxt in _NEWICK_STRUCTURAL):
                return True
    return False


def leaf_depths_text(text):
    """COMMON.leaf_depths' exact convention, on normalised text.

    leaf_depths() takes a PATH and so cannot be handed normalised text.  This
    reproduces its traversal verbatim -- every leaf's label enters the distinct
    count, including unnamed ones, and a NaN branch length contributes zero
    depth and is counted -- over normalize_newick(text).

    It is NOT a second parser: it calls COMMON.parse_newick.  And it carries a
    PROOF OBLIGATION discharged per tree by measure_trees(): where normalisation
    is a no-op, this function must return exactly what COMMON.leaf_depths()
    returns for the same file, or the run fails.
    """
    root = COMMON.parse_newick(normalize_newick(text))
    depths, labels, nan_bl = [], [], 0
    stack = [(root, 0.0)]
    while stack:
        (kids, name, bl), acc = stack.pop()
        if bl != bl:                        # NaN branch length
            nan_bl += 1
            bl = 0.0
        d = acc + bl
        if kids:
            for k in kids:
                stack.append((k, d))
        else:
            depths.append(d)
            labels.append(name)
    return depths, len(depths), len(set(labels)), nan_bl


def tree_facts(text):
    """(named leaf labels, root-to-tip depths, n_leaves, n_nan_branch_lengths).

    The audit's NaN convention: a NaN branch length contributes zero depth and
    is counted.  Unnamed leaves are counted as leaves but contribute no label,
    which is the convention the intersections use; truth_dated.csv's
    `distinct_labels` uses leaf_depths' own convention and is compared against
    that function directly, never against this one."""
    root = COMMON.parse_newick(normalize_newick(text))
    labels, depths, nan_bl = [], [], 0
    stack = [(root, 0.0)]
    while stack:
        (kids, name, bl), acc = stack.pop()
        if bl != bl:
            nan_bl += 1
            bl = 0.0
        d = acc + bl
        if kids:
            for k in kids:
                stack.append((k, d))
        else:
            depths.append(d)
            if name:
                labels.append(name)
    return labels, depths, len(depths), nan_bl


def read_first_tree(path, max_bytes=None):
    """The first tree of a Newick or NEXUS file, with a NEXUS translate table
    applied and [&...] comments removed.  Read line by line, so a large
    multi-tree posterior set costs only the bytes up to the end of its first
    tree; `max_bytes` bounds how far the search for that tree may go."""
    try:
        fh = open(path, "r", encoding="utf-8", errors="replace")
    except OSError:
        return None
    translate, in_tr, buf, seen, nexus = {}, False, None, 0, None
    with fh:
        for line in fh:
            seen += len(line)
            s = re.sub(r"\[[^\]]*\]", "", line).strip()
            if nexus is None and s:
                nexus = s.upper().startswith("#NEXUS")
            if buf is not None:
                buf += s
                if ";" in s:
                    break
                continue
            if nexus:
                low = s.lower()
                if low.startswith("translate"):
                    in_tr = True
                    s = s[len("translate"):].strip()
                    if not s:
                        continue
                if in_tr:
                    for item in s.rstrip(";").split(","):
                        kv = item.split(None, 1)
                        k, v = (kv[0], kv[1]) if len(kv) == 2 else ("", "")
                        if k:
                            translate[k] = v.strip().strip("'")
                    if s.endswith(";"):
                        in_tr = False
                    continue
                m = re.match(r"^u?tree\b[^=]*=\s*(\(.*)$", s, re.I)
                if m:
                    buf = m.group(1)
                    if ";" in buf:
                        break
            elif "(" in s:
                buf = s[s.index("("):]
                if ";" in buf:
                    break
            if max_bytes is not None and seen > max_bytes:
                return None
    if buf is None or ";" not in buf:
        return None
    buf = buf[:buf.index(";") + 1]
    if translate:
        buf = re.sub(r"([(,])\s*([^(),:;\s]+)",
                     lambda m: m.group(1) + translate.get(m.group(2),
                                                          m.group(2)), buf)
    return buf

# ==========================================================================
# CONTEXT -- everything measured, in one object
# ==========================================================================

class Ctx(object):
    def __init__(self, quick=False):
        self.quick = quick
        self.dated_rows = read_csv(TRUTH["dated"])
        self.tips_rows = read_csv(TRUTH["tips"])
        self.species_rows = read_csv(TRUTH["species"])
        self.prov_rows = read_csv(TRUTH["provenance"])

        self.dated = OrderedDict((r["filename"], r) for r in self.dated_rows)
        self.tips = OrderedDict((r["filename"], r) for r in self.tips_rows)
        self.species = OrderedDict((r["partition"], r)
                                   for r in self.species_rows)
        self.canon_prov = OrderedDict(
            (r["partition_name"], r) for r in self.prov_rows
            if truthy(r.get("is_canonical")))
        self.canon_file = OrderedDict(
            (p, (r.get("filename") or "").strip())
            for p, r in self.canon_prov.items())

        self.eukaryotic = [p for p, r in self.species.items()
                           if not truthy(r.get("is_prokaryotic"))]
        euk = set(self.eukaryotic)
        self.canon_euk = [p for p in self.canon_prov if p in euk]
        self.canon_prok = [p for p in self.canon_prov if p not in euk]

        # ---- name-resolution layer -------------------------------------
        self.by_id = {}
        with open(DICTIONARY, newline="", encoding="utf-8",
                  errors="replace") as fh:
            for row in csv.reader(fh):
                if len(row) >= 2 and row[0] != "id":
                    self.by_id[row[0]] = row[1]
        self.by_group = defaultdict(dict)
        for row in read_csv(FULL_MAPPING, required=False):
            g, o, s = (row.get("group"), row.get("original_label"),
                       row.get("standardized_name"))
            if g and o and s:
                self.by_group[g][o] = s

        # ---- measurements, filled by measure_trees() --------------------
        self.measured = {}        # filename -> dict of measurements
        self.labels = {}          # partition -> set(standardised name)
        self.idmap = {}           # partition -> {standardised name -> id}
        self.source_cache = {}
        self.tip_evidence = read_csv(MOL_TIP_EVIDENCE, required=False)
        self.time_scale = {(r.get("filename") or "").strip():
                           (r.get("time_scale") or "").strip()
                           for r in read_csv(TIME_SCALE_CHECKS, required=False)}
        self.col_rows = defaultdict(dict)
        for r in read_csv(COL_RESOLUTION, required=False):
            self.col_rows[r["partition"]][r["standardized_id"]] = r
        self.c0_names = {}        # partition -> {counting unit -> display}
        self.unit_of = {}         # partition -> {standardized id -> unit}
        self.c0_unresolved = {}
        self.c1_sets = {}         # partition -> set of C1 units, or None

        self.crown = OrderedDict(
            (r["partition"], r)
            for r in read_csv(CROWN_INTERVALS, required=False))

        # ---- deposit outputs under test ---------------------------------
        self.basis_rows = read_csv(BASIS)
        self.basis = OrderedDict((r["partition"], r) for r in self.basis_rows)
        with open(HEADLINE, encoding="utf-8") as fh:
            self.head = json.load(fh)

    # -- accessors on ground truth ------------------------------------
    def group(self, part):
        return (self.species[part].get("group") or "").strip()

    def described(self, part):
        return as_int(self.species[part].get("described"), 0)

    def labels_unique(self, part):
        return as_int(self.tips[self.canon_file[part]]["tips_adopted"], 0)

    def labels_raw(self, part):
        return as_int(self.tips[self.canon_file[part]]["raw_leaves"], 0)

    def is_unrepresented(self, part):
        return truthy(self.species[part].get("is_unrepresented"))

    def is_prokaryotic(self, part):
        return truthy(self.species[part].get("is_prokaryotic"))

    def archival_class(self, part):
        return (self.canon_prov[part].get("archival_class") or "").strip()

    def interval_recoverable(self, part):
        """C1b flag: a counted row in Data/crown_intervals.csv (S02)."""
        r = self.crown.get(part)
        return bool(r) and r.get("counts_for_c1b") == "yes"

    def posterior_n(self, part):
        r = self.crown.get(part)
        if not r or r.get("counts_for_c1b") != "yes" \
                or as_int(r.get("n_trees"), 0) <= 1:
            return ""
        return str(as_int(r["n_trees"]))

    # -- MEASURED chronogram verdict, never a ledger flag --------------
    def cv(self, part):
        return self.measured[self.canon_file[part]]["cv"]

    def is_geometric_chronogram(self, part):
        cv = self.cv(part)
        return cv is not None and cv <= CHRONOGRAM_TOLERANCE

    def is_chronogram(self, part):
        """D14: ultrametric AND in absolute time."""
        return (self.is_geometric_chronogram(part)
                and self.time_scale.get(self.canon_file[part]) == "absolute")

    # -- name resolution ------------------------------------------------
    def resolve_ids(self, labels):
        return set(self.by_id.get(lab, lab) for lab in labels)

    def resolve_originals(self, labels, group):
        m = self.by_group.get(group, {})
        return set(m.get(lab, lab) for lab in labels)

    def source_sets(self, path, group):
        """(resolved labels, verbatim labels) of a deposited source tree."""
        key = (path, group)
        if key in self.source_cache:
            return self.source_cache[key]
        empty = (frozenset(), frozenset())
        if path.endswith(TAXA_LIST_EXT):
            with open(path, encoding="utf-8") as fh:
                labels = [ln.strip() for ln in fh
                          if ln.strip() and not ln.startswith("#")]
            self.source_cache[key] = (
                frozenset(self.resolve_originals(labels, group)),
                frozenset(labels))
            return self.source_cache[key]
        text = read_first_tree(path, SOURCE_SCAN_MAX_BYTES)
        if not text:
            self.source_cache[key] = empty
            return empty
        try:
            labels, _d, _n, _nb = tree_facts(text)
        except Exception:                                  # noqa: BLE001
            self.source_cache[key] = empty
            return empty
        self.source_cache[key] = (
            frozenset(self.resolve_originals(labels, group)),
            frozenset(labels))
        return self.source_cache[key]


def classify_assembly(otl, inference):
    """The NOMINAL topology_assembly annotation.  Reimplemented here (it is
    four lines) so the classifier the deposit prints in Methods is verified,
    not inherited."""
    if (otl or "").strip() == ASSEMBLY_OTL_TRIGGER:
        return ASSEMBLY_AUTOMATED
    if (inference or "").strip() in ASSEMBLY_INFERENCE_TRIGGERS:
        return ASSEMBLY_COMPILED
    return ASSEMBLY_PRIMARY


# ==========================================================================
# CHECK: INPUTS
# ==========================================================================

REQUIRED = [
    ("Data/_rebuild/truth_dated.csv", TRUTH["dated"], "measured ground truth"),
    ("Data/_rebuild/truth_tips.csv", TRUTH["tips"], "measured ground truth"),
    ("Data/_rebuild/truth_species.csv", TRUTH["species"], "measured ground truth"),
    ("Data/_rebuild/truth_provenance.csv", TRUTH["provenance"],
     "measured ground truth"),
    ("Paper/audit-2026-08-20/code/common.py", AUDIT_COMMON,
     "the reference Newick parser"),
    ("Data/standardized/trees", STD_TREES, "the shipped trees"),
    ("Data/standardized/dictionary.csv", DICTIONARY, "id -> name"),
    ("Data/standardized/full_mapping.csv", FULL_MAPPING, "original -> name"),
    ("Data/tree_provenance.csv", TREE_PROVENANCE,
     "pre-audit ledger; annotation fields only"),
    ("Data/coverage_basis.csv", BASIS, "run scripts/rebuild_atlas.py first"),
    ("Data/molecular_tips.csv", MOLTIPS, "run scripts/rebuild_atlas.py first"),
    ("Paper/inputs/headline_numbers.json", HEADLINE,
     "run scripts/rebuild_atlas.py first"),
    ("Paper/inputs/coverage_per_partition.csv", COVERAGE_PP,
     "run scripts/rebuild_atlas.py first"),
    ("Paper/inputs/fig1_coverage_bars.csv", FIG1_BARS,
     "run scripts/rebuild_atlas.py first"),
    ("Paper/inputs/fig2_coverage_scatter.csv", FIG2_SCATTER,
     "run scripts/rebuild_atlas.py first"),
    ("Paper/inputs/fig3_c0c1_bars.csv", FIG3_BARS,
     "run scripts/rebuild_atlas.py first"),
    ("Paper/inputs/fig3_dated_undated_counts.csv", FIG3_COUNTS,
     "run scripts/rebuild_atlas.py first"),
    ("Data/site/partition_summary.csv", PARTITION_SUMMARY,
     "run scripts/rebuild_atlas.py first"),
]


def check_inputs():
    missing = []
    for label, path, hint in REQUIRED:
        R.bump()
        if not os.path.exists(path):
            missing.append(label)
            R.fail("missing %s (%s)" % (label, hint))
    if not os.path.isdir(PRE_AUDIT):
        R.note("the pre-audit snapshot %s is absent; the PREAUDIT check falls "
               "back to the live ledgers, which the rebuild may have already "
               "overwritten" % rel(PRE_AUDIT))
    return not missing


# ==========================================================================
# CHECK: TREES -- re-measure every shipped .nwk
# ==========================================================================

def measure_trees(ctx):
    """Re-measure the shipped trees with the audit's own parser and assert
    truth_dated.csv / truth_tips.csv against the files.  This is the check that
    makes every downstream number rest on geometry."""
    on_disk = set(fn for fn in os.listdir(STD_TREES) if fn.endswith(".nwk"))
    listed = set(ctx.dated)
    R.bump()
    if on_disk != listed:
        extra = sorted(on_disk - listed)
        gone = sorted(listed - on_disk)
        R.fail("Data/standardized/trees and truth_dated.csv disagree on the "
               "set of shipped trees",
               ["on disk only: %s" % (extra[:10] or "-"),
                "in ledger only: %s" % (gone[:10] or "-")])

    R.bump()
    if set(ctx.dated) != set(ctx.tips):
        R.fail("truth_dated.csv and truth_tips.csv key different file sets "
               "(%d vs %d rows)" % (len(ctx.dated), len(ctx.tips)))

    ctx.canon_sources = {}
    ctx.source_leaves = {}
    ctx.n_pruned = defaultdict(int)
    for r in read_csv(P("Data", "canonical_sources.csv"), required=False):
        ctx.canon_sources[(r.get("filename") or "").strip()] = \
            (r.get("source") or "").strip()
    for r in read_csv(P("Data", "pruned_tips.csv"), required=False):
        ctx.n_pruned[(r.get("filename") or "").strip()] += 1

    canon_files = set(ctx.canon_file.values())
    targets = sorted(canon_files) if ctx.quick else sorted(listed & on_disk)
    if ctx.quick:
        R.skip("--quick: only the %d partition-canonical trees were "
               "re-measured; %d shipped trees were NOT verified against their "
               "files. NOT a valid acceptance run."
               % (len(canon_files), len(listed) - len(canon_files)))
        R.note("--quick was used. A green run under --quick is not an "
               "acceptance result.")

    for fname in targets:
        path = os.path.join(STD_TREES, fname)
        d = ctx.dated.get(fname) or {}
        t = ctx.tips.get(fname) or {}
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
            depths, n_leaves, n_distinct, n_nan = leaf_depths_text(text)
            # PROOF OBLIGATION, discharged on EVERY tree.  The shim must agree
            # with COMMON.leaf_depths called on the file itself, and where it
            # does not, the disagreement must be attributable to structural
            # whitespace -- otherwise the shim has drifted from the parser it
            # wraps and is a second, divergent parser, which is the defect this
            # script exists to catch.
            ref = COMMON.leaf_depths(path)
            R.bump()
            if (ref[1], ref[2], ref[3]) == (n_leaves, n_distinct, n_nan):
                pass
            elif _has_structural_whitespace(text):
                R.note("%s: carries whitespace at structural positions, which "
                       "the audit's reference parser does not skip (after ',' "
                       "or before '('). It reads %d leaves / %d distinct where "
                       "the file has %d / %d. The file is normalised before "
                       "parsing; truth_tips.csv is correct and the PARSER was "
                       "the defect."
                       % (fname, ref[1], ref[2], n_leaves, n_distinct))
            else:
                R.fail("%s: leaf_depths_text disagrees with the reference "
                       "COMMON.leaf_depths (%s vs %s) on a file that carries "
                       "no structural whitespace -- the normalisation shim has "
                       "drifted from the parser it wraps"
                       % (fname, (n_leaves, n_distinct, n_nan),
                          (ref[1], ref[2], ref[3])))
        except RecursionError:
            R.bump()
            R.fail("%s: the reference parser exhausted the stack; raise "
                   "sys.setrecursionlimit or report the tree" % fname)
            continue
        except Exception as exc:                           # noqa: BLE001
            R.bump()
            R.fail("%s: could not be parsed (%s)" % (fname, exc))
            continue

        cv = coefficient_of_variation(depths)
        chrono = cv is not None and cv <= CHRONOGRAM_TOLERANCE
        # A second traversal gives the label set the intersections need and
        # the NaN count on the tree_facts convention.
        named, depths2, _n2, nan2 = tree_facts(text)

        ctx.measured[fname] = {
            "cv": cv, "chronogram": chrono, "leaves": n_leaves,
            "distinct": n_distinct, "nan": n_nan, "nan_named": nan2,
            "named_labels": named,
            "root_height": max(depths2) if depths2 else 0.0,
        }

        R.eq("%s leaves" % fname, as_int(d.get("leaves")), n_leaves)
        R.eq("%s distinct_labels" % fname, as_int(d.get("distinct_labels")),
             n_distinct)
        R.eq("%s nan_branch_lengths" % fname,
             as_int(d.get("nan_branch_lengths")), n_nan)
        R.eq("%s roottotip_cv" % fname, as_float(d.get("roottotip_cv")), cv,
             atol=CV_RECHECK_TOL)
        R.eq("%s is_chronogram" % fname, tri(d.get("is_chronogram")),
             "yes" if chrono else "no")

        # truth_tips.csv, on its own declared convention
        R.eq("%s raw_leaves" % fname, as_int(t.get("raw_leaves")), n_leaves)
        R.eq("%s tips distinct_labels" % fname,
             as_int(t.get("distinct_labels")), n_distinct)
        # truth_tips.csv's `convention` is "<rule>" optionally followed by
        # "; published=<reason>".  The RULE is the settled tip convention and
        # must be one of the two admissible ones.  The optional REASON explains
        # why the PRE-AUDIT published count differs, and is not a second rule.
        conv_raw = (t.get("convention") or "").strip()
        parts = [s.strip() for s in conv_raw.split(";") if s.strip()]
        conv = parts[0] if parts else ""
        reasons = [s for s in parts[1:] if s.startswith("published=")]
        strays = [s for s in parts[1:] if not s.startswith("published=")]
        adopted = as_int(t.get("tips_adopted"))
        R.bump()
        if conv == "distinct_labels":
            if adopted != n_distinct:
                R.fail("%s: convention='distinct_labels' but tips_adopted=%s "
                       "against %d distinct labels in the file"
                       % (fname, adopted, n_distinct))
        elif conv == "raw_leaves":
            if adopted != n_leaves:
                R.fail("%s: convention='raw_leaves' but tips_adopted=%s "
                       "against %d leaves in the file"
                       % (fname, adopted, n_leaves))
        else:
            R.fail("%s: unknown tip convention %r -- the tip convention must "
                   "be one settled, named rule" % (fname, conv))
        R.bump()
        if strays:
            R.fail("%s: convention %r carries a clause that is neither the "
                   "rule nor a published= reason: %s"
                   % (fname, conv_raw, strays))
        # source_ntips: the leaf count of the source tree this canonical file
        # is built from (Data/canonical_sources.csv), measured here from that
        # file.  It is NOT published_ntips, which keeps the count the published
        # table carried; the two columns answer different questions.
        src_n = as_int(t.get("source_ntips"))
        src_path = ctx.canon_sources.get(fname)
        R.bump()
        if src_path is None:
            if src_n is not None:
                R.fail("%s: truth_tips.csv carries source_ntips=%s but "
                       "Data/canonical_sources.csv names no source for it"
                       % (fname, src_n))
        elif not os.path.exists(P(src_path)):
            R.fail("%s: its source %s (Data/canonical_sources.csv) is missing"
                   % (fname, src_path))
        else:
            with open(P(src_path), encoding="utf-8", errors="replace") as fh:
                _l, _d, src_leaves, _nan = tree_facts(fh.read())
            ctx.source_leaves[fname] = src_leaves
            R.eq("%s source_ntips vs the leaves of %s" % (fname, src_path),
                 src_n, src_leaves)
            n_d7 = ctx.n_pruned.get(fname, 0)
            R.bump()
            if adopted is not None and src_leaves - adopted < n_d7:
                R.fail("%s: the source has %d leaves and %d are pruned by D7 "
                       "(Data/pruned_tips.csv), yet %s tips are shipped"
                       % (fname, src_leaves, n_d7, adopted))
        pub_n = as_int(t.get("published_ntips"))
        delta = as_int(t.get("delta"))
        if pub_n is not None and delta is not None and adopted is not None:
            R.bump()
            if abs(delta) != abs(adopted - pub_n):
                R.fail("%s: delta=%s does not reconcile tips_adopted=%s with "
                       "published_ntips=%s" % (fname, delta, adopted, pub_n))
            # A discrepancy with the pre-audit published count must ALWAYS be
            # explained, and an explanation must never be attached where there
            # is nothing to explain.  This is what lets CROSSLEDGER downgrade a
            # stale ledger tip count to a note: the ground truth has to have
            # anticipated it and named a reason.
            R.bump()
            if delta != 0 and not reasons:
                R.fail("%s: tips_adopted=%s differs from published_ntips=%s "
                       "(delta %s) but convention %r gives no published= "
                       "reason. Every discrepancy with the published count "
                       "must be explained in the ground truth."
                       % (fname, adopted, pub_n, delta, conv_raw))
            elif delta == 0 and reasons:
                R.fail("%s: convention %r explains a discrepancy with the "
                       "published count, but tips_adopted and published_ntips "
                       "agree (%s). An unexplained explanation is a stale note."
                       % (fname, conv_raw, adopted))

    # Label sets for the canonical partitions, resolved through the deposit's
    # own id mapping.  These are the sets every numerator is formed from.
    for part, fname in ctx.canon_file.items():
        m = ctx.measured.get(fname)
        if m is None:
            R.bump()
            R.fail("canonical tree %s of partition %r was not measured"
                   % (fname, part))
            continue
        raw = m["named_labels"]
        ctx.labels[part] = ctx.resolve_ids(raw)
        first = {}
        for lab in raw:
            nm = ctx.by_id.get(lab, lab)
            if nm not in first or lab < first[nm]:
                first[nm] = lab
        ctx.idmap[part] = first
        R.eq("%s labels_unique vs the shipped file" % part,
             ctx.labels_unique(part), len(ctx.labels[part]))
        # C0 counting units: GTDB species clusters for a prokaryotic partition
        # (D12), CoL accepted names for a eukaryotic one (S05).
        if ctx.is_prokaryotic(part):
            ctx.c0_names[part] = {nm: nm for nm in ctx.labels[part]}
            ctx.c0_unresolved[part] = 0
            ctx.unit_of[part] = {sid: nm for nm, sid in first.items()}
            continue
        rows = ctx.col_rows.get(part, {})
        tips = set(first.values())
        R.bump()
        if set(rows) != tips:
            R.fail("%s: col_resolution.csv rows (%d) do not equal the shipped "
                   "tips (%d)" % (part, len(rows), len(tips)))
        acc, unres, unit_of = {}, 0, {}
        for sid in tips:
            r = rows.get(sid) or {}
            if r.get("match_type") in COL_RESOLUTION_RESOLVED and r.get("col_id"):
                acc[r["col_id"]] = r.get("col_name") or r["col_id"]
                unit_of[sid] = r["col_id"]
            else:
                unres += 1
        ctx.c0_names[part] = acc
        ctx.unit_of[part] = unit_of
        ctx.c0_unresolved[part] = unres


# ==========================================================================
# CHECK: BRANCHLEN -- every shipped tree reproduces its declared reference
# ==========================================================================
#
# The shipped trees once carried branch lengths rewritten by a relabelling step
# that read every digit run followed by ',', ')' or ':' as an old tip id, dropped
# its leading zeros and mapped it to the new id (0.007143 -> 0.7146;
# 27.05691095 -> 27.5691095); see scripts/restore_branch_lengths.py.  The CV filter cannot catch that: a
# corrupted ultrametric tree can still pass a CV <= 0.05 test.  So every shipped
# tree must declare, in Data/standardized/tree_references.csv, the file its
# branch lengths come from and the rule that ties the two:
#
#   tokens_equal_after_tip_relabel  the reference is the pre-remap standardized
#       copy; tip ids map old id -> standardised name -> current id, and EVERY
#       other token (branch length, internal label, comment) is identical.
#   edge_length_multiset_equal      the reference is a deposited source tree;
#       the multisets of branch lengths agree to BL_REL_TOL.  Checked when the
#       source is present (Data/raw is not distributed); a note otherwise.

_BL_TOKEN_RE = re.compile(r"[(),;]|:[^,();\[]*|\[[^\]]*\]|'(?:[^']|'')*'"
                          r"|[^(),:;\[']+")


def _bl_tokens(text):
    out = []
    for m in _BL_TOKEN_RE.finditer(text):
        t = m.group(0)
        if len(t) == 1 and t in "(),;":
            out.append(("p", t))
        elif t.startswith(":"):
            out.append(("b", t))
        elif t.startswith("["):
            out.append(("c", t))
        elif t.strip():
            out.append(("l", t.strip()))
    return out


def _edge_lengths(text):
    text = re.sub(r"\[[^\]]*\]", "", text)     # [&...] annotations
    root = COMMON.parse_newick(normalize_newick(text))
    out, stack = [], [(root, True)]
    while stack:
        (kids, _name, bl), is_root = stack.pop()
        if not is_root:
            out.append(bl)
        for k in kids:
            stack.append((k, False))
    return sorted(x for x in out if x == x)


_ZERO_AFTER_POINT = re.compile(r"(:\d+)\.0+(\d)")


def _relabel_diff(a, b, old2new):
    """(number of differing tokens, first difference) between a reference
    token list `a` and a shipped token list `b`, tips compared through the
    old -> current id map and every other token compared verbatim."""
    if len(a) != len(b):
        return 1, ("token count", len(a), len(b))
    bad, first, prev = 0, None, None
    for (ka, va), (kb, vb) in zip(a, b):
        if ka != kb:
            return bad + 1, first or ("structure", va, vb)
        if ka == "l" and prev in ("(", ",", None):
            if old2new.get(va.strip("'")) != vb.strip("'"):
                bad += 1
                first = first or ("tip", va, vb)
        elif va != vb:
            bad += 1
            first = first or (ka, va, vb)
        prev = va if ka == "p" else ka
    return bad, first


# --------------------------------------------------------------------------
# D4 species key and canonical re-derivation, written independently of
# scripts/canonical_trees.py so that the check is a re-derivation.
# --------------------------------------------------------------------------

D4_OPEN = {"sp", "spp", "sp.", "spp.", "cf", "cf.", "aff", "aff.", "af", "nr",
           "nr.", "indet", "indet.", "x", "hybrid", "nov", "n", "undescribed",
           "unidentified", "unknown", "gen", "trans", "genome"}
D4_QUALIFIER = {"x", "hybrid", "cf", "cf.", "aff", "aff.", "af"}
D4_SUPRA = ("idae", "aceae", "inae", "oidea", "ales", "phyta", "mycota",
            "mycetes", "opsida", "ineae")


# A token with a suprageneric suffix is still a genus when the Catalogue of
# Life lists a genus of that name (audit A10: Xerophyta, Elapsoidea,
# Pimephales, Smittoidea ...).  The list is committed; check_names re-derives
# it from the release index when the index is present.
D4_SUPRA_GENERA_FILE = P("Data", "standardized",
                         "col269_suprageneric_genera.txt")
COL_INDEX = P("Data", "raw", "col", "col269_species_index.tsv")


def _read_supra_genera():
    out = set()
    if os.path.exists(D4_SUPRA_GENERA_FILE):
        with open(D4_SUPRA_GENERA_FILE, encoding="utf-8") as fh:
            out = {ln.strip() for ln in fh
                   if ln.strip() and not ln.startswith("#")}
    return frozenset(out)


D4_SUPRA_GENERA = _read_supra_genera()


def d4_is_genus(g):
    if not re.match(r"^[A-Z][a-z]+$", g):
        return False
    return g in D4_SUPRA_GENERA or not g.lower().endswith(D4_SUPRA)


def d4_species_key(name, sid):
    toks = [t for t in re.split(r"[_.\s]+", name or "") if t]
    if {t.lower() for t in toks} & D4_QUALIFIER:
        return None
    for g, e in zip(toks, toks[1:]):
        genus = d4_is_genus(g)
        if genus and re.match(r"^[a-z][a-z-]+$", e):
            if e.lower() in D4_OPEN or re.match(r"^sp\d*$", e):
                return None
            return g + "_" + e
        if genus and e.lower() in D4_OPEN:
            return None
    return None


def d4_keys(labels):
    out = set()
    for nm in labels:
        k = d4_species_key(nm, nm)
        out.add(k if k else "label:" + nm)
    return out


def _leaf_order(text):
    """Leaf labels, left to right."""
    out = []
    stack = [COMMON.parse_newick(normalize_newick(text))]
    while stack:
        kids, name, _bl = stack.pop()
        if kids:
            stack.extend(reversed(kids))
        else:
            out.append(name.strip("'"))
    return out


def _rederive_canonical(src_text, pruned_ids, by_id):
    """(retained leaf order, sorted edge lengths) of the source after removing
    the pruned ids and, left to right, every later tip of a species key
    already seen.  An edge of the result is the sum of the source edges from a
    retained node up to its nearest retained ancestor; the result's root is the
    MRCA of the retained tips and has no edge."""
    root = COMMON.parse_newick(normalize_newick(src_text))
    seen = set()
    nodes, parent, order = [root], {0: None}, []
    stack = [0]
    while stack:                      # pre-order, children left to right
        i = stack.pop()
        order.append(i)
        kids = nodes[i][0]
        base = len(nodes)
        nodes.extend(kids)
        for j in range(len(kids)):
            parent[base + j] = i
        stack.extend(range(base + len(kids) - 1, base - 1, -1))
    keep_leaf = {}
    for i in order:
        kids, name, _bl = nodes[i]
        if kids:
            continue
        sid = name.strip("'")
        if sid in pruned_ids:
            keep_leaf[i] = False
            continue
        key = d4_species_key(by_id.get(sid, ""), sid) or ("id:" + sid)
        keep_leaf[i] = key not in seen
        seen.add(key)
    count = defaultdict(int)          # retained tips below
    nkids = defaultdict(int)          # children with a retained tip
    for i in reversed(order):
        if not nodes[i][0]:
            count[i] = 1 if keep_leaf[i] else 0
        p = parent[i]
        if p is not None and count[i]:
            count[p] += count[i]
            nkids[p] += 1
    retained = {i for i in order
                if (not nodes[i][0] and keep_leaf[i]) or nkids[i] >= 2}
    total = count[0]
    edges, acc = [], {}
    for i in order:                   # sum of edges up to nearest retained
        p = parent[i]
        bl = nodes[i][2]
        bl = 0.0 if bl is None else bl
        acc[i] = bl if p is None or p in retained else acc[p] + bl
        if i in retained and p is not None and count[i] < total:
            edges.append(acc[i])
        elif i in retained and p is not None and count[i] == total:
            pass                      # the new root: its stem is dropped
    leaves_out = [nodes[i][1].strip("'") for i in order
                  if not nodes[i][0] and keep_leaf[i]]
    return leaves_out, sorted(x for x in edges if x == x)


def _multiset_off(ea, eb):
    return [(x, y) for x, y in zip(ea, eb)
            if abs(x - y) > BL_REL_TOL * max(1.0, abs(x))]


def check_branchlen(ctx):
    rows = read_csv(TREE_REFERENCES, required=False)
    R.bump()
    if not rows:
        R.fail("%s is missing or empty: no shipped tree declares the file its "
               "branch lengths come from" % rel(TREE_REFERENCES))
        return
    on_disk = sorted(fn for fn in os.listdir(STD_TREES) if fn.endswith(".nwk"))
    declared = defaultdict(list)
    for r in rows:
        declared[(r.get("filename") or "").strip()].append(r)
    for fn in on_disk:
        R.bump()
        if fn not in declared:
            R.fail("%s declares no reference in %s"
                   % (fn, rel(TREE_REFERENCES)))
    for fn in sorted(declared):
        R.bump()
        if "/" in fn:
            # an archived file that a shipped tree is built from: it must exist
            if not os.path.exists(P(*fn.split("/"))):
                R.fail("%s declares a reference for %s, which does not exist"
                       % (rel(TREE_REFERENCES), fn))
        elif fn not in on_disk:
            R.fail("%s declares a reference for %s, which is not shipped"
                   % (rel(TREE_REFERENCES), fn))

    # old id -> current id through the standardised name
    need_relabel = any((r.get("rule") or "").strip()
                       == "tokens_equal_after_tip_relabel" for r in rows)
    old2new = {}
    if need_relabel:
        name2new = {}
        for k, v in ctx.by_id.items():
            name2new.setdefault(v.replace(" ", "_"), k)
        for r in read_csv(PRE_REMAP_DICT):
            old2new[r["id"]] = name2new.get(
                (r["standardized_name"] or "").replace(" ", "_"))

    self_test = None
    self_test_ms = None
    pruned = defaultdict(set)
    for r in read_csv(PRUNED_TIPS, required=False):
        pruned[(r.get("filename") or "").strip()].add(
            (r.get("standardized_id") or "").strip())
    targets = list(on_disk) + sorted(f for f in declared if "/" in f)
    for fn in targets:
        tpath = P(*fn.split("/")) if "/" in fn else os.path.join(STD_TREES, fn)
        if not os.path.exists(tpath):
            continue
        with open(tpath, encoding="utf-8", errors="replace") as fh:
            shipped = fh.read()
        for r in declared.get(fn, []):
            rule = (r.get("rule") or "").strip()
            ref = P(*(r.get("reference") or "").strip().split("/"))
            R.bump()
            if rule == "tokens_equal_after_tip_relabel":
                if not os.path.exists(ref):
                    R.fail("%s: declared reference %s is missing"
                           % (fn, rel(ref)))
                    continue
                with open(ref, encoding="utf-8", errors="replace") as fh:
                    a = _bl_tokens(fh.read())
                b = _bl_tokens(shipped)
                bad, first = _relabel_diff(a, b, old2new)
                if self_test is None and not bad and any(
                        k == "b" and _ZERO_AFTER_POINT.search(v) for k, v in a):
                    self_test = (fn, a, shipped)
                if bad:
                    R.fail("%s: %d token(s) differ from the reference %s after "
                           "tip relabelling; first: %s %r -> %r"
                           % (fn, bad, rel(ref), first[0], first[1], first[2]))
            elif rule == "edge_length_multiset_equal":
                if not os.path.exists(ref):
                    R.note("%s: declared source %s is not present (Data/raw is "
                           "not distributed); its branch lengths were NOT "
                           "compared on this run" % (fn, rel(ref)))
                    continue
                src = read_first_tree(ref)
                scale = as_float((r.get("scale") or "").strip()) or 1.0
                ea = sorted(x * scale for x in _edge_lengths(src or ""))
                eb = _edge_lengths(shipped)
                if len(ea) != len(eb):
                    R.fail("%s: %d branch lengths against %d in its source %s"
                           % (fn, len(eb), len(ea), rel(ref)))
                    continue
                off = _multiset_off(ea, eb)
                if not off and self_test_ms is None and eb:
                    self_test_ms = (fn, ea, eb)
                if off:
                    R.fail("%s: %d of %d branch lengths differ from its source "
                           "%s; first %r -> %r"
                           % (fn, len(off), len(ea), rel(ref), off[0][0],
                              off[0][1]))
            elif rule == "built_from_source":
                if not os.path.exists(ref):
                    R.fail("%s: declared source %s is missing" % (fn, rel(ref)))
                    continue
                with open(ref, encoding="utf-8", errors="replace") as fh:
                    src = fh.read()
                want_leaves, want_edges = _rederive_canonical(
                    src, pruned.get(fn, set()), ctx.by_id)
                got_leaves = _leaf_order(shipped)
                if got_leaves != want_leaves:
                    R.fail("%s: tips differ from its source %s after D7 pruning "
                           "and D4 collapsing (%d shipped, %d re-derived)"
                           % (fn, rel(ref), len(got_leaves), len(want_leaves)))
                    continue
                off = _multiset_off(want_edges, _edge_lengths(shipped))
                if off or len(want_edges) != len(_edge_lengths(shipped)):
                    R.fail("%s: branch lengths differ from those re-derived from "
                           "%s (%d differ)" % (fn, rel(ref), len(off)))
            else:
                R.fail("%s: unknown reference rule %r in %s"
                       % (fn, rule, rel(TREE_REFERENCES)))

    # SELF-TEST.  The check must prove, on every run, that it detects the
    # corruption it exists for: one passing tree is corrupted in memory the way
    # the relabelling step did it (zeros after the decimal point dropped), and
    # the comparison must flag it.
    R.bump()
    if self_test is None:
        R.fail("BRANCHLEN self-test could not run: no passing tree has a branch "
               "length with a zero after the decimal point")
    else:
        fn, a, shipped = self_test
        corrupted = _ZERO_AFTER_POINT.sub(r"\1.\2", shipped, count=1)
        bad, _first = _relabel_diff(a, _bl_tokens(corrupted), old2new)
        if corrupted == shipped or not bad:
            R.fail("BRANCHLEN self-test: a copy of %s with one branch length's "
                   "post-decimal zeros dropped was NOT flagged -- the check is "
                   "blind to the defect it exists for" % fn)

    # SELF-TEST of the multiset rule: one branch length of a passing tree is
    # perturbed in memory by 1% and the comparison must flag it.
    R.bump()
    if self_test_ms is None:
        R.note("BRANCHLEN multiset self-test not run: no edge_length_multiset_"
               "equal reference is present on this run (Data/raw is not "
               "distributed)")
    else:
        fn, ea, eb = self_test_ms
        k = max(range(len(eb)), key=lambda i: abs(eb[i]))
        bent = sorted(eb[:k] + [eb[k] * 1.01 + 1e-9] + eb[k + 1:])
        if not _multiset_off(ea, bent):
            R.fail("BRANCHLEN multiset self-test: %s with one branch length "
                   "changed by 1%% was NOT flagged" % fn)


# ==========================================================================
# CHECK: DECISIONS -- D4, D7 and D8 re-derived from their registers
# ==========================================================================

def check_decisions(ctx):
    """D4: every shipped tree counts one tip per distinct species label.
    D7: no pruned tip survives in its canonical tree; the per-partition
    checked/pruned record in headline_numbers.json is re-derived.
    D8: the coverage of each named deposited chronogram is re-measured and the
    outcome re-derived from the count."""
    for fn, t in ctx.tips.items():
        conv = (t.get("convention") or "").split(";")[0].strip()
        R.eq("D4 %s tip convention" % fn, conv, "distinct_labels")
    # D4 on the canonical files: no two tips share a species key, and no tip
    # id occurs twice.  The key is re-derived here from the standardized name.
    for part in sorted(ctx.canon_prov):
        fn = ctx.canon_file[part]
        m = ctx.measured.get(fn)
        if m is None:
            continue
        ids = list(m["named_labels"])
        keys = defaultdict(list)
        for sid in ids:
            k = d4_species_key(ctx.by_id.get(sid, ""), sid)
            if k:
                keys[k].append(sid)
        dup = {k: v for k, v in keys.items() if len(v) > 1}
        R.bump()
        if dup or m["leaves"] != m["distinct"]:
            R.fail("D4 %s: %d species key(s) carried by more than one tip "
                   "(e.g. %s) and %d repeated tip id(s)"
                   % (fn, len(dup), list(dup.items())[:2],
                      m["leaves"] - m["distinct"]))
    # SELF-TEST: the re-derivation must catch a tip that should have been
    # pruned and a duplicate species that should have been collapsed.
    R.bump()
    probe = "((1:1,2:1):1,(3:1,4:1):1);"
    by = {"1": "Aus_bus", "2": "Aus_bus_subsp._yus", "3": "Cus_dus",
          "4": "Eus_fus"}
    lv, ed = _rederive_canonical(probe, {"4"}, by)
    if lv != ["1", "3"] or sorted(ed) != [2.0, 2.0]:
        R.fail("D4/D7 self-test: re-derivation of a 4-tip probe gave %r %r "
               "instead of ['1', '3'] [2.0, 2.0]" % (lv, ed))

    head = (ctx.head.get("revision_checks") or {})
    checked = {(r.get("filename") or "").strip()
               for r in read_csv(OFFTARGET_CHECKS, required=False)}
    pruned = defaultdict(list)
    for r in read_csv(PRUNED_TIPS, required=False):
        pruned[(r.get("filename") or "").strip()].append(
            (r.get("standardized_id") or "").strip())
    canon_files = set(ctx.canon_file.values())
    for fn in pruned:
        R.bump()
        if fn not in canon_files:
            R.fail("%s lists tips pruned from %s, which is not a canonical tree"
                   % (rel(PRUNED_TIPS), fn))
    got = (head.get("offtarget_D7") or {}).get("per_partition") or {}
    unchecked = []
    for part in sorted(ctx.canon_prov):
        fn = ctx.canon_file[part]
        ids = pruned.get(fn, [])
        m = ctx.measured.get(fn)
        if m is not None:
            present = set(m["named_labels"])
            still = [i for i in ids if i in present]
            R.bump()
            if still:
                R.fail("D7 %s: %d pruned tip(s) still in the shipped tree (e.g. %s)"
                       % (fn, len(still), still[:5]))
        R.bump()
        if ids and fn not in checked:
            R.fail("D7 %s: pruned tips without a row in %s"
                   % (fn, rel(OFFTARGET_CHECKS)))
        if fn not in checked:
            unchecked.append(part)
        g = got.get(part) or {}
        R.eq("D7 %s checked" % part, g.get("checked"), fn in checked)
        R.eq("D7 %s n_pruned" % part, g.get("n_pruned"), len(ids))
    ot = head.get("offtarget_D7") or {}
    R.eq("D7 n_canonical", ot.get("n_canonical"), len(ctx.canon_prov))
    R.eq("D7 n_checked", ot.get("n_checked"), len(ctx.canon_prov) - len(unchecked))
    R.eq("D7 unchecked_partitions", ot.get("unchecked_partitions"), unchecked)
    if unchecked:
        R.note("D7: %d of %d canonical trees have not been checked for tips "
               "outside their partition: %s"
               % (len(unchecked), len(ctx.canon_prov), ", ".join(unchecked)))

    got = ((head.get("deposited_chronogram_D8") or {}).get("per_partition")
           or {})
    rows = read_csv(DATING_DEPOSIT_CHECKS, required=False)
    R.eq("D8 partitions", sorted(got), sorted((r.get("partition") or "").strip()
                                              for r in rows))
    for r in rows:
        part = (r.get("partition") or "").strip()
        g = got.get(part) or {}
        src = P("Data", (r.get("deposited_chronogram") or "").strip())
        n = ctx.labels_unique(part)
        R.eq("D8 %s shipped_tips" % part, g.get("shipped_tips"), n)
        if not (r.get("deposited_chronogram") or "").strip():
            R.eq("D8 %s outcome" % part, g.get("outcome"), "none_deposited")
        elif os.path.exists(src):
            covered = len(ctx.source_sets(src, ctx.group(part))[0]
                          & ctx.labels.get(part, set()))
            R.eq("D8 %s chronogram_tips_in_shipped" % part,
                 g.get("chronogram_tips_in_shipped"), covered)
            R.eq("D8 %s outcome" % part, g.get("outcome"),
                 "chronogram_covers_shipped_tips" if covered == n
                 else "gap_in_deposit")
            # R1-09: a deposited chronogram that lacks shipped tips is
            # recorded as such for the reader
            R.eq("D8 %s study_dating" % part,
                 (r.get("study_dating") or "").strip(),
                 "absolute_deposited" if covered == n
                 else "absolute_deposited_fewer_tips")
            R.bump()
            if covered == n and ctx.canon_file[part] \
                    and not ctx.is_geometric_chronogram(part):
                R.fail("D8 %s: the deposited chronogram covers every shipped "
                       "tip but the canonical tree is not a chronogram -- D8 "
                       "requires the swap" % part)
        else:
            R.note("D8 %s: %s is not present; its coverage was not re-measured"
                   % (part, rel(src)))


# ==========================================================================
# CHECK: NAMES -- S05, labels resolving to a CoL accepted name, re-derived
# ==========================================================================

def check_crown(ctx):
    """S02.  Crown ages and their intervals.

    1. The register (Data/crown_interval_sources.csv) and its output
       (Data/crown_intervals.csv) list the same rows; every kind is known and
       carries the C1b verdict of the table above.
    2. An 'outside' list names exactly as many labels as decision D7 pruned
       from the canonical file, so the focal node is the crown of the
       retained tips.
    3. Every summary-tree row whose source is present is RE-DERIVED here with
       scripts/crown_intervals.py's own reader and must reproduce the
       committed interval and point; its point must equal the canonical root
       height.  Every sample row whose source is present must have the
       committed sha256 (same input as the committed interval).
    4. Paper/inputs/crown_ages.csv has one row per eukaryotic canonical
       chronogram in absolute time, its crown age is the root height of the
       shipped file re-measured here, and its interval fields are the
       register's.
    5. C1b's flag in the basis follows the register (checked in BASIS)."""
    import importlib.util
    reg = read_csv(CROWN_SOURCES)
    out = OrderedDict((r["partition"], r) for r in read_csv(CROWN_INTERVALS))
    R.eq("crown register rows == crown_intervals.csv rows",
         sorted((r["partition"], r["kind"]) for r in reg),
         sorted((r["partition"], r["kind"]) for r in out.values()))
    pruned = defaultdict(int)
    for r in read_csv(PRUNED_TIPS, required=False):
        pruned[r["filename"].strip()] += 1
    spec = importlib.util.spec_from_file_location(
        "_crown_intervals", P("scripts", "crown_intervals.py"))
    ci = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ci)
    R.eq("crown_intervals.py kind table == verify's",
         {k: ("yes" if v else "no") for k, v in ci.KINDS.items()},
         CROWN_KIND_COUNTS)
    R.eq("crown_intervals.py degenerate width == verify's",
         ci.DEGENERATE_WIDTH_MYR, CROWN_DEGENERATE_MYR)
    R.eq("crown_intervals.py bound tolerance == verify's",
         ci.BOUND_TOL_MYR, CROWN_BOUND_TOL_MYR)
    R.eq("crown_intervals.py D25 containment == verify's",
         ci.PRIOR_CONTAINMENT_MIN, PRIOR_CONTAINMENT_MIN)
    R.eq("crown_intervals.py D25 width ratio == verify's",
         ci.PRIOR_WIDTH_RATIO_MIN, PRIOR_WIDTH_RATIO_MIN)
    # D25 self-test: a posterior that reproduces its prior is caught, one
    # that the data moved past a bound is not, a one-sided bound never is
    R.eq("D25 self-test reproduced",
         d25_prior_test(236.6, 249.6, {"calibration": "two_sided",
                                       "prior_lo_ma": "236.7",
                                       "prior_hi_ma": "251.2"})[5], True)
    R.eq("D25 self-test moved",
         d25_prior_test(166.7, 210.9, {"calibration": "two_sided",
                                       "prior_lo_ma": "165.2",
                                       "prior_hi_ma": "201.5"})[5], False)
    R.eq("D25 self-test narrowed",
         d25_prior_test(180.0, 190.0, {"calibration": "two_sided",
                                       "prior_lo_ma": "165.2",
                                       "prior_hi_ma": "201.5"})[5], False)
    R.eq("D25 self-test one-sided",
         d25_prior_test(498.8, 525.0, {"calibration": "one_sided"})[5], False)
    cals = OrderedDict((r["partition"], r)
                       for r in read_csv(CALIBRATION_REGISTER))
    for part, c in cals.items():
        R.eq("%s calibration class known" % part,
             c.get("calibration") in CAL_CLASSES, True)
        R.eq("%s calibration has a source and evidence" % part,
             bool(c.get("source")) and bool(c.get("evidence")), True)
        if c.get("calibration") == "two_sided":
            lo_, hi_ = as_float(c.get("prior_lo_ma")), as_float(c.get("prior_hi_ma"))
            R.eq("%s two-sided prior range numeric and ordered" % part,
                 lo_ is not None and hi_ is not None and lo_ < hi_, True)
    data_root = P("Data")
    for r in reg:
        part, kind = r["partition"], r["kind"]
        o = out.get(part) or {}
        R.bump()
        if kind not in CROWN_KIND_COUNTS:
            R.fail("%s: unknown interval kind %r" % (part, kind))
            continue
        # D5/D22, recomputed here from the committed interval
        lo, hi = as_float(o.get("interval_lo_ma")), as_float(o.get("interval_hi_ma"))
        c = as_float(o.get("interval_center_ma"))
        if c is None:
            c = as_float(o.get("source_point_height_ma"))
        degenerate = (hi - lo) < CROWN_DEGENERATE_MYR
        at_bound = (not degenerate
                    and (c - lo < CROWN_BOUND_TOL_MYR
                         or hi - c < CROWN_BOUND_TOL_MYR))
        if CROWN_KIND_COUNTS[kind] == "yes":
            R.eq("%s has a calibration-register row (D25)" % part,
                 part in cals, True)
        klass, plo, phi, cont, ratio, repro = d25_prior_test(lo, hi, cals.get(part))
        R.eq("%s calibration class" % part, o.get("calibration", ""), klass)
        R.eq("%s prior_reproduced" % part, o.get("prior_reproduced"),
             "yes" if repro else "no")
        for k, v in (("prior_lo_ma", plo), ("prior_hi_ma", phi),
                     ("prior_containment", cont), ("prior_width_ratio", ratio)):
            R.eq("%s %s" % (part, k), as_float(o.get(k)), v, atol=1e-6)
        # D29: root ages drawn from the calibration prior return the prior
        why = ("prior_reproduced" if kind == "calibration_prior_sample" else
               "kind" if CROWN_KIND_COUNTS[kind] == "no" else
               "degenerate" if degenerate else "at_bound" if at_bound
               else "prior_reproduced" if repro else "")
        R.eq("%s crown degenerate" % part, o.get("degenerate"),
             "yes" if degenerate else "no")
        R.eq("%s crown at_bound" % part, o.get("at_bound"),
             "yes" if at_bound else "no")
        R.eq("%s crown c1b_exclusion" % part, o.get("c1b_exclusion", ""), why)
        R.eq("%s crown counts_for_c1b" % part, o.get("counts_for_c1b"),
             "no" if why else "yes")
        fn = ctx.canon_file.get(part)
        outside = [x for x in r["outside"].split(";") if x.strip()]
        kept = [x for x in (r.get("focal_outgroups") or "").split(";")
                if x.strip()]
        R.bump()
        if not set(kept) <= set(outside):
            R.fail("%s: focal_outgroups is not a subset of outside" % part)
        R.eq("%s crown n_focal_outgroups" % part,
             as_int(o.get("n_focal_outgroups")), len(kept))
        # A source that IS the shipped canonical file is already pruned.
        shipped_src = r["source"].strip() == "standardized/trees/%s" % fn
        R.eq("%s crown outside labels - focal outgroups == D7 pruned tips of %s"
             % (part, fn), len(outside) - len(kept),
             0 if shipped_src else pruned.get(fn, 0))
        paths = []
        for spec_ in r["source"].split("|"):
            paths += sorted(glob.glob(os.path.join(data_root, spec_.strip())))
        if not paths:
            R.note("%s: crown-interval source %s is not present here (raw data "
                   "are not distributed); the committed interval was not "
                   "re-derived" % (part, r["source"]))
            continue
        if kind in CROWN_SUMMARY_KINDS:
            got = ci.measure(r, data_root, 1)
            for k in ("source_point_height_ma", "interval_lo_ma",
                      "interval_center_ma", "interval_hi_ma",
                      "n_focal_labels", "source_sha256",
                      "source_root_height_ma"):
                R.eq("%s crown %s re-derived" % (part, k),
                     ci.fmt(got.get(k, "")), o.get(k, ""))
            h = ctx.measured.get(fn, {}).get("root_height")
            sp = as_float(o.get("source_point_height_ma"))
            R.bump()
            if not kept:
                if h is None or abs(sp - h) > CROWN_POINT_TOL:
                    R.fail("%s: focal node of %s at %s Ma, canonical root at "
                           "%s Ma" % (part, r["source"], sp, h))
            else:
                sr = as_float(o.get("source_root_height_ma"))
                if h is None or sr is None or abs(sr - h) > CROWN_POINT_TOL \
                        or sp >= h:
                    R.fail("%s: with kept outgroups the deposit's root (%s) "
                           "must equal the canonical root (%s) and the focal "
                           "node (%s) lie below it" % (part, sr, h, sp))
        else:
            R.eq("%s crown source sha256" % part, ci.sha256_of(paths),
                 o.get("source_sha256"))
    # 4. crown_ages.csv
    clades = {r["partition"]: r for r in read_csv(CROWN_CLADES)}
    ages = OrderedDict((r["partition"], r) for r in read_csv(CROWN_AGES))
    want = sorted(p for p in ctx.canon_euk
                  if ctx.is_chronogram(p))
    R.eq("crown_ages.csv partitions == eukaryotic absolute chronograms",
         sorted(ages), want)
    for part in want:
        a = ages.get(part) or {}
        h = ctx.measured.get(ctx.canon_file[part], {}).get("root_height")
        R.eq("%s geometry_root_ma == root height of %s"
             % (part, ctx.canon_file[part]),
             as_float(a.get("geometry_root_ma")), h, atol=CROWN_VALUE_TOL)
        R.bump()
        if not (clades.get(part) or {}).get("crown_clade"):
            R.fail("%s has no crown_clade in Data/crown_clades.csv" % part)
        R.eq("%s crown_clade" % part, a.get("crown_clade"),
             (clades.get(part) or {}).get("crown_clade"))
        o = ctx.crown.get(part)
        # D20, recomputed: the deposit's own summary of the focal node
        if o and o["kind"] in CROWN_CENTER_KINDS and o.get("interval_center_ma"):
            want_age = as_float(o["interval_center_ma"])
        elif o and o["kind"] in CROWN_SUMMARY_KINDS:
            want_age = as_float(o["source_point_height_ma"])
        else:
            want_age = h
        R.eq("%s crown_age_ma by D20" % part, as_float(a.get("crown_age_ma")),
             want_age, atol=CROWN_VALUE_TOL)
        R.eq("%s crown interval kind" % part, a.get("interval_kind", ""),
             o["kind"] if o else "")
        if o:
            for k in ("interval_lo_ma", "interval_hi_ma"):
                R.eq("%s crown_ages %s" % (part, k), as_float(a.get(k)),
                     as_float(o.get(k)), atol=CROWN_VALUE_TOL)
            for k in ("degenerate", "at_bound", "counts_for_c1b",
                      "c1b_exclusion", "calibration", "prior_reproduced"):
                R.eq("%s crown_ages %s" % (part, k), a.get(k), o.get(k))
            # a reported crown age never lies outside the deposit's interval
            R.bump()
            if not (as_float(o["interval_lo_ma"]) - CROWN_VALUE_TOL
                    <= as_float(a.get("crown_age_ma"))
                    <= as_float(o["interval_hi_ma"]) + CROWN_VALUE_TOL):
                R.fail("%s: crown age %s Ma lies outside the deposit's "
                       "interval [%s, %s]" % (part, a.get("crown_age_ma"),
                                              o["interval_lo_ma"],
                                              o["interval_hi_ma"]))
    hc = ctx.head.get("crown_ages", {})
    R.eq("headline crown_ages.n_absolute_chronograms",
         hc.get("n_absolute_chronograms"), len(want))
    R.eq("headline crown_ages.n_with_counted_interval",
         hc.get("n_with_counted_interval"),
         sum(1 for p in want if ctx.interval_recoverable(p)))
    R.eq("headline crown_ages.n_degenerate_counted_intervals",
         hc.get("n_degenerate_counted_intervals"),
         sum(1 for p in want if ctx.interval_recoverable(p)
             and (ctx.crown.get(p) or {}).get("degenerate") == "yes"))


def check_union(ctx):
    """S09.  The union bound, recomputed here from the two resolution files.

    Every non-canonical eukaryotic tree has exactly its tips in
    col_resolution_noncanonical.csv, each with the prescribed query; a tip
    that adds to C0 resolves to an accepted name no canonical tree already
    counts; it is assigned to a partition by the family map measured on the
    canonical tips, else by its tree's partition (truth_dated partition_name,
    or the majority partition of its mapped tips)."""
    head = ctx.head.get("union_bound_S09") or {}
    override = {(r.get("standardized_id") or "").strip():
                (r.get("query") or "").strip()
                for r in read_csv(NAME_OVERRIDES, required=False)}
    rows = defaultdict(list)
    for r in read_csv(COL_RESOLUTION_NONCANON):
        rows[r["filename"]].append(r)
    canon_parts = [p for p in ctx.canon_euk
                   if as_int((ctx.basis.get(p) or {}).get("c0_units"), 0) > 0]
    fam = defaultdict(Counter)
    counted = set()
    for p in canon_parts:
        for r in ctx.col_rows.get(p, {}).values():
            if r.get("match_type") in COL_RESOLUTION_RESOLVED and r.get("col_id"):
                counted.add(r["col_id"])
                if r.get("col_family"):
                    fam[r["col_family"]][p] += 1
    fam_part = {f: min(c.items(), key=lambda kv: (-kv[1], kv[0]))[0]
                for f, c in fam.items()}
    prok = {p for p in ctx.species if ctx.is_prokaryotic(p)}
    noncanon = [fn for fn, d in ctx.dated.items()
                if not truthy(d.get("is_canonical"))
                and (d.get("partition_name") or "").strip() not in prok]
    R.eq("S09 non-canonical eukaryotic trees", head.get("n_noncanonical_trees"),
         len(noncanon))
    added = defaultdict(set)
    bad_ids = bad_query = 0
    for fn in sorted(noncanon):
        with open(P("Data", "standardized", "names", fn + ".json"),
                  encoding="utf-8") as fh:
            ids = set(json.load(fh))
        got = rows.get(fn, [])
        if {r["standardized_id"] for r in got} != ids or len(got) != len(ids):
            bad_ids += 1
        for r in got:
            sid = r["standardized_id"]
            want = override.get(sid) or d4_species_key(ctx.by_id.get(sid, ""),
                                                       sid) or ""
            if (r.get("query") or "") != want:
                bad_query += 1
        res = [r for r in got if r.get("match_type") in COL_RESOLUTION_RESOLVED
               and r.get("col_id")]
        tag = (ctx.dated[fn].get("partition_name") or "").strip()
        if tag in canon_parts:
            tp = tag
        else:
            c = Counter(fam_part[r["col_family"]] for r in res
                        if r["col_family"] in fam_part)
            tp = min(c.items(), key=lambda kv: (-kv[1], kv[0]))[0] if c else ""
        R.eq("S09 %s tree partition" % fn,
             ((head.get("per_file") or {}).get(fn) or {}).get("tree_partition"),
             tp)
        for r in res:
            p = fam_part.get(r["col_family"]) or tp
            if p and r["col_id"] not in counted:
                added[p].add(r["col_id"])
    R.bump()
    if bad_ids:
        R.fail("col_resolution_noncanonical.csv: %d tree(s) whose rows are not "
               "exactly their tips" % bad_ids)
    R.bump()
    if bad_query:
        R.fail("col_resolution_noncanonical.csv: %d row(s) carry a query other "
               "than the species key or the documented override" % bad_query)
    # an id assigned to two partitions would be double counted
    allids = [i for v in added.values() for i in v]
    R.eq("S09 no added name in two partitions", len(allids), len(set(allids)))
    n_add = len(set(allids))
    c0 = sum(as_int((ctx.basis.get(p) or {}).get("c0_units"), 0)
             for p in ctx.eukaryotic)
    den = sum(ctx.described(p) for p in ctx.eukaryotic)
    R.eq("S09 names_added", head.get("names_added"), n_add)
    R.eq("S09 C0_union", head.get("C0_union"), c0 + n_add)
    R.eq("S09 C0_union_pct", head.get("C0_union_pct"), pct(c0 + n_add, den),
         atol=1e-9)
    c1 = ((ctx.head.get("rungs") or {}).get("C1") or {}).get("numerator_sum")
    R.eq("S09 C1_union_low = C1", head.get("C1_union_low"), c1, atol=1e-4)
    R.eq("S09 C1_union_high = C1 + names_added", head.get("C1_union_high"),
         (c1 or 0) + n_add, atol=1e-4)
    for p, v in (head.get("per_partition") or {}).items():
        R.eq("S09 %s names added" % p, v.get("names_added_by_other_trees"),
             len(added.get(p, ())))
        R.eq("S09 %s C0_union" % p, v.get("C0_union"),
             as_int((ctx.basis.get(p) or {}).get("c0_units"), 0)
             + len(added.get(p, ())))


def check_colden(ctx):
    """S08 / D9.  The Catalogue of Life single-authority denominators.
    The counts themselves need the ChecklistBank API (scripts/
    col_denominators.py); what is checked here is that the file covers every
    partition, that filling holes never lowers a count, that the completeness
    residual is written, and that every sensitivity figure in
    headline_numbers.json follows from the file and the basis."""
    with open(COL_DENOMINATORS, newline="", encoding="utf-8") as fh:
        lines = fh.readlines()
    col = OrderedDict((r["partition"], r) for r in csv.DictReader(
        ln for ln in lines if not ln.startswith("#")))
    R.eq("col_denominators.csv partitions", sorted(col), sorted(ctx.species))
    for p, r in col.items():
        R.bump()
        if as_int(r["col_filled"], 0) < as_int(r["col_base"], 0):
            R.fail("%s: col_filled %s < col_base %s"
                   % (p, r["col_filled"], r["col_base"]))
        R.bump()
        if (as_int(r["col_filled"], 0) > as_int(r["col_base"], 0)) \
                != bool(r["filled_classes"].strip()):
            R.fail("%s: a filled count must list its filled classes" % p)
    R.bump()
    comp = [ln for ln in lines if ln.startswith("# completeness")]
    if not comp:
        R.fail("col_denominators.csv carries no completeness residual line")
    else:
        R.note(comp[0].strip("# \n"))
    h = (((ctx.head.get("sensitivities") or {}).get("denominator") or {})
         .get("col_single_authority") or {})
    rungs = ctx.head.get("rungs") or {}
    for variant, field in (("filled", "col_filled"), ("base_only", "col_base")):
        v = h.get(variant) or {}
        den, kept = 0, []
        for p in ctx.eukaryotic:
            c = as_int(col[p][field], 0)
            if c <= 0:
                kept.append(p)
            den += c if c > 0 else ctx.described(p)
        R.eq("COL %s denominator" % variant, v.get("denominator"), den)
        R.eq("COL %s partitions keeping primary" % variant,
             sorted(v.get("partitions_keeping_primary") or []), sorted(kept))
        for k, rv in (v.get("rungs") or {}).items():
            num = (rungs.get(k) or {}).get("numerator_sum")
            R.eq("COL %s %s percent" % (variant, k), rv.get("percent"),
                 pct(num, den), atol=1e-9)
            R.eq("COL %s %s 1dp" % (variant, k),
                 (v.get("rungs_pct_1dp") or {}).get(k),
                 round(pct(num, den), 1))
        for p, pv in (v.get("per_partition") or {}).items():
            c = as_int(col[p][field], 0)
            d = c if c > 0 else ctx.described(p)
            n0 = as_int((ctx.basis.get(p) or {}).get("c0_units"), 0)
            R.eq("COL %s %s C0 pct" % (variant, p), pv.get("C0_pct_col"),
                 round(pct(n0, d), 1))


def check_opentree(ctx):
    """P01.  The Open Tree comparison recomputed from Data/external (the
    extract of synthesis v16.1), the canonical tips' S05 queries and
    Data/molecular_tips.csv."""
    import gzip
    cnt = {r["clade"]: r for r in read_csv(OT_COUNTS)}
    R.eq("opentree tips from phylogeny (all life)",
         as_int(cnt["__all"]["ot_phylo"]),
         (ctx.head.get("opentree_P01") or {}).get("ot_phylo_all_life"))
    ot = defaultdict(set)
    with gzip.open(OT_SPECIES, "rt", encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            _u, nm, clade = line.rstrip("\n").split("\t")
            k = d4_species_key(nm.replace(" ", "_"), "ot")
            if k:
                ot[clade].add(k)
    mol = defaultdict(set)
    for r in read_csv(MOLTIPS):
        mol[r["partition"]].add(r["standardized_id"])
    head = ctx.head.get("opentree_P01") or {}
    per = head.get("per_partition") or {}
    csvrows = {r["partition"]: r for r in read_csv(OT_CSV)}
    want_parts = sorted(p for p in ctx.canon_euk if p in cnt
                        and as_int((ctx.basis.get(p) or {}).get("c0_units"), 0))
    R.eq("opentree partitions compared", sorted(per), want_parts)
    R.eq("opentree_comparison.csv partitions", sorted(csvrows), want_parts)
    ot_acc = defaultdict(set)
    for r in read_csv(OT_SPECIES_COL):
        if r["match_type"] in COL_RESOLUTION_RESOLVED and r["col_id"]:
            ot_acc[r["clade"]].add(r["col_id"])
    tot_acc = [0, 0]
    tot_both = 0
    for part in want_parts:
        rows = ctx.col_rows.get(part, {})
        ids = set(ctx.idmap.get(part, {}).values())

        def key(sid):
            q = (rows.get(sid) or {}).get("query")
            if q is not None:
                return q or None
            return d4_species_key(ctx.by_id.get(sid, ""), sid)

        a = {k for k in map(key, ids) if k}
        o = ot.get(part, set())
        v = per.get(part) or {}
        R.eq("opentree %s keys_atlas" % part, v.get("keys_atlas"), len(a))
        R.eq("opentree %s keys_both" % part, v.get("keys_both"), len(a & o))
        R.eq("opentree %s ot_phylo_extant_species" % part,
             v.get("ot_phylo_extant_species"),
             as_int(cnt[part]["ot_phylo_extant_species"]))
        R.eq("opentree %s csv keys_both" % part,
             as_int(csvrows[part]["keys_both"]), len(a & o))
        if part in mol:
            m = {k for k in map(key, mol[part]) if k}
            R.eq("opentree %s keys_molecular_and_ot" % part,
                 v.get("keys_molecular_and_ot"), len(m & o))
            # R1-14: the same comparison on COL26.9 accepted names
            macc = {(rows.get(sid) or {}).get("col_id") for sid in mol[part]
                    if (rows.get(sid) or {}).get("match_type")
                    in COL_RESOLUTION_RESOLVED}
            macc.discard(None)
            macc.discard("")
            oa = ot_acc.get(part, set())
            R.eq("opentree %s accepted_molecular" % part,
                 v.get("accepted_molecular"), len(macc))
            R.eq("opentree %s accepted_molecular_and_ot" % part,
                 v.get("accepted_molecular_and_ot"), len(macc & oa))
            R.eq("opentree %s accepted_molecular_not_in_ot" % part,
                 v.get("accepted_molecular_not_in_ot"), len(macc - oa))
            if ctx.basis[part].get("c1_basis") == "enumerated":
                R.eq("opentree %s accepted molecular names = C1 units" % part,
                     len(macc), as_int(ctx.basis[part]["c1_units"]))
            tot_acc[0] += len(macc)
            tot_acc[1] += len(macc & oa)
        tot_both += len(a & o)
    R.eq("opentree total keys_both", (head.get("keys") or {}).get("keys_both"),
         tot_both)
    me = head.get("molecular_enumerated") or {}
    R.eq("opentree total accepted_molecular", me.get("accepted_molecular"),
         tot_acc[0])
    R.eq("opentree total accepted_molecular_and_ot",
         me.get("accepted_molecular_and_ot"), tot_acc[1])
    tot_ot = sum(as_int(cnt[p]["ot_phylo_extant_species"]) for p in want_parts)
    R.eq("opentree OT phylo species in atlas clades",
         head.get("ot_phylo_extant_species_in_atlas_clades"), tot_ot)
    c1 = ((ctx.head.get("rungs") or {}).get("C1") or {}).get("numerator_sum")
    R.eq("opentree atlas C1 in compared partitions = rung C1",
         head.get("atlas_C1_in_compared_partitions"), c1, atol=1e-3)


WORKED_FAMILY_MIN_TIPS = 20     # D23


WORKED_CI_CUTOFF = 3.841459    # qchisq(0.95, 1): every profile interval (D30)
WORKED_LIKELIHOOD_MAX_TIPS = 15000   # independent likelihood re-derivation


def _worked_branching_times(path, keep, scale):
    """Branching times of the tree at `path` restricted to the tips in `keep`
    (None: all), scaled by `scale`; a node with k kept children contributes
    k - 1 events (worked_example_bd.R branching_events after keep.tip)."""
    with open(path, encoding="utf-8") as fh:
        root = COMMON.parse_newick(normalize_newick(fh.read()))
    out, tip_depth = [], [None]

    def walk(node, depth):
        kids, name, bl = node
        if not kids:
            ok = keep is None or name in keep
            if ok and (tip_depth[0] is None or depth > tip_depth[0]):
                tip_depth[0] = depth
            return ok
        n_kept = 0
        for k in kids:
            kb = k[2] if k[2] == k[2] else 0.0
            n_kept += 1 if walk(k, depth + kb) else 0
        if n_kept >= 2:
            out.extend([depth] * (n_kept - 1))
        return n_kept > 0
    import sys as _sys
    _sys.setrecursionlimit(max(_sys.getrecursionlimit(), 100000))
    walk(root, 0.0)
    h = tip_depth[0]          # the deepest kept tip, as node.depth.edgelength
    return sorted((max(h - d, 0.0) * scale for d in out), reverse=True)


def _worked_polytomy_counts(path, keep, scale):
    """(branching events, events at nodes with more than two children) of the
    tree at `path` restricted to `keep` (None: all), after internal branches
    shorter than ZERO_BRANCH_TOL (scaled) are collapsed, as ape::di2multi does
    in worked_example_bd.R: an independent count."""
    with open(path, encoding="utf-8") as fh:
        root = COMMON.parse_newick(normalize_newick(fh.read()))
    tot = [0, 0]

    def fin(k):
        tot[0] += k - 1
        if k > 2:
            tot[1] += k - 1

    def walk(node):
        kids, name, bl = node
        b = (bl if bl == bl and bl is not None else 0.0) * scale
        if not kids:
            return ("tip", 0, b) if (keep is None or name in keep) else None
        res = [x for x in (walk(k) for k in kids) if x is not None]
        if not res:
            return None
        if len(res) == 1:
            kind, k, ln = res[0]
            return (kind, k, ln + b)
        k_eff = 0
        for kind, k, ln in res:
            if kind == "node" and ln < ZERO_BRANCH_TOL:
                k_eff += k
            else:
                k_eff += 1
                if kind == "node":
                    fin(k)
        return ("node", k_eff, b)
    import sys as _sys
    _sys.setrecursionlimit(max(_sys.getrecursionlimit(), 100000))
    r = walk(root)
    if r and r[0] == "node":
        fin(r[1])
    return tot[0], tot[1]


def _worked_loglik(r, eps, x, rho):
    """Stadler (2009) constant-rate birth-death log-likelihood with uniform
    sampling, conditioned on the crown age and survival of both crown
    lineages: an implementation independent of worked_example_bd.R."""
    lam = r / (1.0 - eps)
    mu = eps * lam
    n = len(x) + 1
    tot = 0.0
    for i, t in enumerate(x):
        e = math.exp(-r * t)
        den = rho * lam + (lam * (1.0 - rho) - mu) * e
        lp1 = math.log(rho) + 2.0 * math.log(r) - r * t - 2.0 * math.log(den)
        tot += 2.0 * lp1 if i == 0 else lp1
        if i == 0:
            l1p0 = math.log(rho) + math.log(r) - math.log(den)
    return (n - 2) * math.log(lam) + tot - 2.0 * l1p0


def _worked_profile_eps(eps, x, rho, r0):
    """max over r of the log-likelihood at fixed epsilon (golden section on
    log r around r0)."""
    a, b = math.log(r0) - 4.0, math.log(r0) + 4.0
    g = (math.sqrt(5.0) - 1.0) / 2.0
    f = lambda q: _worked_loglik(math.exp(q), eps, x, rho)
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(120):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
    return max(fc, fd)


def check_worked_likelihood(ctx, fits, mol, sel):
    """D30.  For every uniform-sampling fit of up to
    WORKED_LIKELIHOOD_MAX_TIPS tips, the log-likelihood is recomputed here
    from the tree file with an independent implementation, and the epsilon
    interval is checked against the cutoff qchisq(0.95, 1) at each end that
    is not the parameter bound: the profile likelihood-ratio statistic at a
    reported end is 3.84 whether or not the estimate is at 0."""
    for f in fits:
        if f["rho_convention"] == "family" \
                or as_int(f["n_tips"]) > WORKED_LIKELIHOOD_MAX_TIPS:
            continue
        p = f["partition"]
        keep = None if f["tree"] == "full" else mol[p]
        x = _worked_branching_times(os.path.join(STD_TREES, f["canonical_tree"]),
                                    keep, as_float(sel[p]["time_scale_factor"]))
        R.eq("S13 %s %s %s branching times = tips - 1"
             % (p, f["tree"], f["rho_convention"]), len(x) + 1,
             as_int(f["n_tips"]))
        ne, npe = _worked_polytomy_counts(
            os.path.join(STD_TREES, f["canonical_tree"]), keep,
            as_float(sel[p]["time_scale_factor"]))
        R.eq("S13 %s %s %s events and events at polytomies, counted here"
             % (p, f["tree"], f["rho_convention"]), (ne, npe),
             (as_int(f["n_events"]), as_int(f["n_polytomy_events"])))
        rho, r0 = as_float(f["rho"]), as_float(f["r"])
        ll = _worked_profile_eps(as_float(f["epsilon"]), x, rho, r0)
        R.eq("S13 %s %s %s log-likelihood, independent implementation"
             % (p, f["tree"], f["rho_convention"]), round(ll, 2),
             round(as_float(f["loglik"]), 2), atol=0.02)
        for end in ("epsilon_lo", "epsilon_hi"):
            v = as_float(f[end])
            if v <= 0.0 or v >= 1.0 - 1e-6:
                continue
            lr = 2.0 * (as_float(f["loglik"]) - _worked_profile_eps(v, x, rho, r0))
            R.eq("S13 %s %s %s %s: profile statistic = qchisq(0.95, 1)"
                 % (p, f["tree"], f["rho_convention"], end), round(lr, 2),
                 round(WORKED_CI_CUTOFF, 2), atol=0.03)


def check_worked(ctx):
    """S13 (D23).  The worked example (Paper/code/analysis/worked_example_bd.R).
    The fits need R; what is checked here is that they were run on the
    CURRENT deposit and say what the headline says: the eligible set follows
    from the basis and the enumerated molecular tips (every eligible partition
    reported, none selected); every fitted tree has the tip count of the
    canonical file or of its molecular tips; the time scale is the reported
    crown age (D20); every sampling fraction is the one its convention
    prescribes; each estimate lies inside its own interval; the likelihood
    agrees with phytools::fit.bd; and the headline contrasts and the
    per-family molecular fractions are recomputed here."""
    mol = defaultdict(set)
    for r in read_csv(MOLTIPS):
        mol[r["partition"]].add(r["standardized_id"])
    ages = {r["partition"]: r for r in read_csv(CROWN_AGES)}
    elig = []
    for p in sorted(ctx.basis):
        b = ctx.basis[p]
        canon = set(ctx.idmap.get(p, {}).values())
        if b.get("is_prokaryotic") == "no" and b.get("is_chronogram") == "yes" \
                and b.get("time_scale") == "absolute" and p in mol \
                and len(mol[p] & canon) < len(canon):
            elig.append(p)
    sel = OrderedDict((r["partition"], r) for r in read_csv(WORKED_SEL))
    R.eq("S13 eligible partitions", list(sel), elig)
    R.bump()
    if any("selected" in r for r in sel.values()):
        R.fail("S13: the selection file still selects a partition; D23 "
               "reports every eligible partition")
    for p in elig:
        a = ages.get(p) or {}
        f = as_float(a.get("crown_age_ma")) / as_float(a.get("geometry_root_ma"))
        R.eq("S13 %s time_scale_factor = crown age / geometry root" % p,
             as_float(sel[p]["time_scale_factor"]), f, atol=1e-9)
        R.eq("S13 %s rho_C0" % p, as_float(sel[p]["rho_C0"]),
             as_int(ctx.basis[p]["c0_units"]) / as_int(ctx.basis[p]["described"]),
             atol=1e-9)
    fits = read_csv(WORKED_FITS)
    want_fits = [(p, t, c) for p in elig for t, c in
                 (("full", "tips"), ("molecular", "tips"),
                  ("molecular", "full_misapplied"), ("molecular", "family"))]
    R.eq("S13 fits: full/tips, molecular/tips, molecular/full_misapplied, "
         "molecular/family per eligible partition",
         sorted((f["partition"], f["tree"], f["rho_convention"]) for f in fits),
         sorted(want_fits))
    per = {}
    for f in fits:
        p, tree, conv = f["partition"], f["tree"], f["rho_convention"]
        per[(p, tree, conv)] = f
        if p not in elig:
            continue
        b = ctx.basis[p]
        canon = set(ctx.idmap.get(p, {}).values())
        n = len(canon) if tree == "full" else len(mol[p] & canon)
        R.eq("S13 %s %s n_tips" % (p, tree), as_int(f["n_tips"]), n)
        d = as_int(b["described"])
        if conv == "family":
            R.eq("S13 %s family fit has no single rho" % p, f["rho"], "")
            R.bump()
            uc = as_float(f.get("family_uniform_check"))
            if uc is None or uc > WORKED_UNIFORM_TOL:
                R.fail("S13 %s: the clade-specific likelihood at equal f "
                       "differs from the uniform one by %s" % (p, uc))
        else:
            rho = n / d if conv == "tips" else len(canon) / d
            R.eq("S13 %s %s %s rho" % (p, tree, conv), as_float(f["rho"]),
                 min(rho, 1.0), atol=WORKED_TOL * max(rho, 1e-12))
        # D30: one cutoff for every interval; an estimate at 0 gets [0, upper]
        bnd = as_float(f["epsilon"]) < WORKED_EPS_BOUNDARY
        R.eq("S13 %s %s %s epsilon_ci" % (p, tree, conv), f.get("epsilon_ci"),
             "estimate at 0" if bnd else "interior")
        if bnd:
            R.eq("S13 %s %s %s epsilon interval at an estimate of 0 starts at 0"
                 % (p, tree, conv), as_float(f["epsilon_lo"]), 0.0)
        if tree == "full":
            R.eq("S13 %s crown age = crown_ages.csv" % p,
                 as_float(f["crown_age_ma"]),
                 as_float((ages.get(p) or {}).get("crown_age_ma")), atol=1e-3)
        for v, lo, hi in (("r", "r_lo", "r_hi"),
                          ("epsilon", "epsilon_lo", "epsilon_hi")):
            R.bump()
            x, a, z = (as_float(f[v]), as_float(f[lo]), as_float(f[hi]))
            if None in (x, a, z) or not (a <= x <= z):
                R.fail("S13 %s %s %s: %s=%s outside its interval [%s, %s]"
                       % (p, tree, conv, v, f[v], f[lo], f[hi]))
    check_worked_likelihood(ctx, fits, mol, sel)
    for c in read_csv(WORKED_CHECK, required=False):
        for a, b in (("this_lambda", "phytools_lambda"),
                     ("this_mu", "phytools_mu")):
            R.eq("S13 %s %s %s likelihood vs phytools::fit.bd %s"
                 % (c["partition"], c["tree"], c.get("rho_convention"), a),
                 as_float(c[a]), as_float(c[b]),
                 atol=1e-3 * max(abs(as_float(c[b]) or 0.0), 1e-6))
    # A15: family sampling fractions, recomputed from the CoL resolution
    fc = {(r["kingdom"], r["family"]): as_int(r["n_accepted_species"])
          for r in read_csv(COL_FAMILY_COUNTS)}
    want = {}
    ntip_fam = {}
    for p in elig:
        tips = defaultdict(int)
        names = defaultdict(set)
        canon = set(ctx.idmap.get(p, {}).values())
        for sid in mol[p] & canon:
            r = ctx.col_rows.get(p, {}).get(sid)
            if not r or r.get("match_type") not in WORKED_RESOLVED \
                    or not (r.get("col_family") or "").strip():
                continue
            k = (r["col_kingdom"], r["col_family"])
            tips[k] += 1
            names[k].add(r["col_id"])
        nt = 0
        for k in tips:
            n_col = fc.get(k)
            fr = min(1.0, len(names[k]) / float(n_col)) if n_col else None
            want[(p,) + k] = (tips[k], len(names[k]), n_col, fr)
            nt += tips[k] if fr is not None else 0
        ntip_fam[p] = nt
    got = {}
    for r in read_csv(WORKED_FAMILY_FRACTIONS):
        got[(r["partition"], r["kingdom"], r["family"])] = (
            as_int(r["n_molecular_tips"]), as_int(r["n_accepted_names"]),
            as_int(r["n_col_species"]) if r["n_col_species"] else None,
            as_float(r["f"]))
    R.eq("S13 family-fraction rows", sorted(got), sorted(want))
    for k in sorted(set(got) & set(want)):
        g_, w_ = got[k], want[k]
        R.eq("S13 %s family counts" % "/".join(k), g_[:3], w_[:3])
        R.eq("S13 %s family f" % "/".join(k), g_[3], w_[3],
             atol=WORKED_TOL * max(w_[3] or 0, 1e-12) if w_[3] is not None else 0)
    for p in elig:
        f = per.get((p, "molecular", "family"))
        if f:
            R.eq("S13 %s tips with a family fraction" % p,
                 as_int(f.get("n_tips_family_fraction")), ntip_fam[p])
    # the headline, recomputed from the fits
    h = ctx.head.get("worked_example_S13") or {}
    R.eq("S13 headline eligible", h.get("eligible"), elig)
    G = lambda p, t, c, k: as_float(per[(p, t, c)][k])
    rhoe = h.get("rho_effect_same_tree") or {}
    imp = (h.get("full_vs_molecular") or {}).get("per_partition") or {}
    for p in elig:
        e = rhoe.get(p) or {}
        R.eq("S13 %s epsilon at rho=tips" % p, e.get("epsilon_rho_tips"),
             G(p, "molecular", "tips", "epsilon"), atol=1e-12)
        R.eq("S13 %s epsilon at rho=C0" % p, e.get("epsilon_rho_full"),
             G(p, "molecular", "full_misapplied", "epsilon"), atol=1e-12)
        R.eq("S13 %s lambda ratio" % p, e.get("lambda_ratio_full_over_tips"),
             round(G(p, "molecular", "full_misapplied", "lambda")
                   / G(p, "molecular", "tips", "lambda"), 4), atol=1e-12)
        # A15: the identity, classified independently
        rt, rc = (G(p, "molecular", "tips", "rho"),
                  G(p, "molecular", "full_misapplied", "rho"))
        lt, lc = (G(p, "molecular", "tips", "loglik"),
                  G(p, "molecular", "full_misapplied", "loglik"))
        ratio = (G(p, "molecular", "full_misapplied", "lambda")
                 / G(p, "molecular", "tips", "lambda"))
        ident = (abs(lt - lc) <= WORKED_IDENTITY_LL_TOL
                 and abs(ratio - rt / rc) <= WORKED_IDENTITY_RATIO_RTOL * rt / rc)
        R.eq("S13 %s identity classification" % p, e.get("classification"),
             "identity" if ident else "boundary")
        R.eq("S13 %s lr_statistic" % p, e.get("lr_statistic"),
             round(2 * (lt - lc), 4), atol=1e-9)
        R.eq("S13 %s rho_C0 rejected" % p, e.get("rho_full_rejected"),
             (not ident) and 2 * (lt - lc) > WORKED_LR_CUTOFF)
        # D30: the fits are nested.  In the parameters of the matching fit,
        # the fit at rho_full is the same model restricted to
        # epsilon >= 1 - rho_tips / rho_full (its mu >= 0), so the pair is a
        # boundary case exactly when the matching estimate lies below that
        # edge, and an interior pair maps epsilon by
        # 1 - eps' = (1 - eps) rho_full / rho_tips.
        et = G(p, "molecular", "tips", "epsilon")
        edge = 1.0 - rt / rc
        R.eq("S13 %s boundary iff matching epsilon below 1 - rho_tips/rho_full"
             % p, e.get("classification"),
             "boundary" if et < edge else "identity")
        R.eq("S13 %s restriction edge inside (0, 1)" % p, 0.0 < edge < 1.0, True)
        if ident:
            R.eq("S13 %s epsilon at rho_full from the nesting map" % p,
                 G(p, "molecular", "full_misapplied", "epsilon"),
                 1.0 - (1.0 - et) * rc / rt, atol=1e-4)
        R.bump()
        if not ident and G(p, "molecular", "full_misapplied", "epsilon") \
                >= WORKED_EPS_BOUNDARY:
            R.fail("S13 %s: the identity fails although epsilon at rho_C0 is "
                   "interior; the fit or the likelihood is wrong" % p)
        i = imp.get(p) or {}
        R.eq("S13 %s r ratio full/molecular" % p,
             i.get("r_ratio_full_over_molecular"),
             round(G(p, "full", "tips", "r") / G(p, "molecular", "tips", "r"),
                   4), atol=1e-12)
        ov = (G(p, "full", "tips", "r_lo") <= G(p, "molecular", "tips", "r_hi")
              and G(p, "molecular", "tips", "r_lo") <= G(p, "full", "tips", "r_hi"))
        R.eq("S13 %s r intervals overlap" % p, i.get("r_intervals_overlap"), ov)
    fvm = h.get("full_vs_molecular") or {}
    ratios = [G(p, "full", "tips", "r") / G(p, "molecular", "tips", "r")
              for p in elig]
    R.eq("S13 n_r_full_greater", fvm.get("n_r_full_greater"),
         sum(1 for x in ratios if x > 1))
    R.eq("S13 n_r_full_smaller", fvm.get("n_r_full_smaller"),
         sum(1 for x in ratios if x < 1))
    # per-family molecular fraction, recomputed from the CoL resolution
    fam_h = h.get("molecular_fraction_by_family") or {}
    for p in elig:
        fam = defaultdict(lambda: [0, 0])
        for sid, r in ctx.col_rows.get(p, {}).items():
            f = (r.get("col_family") or "").strip()
            if f:
                fam[f][0] += 1
                fam[f][1] += 1 if sid in mol[p] else 0
        fr = sorted(m / float(n) for n, m in fam.values()
                    if n >= WORKED_FAMILY_MIN_TIPS)
        got = fam_h.get(p) or {}
        R.eq("S13 %s families with >= %d tips" % (p, WORKED_FAMILY_MIN_TIPS),
             got.get("n_families"), len(fr))
        if fr:
            k = (len(fr) - 1) * 0.5
            lo = int(k)
            hi = min(lo + 1, len(fr) - 1)
            R.eq("S13 %s median family molecular fraction" % p,
                 got.get("median"), fr[lo] + (fr[hi] - fr[lo]) * (k - lo),
                 atol=1e-6)
            R.eq("S13 %s families with no molecular tip" % p,
                 got.get("n_families_no_molecular_tip"),
                 sum(1 for x in fr if x == 0.0))


def check_names(ctx):
    """S05 / D18: the per-partition resolution to CoL accepted names, re-derived
    from Data/standardized/col_resolution.csv.  The match itself needs the CoL
    release (not distributed); what is checked here is that every eukaryotic
    canonical tip has exactly one row, that each row's QUERY is the one the
    rule prescribes (species key, or a documented override), that labels which
    are not species labels are never resolved, and that every count published
    from the file follows from it."""
    head = ctx.head.get("name_resolution_S05") or {}
    got = head.get("per_partition") or {}
    R.eq("S05 partitions", sorted(got), sorted(ctx.canon_euk))
    # A10: genera with a suprageneric suffix are genera when CoL says so.
    # The committed list is re-derived from the release index when present.
    R.bump()
    if not D4_SUPRA_GENERA:
        R.fail("%s is missing or empty; the species key cannot tell a genus "
               "with a suprageneric suffix from a family name"
               % rel(D4_SUPRA_GENERA_FILE))
    if os.path.exists(COL_INDEX):
        mine = set()
        with open(COL_INDEX, encoding="utf-8") as fh:
            for line in fh:
                f = line.split("\t", 6)
                if len(f) > 5 and f[5] and f[4].lower().endswith(D4_SUPRA):
                    mine.add(f[4])
        R.eq("%s re-derived from the CoL index" % rel(D4_SUPRA_GENERA_FILE),
             sorted(D4_SUPRA_GENERA), sorted(mine))
    else:
        R.note("the CoL index is not present; %s was not re-derived"
               % rel(D4_SUPRA_GENERA_FILE))
    for probe, want in (("Xerophyta_suaveolens", "Xerophyta_suaveolens"),
                        ("Elapsoidea_laticincta", "Elapsoidea_laticincta"),
                        ("Pimephales_vigilax", "Pimephales_vigilax"),
                        ("Smittoidea_albula_BLEED1708", "Smittoidea_albula"),
                        ("Rosaceae_sp", None), ("Formicidae_gen", None),
                        ("Smittoidea_aff_reticulata", None)):
        R.eq("D4 species key self-test %s" % probe,
             d4_species_key(probe, probe), want)
    override = {(r.get("standardized_id") or "").strip():
                (r.get("query") or "").strip()
                for r in read_csv(NAME_OVERRIDES, required=False)}
    for r in read_csv(NAME_OVERRIDES, required=False):
        R.bump()
        if not (r.get("evidence") or "").strip():
            R.fail("%s: override for %s cites no evidence"
                   % (rel(NAME_OVERRIDES), r.get("standardized_id")))
    den_euk = sum(ctx.described(p) for p in ctx.eukaryotic)
    t_lab = t_res = t_acc = 0
    union = set()
    csv_rows = {r["partition"]: r for r in read_csv(C0_NAME_RESOLUTION,
                                                    required=False)}
    bad_query = bad_unres = 0
    for part in sorted(ctx.canon_euk):
        rows = ctx.col_rows.get(part, {})
        ids = set(ctx.idmap.get(part, {}).values())
        acc, res = set(), 0
        for sid in ids:
            r = rows.get(sid) or {}
            want = override.get(sid) or d4_species_key(ctx.by_id.get(sid, ""),
                                                       sid) or ""
            if (r.get("query") or "") != want:
                bad_query += 1
            if not want and r.get("match_type") != "unresolved_label":
                bad_unres += 1
            if r.get("match_type") in COL_RESOLUTION_RESOLVED and r.get("col_id"):
                res += 1
                acc.add(r["col_id"])
        n = len(ids)
        g = got.get(part) or {}
        R.eq("S05 %s labels" % part, g.get("labels"), n)
        R.eq("S05 %s resolved_labels" % part, g.get("resolved_labels"), res)
        R.eq("S05 %s distinct_accepted_names" % part,
             g.get("distinct_accepted_names"), len(acc))
        R.eq("S05 %s distinct accepted names = basis c0_units" % part,
             as_int((ctx.basis.get(part) or {}).get("c0_units")), len(acc))
        c = csv_rows.get(part) or {}
        R.eq("c0_name_resolution.csv %s resolved_labels" % part,
             as_int(c.get("resolved_labels")), res)
        R.eq("c0_name_resolution.csv %s unresolved_labels" % part,
             as_int(c.get("unresolved_labels")), n - res)
        t_lab += n
        t_res += res
        t_acc += len(acc)
        union |= acc
    R.bump()
    if bad_query:
        R.fail("col_resolution.csv: %d row(s) carry a query other than the "
               "species key or the documented override" % bad_query)
    R.bump()
    if bad_unres:
        R.fail("col_resolution.csv: %d label(s) that are not species labels "
               "(D13) carry a match" % bad_unres)
    e = head.get("eukaryotic") or {}
    R.eq("S05 eukaryotic unresolved numerator",
         (e.get("unresolved_share") or {}).get("numerator"), t_lab - t_res)
    R.eq("S05 eukaryotic labels", (e.get("unresolved_share") or {}).get(
        "denominator"), t_lab)
    R.eq("S05 C0_resolved_labels", (e.get("C0_resolved_labels") or {}).get(
        "numerator"), t_res)
    R.eq("S05 C0 denominator", (e.get("C0_resolved_labels") or {}).get(
        "denominator"), den_euk)
    R.eq("S05 C0_distinct_accepted_names_sum",
         (e.get("C0_distinct_accepted_names_sum") or {}).get("numerator"), t_acc)
    R.eq("S05 C0_distinct_accepted_names_union",
         (e.get("C0_distinct_accepted_names_union") or {}).get("numerator"),
         len(union))
    R.eq("S05 C0 rung numerator = distinct accepted names (union)",
         ((ctx.head.get("rungs") or {}).get("C0") or {}).get("numerator_union"),
         len(union))


# ==========================================================================
# CHECK: FILE-TRUTH  (the ledgers of record vs geometry)
# ==========================================================================

def ledgers_of_record(ctx):
    """(label, {filename -> 'yes'/'no'/''}) for every ledger that DISTRIBUTES a
    dating claim.  The pre-audit ledgers are handled separately: they disagree
    with geometry on purpose, and that disagreement is the paper's finding."""
    out = []
    out.append(("Data/_rebuild/truth_dated.csv",
                {f: tri(r.get("is_chronogram")) for f, r in ctx.dated.items()}))
    basis_dated = {}
    for r in ctx.basis_rows:
        fn = (r.get("canonical_tree") or "").strip()
        if fn:
            basis_dated[fn] = tri(r.get("is_chronogram"))
    out.append(("Data/coverage_basis.csv", basis_dated))
    ps = {}
    for r in read_csv(PARTITION_SUMMARY, required=False):
        fn = (r.get("filename") or "").strip()
        if fn:
            ps[fn] = tri(r.get("dated"))
    out.append(("Data/site/partition_summary.csv", ps))
    return out


def check_file_truth(ctx):
    canon_files = set(ctx.canon_file.values())
    for label, flags in ledgers_of_record(ctx):
        for fname, flag in sorted(flags.items()):
            m = ctx.measured.get(fname)
            if m is None:
                continue          # --quick, or a file already failed to parse
            if flag == "":
                continue
            R.bump()
            # D14: a ledger of record says dated=yes only for an absolute-time
            # chronogram; a relative-time ultrametric file is not dated.
            want = ("yes" if m["chronogram"]
                    and (label == "Data/_rebuild/truth_dated.csv"
                         or ctx.time_scale.get(fname, "absolute") == "absolute")
                    else "no")
            if flag != want:
                R.fail("%s carries %s as dated=%s, but its root-to-tip CV is "
                       "%.6f (tolerance %.2f) -> %s"
                       % (label, fname, flag, m["cv"] if m["cv"] is not None
                          else float("nan"), CHRONOGRAM_TOLERANCE, want))
            # No NaN branch length is tolerated in a tree a ledger of record
            # carries as dated AND that the atlas distributes as a partition's
            # canonical tree: the parser scores a NaN as zero length, so the
            # ultrametricity verdict would rest on a fabricated branch.
            if flag == "yes" and fname in canon_files and m["nan"]:
                R.bump()
                R.fail("%s carries canonical tree %s as dated, but the file "
                       "has %d NaN branch length(s); the ultrametricity "
                       "verdict rests on branches the parser scored as zero"
                       % (label, fname, m["nan"]))

    # The count of dated canonical trees is a headline in its own right: the
    # ledgers claimed 17/16/18, the files say what they say.  Measured here,
    # never typed.
    measured_chrono = sorted(p for p in ctx.canon_prov
                             if ctx.canon_file[p] in ctx.measured
                             and ctx.is_chronogram(p))
    relative = sorted(p for p in ctx.canon_prov
                      if ctx.canon_file[p] in ctx.measured
                      and ctx.is_geometric_chronogram(p)
                      and ctx.time_scale.get(ctx.canon_file[p]) == "relative")
    R.eq("headline trees.relative_time_partitions",
         sorted((ctx.head.get("trees") or {}).get("relative_time_partitions")
                or []), relative)
    # D14: every ultrametric canonical file declares its time scale; a file
    # declared absolute must be older than MIN_ABSOLUTE_ROOT_MA at its root.
    for p in ctx.canon_prov:
        fn = ctx.canon_file[p]
        m = ctx.measured.get(fn)
        if m is None or not ctx.is_geometric_chronogram(p):
            continue
        ts = ctx.time_scale.get(fn)
        R.bump()
        if ts not in ("absolute", "relative"):
            R.fail("D14 %s is ultrametric but %s gives no absolute/relative "
                   "time scale" % (fn, rel(TIME_SCALE_CHECKS)))
        R.bump()
        if ts == "absolute" and (m.get("root_height") or 0.0) < \
                MIN_ABSOLUTE_ROOT_MA:
            R.fail("D14 %s is declared absolute but its root is at %.4f, below "
                   "%.1f Ma" % (fn, m.get("root_height") or 0.0,
                                MIN_ABSOLUTE_ROOT_MA))
    trees = ctx.head.get("trees", {})
    R.eq("headline trees.n_partition_canonical_chronograms",
         trees.get("n_partition_canonical_chronograms"), len(measured_chrono))
    R.eq("headline trees.chronogram_partitions",
         sorted(trees.get("chronogram_partitions") or []), measured_chrono)
    R.eq("headline trees.n_partition_canonical",
         trees.get("n_partition_canonical"), len(ctx.canon_prov))
    R.eq("headline trees.n_partition_canonical_eukaryotic",
         trees.get("n_partition_canonical_eukaryotic"), len(ctx.canon_euk))
    if not ctx.quick:
        R.eq("headline trees.n_shipped", trees.get("n_shipped"),
             len(ctx.dated))
        R.eq("headline trees.n_chronograms_all_shipped",
             trees.get("n_chronograms_all_shipped"),
             sum(1 for m in ctx.measured.values() if m["chronogram"]))

    check_withdrawn(ctx)


# Places a withdrawn tree must not reach.  Site copies are all named tree.nwk,
# so a copy is recognised by CONTENT (a known sha256, or a Newick file with the
# withdrawn leaf count), never by basename.  Data/raw and Data/curated are not
# distributed (.gitignore) and the withdrawn directory itself is excluded.
WITHDRAWN_SCAN_ROOTS = (("Data", "standardized"), ("Data", "site"),
                        ("Data", "_archive"), ("Paper",))
WITHDRAWN_DIR = P("Data", "_archive", "2026-10-revision",
                  "withdrawn_not_redistributed")
NEWICK_EXT = (".nwk", ".tre", ".tree", ".newick", ".nex", ".trees")


def _sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _newick_leaf_count(path):
    """Leaves of the first tree in a Newick/NEXUS file: commas + 1 inside the
    first '(' ... ';' span.  Approximate for NEXUS translate blocks, exact for
    plain Newick, which is what a copy of a shipped tree is."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    i = text.find("(")
    if i < 0:
        return None
    j = text.find(";", i)
    return text.count(",", i, j if j > 0 else len(text)) + 1


def check_withdrawn(ctx):
    rows = [r for r in read_csv(WITHDRAWN, required=False)
            if (r.get("filename") or "").strip()]
    md_files = {(r.get("filename") or "").strip()
                for r in read_csv(STD_METADATA, required=False)}
    hashes, ntips, groups = set(), set(), set()
    for r in rows:
        hashes |= {h.strip() for h in (r.get("sha256") or "").split(";")
                   if h.strip()}
        n = as_int(r.get("ntips"))
        if n:
            ntips.add(n)
        g = (r.get("group") or "").strip()
        if g:
            groups.add(g)
        R.bump()
        if not hashes or not n or not g:
            R.fail("%s: %s must give group, ntips and sha256 so that copies "
                   "can be recognised by content" % (rel(WITHDRAWN),
                                                     r["filename"]))
    try:
        site = json.load(open(P("Data", "site", "data.json"), encoding="utf-8"))
    except (OSError, ValueError):
        site = {}
    site_trees = site.get("trees") or []
    for r in rows:
        fn = r["filename"].strip()
        g = (r.get("group") or "").strip()
        for label, present in (
                ("Data/standardized/trees", os.path.exists(
                    os.path.join(STD_TREES, fn))),
                ("Data/standardized/names", os.path.exists(
                    P("Data", "standardized", "names", fn + ".json"))),
                ("Data/_rebuild/truth_dated.csv", fn in ctx.dated),
                ("Data/_rebuild/truth_tips.csv", fn in ctx.tips),
                ("Data/standardized/metadata.csv", fn in md_files),
                ("Data/site/data.json", any(
                    t.get("filename") == fn or (g and t.get("group") == g)
                    for t in site_trees))):
            R.bump()
            if present:
                R.fail("%s is withdrawn from distribution (%s) but is still "
                       "present in %s" % (fn, rel(WITHDRAWN), label))
    # The name layer carries only labels of shipped trees.
    n_fm = 0
    for r in read_csv(FULL_MAPPING, required=False):
        if (r.get("group") or "").strip() in groups:
            n_fm += 1
    R.bump()
    if n_fm:
        R.fail("full_mapping.csv still has %d row(s) of a withdrawn group (%s)"
               % (n_fm, ", ".join(sorted(groups))))
    name_ids = set()
    ndir = P("Data", "standardized", "names")
    for fn in os.listdir(ndir):
        if fn.endswith(".json"):
            name_ids |= set(json.load(open(os.path.join(ndir, fn),
                                           encoding="utf-8")).keys())
    # ... plus the labels of the distributed archived trees a canonical file
    # is built from (Data/_archive/2026-10-revision/canonical_sources/) and of
    # the superseded canonical files kept with their name maps.
    rev = P("Data", "_archive", "2026-10-revision")
    for d in (os.path.join(rev, "canonical_sources"),):
        if os.path.isdir(d):
            for fn in os.listdir(d):
                with open(os.path.join(d, fn), encoding="utf-8") as fh:
                    name_ids |= set(_leaf_order(fh.read()))
    sup = os.path.join(rev, "superseded_canonicals", "names")
    if os.path.isdir(sup):
        for fn in os.listdir(sup):
            name_ids |= set(json.load(open(os.path.join(sup, fn),
                                           encoding="utf-8")).keys())
    dict_ids = set(ctx.by_id)
    R.bump()
    if dict_ids != name_ids:
        R.fail("dictionary.csv ids != the labels of the distributed trees "
               "(shipped name maps, canonical sources, superseded canonicals): "
               "%d only in the dictionary, %d only in the trees"
               % (len(dict_ids - name_ids), len(name_ids - dict_ids)))
    col_ids = {r["standardized_id"] for r in read_csv(COL_RESOLUTION,
                                                      required=False)}
    R.bump()
    if col_ids - dict_ids:
        R.fail("col_resolution.csv names %d id(s) absent from dictionary.csv"
               % len(col_ids - dict_ids))
    # Copies anywhere in the distributed tree, recognised by content.
    if not rows:
        return
    for root_parts in WITHDRAWN_SCAN_ROOTS:
        base = P(*root_parts)
        for d, dirs, files in os.walk(base):
            if os.path.abspath(d).startswith(os.path.abspath(WITHDRAWN_DIR)):
                dirs[:] = []
                continue
            for f in files:
                path = os.path.join(d, f)
                if not f.lower().endswith(NEWICK_EXT):
                    continue
                R.bump()
                if _sha256(path) in hashes:
                    R.fail("%s is a byte copy of a withdrawn tree" % rel(path))
                    continue
                if os.path.getsize(path) < 4 * min(ntips):   # >= "(a,"-size tips
                    continue
                if _newick_leaf_count(path) in ntips:
                    R.fail("%s has the leaf count of a withdrawn tree (%s); "
                           "inspect it" % (rel(path),
                                           ", ".join(map(str, sorted(ntips)))))
    # SELF-TEST: the content rule must recognise a renamed site-style copy.
    R.bump()
    probe = P("Data", "_rebuild", ".withdrawn_probe_tree.nwk")
    src = os.path.join(WITHDRAWN_DIR, "Data", "standardized", "trees",
                       rows[0]["filename"].strip())
    if os.path.exists(src):
        # A changed byte defeats the hash; the leaf-count rule must still hit.
        with open(src, "rb") as fi, open(probe, "wb") as fo:
            fo.write(b"\n" + fi.read())
        try:
            hit = (_sha256(probe) not in hashes
                   and _newick_leaf_count(probe) in ntips)
        finally:
            os.remove(probe)
        if not hit:
            R.fail("WITHDRAWN self-test: a byte-modified copy of the withdrawn "
                   "tree was NOT recognised by its leaf count")
    else:
        R.note("WITHDRAWN self-test skipped: the withdrawn file is not present "
               "in this checkout (it is never distributed)")


# ==========================================================================
# CHECK: PREAUDIT -- the ledgers disagree; ground truth must SAY SO correctly
# ==========================================================================

def preaudit_source(relpath):
    """Prefer the frozen pre-audit snapshot; the live file may already have
    been regenerated, in which case comparing truth_dated.csv's record of the
    pre-audit state against it would be meaningless."""
    snap = os.path.join(PRE_AUDIT, relpath)
    if os.path.exists(snap):
        return snap, True
    return P(*relpath.split("/")), False


def check_preaudit(ctx):
    tp_path, tp_frozen = preaudit_source("Data/tree_provenance.csv")
    md_path, md_frozen = preaudit_source("Data/standardized/metadata.csv")
    ps_path, ps_frozen = preaudit_source("Data/site/partition_summary.csv")

    tp = {r["filename"]: tri(r.get("dated"))
          for r in read_csv(tp_path, required=False) if r.get("filename")}
    md = {r["filename"]: tri(r.get("dated"))
          for r in read_csv(md_path, required=False) if r.get("filename")}
    ps = {}
    for r in read_csv(ps_path, required=False):
        fn = (r.get("filename") or "").strip()
        if fn:
            ps[fn] = tri(r.get("dated"))

    if not (tp_frozen and md_frozen and ps_frozen):
        R.note("at least one pre-audit ledger was read from the LIVE tree "
               "rather than the frozen snapshot; if the rebuild has already "
               "regenerated it, truth_dated.csv's flag columns record a state "
               "that no longer exists on disk",
               ["tree_provenance frozen=%s" % tp_frozen,
                "standardized/metadata frozen=%s" % md_frozen,
                "site/partition_summary frozen=%s" % ps_frozen])

    sources = (("flag_tree_provenance", tp),
               ("flag_standardized_metadata", md),
               ("flag_partition_summary", ps))

    disagree = []
    for fname, row in ctx.dated.items():
        seen = set()
        for col, table in sources:
            got = tri(row.get(col))
            want = table.get(fname, "")
            R.bump()
            if got != want:
                R.fail("truth_dated.csv %s for %s says %r; the ledger says %r "
                       "-- ground truth mis-describes the pre-audit state"
                       % (col, fname, got, want))
            if want:
                seen.add(want)
        R.bump()
        want_dis = "yes" if len(seen) > 1 else "no"
        if tri(row.get("ledgers_disagree")) != want_dis:
            R.fail("truth_dated.csv ledgers_disagree for %s says %r; the "
                   "pre-audit flags are %s" % (fname,
                                               tri(row.get("ledgers_disagree")),
                                               sorted(seen) or ["-"]))
        m = ctx.measured.get(fname)
        if m is not None and seen and any(f != ("yes" if m["chronogram"]
                                                else "no") for f in seen):
            disagree.append(fname)

    canon = set(ctx.canon_file.values())
    canon_dis = sorted(f for f in disagree if f in canon)
    R.note("%d shipped tree(s), of which %d partition-canonical, carry a "
           "pre-audit `dated` flag that contradicts the geometry of the file. "
           "THIS IS THE PAPER'S FINDING, not a defect in the deposit: it is "
           "why time calibration is measured and why the ledger's `dated` "
           "column is retired as a source of truth."
           % (len(disagree), len(canon_dis)),
           ["canonical: %s" % ", ".join(canon_dis[:20])] if canon_dis else None)


# ==========================================================================
# CHECK: INVARIANTS
# ==========================================================================

def check_invariants(ctx):
    for part, sp in ctx.species.items():
        d = as_int(sp.get("described"))
        lo = as_int(sp.get("estimated_low"))
        hi = as_int(sp.get("estimated_high"))
        if d is not None and hi is not None:
            R.bump()
            if d > hi:
                R.fail("%s: described %s exceeds estimated_high %s -- more "
                       "species are described than are estimated to exist"
                       % (part, commas(d), commas(hi)))
        if lo is not None and hi is not None:
            R.bump()
            if lo > hi:
                R.fail("%s: estimated_low %s exceeds estimated_high %s"
                       % (part, commas(lo), commas(hi)))
        if not ctx.is_unrepresented(part):
            tips = ctx.labels_unique(part)
            R.bump()
            if d is None or d <= 0:
                R.fail("%s has a canonical tree but described=%r; a coverage "
                       "percentage against it would be undefined"
                       % (part, sp.get("described")))
            elif tips > d:
                R.fail("%s: %s tips exceed %s described species -- coverage "
                       "above 100 %%. The retired min(x,100) cap in "
                       "prep_data.py existed to hide exactly this; it is "
                       "deleted, so this must be fixed at the source"
                       % (part, commas(tips), commas(d)))

    # The basis table's own bounds.
    for r in ctx.basis_rows:
        part = r["partition"]
        f = as_float(r.get("f_mol"))
        lo = as_float(r.get("f_low"))
        hi = as_float(r.get("f_high"))
        if f is None:
            continue
        R.bump()
        if f < 0.0 or f > 1.0:
            R.fail("%s: f_mol=%r is outside [0,1]" % (part, f))
        if lo is not None:
            R.bump()
            if lo > f + ATOL_FRAC:
                R.fail("%s: f_low %r exceeds f_mol %r" % (part, lo, f))
        if hi is not None:
            R.bump()
            if hi < f - ATOL_FRAC:
                R.fail("%s: f_high %r is below f_mol %r" % (part, hi, f))
        cv = as_float(r.get("roottotip_cv"))
        if cv is not None:
            R.bump()
            if cv < 0.0:
                R.fail("%s: negative root-to-tip CV %r" % (part, cv))

    # estimated_low/high sums quoted in the headline must be sums, not typed.
    al = (ctx.head.get("denominators") or {}).get("all_life") or {}
    R.eq("denominators.all_life.estimated_low", al.get("estimated_low"),
         sum(as_int(r.get("estimated_low"), 0) for r in ctx.species_rows))
    R.eq("denominators.all_life.estimated_high", al.get("estimated_high"),
         sum(as_int(r.get("estimated_high"), 0) for r in ctx.species_rows))


# ==========================================================================
# CHECK: BASIS -- the 62-row deposit contract
# ==========================================================================

BASIS_HEADER = [
    "partition", "group", "category", "is_prokaryotic", "described",
    "described_source", "is_unrepresented", "canonical_tree", "labels_raw",
    "labels_unique", "c0_units", "c0_unresolved_labels", "f_mol", "f_basis",
    "f_low", "f_high", "c1_units", "c1_basis",
    "molecular_source_file", "is_chronogram", "time_scale", "roottotip_cv",
    "interval_recoverable", "interval_kind", "posterior_n",
    "topology_assembly", "placement_evidence",
]

TREE_SIDE_FIELDS = ["canonical_tree", "f_mol", "f_basis", "f_low", "f_high",
                    "c1_units", "c1_basis",
                    "molecular_source_file", "is_chronogram", "time_scale",
                    "roottotip_cv",
                    "interval_recoverable", "interval_kind", "posterior_n",
                    "topology_assembly", "placement_evidence"]


def check_basis(ctx):
    R.eq("coverage_basis.csv column order", csv_header(BASIS), BASIS_HEADER)
    R.eq("coverage_basis.csv row count", len(ctx.basis_rows),
         len(ctx.species_rows))
    R.eq("coverage_basis.csv partition set", sorted(ctx.basis),
         sorted(ctx.species))
    R.eq("coverage_basis.csv row order", [r["partition"] for r in ctx.basis_rows],
         sorted(r["partition"] for r in ctx.basis_rows))

    prok = [r for r in ctx.basis_rows if tri(r.get("is_prokaryotic")) == "yes"]
    R.eq("prokaryotic partitions in the basis", len(prok),
         len(ctx.species) - len(ctx.eukaryotic))
    R.eq("eukaryotic partitions in the basis",
         sum(1 for r in ctx.basis_rows
             if tri(r.get("is_prokaryotic")) == "no"), len(ctx.eukaryotic))

    ledger = {r.get("filename"): r for r in read_csv(TREE_PROVENANCE,
                                                     required=False)}

    for part, r in ctx.basis.items():
        sp = ctx.species[part]
        R.eq("%s basis group" % part, (r.get("group") or "").strip(),
             (sp.get("group") or "").strip())
        R.eq("%s basis category" % part, (r.get("category") or "").strip(),
             (sp.get("category") or "").strip())
        R.eq("%s basis described" % part, as_int(r.get("described")),
             ctx.described(part))
        R.eq("%s basis described_source" % part,
             (r.get("described_source") or "").strip(),
             (sp.get("described_source") or "").strip())
        R.eq("%s basis is_prokaryotic" % part, tri(r.get("is_prokaryotic")),
             "yes" if ctx.is_prokaryotic(part) else "no")
        R.eq("%s basis is_unrepresented" % part, tri(r.get("is_unrepresented")),
             "yes" if ctx.is_unrepresented(part) else "no")

        if ctx.is_unrepresented(part):
            R.bump()
            if as_int(r.get("labels_unique"), 0) != 0:
                R.fail("%s is unrepresented but labels_unique=%r; a structural "
                       "zero must be a row with zero tips, not an absent row"
                       % (part, r.get("labels_unique")))
            for f in TREE_SIDE_FIELDS:
                R.bump()
                if (r.get(f) or "").strip():
                    R.fail("%s is unrepresented but carries %s=%r"
                           % (part, f, r.get(f)))
            continue

        fname = ctx.canon_file[part]
        R.eq("%s basis canonical_tree" % part,
             (r.get("canonical_tree") or "").strip(), fname)
        R.eq("%s basis labels_raw" % part, as_int(r.get("labels_raw")),
             ctx.labels_raw(part))
        R.eq("%s basis labels_unique" % part, as_int(r.get("labels_unique")),
             ctx.labels_unique(part))
        R.eq("%s basis c0_units" % part, as_int(r.get("c0_units")),
             len(ctx.c0_names.get(part, {})))
        R.eq("%s basis c0_unresolved_labels" % part,
             as_int(r.get("c0_unresolved_labels")),
             ctx.c0_unresolved.get(part, 0))
        if fname in ctx.measured:
            R.eq("%s basis time_scale" % part,
                 (r.get("time_scale") or "").strip(),
                 ctx.time_scale.get(fname, "")
                 if ctx.is_geometric_chronogram(part) else "not_ultrametric")
        if fname in ctx.measured:
            R.eq("%s basis is_chronogram" % part, tri(r.get("is_chronogram")),
                 "yes" if ctx.is_chronogram(part) else "no")
            # Against the ground truth EXACTLY (both are published at FRAC_DP,
            # so the strings must be identical) and against this script's own
            # re-measurement to within the recheck tolerance.
            R.eq("%s basis roottotip_cv vs truth_dated.csv" % part,
                 (r.get("roottotip_cv") or "").strip(),
                 fnum(as_float(ctx.dated[fname].get("roottotip_cv")), FRAC_DP))
            R.eq("%s basis roottotip_cv vs re-measurement" % part,
                 as_float(r.get("roottotip_cv")), ctx.measured[fname]["cv"],
                 atol=CV_RECHECK_TOL)
        R.eq("%s basis interval_recoverable" % part,
             tri(r.get("interval_recoverable")),
             "yes" if ctx.interval_recoverable(part) else "no")
        R.eq("%s basis posterior_n" % part, (r.get("posterior_n") or "").strip(),
             ctx.posterior_n(part))
        R.eq("%s basis interval_kind" % part,
             (r.get("interval_kind") or "").strip(),
             (ctx.crown[part]["kind"] if ctx.interval_recoverable(part)
              else ""))

        # interval_kind must be empty exactly when nothing is recoverable, so
        # the per-partition table prints "--" and never a misleading zero.
        R.bump()
        has_kind = bool((r.get("interval_kind") or "").strip())
        if has_kind != ctx.interval_recoverable(part):
            R.fail("%s: interval_kind=%r but interval_recoverable=%s"
                   % (part, r.get("interval_kind"),
                      ctx.interval_recoverable(part)))

        lr = ledger.get(fname) or {}
        R.eq("%s topology_assembly" % part,
             (r.get("topology_assembly") or "").strip(),
             classify_assembly(lr.get("otl_pipeline_involvement"),
                               lr.get("inference")))
        # placement_evidence must ship the raw triple so a reader can overturn
        # the classifier in four lines.
        pe = (r.get("placement_evidence") or "")
        for field in ("otl_pipeline_involvement", "inference", "tree_class"):
            val = (lr.get(field) or "").strip()
            if val:
                R.bump()
                if val not in pe:
                    R.fail("%s: placement_evidence=%r does not carry the raw "
                           "%s=%r; the annotation must be overturnable from "
                           "the deposit alone" % (part, pe, field, val))

    # The interval flag is MEASURED (Data/crown_intervals.csv); the ledger's
    # archival_class is prose.  Where the two disagree, say so: the ledger is
    # the one that is wrong, and the reader of truth_provenance.csv must know.
    differ = sorted(p for p in ctx.canon_prov
                    if (ctx.archival_class(p) == "hpd_recoverable")
                    != ctx.interval_recoverable(p))
    if differ:
        R.note("%d canonical partition(s) where truth_provenance.csv "
               "archival_class disagrees with the measured interval register "
               "(Data/crown_intervals.csv), which C1b follows: %s"
               % (len(differ), ", ".join(differ)))


# ==========================================================================
# CHECK: FMOL -- re-measure axis A from the deposited source analyses
# ==========================================================================

def molecular_sentences(note):
    return [s for s in SENTENCE_SPLIT_RE.split(note or "")
            if "molecular_tip_fraction" in s]


def recorded_numerator(sents):
    for m in RATIO_RE.finditer(" ".join(sents)):
        n, d = as_int(m.group(1)), as_int(m.group(2))
        if n is not None and d is not None and 0 < n <= d:
            return n
    return None


def named_candidates(sents):
    named = []
    for s in sents:
        for p in RAW_PATH_RE.findall(s):
            if p not in named:
                named.append(p)
    cands = [P("Data", p) for p in named if p.lower().endswith(SOURCE_TREE_EXT)]
    if named and not cands:
        seen = set()
        for p in named:
            d = os.path.dirname(P("Data", p))
            if d in seen or not os.path.isdir(d):
                continue
            seen.add(d)
            for fn in sorted(os.listdir(d)):
                if fn.lower().endswith(SOURCE_TREE_EXT):
                    cands.append(os.path.join(d, fn))
    return cands


def best_source(ctx, cands, shipped_names, group, recorded,
                excluded=frozenset()):
    """The deposited tree that best evidences the molecular subset, under the
    deposit's own two acceptance rules.  Deterministic; a rerun cannot flip
    the choice."""
    scored = []
    for path in cands:
        if not os.path.exists(path):
            continue
        src, _verbatim = ctx.source_sets(path, group)
        if not src:
            continue
        # D4: a shipped tip stands for its species; compare species keys.
        # D7: tips pruned from the canonical tree are known to lie outside it.
        src = d4_keys(set(src) - excluded)
        shk = d4_keys(shipped_names)
        inter = len(src & shk)
        if inter == 0:
            continue
        if len(src - shk) / float(len(src)) > SOURCE_NONSUBSET_TOL:
            continue
        if recorded is not None:
            dev = abs(inter - recorded) / float(recorded)
            if dev > SOURCE_NUMERATOR_TOL:
                continue
            scored.append((dev, path, inter))
        else:
            scored.append((0.0, path, inter))
    if not scored:
        return None
    scored.sort(key=lambda t: (t[0], t[1]))
    dev, path, inter = scored[0]
    return inter, path


def derive_molecular(ctx, part):
    """Re-derive (f_mol, f_basis, f_low, f_high, source_file, inter, exact)
    for one canonical partition, from the deposit alone.  This is the
    deposit's own four-step ladder, reimplemented so that what the paper calls
    a measurement is verified to be one."""
    pv = ctx.canon_prov[part]
    sents = molecular_sentences(pv.get("conflict_note") or "")
    shipped = ctx.labels.get(part, set())
    n_shipped = ctx.labels_unique(part)
    group = ctx.group(part)
    ledger_f = as_float(pv.get("molecular_tip_fraction"))
    imputation = (pv.get("imputation_status") or "").strip()

    rec = recorded_numerator(sents)
    fn = ctx.canon_file[part]
    excluded = frozenset(ctx.by_id.get((r.get("standardized_id") or "").strip(),
                                       "")
                         for r in read_csv(PRUNED_TIPS, required=False)
                         if (r.get("filename") or "").strip() == fn)
    best = best_source(ctx, named_candidates(sents), shipped, group, rec,
                       excluded)

    f = None
    basis = "unresolved"
    src_file = ""
    inter = exact = None
    contested_f = None

    if best is not None:
        inter, path = best
        exact = len(ctx.source_sets(path, group)[1] & shipped)
        cand_f = inter / float(n_shipped) if n_shipped else None
        asserted_one = (ledger_f is not None and ledger_f >= 1.0
                        and imputation == "none")
        if asserted_one and cand_f is not None and cand_f < 1.0:
            # CONTESTED.  The ground truth ruled for 1.00 and flagged the
            # disagreement; the ruling stands and the implied value becomes the
            # lower end of the band.
            contested_f = cand_f
        else:
            f, src_file = cand_f, rel(path)
            basis = ("source_alignment_intersection"
                     if path.endswith(TAXA_LIST_EXT) else
                     "source_tree_intersection")

    # 1b. per-tip evidence: every row must name a shipped tip; a
    # source_tree_label row must name a label of its source tree when that
    # tree is present.  Rows that fail are not counted.
    ev_names = None
    rows_ev = [r for r in ctx.tip_evidence
               if (r.get("partition") or "").strip() == part]
    if f is None and rows_ev and n_shipped:
        ev_names = set()
        for r in rows_ev:
            nm = ctx.by_id.get((r.get("standardized_id") or "").strip())
            kind = (r.get("evidence_kind") or "").strip()
            if kind not in MOL_EVIDENCE_KINDS or nm is None or nm not in shipped:
                continue
            if kind == "source_tree_label":
                src = P("Data", (r.get("source_file") or "").strip())
                if os.path.exists(src) and (r.get("source_label") or "").strip() \
                        not in ctx.source_sets(src, group)[1]:
                    continue
            ev_names.add(nm)
        f, basis = len(ev_names) / float(n_shipped), "per_tip_evidence"
        src_file = rel(MOL_TIP_EVIDENCE)
        inter = len(ev_names)

    if f is None and rec is not None and n_shipped:
        f, basis = min(1.0, rec / float(n_shipped)), "published_fraction"
    if f is None and imputation == "none":
        f, basis = 1.0, "asserted_no_imputation"
    if f is None:
        f, basis = ledger_f, "unresolved"

    lo = hi = f
    if contested_f is not None and lo is not None:
        lo, hi = min(lo, contested_f), max(hi, contested_f)
    for s in sents:
        if "UNRESOLVED" not in s:
            continue
        for v in (as_float(x) for x in FRACTION_RE.findall(s)):
            if v is not None and 0.0 <= v <= 1.0 and lo is not None:
                lo, hi = min(lo, v), max(hi, v)
    return OrderedDict([("f_mol", f), ("f_basis", basis), ("f_low", lo),
                        ("f_high", hi), ("source_file", src_file),
                        ("measured_intersection", inter),
                        ("exact_label_intersection", exact),
                        ("recorded_numerator", rec),
                        ("contested_f", contested_f),
                        ("evidence_names", ev_names)])


def check_fmol(ctx):
    ev = ctx.head.get("molecular_evidence") or {}
    ctx.derived = OrderedDict()
    for part in ctx.canon_prov:
        if part not in ctx.labels:
            continue                                  # --quick / parse failure
        d = derive_molecular(ctx, part)
        ctx.derived[part] = d
        row = ctx.basis.get(part) or {}
        n = ctx.labels_unique(part)

        R.eq("%s f_basis" % part, (row.get("f_basis") or "").strip(),
             d["f_basis"])
        R.eq("%s f_mol" % part, (row.get("f_mol") or "").strip(),
             fnum(d["f_mol"], FRAC_DP))
        R.eq("%s f_low" % part, (row.get("f_low") or "").strip(),
             fnum(d["f_low"], FRAC_DP))
        R.eq("%s f_high" % part, (row.get("f_high") or "").strip(),
             fnum(d["f_high"], FRAC_DP))
        R.eq("%s molecular_source_file" % part,
             (row.get("molecular_source_file") or "").strip(), d["source_file"])

        # The measured intersection IS the numerator.  Assert the arithmetic
        # closes: f_mol * labels_unique must be the tip count, not a scalar
        # inherited from a ledger.
        if d["f_basis"] in ("source_tree_intersection",
                            "source_alignment_intersection"):
            R.eq("%s f_mol = |shipped n source| / labels_unique" % part,
                 (row.get("f_mol") or "").strip(),
                 fnum(d["measured_intersection"] / float(n), FRAC_DP))
            e = ev.get(part) or {}
            R.eq("%s molecular_evidence.measured_intersection" % part,
                 e.get("measured_intersection"), d["measured_intersection"])
            R.eq("%s molecular_evidence.exact_label_intersection" % part,
                 e.get("exact_label_intersection"),
                 d["exact_label_intersection"])
            if d["recorded_numerator"] is not None \
                    and d["recorded_numerator"] != d["measured_intersection"]:
                R.note("%s: the MEASURED intersection is %s against the %s "
                       "recorded in truth_provenance.csv. The measurement is "
                       "the published value (the recorded ratio is on a "
                       "different tip basis); truth_provenance.csv's "
                       "molecular_tip_fraction for this partition should carry "
                       "a note saying so."
                       % (part, commas(d["measured_intersection"]),
                          commas(d["recorded_numerator"])))
            if d["exact_label_intersection"] is not None \
                    and d["exact_label_intersection"] != d["measured_intersection"]:
                R.note("%s: the adopted rule (both label sets resolved through "
                       "the deposit's own name mapping) measures %s molecular "
                       "tips; naive exact-label matching measures %s. Both are "
                       "published in headline_numbers.json -> "
                       "molecular_evidence, so either can be reproduced."
                       % (part, commas(d["measured_intersection"]),
                          commas(d["exact_label_intersection"])))
        elif d["f_basis"] == "per_tip_evidence":
            rows_ev = [r for r in ctx.tip_evidence
                       if (r.get("partition") or "").strip() == part]
            R.bump()
            bad_kind = sorted({(r.get("evidence_kind") or "").strip()
                               for r in rows_ev} - set(MOL_EVIDENCE_KINDS))
            if bad_kind:
                R.fail("%s: unknown evidence_kind(s) %s in %s"
                       % (part, bad_kind, rel(MOL_TIP_EVIDENCE)))
            R.bump()
            if len(d["evidence_names"]) != len(rows_ev):
                R.fail("%s: %d per-tip evidence rows, %d of them name a shipped "
                       "tip with checkable evidence" % (part, len(rows_ev),
                                                        len(d["evidence_names"])))
            absent = sorted({(r.get("source_file") or "").strip() for r in rows_ev
                             if (r.get("evidence_kind") or "").strip()
                             == "source_tree_label"
                             and not os.path.exists(P("Data", (r.get(
                                 "source_file") or "").strip()))})
            if absent:
                R.note("%s: the source tree(s) %s cited by per-tip evidence are "
                       "not present (Data/raw is not distributed); the cited "
                       "labels were not re-checked on this run"
                       % (part, ", ".join(absent)))
        elif d["f_basis"] == "asserted_no_imputation":
            R.eq("%s asserted_no_imputation => f_mol" % part,
                 as_float(row.get("f_mol")), 1.0, atol=ATOL_FRAC)
            R.eq("%s asserted_no_imputation => imputation_status" % part,
                 (ctx.canon_prov[part].get("imputation_status") or "").strip(),
                 "none")
        elif d["f_basis"] == "published_fraction":
            R.eq("%s published_fraction rebased onto labels_unique" % part,
                 (row.get("f_mol") or "").strip(),
                 fnum(min(1.0, d["recorded_numerator"] / float(n)), FRAC_DP))
            R.note("%s: f_mol is the source publication's fraction rebased "
                   "onto the %s distinct labels the atlas counts. AN UPPER "
                   "BOUND, and the weakest measurement in the ladder."
                   % (part, commas(n)))
        else:
            R.note("%s: f_mol is UNRESOLVED -- the pre-audit ledger scalar is "
                   "carried through. This partition must be named in the "
                   "Results disclosures." % part)

        if d["contested_f"] is not None:
            R.note("%s: a deposited tree named in the ground truth implies "
                   "f_mol=%.6f against the asserted 1.00. The asserted value "
                   "is kept as the point estimate and the implied value is "
                   "published as f_low with its own sensitivity line."
                   % (part, d["contested_f"]))


# ==========================================================================
# THE LADDER, recomputed from measurements alone
# ==========================================================================

def measured_molecular_names(ctx, part, d):
    """The shipped names of `part` that are in its molecular analysis, measured
    here from the source (D4 species keys) or the per-tip evidence; None where
    the partition names no molecular source."""
    if not d["source_file"]:
        return None
    labels = ctx.labels.get(part, set())
    if d["f_basis"] == "per_tip_evidence":
        return {nm for nm in labels if nm in (d["evidence_names"] or set())}
    src_path = P(d["source_file"])
    if not os.path.exists(src_path):
        return None
    src_keys = d4_keys(ctx.source_sets(src_path, ctx.group(part))[0])
    return {nm for nm in labels if d4_keys([nm]) & src_keys}


def measured_c1(ctx, part, d):
    """(c1_units, c1_basis, set or None), measured: distinct counting units
    (CoL accepted names; GTDB clusters for prokaryotes) among the molecular
    tips (audit A01)."""
    names = measured_molecular_names(ctx, part, d)
    c0 = ctx.c0_names.get(part, {})
    if names is not None:
        idmap = ctx.idmap.get(part, {})
        unit_of = ctx.unit_of.get(part, {})
        units = {unit_of[idmap[nm]] for nm in names
                 if nm in idmap and idmap[nm] in unit_of}
        return float(len(units)), "enumerated", units
    if d["f_mol"] is not None and d["f_mol"] == 1.0:
        return float(len(c0)), "all_tips", set(c0)
    return round(len(c0) * (pub(d["f_mol"]) or 0.0), 4), "expectation", None


def my_basis(ctx):
    """The basis this script measured for itself.  It is built from the truth
    files, the shipped trees and the deposited source analyses -- never from
    Data/coverage_basis.csv, which is the file under test."""
    rows = []
    for part in sorted(ctx.species):
        unrep = ctx.is_unrepresented(part)
        row = OrderedDict([
            ("partition", part),
            ("category", (ctx.species[part].get("category") or "").strip()),
            ("described", ctx.described(part)),
            ("is_prokaryotic", "yes" if ctx.is_prokaryotic(part) else "no"),
            ("labels_unique", 0), ("c0_units", 0), ("f_mol", None),
            ("f_low", None),
            ("f_high", None), ("f_basis", ""), ("c1_units", None),
            ("c1_basis", ""), ("is_chronogram", "no"),
            ("interval_recoverable", "no"), ("topology_assembly", ""),
        ])
        if not unrep and part in ctx.derived:
            d = ctx.derived[part]
            lr = ctx.ledger_rows.get(ctx.canon_file[part]) or {}
            row["labels_unique"] = ctx.labels_unique(part)
            row["c0_units"] = len(ctx.c0_names.get(part, {}))
            # The PUBLISHED precision is the operative value: coverage_basis.csv
            # ships f_mol at FRAC_DP places, so the deposited file must
            # reproduce the published totals exactly.
            row["f_mol"] = pub(d["f_mol"])
            row["f_low"] = pub(d["f_low"])
            row["f_high"] = pub(d["f_high"])
            row["f_basis"] = d["f_basis"]
            c1, c1_basis, c1_set = measured_c1(ctx, part, d)
            row["c1_units"], row["c1_basis"] = c1, c1_basis
            ctx.c1_sets[part] = c1_set
            row["is_chronogram"] = "yes" if ctx.is_chronogram(part) else "no"
            row["interval_recoverable"] = ("yes" if ctx.interval_recoverable(part)
                                           else "no")
            row["topology_assembly"] = classify_assembly(
                lr.get("otl_pipeline_involvement"), lr.get("inference"))
        rows.append(row)
    return rows


def weight(row, rung):
    """w(rung).  Nested conjuncts: each rung is the previous one AND one more
    measured condition, so the weights are non-increasing BY CONSTRUCTION."""
    if as_int(row.get("c0_units"), 0) <= 0:
        return 0.0
    w = 1.0
    if rung in ("C1", "C1d", "C1b"):
        c0, c1 = as_int(row.get("c0_units"), 0), as_float(row.get("c1_units"))
        # C1 is the count of distinct accepted names among molecular tips;
        # c1_units=None selects the superseded product c0_units * f_mol
        w *= (c1 / float(c0)) if c1 is not None \
            else (as_float(row.get("f_mol")) or 0.0)
    if rung in ("C1d", "C1b"):
        w *= 1.0 if row.get("is_chronogram") == "yes" else 0.0
    if rung == "C1b":
        w *= 1.0 if row.get("interval_recoverable") == "yes" else 0.0
    return w


def is_euk(row):
    return row.get("is_prokaryotic") == "no"


def rungs_of(ctx, rows, denominator=None, keep=None):
    euk = [r for r in rows if is_euk(r)]
    den = denominator if denominator is not None \
        else sum(as_int(r["described"], 0) for r in euk)
    contributing = [r for r in euk if as_int(r.get("c0_units"), 0) > 0
                    and (keep is None or keep(r))]
    out = OrderedDict()
    for rung in RUNGS:
        total = round(sum(as_int(r["c0_units"], 0) * weight(r, rung)
                          for r in contributing), 4)
        e = OrderedDict([
            ("numerator_sum", total),
            ("numerator_sum_rounded", int(round(total))),
            ("denominator", den),
            ("percent", pct(total, den)),
            ("n_partitions_contributing",
             sum(1 for r in contributing if weight(r, rung) > 0)),
        ])
        if rung != "C0":
            live = [r for r in contributing if weight(r, rung) > 0]
            sets = [ctx.c1_sets.get(r["partition"]) for r in live]
            if live and all(x is not None
                            and as_float(r.get("c1_units")) == float(len(x))
                            for r, x in zip(live, sets)):
                union = set().union(*sets)
                e["numerator_union"] = len(union)
                e["percent_union"] = pct(len(union), den)
                e["shared_labels"] = int(round(total)) - len(union)
        if rung == "C0":
            sets = [set(ctx.c0_names[r["partition"]]) for r in contributing
                    if r["partition"] in ctx.c0_names]
            union = set().union(*sets) if sets else set()
            e["numerator_union"] = len(union)
            e["percent_union"] = pct(len(union), den)
            e["shared_labels"] = int(round(total)) - len(union)
        out[rung] = e
    return out


EXACT_KEYS = ("numerator_sum_rounded", "numerator_union", "shared_labels",
              "n_partitions_contributing", "denominator")

# The sensitivity blocks publish a reduced shape: `numerator` there is the
# ROUNDED numerator sum.  Mapping it explicitly is what lets one comparator
# serve both shapes without either side guessing.
SENS_KEYMAP = {"numerator": "numerator_sum_rounded",
               "percent": "percent",
               "numerator_union": "numerator_union",
               "percent_union": "percent_union"}


def cmp_rungs(label, got, want, keys=("numerator_sum", "percent"), keymap=None):
    for rung in RUNGS:
        g = (got or {}).get(rung) or {}
        w = want[rung]
        for k in keys:
            wk = (keymap or {}).get(k, k)
            if wk not in w:
                continue
            atol = ATOL_PCT if "percent" in k else ATOL_SUM
            if wk in EXACT_KEYS:
                atol = 0.0
            R.eq("%s %s.%s" % (label, rung, k), g.get(k), w[wk], atol=atol)


# ==========================================================================
# CHECK: HEADLINES
# ==========================================================================

def check_headlines(ctx):
    head = ctx.head
    rows = ctx.my_rows
    euk = [r for r in rows if is_euk(r)]
    prok_rows = [r for r in rows if not is_euk(r)]
    den_euk = sum(as_int(r["described"], 0) for r in euk)
    den_prok = sum(as_int(r["described"], 0) for r in prok_rows)

    d = head.get("denominators") or {}
    R.eq("denominators.eukaryotic.described",
         (d.get("eukaryotic") or {}).get("described"), den_euk)
    R.eq("denominators.eukaryotic.n_partitions",
         (d.get("eukaryotic") or {}).get("n_partitions"), len(euk))
    R.eq("denominators.prokaryotic.described",
         (d.get("prokaryotic") or {}).get("described"), den_prok)
    R.eq("denominators.all_life.described",
         (d.get("all_life") or {}).get("described"), den_euk + den_prok)
    R.eq("denominators.all_life.n_partitions",
         (d.get("all_life") or {}).get("n_partitions"), len(rows))

    mine = rungs_of(ctx, rows)
    cmp_rungs("rungs", head.get("rungs"), mine,
              keys=("numerator_sum", "numerator_sum_rounded", "denominator",
                    "percent", "n_partitions_contributing", "numerator_union",
                    "percent_union", "shared_labels"))

    # Every rung's percent must be reproducible from the numerator and
    # denominator printed beside it.  This is the property the retired 23.7 %
    # literal did not have.
    for rung in RUNGS:
        e = (head.get("rungs") or {}).get(rung) or {}
        n, dd = e.get("numerator_sum"), e.get("denominator")
        R.bump()
        if n is None or dd in (None, 0):
            R.fail("rungs.%s has no numerator/denominator pair; a coverage "
                   "figure that is reachable from no pair is exactly the "
                   "defect this deposit exists to retire" % rung)
        elif not close(e.get("percent"), pct(n, dd), ATOL_PCT):
            R.fail("rungs.%s percent=%r is not 100*%r/%r"
                   % (rung, e.get("percent"), n, dd))

    # ---- per category, eukaryotic only -------------------------------
    cats = head.get("per_category_eukaryotic") or {}
    order = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    want_cats = sorted({r["category"] for r in euk},
                       key=lambda c: (order.get(c, len(order)), c))
    R.eq("per_category_eukaryotic keys", list(cats.keys()), want_cats)
    for cat in want_cats:
        rows_c = [r for r in euk if r["category"] == cat]
        dd = sum(as_int(r["described"], 0) for r in rows_c)
        got = cats.get(cat) or {}
        R.eq("category %s described" % cat, got.get("described"), dd)
        R.eq("category %s n_partitions" % cat, got.get("n_partitions"),
             len(rows_c))
        for rung, sk, pk in (("C0", "C0_sum", "C0_pct"),
                             ("C1", "C1_sum", "C1_pct"),
                             ("C1d", "C1d_sum", "C1d_pct")):
            v = sum(as_int(r["c0_units"], 0) * weight(r, rung)
                    for r in rows_c)
            R.eq("category %s %s" % (cat, sk), got.get(sk), int(round(v)))
            R.eq("category %s %s" % (cat, pk), got.get(pk), pct(v, dd),
                 atol=ATOL_PCT)

    # ---- prokaryotes, on their own line, never mixed ------------------
    prok = head.get("prokaryotic_line") or {}
    p_named = sum(as_int(r["c0_units"], 0) for r in prok_rows)
    R.eq("prokaryotic_line.described", prok.get("described"), den_prok)
    R.eq("prokaryotic_line.n_partitions", prok.get("n_partitions"),
         len(prok_rows))
    R.eq("prokaryotic_line.C0_sum", prok.get("C0_sum"), p_named)
    R.eq("prokaryotic_line.C0_pct", prok.get("C0_pct"), pct(p_named, den_prok),
         atol=ATOL_PCT)
    for r in prok_rows:
        pp = (prok.get("per_partition") or {}).get(r["partition"]) or {}
        R.eq("prokaryote %s c0_units" % r["partition"],
             pp.get("c0_units"), as_int(r["c0_units"], 0))
        R.eq("prokaryote %s is_chronogram" % r["partition"],
             tri(pp.get("is_chronogram")), r["is_chronogram"])
        R.eq("prokaryote %s C1d_sum" % r["partition"], pp.get("C1d_sum"),
             int(round(as_int(r["c0_units"], 0) * weight(r, "C1d"))))
    # A prokaryote must never appear in a eukaryotic numerator.
    R.bump()
    if any(r["partition"] in {x["partition"] for x in euk} for r in prok_rows):
        R.fail("a prokaryotic partition appears among the eukaryotic rows")

    # ---- archival classes and the assembly annotation -----------------
    arch = defaultdict(int)
    for part in ctx.canon_prov:
        arch[ctx.archival_class(part)] += 1
    R.eq("archival_classes_canonical",
         OrderedDict(sorted((head.get("archival_classes_canonical")
                             or {}).items())),
         OrderedDict(sorted((k, v) for k, v in arch.items() if k)))
    ta = head.get("topology_assembly") or {}
    want_counts = OrderedDict()
    for c in sorted({r["topology_assembly"] for r in euk
                     if r["topology_assembly"]}):
        want_counts[c] = sum(1 for r in euk if r["topology_assembly"] == c)
    R.eq("topology_assembly.counts",
         OrderedDict(sorted((ta.get("counts") or {}).items())),
         OrderedDict(sorted(want_counts.items())))
    R.bump()
    if ">=" in str(ta.get("classifier", "")) \
            or "strong" in str(ta.get("classifier", "")).lower():
        R.fail("topology_assembly is described with an ordering; it is a "
               "NOMINAL annotation and must carry no >= and no strong/weak "
               "reading")
    R.bump()
    if "annotation" not in str(ta.get("source", "")).lower():
        R.fail("topology_assembly.source does not state that it is an "
               "ANNOTATION derived from unreconciled ledger fields")

    # ---- disclosures ---------------------------------------------------
    check_disclosures(ctx, rows, mine)
    # ---- sensitivities --------------------------------------------------
    check_sensitivities(ctx, rows, mine)
    # ---- figure inputs --------------------------------------------------
    check_figure_inputs(ctx, rows)


def check_disclosures(ctx, rows, mine):
    disc = ctx.head.get("disclosures") or {}
    euk = [r for r in rows if is_euk(r)]
    represented = [r for r in euk if as_int(r.get("c0_units"), 0) > 0]
    den = sum(as_int(r["described"], 0) for r in euk)

    def contrib(row, rung):
        return as_int(row["c0_units"], 0) * weight(row, rung)

    # dominance of the single largest contributor
    top = max(represented, key=lambda r: (contrib(r, "C0"), r["partition"]))
    dom = disc.get("dominance") or {}
    R.eq("disclosures.dominance.partition", dom.get("partition"),
         top["partition"])
    for rung in RUNGS:
        R.eq("dominance.share_of_%s_pct" % rung,
             dom.get("share_of_%s_pct" % rung),
             pct(contrib(top, rung), mine[rung]["numerator_sum"]),
             atol=ATOL_PCT)

    # f_basis accounting of C1
    by = defaultdict(float)
    for r in represented:
        by[r["f_basis"]] += contrib(r, "C1")
    acc = disc.get("f_basis_accounting_of_C1") or {}
    R.eq("f_basis_accounting keys", sorted(acc.keys()), sorted(by.keys()))
    c1 = mine["C1"]["numerator_sum"]
    for b, v in by.items():
        got = acc.get(b) or {}
        R.eq("f_basis %s tip_equivalents" % b, got.get("tip_equivalents"),
             int(round(v)))
        R.eq("f_basis %s share_of_C1_pct" % b, got.get("share_of_C1_pct"),
             pct(v, c1), atol=ATOL_PCT)
    R.bump()
    if abs(sum(by.values()) - c1) > ATOL_SUM:
        R.fail("the f_basis accounting sums to %r but C1 is %r; the "
               "accounting must partition C1 exactly"
               % (sum(by.values()), c1))

    # overlapping labels -- the sum/union difference
    where = defaultdict(list)
    shown = {}
    for r in represented:
        for lab, disp in ctx.c0_names.get(r["partition"], {}).items():
            where[lab].append(r["partition"])
            shown[lab] = disp
    dupes = sorted(shown[l] for l, ps in where.items() if len(ps) > 1)
    ov = disc.get("overlapping_labels") or {}
    R.eq("overlapping_labels.n_labels_in_two_or_more_partitions",
         ov.get("n_labels_in_two_or_more_partitions"), len(dupes))
    R.eq("overlapping_labels examples", sorted((ov.get("examples") or {}).keys()),
         dupes)
    R.eq("overlapping_labels.share_of_C0_union_pct",
         ov.get("share_of_C0_union_pct"),
         pct(len(dupes), mine["C0"]["numerator_union"]), atol=ATOL_PCT)
    # the arithmetic that guarantees sum and union can never drift
    R.eq("C0 sum - union == shared labels",
         mine["C0"]["numerator_sum_rounded"] - mine["C0"]["numerator_union"],
         len(dupes))

    # structural zeros
    zeros = [r for r in euk if as_int(r.get("c0_units"), 0) == 0]
    sz = disc.get("structural_zeros") or {}
    R.eq("structural_zeros.n_partitions", sz.get("n_partitions"), len(zeros))
    R.eq("structural_zeros.of_eukaryotic_partitions",
         sz.get("of_eukaryotic_partitions"), len(euk))
    R.eq("structural_zeros.described", sz.get("described"),
         sum(as_int(r["described"], 0) for r in zeros))
    R.eq("structural_zeros.share_of_denominator_pct",
         sz.get("share_of_denominator_pct"),
         pct(sum(as_int(r["described"], 0) for r in zeros), den), atol=ATOL_PCT)
    R.eq("structural_zeros.partitions", sorted(sz.get("partitions") or []),
         sorted(r["partition"] for r in zeros))

    # dark matter -- formed on the UNION numerator throughout
    dm = disc.get("dark_matter") or {}
    c0u = mine["C0"]["numerator_union"]
    R.eq("dark_matter.no_tip_anywhere", dm.get("no_tip_anywhere"), den - c0u)
    R.eq("dark_matter.no_tip_anywhere_pct", dm.get("no_tip_anywhere_pct"),
         pct(den - c0u, den), atol=ATOL_PCT)
    R.eq("dark_matter.no_molecular_placement", dm.get("no_molecular_placement"),
         den - int(round(c1)))
    R.eq("dark_matter.no_molecular_placement_pct",
         dm.get("no_molecular_placement_pct"), pct(den - c1, den), atol=ATOL_PCT)
    R.bump()
    if "union" not in str(dm.get("convention", "")).lower():
        R.fail("dark_matter does not declare the UNION convention; the sum "
               "and the union differ by %d labels and every derived quantity "
               "must use one of them consistently" % len(dupes))

    # the ordering fix
    ofix = {o.get("partition"): o for o in (disc.get("ordering_fix") or [])}
    R.eq("ordering_fix partition set", sorted(ofix.keys()),
         sorted(r["partition"] for r in represented))
    for r in represented:
        o = ofix.get(r["partition"]) or {}
        dd = as_int(r["described"], 0)
        R.eq("ordering_fix %s f_mol" % r["partition"], o.get("f_mol"),
             r["f_mol"], atol=ATOL_FRAC)
        R.eq("ordering_fix %s C0_pct" % r["partition"], o.get("C0_pct"),
             pct(as_int(r["c0_units"], 0), dd), atol=ATOL_PCT)
        R.eq("ordering_fix %s C1_pct" % r["partition"], o.get("C1_pct"),
             pct(contrib(r, "C1"), dd), atol=ATOL_PCT)
        R.eq("ordering_fix %s old_strict_pct" % r["partition"],
             o.get("old_strict_pct"),
             0.0 if r["topology_assembly"] != ASSEMBLY_PRIMARY
             else pct(contrib(r, "C1"), dd), atol=ATOL_PCT)


def check_sensitivities(ctx, rows, mine):
    sens = ctx.head.get("sensitivities") or {}
    euk = [r for r in rows if is_euk(r)]

    # -- 1/4. f_mol bands -------------------------------------------------
    bands = sens.get("f_mol_bands") or {}
    want_banded = sorted(r["partition"] for r in euk
                         if as_int(r.get("c0_units"), 0) > 0
                         and r["f_mol"] is not None and r["f_low"] is not None
                         and r["f_high"] is not None
                         and r["f_high"] > r["f_low"])
    R.eq("sensitivities.f_mol_bands partition set", sorted(bands.keys()),
         want_banded)
    R.bump()
    if not want_banded:
        R.fail("no banded partition was found. The deposit's largest open "
               "uncertainty is an UNRESOLVED molecular fraction; if nothing "
               "bands, the band parser has stopped matching it")
    for part in want_banded:
        row = next(r for r in euk if r["partition"] == part)
        got = bands.get(part) or {}
        R.eq("band %s f_point" % part, got.get("f_point"), row["f_mol"],
             atol=ATOL_FRAC)
        R.eq("band %s f_low" % part, got.get("f_low"), row["f_low"],
             atol=ATOL_FRAC)
        R.eq("band %s f_high" % part, got.get("f_high"), row["f_high"],
             atol=ATOL_FRAC)
        for tag, val in (("at_f_low", row["f_low"]), ("at_f_high", row["f_high"])):
            alt = [dict(r) for r in rows]
            for r in alt:
                if r["partition"] == part:
                    if r.get("c1_units") is not None and r["f_mol"]:
                        r["c1_units"] = round(r["c1_units"] * pub(val)
                                              / r["f_mol"], 4)
                    r["f_mol"] = pub(val)
            cmp_rungs("band %s %s" % (part, tag), got.get(tag),
                      rungs_of(ctx, alt),
                      keys=("numerator", "percent", "numerator_union",
                            "percent_union"), keymap=SENS_KEYMAP)
        R.eq("band %s propagates_to_C1d" % part, got.get("propagates_to_C1d"),
             row["is_chronogram"] == "yes")
        if row["is_chronogram"] != "yes":
            R.note("%s is banded but is NOT a measured chronogram, so the band "
                   "moves C1 and cannot move C1d or C1b. The figures that "
                   "carry the paper's argument are immune to the deposit's "
                   "largest open uncertainty -- state this explicitly."
                   % part)

    # -- 2. leave-assembly-class-out --------------------------------------
    la = sens.get("leave_assembly_class_out") or {}
    classes = sorted({r["topology_assembly"] for r in euk
                      if r["topology_assembly"]})
    R.eq("leave_assembly_class_out.classes_present",
         OrderedDict(sorted((la.get("classes_present") or {}).items())),
         OrderedDict(sorted((c, sum(1 for r in euk
                                    if r["topology_assembly"] == c))
                            for c in classes)))
    for c in classes:
        # Both variants NAME the partitions of class c: retain_only_ keeps
        # exactly them, drop_ removes exactly them.
        in_class = sorted(r["partition"] for r in euk
                          if as_int(r.get("c0_units"), 0) > 0
                          and r["topology_assembly"] == c)
        variants = (
            ("retain_only_", "retained_partitions",
             (lambda cc: (lambda r: r["topology_assembly"] == cc))(c)),
            ("drop_", "dropped_partitions",
             (lambda cc: (lambda r: r["topology_assembly"] != cc))(c)),
        )
        for prefix, key, keep in variants:
            got = la.get(prefix + c) or {}
            R.eq("%s%s %s" % (prefix, c, key), sorted(got.get(key) or []),
                 in_class)
            cmp_rungs("%s%s" % (prefix, c), got.get("rungs"),
                      rungs_of(ctx, rows, keep=keep),
                      keys=("numerator", "percent", "numerator_union",
                            "percent_union"), keymap=SENS_KEYMAP)

    # -- 3. denominator sensitivity and elasticity ------------------------
    ds = sens.get("denominator") or {}
    den_full = sum(as_int(r["described"], 0) for r in euk)
    represented = [r for r in euk if as_int(r.get("c0_units"), 0) > 0]
    den_repr = sum(as_int(r["described"], 0) for r in represented)

    R.eq("denominator.full_denominator.described",
         (ds.get("full_denominator") or {}).get("described"), den_full)
    cmp_rungs("denominator.full", (ds.get("full_denominator") or {}).get("rungs"),
              rungs_of(ctx, rows),
              keys=("numerator", "percent", "numerator_union",
                    "percent_union"), keymap=SENS_KEYMAP)
    R.eq("denominator.represented_partitions_only.described",
         (ds.get("represented_partitions_only") or {}).get("described"), den_repr)
    R.eq("denominator.represented_partitions_only.n_partitions",
         (ds.get("represented_partitions_only") or {}).get("n_partitions"),
         len(represented))
    cmp_rungs("denominator.represented_only",
              (ds.get("represented_partitions_only") or {}).get("rungs"),
              rungs_of(ctx, rows, denominator=den_repr),
              keys=("numerator", "percent", "numerator_union",
                    "percent_union"), keymap=SENS_KEYMAP)

    biggest = max(euk, key=lambda r: (as_int(r["described"], 0), r["partition"]))
    d_big = as_int(biggest["described"], 0)
    lsc = ds.get("largest_single_count") or {}
    R.eq("largest_single_count.partition", lsc.get("partition"),
         biggest["partition"])
    R.eq("largest_single_count.described", lsc.get("described"), d_big)
    R.eq("largest_single_count.share_of_denominator_pct",
         lsc.get("share_of_denominator_pct"), pct(d_big, den_full),
         atol=ATOL_PCT)

    el = ds.get("elasticity_demonstration") or {}
    R.eq("elasticity perturbed_partition", el.get("perturbed_partition"),
         biggest["partition"])
    R.eq("elasticity fraction", el.get("fraction"), ELASTICITY_FRACTION,
         atol=1e-12)
    R.bump()
    if "NOT A RICHNESS CLAIM" not in str(el.get("warning", "")).upper():
        R.fail("the elasticity demonstration does not label itself as a "
               "demonstration rather than a claim about that group's richness")
    for sign, tag in ((-1.0, "minus"), (1.0, "plus")):
        key = "%s_%d_pct" % (tag, int(round(100 * ELASTICITY_FRACTION)))
        got = (el.get("results") or {}).get(key) or {}
        den_alt = den_full + sign * ELASTICITY_FRACTION * d_big
        R.eq("elasticity %s denominator" % key, got.get("denominator"),
             int(round(den_alt)))
        cmp_rungs("elasticity %s" % key, got.get("rungs"),
                  rungs_of(ctx, rows, denominator=den_alt),
                  keys=("numerator", "percent", "numerator_union",
                        "percent_union"), keymap=SENS_KEYMAP)


def check_figure_inputs(ctx, rows):
    """The figure inputs must be projections of the same basis.  A figure that
    can disagree with the headline is a figure that will."""
    euk = {r["partition"]: r for r in rows if is_euk(r)}

    def c0(r):
        return pct(as_int(r["c0_units"], 0), as_int(r["described"], 0))

    def c1(r):
        return pct(as_int(r["c0_units"], 0) * weight(r, "C1"),
                   as_int(r["described"], 0))

    # coverage_per_partition.csv and fig2_coverage_scatter.csv: represented
    # eukaryotic partitions only.
    want_rep = sorted(p for p, r in euk.items()
                      if as_int(r["c0_units"], 0) > 0)
    for path, label in ((COVERAGE_PP, "coverage_per_partition.csv"),
                        (FIG2_SCATTER, "fig2_coverage_scatter.csv")):
        got = read_csv(path, required=False)
        R.eq("%s partition set" % label, sorted(r["partition"] for r in got),
             want_rep)
        for r in got:
            src = euk.get(r["partition"])
            if src is None:
                continue
            R.eq("%s %s coverage_C0_pct" % (label, r["partition"]),
                 (r.get("coverage_C0_pct") or "").strip(),
                 fnum(c0(src), PCT_DP))
            R.eq("%s %s coverage_C1_pct" % (label, r["partition"]),
                 (r.get("coverage_C1_pct") or "").strip(),
                 fnum(c1(src), PCT_DP))
            R.eq("%s %s described" % (label, r["partition"]),
                 as_int(r.get("described")), as_int(src["described"], 0))
            R.eq("%s %s c0_units" % (label, r["partition"]),
                 as_int(r.get("c0_units")),
                 as_int(src["c0_units"], 0))
            R.eq("%s %s f_mol" % (label, r["partition"]),
                 as_float(r.get("f_mol")), src["f_mol"], atol=ATOL_FRAC)
            R.eq("%s %s f_basis" % (label, r["partition"]),
                 (r.get("f_basis") or "").strip(), src["f_basis"])

    # THE DATING FIX: dating is a FILTER, never a coverage column here.
    hdr = csv_header(COVERAGE_PP)
    R.bump()
    for banned in ("coverage_C1d_pct", "coverage_C1b_pct", "C1d", "C1b"):
        if banned in hdr:
            R.fail("coverage_per_partition.csv has a %r column. Dating is a "
                   "FILTER, not a coverage column in the per-partition table: "
                   "a completely sampled but uncalibrated partition must read "
                   "'not ultrametric (CV ...)', never '0'" % banned)
    for need in ("is_chronogram", "roottotip_cv", "interval_kind"):
        R.bump()
        if need not in hdr:
            R.fail("coverage_per_partition.csv is missing the %r column; the "
                   "reader cannot set their own tolerance without it" % need)

    # fig3_c0c1_bars.csv: ALL eukaryotic partitions, structural zeros shown.
    got = read_csv(FIG3_BARS, required=False)
    R.eq("fig3_c0c1_bars.csv partition set",
         sorted(r["partition"] for r in got), sorted(euk))
    n_zero = sum(1 for r in got if tri(r.get("is_structural_zero")) == "yes")
    R.eq("fig3_c0c1_bars.csv structural zeros shown", n_zero,
         sum(1 for r in euk.values() if as_int(r["c0_units"], 0) == 0))

    # fig1_coverage_bars.csv: the ALL-LIFE inventory view.
    got = read_csv(FIG1_BARS, required=False)
    R.eq("fig1_coverage_bars.csv partition set",
         sorted(r["group"] for r in got), sorted(ctx.species))
    for r in got:
        part = r["group"]
        R.eq("fig1 %s described" % part, as_int(r.get("described")),
             ctx.described(part))
        if ctx.is_unrepresented(part):
            continue
        tips = ctx.labels_unique(part)
        chrono = ctx.is_chronogram(part) if ctx.canon_file[part] in ctx.measured \
            else None
        if chrono is None:
            continue
        R.eq("fig1 %s tips" % part, as_int(r.get("tips")), tips)
        R.eq("fig1 %s in_dated_tree" % part, as_int(r.get("in_dated_tree")),
             tips if chrono else 0)
        R.eq("fig1 %s in_undated_tree" % part, as_int(r.get("in_undated_tree")),
             0 if chrono else tips)
        R.eq("fig1 %s dated" % part, tri(r.get("dated")),
             "yes" if chrono else "no")

    # fig3_dated_undated_counts.csv -- counts of canonical trees by category
    # and by the SHIPPED FILE's ultrametricity verdict.
    if any(ctx.canon_file[p] not in ctx.measured for p in ctx.canon_prov):
        R.note("some canonical tree was not measured, so "
               "fig3_dated_undated_counts.csv was not cross-checked")
        return
    got = read_csv(FIG3_COUNTS, required=False)
    tally = defaultdict(lambda: {"Dated": 0, "Undated": 0})
    for part in ctx.canon_prov:
        cat = (ctx.species[part].get("category") or "").strip()
        tally[cat]["Dated" if ctx.is_chronogram(part) else "Undated"] += 1
    for r in got:
        R.eq("fig3 counts %s/%s" % (r.get("category"), r.get("type")),
             as_int(r.get("count")),
             tally.get(r.get("category"), {}).get(r.get("type")))
    R.eq("fig3 counts total", sum(as_int(r.get("count"), 0) for r in got),
         len(ctx.canon_prov))


# ==========================================================================
# CHECK: PERCENTS -- every percentage must be reachable from a pair
# ==========================================================================

PCT_KEY_RE = re.compile(r"(_pct$|^percent(_union)?$)")
NUM_KEYS = ("numerator", "numerator_sum", "numerator_sum_rounded",
            "numerator_union", "shared_labels", "tip_equivalents", "C0_sum",
            "C1_sum", "C1d_sum", "C1b_sum", "c0_units", "described",
            "no_tip_anywhere", "no_molecular_placement", "n_partitions",
            "n_labels_in_two_or_more_partitions", "of_eukaryotic_partitions",
            # EXACT weighted-sum numerators.  A rung numerator is fractional;
            # the rounded integer alone cannot reproduce its own percentage,
            # so the exact value is published beside it and is what makes the
            # percentage reachable.
            "C1_sum_exact", "C1d_sum_exact", "C1b_sum_exact",
            "tip_equivalents_exact", "no_molecular_placement_exact",
            "C1_tip_equivalents",
            "contribution_to_C0", "contribution_to_C1",
            "contribution_to_C1d", "contribution_to_C1b")

# Rule 1 tests a percentage against the EXACT numerator where one is published,
# because that is the value the percentage was formed from.  The rounded
# integer is checked separately, against the exact value it displays.
RULE1_NUM_PREFERENCE = ("numerator_sum", "numerator", "numerator_sum_rounded",
                        "tip_equivalents_exact", "tip_equivalents",
                        "C0_sum", "C1_sum_exact", "C1_sum",
                        "c0_units", "described")
DEN_KEYS = ("denominator", "described", "estimated_total", "estimated_high",
            "numerator_sum", "numerator_sum_rounded", "numerator_union",
            "of_eukaryotic_partitions", "n_partitions",
            # rung totals published as the base of a share_of_* percentage
            "C0_total", "C1_total", "C1d_total", "C1b_total")


def _walk_json(node, path):
    """(path, key, value) for every scalar in a JSON document."""
    if isinstance(node, dict):
        for k, v in node.items():
            if isinstance(v, (dict, list)):
                for t in _walk_json(v, "%s.%s" % (path, k)):
                    yield t
            else:
                yield path, str(k), v
    elif isinstance(node, list):
        for i, v in enumerate(node):
            if isinstance(v, (dict, list)):
                for t in _walk_json(v, "%s[%d]" % (path, i)):
                    yield t
            else:
                yield path, "[%d]" % i, v


def check_percents(ctx):
    """Account for EVERY percentage in headline_numbers.json.

    The retired 23.7 % figure was reachable from no numerator/denominator pair
    anywhere in the deposit -- its implied denominator, 2,084,409, exists in no
    file.  Two rules follow.

    RULE 1 (fatal).  Where a percentage sits in the same object as a numerator
    AND a denominator, it must be exactly 100*n/d.  A percentage that
    contradicts the pair printed beside it is worse than one with no pair.

    RULE 2 (reported).  Every percentage must be REACHABLE: some numerator and
    some denominator that the file actually carries must reproduce it.  A
    percentage that no pair in the whole document can produce is exactly what
    23.7 % was, and it is named here by its full JSON path."""
    nums, dens = set(), set()
    for _path, k, v in _walk_json(ctx.head, "$"):
        f = as_float(v)
        if f is None or isinstance(v, bool):
            continue
        if k in NUM_KEYS:
            nums.add(f)
        if k in DEN_KEYS and f > 0:
            dens.add(f)
    pool = sorted(dens)

    # RULE 1 -- the pair printed beside it
    def rule1(node, path):
        if isinstance(node, dict):
            n = d = None
            # A rung numerator is a WEIGHTED sum and is fractional; the integer
            # beside it is a display rounding.  The percentage is defined on the
            # EXACT value, so that is what rule 1 tests against -- and the
            # rounded value is separately asserted to be its correct rounding,
            # which is a stricter test than either alone.
            for k in RULE1_NUM_PREFERENCE:
                if as_float(node.get(k)) is not None:
                    n = as_float(node[k])
                    break
            if as_float(node.get("denominator")) not in (None, 0):
                d = as_float(node["denominator"])
            exact = as_float(node.get("numerator_sum"))
            shown = as_float(node.get("numerator"))
            if exact is not None and shown is not None:
                R.bump()
                if int(round(exact)) != int(round(shown)):
                    R.fail("%s: numerator=%r is not the rounding of "
                           "numerator_sum=%r" % (path, shown, exact))
            for k, v in node.items():
                if isinstance(v, (dict, list)):
                    rule1(v, "%s.%s" % (path, k))
                    continue
                if not PCT_KEY_RE.search(str(k)) or as_float(v) is None:
                    continue
                if k == "percent" and n is not None and d:
                    R.bump()
                    if not close(as_float(v), pct(n, d), ATOL_PCT):
                        R.fail("%s.%s = %r is not 100*%r/%r, the pair printed "
                               "beside it" % (path, k, v, n, d))
                if k == "percent_union" and d \
                        and as_float(node.get("numerator_union")) is not None:
                    R.bump()
                    if not close(as_float(v),
                                 pct(as_float(node["numerator_union"]), d),
                                 ATOL_PCT):
                        R.fail("%s.%s = %r is not 100*numerator_union/"
                               "denominator" % (path, k, v))
        elif isinstance(node, list):
            for i, v in enumerate(node):
                rule1(v, "%s[%d]" % (path, i))

    rule1(ctx.head, "$")

    # RULE 2 -- reachability from SOME pair the file carries
    orphans, total = [], 0
    for path, k, v in _walk_json(ctx.head, "$"):
        if not PCT_KEY_RE.search(k):
            continue
        f = as_float(v)
        if f is None:
            continue
        total += 1
        if f == 0.0:
            continue
        ok = False
        for d in pool:
            n = f * d / 100.0
            for cand in nums:
                if abs(cand - n) <= max(ATOL_SUM, abs(n) * 1e-9):
                    ok = True
                    break
            if ok:
                break
        if not ok:
            orphans.append("%s.%s = %r" % (path, k, v))
    R.bump(total)
    R.bump()
    if total == 0:
        R.fail("headline_numbers.json contains no percentage at all; the "
               "headline file is not being generated")
    if orphans:
        R.note("%d of %d percentages in headline_numbers.json are not "
               "reproducible from any numerator/denominator the file itself "
               "carries. Each must be wired to its pair or deleted -- an "
               "unreachable percentage is precisely what 23.7 %% was."
               % (len(orphans), total), orphans[:25])
    else:
        R.note("all %d percentages in headline_numbers.json are reachable "
               "from a numerator and a denominator the file carries." % total)


# ==========================================================================
# CHECK: MONOTONICITY
# ==========================================================================

def check_monotonicity(ctx):
    rows = ctx.my_rows
    euk = [r for r in rows if is_euk(r)]

    # 1. per partition: nested conjuncts, so weights are non-increasing.
    for r in euk:
        ws = [weight(r, k) for k in RUNGS]
        for i in range(len(ws) - 1):
            R.bump()
            if ws[i + 1] > ws[i] + 1e-12:
                R.fail("%s: w(%s)=%r exceeds w(%s)=%r -- the rungs are not "
                       "nested" % (r["partition"], RUNGS[i + 1], ws[i + 1],
                                   RUNGS[i], ws[i]))
        n = as_int(r["c0_units"], 0)
        cs = [n * w for w in ws]
        for i in range(len(cs) - 1):
            R.bump()
            if cs[i + 1] > cs[i] + 1e-9:
                R.fail("%s: contribution to %s exceeds its contribution to %s"
                       % (r["partition"], RUNGS[i + 1], RUNGS[i]))

    # 2. totals, on the published file as well as on the measurement.
    for label, src in (("measured", rungs_of(ctx, rows)),
                       ("headline_numbers.json", ctx.head.get("rungs") or {})):
        for i in range(len(RUNGS) - 1):
            a = as_float((src.get(RUNGS[i]) or {}).get("numerator_sum"))
            b = as_float((src.get(RUNGS[i + 1]) or {}).get("numerator_sum"))
            R.bump()
            if a is None or b is None:
                R.fail("%s: rung %s or %s has no numerator_sum"
                       % (label, RUNGS[i], RUNGS[i + 1]))
            elif b > a + ATOL_SUM:
                R.fail("%s: %s numerator %r exceeds %s numerator %r"
                       % (label, RUNGS[i + 1], b, RUNGS[i], a))

    # 3. the per-partition table: C1 can never exceed C0.
    for r in read_csv(COVERAGE_PP, required=False):
        a, b = as_float(r.get("coverage_C0_pct")), as_float(r.get("coverage_C1_pct"))
        R.bump()
        if a is not None and b is not None and b > a + 1e-9:
            R.fail("coverage_per_partition.csv %s: C1 %.4f exceeds C0 %.4f"
                   % (r["partition"], b, a))

    # 4. every restricted sensitivity is a restriction: numerators may only go
    #    down when partitions are removed from the numerator side.
    full = rungs_of(ctx, rows)
    la = (ctx.head.get("sensitivities") or {}).get("leave_assembly_class_out") or {}
    for key, entry in la.items():
        if not isinstance(entry, dict) or "rungs" not in entry:
            continue
        for rung in RUNGS:
            got = as_float(((entry["rungs"] or {}).get(rung) or {}).get("numerator"))
            R.bump()
            if got is not None and got > full[rung]["numerator_sum_rounded"] + 1:
                R.fail("sensitivity %s raises the %s numerator to %r above the "
                       "unrestricted %r; dropping partitions from a numerator "
                       "cannot increase it"
                       % (key, rung, got, full[rung]["numerator_sum_rounded"]))

    # 5. bands: the low end can never exceed the high end.
    bands = (ctx.head.get("sensitivities") or {}).get("f_mol_bands") or {}
    for part, e in bands.items():
        for rung in RUNGS:
            lo = as_float(((e.get("at_f_low") or {}).get(rung) or {}).get("percent"))
            hi = as_float(((e.get("at_f_high") or {}).get(rung) or {}).get("percent"))
            R.bump()
            if lo is not None and hi is not None and lo > hi + ATOL_PCT:
                R.fail("band %s: %s at f_low (%r) exceeds %s at f_high (%r)"
                       % (part, rung, lo, rung, hi))

    # 6. denominator monotonicity: a larger denominator cannot raise a percent.
    ds = (ctx.head.get("sensitivities") or {}).get("denominator") or {}
    series = []
    for key in ("represented_partitions_only", "full_denominator"):
        e = ds.get(key) or {}
        if e.get("described") and e.get("rungs"):
            series.append((as_float(e["described"]), key, e["rungs"]))
    for key, e in ((k, v) for k, v in
                   ((ds.get("elasticity_demonstration") or {}).get("results")
                    or {}).items()):
        if e.get("denominator") and e.get("rungs"):
            series.append((as_float(e["denominator"]), key, e["rungs"]))
    series.sort(key=lambda t: t[0])
    for i in range(len(series) - 1):
        (d1, k1, r1), (d2, k2, r2) = series[i], series[i + 1]
        for rung in RUNGS:
            p1 = as_float((r1.get(rung) or {}).get("percent"))
            p2 = as_float((r2.get(rung) or {}).get("percent"))
            R.bump()
            if p1 is not None and p2 is not None and p2 > p1 + ATOL_PCT:
                R.fail("denominator %s (%s) gives %s = %r, but the larger "
                       "denominator %s (%s) gives %r -- a larger denominator "
                       "cannot raise a coverage percentage"
                       % (commas(d1), k1, rung, p1, commas(d2), k2, p2))

    # 7. the chronogram tolerance: relaxing it can only admit more trees.
    counts = []
    for tol in (0.001, 0.01, CHRONOGRAM_TOLERANCE, 0.10, 0.50):
        counts.append((tol, sum(1 for p in ctx.canon_prov
                                if ctx.canon_file[p] in ctx.measured
                                and (ctx.measured[ctx.canon_file[p]]["cv"] or 9e9)
                                <= tol)))
    for i in range(len(counts) - 1):
        R.bump()
        if counts[i + 1][1] < counts[i][1]:
            R.fail("raising the ultrametricity tolerance from %r to %r lowers "
                   "the chronogram count from %d to %d"
                   % (counts[i][0], counts[i + 1][0], counts[i][1],
                      counts[i + 1][1]))
    R.note("canonical chronogram count by ultrametricity tolerance: %s. "
           "roottotip_cv ships alongside the bit precisely so a reader can "
           "choose their own row of this table."
           % ", ".join("CV<=%g -> %d" % c for c in counts))


# ==========================================================================
# CHECK: MOLTIPS
# ==========================================================================

def check_c1count(ctx):
    """Audit A01.  C1 counts distinct accepted names among the molecular tips
    (basis c1_units), re-measured here; the superseded product c0_units * f_mol
    is published beside it as a sensitivity, and the gap is always printed."""
    mine = {r["partition"]: r for r in ctx.my_rows}
    for part, r in ctx.basis.items():
        m = mine.get(part)
        if m is None or ctx.is_unrepresented(part) or part not in ctx.derived:
            continue
        R.eq("%s basis c1_units" % part, as_float(r.get("c1_units")),
             m["c1_units"], atol=1e-4)
        R.eq("%s basis c1_basis" % part, (r.get("c1_basis") or "").strip(),
             m["c1_basis"])
        R.bump()
        if m["c1_units"] is not None and m["c1_units"] > as_int(m["c0_units"], 0):
            R.fail("%s: C1 (%s) exceeds C0 (%s); C1 is a subset of C0"
                   % (part, m["c1_units"], m["c0_units"]))
    sens = (ctx.head.get("sensitivities") or {}).get("c1_definition")
    R.bump()
    if not sens:
        R.fail("headline_numbers.json carries no sensitivities.c1_definition: "
               "the superseded C1 product must be published beside the count")
        return
    alt = [dict(r) for r in ctx.my_rows]
    for r in alt:
        r["c1_units"] = None
    prod = rungs_of(ctx, alt)
    count = rungs_of(ctx, ctx.my_rows)
    cmp_rungs("c1_definition.product", sens.get("product"), prod,
              keys=("numerator", "percent"), keymap=SENS_KEYMAP)
    diff = sens.get("difference_pp") or {}
    for k in ("C1", "C1d", "C1b"):
        # published to 4 decimals (PCT_DP of the rebuild)
        R.eq("c1_definition.difference_pp.%s" % k, diff.get(k),
             round(count[k]["percent"] - prod[k]["percent"], 4), atol=1e-9)
    per = sens.get("per_partition") or {}
    for r in ctx.my_rows:
        if not is_euk(r) or as_int(r.get("c0_units"), 0) <= 0:
            continue
        g = per.get(r["partition"]) or {}
        R.eq("c1_definition %s c1_count" % r["partition"], g.get("c1_count"),
             r["c1_units"], atol=1e-4)
        R.eq("c1_definition %s c1_product" % r["partition"],
             g.get("c1_product"),
             as_int(r["c0_units"], 0) * (r["f_mol"] or 0.0), atol=1e-3)
    R.note("C1 is counted as distinct accepted names among molecular tips: "
           "C1 %.4f %%, C1d %.4f %%, C1b %.4f %%. The superseded product "
           "c0_units * f_mol gives C1 %.4f %%, C1d %.4f %%, C1b %.4f %% "
           "(difference %+.4f / %+.4f / %+.4f pp)."
           % (count["C1"]["percent"], count["C1d"]["percent"],
              count["C1b"]["percent"], prod["C1"]["percent"],
              prod["C1d"]["percent"], prod["C1b"]["percent"],
              count["C1"]["percent"] - prod["C1"]["percent"],
              count["C1d"]["percent"] - prod["C1d"]["percent"],
              count["C1b"]["percent"] - prod["C1b"]["percent"]))
    check_c1b_sensitivity(ctx)


GBIF_UNMATCHED = P("Data", "external", "gbif_unmatched.csv")


def check_gbif(ctx):
    """A10.  C0 with the COL26.9-unmatched tips that the GBIF backbone
    resolves, recomputed from Data/external/gbif_unmatched.csv: a strict
    EXACT or FUZZY match at rank SPECIES with an accepted, synonym or doubtful
    status adds its accepted species when the partition's C0 has no CoL name
    with the same 'Genus epithet'."""
    sens = (ctx.head.get("sensitivities") or {}).get("c0_gbif_unmatched")
    R.bump()
    if not sens:
        R.fail("headline_numbers.json carries no "
               "sensitivities.c0_gbif_unmatched (A10)")
        return
    gb = {(r["query"], r["kingdom"]): r for r in read_csv(GBIF_UNMATCHED)}
    ok_status = {"ACCEPTED", "SYNONYM", "HETEROTYPIC_SYNONYM",
                 "HOMOTYPIC_SYNONYM", "PROPARTE_SYNONYM", "DOUBTFUL"}
    added, n_q, n_res, c0 = 0, 0, 0, 0
    per = sens.get("per_partition") or {}
    for r in ctx.my_rows:
        if not is_euk(r) or as_int(r["c0_units"], 0) <= 0:
            continue
        p = r["partition"]
        c0 += as_int(r["c0_units"], 0)
        rows = ctx.col_rows.get(p, {})
        kc = Counter(x.get("col_kingdom") for x in rows.values()
                     if x.get("col_kingdom"))
        king = kc.most_common(1)[0][0] if kc else ""
        have = {" ".join((x.get("col_name") or "").split()[:2])
                for x in rows.values()
                if x.get("match_type") in COL_RESOLUTION_RESOLVED
                and x.get("col_id")}
        add = set()
        for x in rows.values():
            if x.get("match_type") != "unmatched" or not x.get("query"):
                continue
            n_q += 1
            g = gb.get((x["query"], king))
            R.bump()
            if g is None:
                R.fail("%s: unmatched query %r has no row in %s"
                       % (p, x["query"], rel(GBIF_UNMATCHED)))
                continue
            if g["matchType"] in ("EXACT", "FUZZY") and g["rank"] == "SPECIES" \
                    and g["status"] in ok_status and g["species"]:
                n_res += 1
                if g["species"] not in have:
                    add.add(g["species"])
        added += len(add)
        if p in per or add:
            R.eq("c0_gbif %s species added" % p,
                 (per.get(p) or {}).get("gbif_species_added"), len(add))
    den = sum(as_int(r["described"], 0) for r in ctx.my_rows if is_euk(r))
    R.eq("c0_gbif unmatched tips", sens.get("unmatched_tips"), n_q)
    R.eq("c0_gbif resolved by GBIF", sens.get("unmatched_tips_resolved_by_gbif"),
         n_res)
    R.eq("c0_gbif species added", sens.get("species_added"), added)
    R.eq("c0_gbif numerator", sens.get("numerator"), c0 + added)
    R.eq("c0_gbif percent", sens.get("percent"), pct(c0 + added, den),
         atol=ATOL_PCT)
    R.note("C0 under COL26.9 alone %.4f %%; with the %d tips that only the "
           "GBIF backbone resolves (%d species added) %.4f %%"
           % (pct(c0, den), n_res, added, pct(c0 + added, den)))


def check_c1b_sensitivity(ctx):
    # D22 sensitivity: C1b as if the calibration-fixed intervals counted
    sens = (ctx.head.get("sensitivities") or {}).get("c1b_calibration_fixed")
    R.bump()
    if not sens:
        R.fail("headline_numbers.json carries no "
               "sensitivities.c1b_calibration_fixed (D22)")
        return
    added = sorted(p for p, o in ctx.crown.items()
                   if o.get("c1b_exclusion") in ("degenerate", "at_bound",
                                                 "prior_reproduced")
                   or o.get("kind") == "calibration_prior_sample")
    R.eq("c1b_calibration_fixed.partitions_added",
         sens.get("partitions_added"), added)
    alt = [dict(r) for r in ctx.my_rows]
    for r in alt:
        if r["partition"] in added and r["c0_units"]:
            r["interval_recoverable"] = "yes"
    cmp_rungs("c1b_calibration_fixed", sens.get("rungs"), rungs_of(ctx, alt),
              keys=("numerator", "percent"), keymap=SENS_KEYMAP)
    # R2-06: every archived interval counted, whatever its kind
    sa = (ctx.head.get("sensitivities") or {}).get("c1b_all_intervals")
    R.bump()
    if not sa:
        R.fail("headline_numbers.json carries no sensitivities.c1b_all_intervals")
        return
    added = sorted(p for p, o in ctx.crown.items()
                   if o.get("counts_for_c1b") != "yes")
    R.eq("c1b_all_intervals.partitions_added", sa.get("partitions_added"), added)
    alt = [dict(r) for r in ctx.my_rows]
    for r in alt:
        if r["partition"] in added and r["c0_units"]:
            r["interval_recoverable"] = "yes"
    cmp_rungs("c1b_all_intervals", sa.get("rungs"), rungs_of(ctx, alt),
              keys=("numerator", "percent"), keymap=SENS_KEYMAP)


def check_moltips(ctx):
    R.eq("molecular_tips.csv columns", csv_header(MOLTIPS),
         ["partition", "standardized_id"])
    rows = read_csv(MOLTIPS, required=False)
    by_part = defaultdict(list)
    for r in rows:
        by_part[r.get("partition")].append(r.get("standardized_id"))

    want = sorted(p for p, d in ctx.derived.items() if d["source_file"])
    R.eq("molecular_tips.csv partition set", sorted(by_part.keys()), want)

    for part in want:
        d = ctx.derived[part]
        ids = by_part.get(part, [])
        R.eq("molecular_tips %s row count" % part, len(ids),
             d["measured_intersection"])
        R.bump()
        if len(set(ids)) != len(ids):
            R.fail("molecular_tips.csv: %s has %d duplicate (partition, id) "
                   "rows" % (part, len(ids) - len(set(ids))))
        if d["f_basis"] == "per_tip_evidence":
            names = d["evidence_names"] or set()
            bad = [i for i in ids if ctx.by_id.get(i, i) not in names]
            R.bump()
            if bad:
                R.fail("molecular_tips.csv: %d id(s) for %s carry no per-tip "
                       "evidence (e.g. %s)" % (len(bad), part, bad[:5]))
            continue
        # every id must resolve, through the deposit's own mapping, to a name
        # that is genuinely in the intersection.
        src_path = P(d["source_file"]) if d["source_file"] else None
        if not src_path or not os.path.exists(src_path):
            R.bump()
            R.fail("the molecular source recorded for %s (%r) does not exist; "
                   "C1 for that partition rests on nothing" % (part, src_path))
            continue
        # D4: a shipped tip is molecular when its species key is in the source.
        src_keys = d4_keys(ctx.source_sets(src_path, ctx.group(part))[0])
        inter = {nm for nm in ctx.labels.get(part, set())
                 if d4_keys([nm]) & src_keys}
        bad = [i for i in ids if ctx.by_id.get(i, i) not in inter]
        R.bump()
        if bad:
            R.fail("molecular_tips.csv: %d id(s) for %s do not resolve into "
                   "the measured molecular tip set (e.g. %s)"
                   % (len(bad), part, bad[:5]))
        idmap = ctx.idmap.get(part, {})
        missing = [nm for nm in inter if nm not in idmap]
        R.bump()
        if missing:
            R.fail("%d name(s) of %s's molecular tip set have no id in the "
                   "shipped tree (e.g. %s)" % (len(missing), part, missing[:5]))

    unenumerated = sorted(p for p in ctx.canon_euk
                          if p in ctx.derived and not ctx.derived[p]["source_file"])
    if unenumerated:
        R.note("%d of %d eukaryotic canonical partitions are NOT per-tip "
               "enumerated in molecular_tips.csv; their C1 equals C0 on an "
               "f_mol = 1 that is asserted, not measured per tip. Until every "
               "partition is enumerated the paper writes TIP-EQUIVALENTS, "
               "never species, at C1 and below: %s"
               % (len(unenumerated), len(ctx.canon_euk),
                  ", ".join(unenumerated)))


# ==========================================================================
# CHECK: UNITS
# ==========================================================================

def check_units(ctx):
    rungs = ctx.head.get("rungs") or {}
    R.eq("rungs.C0.unit", (rungs.get("C0") or {}).get("unit"), "species")
    for k in ("C1", "C1d", "C1b"):
        R.eq("rungs.%s.unit" % k, (rungs.get(k) or {}).get("unit"),
             "tip-equivalents")
    R.bump()
    if "tip-equivalent" not in str(ctx.head.get("units", "")).lower():
        R.fail("headline_numbers.json `units` does not state the "
               "tip-equivalents rule for C1 and below")
    R.bump()
    mono = str(ctx.head.get("monotonicity", "")).lower()
    if "nested tip sets are not claimed" not in mono \
            and "not claimed" not in mono:
        R.fail("the monotonicity statement claims more than the computed "
               "object supports: the rungs are nested WEIGHTED SUMS, and "
               "nested tip SETS must not be claimed until every partition is "
               "enumerated")
    R.bump()
    nd = str(ctx.head.get("naming_discipline", "")).lower()
    if "placement traceability" not in nd:
        R.fail("the naming discipline for C1 is missing or incomplete: C1 is "
               "MOLECULARLY SUPPORTED PLACEMENT and must be stated NOT to be "
               "a predicate on placement traceability -- f_mol is an axis-(a) "
               "quantity, and claiming the two-axis conflict was dissolved is "
               "an overclaim a referee will catch")
    # "dissolved" must not be ASSERTED -- but the deposit is expected to warn
    # against the overclaim in so many words, and a bare substring test cannot
    # tell the claim from the prohibition.  Scan sentence by sentence: a
    # sentence that mentions dissolution must also carry a prohibition marker.
    R.bump()
    _PROHIBIT = ("overclaim", "must not", "never", "is not", "not a",
                 "do not", "cannot", "would be wrong", "refute", "rather than")
    claimed = [s.strip() for s in re.split(r"[.;]", nd)
               if "dissolve" in s and not any(p in s for p in _PROHIBIT)]
    if claimed:
        R.fail("the naming discipline ASSERTS the two-axis conflict was "
               "'dissolved'. It was not: the genuine axis-(b) distinction is "
               "reported as the NOMINAL topology_assembly annotation because "
               "the deposit carries no measured per-tip placement-warrant flag",
               claimed[:3])


# ==========================================================================
# CHECK: CROSSLEDGER -- no field may hold two different values
# ==========================================================================

def note_documents(note, value):
    """Does the ground truth's own conflict note record this competing value?

    A disagreement the ground truth NAMES is reconciliation; a disagreement it
    does not name is rot."""
    if value is None:
        return False
    blob = re.sub(r"(?<=\d),(?=\d{3}\b)", "", str(note or ""))
    return re.search(r"(?<![\d.])%d(?![\d.])" % int(value), blob) is not None


def check_crossledger(ctx):
    # ---- described species, per partition -----------------------------
    ps = {r.get("partition"): r for r in read_csv(PARTITION_SUMMARY,
                                                  required=False)}
    de = {r.get("group"): r for r in read_csv(DATA_ESTIMATES, required=False)}
    f1 = {r.get("group"): r for r in read_csv(FIG1_BARS, required=False)}

    unjoined = sorted(k for k in de if k and k not in ctx.species)
    if unjoined:
        R.note("%d row(s) of Data/data_estimates.csv carry a `group` that is "
               "no partition of truth_species.csv and were not cross-checked: "
               "%s" % (len(unjoined), ", ".join(unjoined[:12])))

    for part, sp in ctx.species.items():
        truth_v = ctx.described(part)
        note = sp.get("conflict_note") or ""
        claims = [
            ("Data/coverage_basis.csv",
             as_int((ctx.basis.get(part) or {}).get("described"))),
            ("Data/site/partition_summary.csv",
             as_int((ps.get(part) or {}).get("described_species"))),
            ("Paper/inputs/fig1_coverage_bars.csv",
             as_int((f1.get(part) or {}).get("described"))),
            ("Data/data_estimates.csv",
             as_int((de.get(part) or {}).get("described"))),
        ]
        for src, v in claims:
            if v is None:
                continue
            R.bump()
            if v == truth_v:
                continue
            if note_documents(note, v):
                R.note("%s: %s carries described=%s against ground truth %s; "
                       "the disagreement is DOCUMENTED in "
                       "truth_species.csv's conflict_note"
                       % (part, src, commas(v), commas(truth_v)))
            else:
                R.fail("%s: %s says described=%s, truth_species.csv says %s, "
                       "and nothing in the ground truth documents the "
                       "difference. One field, two values."
                       % (part, src, commas(v), commas(truth_v)))

        for field, col_ps, col_de in (("estimated_low", "estimated_low",
                                       "estimated_low"),
                                      ("estimated_high", "estimated_high",
                                       "estimated_high")):
            tv = as_int(sp.get(field))
            for src, v in (("Data/site/partition_summary.csv",
                            as_int((ps.get(part) or {}).get(col_ps))),
                           ("Data/data_estimates.csv",
                            as_int((de.get(part) or {}).get(col_de)))):
                if v is None or tv is None:
                    continue
                R.bump()
                if v != tv and not note_documents(note, v):
                    R.fail("%s: %s says %s=%s, truth_species.csv says %s"
                           % (part, src, field, commas(v), commas(tv)))

    # ---- tip counts, per shipped tree -----------------------------------
    md = {r.get("filename"): r for r in read_csv(STD_METADATA, required=False)}
    tp = {r.get("filename"): r for r in read_csv(TREE_PROVENANCE,
                                                 required=False)}
    for fname, t in ctx.tips.items():
        raw = as_int(t.get("raw_leaves"))
        uniq = as_int(t.get("distinct_labels"))
        admissible = {v for v in (raw, uniq) if v is not None}
        # A pre-audit ledger may legitimately carry a STALE tip count -- that
        # the pre-audit ledgers disagree with the files is the paper's finding,
        # and Data/tree_provenance.csv is retained deliberately as the
        # pre-audit record.  It is downgraded to a note ONLY when the ground
        # truth ANTICIPATED that exact value: truth_tips.csv must carry it as
        # `published_ntips` AND explain it with a `published=<reason>` clause
        # in `convention`.  Any other mismatch is still a hard failure, and the
        # convention check above independently fails a discrepancy with no
        # reason -- so this cannot become a blanket amnesty.
        pub_n = as_int(t.get("published_ntips"))
        conv_raw = (t.get("convention") or "")
        pub_reason = None
        for clause in conv_raw.split(";"):
            clause = clause.strip()
            if clause.startswith("published="):
                pub_reason = clause.split("=", 1)[1].strip()
        documented = (pub_n is not None and pub_n not in admissible
                      and pub_reason)
        src_n = ctx.source_leaves.get(fname)
        for src, v in (("Data/standardized/metadata.csv",
                        as_int((md.get(fname) or {}).get("ntips"))),
                       ("Data/tree_provenance.csv",
                        as_int((tp.get(fname) or {}).get("ntips_effective"))),
                       ("Data/tree_provenance.csv:ntips_raw",
                        as_int((tp.get(fname) or {}).get("ntips_raw")))):
            if v is None:
                continue
            R.bump()
            if v in admissible:
                continue
            # Data/tree_provenance.csv:ntips_raw is DEFINED as the leaf count
            # of the source tree before D7 pruning and D4 collapsing.  It is
            # accepted only in that column, and only when it equals the leaf
            # count this script measured from the source file itself.
            if (src == "Data/tree_provenance.csv:ntips_raw"
                    and src_n is not None and v == src_n):
                R.note("%s: %s carries the SOURCE-TREE leaf count %s, "
                       "measured from %s, before D7 pruning (%d tips, "
                       "Data/pruned_tips.csv) and D4 collapsing; the shipped "
                       "file has %s leaves."
                       % (fname, src, commas(v), ctx.canon_sources.get(fname),
                          ctx.n_pruned.get(fname, 0), commas(raw)))
                continue
            if documented and v == pub_n:
                R.note("%s: %s carries the PRE-AUDIT tip count %s (file: %s "
                       "leaves, %s distinct labels). truth_tips.csv records it "
                       "as published_ntips and explains it as '%s'. This is "
                       "the paper's finding, not a defect -- the ledger is "
                       "retained as the pre-audit record and demoted in "
                       "Methods." % (fname, src, commas(v), commas(raw),
                                     commas(uniq), pub_reason))
                continue
            R.fail("%s: %s says %s tips; the file has %s leaves and %s "
                   "distinct labels. The value matches NEITHER admissible "
                   "convention, and the ground truth does not document it "
                   "(published_ntips=%s, reason=%s)."
                   % (fname, src, commas(v), commas(raw), commas(uniq),
                      commas(pub_n) if pub_n is not None else "-",
                      pub_reason or "-"))

    # ---- the canonical tree of each partition ---------------------------
    for part in ctx.canon_prov:
        want = ctx.canon_file[part]
        for src, v in (("Data/coverage_basis.csv",
                        (ctx.basis.get(part) or {}).get("canonical_tree")),
                       ("Data/site/partition_summary.csv",
                        (ps.get(part) or {}).get("filename"))):
            v = (v or "").strip()
            if not v:
                continue
            R.eq("%s canonical tree per %s" % (part, src), v, want)

    # ---- representation ---------------------------------------------------
    for part in ctx.species:
        want = "yes" if ctx.is_unrepresented(part) else "no"
        for src, v in (("Data/coverage_basis.csv",
                        (ctx.basis.get(part) or {}).get("is_unrepresented")),
                       ("Data/site/partition_summary.csv",
                        (ps.get(part) or {}).get("is_unrepresented")),
                       ("Paper/inputs/fig1_coverage_bars.csv",
                        (f1.get(part) or {}).get("is_unrepresented"))):
            if v is None or str(v).strip() == "":
                continue
            R.eq("%s is_unrepresented per %s" % (part, src), tri(v), want)

    # ---- the denominators quoted in the headline -------------------------
    hd = ctx.head.get("denominators") or {}
    euk_d = sum(ctx.described(p) for p in ctx.eukaryotic)
    R.eq("headline eukaryotic denominator == column sum of the basis file",
         (hd.get("eukaryotic") or {}).get("described"),
         sum(as_int(r.get("described"), 0) for r in ctx.basis_rows
             if tri(r.get("is_prokaryotic")) == "no"))
    R.eq("headline eukaryotic denominator == truth_species.csv column sum",
         (hd.get("eukaryotic") or {}).get("described"), euk_d)


# ==========================================================================
# CHECK: NOLITERALS / NOABSPATHS / RETIRED  (source hygiene)
# ==========================================================================

BENIGN_CONTEXT = (
    "figsize", "dpi", "linewidth", "lw=", "alpha", "fontsize", "markersize",
    "ms=", "pad", "wspace", "hspace", "bbox", "xlim", "ylim", "set_size",
    "subplots_adjust", "rcparams", "zorder", "width=", "height=", "color",
    "labelsize", "titlesize", "elinewidth", "capsize", "rotation", "vmin",
    "vmax", "aspect", "gridspec", "margins", "tick_params", "cm(", "alpha=",
)

def grouped(digits):
    """A regex for an integer however the deposit spells it: 2084409,
    2,084,409, or the supplement's LaTeX house style 2{,}084{,}409."""
    parts = []
    i = len(digits) % 3 or 3
    parts.append(digits[:i])
    while i < len(digits):
        parts.append(digits[i:i + 3])
        i += 3
    sep = r"(?:\s|,|\{,\}|\\,|\.)?"
    return r"(?<![\d.])" + sep.join(parts) + r"(?![\d])"


RETIRED_PATTERNS = [
    (r"(?<![\d.])23\.7(?![\d])",
     "the retired 23.7 % coverage headline -- a hard-coded literal reachable "
     "from no numerator/denominator pair in the deposit"),
    (r"(?<![\d.])8\.65(?![\d])",
     "the retired middle rung 8.65 %, which never had a principled definition"),
    (grouped("2084409"),
     "the phantom denominator 2,084,409, which exists in no file"),
    (r"\b(1[678])\s+(dated|time-calibrated|chronogram|ultrametric)",
     "the ledger claim of 17/16/18 dated canonical trees; the files say what "
     "this script measures"),
]

# Value-plus-context patterns: precise enough not to fire on an unrelated
# number that happens to collide with a retired one.
RETIRED_CONTEXTUAL = [
    (grouped("1633"), r"cnidari|coral|scleractini",
     "the Cnidaria described-species regression 1,633 (a count scoped to "
     "Scleractinia); the adopted figure is the Cnidaria-wide one"),
    (grouped("25000"), r"nematod",
     "the Nematode described-species regression 25,000 (Catalogue of Life); "
     "the adopted figure is Hodda 2022, and the attribution must move with it"),
    (grouped("494005"), r"table\s*1|denominator|numerator|coverage|canonical",
     "494,005 -- the label count across ALL 264 trees, used where Table 1's "
     "own rule (canonical trees only) gives 466,867 union / 466,895 sum"),
]

# In PROSE ONLY, a retired figure inside an explicit retirement disclosure is
# allowed: the paper is expected to volunteer the correction, not hide it.
RETIREMENT_MARKERS = re.compile(
    r"withdraw|retire|superseded|previously|formerly|erratum|no longer|"
    r"replaced|corrigend|earlier version|used to|incorrect|excluded|"
    r"non-canonical|all 264|canonical only|difference|we reported|"
    r"published version|earlier draft", re.I)


def _blank_docstring_lines(text):
    """Blank out DOCSTRING lines only -- bare string-expression statements.

    A string that is merely documentation cannot be USED as a value, so a
    retired figure quoted inside one is prose, not a live constant, and the
    hygiene checks must not accuse it.  Crucially this blanks ONLY bare
    string-expression statements: a string that is assigned to a name, passed
    as an argument, or interpolated is still a value that could be rendered
    onto a figure, and remains fully in scope for every check.  Narrowing the
    check to documentation is not the same as whitelisting a file.
    """
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return text                      # unparseable: check it as-is, verbatim
    lines = text.splitlines()
    blank = set()
    for node in ast.walk(tree):
        body = getattr(node, "body", None)
        if not isinstance(body, list):
            continue
        for stmt in body:
            if (isinstance(stmt, ast.Expr)
                    and isinstance(getattr(stmt, "value", None), ast.Constant)
                    and isinstance(stmt.value.value, str)):
                lo = stmt.lineno
                hi = getattr(stmt, "end_lineno", lo) or lo
                blank.update(range(lo, hi + 1))
    for i in sorted(blank):
        if 1 <= i <= len(lines):
            lines[i - 1] = ""
    return "\n".join(lines)


def strip_python_comments(text):
    text = _blank_docstring_lines(text)
    out = []
    for line in text.splitlines():
        s, i, q = [], 0, None
        while i < len(line):
            ch = line[i]
            if q:
                s.append(ch)
                if ch == "\\":
                    if i + 1 < len(line):
                        s.append(line[i + 1])
                    i += 2
                    continue
                if ch == q:
                    q = None
            elif ch in "\"'":
                q = ch
                s.append(ch)
            elif ch == "#":
                break
            else:
                s.append(ch)
            i += 1
        out.append("".join(s))
    return "\n".join(out)


def strip_tex_comments(text):
    return "\n".join(re.sub(r"(?<!\\)%.*$", "", ln) for ln in text.splitlines())


def strip_comments(path, text):
    low = path.lower()
    if low.endswith((".py", ".r")):
        return strip_python_comments(text)
    if low.endswith((".tex", ".sty")):
        return strip_tex_comments(text)
    return text


def figure_scripts():
    out = []
    for scope in FIGURE_SCOPES:
        out.extend(walk_files(scope, (".py", ".r", ".R")))
    return [p for p in out if os.path.abspath(p) not in TEXT_SCAN_EXEMPT]


def literal_hits(line, forbidden):
    """The literals on one (comment-stripped) source line that reproduce a
    computed coverage quantity: a list of (literal, what it is)."""
    hits = []
    for m in re.finditer(r"(?<![\w.])(\d+(?:\.\d+)?)(?![\w.])", line):
        tok = m.group(1)
        val = as_float(tok)
        if val is None:
            continue
        # The percent scale factor, and only it: a literal 100 that is
        # an operand of a multiplication ("100.0 * n / d").  A reported
        # 100.0% (Archaea, Bacteria C0) written as a value still fails.
        if val == 100.0 and (re.match(r"\s*\*", line[m.end():])
                             or re.search(r"\*\s*$", line[:m.start()])):
            continue
        # PRECISION-AWARE MATCHING.  A literal is compared only
        # against a computed value rounded to the SAME number of
        # decimals the literal itself is written to.  Comparing a
        # 2-decimal literal against a 1-decimal rounding is what made
        # the plot constants 1.30 and 1.10 collide with Sponges C0
        # (1.2566 -> 1.3) and Fungi C0 (1.0645 -> 1.1).  This is
        # strictly tighter arithmetic, not a relaxed rule: a literal
        # written 21.6 still matches the 1-decimal rounding of 21.62.
        if "." not in tok:
            if int(val) < 1000:
                continue
            key = ("int", int(val))
        else:
            if val < 1.0:
                continue
            dp = len(tok.split(".")[1])
            key = (dp, round(val, dp))
        if key in forbidden:
            hits.append((tok, forbidden[key]))
    return hits


def check_noliterals(ctx):
    """No coverage number may be hard-coded in a figure script.

    Two rules.  RULE A is a blacklist of the retired figures -- unconditional,
    because none of them may reappear anywhere.  RULE B is derived from what
    this script measured: any literal in a figure script that reproduces a
    coverage quantity the deposit computes is a literal that will one day
    disagree with it.  Benign geometry (figure sizes, alphas, axis limits) is
    exempted by an explicit, printed context list, and integers below 1000 are
    ignored as too collision-prone to accuse."""
    # the set of numbers a figure script must not contain
    forbidden = {}
    rows = ctx.my_rows
    for r in rows:
        d, n = as_int(r["described"], 0), as_int(r["c0_units"], 0)
        if not d:
            continue
        for rung in RUNGS:
            v = pct(n * weight(r, rung), d)
            if v is None or v < 1.0:
                continue
            for dp in (1, 2, 3, 4):
                forbidden.setdefault((dp, round(v, dp)), "%s %s coverage"
                                     % (r["partition"], rung))
    for rung, e in rungs_of(ctx, rows).items():
        for dp in (1, 2, 3, 4):
            forbidden[(dp, round(e["percent"], dp))] = (
                "the %s headline percentage" % rung)
        forbidden[("int", e["numerator_sum_rounded"])] = "the %s numerator" % rung
        if "numerator_union" in e:
            forbidden[("int", e["numerator_union"])] = ("the %s union numerator"
                                                        % rung)
    forbidden[("int", sum(as_int(r["described"], 0) for r in rows if is_euk(r)))] = \
        "the eukaryotic denominator"

    # SELF-TEST.  A coverage value that rounds to a whole number at one
    # decimal ("86.0", "10.0", "100.0") is still a reported value and must be
    # caught; the percent scale factor "100.0 * n / d" must not be.
    whole = sorted((v, why) for (dp, v), why in
                   ((k, w) for k, w in forbidden.items() if k[0] == 1)
                   if v == int(v))
    R.bump()
    if not whole:
        R.note("NOLITERALS self-test: no coverage value rounds to a whole "
               "number at one decimal, so the whole-number probe has no target")
    else:
        v, why = whole[0]
        probe = "X_PROBE = %.1f" % v
        if not literal_hits(probe, forbidden):
            R.fail("NOLITERALS self-test: the probe %r (%s) was NOT caught; a "
                   "whole-number-valued coverage literal escapes the check"
                   % (probe, why))
    R.bump()
    if literal_hits("pct = 100.0 * num / den", forbidden):
        R.fail("NOLITERALS self-test: the percent scale factor in "
               "'100.0 * num / den' was reported as a coverage literal")

    scripts = figure_scripts()
    R.bump(len(scripts))
    if not scripts:
        R.note("no figure script was found under %s; nothing to scan"
               % ", ".join(rel(s) for s in FIGURE_SCOPES))
    for path in scripts:
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                raw = fh.read()
        except OSError as exc:
            R.fail("could not read %s (%s)" % (rel(path), exc))
            continue
        code = strip_comments(path, raw)
        for ln, line in enumerate(code.splitlines(), 1):
            low = line.lower()
            # RULE A -- the retired figures, unconditional
            for pat, why in RETIRED_PATTERNS:
                if re.search(pat, line):
                    R.fail("%s:%d hard-codes %s" % (rel(path), ln, why),
                           [line.strip()[:160]])
            for pat, ctxpat, why in RETIRED_CONTEXTUAL:
                if re.search(pat, line) and re.search(ctxpat, low):
                    R.fail("%s:%d hard-codes %s" % (rel(path), ln, why),
                           [line.strip()[:160]])
            # RULE B -- anything the deposit computes
            if any(tok in low for tok in BENIGN_CONTEXT):
                continue
            for tok, why in literal_hits(line, forbidden):
                R.fail("%s:%d hard-codes %r, which is %s. Every published "
                       "figure must be one call of atlas_coverage() over "
                       "Data/coverage_basis.csv, or read "
                       "Paper/inputs/headline_numbers.json."
                       % (rel(path), ln, tok, why),
                       [line.strip()[:160]])

        # A figure script that plots coverage but reads no generated input is
        # a figure script that will drift.
        if re.search(r"coverage|\bC0\b|\bC1\b|dark matter", raw, re.I) \
                and not re.search(r"Paper/inputs|headline_numbers\.json|"
                                  r"coverage_basis\.csv|atlas_coverage", raw):
            R.note("%s mentions coverage but reads none of "
                   "Paper/inputs/, headline_numbers.json, coverage_basis.csv "
                   "or atlas_coverage(); check where its numbers come from"
                   % rel(path))


def check_noabspaths(ctx):
    files = []
    for scope in CODE_SCOPES:
        files.extend(walk_files(scope, (".py", ".r", ".R", ".sh", ".tex",
                                        ".yaml", ".yml", ".json", ".md")))
    R.bump(len(files))
    if not files:
        R.fail("no source file was found under %s; the scan scope is wrong"
               % ", ".join(rel(s) for s in CODE_SCOPES))
    for path in files:
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except OSError:
            continue
        # Strip comments and docstrings first, for the same reason as the
        # literal checks: a path quoted in prose that DOCUMENTS its own removal
        # is not a path the code can follow.  Strings that are values are still
        # scanned, so a live os.path.join("/Users/...") is still caught.
        scan_text = strip_comments(path, text)
        raw_lines = text.splitlines()
        for ln, line in enumerate(scan_text.splitlines(), 1):
            if any(m in line for m in HOME_PATH_MARKERS):
                shown = raw_lines[ln - 1] if ln <= len(raw_lines) else line
                R.fail("%s:%d contains an absolute home path; the deposit must "
                       "run from a clone at any location"
                       % (rel(path), ln), [shown.strip()[:160]])


def check_retired(ctx):
    """No retired figure may survive in live figure code or live manuscript
    sources.  In prose ONLY, an occurrence inside an explicit retirement
    disclosure is allowed and reported -- the paper is expected to volunteer
    the correction rather than hide it."""
    files = []
    for scope in LIVE_TEX_SCOPES:
        files.extend(walk_files(scope, (".tex", ".md")))
    files.extend(figure_scripts())
    files = [f for f in files if os.path.abspath(f) not in TEXT_SCAN_EXEMPT]
    R.bump(len(files))
    if not files:
        R.note("no live manuscript source was found under %s"
               % ", ".join(rel(s) for s in LIVE_TEX_SCOPES))
    disclosed = 0
    for path in files:
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                raw = fh.read()
        except OSError:
            continue
        text = strip_comments(path, raw)
        is_prose = path.lower().endswith((".tex", ".md"))
        lines = text.splitlines()
        for ln, line in enumerate(lines, 1):
            hits = [why for pat, why in RETIRED_PATTERNS if re.search(pat, line)]
            hits += [why for pat, cpat, why in RETIRED_CONTEXTUAL
                     if re.search(pat, line) and re.search(cpat, line.lower())]
            if not hits:
                continue
            window = " ".join(lines[max(0, ln - 3):ln + 2])
            if is_prose and RETIREMENT_MARKERS.search(window):
                disclosed += 1
                continue
            for why in hits:
                R.fail("%s:%d carries %s" % (rel(path), ln, why),
                       [line.strip()[:160]])
    if disclosed:
        R.note("%d occurrence(s) of a retired figure appear inside an explicit "
               "retirement disclosure and are allowed. Volunteering the "
               "correction is the right move; make sure each one says "
               "WITHDRAWN, not adjusted." % disclosed)


# ==========================================================================
# MAIN
# ==========================================================================


# ------------------------------------------------------------ DATASETS ---
DATASET_TABLE_CSV = P("Paper", "inputs", "table_datasets.csv")
DATASET_TABLE_TEX = P("Paper", "inputs", "table_datasets.tex")
DATASET_REGISTRY = P("Data", "dataset_table_registry.csv")
TRUTH_DATED_CSV = P("Data", "_rebuild", "truth_dated.csv")
MEE_WRAPPERS = [P("Paper", "MEE", "main.tex"), P("Paper", "MEE", "main_anon.tex")]


def check_datasets(ctx):
    """F01-F07.  The main-text dataset table is generated; re-derive every
    count in it from the basis table, the measured truth and the S09
    resolution, independently of the rebuild, and check that the MEE
    wrappers input the generated table and not a hand-maintained one."""
    rows = read_csv(DATASET_TABLE_CSV)
    truth = OrderedDict((r["filename"], r) for r in read_csv(TRUTH_DATED_CSV))
    scales = {r["filename"]: r["time_scale"] for r in read_csv(TIME_SCALE_CHECKS)}
    nonc = defaultdict(set)
    seen_files = set()
    for r in read_csv(COL_RESOLUTION_NONCANON):
        seen_files.add(r["filename"])
        if r.get("match_type") in COL_RESOLUTION_RESOLVED and r.get("col_id"):
            nonc[r["filename"]].add(r["col_id"])
    by_file = {}
    for r in rows:
        if r["filename"] in truth:
            R.eq("dataset table lists %s once" % r["filename"],
                 r["filename"] in by_file, False)
            by_file[r["filename"]] = r
    # every shipped tree is either its own row or covered by an aggregate row
    agg = [r for r in rows if "aggregate" in r["flags"].split(";")]
    covered = set(by_file)
    for a in agg:
        pat = re.compile("^" + re.escape(a["filename"]).replace(r"\*", ".*") + "$")
        files = [f for f in truth if pat.match(f)]
        covered |= set(files)
        n_chrono = sum(1 for f in files if truth[f]["is_chronogram"] == "yes")
        R.eq("aggregate %s tree count" % a["filename"], a["tips"],
             "%d trees" % len(files))
        R.eq("aggregate %s dated count" % a["filename"], a["dated"],
             "%d of %d" % (n_chrono, len(files)))
    R.eq("every shipped tree appears in the dataset table",
         sorted(set(truth) - covered), [])
    for part, b in ctx.basis.items():
        fn = b.get("canonical_tree") or ""
        if not fn:
            continue
        r = by_file.get(fn)
        if r is None:
            R.fail("canonical tree %s (%s) is missing from the dataset table"
                   % (fn, part))
            continue
        R.eq("%s canonical" % fn, r["is_canonical"], "yes")
        R.eq("%s tips = distinct labels" % fn, as_int(r["tips"]),
             as_int(truth[fn]["distinct_labels"]))
        R.eq("%s names = c0_units" % fn, as_int(r["accepted_names"]),
             as_int(b["c0_units"]))
        R.eq("%s described = basis" % fn, as_int(r["described"]),
             as_int(b["described"]))
        R.eq("%s described source = basis" % fn, r["described_source"],
             b["described_source"])
        R.eq("%s coverage" % fn, r["coverage_pct"],
             fnum(100.0 * as_int(b["c0_units"]) / as_int(b["described"]), 1))
        want = ("yes" if b["is_chronogram"] == "yes" else
                "relative" if scales.get(fn) == "relative" else "no")
        R.eq("%s dated" % fn, r["dated"], want)
    canon_files = {b.get("canonical_tree") for b in ctx.basis.values()}
    for fn, r in by_file.items():
        if fn in canon_files:
            continue
        R.eq("%s not canonical" % fn, r["is_canonical"], "no")
        R.eq("%s tips" % fn, as_int(r["tips"]),
             as_int(truth[fn]["distinct_labels"]))
        want_names = len(nonc[fn]) if fn in seen_files else ""
        R.eq("%s accepted names" % fn,
             as_int(r["accepted_names"]) if r["accepted_names"] != "" else "",
             want_names)
        R.eq("%s dated" % fn, r["dated"],
             "ultrametric" if truth[fn]["is_chronogram"] == "yes" else "no")
        if r["described"]:
            R.eq("%s described has a source" % fn, bool(r["described_source"]),
                 True)
    unrep = sorted(r["group"] for r in rows if r["section"] == "Not yet represented")
    R.eq("unrepresented rows = basis is_unrepresented (eukaryotic)", unrep,
         sorted(p for p, b in ctx.basis.items()
                if b.get("is_prokaryotic") == "no"
                and b.get("is_unrepresented") == "yes"))
    for r in rows:
        if r["section"] == "Not yet represented":
            R.eq("%s unrepresented described" % r["group"],
                 as_int(r["described"]),
                 as_int(ctx.basis[r["group"]]["described"]))
    tex = open(DATASET_TABLE_TEX, encoding="utf-8").read()
    body = [ln for ln in tex.splitlines() if ln.endswith(r"\\")
            and not ln.startswith(r"\textbf{\textit{")]
    R.eq("dataset table .tex has one line per csv row", len(body), len(rows))
    # The revised Research Article carries the dataset table in Supplement S2
    # (OUTLINE.md section 8); its main-text tables are Tables 1-3, generated
    # and checked by TABLES.  The MEE wrappers must input the revised body,
    # the body must input the three generated tables, and no MEE file may
    # input a hand-maintained table.
    body = open(MANUSCRIPT_BODY, encoding="utf-8").read()
    for t in ("table_which_tree.tex", "table_rungs_sensitivity.tex",
              "table_worked_example.tex"):
        R.eq("the manuscript body inputs the generated %s" % t,
             ("{../inputs/%s}" % t) in body, True)
    R.eq("the manuscript body does not input a hand-maintained table",
         "main_table_rows" in body, False)
    for w in MEE_WRAPPERS:
        s = open(w, encoding="utf-8").read()
        R.eq("%s inputs the revised body" % rel(w),
             r"\input{../shared/manuscript_body_v2.tex}" in s, True)
        R.eq("%s does not input a hand-maintained table" % rel(w),
             "main_table_rows" in s, False)
    sbody = P("Paper", "shared", "supplementary_body_v2.tex")
    R.eq("the Supplement inputs the generated dataset table",
         os.path.exists(sbody) and "{../inputs/table_datasets.tex}"
         in open(sbody, encoding="utf-8").read(), True)


# ------------------------------------------------------------- FIGURES ---
FIGURE_DIR = P("Paper", "figures")
FIGURE_CODE = P("Paper", "code", "figures")
MANUSCRIPT_BODY = P("Paper", "shared", "manuscript_body_v2.tex")
MANUSCRIPT_META = P("Paper", "shared", "metadata_v2.tex")
MANUSCRIPT_TEX = [MANUSCRIPT_BODY, MANUSCRIPT_META] + MEE_WRAPPERS + [P("Paper", "bioRxiv", "main.tex")]


def check_figures(ctx):
    """F01.  Every figure the MEE manuscript draws exists and is written by a
    script under Paper/code/figures/; the figure macro carries no absolute
    path and no placeholder survives."""
    used = set()
    for path in MANUSCRIPT_TEX:
        s = open(path, encoding="utf-8").read()
        s = re.sub(r"(?m)(?<!\\)%.*$", "", s)
        for m in re.finditer(r"\\(?:AtlasFigure|inlinefigure|includegraphics)"
                             r"(?:\[[^\]]*\])?\{([^}#]+\.(?:pdf|png))\}", s):
            used.add(m.group(1))
        R.eq("%s has no placeholder figure text" % rel(path),
             "pending regeneration" in s, False)
    blocks = re.findall(r"\\newcommand\{\\(\w+FigureBlock)\}",
                        open(MEE_WRAPPERS[0], encoding="utf-8").read())
    body = open(MANUSCRIPT_BODY, encoding="utf-8").read()
    for blk in blocks:
        if ("\\" + blk) not in body:
            # defined but never placed: its graphic is not part of the paper
            m = re.search(r"\\newcommand\{\\%s\}\{\\inlinefigure\{([^}]+)\}" % blk,
                          open(MEE_WRAPPERS[0], encoding="utf-8").read())
            if m:
                used.discard(m.group(1))
    R.eq("the manuscript draws at least one figure", bool(used), True)
    scripts = {}
    for fn in os.listdir(FIGURE_CODE):
        if fn.endswith((".py", ".R", ".r")):
            scripts[fn] = open(os.path.join(FIGURE_CODE, fn), encoding="utf-8").read()
    for g in sorted(used):
        R.eq("figure %s exists" % g,
             os.path.exists(os.path.join(FIGURE_DIR, g)), True)
        stem = os.path.splitext(g)[0]
        writers = [n for n, t in scripts.items() if '"%s"' % stem in t
                   or "'%s'" % stem in t or "%s." % stem in t]
        R.eq("figure %s is written by a script in %s" % (g, rel(FIGURE_CODE)),
             bool(writers), True)
    meta = open(MANUSCRIPT_META, encoding="utf-8").read()
    m = re.search(r"\\providecommand\{\\AtlasFigDir\}\{([^}]*)\}", meta)
    R.eq("AtlasFigDir is defined", bool(m), True)
    if m:
        R.eq("AtlasFigDir is relative", m.group(1).startswith(("/", "~")), False)


# --------------------------------------------------------------------------
# CROWNREG (audit A13): one crown-age register
# --------------------------------------------------------------------------
CROWN_AGES_CSV = P("Paper", "inputs", "crown_ages.csv")
CROWN_MACROS = P("Paper", "inputs", "crown_ages_macros.tex")
CROWN_POINTS = P("Paper", "inputs", "fig_crown_age_points.csv")
CROWN_MANIFEST = P("Paper", "inputs", "fig2_tree_manifest.csv")
SUPERSEDED_CROWN = P("Data", "_archive", "2026-10-revision",
                     "superseded_crown_register")
# Directories whose CSV files are scanned for a crown-age column.
CROWNREG_SCAN_DIRS = (P("Paper", "inputs"), P("Data"), P("Data", "standardized"),
                      P("Data", "_rebuild"))
# Text files that may not cite a superseded register.
CROWNREG_TEXT_GLOBS = ("scripts/*.py", "Paper/code/figures/*", "Paper/code/analysis/*",
                       "Paper/shared/*.tex", "Paper/MEE/*.tex")
CROWN_COL_RE = re.compile(r"crown", re.I)
CROWN_NONAGE_COLS = ("crown_clade",)


def crown_age_columns(header):
    """The columns of a CSV header that can hold a crown age."""
    return [c for c in header if CROWN_COL_RE.search(c or "")
            and c not in CROWN_NONAGE_COLS]


def crown_macro_name(partition):
    return "".join(w[:1].upper() + w[1:]
                   for w in re.split(r"[^A-Za-z]+", partition) if w)


def check_crownreg(ctx):
    """A13.  Paper/inputs/crown_ages.csv is the only crown-age register.

    1. Every CSV in Paper/inputs/, Data/, Data/standardized/ and
       Data/_rebuild/ with a crown-age column is one of the files below, and
       each of those agrees with crown_ages.csv for the same partition within
       CROWN_VALUE_TOL.  A file with a crown-age column that is not listed is a
       second register and fails.
    2. The generated macros (crown_ages_macros.tex) carry the register's
       values.
    3. No live script, figure code or manuscript source cites a file retired
       to Data/_archive/2026-10-revision/superseded_crown_register/."""
    ages = OrderedDict((r["partition"], r) for r in read_csv(CROWN_AGES_CSV))
    tol = CROWN_VALUE_TOL + 1e-9      # both sides written at 4 dp

    def age(part, k):
        return as_float(ages[part].get(k))

    def drawn_interval(a):
        lo, hi = as_float(a.get("interval_lo_ma")), as_float(a.get("interval_hi_ma"))
        if lo is None or hi is None or hi - lo <= CROWN_DEGENERATE_MYR:
            return None
        return lo, hi

    def h_self(path, rows):
        R.eq("crown_ages.csv partitions unique", len(rows), len(ages))

    def h_points(path, rows):
        canon = [r for r in rows if r.get("partition")]
        R.eq("fig_crown_age_points.csv partitions == crown_ages.csv",
             sorted(r["partition"] for r in canon), sorted(ages))
        for r in canon:
            a = ages.get(r["partition"])
            if not a:
                continue
            R.eq("%s figure crown age" % r["partition"], r.get("crown_age"),
                 a["crown_age_ma"], atol=tol)
            iv = drawn_interval(a)
            if iv:
                R.eq("%s figure interval lo" % r["partition"],
                     r.get("crown_age_low"), iv[0], atol=tol)
                R.eq("%s figure interval hi" % r["partition"],
                     r.get("crown_age_high"), iv[1], atol=tol)
                R.eq("%s figure interval counts_for_c1b" % r["partition"],
                     r.get("counts_for_c1b"), a.get("counts_for_c1b"))
            else:
                R.eq("%s figure draws no interval" % r["partition"],
                     (r.get("crown_age_low") or "") + (r.get("crown_age_high") or ""), "")
        for r in rows:
            if not r.get("partition"):
                R.eq("%s (context tree) draws no interval" % r.get("name"),
                     (r.get("crown_age_low") or "") + (r.get("crown_age_high") or ""), "")

    def h_manifest(path, rows):
        canon = [r for r in rows if r.get("source_kind") == "canonical"]
        R.eq("crown manifest canonical rows == crown_ages.csv",
             sorted(r.get("partition") for r in canon), sorted(ages))
        for r in canon:
            a = ages.get(r.get("partition"))
            if a:
                R.eq("%s manifest file" % r["partition"], r.get("file"),
                     a["canonical_tree"])
            R.eq("%s manifest asserts no crown age" % r.get("partition"),
                 r.get("crown_age_fallback", ""), "")
        for r in rows:
            if r.get("source_kind") != "canonical":
                R.eq("%s (context tree) names no partition" % r.get("name"),
                     r.get("partition", ""), "")

    def h_worked_sel(path, rows):
        for r in rows:
            if r.get("partition") in ages:
                R.eq("%s worked-example crown age" % r["partition"],
                     r.get("crown_age_ma"), ages[r["partition"]]["crown_age_ma"],
                     atol=tol)
            else:
                R.fail("%s: partition %r has no crown_ages.csv row"
                       % (rel(path), r.get("partition")))

    def h_worked(path, rows):
        for r in rows:
            if r.get("partition") not in ages:
                R.fail("%s: partition %r has no crown_ages.csv row"
                       % (rel(path), r.get("partition")))
            elif r.get("tree") == "full":
                R.eq("%s worked-example full-tree crown age" % r["partition"],
                     r.get("crown_age_ma"), ages[r["partition"]]["crown_age_ma"],
                     atol=tol)

    def h_whichtree(path, rows):
        # Table 1 prints each chronogram's crown age and interval at one
        # decimal: it must be crown_ages.csv rounded, never a second source
        reg = {r["partition"]: r for r in read_csv(CROWN_AGES_CSV)}
        for r in rows:
            c = reg.get(r["partition"])
            if not r["crown_age_ma"]:
                R.eq("%s: %s has no crown age only when it is not a "
                     "chronogram" % (rel(path), r["partition"]), c, None)
                continue
            R.eq("%s: %s listed in crown_ages.csv" % (rel(path), r["partition"]),
                 c is not None, True)
            if c is None:
                continue
            for k in ("crown_age_ma", "interval_lo_ma", "interval_hi_ma"):
                want = "%.1f" % as_float(c[k]) if c[k] else ""
                R.eq("%s: %s %s" % (rel(path), r["partition"], k), r[k], want)

    handlers = {
        TABLE_WHICH: h_whichtree,
        P("Paper", "inputs", "supp_crown.csv"): h_whichtree,
        CROWN_AGES_CSV: h_self,
        CROWN_POINTS: h_points,
        CROWN_MANIFEST: h_manifest,
        P("Paper", "inputs", "worked_example_bd_selection.csv"): h_worked_sel,
        P("Paper", "inputs", "worked_example_bd.csv"): h_worked,
    }
    # self-test: the classifier sees the retired register's age columns
    R.eq("crown_age_columns self-test",
         crown_age_columns(["dataset", "crown_age_point", "crown_age_low",
                            "crown_ma", "crown_clade"]),
         ["crown_age_point", "crown_age_low", "crown_ma"])
    seen = set()
    for d in CROWNREG_SCAN_DIRS:
        for path in sorted(glob.glob(os.path.join(d, "*.csv"))):
            with open(path, newline="", encoding="utf-8", errors="replace") as fh:
                header = next(csv.reader(fh), [])
            cols = crown_age_columns(header)
            if not cols:
                continue
            seen.add(path)
            fn = handlers.get(path)
            R.bump()
            if fn is None:
                R.fail("%s carries crown-age column(s) %s but is not "
                       "crown_ages.csv or a file checked against it: a second "
                       "crown-age register" % (rel(path), ", ".join(cols)))
                continue
            fn(path, read_csv(path))
    for path in handlers:
        R.eq("%s present and scanned" % rel(path), path in seen, True)
    # 2. macros
    want = {}
    for part, a in ages.items():
        n = crown_macro_name(part)
        g = as_float(a["crown_age_ma"])
        want["CrownAge" + n] = "%.1f" % g
        lo, hi = as_float(a.get("interval_lo_ma")), as_float(a.get("interval_hi_ma"))
        if lo is not None and hi is not None:
            want["CrownLo" + n] = "%.1f" % lo
            want["CrownHi" + n] = "%.1f" % hi
            want["CrownWidth" + n] = "%.1f" % (hi - lo)
            want["CrownWidthPct" + n] = "%.1f" % (100.0 * (hi - lo) / g)
    got = {}
    if os.path.exists(CROWN_MACROS):
        got = dict(re.findall(r"\\newcommand\{\\(\w+)\}\{([^}]*)\}",
                              open(CROWN_MACROS, encoding="utf-8").read()))
    R.eq("crown_ages_macros.tex == crown_ages.csv", got, want)
    # 3. no live citation of a superseded register
    retired = sorted({os.path.basename(f) for f in
                      glob.glob(os.path.join(SUPERSEDED_CROWN, "**", "*"),
                                recursive=True) if os.path.isfile(f)}
                     - {"data_provenance.csv", "fig2_tree_manifest.csv"})
    R.eq("superseded crown registers archived", bool(retired), True)
    for g in CROWNREG_TEXT_GLOBS:
        for path in sorted(glob.glob(P(g))):
            if not os.path.isfile(path) or path.endswith((".pdf", ".png", ".log",
                                                          ".aux", ".bbl")):
                continue
            text = open(path, encoding="utf-8", errors="replace").read()
            for name in retired:
                R.bump()
                if name in text:
                    R.fail("%s cites the retired crown register %s"
                           % (rel(path), name))


# --------------------------------------------------------------------------
# NUMMACROS: the prose number macros re-derived
# --------------------------------------------------------------------------
NUMBER_MACROS = P("Paper", "inputs", "numbers_macros.tex")
GENERATED_MACRO_FILES = (NUMBER_MACROS, CROWN_MACROS,
                         P("Paper", "inputs", "table_datasets_macros.tex"))


def _mpct(x):
    return "%.1f" % x if abs(x) >= 1.0 else "%.2f" % x


def _mint(n):
    return "{:,}".format(int(round(float(n)))).replace(",", "{,}")


def _macros_in(path):
    if not os.path.exists(path):
        return []
    # the value runs to the last brace on the line: it may hold {,}
    return re.findall(r"^\\newcommand\{\\(\w+)\}\{(.*)\}\s*$",
                      open(path, encoding="utf-8").read(), re.M)


UNREP_NAMED = 4      # largest groups without a tree named in the text


def union_losses(ctx, canon_e, crows):
    """Levels by category and losses by reason on the union convention,
    from the measured name sets (ctx.c0_names, ctx.c1_sets)."""
    crown = {r["partition"]: r for r in crows}
    parts = sorted(r["partition"] for r in canon_e)
    row = {r["partition"]: r for r in canon_e}

    def reason(p):
        r, c = row[p], crown.get(p)
        if r["is_chronogram"] != "yes":
            return "RelTime" if r["time_scale"] == "relative" else "NotChron"
        if c is None or not c["interval_kind"]:
            return "NoInterval"
        if c["counts_for_c1b"] == "yes":
            return ""
        return "Kind" if c["c1b_exclusion"] == "kind" else "Calib"
    u = {k: set() for k in ("C0", "C1", "C1d", "C1b")}
    for p in parts:
        c0, c1 = set(ctx.c0_names[p]), ctx.c1_sets.get(p)
        R.eq("%s C1 set enumerated for the union losses" % p,
             c1 is not None, True)
        c1 = set(c1 or ())
        u["C0"] |= c0
        u["C1"] |= c1
        if row[p]["is_chronogram"] == "yes":
            u["C1d"] |= c1
            if reason(p) == "":
                u["C1b"] |= c1
    lost = Counter()
    holder = {}
    for p in parts:
        c1 = set(ctx.c1_sets.get(p) or ())
        for sp in set(ctx.c0_names[p]):
            lvl = ("C1b" if sp in u["C1b"] else "C1d" if sp in u["C1d"]
                   else "C1" if sp in u["C1"] else "C0")
            here = ("C0" if sp not in c1 else "C1" if row[p]["is_chronogram"]
                    != "yes" else "C1b" if reason(p) == "" else "C1d")
            if here == lvl and sp not in holder:
                holder[sp] = p
    levels = {k: Counter() for k in u}
    for sp in u["C0"]:
        p = holder[sp]
        cat = row[p]["category"]
        for k in ("C0", "C1", "C1d", "C1b"):
            if sp in u[k]:
                levels[k][cat] += 1
        if sp not in u["C1"]:
            lost["TaxPlaced"] += 1
        elif sp not in u["C1b"]:
            lost[reason(p)] += 1
    return lost, levels


def _names(ps):
    ps = sorted(ps)
    if not ps:
        return ""
    return ps[0] if len(ps) == 1 else ", ".join(ps[:-1]) + " and " + ps[-1]


def _worked_poly_macros(eligible, fits):
    """Polytomy macros, formed here from the fit files of worked_example_bd.R."""
    w = OrderedDict()
    poly = read_csv(WORKED_POLY)
    res = defaultdict(lambda: defaultdict(dict))
    for r in poly:
        res[r["partition"]][as_int(r["seed"])][(r["tree"], r["rho_convention"])] = r
    with_poly, allr = [], []
    nb = nrej = 0
    for p in eligible:
        ff, fm = fits[p + "|full|tips"], fits[p + "|molecular|tips"]
        if as_int(ff["n_polytomy_events"]) or as_int(fm["n_polytomy_events"]):
            with_poly.append(p)
        lrs, rr, bnd, rej = [], [], [], []
        for sd, f in sorted(res[p].items()):
            lr = 2.0 * (as_float(f[("molecular", "tips")]["loglik"])
                        - as_float(f[("molecular", "full_misapplied")]["loglik"]))
            em = as_float(f[("molecular", "full_misapplied")]["epsilon"])
            lrs.append(round(lr, 4))
            rr.append(round(as_float(f[("full", "tips")]["r"])
                            / as_float(f[("molecular", "tips")]["r"]), 4))
            bnd.append(em < 1e-6 and lr > 2e-3)
            rej.append(em < 1e-6 and lr > WORKED_LR_CUTOFF)
        R.eq("worked polytomy %s: one resolved fit per seed and convention" % p,
             bool(res[p]) and all(len(v) == 3 for v in res[p].values()), True)
        allr += rr
        nb += all(bnd)
        nrej += all(rej)
        if p in with_poly:
            n = crown_macro_name(p)
            for tree, nm in (("full", "Full"), ("molecular", "Pruned")):
                f = fits["%s|%s|tips" % (p, tree)]
                w["WEPolyPct%s%s" % (nm, n)] = _mpct(
                    100.0 * as_int(f["n_polytomy_events"]) / as_int(f["n_events"]))
            ft_, fc_ = fits[p + "|molecular|tips"], fits[p + "|molecular|full_misapplied"]
            if as_float(fc_["epsilon"]) < 1e-6 and abs(as_float(ft_["loglik"])
                                                       - as_float(fc_["loglik"])) > 1e-3:
                w["WEPolyLRMin" + n] = "%.1f" % min(lrs)
                w["WEPolyLRMax" + n] = "%.1f" % max(lrs)
            w["WEPolyRRatioMin" + n] = "%.2f" % min(rr)
            w["WEPolyRRatioMax" + n] = "%.2f" % max(rr)
    w["WEPolyNGroups"] = _mint(len(with_poly))
    w["WEPolyNames"] = _names(with_poly)
    w["WEPolyRRatioMin"] = "%.2f" % min(allr)
    w["WEPolyRRatioMax"] = "%.2f" % max(allr)
    w["WEPolyNBoundary"] = _mint(nb)
    w["WEPolyNRejected"] = _mint(nrej)
    return w


def _check_supp_polytomy(byk):
    """Supplement table of the polytomy-resolution sensitivity, re-derived
    from the fit files."""
    poly = read_csv(WORKED_POLY)
    got = {r["partition"]: r for r in _supp_csv("supp_worked_polytomy")}
    parts = sorted({k[0] for k in byk})
    R.eq("Supplement S9 polytomy table: one row per eligible group",
         sorted(got), parts)

    def rng(vals, fmt):
        lo, hi = fmt % min(vals), fmt % max(vals)
        return lo if lo == hi else "%s--%s" % (lo, hi)
    for p in parts:
        g = got.get(p)
        if g is None:
            continue
        seeds = defaultdict(dict)
        for r in poly:
            if r["partition"] == p:
                seeds[r["seed"]][(r["tree"], r["rho_convention"])] = r
        rr, em, lr, bnd, rej = [], [], [], [], []
        for f in seeds.values():
            rr.append(round(as_float(f[("full", "tips")]["r"])
                            / as_float(f[("molecular", "tips")]["r"]), 4))
            em.append(as_float(f[("molecular", "tips")]["epsilon"]))
            x = round(2.0 * (as_float(f[("molecular", "tips")]["loglik"])
                             - as_float(f[("molecular", "full_misapplied")]["loglik"])), 4)
            lr.append(x)
            e_ = as_float(f[("molecular", "full_misapplied")]["epsilon"])
            bnd.append(e_ < 1e-6 and x > 2e-3)
            rej.append(e_ < 1e-6 and x > WORKED_LR_CUTOFF)
        ff, fm = byk[(p, "full", "tips")], byk[(p, "molecular", "tips")]
        fc = byk[(p, "molecular", "full_misapplied")]
        want = OrderedDict([
            ("events_full", ff["n_events"]),
            ("polytomy_events_full", ff["n_polytomy_events"]),
            ("events_pruned", fm["n_events"]),
            ("polytomy_events_pruned", fm["n_polytomy_events"]),
            ("r_ratio", "%.2f" % round(as_float(ff["r"]) / as_float(fm["r"]), 4)),
            ("r_ratio_resolved", rng(rr, "%.2f")),
            ("epsilon_matching", "%.2f" % as_float(fm["epsilon"])),
            ("epsilon_matching_resolved", rng(em, "%.2f")),
            ("lr", "%.2f" % round(2.0 * (as_float(fm["loglik"])
                                         - as_float(fc["loglik"])), 4)),
            ("lr_resolved", rng(lr, "%.2f")),
            ("outcome_resolved", "mu0_rejected" if all(rej) else "mu0"
             if all(bnd) else "rescaled" if not any(bnd) else "mixed"),
        ])
        for k, v in want.items():
            R.eq("Supplement S9 polytomy table %s %s" % (p, k), g[k], v)


def check_numbermacros(ctx):
    """Every number the manuscript prose quotes is a macro in
    Paper/inputs/numbers_macros.tex.  This check re-derives each one: counts
    from the basis table and the headline numerators, percentages formed here
    from numerator and denominator (never copied from a 'percent' field),
    worked-example estimates from the fitted values.  It also requires that no
    macro name is defined twice across the generated macro files."""
    h, basis = ctx.head, ctx.basis_rows
    euk = [r for r in basis if r["is_prokaryotic"] == "no"]
    D = sum(as_int(r["described"], 0) for r in euk)
    w = OrderedDict()
    w["NumPartitionsAll"] = _mint(len(basis))
    w["NumPartitionsEuk"] = _mint(len(euk))
    w["NumPartitionsProk"] = _mint(len(basis) - len(euk))
    tr = h["trees"]
    w["NumTreesShipped"] = _mint(tr["n_shipped"])
    w["NumTreesChron"] = _mint(tr["n_chronograms_all_shipped"])
    canon = [r for r in basis if as_int(r["c0_units"], 0) > 0]
    canon_e = [r for r in canon if r["is_prokaryotic"] == "no"]
    w["NumCanonical"] = _mint(len(canon))
    w["NumCanonicalEuk"] = _mint(len(canon_e))
    w["NumCanonicalChron"] = _mint(sum(r["is_chronogram"] == "yes"
                                       for r in canon))
    w["NumCanonicalRelTime"] = _mint(len(tr["relative_time_partitions"]))
    w["DescribedEuk"] = _mint(D)
    w["DescribedProk"] = _mint(sum(as_int(r["described"], 0) for r in basis
                                   if r["is_prokaryotic"] != "no"))
    Drep = sum(as_int(r["described"], 0) for r in canon_e)
    w["DescribedRepresented"] = _mint(Drep)
    names = (("C0", "CZero"), ("C1", "COne"), ("C1d", "COneD"),
             ("C1b", "COneB"))
    rg = h["rungs"]
    for key, nm in names:
        num = rg[key]["numerator_union"]
        w[nm + "N"] = _mint(num)
        w[nm + "Pct"] = _mpct(100.0 * num / D)
        w[nm + "NPart"] = _mint(rg[key]["n_partitions_contributing"])
    bases = [r["c1_basis"] for r in canon_e]
    w["COneNEnumerated"] = _mint(bases.count("enumerated"))
    w["COneNAllTips"] = _mint(bases.count("all_tips"))
    asserted = sum(as_float(r["c1_units"]) for r in canon_e
                   if r["c1_basis"] == "all_tips")
    w["COneAssertedN"] = _mint(asserted)
    w["COneAssertedPct"] = _mpct(100.0 * asserted / rg["C1"]["numerator_sum"])
    w["UnrepresentedPct"] = _mpct(100.0 * (D - Drep) / D)
    for (a, an), (b, bn) in zip(names[1:], names[:-1]):
        w[an + "Of" + bn] = _mpct(100.0 * rg[a]["numerator_union"]
                                  / rg[b]["numerator_union"])
    sp = OrderedDict((r["partition"], r) for r in basis)["Seed plants"]
    c0sp, c1sp = as_int(sp["c0_units"], 0), as_float(sp["c1_units"])
    sp_vals = {"C0": c0sp, "C1": c1sp,
               "C1d": c1sp if sp["is_chronogram"] == "yes" else 0.0,
               "C1b": c1sp if (sp["is_chronogram"] == "yes" and
                               sp["interval_recoverable"] == "yes") else 0.0}
    for key, nm in names:
        w["SeedPlantsShare" + nm] = _mpct(100.0 * sp_vals[key]
                                          / rg[key]["numerator_sum"])
    # the levels without the seed-plant tree and its described species
    for key, nm in names[:3]:
        w["ExclSeedPlants" + nm + "Pct"] = _mpct(
            100.0 * (rg[key]["numerator_sum"] - sp_vals[key])
            / (D - as_int(sp["described"], 0)))
    for cat, v in h["per_category_eukaryotic"].items():
        cn = crown_macro_name(cat)
        for key, nm in names[:3]:
            w["Cat%s%sPct" % (cn, nm)] = _mpct(100.0 * v["%s_sum" % key]
                                               / v["described"])
    # R1-11, R1-06, R1-09
    w["TipsNotNamedN"] = _mint(sum(as_int(r["labels_unique"], 0)
                                   - as_int(r["c0_units"], 0) for r in canon_e))
    w["SensTipsCZeroPct"] = _mpct(100.0 * sum(
        min(as_int(r["labels_unique"], 0), as_int(r["described"], 0))
        for r in canon_e) / D)
    for r in canon_e:
        w["TipsNotNamed" + crown_macro_name(r["partition"])] = _mint(
            as_int(r["labels_unique"], 0) - as_int(r["c0_units"], 0))
    nmol_ = Counter(r["partition"] for r in read_csv(MOLTIPS))
    taxo = sorted(r["partition"] for r in canon_e
                  if r["c1_basis"] == "enumerated"
                  and nmol_[r["partition"]] < as_int(r["labels_unique"], 0))
    mo = [r["partition"] for r in read_csv(MOLECULAR_ONLY_TREES)]
    R.eq("R1-06 molecular-only trees listed for exactly the groups with "
         "taxonomically placed tips", sorted(mo), taxo)
    w["NGroupsTaxonomicTips"] = _mint(len(taxo))
    w["NMolOnlyTrees"] = _mint(len(mo))
    sdv = {r["partition"]: (r.get("study_dating") or "").strip()
           for r in read_csv(DATING_DEPOSIT_CHECKS)}
    for key, nm in STUDY_DATING_MACRO.items():
        w["StudyDating%sN" % nm] = _mint(sum(
            1 for r in canon_e if sdv.get(r["partition"]) == key
            and not (r["is_chronogram"] == "yes"
                     and r["time_scale"] == "absolute")))
    nr = h["name_resolution_S05"]["eukaryotic"]["unresolved_share"]
    w["LabelsEukN"] = _mint(nr["denominator"])
    w["UnresolvedLabelsN"] = _mint(nr["numerator"])
    w["UnresolvedLabelsPct"] = _mpct(100.0 * nr["numerator"] / nr["denominator"])
    s = h["sensitivities"]
    g = s["c0_gbif_unmatched"]
    w["SensGbifSpeciesAdded"] = _mint(g["species_added"])
    w["SensGbifCZeroPct"] = _mpct(100.0 * g["numerator"] / D)
    col = s["denominator"]["col_single_authority"]["filled"]
    w["SensColDenominator"] = _mint(col["denominator"])
    for key, nm in names:
        cr = col["rungs"][key]
        w["SensCol%sPct" % nm] = _mpct(100.0 * cr.get("numerator_union",
                                                      cr["numerator"])
                                       / col["denominator"])
        w["SensRep%sPct" % nm] = _mpct(100.0 * rg[key]["numerator_union"]
                                       / Drep)
    u = h["union_bound_S09"]
    w["SensUnionNTrees"] = _mint(u["n_noncanonical_trees"])
    w["SensUnionNamesAdded"] = _mint(u["C0_union"] - u["C0_canonical"])
    w["SensUnionCZeroPct"] = _mpct(100.0 * u["C0_union"] / D)
    w["SensUnionCOneHighPct"] = _mpct(100.0 * u["C1_union_high"] / D)
    prod = s["c1_definition"]["product"]
    for key, nm in names[1:3]:
        w["SensProduct%sPct" % nm] = _mpct(100.0 * prod[key]["numerator_sum"]
                                           / D)
    alli = s["c1b_all_intervals"]["rungs"]["C1b"]
    w["SensAllIntCOneBN"] = _mint(alli["numerator_union"])
    w["SensAllIntCOneBPct"] = _mpct(100.0 * alli["numerator_union"] / D)
    cal = s["c1b_calibration_fixed"]["rungs"]["C1b"]
    w["SensCalibCOneBN"] = _mint(cal["numerator_union"])
    w["SensCalibCOneBPct"] = _mpct(100.0 * cal["numerator_union"] / D)
    # round F1: stricter interval classifications, from the register
    civ = {r["partition"]: r for r in read_csv(CROWN_INTERVALS)}
    wide = sorted(p_ for p_, r in civ.items() if r["counts_for_c1b"] == "yes"
                  and r["calibration"] == "two_sided" and r["prior_width_ratio"]
                  and as_float(r["prior_width_ratio"]) >= C1B_WIDTH_RULE_RATIO)
    unch = sorted(p_ for p_, r in civ.items() if r["counts_for_c1b"] == "yes"
                  and r["calibration"] == "not_checked")
    for nm, drop in (("Width", set(wide)), ("Strict", set(wide) | set(unch))):
        num = sum(as_float(r["c1_units"]) for r in canon_e
                  if r["is_chronogram"] == "yes" and r["time_scale"] == "absolute"
                  and r["interval_recoverable"] == "yes"
                  and r["partition"] not in drop)
        w["Sens%sCOneBN" % nm] = _mint(num)
        w["Sens%sCOneBPct" % nm] = _mpct(100.0 * num / D)
    w["SensWidthNames"] = _names(wide)
    w["SensUncheckedNames"] = _names([x for x in unch if x not in wide])
    w["RuleWidthRatio"] = "%g" % C1B_WIDTH_RULE_RATIO
    w["RuleTipsOverNamedPct"] = "%.0f" % (100 * TABLE1_TIPS_OVER_NAMED)
    w["RuleElasticityPct"] = "%.0f" % (100 * ELASTICITY_FRACTION)
    for key, nm in names[:2]:
        n_ = rg[key]["numerator_union"]
        w["SensAllMinus%sPct" % nm] = _mpct(100.0 * n_ / (D * (1 - ELASTICITY_FRACTION)))
        w["SensAllPlus%sPct" % nm] = _mpct(100.0 * n_ / (D * (1 + ELASTICITY_FRACTION)))
    nob = s["leave_assembly_class_out"]["drop_automated_synthesis"]["rungs"]
    for key, nm in names[:2]:
        w["SensNoBirds%sPct" % nm] = _mpct(100.0 * nob[key]["numerator_union"]
                                           / D)
    ot = h["opentree_P01"]
    w["OTPhyloEukN"] = _mint(ot["ot_phylo_extant_species_eukaryota"])
    w["OTPhyloEukPct"] = _mpct(100.0 * ot["ot_phylo_extant_species_eukaryota"]
                               / D)
    w["OTCOneRatio"] = "%.2f" % (ot["atlas_C1_in_compared_partitions"]
                                 / ot["ot_phylo_extant_species_in_atlas_clades"])
    w["OTKeysBoth"] = _mint(ot["keys"]["keys_both"])
    me = ot["molecular_enumerated"]
    w["OTKeysMolecular"] = _mint(me["keys_atlas_molecular"])
    w["OTKeysMolecularAndOT"] = _mint(me["keys_molecular_and_ot"])
    w["OTAccMolecular"] = _mint(me["accepted_molecular"])
    w["OTAccMolecularAndOT"] = _mint(me["accepted_molecular_and_ot"])
    w["OTAccMolecularNotOT"] = _mint(me["accepted_molecular"]
                                     - me["accepted_molecular_and_ot"])
    w["OTAccMolecularAndOTPct"] = _mpct(100.0 * me["accepted_molecular_and_ot"]
                                        / me["accepted_molecular"])
    w["OTNPartitionsCompared"] = _mint(len(me["partitions"]))
    crows = read_csv(CROWN_AGES_CSV)
    w["CrownNChron"] = _mint(len(crows))
    w["CrownNCounted"] = _mint(sum(r["counts_for_c1b"] == "yes" for r in crows))
    cn = sorted(r["partition"] for r in crows if r["counts_for_c1b"] == "yes")
    w["CrownNamesCounted"] = ", ".join(cn[:-1]) + " and " + cn[-1]
    excl = Counter(r["c1b_exclusion"] for r in crows if r["c1b_exclusion"])
    for k in sorted(excl):
        w["CrownExcl" + crown_macro_name(k)] = _mint(excl[k])
    w["CrownNWithInterval"] = _mint(sum(bool(r["interval_kind"]) for r in crows))
    w["CrownNNoInterval"] = _mint(sum(not r["interval_kind"] for r in crows))
    cset = [r["partition"] for r in crows if r["counts_for_c1b"] == "yes"
            and r["interval_kind"] in CROWN_SET_KINDS]
    csum = [r["partition"] for r in crows if r["counts_for_c1b"] == "yes"
            and r["interval_kind"] not in CROWN_SET_KINDS]
    oset = [r["partition"] for r in crows if r["counts_for_c1b"] != "yes"
            and r["interval_kind"] in CROWN_SET_KINDS]
    for nm, ps in (("CountedSet", cset), ("CountedSummary", csum),
                   ("SetNotCounted", oset)):
        w["CrownN" + nm] = _mint(len(ps))
        w["CrownNames" + nm.replace("NotCounted", "NotCounted")] = _names(ps)
    # declared rule parameters, read from the script that applies them
    spec = importlib.util.spec_from_file_location(
        "crown_intervals_rules", P("scripts", "crown_intervals.py"))
    ci = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ci)
    w["RuleChronCV"] = "%g" % CHRONOGRAM_TOLERANCE
    w["RuleDegenerateMyr"] = "%g" % ci.DEGENERATE_WIDTH_MYR
    w["RuleBoundMyr"] = "%g" % ci.BOUND_TOL_MYR
    w["RulePriorContainPct"] = "%g" % (100 * ci.PRIOR_CONTAINMENT_MIN)
    w["RulePriorWidthRatio"] = "%g" % ci.PRIOR_WIDTH_RATIO_MIN
    w["RuleIntervalPct"] = "%g" % round(100 * (ci.QUANTILES[2] - ci.QUANTILES[0]))
    # C1 lost before C1b, by reason, on the union convention (R2-15): the
    # losses are set differences of the unions of consecutive levels; a lost
    # species is charged to the reason of the first group, in name order,
    # whose tree holds it at the level it reached.
    lost, levels_by_cat = union_losses(ctx, canon_e, crows)
    for key, nm in (("C0", "C0"), ("C1", "C1"), ("C1d", "C1d"), ("C1b", "C1b")):
        R.eq("union level %s = headline union" % key,
             sum(levels_by_cat[key].values()), rg[key]["numerator_union"])
    R.eq("C1 = C1b + every loss (union)",
         sum(lost[k] for k in ("NotChron", "RelTime", "NoInterval", "Kind",
                               "Calib")) + rg["C1b"]["numerator_union"],
         rg["C1"]["numerator_union"])
    ctx.union_losses = (lost, levels_by_cat)
    for k in ("NotChron", "RelTime", "NoInterval", "Kind", "Calib"):
        w["Lost%sN" % k] = _mint(lost[k])
    und = lost["NotChron"] + lost["RelTime"]
    w["LostUndatedOfCOnePct"] = _mpct(100.0 * und / rg["C1"]["numerator_union"])
    w["LostRelTimeOfUndatedPct"] = _mpct(100.0 * lost["RelTime"] / und)
    for k in ("NoInterval", "Kind", "Calib"):
        w["Lost%sOfCOneDPct" % k] = _mpct(100.0 * lost[k]
                                          / rg["C1d"]["numerator_union"])
    # the largest groups without a reference tree, named (R2-24)
    unrep = sorted((r for r in euk if not as_int(r["c0_units"], 0)),
                   key=lambda r: (-as_int(r["described"], 0), r["partition"]))
    big = [r["partition"] for r in unrep[:UNREP_NAMED]]
    w["UnrepLargestNames"] = ", ".join(big[:-1]) + " and " + big[-1]
    # Open Tree: share of the molecular species absent from its
    # phylogeny-placed set, by group (R2-10)
    me_ = h["opentree_P01"]["molecular_enumerated"]
    per_ = h["opentree_P01"]["per_partition"]
    R.eq("Open Tree absent-species shares sum to the total",
         sum(per_[p]["accepted_molecular_not_in_ot"] for p in me_["partitions"]),
         me_["accepted_molecular_not_in_ot"])
    for p in me_["partitions"]:
        w["OTNotOTSharePct" + crown_macro_name(p)] = _mpct(
            100.0 * per_[p]["accepted_molecular_not_in_ot"]
            / me_["accepted_molecular_not_in_ot"])
    # crown-age scaling of the fitted trees (R2-30), recomputed from ages
    sel = read_csv(WORKED_SEL)
    sf = []
    for r in sel:
        sf.append(as_float(r["crown_age_ma"]) / as_float(r["geometry_root_ma"]))
        R.eq("worked example %s scale factor" % r["partition"],
             as_float(r["time_scale_factor"]), sf[-1], atol=1e-9)
    w["WEScaleMin"] = "%.2f" % min(sf)
    w["WEScaleMax"] = "%.2f" % max(sf)
    we = h["worked_example_S13"]
    fits = we["fits"]

    eff = we["rho_effect_same_tree"]
    w["WENEligible"] = _mint(len(we["eligible"]))
    el = sorted(we["eligible"])
    w["WENamesEligible"] = ", ".join(el[:-1]) + " and " + el[-1]
    w["WENIdentity"] = _mint(sum(v["classification"] == "identity"
                                 for v in eff.values()))
    w["WENBoundary"] = _mint(sum(v["classification"] == "boundary"
                                 for v in eff.values()))
    w["WENRejected"] = _mint(sum(v["lr_statistic"] > WORKED_LR_CUTOFF
                                 for v in eff.values()))
    # D30: family and uniform fits side by side, no likelihood comparison
    epsd = {p: (as_float(fits[p + "|molecular|family"]["epsilon"]),
                as_float(fits[p + "|molecular|tips"]["epsilon"]))
            for p in we["eligible"]}
    w["WENFamilyEpsLower"] = _mint(sum(a < b for a, b in epsd.values()))
    w["WENFamilyEpsHigher"] = _mint(sum(a > b for a, b in epsd.values()))
    R.bump()
    if "loglik_family_minus_uniform" in json.dumps(we):
        R.fail("S13: the headline still carries a family-minus-uniform "
               "log-likelihood difference, which is not a model comparison "
               "(D30)")
    w["RuleLRCritical"] = "%.2f" % WORKED_LR_CUTOFF
    # the cutoff is the critical value of the equal chi2(0):chi2(1) mixture at
    # the stated level: P(chi2_1 > c) = 2 * level
    import statistics as _st
    R.eq("LR cutoff = critical value of the mixture at the test level",
         round(WORKED_LR_CUTOFF, 5),
         round(_st.NormalDist().inv_cdf(1.0 - WORKED_TEST_LEVEL) ** 2, 5))
    w["RuleTestLevelPct"] = "%g" % (100 * WORKED_TEST_LEVEL)
    rf = {p: (as_float(fits[p + "|full|tips"]["r"]),
              as_float(fits[p + "|molecular|tips"]["r"])) for p in we["eligible"]}
    w["WENFullGreater"] = _mint(sum(a > b for a, b in rf.values()))
    w["WENFullSmaller"] = _mint(sum(a < b for a, b in rf.values()))
    w["WENFullOverlap"] = _mint(sum(
        as_float(fits[p + "|full|tips"]["r_lo"]) <=
        as_float(fits[p + "|molecular|tips"]["r_hi"]) and
        as_float(fits[p + "|molecular|tips"]["r_lo"]) <=
        as_float(fits[p + "|full|tips"]["r_hi"]) for p in we["eligible"]))
    fam = we["family_sampling_fit"]["per_partition"]
    w["WENFamilyOverlap"] = _mint(sum(
        as_float(fits[p + "|molecular|family"]["epsilon_lo"]) <=
        as_float(fits[p + "|molecular|tips"]["epsilon_hi"]) and
        as_float(fits[p + "|molecular|tips"]["epsilon_lo"]) <=
        as_float(fits[p + "|molecular|family"]["epsilon_hi"])
        for p in we["eligible"]))
    for p in we["eligible"]:
        n = crown_macro_name(p)
        ft, fc = fits[p + "|molecular|tips"], fits[p + "|molecular|full_misapplied"]
        ff, fm = fits[p + "|full|tips"], fits[p + "|molecular|family"]
        w["WERhoTips" + n] = "%.2f" % as_float(ft["rho"])
        w["WERhoFull" + n] = "%.2f" % as_float(fc["rho"])
        w["WEEpsTips" + n] = "%.2f" % as_float(ft["epsilon"])
        w["WEEpsFull" + n] = "%.2f" % as_float(fc["epsilon"])
        w["WELambdaRatio" + n] = "%.2f" % (as_float(fc["lambda"])
                                           / as_float(ft["lambda"]))
        w["WELR" + n] = "%.2f" % max(0.0, 2.0 * (as_float(ft["loglik"])
                                                 - as_float(fc["loglik"])))
        w["WERFull" + n] = "%.3f" % as_float(ff["r"])
        w["WERMol" + n] = "%.3f" % as_float(ft["r"])
        w["WERRatio" + n] = "%.2f" % (as_float(ff["r"]) / as_float(ft["r"]))
        w["WEEpsFamily" + n] = "%.2f" % as_float(fm["epsilon"])
        w["WERFamily" + n] = "%.3f" % as_float(fm["r"])
    # round F1 (F1-14): polytomies and the fits with them resolved
    w.update(_worked_poly_macros(we["eligible"], fits))
    for r in canon_e:
        n = crown_macro_name(r["partition"])
        d = as_int(r["described"], 0)
        c0, c1 = as_int(r["c0_units"], 0), as_float(r["c1_units"])
        w["Described" + n] = _mint(d)
        w["CZeroN" + n] = _mint(c0)
        w["CZeroPct" + n] = _mpct(100.0 * c0 / d)
        w["COneN" + n] = _mint(c1)
        w["COnePct" + n] = _mpct(100.0 * c1 / d)
        w["FMol" + n] = "%.2f" % as_float(r["f_mol"])
    got = OrderedDict(_macros_in(NUMBER_MACROS))
    R.eq("numbers_macros.tex present", os.path.exists(NUMBER_MACROS), True)
    R.eq("numbers_macros.tex macro names", sorted(got), sorted(w))
    for k, v in w.items():
        if k in got:
            R.eq("\\%s" % k, got[k], v)
    seen = Counter(k for f in GENERATED_MACRO_FILES for k, _ in _macros_in(f))
    for k, c in sorted(seen.items()):
        if c > 1:
            R.bump()
            R.fail("macro \\%s is defined %d times across the generated "
                   "macro files" % (k, c))


# --------------------------------------------------------------------------
# TABLES: the three tables and the Figure 1 projection of the revised text
# --------------------------------------------------------------------------
TABLE_WHICH = P("Paper", "inputs", "table_which_tree.csv")
TABLE_SENS = P("Paper", "inputs", "table_rungs_sensitivity.csv")
TABLE_WORKED = P("Paper", "inputs", "table_worked_example.csv")
FIG_CASCADE = P("Paper", "inputs", "fig1_coverage_cascade.csv")
MOLECULAR_TIPS_CSV = P("Data", "molecular_tips.csv")


def _rho(x):
    return "%.2f" % x if x >= 0.1 else "%#.2g" % x


def check_tables(ctx):
    """Table 1 (which tree for which group), Table 2 (four levels and their
    sensitivity), Table 3 (worked example) and the Figure 1 cascade, each
    re-derived here from the basis table, the crown-age register, the
    molecular tip list and the fits, and compared cell by cell with the
    generated CSV that the .tex twin is rendered from."""
    basis = ctx.basis_rows
    canon_e = [r for r in basis if r["is_prokaryotic"] == "no"
               and as_int(r["c0_units"], 0) > 0]
    crown = {r["partition"]: r for r in read_csv(CROWN_AGES_CSV)}
    nmol = Counter(r["partition"] for r in read_csv(MOLECULAR_TIPS_CSV))
    reg = {r["filename"]: r for r in read_csv(DATASET_REGISTRY)}
    got = {r["partition"]: r for r in read_csv(TABLE_WHICH)}
    R.eq("Table 1 lists exactly the eukaryotic groups with a tree",
         sorted(got), sorted(r["partition"] for r in canon_e))
    for b in canon_e:
        p, g = b["partition"], got.get(b["partition"])
        if g is None:
            continue
        d, tips = as_int(b["described"]), as_int(b["labels_unique"])
        mt = nmol[p] if b["c1_basis"] == "enumerated" else tips
        c = crown.get(p)
        R.eq("Table 1 %s source" % p, g["bibkey"],
             reg[b["canonical_tree"]]["bibkey"])
        R.eq("Table 1 %s tips" % p, as_int(g["tips"]), tips)
        R.eq("Table 1 %s named %%" % p, g["named_pct"],
             _mpct(100.0 * as_int(b["c0_units"]) / d))
        R.eq("Table 1 %s molecular tips" % p, as_int(g["molecular_tips"]), mt)
        R.eq("Table 1 %s molecular share of tips" % p, g["molecular_share"],
             "%.2f" % (mt / tips))
        R.eq("Table 1 %s molecular %%" % p, g["molecular_pct"],
             _mpct(100.0 * as_float(b["c1_units"]) / d))
        R.eq("Table 1 %s time scale" % p, g["time_scale"], b["time_scale"])
        R.eq("Table 1 %s chronogram has a crown age" % p,
             bool(g["crown_age_ma"]), b["is_chronogram"] == "yes")
        if c is not None and b["is_chronogram"] == "yes":
            R.eq("Table 1 %s crown age" % p, g["crown_age_ma"],
                 "%.1f" % as_float(c["crown_age_ma"]))
            R.eq("Table 1 %s interval kind" % p, g["interval_kind"],
                 c["interval_kind"])
            want = ("none" if not c["interval_kind"] else "counted"
                    if c["counts_for_c1b"] == "yes" else c["c1b_exclusion"])
            R.eq("Table 1 %s interval status" % p, g["interval_status"], want)
        R.eq("Table 1 %s described" % p, as_int(g["described"]), d)
        sdv = {r["partition"]: (r.get("study_dating") or "").strip()
               for r in read_csv(DATING_DEPOSIT_CHECKS)}
        R.eq("Table 1 %s study dating" % p, g["study_dating"], sdv.get(p) or (
            "absolute_deposited" if b["is_chronogram"] == "yes"
            and b["time_scale"] == "absolute" else "relative_only"
            if b["time_scale"] == "relative" else None))
        R.eq("Table 1 %s study dating value known" % p,
             g["study_dating"] in STUDY_DATING_VALUES, True)
        if b["is_chronogram"] == "yes" and b["time_scale"] == "absolute":
            R.eq("Table 1 %s absolute chronogram -> study deposited it" % p,
                 g["study_dating"], "absolute_deposited")
        elif b["time_scale"] == "relative":
            R.eq("Table 1 %s relative time -> relative_only" % p,
                 g["study_dating"], "relative_only")
        else:
            R.eq("Table 1 %s not a chronogram -> not absolute_deposited" % p,
                 g["study_dating"] != "absolute_deposited", True)
        mo = {r["partition"] for r in read_csv(MOLECULAR_ONLY_TREES)}
        R.eq("Table 1 %s molecular-only tree" % p, g["molecular_only_tree"],
             "yes" if p in mo else "no")
        R.eq("Table 1 %s rho full" % p, g["rho_full"], _rho(tips / d))
        R.eq("Table 1 %s rho molecular" % p, g["rho_molecular"], _rho(mt / d))
        R.eq("Table 1 %s tips exceed named species by the flag fraction" % p,
             g["tips_exceed_named"],
             "yes" if tips > (1 + TABLE1_TIPS_OVER_NAMED) * as_int(b["c0_units"])
             else "no")
    # the molecular tip counts behind rho equal the worked-example fits
    fits = ctx.head["worked_example_S13"]["fits"]
    for p in ctx.head["worked_example_S13"]["eligible"]:
        R.eq("Table 1 %s molecular tips = worked-example molecular tree" % p,
             as_int(got[p]["molecular_tips"]),
             as_int(fits[p + "|molecular|tips"]["n_tips"]))
        R.eq("Table 1 %s tips = worked-example full tree" % p,
             as_int(got[p]["tips"]), as_int(fits[p + "|full|tips"]["n_tips"]))

    # Table 2: percentages formed here from numerator and denominator
    h, D = ctx.head, sum(as_int(r["described"], 0) for r in basis
                         if r["is_prokaryotic"] == "no")
    Drep = sum(as_int(r["described"], 0) for r in canon_e)
    s = h["sensitivities"]
    col = s["denominator"]["col_single_authority"]["filled"]
    nob = s["leave_assembly_class_out"]["drop_automated_synthesis"]["rungs"]
    u = h["union_bound_S09"]
    rows = {r["level"]: r for r in read_csv(TABLE_SENS)}
    for key, lvl in (("C0", "named"), ("C1", "molecular"), ("C1d", "dated"),
                     ("C1b", "interval")):
        r, e = rows[lvl], h["rungs"][key]
        n = e["numerator_union"]
        R.eq("Table 2 %s species" % lvl, r["species"].replace("{,}", ""),
             str(n))
        R.eq("Table 2 %s %%" % lvl, r["pct"], _mpct(100.0 * n / D))
        R.eq("Table 2 %s represented only" % lvl, r["rep"],
             _mpct(100.0 * n / Drep))
        cr = col["rungs"][key]
        R.eq("Table 2 %s CoL" % lvl, r["col"],
             _mpct(100.0 * cr.get("numerator_union", cr["numerator"])
                   / col["denominator"]))
        nb = nob[key]
        R.eq("Table 2 %s Birds removed" % lvl, r["nobirds"],
             _mpct(100.0 * nb.get("numerator_union", nb["numerator"]) / D))
    R.eq("Table 2 named union", rows["named"]["union"],
         _mpct(100.0 * u["C0_union"] / D))
    R.eq("Table 2 molecular union bound", rows["molecular"]["union"],
         "<= " + _mpct(100.0 * u["C1_union_high"] / D))
    R.eq("Table 2 named GBIF", rows["named"]["gbif"],
         _mpct(100.0 * s["c0_gbif_unmatched"]["numerator"] / D))
    for key, lvl in (("C1", "molecular"), ("C1d", "dated")):
        R.eq("Table 2 %s product" % lvl, rows[lvl]["product"],
             _mpct(100.0 * s["c1_definition"]["product"][key]["numerator_sum"]
                   / D))
    cal = s["c1b_calibration_fixed"]["rungs"]["C1b"]
    R.eq("Table 2 interval calibration-fixed counted", rows["interval"]["calib"],
         _mpct(100.0 * cal["numerator_union"] / D))
    alli = s["c1b_all_intervals"]["rungs"]["C1b"]
    R.eq("Table 2 interval every archived interval counted",
         rows["interval"]["allint"], _mpct(100.0 * alli["numerator_union"] / D))

    # Table 3: from the fits
    we = h["worked_example_S13"]
    t3 = {r["partition"]: r for r in read_csv(TABLE_WORKED)}
    R.eq("Table 3 lists every eligible group", sorted(t3), sorted(we["eligible"]))
    for p in we["eligible"]:
        r = t3[p]
        ft, fc = fits[p + "|molecular|tips"], fits[p + "|molecular|full_misapplied"]
        ff, fm = fits[p + "|full|tips"], fits[p + "|molecular|family"]
        for col_, f, k, fmt in (("rho_molecular", ft, "rho", "%.2f"),
                                ("rho_fulltree", fc, "rho", "%.2f"),
                                ("eps_molecular", ft, "epsilon", "%.2f"),
                                ("eps_fullrho", fc, "epsilon", "%.2f"),
                                ("r_full", ff, "r", "%.3f"),
                                ("r_molecular", ft, "r", "%.3f"),
                                ("eps_family", fm, "epsilon", "%.2f")):
            R.eq("Table 3 %s %s" % (p, col_), r[col_], fmt % as_float(f[k]))
        R.eq("Table 3 %s lambda ratio" % p, r["lambda_ratio"],
             "%.2f" % (as_float(fc["lambda"]) / as_float(ft["lambda"])))
        lr = max(0.0, 2.0 * (as_float(ft["loglik"]) - as_float(fc["loglik"])))
        R.eq("Table 3 %s LR" % p, r["lr"], "%.2f" % lr)
        for col_, k in (("r_family", "r"), ("r_family_lo", "r_lo"),
                        ("r_family_hi", "r_hi"), ("eps_family_lo", "epsilon_lo"),
                        ("eps_family_hi", "epsilon_hi")):
            R.eq("Table 3 %s %s" % (p, col_), r[col_],
                 ("%.3f" if col_.startswith("r_") else "%.2f")
                 % as_float(fm[k]))
        R.eq("Table 3 %s has no log-likelihood comparison column" % p,
             "dll_family" in r, False)
        R.eq("Table 3 %s rejected" % p, r["classification"].endswith(
            "(rejected)"), lr > WORKED_LR_CUTOFF)

    # Figure 1 cascade: per-category levels and the losses, on the union
    # convention, cell by cell against union_losses() (R2-15)
    cas = read_csv(FIG_CASCADE)
    ul, lc = getattr(ctx, "union_losses", (None, None))
    R.eq("Figure 1 union losses computed", ul is not None, True)
    if ul is None:
        return
    for r in cas:
        if r["kind"] == "level":
            R.eq("Figure 1 %s %s" % (r["level"], r["group"]),
                 int(r["species"]), lc[r["level"]][r["group"]])
        else:
            R.eq("Figure 1 loss %s" % r["group"], int(r["species"]),
                 ul[r["group"]])
    lv = defaultdict(int)
    for r in cas:
        if r["kind"] == "level":
            lv[r["level"]] += int(r["species"])
    for key in ("C0", "C1", "C1d", "C1b"):
        R.eq("Figure 1 %s total = union" % key, lv[key],
             h["rungs"][key]["numerator_union"])
    loss = defaultdict(int)
    for r in cas:
        if r["kind"] == "loss":
            loss[r["level"]] += int(r["species"])
    for a, b, step in (("C0", "C1", "C0>C1"), ("C1", "C1d", "C1>C1d"),
                       ("C1d", "C1b", "C1d>C1b")):
        R.eq("Figure 1 losses %s" % step, loss[step],
             h["rungs"][a]["numerator_union"] - h["rungs"][b]["numerator_union"])


# --------------------------------------------------------------------------
# SUPPLEMENT: the Supporting Information, its numbering and its tables
# --------------------------------------------------------------------------
SUPP_BODY = P("Paper", "shared", "supplementary_body_v2.tex")
SUPP_META = P("Paper", "shared", "supplementary_metadata_v2.tex")
SUPP_MACROS = P("Paper", "inputs", "supp_macros.tex")
SUPP_IN = lambda name: P("Paper", "inputs", name)          # noqa: E731
REFERENCE_TREE_SOURCES = P("Data", "reference_tree_sources.csv")
CANONICAL_SOURCE_PREP = P("Data", "canonical_source_prep.csv")
SUPP_SELECTION = {"as_archived", "release", "study_chronogram",
                  "more_branch_lengths", "first_of_set", "mcc_of_set"}
# the decision a prepared source cites -> the selection it implies
SUPP_PREP_DECISION = {"D8": "study_chronogram", "D15": "more_branch_lengths",
                      "D16": "first_of_set", "D17": "mcc_of_set"}
SUPP_EDGE_RTOL = 1e-9
SUPP_WRAPPERS = [P("Paper", "MEE", "supplementary.tex"),
                 P("Paper", "MEE", "supplementary_anon.tex"),
                 P("Paper", "bioRxiv", "supplementary.tex")]


def _supp_first_tree_text(path, index=1):
    """The index-th tree of a Newick or NEXUS file, comments removed."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    text = re.sub(r"\[[^\]]*\]", "", text)
    if text.lstrip().upper().startswith("#NEXUS"):
        trees = re.findall(r"(?im)^\s*tree\s+[^=]+=\s*(\(.*?;)", text, re.S)
    else:
        trees = [t + ";" for t in text.split(";") if "(" in t]
    if len(trees) < index:
        return None
    return trees[index - 1]


def _supp_edges(text):
    return sorted(float(x) for x in
                  re.findall(r":\s*([-+0-9.eE]+|nan|NaN)", text)
                  if x.lower() != "nan")


def _supp_csv(name):
    return read_csv(SUPP_IN(name + ".csv"))


def _supp_cmp_table(name, want_rows, keys):
    got = _supp_csv(name)
    R.eq("Supplement %s rows" % name, [r[keys[0]] for r in got],
         [str(w[keys[0]]) for w in want_rows])
    gotk = {r[keys[0]]: r for r in got}
    for w in want_rows:
        g = gotk.get(str(w[keys[0]]))
        if g is None:
            continue
        for k in keys[1:]:
            R.eq("Supplement %s %s %s" % (name, w[keys[0]], k), g.get(k),
                 str(w[k]))


def _supp_mpct(x):
    return "--" if x is None else ("%.1f" % x if abs(x) >= 1.0 else "%.2f" % x)


def check_supplement(ctx):
    """The Supporting Information: (1) it has exactly one section per
    Supplement number the main text cites, numbered in order of first
    citation, and the wrappers read the revised body; (2) every generated
    table it inputs exists, and each supp_*.csv is re-derived here, cell by
    cell, from the registers, the shipped trees and the source files; (3) the
    reference trees' archived sources reproduce the branch lengths of the
    source files the register names; (4) its macros are re-derived."""
    meta = open(MANUSCRIPT_META, encoding="utf-8").read()
    body = strip_tex_comments(open(MANUSCRIPT_BODY, encoding="utf-8").read())
    ids = re.findall(r"^\\newcommand\{\\(IdSupp\w+)\}\{S(\d+)\}", meta, re.M)
    nums = [int(n) for _k, n in ids]
    R.eq("Supplement numbers are S1..Sn in definition order", nums,
         list(range(1, len(ids) + 1)))
    text = body + "\n" + strip_tex_comments(meta)
    first = {}
    for k, _n in ids:
        m = re.search(r"\\%s(?![A-Za-z])" % k,
                      re.sub(r"^\\newcommand\{\\IdSupp\w+\}\{S\d+\}\s*$", "",
                             text, flags=re.M))
        R.eq("Supplement %s is cited in the article" % k, m is not None, True)
        if m is not None:
            # first citation: the body counts in reading order; a citation
            # only in a caption or statement counts after the body
            mb = re.search(r"\\%s(?![A-Za-z])" % k, body)
            first[k] = mb.start() if mb else len(body) + m.start()
    order = [k for k, _n in sorted(first.items(), key=lambda kv: kv[1])]
    R.eq("Supplement numbers follow the order of first citation", order,
         [k for k, _n in ids])
    sb = strip_tex_comments(open(SUPP_BODY, encoding="utf-8").read())
    secs = re.findall(r"\\section\{Supplement \\(IdSupp\w+)", sb)
    R.eq("Supplement sections = the cited numbers, in order", secs,
         [k for k, _n in ids])
    for w in SUPP_WRAPPERS:
        t = strip_tex_comments(open(w, encoding="utf-8").read())
        R.eq("%s inputs the revised supplement body" % rel(w),
             "supplementary_body_v2" in t and "supplementary_body}" not in t
             and "supplementary_body.tex" not in t, True)
    inputs = re.findall(r"\\inputtable\{\.\./inputs/([\w.]+)\}", sb)
    for f in inputs:
        R.eq("Supplement input %s exists" % f, os.path.exists(SUPP_IN(f)), True)
    for f in sorted(os.listdir(P("Paper", "inputs"))):
        if f.startswith("supp_") and f.endswith(".tex") and f != "supp_macros.tex":
            R.eq("generated %s is used by the Supplement" % f, f in inputs, True)
    # every macro the Supplement uses is defined
    defined = set()
    for f in GENERATED_MACRO_FILES + (MANUSCRIPT_META, SUPP_META, SUPP_MACROS):
        if os.path.exists(f):
            defined |= set(re.findall(r"\\(?:newcommand|providecommand)"
                                      r"\{\\([A-Za-z]+)\}", open(f).read()))
    for name in sorted(set(re.findall(r"\\((?:Supp|Id|C[A-Z]|Sens|WE|OT|"
                                      r"Crown|Rule|Num|Described|DT|Labels|"
                                      r"Unresolved)[A-Za-z]*)", sb))):
        R.eq("Supplement macro \\%s is defined" % name, name in defined, True)

    basis = ctx.basis_rows
    canon_e = [r for r in basis if r["is_prokaryotic"] == "no"
               and as_int(r["c0_units"], 0) > 0]
    h = ctx.head

    # ---- reference trees: register and source files -------------------
    reg = {r["partition"]: r for r in read_csv(REFERENCE_TREE_SOURCES)}
    prep = {r["filename"]: r for r in read_csv(CANONICAL_SOURCE_PREP)}
    dsreg = {r["filename"]: r for r in read_csv(DATASET_REGISTRY)}
    R.eq("reference_tree_sources.csv covers the reference trees",
         sorted(reg), sorted(r["partition"] for r in canon_e))
    for b in canon_e:
        p, fn = b["partition"], b["canonical_tree"]
        r = reg.get(p)
        if r is None:
            continue
        R.eq("%s register filename" % p, r["filename"], fn)
        R.eq("%s selection known" % p, r["selection"] in SUPP_SELECTION, True)
        R.eq("%s register evidence present" % p,
             bool((r.get("evidence") or "").strip()), True)
        if fn in prep:
            R.eq("%s source = prepared source" % p, r["source_file"],
                 prep[fn]["raw_source"])
            dec = re.match(r"\s*(D\d+)", prep[fn]["evidence"])
            dec = dec.group(1) if dec else None
            if dec == "D6":   # 'D6/D16'
                dec = "D16"
            want = SUPP_PREP_DECISION.get(dec, "release")
            R.eq("%s selection follows the prepared source's rule" % p,
                 r["selection"], want)
            continue
        R.eq("%s is not a prepared source -> as archived or release" % p,
             r["selection"] in ("as_archived", "release"), True)
        sp = P("Data", *r["source_file"].split("/"))
        cp = P(*ctx.canon_sources.get(fn, "").split("/"))
        R.eq("%s source file exists" % p, os.path.exists(sp), True)
        if not (os.path.exists(sp) and os.path.exists(cp)):
            continue
        idx = as_int(r.get("tree_index"), 1) or 1
        a = _supp_edges(_supp_first_tree_text(sp, idx) or "")
        c = _supp_edges(_supp_first_tree_text(cp) or "")
        R.bump()
        if len(a) != len(c) or any(abs(x - y) > SUPP_EDGE_RTOL * max(1.0, abs(x))
                                   for x, y in zip(a, c)):
            R.fail("%s: the archived source %s does not reproduce the branch "
                   "lengths of %s (%d vs %d edges)"
                   % (p, rel(cp), r["source_file"], len(c), len(a)))
    order = lambda rows: sorted(rows, key=lambda x: (            # noqa: E731
        CATEGORY_ORDER.index(x["category"]), x["partition"]))
    want = order([OrderedDict([
        ("partition", b["partition"]), ("category", b["category"]),
        ("filename", b["canonical_tree"]),
        ("bibkey", dsreg[b["canonical_tree"]]["bibkey"]),
        ("source_file", os.path.basename(reg[b["partition"]]["source_file"])),
        ("selection", reg[b["partition"]]["selection"]),
        ("tree_used", reg[b["partition"]]["tree_used"])])
        for b in canon_e if b["partition"] in reg])
    _supp_cmp_table("supp_reference_trees", want,
                    ["partition", "filename", "bibkey", "source_file",
                     "selection", "tree_used"])

    # ---- removed and merged tips ---------------------------------------
    rows = []
    for b in canon_e:
        p, fn = b["partition"], b["canonical_tree"]
        src = ctx.source_leaves.get(fn)
        tips = as_int(b["labels_unique"])
        nokey = sum(1 for nm, sid in ctx.idmap[p].items()
                    if d4_species_key(nm, sid) is None)
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("source_tips", src), ("removed_outside", ctx.n_pruned.get(fn, 0)),
            ("merged_same_species", None if src is None
             else src - ctx.n_pruned.get(fn, 0) - tips),
            ("tips", tips), ("unresolved_labels_kept", nokey)]))
    _supp_cmp_table("supp_pruned", order(rows),
                    ["partition", "source_tips", "removed_outside",
                     "merged_same_species", "tips", "unresolved_labels_kept"])

    # ---- name resolution, from col_resolution.csv ----------------------
    gb = h["sensitivities"]["c0_gbif_unmatched"]["per_partition"]
    rows = []
    for b in canon_e:
        p = b["partition"]
        cr = ctx.col_rows.get(p, {})
        mt = Counter(r["match_type"] for r in cr.values())
        res = [r for r in cr.values() if r["match_type"] in
               COL_RESOLUTION_RESOLVED and r["col_id"]]
        acc = sum(n for k, n in mt.items() if k.startswith(("accepted",
                                                             "provisional")))
        syn = sum(n for k, n in mt.items() if k.startswith("synonym"))
        n = len(cr)
        named = len({r["col_id"] for r in res})
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]), ("labels", n),
            ("accepted", acc), ("synonym", syn), ("unresolved", n - len(res)),
            ("unresolved_pct", "%.1f" % (100.0 * (n - len(res)) / n)),
            ("named_species", named), ("merged", len(res) - named),
            ("gbif_added", (gb.get(p) or {}).get("gbif_species_added", 0)),
            ("override", sum(1 for r in res
                             if r.get("query_source") == "override"))]))
        R.eq("%s named species = c0_units" % p, named, as_int(b["c0_units"]))
    _supp_cmp_table("supp_names", order(rows),
                    ["partition", "labels", "accepted", "synonym", "unresolved",
                     "unresolved_pct", "named_species", "merged", "gbif_added",
                     "override"])
    ovr_ids = {r["standardized_id"] for r in read_csv(NAME_OVERRIDES)}
    R.eq("every override query comes from a row of name_overrides.csv",
         all(sid in ovr_ids for rows_ in ctx.col_rows.values()
             for sid, r in rows_.items() if r.get("query_source") == "override"),
         True)

    # ---- molecular tips --------------------------------------------------
    nmol = Counter(r["partition"] for r in read_csv(MOLECULAR_TIPS_CSV))
    ev = h["molecular_evidence"]
    mo = {r["partition"] for r in read_csv(MOLECULAR_ONLY_TREES)}
    rows = []
    for b in canon_e:
        p = b["partition"]
        tips = as_int(b["labels_unique"])
        mt = nmol[p] if b["c1_basis"] == "enumerated" else tips
        ex = ev[p].get("exact_label_intersection")
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("f_basis", b["f_basis"]),
            ("source_file", os.path.basename(b["molecular_source_file"] or "")),
            ("tips", tips), ("molecular_tips", mt),
            ("molecular_fraction", "%.3f" % (mt / tips)),
            ("exact_label_tips", "" if ex is None else int(ex)),
            ("molecular_species", int(round(as_float(b["c1_units"])))),
            ("product", "%.1f" % (as_int(b["c0_units"]) * as_float(b["f_mol"]))),
            ("molecular_only_tree", "yes" if p in mo else "no")]))
    _supp_cmp_table("supp_molecular", order(rows),
                    ["partition", "f_basis", "source_file", "tips",
                     "molecular_tips", "molecular_fraction", "exact_label_tips",
                     "molecular_species", "product", "molecular_only_tree"])

    # ---- time scale ------------------------------------------------------
    sdv = {r["partition"]: (r.get("study_dating") or "").strip()
           for r in read_csv(DATING_DEPOSIT_CHECKS)}
    rows = []
    for b in canon_e:
        p, fn = b["partition"], b["canonical_tree"]
        m = ctx.measured.get(fn, {})
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]),
            ("roottotip_cv", "%.4f" % m.get("cv", float("nan"))),
            ("ultrametric", "yes" if m.get("cv", 1.0) <= CHRONOGRAM_TOLERANCE
             else "no"),
            ("root_height", "%.4g" % m.get("root_height", float("nan"))),
            ("time_scale", b["time_scale"]),
            ("study_dating", sdv.get(p) or (
                "absolute_deposited" if b["is_chronogram"] == "yes"
                else "relative_only"))]))
    _supp_cmp_table("supp_time_scale", order(rows),
                    ["partition", "roottotip_cv", "ultrametric", "root_height",
                     "time_scale", "study_dating"])

    # ---- crown ages and intervals ----------------------------------------
    cat = {r["partition"]: r["category"] for r in basis}
    rows = []
    for c in read_csv(CROWN_AGES_CSV):
        p, kind = c["partition"], c["interval_kind"]
        iv = ctx.crown.get(p, {})
        st = ("none" if not kind else "counted"
              if c["counts_for_c1b"] == "yes" else c["c1b_exclusion"])
        f1 = lambda x: "%.1f" % as_float(x) if x else ""     # noqa: E731
        rows.append(OrderedDict([
            ("partition", p), ("category", cat[p]),
            ("crown_clade", c["crown_clade"]),
            ("crown_age_ma", f1(c["crown_age_ma"])),
            ("interval_lo_ma", f1(c["interval_lo_ma"]) if kind else ""),
            ("interval_hi_ma", f1(c["interval_hi_ma"]) if kind else ""),
            ("interval_kind", kind),
            ("calibration", c["calibration"] if kind else ""),
            ("prior_lo_ma", f1(iv.get("prior_lo_ma"))),
            ("prior_hi_ma", f1(iv.get("prior_hi_ma"))),
            ("prior_containment_pct", "%.0f" % (100 * as_float(
                iv["prior_containment"])) if iv.get("prior_containment")
             else ""),
            ("status", st)]))
        if st == "prior_reproduced" and iv.get("prior_containment") \
                and iv.get("kind") != "calibration_prior_sample":
            R.eq("%s prior rule: containment and width ratio" % p,
                 as_float(iv["prior_containment"]) >= PRIOR_CONTAINMENT_MIN
                 and as_float(iv["prior_width_ratio"]) >= PRIOR_WIDTH_RATIO_MIN,
                 True)
    _supp_cmp_table("supp_crown", order(rows),
                    ["partition", "crown_clade", "crown_age_ma",
                     "interval_lo_ma", "interval_hi_ma", "interval_kind",
                     "calibration", "prior_lo_ma", "prior_hi_ma",
                     "prior_containment_pct", "status"])

    # ---- sensitivity: the alternatives recomputed from the basis --------
    D = sum(as_int(r["described"], 0) for r in basis
            if r["is_prokaryotic"] == "no")
    crown_ok = {r["partition"]: r["counts_for_c1b"] == "yes"
                for r in read_csv(CROWN_AGES_CSV)}

    def levels(rows_, keep=lambda r: True):
        out = OrderedDict((k, 0.0) for k in RUNGS)
        for r in rows_:
            if not keep(r):
                continue
            out["C0"] += as_int(r["c0_units"], 0)
            c1 = as_float(r["c1_units"]) or 0.0
            out["C1"] += c1
            if r["is_chronogram"] == "yes" and r["time_scale"] == "absolute":
                out["C1d"] += c1
                if crown_ok.get(r["partition"]):
                    out["C1b"] += c1
        return out

    L = levels(canon_e)
    sens = {r["alternative"]: r for r in _supp_csv("supp_sensitivity")}

    def check_row(key, den, lv):
        r = sens.get(key)
        R.eq("Supplement S7 row %s present" % key, r is not None, True)
        if r is None:
            return
        R.eq("Supplement S7 %s denominator" % key, as_int(r["denominator"]),
             int(round(den)))
        for k in RUNGS:
            if lv.get(k) is not None:
                R.eq("Supplement S7 %s %s" % (key, k), r[k],
                     _supp_mpct(100.0 * lv[k] / den))

    check_row("reported", D, L)
    f = ELASTICITY_FRACTION
    ins = [r for r in basis if r["partition"] ==
           h["sensitivities"]["denominator"]["largest_single_count"]["partition"]]
    di = as_int(ins[0]["described"])
    check_row("insects_minus", D - f * di, L)
    check_row("insects_plus", D + f * di, L)
    check_row("all_minus", D * (1 - f), L)
    check_row("all_plus", D * (1 + f), L)
    col = {r["partition"]: r for r in read_csv(COL_DENOMINATORS)}
    dcol = sum(as_int(col[r["partition"]]["col_filled"], 0) or
               as_int(r["described"], 0) for r in basis
               if r["is_prokaryotic"] == "no")
    check_row("col_filled", dcol, L)
    check_row("rep", sum(as_int(r["described"], 0) for r in canon_e), L)
    dom = max(canon_e, key=lambda r: as_int(r["c0_units"], 0))
    check_row("no_seed", D - as_int(dom["described"]),
              levels(canon_e, lambda r: r is not dom))
    asm = {r["partition"]: r["topology_assembly"] for r in canon_e}
    for cls in ("automated_synthesis", "compiled", "primary_analysis"):
        check_row("drop_" + cls, D,
                  levels(canon_e, lambda r, c=cls: asm[r["partition"]] != c))
        check_row("retain_only_" + cls, D,
                  levels(canon_e, lambda r, c=cls: asm[r["partition"]] == c))
    u = h["union_bound_S09"]
    check_row("union", D, OrderedDict([("C0", u["C0_union"]),
                                       ("C1", u["C1_union_high"])]))
    civ = {r["partition"]: r for r in read_csv(CROWN_INTERVALS)}
    wide = {p_ for p_, r in civ.items() if r["counts_for_c1b"] == "yes"
            and r["calibration"] == "two_sided" and r["prior_width_ratio"]
            and as_float(r["prior_width_ratio"]) >= C1B_WIDTH_RULE_RATIO}
    unch = {p_ for p_, r in civ.items() if r["counts_for_c1b"] == "yes"
            and r["calibration"] == "not_checked"}
    for key, drop in (("width_rule", wide), ("strict", wide | unch)):
        lv = levels(canon_e, lambda r, d_=drop: r["partition"] not in d_)
        check_row(key, D, OrderedDict([("C1b", lv["C1b"])]))
        R.eq("Supplement S7 %s prints the interval level only" % key,
             [sens[key][k] for k in ("C0", "C1", "C1d")], ["", "", ""])

    # ---- denominators by group -------------------------------------------
    uni = u["per_partition"]
    rows = []
    for b in basis:
        if b["is_prokaryotic"] != "no":
            continue
        p, d = b["partition"], as_int(b["described"], 0)
        n0 = as_int(b["c0_units"], 0)
        cf = as_int(col[p]["col_filled"], 0)
        rows.append(OrderedDict([
            ("partition", p), ("category", b["category"]), ("described", d),
            ("col_count", cf),
            ("named_pct", _supp_mpct(100.0 * n0 / d) if n0 else ""),
            ("named_col_pct", _supp_mpct(100.0 * n0 / cf) if n0 and cf else ""),
            ("named_all_trees_pct", _supp_mpct(100.0 * uni[p]["C0_union"] / d)
             if p in uni else "")]))
    rows.sort(key=lambda x: (CATEGORY_ORDER.index(x["category"]),
                             -x["described"], x["partition"]))
    _supp_cmp_table("supp_denominators", rows,
                    ["partition", "described", "col_count", "named_pct",
                     "named_col_pct", "named_all_trees_pct"])

    # ---- Open Tree --------------------------------------------------------
    rows = []
    for r in read_csv(OT_CSV):
        if r["atlas_C1_tips_enumerated"] != "yes":
            continue
        am, both = as_int(r["accepted_molecular"]), as_int(
            r["accepted_molecular_and_ot"])
        rows.append(OrderedDict([
            ("partition", r["partition"]), ("category", cat[r["partition"]]),
            ("molecular_species", am), ("both", both),
            ("molecular_not_ot", am - both),
            ("both_pct", "%.1f" % (100.0 * both / am) if am else "")]))
    _supp_cmp_table("supp_opentree", order(rows),
                    ["partition", "molecular_species", "both",
                     "molecular_not_ot", "both_pct"])
    R.eq("Supplement S8 groups = groups with listed molecular tips",
         len(rows), sum(1 for b in canon_e if b["c1_basis"] == "enumerated"))

    # ---- worked example ----------------------------------------------------
    fits = read_csv(WORKED_FITS)
    got = _supp_csv("supp_worked_fits")
    R.eq("Supplement S9 one row per fit", len(got), len(fits))
    byk = {(r["partition"], r["tree"], r["rho_convention"]): r for r in fits}
    for g in got:
        fr = byk.get((g["partition"], g["tree"], g["rho_convention"]))
        R.eq("Supplement S9 fit %s/%s/%s exists" % (g["partition"], g["tree"],
                                                 g["rho_convention"]),
             fr is not None, True)
        if fr is None:
            continue
        for k, fmt in (("r", "%.4f"), ("r_lo", "%.4f"), ("r_hi", "%.4f"),
                       ("epsilon", "%.3f"), ("epsilon_lo", "%.3f"),
                       ("epsilon_hi", "%.3f"), ("loglik", "%.2f")):
            R.eq("Supplement S9 %s %s %s %s" % (g["partition"], g["tree"],
                                                 g["rho_convention"], k),
                 g[k], fmt % as_float(fr[k]))
        R.eq("Supplement S9 %s %s tips" % (g["partition"], g["tree"]),
             as_int(g["n_tips"]), as_int(fr["n_tips"]))
        R.eq("Supplement S9 %s %s %s events, at polytomies" % (
            g["partition"], g["tree"], g["rho_convention"]),
            (g["n_events"], g["n_polytomy_events"]),
            (fr["n_events"], fr["n_polytomy_events"]))
    _check_supp_polytomy(byk)
    fam = read_csv(WORKED_FAMILY_FRACTIONS)
    for g in _supp_csv("supp_worked_families"):
        fs = sorted(as_float(r["f"]) for r in fam
                    if r["partition"] == g["partition"])
        R.eq("Supplement S9 %s families" % g["partition"], as_int(g["families"]),
             len(fs))
        if fs:
            R.eq("Supplement S9 %s min/max fraction" % g["partition"],
                 (g["f_min"], g["f_max"]), ("%.2f" % fs[0], "%.2f" % fs[-1]))
        fr_ = [r for r in fam if r["partition"] == g["partition"]
               and r["n_col_species"]]
        R.eq("Supplement S9 %s pooled family fraction" % g["partition"],
             g["f_pooled"], _rho(sum(as_int(r["n_accepted_names"]) for r in fr_)
                                 / float(sum(as_int(r["n_col_species"])
                                             for r in fr_))))
        fo = {r["child"]: r for r in read_csv(WORKED_FORDER)
              if r["partition"] == g["partition"]}
        R.eq("Supplement S9 %s family fit, u from the first and the last child"
             % g["partition"],
             (g["r_first"], g["r_last"], g["epsilon_first"], g["epsilon_last"]),
             ("%.4f" % as_float(fo["first"]["r"]), "%.4f" % as_float(fo["last"]["r"]),
              "%.3f" % as_float(fo["first"]["epsilon"]),
              "%.3f" % as_float(fo["last"]["epsilon"])))
        R.eq("Supplement S9 %s family fit with u from the first child = the "
             "reported family fit" % g["partition"],
             (as_float(fo["first"]["r"]), as_float(fo["first"]["loglik"])),
             (as_float(byk[(g["partition"], "molecular", "family")]["r"]),
              as_float(byk[(g["partition"], "molecular", "family")]["loglik"])))
    for f in read_csv(WORKED_FITS):
        if f["rho_convention"] == "family":
            R.eq("Supplement S9 %s family likelihood = uniform at equal f, "
                 "within the stated tolerance" % f["partition"],
                 abs(as_float(f["family_uniform_check"])) <= WORKED_UNIFORM_TOL,
                 True)

    # ---- data and code ----------------------------------------------------
    R.eq("Supplement check list = verify's CHECKS",
         [(r["check"], r["title"]) for r in _supp_csv("supp_checks")],
         [tuple(c) for c in CHECKS])
    hdr = set(basis[0].keys())
    for r in _supp_csv("supp_terms"):
        for c in [x.strip() for x in r["columns"].split(";") if x.strip()]:
            R.eq("Supplement term column %s is in the basis" % c, c in hdr, True)
        if r["file"]:
            R.eq("Supplement term file %s is in the release" % r["file"],
                 release_path_known(r["file"]), True)

    # ---- macros -------------------------------------------------------------
    mac = OrderedDict(_macros_in(SUPP_MACROS))
    num = lambda k: int(mac[k].replace("{,}", ""))           # noqa: E731
    pr = _supp_csv("supp_pruned")
    for k, col_ in (("SuppRemovedN", "removed_outside"),
                    ("SuppMergedN", "merged_same_species"),
                    ("SuppOpenLabelsN", "unresolved_labels_kept"),
                    ("SuppSourceTipsN", "source_tips"), ("SuppTipsN", "tips")):
        R.eq("macro %s = sum of Table S3" % k, num(k),
             sum(as_int(r[col_]) for r in pr))
    R.eq("macro SuppTipsN = LabelsEukN", num("SuppTipsN"),
         sum(as_int(r["labels_unique"]) for r in canon_e))
    R.eq("macro SuppOverridesN", num("SuppOverridesN"),
         len(read_csv(NAME_OVERRIDES)))
    R.eq("macro SuppRejectedTargetsN", num("SuppRejectedTargetsN"),
         len(read_csv(P("Data", "standardized", "name_rejected_targets.csv"))))
    R.eq("macro SuppSupraGeneraN", num("SuppSupraGeneraN"),
         len(D4_SUPRA_GENERA))
    R.eq("macro SuppProfileCutoff", mac["SuppProfileCutoff"],
         "%.2f" % (WORKED_CI_CUTOFF / 2.0))
    R.eq("macro SuppVerifyChecksN", num("SuppVerifyChecksN"), len(CHECKS))
    R.eq("macro SuppUniformCheckTol", mac["SuppUniformCheckTol"],
         r"$10^{%d}$" % round(math.log10(WORKED_UNIFORM_TOL)))
    R.eq("macro SuppCrossCheckTolPct", mac["SuppCrossCheckTolPct"],
         "%g" % (100 * WORKED_CROSSCHECK_RTOL))
    rd = 0.0
    for r in read_csv(WORKED_CHECK):
        lam = max(as_float(r["this_lambda"]), as_float(r["phytools_lambda"]))
        for a_, b_ in (("this_lambda", "phytools_lambda"), ("this_mu", "phytools_mu")):
            x, y = as_float(r[a_]), as_float(r[b_])
            if lam > 1e-6:
                rd = max(rd, abs(x - y) / lam)
    R.eq("phytools cross-check within the stated tolerance",
         rd <= WORKED_CROSSCHECK_RTOL, True)
    sel = Counter(r["selection"] for r in reg.values())
    for code in SUPP_SELECTION:
        nm = "SuppSel" + "".join(w.capitalize() for w in code.split("_"))
        R.eq("macro %s" % nm, num(nm), sel.get(code, 0))
    ts = _supp_csv("supp_time_scale")
    ul = [as_float(r["roottotip_cv"]) for r in ts if r["ultrametric"] == "yes"]
    ot = [as_float(r["roottotip_cv"]) for r in ts if r["ultrametric"] == "no"]
    R.eq("macro SuppMaxCVUltrametric", mac["SuppMaxCVUltrametric"],
         "%.4f" % max(ul))
    R.eq("macro SuppMinCVOther", mac["SuppMinCVOther"], "%.4f" % min(ot))
    R.eq("ultrametric and other CVs are separated by the rule",
         max(ul) <= CHRONOGRAM_TOLERANCE < min(ot), True)
    shipped = {f for f in os.listdir(STD_TREES) if f.endswith(".nwk")}
    rest = [r for r in read_csv(P("Data", "_archive", "2026-10-revision",
                                  "corrupted_branch_lengths", "RESTORED.csv"))
            if r["file"].startswith("Data/standardized/trees/")
            and r["changed"] == "yes"
            and os.path.basename(r["file"]) in shipped]
    R.eq("macro SuppRestoredTreesN", num("SuppRestoredTreesN"), len(rest))
    R.eq("macro SuppRestoredBranchLengthsN", num("SuppRestoredBranchLengthsN"),
         sum(as_int(r["branch_lengths_restored"], 0) for r in rest))
    R.eq("macro SuppShippedTreesN", num("SuppShippedTreesN"), len(shipped))
    R.eq("macro SuppPolySeedsN", num("SuppPolySeedsN"),
         len({r["seed"] for r in read_csv(WORKED_POLY)}))
    fo = defaultdict(dict)
    for r in read_csv(WORKED_FORDER):
        fo[r["partition"]][r["child"]] = r
    R.eq("macro SuppFamOrderMaxRDiffPct", mac["SuppFamOrderMaxRDiffPct"],
         "%.1f" % (100 * max(abs(as_float(v["last"]["r"]) - as_float(v["first"]["r"]))
                             / as_float(v["first"]["r"]) for v in fo.values())))
    R.eq("macro SuppFamOrderMaxEpsDiff", mac["SuppFamOrderMaxEpsDiff"],
         "%.3f" % max(abs(as_float(v["last"]["epsilon"]) - as_float(v["first"]["epsilon"]))
                      for v in fo.values()))
    R.eq("macro SuppFamOrderMaxDLL", mac["SuppFamOrderMaxDLL"],
         "%.1f" % max(abs(as_float(v["first"]["loglik_last_at_first_mle"])
                          - as_float(v["first"]["loglik"])) for v in fo.values()))



CATEGORY_NAMES_FIG = ("Vertebrates", "Plants", "Arthropods", "Other animals",
                      "Microbes & protists")


# --------------------------------------------------------------------------
# PROSELIT: no numeric literal in the revised manuscript text
# --------------------------------------------------------------------------
PROSE_V2 = [P("Paper", "shared", "manuscript_body_v2.tex"),
            P("Paper", "shared", "metadata_v2.tex"),
            P("Paper", "shared", "supplementary_body_v2.tex"),
            P("Paper", "shared", "supplementary_metadata_v2.tex")]
# Constructs whose digits are not reported quantities: citations and
# cross-references, file inputs, URLs, layout lengths and column specs.
PROSE_STRIP = [
    r"\\(?:cite[tp]?\*?|citeauthor|citeyear|ref|eqref|label|inputtable|input|AtlasFigure|"
    r"includegraphics|url|href|bibliography|bibliographystyle|"
    r"graphicspath)\s*(?:\[[^\]]*\]\s*)*\{[^{}]*\}",
    r"\\cmidrule\s*(?:\([^)]*\))?\s*\{[^{}]*\}",
    r"\\(?:newcommand|renewcommand|providecommand)\{\\[A-Za-z]+\}\[\d\]",
    r"#\d",
    r"\\(?:hspace|vspace|setlength|addtolength|setstretch|multicolumn|"
    r"arraystretch|resizebox|parbox|makebox|rule)\*?\s*(?:\{[^{}]*\}\s*){0,3}",
    r"\\renewcommand\{\\arraystretch\}\{[^{}]*\}",
    r"\\begin\{(?:longtable|tabular\*?|tabularx|table|figure|landscape)\}"
    r"\s*(?:\[[^\]]*\])?\s*(?:\{(?:[^{}]|\{[^{}]*\})*\})?",
    r"\[(?:width|height|scale|angle)=[^\]]*\]",
    r"\\(?:p|m|b)\{[^{}]*\}",
    r"\\\\\[[^\]]*\]",
    r"\b\d+(?:\.\d+)?\s*(?:pt|cm|mm|em|ex)\b",
    r"\\(?:[a-zA-Z]+)",            # macro names (letters only)
]
ID_MACRO_RE = re.compile(r"^\\newcommand\{\\Id[A-Za-z]+\}\{(.*)\}\s*$", re.M)


def _prose_scrub(text):
    text = strip_tex_comments(text)
    ids = ID_MACRO_RE.findall(text)
    text = ID_MACRO_RE.sub("", text)
    # math: only the structural 0 and 1 are allowed (mu = 0, 1 - epsilon)
    def math(m):
        return re.sub(r"(?<![\d.])[01](?![\d.])", "", m.group(0))
    text = re.sub(r"\$[^$]*\$", math, text)
    for pat in PROSE_STRIP:
        text = re.sub(pat, " ", text)
    return text, ids


def check_proselit(ctx):
    """The revised manuscript files write no number by hand: every digit that
    survives the removal of citations, references, file names, layout lengths,
    macro names and the structural 0 and 1 in mathematics is a literal and
    fails.  Identifiers that are not quantities (release and version names,
    supplement section labels) live only in \\Id... macros of metadata_v2.tex,
    each printed here for review.  A self-test plants a literal and requires
    that it is caught."""
    R.bump()
    probe, _ = _prose_scrub(r"Named coverage is \CZeroPct{} \% of 2{,}159{,}542 "
                            r"species \citep{x2019}.")
    if not re.search(r"\d", probe):
        R.fail("PROSELIT self-test: a planted literal was not caught")
    R.bump()
    clean, _ = _prose_scrub(r"Named coverage is \CZeroPct{} \% "
                            r"\citep{x2019} (Table~\ref{tab:a1}), $\mu = 0$.")
    if re.search(r"\d", clean):
        R.fail("PROSELIT self-test: a macro, citation, reference or the "
               "structural zero was reported as a literal: %r" % clean)
    for path in PROSE_V2:
        R.bump()
        if not os.path.exists(path):
            R.fail("%s does not exist" % rel(path))
            continue
        raw = open(path, encoding="utf-8").read()
        text, ids = _prose_scrub(raw)
        for i, ln in enumerate(text.splitlines(), 1):
            for m in re.finditer(r"\d[\d.,]*", ln):
                R.fail("%s: numeric literal %r in the prose (scrubbed line %d)"
                       % (rel(path), m.group(0), i), [ln.strip()[:160]])
        for v in ids:
            R.note("%s identifier (not a quantity): %s" % (rel(path), v))
        # every macro the prose uses is defined somewhere it can be read from
    body = open(PROSE_V2[0], encoding="utf-8").read() if os.path.exists(
        PROSE_V2[0]) else ""
    used = set(re.findall(r"\\([A-Z][A-Za-z]+)\b", strip_tex_comments(body)))
    defined = set()
    for f in GENERATED_MACRO_FILES + (PROSE_V2[1],):
        if os.path.exists(f):
            defined |= set(re.findall(r"\\(?:newcommand|providecommand)"
                                      r"\{\\([A-Za-z]+)\}", open(f).read()))
    generated_like = re.compile(r"^(?:C(?:Zero|One)|COneD|COneB|Num|Described|"
                                r"Sens|Lost|WE|OT|Crown|Cat|FMol|Rule|"
                                r"Unresolved|Labels|DT|SeedPlantsShare|Excl|"
                                r"TipsNotNamed|Unrep|Study|NMol|NGroups)")
    for name in sorted(used):
        if generated_like.match(name):
            R.bump()
            if name not in defined:
                R.fail("the manuscript uses \\%s, which no generated macro "
                       "file defines" % name)
    check_text_facts()


AI_DISCLOSURE = P("Paper", "MEE", "ai_disclosure.txt")
GBIF_UNMATCHED = P("Data", "external", "gbif_unmatched.csv")


def _expand_ids(text, ids):
    for k, v in sorted(ids.items(), key=lambda kv: -len(kv[0])):
        text = re.sub(r"\\%s(?:\{\})?(?![A-Za-z])" % k, lambda _m: v, text)
    return text


def check_text_facts():
    """Statements the text makes about its own numbers and companion files
    (round-3 findings): the label arithmetic of the Names paragraph closes;
    the AI statement filed with the submission is the Methods statement word
    for word; every code printed in Table 1 is decoded in its caption; the
    GBIF query date names every date on which a query was made."""
    mac = OrderedDict(_macros_in(NUMBER_MACROS))
    num = lambda k: int(mac[k].replace("{,}", ""))
    R.eq("Names: tip labels - tips not named = named species",
         num("LabelsEukN") - num("TipsNotNamedN"), num("CZeroN"))
    meta = open(PROSE_V2[1], encoding="utf-8").read()
    ids = dict(re.findall(r"^\\newcommand\{\\(Id\w+)\}\{(.*)\}\s*$", meta,
                          re.M))
    m = re.search(r"^\\newcommand\{\\PaperAIStatement\}\{(.*)\}\s*$", meta,
                  re.M)
    R.eq("AI statement macro present", bool(m), True)
    if m:
        want = _expand_ids(m.group(1), ids).replace("~", " ")
        got = ""
        if os.path.exists(AI_DISCLOSURE):
            lines = open(AI_DISCLOSURE, encoding="utf-8").read().split("\n")
            got = " ".join(x.strip() for x in lines[1:] if x.strip())
        R.eq("Paper/MEE/ai_disclosure.txt = the Methods AI statement", got, want)
    k = meta.find("\\newcommand{\\TabWhichTreeCaption}")
    cap = meta[k:meta.find("\n\\newcommand", k + 1)] if k >= 0 else ""
    codes = set()
    with open(P("Paper", "inputs", "table_which_tree.tex"), encoding="utf-8") as fh:
        for ln in fh:
            c = [x.strip() for x in ln.split("&")]
            if len(c) == 12:
                codes |= {c[8], c[9]}
    for code in sorted(codes - {"---"}):
        R.eq("Table 1 code %s decoded in the caption" % code,
             bool(re.search(r"(?<![A-Za-z])%s," % re.escape(code), cap)), True)
    dates = sorted({r["queried_on"] for r in read_csv(GBIF_UNMATCHED)})
    gd = ids.get("IdGbifDate", "")
    months = ["January", "February", "March", "April", "May", "June", "July",
              "August", "September", "October", "November", "December"]
    for d in dates:
        y, mo, dd = (int(x) for x in d.split("-"))
        R.eq("IdGbifDate names the query date %s" % d,
             bool(re.search(r"(?<!\d)%d(?!\d)" % dd, gd)) and months[mo - 1] in gd
             and str(y) in gd, True)


# --------------------------------------------------------------------------
# RELEASE: the public data release the paper names
# --------------------------------------------------------------------------
RELEASE_MAP = P("release", "path_map.csv")
RELEASE_FILES = P("release", "files")
RELEASE_TEX = PROSE_V2 + [P("Paper", "bioRxiv", "main.tex"),
                          P("Paper", "bioRxiv", "supplementary.tex")] + MEE_WRAPPERS
_RELEASE_PATHS = None


def release_paths():
    """{release path: (working path or None for a release-only file)}."""
    global _RELEASE_PATHS
    if _RELEASE_PATHS is None:
        out = OrderedDict()
        for r in read_csv(RELEASE_MAP):
            out[r["release_path"]] = r["working_path"]
        for d, dirs, files in os.walk(RELEASE_FILES):
            dirs[:] = [x for x in dirs if x != "__pycache__"]
            for f in files:
                if f != ".DS_Store" and not f.endswith(".pyc"):
                    out[os.path.relpath(os.path.join(d, f), RELEASE_FILES)
                        .replace(os.sep, "/")] = None
        _RELEASE_PATHS = out
    return _RELEASE_PATHS


def release_path_known(path):
    """A file of the release, or a directory holding at least one."""
    rp = release_paths()
    if path.endswith("/"):
        return any(k.startswith(path) for k in rp)
    return path in rp


def _check_release_readme():
    """Round F1: the numbers the release README prints equal the macros of the
    paper (the README is a release-only file, written by hand)."""
    path = os.path.join(RELEASE_FILES, "README.md")
    t = open(path, encoding="utf-8").read()
    mac = OrderedDict(_macros_in(NUMBER_MACROS))
    plain = lambda v: v.replace("{,}", ",")                    # noqa: E731
    m = re.search(r"The ([\d,]+) described eukaryotic species are divided into "
                  r"(\d+) groups; (\d+) of them", t)
    R.eq("release README: described species, groups, groups with a tree",
         m.groups() if m else None,
         (plain(mac["DescribedEuk"]), plain(mac["NumPartitionsEuk"]),
          plain(mac["NumCanonicalEuk"])))
    for lv, nm in (("Named", "CZero"), ("Molecular", "COne"), ("Dated", "COneD"),
                   ("Interval", "COneB")):
        m = re.search(r"^\| %s \|[^\n]*\| ([\d,]+) \| ([\d.]+) \|$" % lv, t, re.M)
        R.eq("release README level %s" % lv, m.groups() if m else None,
             (plain(mac[nm + "N"]), mac[nm + "Pct"]))
    m = re.search(r"The interval level is ([\d.]+)% when every archived", t)
    R.eq("release README: every archived interval counted",
         m.group(1) if m else None, mac["SensAllIntCOneBPct"])


def check_release(ctx):
    """The paper names data release \\IdDataRelease and prints paths in it.
    scripts/build_release.py writes the release from release/path_map.csv
    (one row per mapped file) and release/files/ (files only the release
    holds).  Every mapped working file exists and is listed once; every path
    the paper prints (\\nolinkurl, the terms table) is a file or directory of
    the release; every generated file the wrappers or bodies input from
    Paper/inputs/ is mapped into the release; the release's CITATION.cff
    carries the paper's title and the release version printed in the text."""
    rows = read_csv(RELEASE_MAP)
    _check_release_readme()
    seen = Counter(r["release_path"] for r in rows)
    R.eq("path_map.csv lists each release path once",
         sorted(k for k, n in seen.items() if n > 1), [])
    missing = [r["working_path"] for r in rows
               if not os.path.isfile(P(*r["working_path"].split("/")))]
    R.eq("every mapped working file exists", missing, [])
    R.eq("path_map.csv transforms are copy or gzip",
         sorted({r["transform"] for r in rows} - {"copy", "gzip"}), [])
    only = [k for k, v in release_paths().items() if v is None]
    R.eq("no release path is both mapped and release-only",
         sorted(set(only) & set(seen)), [])
    by_working = {r["working_path"]: r["release_path"] for r in rows}

    printed, inputs = set(), set()
    for path in RELEASE_TEX:
        if not os.path.exists(path):
            continue
        t = strip_tex_comments(open(path, encoding="utf-8").read())
        printed.update(re.findall(r"\\nolinkurl\{([^}]*)\}", t))
        for m in re.finditer(r"\\(?:input|inputtable)\{\.\./inputs/([^}]+)\}", t):
            inputs.add("Paper/inputs/" + m.group(1))
    for r in read_csv(SUPP_IN("supp_terms.csv")):
        if r["file"]:
            printed.add(r["file"])
    R.eq("the paper prints at least one release path", bool(printed), True)
    for p_ in sorted(printed):
        R.eq("printed path %s is in the release" % p_, release_path_known(p_), True)
    R.eq("the wrappers input generated files", bool(inputs), True)
    for w in sorted(inputs):
        R.eq("generated input %s is in the release" % w, w in by_working, True)

    cff = os.path.join(RELEASE_FILES, "CITATION.cff")
    meta = open(MANUSCRIPT_META, encoding="utf-8").read()
    title = re.search(r"\\newcommand\{\\PaperTitle\}\{(.*)\}\s*$", meta, re.M)
    rel_id = re.search(r"\\newcommand\{\\IdDataRelease\}\{v?([^}]*)\}", meta)
    R.eq("metadata_v2.tex defines PaperTitle and IdDataRelease",
         bool(title and rel_id), True)
    if title and rel_id and os.path.exists(cff):
        c = open(cff, encoding="utf-8").read()
        m = re.search(r"preferred-citation:.*?\n\s+title:\s*\"([^\"]*)\"", c, re.S)
        R.eq("CITATION.cff preferred-citation title = PaperTitle",
             m.group(1) if m else None, title.group(1))
        m = re.search(r"^version:\s*\"?([^\"\s]+)", c, re.M)
        R.eq("CITATION.cff version = IdDataRelease",
             m.group(1) if m else None, rel_id.group(1))


# --------------------------------------------------------------------------
# BIORXIV: the preprint wrappers
# --------------------------------------------------------------------------
BIORXIV_MAIN = P("Paper", "bioRxiv", "main.tex")
BIORXIV_SUPP = P("Paper", "bioRxiv", "supplementary.tex")
BIORXIV_REVIEW_ONLY = (r"\paperanontrue", r"\AnonDataLink", r"\PaperRunningHead",
                       r"\PaperDataReviewStatement", r"\linenumberstrue",
                       r"\PaperAbstractPoints", "main_table_rows")


def check_biorxiv(ctx):
    """The bioRxiv preprint is the named build of the revised article: the
    main wrapper reads manuscript_body_v2.tex and metadata_v2.tex, the
    Supplement wrapper supplementary_body_v2.tex; neither carries an element
    that only the journal review needs (anonymisation, running head, review
    link statement, line numbers switched on, the numbered abstract, the
    hand-maintained table); the main wrapper prints the abstract, the data
    and code statements and the author block.  The data statement names the
    repository, the release, its version DOI and the concept DOI; while the
    version DOI is not yet minted, a note says so."""
    m = strip_tex_comments(open(BIORXIV_MAIN, encoding="utf-8").read())
    s = strip_tex_comments(open(BIORXIV_SUPP, encoding="utf-8").read())
    for need in (r"\input{../shared/manuscript_body_v2.tex}",
                 r"\input{../shared/metadata_v2.tex}", r"\PaperAbstract",
                 r"\PaperDataAvailability", r"\PaperCodeAvailability",
                 r"\PaperAuthorName", r"\PaperAffiliation", r"\PaperORCID",
                 r"\PaperCompetingInterests", r"\PaperFunding"):
        R.eq("bioRxiv main.tex has %s" % need, need in m, True)
    for need in (r"\input{../shared/supplementary_body_v2.tex}",
                 r"\input{../shared/metadata_v2.tex}",
                 r"\input{../shared/supplementary_metadata_v2.tex}",
                 r"\PaperAuthorName"):
        R.eq("bioRxiv supplementary.tex has %s" % need, need in s, True)
    for bad in BIORXIV_REVIEW_ONLY:
        R.eq("bioRxiv main.tex has no %s" % bad, bad in m, False)
        R.eq("bioRxiv supplementary.tex has no %s" % bad, bad in s, False)
    meta = open(MANUSCRIPT_META, encoding="utf-8").read()
    da = re.search(r"\\newcommand\{\\PaperDataAvailability\}\{(.*)\}\s*$", meta, re.M)
    R.eq("metadata_v2.tex defines PaperDataAvailability", bool(da), True)
    if da:
        for need in ("https://github.com/franciscorichter/phylo-species-atlas",
                     "https://doi.org/10.5281/zenodo.20127157",
                     r"\IdDataRelease", r"\DataVersionDOI"):
            R.eq("the data statement names %s" % need, need in da.group(1), True)
    v = re.search(r"\\newcommand\{\\DataVersionDOI\}\{(.*)\}\s*$", meta, re.M)
    R.eq("metadata_v2.tex defines DataVersionDOI", bool(v), True)
    minted = re.search(r"^\\newif\\ifDataVersionMinted\\DataVersionMintedtrue",
                       meta, re.M) is not None
    R.eq("a version DOI is cited only once minted, and then is a Zenodo DOI",
         (not minted) or (v is not None and "doi.org/10.5281/zenodo." in v.group(1)),
         True)
    if not minted:
        R.note("the data statement cites tag \\IdDataRelease and the concept "
               "DOI; the version DOI is added (\\DataVersionMintedtrue) once "
               "Zenodo mints it")


CHECKS = [
    ("INPUTS", "every required file exists"),
    ("TREES", "re-measure the shipped .nwk files"),
    ("BRANCHLEN", "every shipped tree reproduces its declared reference"),
    ("FILE-TRUTH", "dated if and only if ultrametric, in every metadata file"),
    ("DECISIONS", "tips removed and merged, and chronograms chosen, re-derived from their registers"),
    ("NAMES", "Catalogue of Life name resolution re-derived"),
    ("UNION", "named species of all redistributed trees, recomputed"),
    ("COLDEN", "Catalogue of Life counts used as denominators"),
    ("OPENTREE", "Open Tree of Life comparison recomputed"),
    ("CROWN", "crown ages and their archived intervals"),
    ("CROWNREG", "one register of crown ages for text, tables and figures"),
    ("PREAUDIT", "superseded tree metadata described correctly"),
    ("INVARIANTS", "described/estimated/tips bounds, per partition"),
    ("BASIS", "the per-group table: one row per group, every column"),
    ("FMOL", "every molecular fraction re-measured from the source analysis"),
    ("HEADLINES", "every headline recomputed independently"),
    ("PERCENTS", "every percentage reaches a numerator and a denominator"),
    ("MONOTONICITY", "a stricter setting never increases coverage"),
    ("MOLTIPS", "molecular_tips.csv enumerates the measured tip sets"),
    ("C1COUNT", "molecular coverage counted as accepted names among molecular tips"),
    ("WORKED", "worked example fitted on the current trees"),
    ("GBIF", "named coverage with GBIF-resolved labels"),
    ("UNITS", "levels below named coverage counted in tips where tips are not listed"),
    ("CROSSLEDGER", "no field holds two different values"),
    ("NOLITERALS", "no coverage number hard-coded in a figure script"),
    ("NOABSPATHS", "no absolute home-directory path in code"),
    ("RETIRED", "no retired figure survives outside a disclosure"),
    ("DATASETS", "the dataset table re-derived from the basis"),
    ("FIGURES", "every manuscript figure exists and has a script"),
    ("NUMMACROS", "every prose number macro re-derived"),
    ("TABLES", "Tables 1-3 and the Figure 1 cascade re-derived"),
    ("PROSELIT", "no numeric literal in the revised manuscript text"),
    ("SUPPLEMENT", "the Supplement's sections, tables and macros re-derived"),
    ("RELEASE", "every path the paper prints is in the data release"),
    ("BIORXIV", "the bioRxiv wrappers build the revised, named preprint"),
]


def run(args):
    global COMMON

    sys.stdout.write("verify_atlas -- acceptance test for %s\n" % rel(ROOT))
    sys.stdout.write("%s\n" % ("-" * 78))

    if not run_check("INPUTS", "every required file exists", check_inputs):
        sys.stdout.write("\nRequired inputs are missing; the remaining checks "
                         "cannot run.\n")
        return 1

    COMMON = load_audit_common()
    ctx = Ctx(quick=args.quick)
    ctx.ledger_rows = {r.get("filename"): r
                       for r in read_csv(TREE_PROVENANCE, required=False)}
    ctx.derived = OrderedDict()

    run_check("TREES", "re-measure the shipped .nwk files", measure_trees, ctx)
    run_check("BRANCHLEN", "every shipped tree reproduces its declared reference",
              check_branchlen, ctx)
    run_check("FILE-TRUTH", "dated if and only if ultrametric, in every metadata file",
              check_file_truth, ctx)
    run_check("DECISIONS", "tips removed and merged, and chronograms chosen, re-derived from their registers",
              check_decisions, ctx)
    run_check("NAMES", "Catalogue of Life name resolution re-derived",
              check_names, ctx)
    run_check("UNION", "named species of all redistributed trees, recomputed",
              check_union, ctx)
    run_check("COLDEN", "Catalogue of Life counts used as denominators",
              check_colden, ctx)
    run_check("OPENTREE", "Open Tree of Life comparison recomputed",
              check_opentree, ctx)
    run_check("CROWN", "crown ages and their archived intervals",
              check_crown, ctx)
    run_check("CROWNREG", "one register of crown ages for text, tables and figures", check_crownreg, ctx)
    run_check("PREAUDIT", "ground truth describes the pre-audit ledgers",
              check_preaudit, ctx)
    run_check("INVARIANTS", "described/estimated/tips bounds, per partition",
              check_invariants, ctx)
    run_check("BASIS", "the per-group table: one row per group, every column",
              check_basis, ctx)
    run_check("FMOL", "every molecular fraction re-measured from the source analysis",
              check_fmol, ctx)

    ctx.my_rows = my_basis(ctx)

    run_check("HEADLINES", "every headline recomputed independently",
              check_headlines, ctx)
    run_check("PERCENTS", "every percentage reaches a numerator and denominator",
              check_percents, ctx)
    run_check("MONOTONICITY", "a stricter setting never increases coverage",
              check_monotonicity, ctx)
    run_check("MOLTIPS", "molecular_tips.csv enumerates the measured tip sets",
              check_moltips, ctx)
    run_check("C1COUNT", "molecular coverage counted as accepted names among molecular tips",
              check_c1count, ctx)
    run_check("WORKED", "worked example fitted on the current trees",
              check_worked, ctx)
    run_check("GBIF", "named coverage with GBIF-resolved labels",
              check_gbif, ctx)
    run_check("UNITS", "levels below named coverage counted in tips where tips are not listed",
              check_units, ctx)
    run_check("CROSSLEDGER", "no field holds two different values",
              check_crossledger, ctx)
    run_check("NOLITERALS", "no coverage number hard-coded in a figure script",
              check_noliterals, ctx)
    run_check("NOABSPATHS", "no absolute home-directory path in code",
              check_noabspaths, ctx)
    run_check("RETIRED", "no retired figure survives outside a disclosure",
              check_retired, ctx)
    run_check("DATASETS", "the dataset table re-derived from the basis",
              check_datasets, ctx)
    run_check("FIGURES", "every manuscript figure exists and has a script",
              check_figures, ctx)
    run_check("NUMMACROS", "every prose number macro re-derived",
              check_numbermacros, ctx)
    run_check("TABLES", "Tables 1-3 and the Figure 1 cascade re-derived",
              check_tables, ctx)
    run_check("PROSELIT", "no numeric literal in the revised manuscript text",
              check_proselit, ctx)
    run_check("SUPPLEMENT", "the Supplement's sections, tables and macros "
              "re-derived", check_supplement, ctx)
    run_check("RELEASE", "every path the paper prints is in the data release",
              check_release, ctx)
    run_check("BIORXIV", "the bioRxiv wrappers build the revised, named preprint",
              check_biorxiv, ctx)

    # ---- the report ---------------------------------------------------
    sys.stdout.write("%s\n" % ("-" * 78))
    notes = [(c["id"], n) for c in R.checks.values() for n in c["notes"]]
    if notes:
        sys.stdout.write("\nNOTES -- not failures; things a reader must see\n")
        for cid, n in notes:
            sys.stdout.write("  [%s] %s\n" % (cid, n["message"]))
            for d in n.get("detail", [])[:8]:
                sys.stdout.write("        %s\n" % d)

    total = sum(c["assertions"] for c in R.checks.values())
    sys.stdout.write("\n%d assertions, %d failures, %d notes\n"
                     % (total, R.n_failures, R.n_notes))

    if args.json:
        payload = OrderedDict([
            ("root", rel(ROOT)),
            ("quick", bool(args.quick)),
            ("assertions", total),
            ("failures", R.n_failures),
            ("notes", R.n_notes),
            ("checks", list(R.checks.values())),
        ])
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, ensure_ascii=False, default=str)
            fh.write("\n")
        sys.stdout.write("report written to %s\n" % args.json)

    if R.n_failures:
        sys.stdout.write("\nFAILED. The deposit is NOT acceptable as it "
                         "stands. Every failure above names the file, the "
                         "value the deposit carries and the value this script "
                         "measured.\n")
        return 1
    if args.quick:
        sys.stdout.write("\nPASSED under --quick, which does NOT re-measure "
                         "every shipped tree. Run without --quick before "
                         "treating this as an acceptance result.\n")
        return 0
    sys.stdout.write("\nPASSED. Every headline in "
                     "Paper/inputs/headline_numbers.json was recomputed from "
                     "the truth files and the shipped trees and agrees.\n")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Acceptance test for the phylo-atlas deposit. "
                    "Non-zero exit on any failure.")
    ap.add_argument("--quick", action="store_true",
                    help="re-measure only the partition-canonical trees "
                         "(not a valid acceptance run)")
    ap.add_argument("--json", metavar="PATH",
                    help="write a machine-readable report")
    ap.add_argument("--list", action="store_true",
                    help="list the checks and exit")
    args = ap.parse_args(argv)

    if args.list:
        for cid, title in CHECKS:
            sys.stdout.write("  %-13s %s\n" % (cid, title))
        return 0

    try:
        return run(args)
    except Exception:                                      # noqa: BLE001
        sys.stderr.write("verify_atlas aborted:\n%s\n" % traceback.format_exc())
        return 2


if __name__ == "__main__":
    # The reference parser is recursive; the largest shipped tree has ~356k
    # leaves.  Run on a thread with a large stack where the platform allows it.
    _status = [2]

    def _target():
        _status[0] = main()

    try:
        threading.stack_size(256 * 1024 * 1024)
        t = threading.Thread(target=_target)
        t.start()
        t.join()
    except (ValueError, RuntimeError, threading.ThreadError):
        _status[0] = main()
    sys.exit(_status[0])
